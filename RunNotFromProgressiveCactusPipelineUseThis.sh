#!/bin/bash
set -euo pipefail
shopt -s nullglob

usage() {
    cat <<'EOF'
Usage:
  bash run_pipeline2.sh --input aln1.maf aln2.maf ... --reference REF --output consensus.maf --threads 60 --common_workers 15 [options]

Required arguments:
  --input                 Input MAF files (space-separated).
  --reference             Reference taxon name used to reorder blocks and split files.
  --output                Final consensus MAF file.
  --threads               Total CPU budget available on the machine.
  --common_workers        CPUs used by one consensus-extraction program.

Optional arguments:
  --sep_length            Chunk length for long chromosome files. Default: 218505
  --chrom_length_threshold
                          Chromosome files longer than this threshold will be split. Default: 15819469
  --global_num            Number of internal subsegments used by one consensus-extraction program. Default: 40
  --separate_workers      CPUs used by the internal separate() stage of one consensus-extraction program. Default: 4
  -h, --help              Show this help message and exit.

Notes:
  * This script is expected to be run in the same directory as:
      - reorder_maf_by_reference.py
      - split_maf_by_refseq.py
      - Extract_consensus_multi_alns.py
      - sep_file/1.sh
  * threads must be >= common_workers.
  * threads is the total CPU budget on the machine.
  * common_workers is the CPU usage of one consensus-extraction program.
  * If threads=20, a practical choice for common_workers is often 10 or 15.
  * The maximum number of consensus jobs run in parallel is floor(threads / common_workers).
  * The long-chromosome splitting scripts are treated as ~1 CPU jobs, so up to threads of them may run in parallel.
EOF
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

run_parallel_commands() {
    local max_jobs="$1"
    shift
    local -a cmds=("$@")
    local -a pids=()
    local fail=0

    if (( ${#cmds[@]} == 0 )); then
        return 0
    fi

    if ! is_positive_integer "$max_jobs"; then
        max_jobs=1
    fi

    for cmd in "${cmds[@]}"; do
        while (( $(jobs -rp | wc -l) >= max_jobs )); do
            sleep 1
        done

        bash -euo pipefail -c "$cmd" &
        pids+=("$!")
    done

    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            fail=1
        fi
    done

    if (( fail != 0 )); then
        echo "Error: at least one parallel job failed."
        return 1
    fi
}

preprocess_one_input() {
    local file="$1"
    local reference="$2"

    local pre dir p
    pre=$(basename "$file" .maf)
    dir="$pre"

    # Step 5: remove duplicate sequences
    conda run -n py2 mafDuplicateFilter --maf "$file" > "${pre}_noDup.maf"

    # Step 6: reorder taxa and discard blocks without the reference
    python reorder_maf_by_reference.py \
        --input "${pre}_noDup.maf" \
        --output "${pre}_noDup_ordered.maf" \
        --reference "$reference"

    # Step 7: split by reference sequence ID
    mkdir -p "$dir"
    python split_maf_by_refseq.py \
        --input "${pre}_noDup_ordered.maf" \
        --output_dir "$dir" \
        --reference "$reference"

    # Step 8: normalize reference strand and coordinate order
    (
        cd "$dir" || exit 1
        for f in *.maf
        do
            [ -e "$f" ] || continue
            p=$(basename "$f" .maf)

            conda run -n py2 mafStrander --maf "$f" --seq "$p" --strand + > "${p}_positive.maf"
            conda run -n py2 mafSorter --maf "${p}_positive.maf" --seq "$p" > "${p}_positive_sorted.maf"

            rm -f "${p}_positive.maf" "$f"
            mv "${p}_positive_sorted.maf" "$f"
        done
    )
}
export -f preprocess_one_input

############################
# 0. 参数默认值
# 0. Default parameter values
############################
inputs=()
threads=""
output=""
reference=""
sep_length=218505
chrom_length_threshold=15819469
global_num=40
separate_workers=4
common_workers=""
last_dir=""

############################
# 1. 解析命令行参数
# 1. Parse command-line arguments
############################
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                inputs+=("$1")
                shift
            done
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        --output)
            output="$2"
            shift 2
            ;;
        --reference)
            reference="$2"
            shift 2
            ;;
        --sep_length)
            sep_length="$2"
            shift 2
            ;;
        --chrom_length_threshold)
            chrom_length_threshold="$2"
            shift 2
            ;;
        --global_num)
            global_num="$2"
            shift 2
            ;;
        --separate_workers)
            separate_workers="$2"
            shift 2
            ;;
        --common_workers)
            common_workers="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown parameter: $1"
            echo
            usage
            exit 1
            ;;
    esac
done

############################
# 2. 参数检查
# 2. Check parameters
############################
if [[ ${#inputs[@]} -eq 0 ]]; then
    echo "Error: --input requires at least one MAF file."
    exit 1
fi

if [[ -z "$output" ]]; then
    echo "Error: --output is required."
    exit 1
fi

if [[ -z "$reference" ]]; then
    echo "Error: --reference is required."
    exit 1
fi

if [[ -z "$threads" ]]; then
    echo "Error: --threads is required."
    exit 1
fi

if [[ -z "$common_workers" ]]; then
    echo "Error: --common_workers is required."
    exit 1
fi

for value_name in threads common_workers sep_length chrom_length_threshold global_num separate_workers; do
    value="${!value_name}"
    if ! is_positive_integer "$value"; then
        echo "Error: --${value_name} must be a positive integer."
        exit 1
    fi
done

if (( threads < common_workers )); then
    echo "Error: --threads must be >= --common_workers."
    exit 1
fi

############################
# 3. 基础文件检查
# 3. Check required files
############################
if [[ ! -f "reorder_maf_by_reference.py" ]]; then
    echo "Error: 当前目录下找不到 reorder_maf_by_reference.py"
    exit 1
fi

if [[ ! -f "split_maf_by_refseq.py" ]]; then
    echo "Error: 当前目录下找不到 split_maf_by_refseq.py"
    exit 1
fi

if [[ ! -f "Extract_consensus_multi_alns.py" ]]; then
    echo "Error: 当前目录下找不到 Extract_consensus_multi_alns.py"
    exit 1
fi

if [[ ! -f "sep_file/1.sh" ]]; then
    echo "Error: 当前目录下找不到 sep_file/1.sh"
    exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "Error: 找不到 conda"
    exit 1
fi

if ! conda env list | awk 'NR>2 {print $1}' | grep -qx "py2"; then
    echo "Error: conda environment 'py2' does not exist."
    exit 1
fi

if ! conda run -n py2 bash -c 'command -v mafDuplicateFilter >/dev/null 2>&1'; then
    echo "Error: mafDuplicateFilter not found in conda environment 'py2'."
    exit 1
fi

if ! conda run -n py2 bash -c 'command -v mafStrander >/dev/null 2>&1'; then
    echo "Error: mafStrander not found in conda environment 'py2'."
    exit 1
fi

if ! conda run -n py2 bash -c 'command -v mafSorter >/dev/null 2>&1'; then
    echo "Error: mafSorter not found in conda environment 'py2'."
    exit 1
fi

############################
# 4. Preprocess each input MAF in parallel
#    For each input file:
#    1) remove duplicate sequences
#    2) reorder taxa so that the reference is first
#    3) split by reference sequence ID
#    4) normalize reference strand and coordinate order
############################
last_dir=$(basename "${inputs[${#inputs[@]}-1]}" .maf)

preprocess_cmds=()

for file in "${inputs[@]}"
do
    preprocess_cmds+=("preprocess_one_input \"$file\" \"$reference\"")
done

run_parallel_commands "$threads" "${preprocess_cmds[@]}"
############################
# 5. 准备中间目录
# 5. Prepare intermediate directories
############################
mkdir -p common_command_files
mkdir -p common_output_files
mkdir -p sep

############################
# 6. 用最后一个输入拆出来的目录作为模板，
#     遍历其中每条染色体 maf，判断是否需要再切块
# 6. Use the directory generated from the last input file as the template,
#     iterate through each chromosome MAF in that directory,
#     and determine whether further splitting is needed
############################
cd sep || exit 1

for i in ../"${last_dir}"/*.maf
do
    [ -e "$i" ] || continue

    head_chrom=$(basename "$i" .maf)
    length=$(awk '$1=="s"{print $6; exit}' "$i")

    mkdir -p ../common_command_files/"${head_chrom}"
    mkdir -p ../common_output_files/"${head_chrom}"

    if (( length > chrom_length_threshold )); then
        subdir="sep_${head_chrom}"
        mkdir -p "$subdir" "$head_chrom"
        cd "$subdir" || exit 1

        num=$(( (length + sep_length - 1) / sep_length ))

        for n in $(seq 1 "$num")
        do
            cp ../../sep_file/1.sh "${n}.sh"

            start=$(( (n - 1) * sep_length ))
            end=$(( n * sep_length ))
            if (( end > length )); then
                end=$length
            fi

            sed -i \
                -e "s#mkdir ../ref1#mkdir -p ../${head_chrom}#" \
                -e "s#../GALGA\.1\.maf 1 0 10925261 ../ref1#../${head_chrom}.maf ${n} ${start} ${end} ../${head_chrom}#" \
                "${n}.sh"
        done

        cd ..

        for n in $(seq 1 "$num")
        do
            start=$(( (n - 1) * sep_length ))
            end=$(( n * sep_length ))
            if (( end > length )); then
                end=$length
            fi

            input_args=""
            for file in "${inputs[@]}"
            do
                d=$(basename "$file" .maf)
                input_args="${input_args} ../../${d}/${head_chrom}/${n}_${start}_${end}.maf"
            done

            echo "python ../../Extract_consensus_multi_alns.py --input${input_args} --output ../../common_output_files/${head_chrom}/part${n}.maf --global_num ${global_num} --global_start ${start} --global_end ${end} --separate_workers ${separate_workers} --common_workers ${common_workers}" \
                > ../common_command_files/"${head_chrom}"/"${n}.sh"
        done

    else
        input_args=""
        for file in "${inputs[@]}"
        do
            d=$(basename "$file" .maf)
            input_args="${input_args} ../../${d}/${head_chrom}.maf"
        done

        echo "python ../../Extract_consensus_multi_alns.py --input${input_args} --output ../../common_output_files/${head_chrom}/part1.maf --global_num ${global_num} --global_start 0 --global_end ${length} --separate_workers ${separate_workers} --common_workers ${common_workers}" \
            > ../common_command_files/"${head_chrom}"/1.sh
    fi
done

cd ..

############################
# 7. 把 sep/ 下生成的 sep_<chrom> 目录复制到每个输入目录中
# 7. Copy the generated sep_<chrom> directories under sep/
#     into each input directory
############################
for file in "${inputs[@]}"
do
    dir=$(basename "$file" .maf)

    for sep_dir in sep/*
    do
        [ -d "$sep_dir" ] || continue
        base_sep=$(basename "$sep_dir")
        rm -rf "$dir/$base_sep"
        cp -r "$sep_dir" "$dir"/
    done
done

############################
# 8. 在每个输入目录里执行切块脚本
# 8. Execute the chunking scripts in each input directory
############################
split_job_cmds=()

for file in "${inputs[@]}"
do
    dir=$(basename "$file" .maf)

    for sep_dir in "$dir"/sep_*
    do
        [ -d "$sep_dir" ] || continue
        for f in "$sep_dir"/*.sh
        do
            [ -e "$f" ] || continue
            fname=$(basename "$f")
            split_job_cmds+=("cd \"$sep_dir\" && bash \"$fname\"")
        done
    done
done

run_parallel_commands "$threads" "${split_job_cmds[@]}"

############################
# 9. 执行 common_command_files 里的命令
# 9. Execute the commands in common_command_files
############################
max_common_jobs=$(( threads / common_workers ))
if (( max_common_jobs < 1 )); then
    max_common_jobs=1
fi

common_job_cmds=()

for d in common_command_files/*/
do
    [ -d "$d" ] || continue
    for f in "$d"*.sh
    do
        [ -e "$f" ] || continue
        fname=$(basename "$f")
        common_job_cmds+=("cd \"$d\" && bash \"$fname\"")
    done
done

run_parallel_commands "$max_common_jobs" "${common_job_cmds[@]}"

############################
# 10. 先把每条染色体的 part*.maf 合并成 <chrom>.maf
# 10. First merge part*.maf for each chromosome into <chrom>.maf
############################
cd common_output_files || exit 1

for d in */
do
    dir="${d%/}"
    out_file="${dir}.maf"
    : > "$out_file"

    part_files=( "$dir"/part*.maf )
    [ ${#part_files[@]} -gt 0 ] || continue

    first=1
    while IFS= read -r f
    do
        if (( first == 1 )); then
            cat "$f" >> "$out_file"
            first=0
        else
            awk 'BEGIN{p=0} !p && /^#/ {next} {p=1; print}' "$f" >> "$out_file"
        fi
    done < <(printf '%s\n' "${part_files[@]}" | sort -V)
done

cd ..

############################
# 11. 再把所有染色体 maf 合并成最终输出
# 11. Then merge all chromosome MAFs into the final output
############################
: > "$output"

merged_chrom_files=( common_output_files/*.maf )
if [ ${#merged_chrom_files[@]} -gt 0 ]; then
    first=1
    while IFS= read -r f
    do
        if (( first == 1 )); then
            cat "$f" >> "$output"
            first=0
        else
            awk 'BEGIN{p=0} !p && /^#/ {next} {p=1; print}' "$f" >> "$output"
        fi
    done < <(printf '%s\n' "${merged_chrom_files[@]}" | sort -V)
fi

############################
# 12. 删除中间文件夹和中间 maf，只保留输入 maf 和最终 output
# 12. Delete intermediate directories and intermediate MAF files,
#     keeping only the input MAF files and the final output
############################
for file in "${inputs[@]}"
do
    dir=$(basename "$file" .maf)
    rm -rf "$dir"
    rm -f "${dir}_noDup.maf"
    rm -f "${dir}_noDup_ordered.maf"
done

rm -rf sep common_command_files common_output_files

echo "Pipeline finished. Final output: $output"