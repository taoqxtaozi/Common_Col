#!/bin/bash
set -euo pipefail
shopt -s nullglob

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage() {
    cat <<'EOF'
Usage:
  bash RunPipelineUseThis.sh --input /path/to/aln1.maf /path/to/aln2.maf ... --output consensus.maf --threads 60 --common_workers 15 [options]

Required arguments:
  --input                 Input MAF files (space-separated). All input MAF files must be in the same directory.
  --output                Final consensus MAF file. If a relative path is given, it will be created in the input MAF directory.
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
  * This script may be run from any directory.
  * The code files are expected to be located in the same directory as this script:
      - seperate_maffile_by_refchrom.py
      - seperate_maffile_by_index.py
      - Extract_consensus_multi_alns.py
      - sep_file/1.sh
  * All intermediate files and directories will be created in the directory containing the input MAF files.
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

############################
# 0. Default parameter values
############################
raw_inputs=()
inputs=()
threads=""
output=""
output_abs=""
work_dir=""
sep_length=218505
chrom_length_threshold=15819469
global_num=40
separate_workers=4
common_workers=""
last_dir=""

############################
# 1. Parse command-line arguments
############################
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                raw_inputs+=("$1")
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
# 2. Check parameters
############################
if [[ ${#raw_inputs[@]} -eq 0 ]]; then
    echo "Error: --input requires at least one MAF file."
    exit 1
fi

if [[ -z "$output" ]]; then
    echo "Error: --output is required."
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
# 3. Resolve input paths and working directory
############################
for file in "${raw_inputs[@]}"
do
    abs_file=$(readlink -f "$file")
    if [[ ! -f "$abs_file" ]]; then
        echo "Error: input file not found: $file"
        exit 1
    fi

    file_dir=$(dirname "$abs_file")

    if [[ -z "$work_dir" ]]; then
        work_dir="$file_dir"
    elif [[ "$file_dir" != "$work_dir" ]]; then
        echo "Error: all input MAF files must be in the same directory."
        echo "  First input directory: $work_dir"
        echo "  Conflicting directory: $file_dir"
        exit 1
    fi

    inputs+=("$abs_file")
done

if [[ "$output" = /* ]]; then
    output_abs="$output"
else
    output_abs="$work_dir/$output"
fi

mkdir -p "$(dirname "$output_abs")"

############################
# 4. Check required files
############################
if [[ ! -f "$script_dir/seperate_maffile_by_refchrom.py" ]]; then
    echo "Error: seperate_maffile_by_refchrom.py not found in script directory."
    exit 1
fi

if [[ ! -f "$script_dir/seperate_maffile_by_index.py" ]]; then
    echo "Error: seperate_maffile_by_index.py not found in script directory."
    exit 1
fi

if [[ ! -f "$script_dir/Extract_consensus_multi_alns.py" ]]; then
    echo "Error: Extract_consensus_multi_alns.py not found in script directory."
    exit 1
fi

if [[ ! -f "$script_dir/sep_file/1.sh" ]]; then
    echo "Error: sep_file/1.sh not found in script directory."
    exit 1
fi

############################
# 5. Split each input MAF by reference sequence ID in parallel
############################
last_dir=$(basename "${inputs[${#inputs[@]}-1]}" .maf)

split_input_cmds=()

for file in "${inputs[@]}"
do
    dir=$(basename "$file" .maf)
    split_input_cmds+=("mkdir -p \"$work_dir/$dir\" && python \"$script_dir/seperate_maffile_by_refchrom.py\" \"$file\" \"$work_dir/$dir\"")
done

run_parallel_commands "$threads" "${split_input_cmds[@]}"

############################
# 6. Prepare intermediate directories in the working directory
############################
mkdir -p "$work_dir/common_command_files"
mkdir -p "$work_dir/common_output_files"
mkdir -p "$work_dir/sep"

############################
# 7. Use the directory generated from the last input file as the template,
#    iterate through each chromosome MAF in that directory,
#    and determine whether further splitting is needed
############################
cd "$work_dir/sep" || exit 1

for i in "$work_dir/${last_dir}"/*.maf
do
    [ -e "$i" ] || continue

    head_chrom=$(basename "$i" .maf)
    length=$(awk '$1=="s"{print $6; exit}' "$i")

    mkdir -p "$work_dir/common_command_files/${head_chrom}"
    mkdir -p "$work_dir/common_output_files/${head_chrom}"

    if (( length > chrom_length_threshold )); then
        subdir="$work_dir/sep/sep_${head_chrom}"
        mkdir -p "$subdir" "$work_dir/sep/${head_chrom}"
        cd "$subdir" || exit 1

        num=$(( (length + sep_length - 1) / sep_length ))

        for n in $(seq 1 "$num")
        do
            start=$(( (n - 1) * sep_length ))
            end=$(( n * sep_length ))
            if (( end > length )); then
                end=$length
            fi

            cat > "${n}.sh" <<EOF
#!/bin/bash
set -euo pipefail
mkdir -p "../${head_chrom}"
python "$script_dir/seperate_maffile_by_index.py" "../${head_chrom}.maf" "${n}" "${start}" "${end}" "../${head_chrom}"
EOF
            chmod +x "${n}.sh"
        done

        cd "$work_dir/sep" || exit 1

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
                input_args="${input_args} \"$work_dir/${d}/${head_chrom}/${n}_${start}_${end}.maf\""
            done

            echo "python \"$script_dir/Extract_consensus_multi_alns.py\" --input${input_args} --output \"$work_dir/common_output_files/${head_chrom}/part${n}.maf\" --global_num ${global_num} --global_start ${start} --global_end ${end} --separate_workers ${separate_workers} --common_workers ${common_workers}" \
                > "$work_dir/common_command_files/${head_chrom}/${n}.sh"
        done

    else
        input_args=""
        for file in "${inputs[@]}"
        do
            d=$(basename "$file" .maf)
            input_args="${input_args} \"$work_dir/${d}/${head_chrom}.maf\""
        done

        echo "python \"$script_dir/Extract_consensus_multi_alns.py\" --input${input_args} --output \"$work_dir/common_output_files/${head_chrom}/part1.maf\" --global_num ${global_num} --global_start 0 --global_end ${length} --separate_workers ${separate_workers} --common_workers ${common_workers}" \
            > "$work_dir/common_command_files/${head_chrom}/1.sh"
    fi
done

cd "$work_dir" || exit 1

############################
# 8. Copy the generated sep_<chrom> directories under sep/
#    into each input directory
############################
for file in "${inputs[@]}"
do
    dir=$(basename "$file" .maf)

    for sep_dir in "$work_dir"/sep/*
    do
        [ -d "$sep_dir" ] || continue
        base_sep=$(basename "$sep_dir")
        rm -rf "$work_dir/$dir/$base_sep"
        cp -r "$sep_dir" "$work_dir/$dir"/
    done
done

############################
# 9. Execute the chunking scripts in parallel
#    Split jobs are treated as ~1 CPU jobs, so up to threads jobs may run in parallel.
############################
split_job_cmds=()

for file in "${inputs[@]}"
do
    dir=$(basename "$file" .maf)

    for sep_dir in "$work_dir/$dir"/sep_*
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
# 10. Execute the commands in common_command_files in parallel
#     Maximum parallel consensus jobs = floor(threads / common_workers)
############################
max_common_jobs=$(( threads / common_workers ))
if (( max_common_jobs < 1 )); then
    max_common_jobs=1
fi

common_job_cmds=()

for d in "$work_dir"/common_command_files/*/
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
# 11. First merge part*.maf for each chromosome into <chrom>.maf
############################
cd "$work_dir/common_output_files" || exit 1

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

cd "$work_dir" || exit 1

############################
# 12. Then merge all chromosome MAFs into the final output
############################
: > "$output_abs"

merged_chrom_files=( "$work_dir"/common_output_files/*.maf )
if [ ${#merged_chrom_files[@]} -gt 0 ]; then
    first=1
    while IFS= read -r f
    do
        if (( first == 1 )); then
            cat "$f" >> "$output_abs"
            first=0
        else
            awk 'BEGIN{p=0} !p && /^#/ {next} {p=1; print}' "$f" >> "$output_abs"
        fi
    done < <(printf '%s\n' "${merged_chrom_files[@]}" | sort -V)
fi

############################
# 13. Delete intermediate directories, keeping only the input MAF files
#     and the final output
############################
for file in "${inputs[@]}"
do
    dir=$(basename "$file" .maf)
    rm -rf "$work_dir/$dir"
done

rm -rf "$work_dir/sep" "$work_dir/common_command_files" "$work_dir/common_output_files"

echo "Pipeline finished. Final output: $output_abs"