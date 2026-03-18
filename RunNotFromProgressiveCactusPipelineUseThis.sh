#!/bin/bash
set -euo pipefail
shopt -s nullglob

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage() {
    cat <<'EOF'
Usage:
  bash RunNotFromProgressiveCactusPipelineUseThis.sh --input /path/to/aln1.maf /path/to/aln2.maf ... --reference REF --pre output_prefix --threads 60 --common_workers 15 [options]

Required arguments:
  --input                 Input MAF files (space-separated). All input MAF files must be in the same directory.
  --reference             Reference taxon name used to reorder blocks and used as the first taxon in FASTA output.
  --pre                   Output prefix. Two files will be generated: <pre>.maf and <pre>.fasta
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
      - reorder_maf_by_reference.py
      - split_maf_by_refseq.py
      - seperate_maffile_by_index.py
      - Extract_consensus_multi_alns.py
      - maf_to_concat_fasta.py
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

preprocess_one_input() {
    local file="$1"
    local reference="$2"
    local script_dir="$3"
    local work_dir="$4"

    local pre dir p
    pre=$(basename "$file" .maf)
    dir="$work_dir/$pre"

    # Step 5: remove duplicate sequences
    conda run -n py2 mafDuplicateFilter --maf "$file" > "$work_dir/${pre}_noDup.maf"

    # Step 6: reorder taxa and discard blocks without the reference
    python "$script_dir/reorder_maf_by_reference.py" \
        --input "$work_dir/${pre}_noDup.maf" \
        --output "$work_dir/${pre}_noDup_ordered.maf" \
        --reference "$reference"

    # Step 7: split by reference sequence ID
    mkdir -p "$dir"
    python "$script_dir/split_maf_by_refseq.py" \
        --input "$work_dir/${pre}_noDup_ordered.maf" \
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
# 0. Default parameter values
############################
raw_inputs=()
inputs=()
threads=""
pre=""
output_prefix_abs=""
output_maf_abs=""
output_fasta_abs=""
reference=""
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
        --pre)
            pre="$2"
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
# 2. Check parameters
############################
if [[ ${#raw_inputs[@]} -eq 0 ]]; then
    echo "Error: --input requires at least one MAF file."
    exit 1
fi

if [[ -z "$pre" ]]; then
    echo "Error: --pre is required."
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

if [[ "$pre" = /* ]]; then
    output_prefix_abs="$pre"
else
    output_prefix_abs="$work_dir/$pre"
fi

output_maf_abs="${output_prefix_abs}.maf"
output_fasta_abs="${output_prefix_abs}.fasta"

mkdir -p "$(dirname "$output_prefix_abs")"

############################
# 4. Check required files
############################
if [[ ! -f "$script_dir/reorder_maf_by_reference.py" ]]; then
    echo "Error: reorder_maf_by_reference.py not found in script directory."
    exit 1
fi

if [[ ! -f "$script_dir/split_maf_by_refseq.py" ]]; then
    echo "Error: split_maf_by_refseq.py not found in script directory."
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

if [[ ! -f "$script_dir/maf_to_concat_fasta.py" ]]; then
    echo "Error: maf_to_concat_fasta.py not found in script directory."
    exit 1
fi

if [[ ! -f "$script_dir/sep_file/1.sh" ]]; then
    echo "Error: sep_file/1.sh not found in script directory."
    exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda not found."
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
# 5. Preprocess each input MAF in parallel
#    1) remove duplicate sequences
#    2) reorder taxa so that the reference is first
#    3) split by reference sequence ID
#    4) normalize reference strand and coordinate order
############################
last_dir=$(basename "${inputs[${#inputs[@]}-1]}" .maf)

preprocess_cmds=()

for file in "${inputs[@]}"
do
    preprocess_cmds+=("preprocess_one_input \"$file\" \"$reference\" \"$script_dir\" \"$work_dir\"")
done

run_parallel_commands "$threads" "${preprocess_cmds[@]}"

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
# 12. Then merge all chromosome MAFs into the final MAF output
############################
: > "$output_maf_abs"

merged_chrom_files=( "$work_dir"/common_output_files/*.maf )
if [ ${#merged_chrom_files[@]} -gt 0 ]; then
    first=1
    while IFS= read -r f
    do
        if (( first == 1 )); then
            cat "$f" >> "$output_maf_abs"
            first=0
        else
            awk 'BEGIN{p=0} !p && /^#/ {next} {p=1; print}' "$f" >> "$output_maf_abs"
        fi
    done < <(printf '%s\n' "${merged_chrom_files[@]}" | sort -V)
fi

############################
# 13. Convert the final MAF to concatenated FASTA
############################
taxa_tmp=$(mktemp "$work_dir/taxa_list.XXXXXX")

if ! LC_ALL=C awk -v ref="$reference" '
$1=="s"{
    split($2,a,".")
    sp=a[1]
    if (sp == ref) {
        seen_ref=1
    } else {
        seen[sp]=1
    }
}
END{
    if (!seen_ref) exit 2
    for (sp in seen) print sp
}' "$output_maf_abs" > "$taxa_tmp"; then
    rm -f "$taxa_tmp"
    echo "Error: reference '$reference' was not found in the final MAF."
    exit 1
fi

mapfile -t other_taxa < <(sort "$taxa_tmp")
rm -f "$taxa_tmp"

taxa_args=("$reference")
for tax in "${other_taxa[@]}"; do
    taxa_args+=("$tax")
done

python "$script_dir/maf_to_concat_fasta.py" "${taxa_args[@]}" < "$output_maf_abs" > "$output_fasta_abs"

############################
# 14. Delete intermediate directories and intermediate MAF files,
#     keeping only the input MAF files and the final outputs
############################
for file in "${inputs[@]}"
do
    dir=$(basename "$file" .maf)
    rm -rf "$work_dir/$dir"
    rm -f "$work_dir/${dir}_noDup.maf"
    rm -f "$work_dir/${dir}_noDup_ordered.maf"
done

rm -rf "$work_dir/sep" "$work_dir/common_command_files" "$work_dir/common_output_files"

echo "Pipeline finished."
echo "Final MAF output: $output_maf_abs"
echo "Final FASTA output: $output_fasta_abs"