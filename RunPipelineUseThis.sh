#!/bin/bash
set -euo pipefail
shopt -s nullglob

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

patched_halWriteNucleotides="$script_dir/tool_used/bin/halWriteNucleotides"
patched_phyloFit="$script_dir/tool_used/bin/phyloFit"
patched_ancestorsML="$script_dir/tool_used/bin/ancestorsML"
patched_maf2hal="$script_dir/tool_used/bin/maf2hal"
patched_hal2fasta="$script_dir/tool_used/bin/hal2fasta"
patched_halAppendSubtree="$script_dir/tool_used/bin/halAppendSubtree"
patched_seqkit="$script_dir/tool_used/bin/seqkit"

for exe in \
    "$patched_halWriteNucleotides" \
    "$patched_phyloFit" \
    "$patched_ancestorsML" \
    "$patched_maf2hal" \
    "$patched_hal2fasta" \
    "$patched_halAppendSubtree" \
    "$patched_seqkit"
do
    if [[ ! -f "$exe" ]]; then
        echo "Error: required bundled executable not found: $exe"
        exit 1
    fi

    chmod +x "$exe"

    if [[ ! -x "$exe" ]]; then
        echo "Error: bundled executable is still not executable: $exe"
        exit 1
    fi
done

usage() {
    cat <<'USAGEEOF'
Usage:
  bash RunPipelineUseThis.sh --input /path/to/aln1.maf /path/to/aln2.maf ... --reference REF --pre output_prefix --threads 60 --common_workers 15 --genome_path genome.txt [options]

Required arguments:
  --input                 Input MAF files (space-separated).
  --reference             Reference taxon name. It will be placed first in the final FASTA output.
  --pre                   Output prefix. Two files will be generated first: <pre>.maf and <pre>.fasta
  --threads               Total CPU budget available on the machine.
  --common_workers        CPUs used by one consensus-extraction program.
                          The maximum number of chromosome-level ancestor-inference
                          jobs run in parallel is floor(threads / common_workers).
                          common_workers must be <= threads.
  --genome_path           Two-column file: TAXON GENOME_PATH. Extra taxa not present in the final alignment will be ignored.

Optional arguments:
  --Anc_name              Name of the consensus ancestor. Default: AncOfAllTaxa.
                          When --ifRefNonOutgroup 0, this is the ancestor of all
                          aligned taxa, with the reference treated as a direct
                          alignment element rather than an outgroup.
                          When --ifRefNonOutgroup 1, this is the top ancestor of
                          the reference and the ingroup ancestor.
                          Do not set it to the same value as --Anc_name_ofIngroup;
                          avoid using the reserved default name AncOfIngroup.
  --ModelFile             Path to the base model file used to generate <model>_useThis.mod
                          for phyloFit/ancestorsML.
                          Default: tryMLstartree.mod in the same directory as this script.
  --sep_length            Chunk length for long chromosome files. Default: 218505
  --chrom_length_threshold
                          Chromosome files longer than this threshold will be split. Default: 15819469
  --global_num            Number of internal subsegments used by one consensus-extraction program. Default: 40
  --separate_workers      CPUs used by the internal separate() stage of one consensus-extraction program. Default: 4
  --ifRefNonOutgroup      Whether the reference should be treated as a non-outgroup
                          direct alignment element for ancestor output.
                          0 = reference is not an outgroup; infer one ancestor
                              from reference + non-reference taxa. Default: 0.
                          1 = reference is an outgroup; infer an ingroup ancestor
                              first, then infer the top ancestor:
                              (reference,(non-reference taxa)Anc_name_ofIngroup)Anc_name.
  --ifDeleteImmdiFiles    Whether to keep intermediate files and directories.
                          0 = delete all intermediate files and directories. Default: 0.
                          1 = keep all intermediate files and directories.
  --Anc_name_ofIngroup    Ingroup ancestor name used only when --ifRefNonOutgroup 1,
                          i.e. when the reference is treated as an outgroup.
                          It is the ancestor of non-reference taxa and must differ
                          from --Anc_name.
                          Default: AncOfIngroup.
                          If it equals --Anc_name, both names will be automatically
                          reset to AncOfAllTaxa and AncOfAllTaxaExceptOutgroup.
  -h, --help              Show this help message and exit.

Notes:
  * This script may be run from any directory.
  * The code files are expected to be located in the same directory as this script:
      - seperate_maffile_by_refchrom.py
      - seperate_maffile_by_index.py
      - Extract_consensus_multi_alns.py
      - maf_to_concat_fasta.py
      - add_ance2mafwithN.py
      - check_if_maf_contain_all_taxachrom.py
      - replace_ancestor_seq_in_maf.py
      - sep_file/1.sh
  * add_ance2mafwithN.py is expected to support the fifth argument:
      python add_ance2mafwithN.py input.maf output.maf AncName SeqName RefGenome.fa
    so that per-chromosome ancestor rows are named AncName.SeqName and
    missing intervals are filled from the corresponding reference sequence.
  * All intermediate files and directories will be created under a tmp_XXXXXX directory
    inside the directory containing --pre.
  * By default, the temporary workspace is deleted after completion.
    Use --ifDeleteImmdiFiles 1 to retain all intermediate files and directories.
  * threads must be >= common_workers.
  * threads is the total CPU budget on the machine.
  * common_workers is the CPU usage of one consensus-extraction program.
  * The maximum number of chromosome-level ancestor-inference jobs run in parallel is also floor(threads / common_workers).
  * If threads=20, a practical choice for common_workers is often 10 or 15.
  * The maximum number of consensus jobs run in parallel is floor(threads / common_workers).
  * The long-chromosome splitting scripts are treated as ~1 CPU jobs, so up to threads of them may run in parallel.
USAGEEOF
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

run_parallel_commands() {
    local max_jobs="$1"
    shift
    local -a cmds=("$@")
    local fail=0
    local status_dir
    local idx=0

    if (( ${#cmds[@]} == 0 )); then
        return 0
    fi

    if ! is_positive_integer "$max_jobs"; then
        max_jobs=1
    fi

    status_dir=$(mktemp -d)

    for cmd in "${cmds[@]}"; do
        while (( $(jobs -rp | wc -l) >= max_jobs )); do
            sleep 1
        done

        idx=$((idx + 1))
        (
            set +e
            bash -euo pipefail -c "$cmd"
            rc=$?
            echo "$rc" > "$status_dir/$idx.rc"
            exit 0
        ) &
    done

    wait

    for rc_file in "$status_dir"/*.rc; do
        [[ -e "$rc_file" ]] || continue
        rc=$(cat "$rc_file")
        if [[ "$rc" != "0" ]]; then
            fail=1
        fi
    done

    rm -rf "$status_dir"

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
input_tags=()
threads=""
pre=""
output_prefix_abs=""
output_maf_abs=""
output_fasta_abs=""
output_hal_abs=""
reference=""
work_dir=""
sep_length=218505
chrom_length_threshold=15819469
global_num=40
separate_workers=4
common_workers=""
last_dir=""
genome_path_file=""
anc_name="AncOfAllTaxa"
anc_name_of_ingroup="AncOfIngroup"
ifRefNonOutgroup=0
ifDeleteImmdiFiles=0
model_file="$script_dir/tryMLstartree.mod"
output_dir=""
output_base=""
ancestor_prefix_abs=""
ancestor_maf_abs=""
ancestor_hal_abs=""
ancestor_tsv_abs=""
ancestor_fasta_abs=""
model_use_abs=""
jc69_root_abs=""
filtered_genome_txt_abs=""
info_tmp_dir=""
reference_genome_abs=""

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
        --reference)
            reference="$2"
            shift 2
            ;;
        --pre)
            pre="$2"
            shift 2
            ;;
        --genome_path)
            genome_path_file="$2"
            shift 2
            ;;
        --Anc_name)
            anc_name="$2"
            shift 2
            ;;
        --Anc_name_ofIngroup)
            anc_name_of_ingroup="$2"
            shift 2
            ;;
        --ifRefNonOutgroup)
            ifRefNonOutgroup="$2"
            shift 2
            ;;
        --ifDeleteImmdiFiles)
            ifDeleteImmdiFiles="$2"
            shift 2
            ;;
        --ModelFile)
            model_file="$2"
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

if [[ -z "$reference" ]]; then
    echo "Error: --reference is required."
    exit 1
fi

if [[ -z "$pre" ]]; then
    echo "Error: --pre is required."
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

if [[ -z "$genome_path_file" ]]; then
    echo "Error: --genome_path is required."
    exit 1
fi


for value_name in threads common_workers sep_length chrom_length_threshold global_num separate_workers; do
    value="${!value_name}"
    if ! is_positive_integer "$value"; then
        echo "Error: --${value_name} must be a positive integer."
        exit 1
    fi
done

if [[ ! "$ifRefNonOutgroup" =~ ^[01]$ ]]; then
    echo "Error: --ifRefNonOutgroup must be 0 or 1."
    exit 1
fi

if [[ ! "$ifDeleteImmdiFiles" =~ ^[01]$ ]]; then
    echo "Error: --ifDeleteImmdiFiles must be 0 or 1."
    exit 1
fi


if [[ "$ifRefNonOutgroup" == "1" && "$anc_name_of_ingroup" == "$anc_name" ]]; then
    echo "Warning: --Anc_name and --Anc_name_ofIngroup are identical when --ifRefNonOutgroup is 1." >&2
    echo "Warning: resetting --Anc_name to AncOfAllTaxa and --Anc_name_ofIngroup to AncOfAllTaxaExceptOutgroup." >&2
    anc_name="AncOfAllTaxa"
    anc_name_of_ingroup="AncOfAllTaxaExceptOutgroup"
fi

if (( threads < common_workers )); then
    echo "Error: --threads must be >= --common_workers."
    exit 1
fi

if [[ ! -f "$genome_path_file" ]]; then
    echo "Error: genome path file not found: $genome_path_file"
    exit 1
fi

if [[ ! -f "$model_file" ]]; then
    echo "Error: model file not found: $model_file"
    echo "Provide --ModelFile, or place tryMLstartree.mod in the same directory as this script: $script_dir"
    exit 1
fi

genome_path_file=$(readlink -f "$genome_path_file")
genome_path_dir=$(dirname "$genome_path_file")
model_file=$(readlink -f "$model_file")

############################
# 3. Resolve output prefix and create temporary working directory
############################
pre_dir_arg=$(dirname "$pre")
pre_base=$(basename "$pre")
mkdir -p "$pre_dir_arg"
output_dir=$(cd "$pre_dir_arg" && pwd)
output_prefix_abs="$output_dir/$pre_base"
output_maf_abs="${output_prefix_abs}.maf"
output_fasta_abs="${output_prefix_abs}.fasta"
output_hal_abs="${output_prefix_abs}.hal"
output_base=$(basename "$output_prefix_abs")
ancestor_prefix_abs="${output_prefix_abs}_${anc_name}"
ancestor_maf_abs="${output_prefix_abs}_${anc_name}.maf"
ancestor_hal_abs="${output_prefix_abs}_${anc_name}.hal"
ancestor_tsv_abs="${output_prefix_abs}_${anc_name}_JC69modle.tsv"
ancestor_fasta_abs="${output_dir}/${anc_name}.fa"
model_use_abs="${output_dir}/$(basename "${model_file%.*}")_useThis.mod"
jc69_root_abs="${output_dir}/JC69modle"

mkdir -p "$output_dir"
work_dir=$(mktemp -d "$output_dir/tmp_XXXXXX")
filtered_genome_txt_abs="$work_dir/${output_base}_${anc_name}.genome.txt"
info_tmp_dir="$work_dir/genome_info"
mkdir -p "$info_tmp_dir"

############################
# 4. Resolve input paths (inputs may be in different directories)
############################
idx=0
for file in "${raw_inputs[@]}"
do
    abs_file=$(readlink -f "$file")
    if [[ ! -f "$abs_file" ]]; then
        echo "Error: input file not found: $file"
        rm -rf "$work_dir"
        exit 1
    fi
    inputs+=("$abs_file")
    idx=$((idx + 1))
    tag=$(printf 'in%03d_%s' "$idx" "$(basename "$abs_file" .maf)")
    input_tags+=("$tag")
done

last_dir="${input_tags[${#input_tags[@]}-1]}"

############################
# 5. Check required files
############################
if [[ ! -f "$script_dir/seperate_maffile_by_refchrom.py" ]]; then
    echo "Error: seperate_maffile_by_refchrom.py not found in script directory."
    rm -rf "$work_dir"
    exit 1
fi

if [[ ! -f "$script_dir/seperate_maffile_by_index.py" ]]; then
    echo "Error: seperate_maffile_by_index.py not found in script directory."
    rm -rf "$work_dir"
    exit 1
fi

if [[ ! -f "$script_dir/Extract_consensus_multi_alns.py" ]]; then
    echo "Error: Extract_consensus_multi_alns.py not found in script directory."
    rm -rf "$work_dir"
    exit 1
fi

if [[ ! -f "$script_dir/maf_to_concat_fasta.py" ]]; then
    echo "Error: maf_to_concat_fasta.py not found in script directory."
    rm -rf "$work_dir"
    exit 1
fi

if [[ ! -f "$script_dir/add_ance2mafwithN.py" ]]; then
    echo "Error: add_ance2mafwithN.py not found in script directory."
    rm -rf "$work_dir"
    exit 1
fi

if [[ ! -f "$script_dir/check_if_maf_contain_all_taxachrom.py" ]]; then
    echo "Error: check_if_maf_contain_all_taxachrom.py not found in script directory."
    rm -rf "$work_dir"
    exit 1
fi

if [[ ! -f "$script_dir/replace_ancestor_seq_in_maf.py" ]]; then
    echo "Error: replace_ancestor_seq_in_maf.py not found in script directory."
    rm -rf "$work_dir"
    exit 1
fi

############################
# 6. Split each input MAF by reference sequence ID in parallel
############################
split_input_cmds=()
for i in "${!inputs[@]}"
do
    file="${inputs[$i]}"
    dir="${input_tags[$i]}"
    split_input_cmds+=("mkdir -p \"$work_dir/$dir\" && python \"$script_dir/seperate_maffile_by_refchrom.py\" \"$file\" \"$work_dir/$dir\"")
done

run_parallel_commands "$threads" "${split_input_cmds[@]}"

############################
# 7. Prepare intermediate directories in the working directory
############################
mkdir -p "$work_dir/common_command_files"
mkdir -p "$work_dir/common_output_files"

# Record chromosome MAFs that need coordinate splitting.
# Keep an explicit count because Bash 4.2 + `set -u` can treat expansion of
# an empty array (e.g. "${long_chroms[@]}") as an unbound-variable error.
long_chroms=()
long_chrom_count=0

############################
# 8. Use the directory generated from the last input file as the template.
#    For long chromosomes, prepare consensus commands for all expected chunks.
#    The actual chunk files will be generated later by
#    seperate_maffile_by_index.py, which scans each chromosome MAF
#    only once and writes all chunks in one run.
############################

for i in "$work_dir/${last_dir}"/*.maf
do
    [ -e "$i" ] || continue

    head_chrom=$(basename "$i" .maf)

    # srcSize of the first reference s-line.
    # This is the total reference chromosome/sequence length.
    length=$(awk '$1=="s"{print $6; exit}' "$i")

    mkdir -p "$work_dir/common_command_files/${head_chrom}"
    mkdir -p "$work_dir/common_output_files/${head_chrom}"

    if (( length > chrom_length_threshold )); then

        # Remember this chromosome so that every input alignment can be
        # split once in Step 10.
        long_chroms[$long_chrom_count]="$head_chrom"
        long_chrom_count=$((long_chrom_count + 1))

        num=$(( (length + sep_length - 1) / sep_length ))

        echo "Long chromosome detected: ${head_chrom}"
        echo "  Reference length: ${length}"
        echo "  Chunk length: ${sep_length}"
        echo "  Number of chunks: ${num}"

        # Prepare the downstream consensus command corresponding to
        # every chunk that the faster splitter will generate.
        for n in $(seq 1 "$num")
        do
            start=$(( (n - 1) * sep_length ))
            end=$(( n * sep_length ))

            if (( end > length )); then
                end=$length
            fi

            input_args=""

            for j in "${!inputs[@]}"
            do
                d="${input_tags[$j]}"

                input_args="${input_args} \"$work_dir/${d}/${head_chrom}/${n}_${start}_${end}.maf\""
            done

            echo "python \"$script_dir/Extract_consensus_multi_alns.py\" --input${input_args} --output \"$work_dir/common_output_files/${head_chrom}/part${n}.maf\" --global_num ${global_num} --global_start ${start} --global_end ${end} --separate_workers ${separate_workers} --common_workers ${common_workers}" \
                > "$work_dir/common_command_files/${head_chrom}/${n}.sh"
        done

    else

        # Short chromosomes are not split.
        input_args=""

        for j in "${!inputs[@]}"
        do
            d="${input_tags[$j]}"
            input_args="${input_args} \"$work_dir/${d}/${head_chrom}.maf\""
        done

        echo "python \"$script_dir/Extract_consensus_multi_alns.py\" --input${input_args} --output \"$work_dir/common_output_files/${head_chrom}/part1.maf\" --global_num ${global_num} --global_start 0 --global_end ${length} --separate_workers ${separate_workers} --common_workers ${common_workers}" \
            > "$work_dir/common_command_files/${head_chrom}/1.sh"
    fi
done


############################
# 9. No separate script-copying stage is needed anymore.
#
# The old pipeline generated hundreds of small shell scripts under sep_<chrom>
# and copied them into every input directory.
#
# The faster splitter instead handles one whole chromosome MAF per invocation.
############################


############################
# 10. Split each long chromosome MAF in parallel.
#
# IMPORTANT:
# One Python process now handles ONE chromosome MAF:
#
#     input chromosome MAF
#         -> scan once
#         -> generate all coordinate chunks
#
# Example:
#
#     GALGA.1.maf
#         -> 1_0_218505.maf
#         -> 2_218505_437010.maf
#         -> ...
#         -> 920_200806095_200994015.maf
#
# Up to --threads chromosome-level splitting jobs may run simultaneously.
############################

split_job_cmds=()

# If every reference sequence is shorter than --chrom_length_threshold,
# long_chrom_count remains 0. In that case Step 10 must be skipped entirely.
# This avoids expanding an empty long_chroms array under Bash 4.2 + `set -u`.
if (( long_chrom_count > 0 )); then
    for j in "${!inputs[@]}"
    do
        dir="${input_tags[$j]}"

        for (( k=0; k<long_chrom_count; k++ ))
        do
            head_chrom="${long_chroms[$k]}"
            input_maf="$work_dir/$dir/${head_chrom}.maf"
            output_chunk_dir="$work_dir/$dir/${head_chrom}"

            if [[ ! -f "$input_maf" ]]; then
                echo "Error: chromosome MAF not found: $input_maf"
                exit 1
            fi

            split_job_cmds+=(
                "mkdir -p \"$output_chunk_dir\" && python \"$script_dir/seperate_maffile_by_index.py\" \"$input_maf\" \"$sep_length\" \"$output_chunk_dir\""
            )
        done
    done
fi

if (( ${#split_job_cmds[@]} > 0 )); then
    echo "Starting faster chromosome splitting."
    echo "Number of chromosome-level splitting jobs: ${#split_job_cmds[@]}"
    echo "Maximum simultaneous splitting jobs: ${threads}"

    run_parallel_commands "$threads" "${split_job_cmds[@]}"
fi

############################
# 11. Execute the commands in common_command_files in parallel
#     Maximum parallel consensus jobs = floor(threads / common_workers)
############################
max_common_jobs=$(( threads / common_workers ))
if (( max_common_jobs < 1 )); then
    max_common_jobs=1
fi

# Use the same CPU-budget logic for chromosome-level ancestor inference.
# Example: threads=80 and common_workers=5 -> at most 16 chromosome jobs
# run simultaneously in each ancestor-inference round.
max_ancestor_jobs=$(( threads / common_workers ))
if (( max_ancestor_jobs < 1 )); then
    max_ancestor_jobs=1
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

if (( ${#common_job_cmds[@]} > 0 )); then
    run_parallel_commands "$max_common_jobs" "${common_job_cmds[@]}"
fi

############################
# 12. First merge part*.maf for each chromosome into <chrom>.maf
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
# 13. Then merge all chromosome MAFs into the final MAF output
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
# 14. Convert the final MAF to concatenated FASTA
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
    rm -rf "$work_dir"
    exit 1
fi

mapfile -t other_taxa < <(LC_ALL=C sort -f -s "$taxa_tmp")
rm -f "$taxa_tmp"

taxa_args=("$reference")
for tax in "${other_taxa[@]}"; do
    taxa_args+=("$tax")
done

python "$script_dir/maf_to_concat_fasta.py" "${taxa_args[@]}" < "$output_maf_abs" > "$output_fasta_abs"

############################
# 15. Filter genome_path to taxa present in the final alignment and build .info.txt files
############################
declare -A wanted_taxa=()
declare -A genome_paths=()
for tax in "${taxa_args[@]}"; do
    wanted_taxa["$tax"]=1
done

: > "$filtered_genome_txt_abs"
while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "$line" ]] && continue
    [[ "$line" =~ ^[[:space:]]*# ]] && continue
    read -r tax path _ <<< "$line"
    [[ -z "${tax:-}" || -z "${path:-}" ]] && continue
    if [[ -n "${wanted_taxa[$tax]+x}" ]]; then
        if [[ "$path" = /* ]]; then
            abs_path=$(readlink -f "$path")
        else
            abs_path=$(readlink -f "$genome_path_dir/$path")
        fi
        genome_paths["$tax"]="$abs_path"
        printf '%s %s\n' "$tax" "$abs_path" >> "$filtered_genome_txt_abs"
    fi
done < "$genome_path_file"

for tax in "${taxa_args[@]}"; do
    if [[ -z "${genome_paths[$tax]+x}" ]]; then
        echo "Error: taxon '$tax' is present in the final MAF but missing from --genome_path file."
        rm -rf "$work_dir"
        exit 1
    fi
done

reference_genome_abs="${genome_paths[$reference]}"
if [[ ! -f "$reference_genome_abs" ]]; then
    echo "Error: reference genome FASTA not found: $reference_genome_abs"
    rm -rf "$work_dir"
    exit 1
fi

for tax in "${taxa_args[@]}"; do
    genome_abs="${genome_paths[$tax]}"
    info_abs="$info_tmp_dir/${tax}.info.txt"
    "$patched_seqkit" fx2tab -n -i -l "$genome_abs" > "$info_abs"
done

############################
# 16. Infer ancestor HAL/FASTA branches
#     In this version, ancestor inference is done per reference chromosome.
#     Each chromosome consensus MAF first gets chromosome-specific ancestor rows
#     such as AncName.chr1. After per-chromosome ancestor sequences are inferred
#     and written back into MAFs, the chromosome-level MAFs/HALs/FASTAs are merged
#     into the final outputs.
############################
if (( ${#other_taxa[@]} == 0 )); then
    echo "Error: no non-reference taxa found in the final MAF; cannot build ancestor model."
    exit 1
fi

join_by_comma() {
    local IFS=,
    echo "$*"
}

get_reference_lookup_seq_name() {
    local seq_id="$1"
    local label="$2"
    local seq_part="$seq_id"
    local prefix="${label}."

    if [[ -n "$label" && "$seq_part" == "$prefix"* ]]; then
        seq_part="${seq_part#${prefix}}"
    fi

    echo "$seq_part"
}

clean_reference_seq_name() {
    local seq_id="$1"
    local label="$2"

    python - "$seq_id" "$label" <<'PYCLEAN'
import sys

seq_id = sys.argv[1]
label = sys.argv[2]


def is_connector(ch: str) -> bool:
    return not ch.isalnum()


def clean_id(seq_id: str, label: str) -> str:
    if not label or label not in seq_id:
        return seq_id

    idx = seq_id.find(label)
    before = seq_id[:idx]
    after = seq_id[idx + len(label):]

    # Case 1:
    #   GALGAchr1   -> chr1
    #   GALGA.chr1  -> chr1
    #   GALGA_chr1  -> chr1
    #   GALGA/chr1  -> chr1
    if before == "":
        if after[:1] and is_connector(after[0]):
            after = after[1:]
        return after if after else seq_id

    # Case 2:
    #   a_GALGA_chr1 -> a_chr1
    #   a/GALGA/chr1 -> a/chr1
    #   a.GALGA.chr1 -> a.chr1
    #
    # If both sides of label are connector-like, keep the left connector
    # and remove the right connector, so we do not produce a__chr1 or a//chr1.
    if before[-1:] and after[:1] and is_connector(before[-1]) and is_connector(after[0]):
        after = after[1:]

    new_id = before + after
    return new_id if new_id else seq_id

# If the MAF-derived chromosome name starts with the standard MAF prefix
# TAXON.SEQ, strip only that leading TAXON. first, then clean the remaining
# sequence ID. This handles GALGA.GALGA.chr20 -> GALGA.chr20 -> chr20.
seq_part = seq_id
prefix = label + "."
if label and seq_part.startswith(prefix):
    seq_part = seq_part[len(prefix):]

print(clean_id(seq_part, label))
PYCLEAN
}

build_ladder_tree_from_taxa() {
    local taxa_space="$*"
    awk -v t="$taxa_space" 'BEGIN{
        n=split(t,a," ")
        if(n < 1){exit 2}
        s=a[1]":1"
        for(i=2;i<=n;i++){
            edge=(i<n?":0!":":0")
            s="(" s "," a[i]":1)" edge
        }
        print s
    }'
}

remove_reference_rows_from_maf() {
    local maf_abs="$1"
    local ref="$2"
    local tmp_maf="${maf_abs%.maf}_2.maf"

    LC_ALL=C awk -v ref="$ref" '
    $1 ~ /^(s|i|e|q)$/ {
        split($2,a,".")
        if (a[1] == ref) next
    }
    {print}
    ' "$maf_abs" > "$tmp_maf"

    mv "$tmp_maf" "$maf_abs"
}

remove_taxa_rows_from_maf() {
    local input_maf="$1"
    local output_maf="$2"
    shift 2
    local taxa_to_remove="$*"

    LC_ALL=C awk -v taxa="$taxa_to_remove" '
    BEGIN{
        n=split(taxa,a," ")
        for(i=1;i<=n;i++){
            if(a[i] != "") remove[a[i]]=1
        }
    }
    $1 ~ /^(s|i|e|q)$/ {
        split($2,b,".")
        if (remove[b[1]]) next
    }
    {print}
    ' "$input_maf" > "$output_maf"
}

replace_reference_sequence_id_in_maf() {
    local input_maf="$1"
    local output_maf="$2"
    local ref="$3"
    local old_seq="$4"
    local new_seq="$5"

    LC_ALL=C awk -v ref="$ref" -v old_seq="$old_seq" -v new_seq="$new_seq" '
    BEGIN{
        old_id=ref "." old_seq
        new_id=ref "." new_seq
    }
    $1 ~ /^(s|i|e|q)$/ {
        if ($2 == old_id) $2 = new_id
    }
    {print}
    ' "$input_maf" > "$output_maf"
}


merge_maf_files() {
    local output_maf="$1"
    shift
    local -a maf_files=("$@")

    : > "$output_maf"

    if (( ${#maf_files[@]} == 0 )); then
        echo "Error: no MAF files provided for merge: $output_maf"
        exit 1
    fi

    local first=1
    local f=""

    while IFS= read -r f
    do
        [[ -n "$f" ]] || continue
        if (( first == 1 )); then
            cat "$f" >> "$output_maf"
            first=0
        else
            awk 'BEGIN{p=0} !p && /^#/ {next} {p=1; print}' "$f" >> "$output_maf"
        fi
    done < <(printf '%s\n' "${maf_files[@]}" | sort -V)
}

merge_fasta_files() {
    local output_fasta="$1"
    shift
    local -a fasta_files=("$@")
    local f=""

    : > "$output_fasta"

    if (( ${#fasta_files[@]} == 0 )); then
        echo "Error: no FASTA files provided for merge: $output_fasta"
        exit 1
    fi

    while IFS= read -r f
    do
        [[ -n "$f" ]] || continue
        cat "$f" >> "$output_fasta"
    done < <(printf '%s\n' "${fasta_files[@]}" | sort -V)
}

write_reference_sequence_as_fasta() {
    local reference_fasta="$1"
    local lookup_seq_name="$2"
    local output_seq_name="$3"
    local output_fasta="$4"

    python - "$reference_fasta" "$lookup_seq_name" "$output_seq_name" "$output_fasta" <<'PYREF'
import sys

reference_fasta = sys.argv[1]
lookup_seq_name = sys.argv[2]
output_seq_name = sys.argv[3]
output_fasta = sys.argv[4]

found = False
collect = False
seq_parts = []

with open(reference_fasta, "r", encoding="utf-8") as fin:
    for line in fin:
        if line.startswith(">"):
            current_id = line[1:].strip().split()[0]
            if collect:
                break
            collect = (current_id == lookup_seq_name)
            if collect:
                found = True
        else:
            if collect:
                seq_parts.append(line.strip())

if not found:
    sys.stderr.write(
        f"ERROR: sequence '{lookup_seq_name}' was not found in reference FASTA: {reference_fasta}\n"
    )
    sys.exit(1)

seq = "".join(seq_parts)

with open(output_fasta, "w", encoding="utf-8") as fout:
    fout.write(f">{output_seq_name}\n")
    for i in range(0, len(seq), 60):
        fout.write(seq[i:i+60] + "\n")
PYREF
}

run_task_level_phylofit_model() {
    local branch_model_tree="$1"
    local branch_maf_abs="$2"
    local branch_model_use_abs="$3"
    local branch_jc69_root_abs="$4"

    cp "$model_file" "$branch_model_use_abs"
    printf 'TREE: %s;\n' "$branch_model_tree" >> "$branch_model_use_abs"

    "$patched_phyloFit" --init-model "$branch_model_use_abs" --subst-mod JC69 --out-root "$branch_jc69_root_abs" --msa-format MAF "$branch_maf_abs"

    if [[ ! -s "${branch_jc69_root_abs}.mod" ]]; then
        echo "Error: task-level phyloFit did not produce model file: ${branch_jc69_root_abs}.mod" >&2
        exit 1
    fi
}

run_ancestors_and_hal2fasta_with_model() {
    local branch_hal_abs="$1"
    local branch_ancestor_name="$2"
    local branch_model_mod_abs="$3"
    local branch_tsv_abs="$4"
    local branch_fasta_abs="$5"
    local branch_log_prefix_abs="$6"

    if [[ ! -s "$branch_model_mod_abs" ]]; then
        echo "Error: ancestor model file not found or empty: $branch_model_mod_abs" >&2
        exit 1
    fi

    "$patched_ancestorsML" --printWrites "$branch_hal_abs" "$branch_ancestor_name" "$branch_model_mod_abs" > "$branch_tsv_abs"
    env -u LD_LIBRARY_PATH "$patched_halWriteNucleotides" "$branch_hal_abs" "$branch_tsv_abs" > "${branch_log_prefix_abs}_halWriteNucleotides.log"
    "$patched_hal2fasta" "$branch_hal_abs" "$branch_ancestor_name" > "$branch_fasta_abs"
}

run_per_chrom_direct_ancestor_branch() {
    local branch_ancestor_name="$anc_name"
    local branch_label="reference is not an outgroup; one ancestor inferred from reference + non-reference taxa"
    local ingroup_tree=""
    local model_tree=""
    local chrom_maf_abs=""
    local chrom_name=""
    local ref_seq_name=""
    local ref_seq_lookup_name=""
    local chrom_dir=""
    local chrom_prefix_abs=""
    local chrom_added_maf_abs=""
    local chrom_reference_fasta_abs=""
    local chrom_maf_for_add_abs=""
    local chrom_hal_abs=""
    local chrom_jc69_root_abs=""
    local chrom_tsv_abs=""
    local chrom_fasta_abs=""
    local chrom_new_maf_abs=""
    local tax=""
    local info_abs=""
    local branch_task_model_input_maf_abs=""
    local branch_task_model_use_abs=""
    local branch_task_jc69_root_abs=""

    local -a chrom_added_maf_files=()
    local -a chrom_new_maf_files=()
    local -a chrom_fasta_files=()

    echo "Running ancestor branch: ${branch_label}"
    echo "  Ancestor name: ${branch_ancestor_name}"
    echo "  Ancestor inference unit: per reference chromosome"

    ingroup_tree=$(build_ladder_tree_from_taxa "${other_taxa[@]}")
    model_tree="(${reference}:1,${ingroup_tree})"

    local -a chrom_maf_files=( "$work_dir"/common_output_files/*.maf )
    if (( ${#chrom_maf_files[@]} == 0 )); then
        echo "Error: no chromosome consensus MAF files found under $work_dir/common_output_files"
        exit 1
    fi

    while IFS= read -r chrom_maf_abs
    do
        [[ -n "$chrom_maf_abs" ]] || continue

        chrom_name=$(basename "$chrom_maf_abs" .maf)
        ref_seq_lookup_name=$(get_reference_lookup_seq_name "$chrom_name" "$reference")
        ref_seq_name=$(clean_reference_seq_name "$chrom_name" "$reference")
        chrom_dir="$work_dir/common_output_files/${chrom_name}"
        mkdir -p "$chrom_dir"

        chrom_prefix_abs="$chrom_dir/${output_base}_${branch_ancestor_name}"
        chrom_added_maf_abs="${chrom_prefix_abs}.maf"

        echo "  Preparing chromosome for task-level model: ${chrom_name}"

        if ! awk '$1=="s"{found=1; exit} END{exit found ? 0 : 1}' "$chrom_maf_abs"; then
            chrom_fasta_abs="${chrom_dir}/${branch_ancestor_name}.${ref_seq_name}.fa"

            echo "  No consensus alignment blocks found for ${chrom_name}; copying reference sequence to ${branch_ancestor_name}.${ref_seq_name}.fa."

            write_reference_sequence_as_fasta \
                "$reference_genome_abs" \
                "$ref_seq_lookup_name" \
                "$ref_seq_name" \
                "$chrom_fasta_abs"

            chrom_fasta_files+=( "$chrom_fasta_abs" )
            continue
        fi

        chrom_maf_for_add_abs="${chrom_prefix_abs}_refseq_for_add.maf"
        replace_reference_sequence_id_in_maf \
            "$chrom_maf_abs" \
            "$chrom_maf_for_add_abs" \
            "$reference" \
            "$ref_seq_lookup_name" \
            "$ref_seq_name"

        chrom_reference_fasta_abs="${chrom_dir}/${reference}.${ref_seq_name}.reference_for_${branch_ancestor_name}.fa"
        write_reference_sequence_as_fasta \
            "$reference_genome_abs" \
            "$ref_seq_lookup_name" \
            "$ref_seq_name" \
            "$chrom_reference_fasta_abs"

        python "$script_dir/add_ance2mafwithN.py" \
            "$chrom_maf_for_add_abs" \
            "$chrom_added_maf_abs" \
            "$branch_ancestor_name" \
            "$ref_seq_name" \
            "$chrom_reference_fasta_abs"

        replace_reference_sequence_id_in_maf \
            "$chrom_added_maf_abs" \
            "${chrom_added_maf_abs}.tmp_restore" \
            "$reference" \
            "$ref_seq_name" \
            "$ref_seq_lookup_name"
        mv "${chrom_added_maf_abs}.tmp_restore" "$chrom_added_maf_abs"

        chrom_added_maf_files+=( "$chrom_added_maf_abs" )

    done < <(printf '%s\n' "${chrom_maf_files[@]}" | sort -V)

    branch_task_model_input_maf_abs="${output_prefix_abs}_${branch_ancestor_name}.maf"
    merge_maf_files "$branch_task_model_input_maf_abs" "${chrom_added_maf_files[@]}"

    branch_task_model_use_abs="${output_prefix_abs}_${branch_ancestor_name}_$(basename "${model_file%.*}")_useThis.mod"
    branch_task_jc69_root_abs="${output_prefix_abs}_${branch_ancestor_name}_JC69modle"
    run_task_level_phylofit_model \
        "$model_tree" \
        "$branch_task_model_input_maf_abs" \
        "$branch_task_model_use_abs" \
        "$branch_task_jc69_root_abs"

    ############################
    # Parallel per-chromosome ancestor inference.
    #
    # All chromosomes use the same task-level model generated above.
    # Each chromosome writes only to its own chromosome directory, so these
    # jobs are independent and can safely run in parallel.
    ############################

    local direct_job_dir="$work_dir/ancestor_jobs/direct_${branch_ancestor_name}"
    local -a direct_job_cmds=()
    local job_script=""
    local job_idx=0

    mkdir -p "$direct_job_dir"

    while IFS= read -r chrom_added_maf_abs
    do
        [[ -n "$chrom_added_maf_abs" ]] || continue

        chrom_dir=$(dirname "$chrom_added_maf_abs")
        chrom_name=$(basename "$chrom_dir")
        ref_seq_lookup_name=$(get_reference_lookup_seq_name "$chrom_name" "$reference")
        ref_seq_name=$(clean_reference_seq_name "$chrom_name" "$reference")

        chrom_prefix_abs="$chrom_dir/${output_base}_${branch_ancestor_name}"
        chrom_hal_abs="${chrom_prefix_abs}.hal"
        chrom_jc69_root_abs="${chrom_prefix_abs}_JC69modle"
        chrom_tsv_abs="${chrom_prefix_abs}_JC69modle.tsv"
        chrom_fasta_abs="${chrom_dir}/${branch_ancestor_name}.${ref_seq_name}.fa"
        chrom_new_maf_abs="${chrom_prefix_abs}_new.maf"

        # Record expected outputs in the parent shell before launching child jobs.
        # Changes to arrays inside child processes would not propagate back.
        chrom_new_maf_files+=( "$chrom_new_maf_abs" )
        chrom_fasta_files+=( "$chrom_fasta_abs" )

        job_idx=$((job_idx + 1))
        job_script=$(printf '%s/job_%04d.sh' "$direct_job_dir" "$job_idx")

        cat > "$job_script" <<EOF
#!/bin/bash
set -euo pipefail

echo "  [ancestor job ${job_idx}] Processing chromosome: ${chrom_name}"

if [[ ! -s "${branch_task_jc69_root_abs}.mod" ]]; then
    echo "Error: ancestor model file not found or empty: ${branch_task_jc69_root_abs}.mod" >&2
    exit 1
fi

"$patched_maf2hal" \
    --refGenome "$branch_ancestor_name" \
    "$chrom_added_maf_abs" \
    "$chrom_hal_abs"

"$patched_ancestorsML" \
    --printWrites \
    "$chrom_hal_abs" \
    "$branch_ancestor_name" \
    "${branch_task_jc69_root_abs}.mod" \
    > "$chrom_tsv_abs"

env -u LD_LIBRARY_PATH \
    "$patched_halWriteNucleotides" \
    "$chrom_hal_abs" \
    "$chrom_tsv_abs" \
    > "${chrom_jc69_root_abs}_halWriteNucleotides.log"

"$patched_hal2fasta" \
    "$chrom_hal_abs" \
    "$branch_ancestor_name" \
    > "$chrom_fasta_abs"

python "$script_dir/replace_ancestor_seq_in_maf.py" \
    "$chrom_added_maf_abs" \
    "$branch_ancestor_name" \
    "$chrom_fasta_abs" \
    "$chrom_new_maf_abs" \
    "$ref_seq_name"

echo "  [ancestor job ${job_idx}] Finished chromosome: ${chrom_name}"
EOF

        chmod +x "$job_script"
        direct_job_cmds+=( "bash \"$job_script\"" )

    done < <(printf '%s\n' "${chrom_added_maf_files[@]}" | sort -V)

    echo "Starting parallel chromosome-level ancestor inference."
    echo "  Number of chromosome jobs: ${#direct_job_cmds[@]}"
    echo "  Maximum simultaneous ancestor jobs: ${max_ancestor_jobs}"

    if (( ${#direct_job_cmds[@]} > 0 )); then
        run_parallel_commands "$max_ancestor_jobs" "${direct_job_cmds[@]}"
    fi

    echo "All chromosome-level ancestor inference jobs finished."

    local final_ancestor_maf_abs="${output_prefix_abs}_${branch_ancestor_name}_new.maf"
    local final_ancestor_checked_maf_abs="${output_prefix_abs}_${branch_ancestor_name}_new2.maf"
    local final_ancestor_fasta_abs="${output_dir}/${branch_ancestor_name}.fa"

    merge_maf_files "$final_ancestor_maf_abs" "${chrom_new_maf_files[@]}"

    merge_fasta_files "$final_ancestor_fasta_abs" "${chrom_fasta_files[@]}"
    "$patched_seqkit" fx2tab -n -i -l "$final_ancestor_fasta_abs" > "$info_tmp_dir/${branch_ancestor_name}.info.txt"

    cp "$final_ancestor_maf_abs" "$final_ancestor_checked_maf_abs"
    for tax in "${taxa_args[@]}" "$branch_ancestor_name"; do
        info_abs="$info_tmp_dir/${tax}.info.txt"
        python "$script_dir/check_if_maf_contain_all_taxachrom.py" "$info_abs" "$final_ancestor_checked_maf_abs" "$tax"
    done

    "$patched_maf2hal" \
        --refGenome "$branch_ancestor_name" \
        "$final_ancestor_checked_maf_abs" \
        "$output_hal_abs"

    generated_branch_summaries+=(
        "${branch_label}|${branch_ancestor_name}|${output_hal_abs}|NA|NA|NA|${final_ancestor_fasta_abs}"
    )
}

run_per_chrom_ref_outgroup_branch() {
    local ingroup_ancestor_name="$anc_name_of_ingroup"
    local top_ancestor_name="$anc_name"

    local target_genomes_csv=""
    local ingroup_model_tree=""
    local top_model_tree=""

    local chrom_maf_abs=""
    local chrom_name=""
    local ref_seq_name=""
    local ref_seq_lookup_name=""
    local chrom_dir=""
    local chrom_prefix_abs=""
    local chrom_reference_fasta_abs=""
    local chrom_maf_for_add_abs=""

    local ingroup_added_maf_abs=""
    local ingroup_noref_maf_abs=""
    local ingroup_hal_abs=""
    local ingroup_jc69_root_abs=""
    local ingroup_tsv_abs=""
    local ingroup_fasta_abs=""
    local ingroup_with_sequence_maf_abs=""
    local ingroup_noref_new_maf_abs=""

    local ref_plus_ingroup_maf_abs=""
    local top_added_maf_abs=""
    local top_hal_abs=""
    local top_jc69_root_abs=""
    local top_tsv_abs=""
    local top_fasta_abs=""
    local top_noingroup_new_maf_abs=""

    local ingroup_task_noref_maf_abs=""
    local ingroup_task_model_use_abs=""
    local ingroup_task_jc69_root_abs=""
    local top_task_added_maf_abs=""
    local top_task_model_use_abs=""
    local top_task_jc69_root_abs=""

    local tax=""
    local info_abs=""

    local -a ingroup_noref_maf_files=()
    local -a ingroup_noref_new_maf_files=()
    local -a ref_plus_ingroup_maf_files=()
    local -a top_added_maf_files=()
    local -a top_noingroup_new_maf_files=()
    local -a ingroup_fasta_files=()
    local -a top_fasta_files=()

    echo "Running reference-outgroup ancestor branch"
    echo "  Ingroup ancestor name: ${ingroup_ancestor_name}"
    echo "  Top ancestor name: ${top_ancestor_name}"
    echo "  Ancestor inference unit: per reference chromosome"

    target_genomes_csv=$(join_by_comma "${other_taxa[@]}")
    ingroup_model_tree=$(build_ladder_tree_from_taxa "${other_taxa[@]}")
    top_model_tree="(${reference}:1,${ingroup_ancestor_name}:1)"

    local -a chrom_maf_files=( "$work_dir"/common_output_files/*.maf )
    if (( ${#chrom_maf_files[@]} == 0 )); then
        echo "Error: no chromosome consensus MAF files found under $work_dir/common_output_files"
        exit 1
    fi

    while IFS= read -r chrom_maf_abs
    do
        [[ -n "$chrom_maf_abs" ]] || continue

        chrom_name=$(basename "$chrom_maf_abs" .maf)
        ref_seq_lookup_name=$(get_reference_lookup_seq_name "$chrom_name" "$reference")
        ref_seq_name=$(clean_reference_seq_name "$chrom_name" "$reference")
        chrom_dir="$work_dir/common_output_files/${chrom_name}"
        mkdir -p "$chrom_dir"

        chrom_prefix_abs="$chrom_dir/${output_base}_${ingroup_ancestor_name}"

        ingroup_added_maf_abs="${chrom_prefix_abs}.maf"
        ingroup_noref_maf_abs="${chrom_prefix_abs}_noref.maf"

        echo "  Preparing chromosome for ingroup task-level model: ${chrom_name}"

        if ! awk '$1=="s"{found=1; exit} END{exit found ? 0 : 1}' "$chrom_maf_abs"; then
            ingroup_fasta_abs="${chrom_dir}/${ingroup_ancestor_name}.${ref_seq_name}.fa"
            top_fasta_abs="${chrom_dir}/${top_ancestor_name}.${ref_seq_name}.fa"

            echo "  No consensus alignment blocks found for ${chrom_name}; copying reference sequence to ${ingroup_ancestor_name}.${ref_seq_name}.fa and ${top_ancestor_name}.${ref_seq_name}.fa."

            write_reference_sequence_as_fasta \
                "$reference_genome_abs" \
                "$ref_seq_lookup_name" \
                "$ref_seq_name" \
                "$ingroup_fasta_abs"

            write_reference_sequence_as_fasta \
                "$reference_genome_abs" \
                "$ref_seq_lookup_name" \
                "$ref_seq_name" \
                "$top_fasta_abs"

            ingroup_fasta_files+=( "$ingroup_fasta_abs" )
            top_fasta_files+=( "$top_fasta_abs" )

            continue
        fi

        # Step 1. Add the ingroup ancestor row as IngroupAncestor.chromName.
        chrom_maf_for_add_abs="${chrom_prefix_abs}_refseq_for_add.maf"
        replace_reference_sequence_id_in_maf \
            "$chrom_maf_abs" \
            "$chrom_maf_for_add_abs" \
            "$reference" \
            "$ref_seq_lookup_name" \
            "$ref_seq_name"

        chrom_reference_fasta_abs="${chrom_dir}/${reference}.${ref_seq_name}.reference_for_${ingroup_ancestor_name}.fa"
        write_reference_sequence_as_fasta \
            "$reference_genome_abs" \
            "$ref_seq_lookup_name" \
            "$ref_seq_name" \
            "$chrom_reference_fasta_abs"

        python "$script_dir/add_ance2mafwithN.py" \
            "$chrom_maf_for_add_abs" \
            "$ingroup_added_maf_abs" \
            "$ingroup_ancestor_name" \
            "$ref_seq_name" \
            "$chrom_reference_fasta_abs"

        replace_reference_sequence_id_in_maf \
            "$ingroup_added_maf_abs" \
            "${ingroup_added_maf_abs}.tmp_restore" \
            "$reference" \
            "$ref_seq_name" \
            "$ref_seq_lookup_name"
        mv "${ingroup_added_maf_abs}.tmp_restore" "$ingroup_added_maf_abs"

        # Step 2. Create the reference-excluded copy used both for
        # task-level ingroup model fitting and per-chromosome ingroup inference.
        cp "$ingroup_added_maf_abs" "$ingroup_noref_maf_abs"
        remove_reference_rows_from_maf "$ingroup_noref_maf_abs" "$reference"

        ingroup_noref_maf_files+=( "$ingroup_noref_maf_abs" )

    done < <(printf '%s\n' "${chrom_maf_files[@]}" | sort -V)

    ingroup_task_noref_maf_abs="${output_prefix_abs}_${ingroup_ancestor_name}_noref.maf"
    merge_maf_files "$ingroup_task_noref_maf_abs" "${ingroup_noref_maf_files[@]}"

    ingroup_task_model_use_abs="${output_prefix_abs}_${ingroup_ancestor_name}_$(basename "${model_file%.*}")_useThis.mod"
    ingroup_task_jc69_root_abs="${output_prefix_abs}_${ingroup_ancestor_name}_JC69modle"
    run_task_level_phylofit_model \
        "$ingroup_model_tree" \
        "$ingroup_task_noref_maf_abs" \
        "$ingroup_task_model_use_abs" \
        "$ingroup_task_jc69_root_abs"

    ############################
    # Round 1:
    # Parallel per-chromosome inference of the ingroup ancestor.
    #
    # Barrier rule:
    # all ingroup chromosome jobs must finish before the top-level MAF is
    # merged and the top-ancestor model is fitted.
    ############################

    local ingroup_job_dir="$work_dir/ancestor_jobs/ingroup_${ingroup_ancestor_name}"
    local -a ingroup_job_cmds=()
    local ingroup_job_script=""
    local ingroup_job_idx=0
    local other_taxa_string="${other_taxa[*]}"

    mkdir -p "$ingroup_job_dir"

    while IFS= read -r ingroup_noref_maf_abs
    do
        [[ -n "$ingroup_noref_maf_abs" ]] || continue

        chrom_dir=$(dirname "$ingroup_noref_maf_abs")
        chrom_name=$(basename "$chrom_dir")
        ref_seq_lookup_name=$(get_reference_lookup_seq_name "$chrom_name" "$reference")
        ref_seq_name=$(clean_reference_seq_name "$chrom_name" "$reference")

        chrom_prefix_abs="$chrom_dir/${output_base}_${ingroup_ancestor_name}"

        ingroup_added_maf_abs="${chrom_prefix_abs}.maf"
        ingroup_hal_abs="${chrom_prefix_abs}_noref.hal"
        ingroup_jc69_root_abs="${chrom_prefix_abs}_JC69modle"
        ingroup_tsv_abs="${chrom_prefix_abs}_JC69modle.tsv"
        ingroup_fasta_abs="${chrom_dir}/${ingroup_ancestor_name}.${ref_seq_name}.fa"
        ingroup_with_sequence_maf_abs="${chrom_prefix_abs}_new.maf"
        ingroup_noref_new_maf_abs="${chrom_prefix_abs}_noref_new.maf"
        ref_plus_ingroup_maf_abs="${chrom_prefix_abs}_new_noIngroup.maf"
        top_added_maf_abs="${chrom_prefix_abs}_new_noIngroup_${top_ancestor_name}.maf"

        # Record expected outputs in the parent shell before launching jobs.
        ingroup_noref_new_maf_files+=( "$ingroup_noref_new_maf_abs" )
        ref_plus_ingroup_maf_files+=( "$ref_plus_ingroup_maf_abs" )
        top_added_maf_files+=( "$top_added_maf_abs" )
        ingroup_fasta_files+=( "$ingroup_fasta_abs" )

        ingroup_job_idx=$((ingroup_job_idx + 1))
        ingroup_job_script=$(printf '%s/job_%04d.sh' "$ingroup_job_dir" "$ingroup_job_idx")

        cat > "$ingroup_job_script" <<EOF
#!/bin/bash
set -euo pipefail

echo "  [ingroup ancestor job ${ingroup_job_idx}] Processing chromosome: ${chrom_name}"

if [[ ! -s "${ingroup_task_jc69_root_abs}.mod" ]]; then
    echo "Error: ingroup ancestor model file not found or empty: ${ingroup_task_jc69_root_abs}.mod" >&2
    exit 1
fi

"$patched_maf2hal" \
    --refGenome "$ingroup_ancestor_name" \
    --targetGenomes "$target_genomes_csv" \
    "$ingroup_noref_maf_abs" \
    "$ingroup_hal_abs"

"$patched_ancestorsML" \
    --printWrites \
    "$ingroup_hal_abs" \
    "$ingroup_ancestor_name" \
    "${ingroup_task_jc69_root_abs}.mod" \
    > "$ingroup_tsv_abs"

env -u LD_LIBRARY_PATH \
    "$patched_halWriteNucleotides" \
    "$ingroup_hal_abs" \
    "$ingroup_tsv_abs" \
    > "${ingroup_jc69_root_abs}_halWriteNucleotides.log"

"$patched_hal2fasta" \
    "$ingroup_hal_abs" \
    "$ingroup_ancestor_name" \
    > "$ingroup_fasta_abs"

# Write the inferred ingroup ancestor sequence back into the
# reference-containing chromosome MAF.
python "$script_dir/replace_ancestor_seq_in_maf.py" \
    "$ingroup_added_maf_abs" \
    "$ingroup_ancestor_name" \
    "$ingroup_fasta_abs" \
    "$ingroup_with_sequence_maf_abs" \
    "$ref_seq_name"

# Also write it back into the reference-excluded chromosome MAF.
python "$script_dir/replace_ancestor_seq_in_maf.py" \
    "$ingroup_noref_maf_abs" \
    "$ingroup_ancestor_name" \
    "$ingroup_fasta_abs" \
    "$ingroup_noref_new_maf_abs" \
    "$ref_seq_name"

# Keep reference + ingroup ancestor only.
LC_ALL=C awk -v taxa="$other_taxa_string" '
BEGIN {
    n=split(taxa,a," ")
    for(i=1;i<=n;i++) {
        if(a[i] != "") remove[a[i]]=1
    }
}
\$1 ~ /^(s|i|e|q)\$/ {
    split(\$2,b,".")
    if(remove[b[1]]) next
}
{print}
' "$ingroup_with_sequence_maf_abs" > "$ref_plus_ingroup_maf_abs"

# Add the top ancestor placeholder row.
python "$script_dir/add_ance2mafwithN.py" \
    "$ref_plus_ingroup_maf_abs" \
    "$top_added_maf_abs" \
    "$top_ancestor_name" \
    "$ref_seq_name"

echo "  [ingroup ancestor job ${ingroup_job_idx}] Finished chromosome: ${chrom_name}"
EOF

        chmod +x "$ingroup_job_script"
        ingroup_job_cmds+=( "bash \"$ingroup_job_script\"" )

    done < <(printf '%s\n' "${ingroup_noref_maf_files[@]}" | sort -V)

    echo "Starting parallel ingroup-ancestor inference."
    echo "  Number of chromosome jobs: ${#ingroup_job_cmds[@]}"
    echo "  Maximum simultaneous ancestor jobs: ${max_ancestor_jobs}"

    if (( ${#ingroup_job_cmds[@]} > 0 )); then
        run_parallel_commands "$max_ancestor_jobs" "${ingroup_job_cmds[@]}"
    fi

    echo "All ingroup-ancestor chromosome jobs finished."

    top_task_added_maf_abs="${output_prefix_abs}_${ingroup_ancestor_name}_new_noIngroup_${top_ancestor_name}.maf"
    merge_maf_files "$top_task_added_maf_abs" "${top_added_maf_files[@]}"

    top_task_model_use_abs="${output_prefix_abs}_${ingroup_ancestor_name}_new_noIngroup_${top_ancestor_name}_$(basename "${model_file%.*}")_useThis.mod"
    top_task_jc69_root_abs="${output_prefix_abs}_${ingroup_ancestor_name}_new_noIngroup_${top_ancestor_name}_JC69modle"
    run_task_level_phylofit_model \
        "$top_model_tree" \
        "$top_task_added_maf_abs" \
        "$top_task_model_use_abs" \
        "$top_task_jc69_root_abs"

    ############################
    # Round 2:
    # Parallel per-chromosome inference of the top ancestor.
    #
    # This round starts only after every Round-1 chromosome job has finished
    # and the shared top-ancestor model has been generated.
    ############################

    local top_job_dir="$work_dir/ancestor_jobs/top_${top_ancestor_name}"
    local -a top_job_cmds=()
    local top_job_script=""
    local top_job_idx=0

    mkdir -p "$top_job_dir"

    while IFS= read -r top_added_maf_abs
    do
        [[ -n "$top_added_maf_abs" ]] || continue

        chrom_dir=$(dirname "$top_added_maf_abs")
        chrom_name=$(basename "$chrom_dir")
        ref_seq_lookup_name=$(get_reference_lookup_seq_name "$chrom_name" "$reference")
        ref_seq_name=$(clean_reference_seq_name "$chrom_name" "$reference")

        chrom_prefix_abs="$chrom_dir/${output_base}_${ingroup_ancestor_name}"

        top_hal_abs="${chrom_prefix_abs}_new_noIngroup_${top_ancestor_name}.hal"
        top_jc69_root_abs="${chrom_prefix_abs}_new_noIngroup_${top_ancestor_name}_JC69modle"
        top_tsv_abs="${chrom_prefix_abs}_new_noIngroup_${top_ancestor_name}_JC69modle.tsv"
        top_fasta_abs="${chrom_dir}/${top_ancestor_name}.${ref_seq_name}.fa"
        top_noingroup_new_maf_abs="${chrom_prefix_abs}_new_noIngroup_${top_ancestor_name}_new.maf"

        # Record expected outputs in the parent shell before launching jobs.
        top_noingroup_new_maf_files+=( "$top_noingroup_new_maf_abs" )
        top_fasta_files+=( "$top_fasta_abs" )

        top_job_idx=$((top_job_idx + 1))
        top_job_script=$(printf '%s/job_%04d.sh' "$top_job_dir" "$top_job_idx")

        cat > "$top_job_script" <<EOF
#!/bin/bash
set -euo pipefail

echo "  [top ancestor job ${top_job_idx}] Processing chromosome: ${chrom_name}"

if [[ ! -s "${top_task_jc69_root_abs}.mod" ]]; then
    echo "Error: top ancestor model file not found or empty: ${top_task_jc69_root_abs}.mod" >&2
    exit 1
fi

"$patched_maf2hal" \
    --refGenome "$top_ancestor_name" \
    "$top_added_maf_abs" \
    "$top_hal_abs"

"$patched_ancestorsML" \
    --printWrites \
    "$top_hal_abs" \
    "$top_ancestor_name" \
    "${top_task_jc69_root_abs}.mod" \
    > "$top_tsv_abs"

env -u LD_LIBRARY_PATH \
    "$patched_halWriteNucleotides" \
    "$top_hal_abs" \
    "$top_tsv_abs" \
    > "${top_jc69_root_abs}_halWriteNucleotides.log"

"$patched_hal2fasta" \
    "$top_hal_abs" \
    "$top_ancestor_name" \
    > "$top_fasta_abs"

python "$script_dir/replace_ancestor_seq_in_maf.py" \
    "$top_added_maf_abs" \
    "$top_ancestor_name" \
    "$top_fasta_abs" \
    "$top_noingroup_new_maf_abs" \
    "$ref_seq_name"

echo "  [top ancestor job ${top_job_idx}] Finished chromosome: ${chrom_name}"
EOF

        chmod +x "$top_job_script"
        top_job_cmds+=( "bash \"$top_job_script\"" )

    done < <(printf '%s\n' "${top_added_maf_files[@]}" | sort -V)

    echo "Starting parallel top-ancestor inference."
    echo "  Number of chromosome jobs: ${#top_job_cmds[@]}"
    echo "  Maximum simultaneous ancestor jobs: ${max_ancestor_jobs}"

    if (( ${#top_job_cmds[@]} > 0 )); then
        run_parallel_commands "$max_ancestor_jobs" "${top_job_cmds[@]}"
    fi

    echo "All top-ancestor chromosome jobs finished."

    local final_ingroup_maf_abs="${output_prefix_abs}_${ingroup_ancestor_name}_noref_new.maf"
    local final_ingroup_checked_maf_abs="${output_prefix_abs}_${ingroup_ancestor_name}_noref_new2.maf"
    local final_ingroup_hal_abs="${output_prefix_abs}_${ingroup_ancestor_name}_noref_new2.hal"
    local final_top_maf_abs="${output_prefix_abs}_${ingroup_ancestor_name}_new_noIngroup_${top_ancestor_name}_new.maf"
    local final_top_checked_maf_abs="${output_prefix_abs}_${ingroup_ancestor_name}_new_noIngroup_${top_ancestor_name}_new2.maf"
    local final_top_hal_abs="${output_prefix_abs}_${ingroup_ancestor_name}_new_noIngroup_${top_ancestor_name}_new.hal"
    local final_ingroup_fasta_abs="${output_dir}/${ingroup_ancestor_name}.fa"
    local final_top_fasta_abs="${output_dir}/${top_ancestor_name}.fa"

    merge_maf_files "$final_ingroup_maf_abs" "${ingroup_noref_new_maf_files[@]}"

    merge_fasta_files "$final_ingroup_fasta_abs" "${ingroup_fasta_files[@]}"
    "$patched_seqkit" fx2tab -n -i -l "$final_ingroup_fasta_abs" > "$info_tmp_dir/${ingroup_ancestor_name}.info.txt"

    cp "$final_ingroup_maf_abs" "$final_ingroup_checked_maf_abs"
    for tax in "${taxa_args[@]}" "$ingroup_ancestor_name"; do
        info_abs="$info_tmp_dir/${tax}.info.txt"
        python "$script_dir/check_if_maf_contain_all_taxachrom.py" "$info_abs" "$final_ingroup_checked_maf_abs" "$tax"
    done

    "$patched_maf2hal" \
        --refGenome "$ingroup_ancestor_name" \
        --targetGenomes "$target_genomes_csv" \
        "$final_ingroup_checked_maf_abs" \
        "$final_ingroup_hal_abs"

    merge_maf_files "$final_top_maf_abs" "${top_noingroup_new_maf_files[@]}"

    merge_fasta_files "$final_top_fasta_abs" "${top_fasta_files[@]}"
    "$patched_seqkit" fx2tab -n -i -l "$final_top_fasta_abs" > "$info_tmp_dir/${top_ancestor_name}.info.txt"

    cp "$final_top_maf_abs" "$final_top_checked_maf_abs"
    for tax in "$reference" "$ingroup_ancestor_name" "$top_ancestor_name"; do
    	info_abs="$info_tmp_dir/${tax}.info.txt"
    	python "$script_dir/check_if_maf_contain_all_taxachrom.py" "$info_abs" "$final_top_checked_maf_abs" "$tax"
    done

    "$patched_maf2hal" \
        --refGenome "$top_ancestor_name" \
        "$final_top_checked_maf_abs" \
        "$final_top_hal_abs"

    cp "$final_top_hal_abs" "$output_hal_abs"

    "$patched_halAppendSubtree" \
        "$output_hal_abs" \
        "$final_ingroup_hal_abs" \
        "$ingroup_ancestor_name" \
        "$ingroup_ancestor_name" \
        --merge

    generated_branch_summaries+=(
        "reference outgroup; final HAL contains top and ingroup ancestor levels|${top_ancestor_name} + ${ingroup_ancestor_name}|${output_hal_abs}|NA|NA|NA|${final_top_fasta_abs},${final_ingroup_fasta_abs}"
    )
}

generated_branch_summaries=()

if [[ "$ifRefNonOutgroup" == "1" ]]; then
    # Reference is treated as an outgroup:
    #   (reference,(non-reference taxa)Anc_name_ofIngroup)Anc_name
    run_per_chrom_ref_outgroup_branch
else
    # Reference is a true direct alignment element, not an outgroup.
    run_per_chrom_direct_ancestor_branch
fi

############################
# 17. Delete or preserve intermediate temporary workspace
############################
echo "Pipeline finished."
echo "Final consensus MAF output: $output_maf_abs"
echo "Final consensus FASTA output: $output_fasta_abs"
if [[ -f "$output_hal_abs" ]]; then
    echo "Final consensus HAL output: $output_hal_abs"
fi

if [[ "$ifRefNonOutgroup" == "1" ]]; then
    if [[ -f "${output_dir}/${anc_name_of_ingroup}.fa" ]]; then
        echo "Final ingroup ancestor FASTA output: ${output_dir}/${anc_name_of_ingroup}.fa"
    fi
    if [[ -f "${output_dir}/${anc_name}.fa" ]]; then
        echo "Final top ancestor FASTA output: ${output_dir}/${anc_name}.fa"
    fi
else
    if [[ -f "${output_dir}/${anc_name}.fa" ]]; then
        echo "Final ancestor FASTA output: ${output_dir}/${anc_name}.fa"
    fi
fi

if [[ "$ifDeleteImmdiFiles" == "0" ]]; then
    rm -rf "$work_dir"
    rm -f "${output_prefix_abs}"_*
    echo "Temporary workspace deleted: $work_dir"
else
    echo "Temporary workspace retained: $work_dir"

    if [[ "$ifRefNonOutgroup" == "1" ]]; then
        echo "Except for ${output_base}.maf ${output_base}.fasta ${output_base}.hal ${anc_name_of_ingroup}.fa ${anc_name}.fa, all other files and directories are intermediate files."
    else
        echo "Except for ${output_base}.maf ${output_base}.fasta ${output_base}.hal ${anc_name}.fa, all other files and directories are intermediate files."
    fi
fi
