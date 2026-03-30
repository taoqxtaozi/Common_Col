#!/bin/bash
set -euo pipefail
shopt -s nullglob

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

patched_halWriteNucleotides="$script_dir/tool_used/bin/halWriteNucleotides"
patched_phyloFit="$script_dir/tool_used/bin/phyloFit"
patched_ancestorsML="$script_dir/tool_used/bin/ancestorsML"
patched_maf2hal="$script_dir/tool_used/bin/maf2hal"
patched_hal2fasta="$script_dir/tool_used/bin/hal2fasta"
patched_seqkit="$script_dir/tool_used/bin/seqkit"

for exe in \
    "$patched_halWriteNucleotides" \
    "$patched_phyloFit" \
    "$patched_ancestorsML" \
    "$patched_maf2hal" \
    "$patched_hal2fasta" \
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
  bash RunPipelineUseThis.sh --input /path/to/aln1.maf /path/to/aln2.maf ... --reference REF --pre output_prefix --threads 60 --common_workers 15 --genome_path genome.txt --Anc_name AncX --ModelFile tryMLstartree.mod [options]

Required arguments:
  --input                 Input MAF files (space-separated).
  --reference             Reference taxon name. It will be placed first in the final FASTA output.
  --pre                   Output prefix. Two files will be generated first: <pre>.maf and <pre>.fasta
  --threads               Total CPU budget available on the machine.
  --common_workers        CPUs used by one consensus-extraction program.
  --genome_path           Two-column file: TAXON GENOME_PATH. Extra taxa not present in the final alignment will be ignored.
  --Anc_name              Name of the consensus ancestor to be inferred, e.g. AncX.
  --ModelFile             Path to the base model file used to generate <model>_useThis.mod for phyloFit/ancestorsML.

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
      - maf_to_concat_fasta.py
      - add_ance2mafwithN.py
      - check_if_maf_contain_all_taxachrom.py
      - sep_file/1.sh
  * All intermediate files and directories will be created under a tmp_XXXXXX directory
    inside the directory containing --pre.
  * threads must be >= common_workers.
  * threads is the total CPU budget on the machine.
  * common_workers is the CPU usage of one consensus-extraction program.
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
reference=""
work_dir=""
sep_length=218505
chrom_length_threshold=15819469
global_num=40
separate_workers=4
common_workers=""
last_dir=""
genome_path_file=""
anc_name=""
model_file=""
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

if [[ -z "$anc_name" ]]; then
    echo "Error: --Anc_name is required."
    exit 1
fi

if [[ -z "$model_file" ]]; then
    echo "Error: --ModelFile is required."
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

if [[ ! -f "$genome_path_file" ]]; then
    echo "Error: genome path file not found: $genome_path_file"
    exit 1
fi

if [[ ! -f "$model_file" ]]; then
    echo "Error: model file not found: $model_file"
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

if [[ ! -f "$script_dir/sep_file/1.sh" ]]; then
    echo "Error: sep_file/1.sh not found in script directory."
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
mkdir -p "$work_dir/sep"

############################
# 8. Use the directory generated from the last input file as the template,
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

            cat > "${n}.sh" <<EOF2
#!/bin/bash
set -euo pipefail
mkdir -p "../${head_chrom}"
python "$script_dir/seperate_maffile_by_index.py" "../${head_chrom}.maf" "${n}" "${start}" "${end}" "../${head_chrom}"
EOF2
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
            for j in "${!inputs[@]}"
            do
                d="${input_tags[$j]}"
                input_args="${input_args} \"$work_dir/${d}/${head_chrom}/${n}_${start}_${end}.maf\""
            done

            echo "python \"$script_dir/Extract_consensus_multi_alns.py\" --input${input_args} --output \"$work_dir/common_output_files/${head_chrom}/part${n}.maf\" --global_num ${global_num} --global_start ${start} --global_end ${end} --separate_workers ${separate_workers} --common_workers ${common_workers}" \
                > "$work_dir/common_command_files/${head_chrom}/${n}.sh"
        done

    else
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

cd "$work_dir" || exit 1

############################
# 9. Copy the generated sep_<chrom> directories under sep/
#    into each input directory
############################
for j in "${!inputs[@]}"
do
    dir="${input_tags[$j]}"

    for sep_dir in "$work_dir"/sep/*
    do
        [ -d "$sep_dir" ] || continue
        base_sep=$(basename "$sep_dir")
        rm -rf "$work_dir/$dir/$base_sep"
        cp -r "$sep_dir" "$work_dir/$dir"/
    done
done

############################
# 10. Execute the chunking scripts in parallel
#     Split jobs are treated as ~1 CPU jobs, so up to threads jobs may run in parallel.
############################
split_job_cmds=()

for j in "${!inputs[@]}"
do
    dir="${input_tags[$j]}"

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

if (( ${#split_job_cmds[@]} > 0 )); then
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

mapfile -t other_taxa < <(sort "$taxa_tmp")
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

for tax in "${taxa_args[@]}"; do
    genome_abs="${genome_paths[$tax]}"
    info_abs="$info_tmp_dir/${tax}.info.txt"
    "$patched_seqkit" fx2tab -n -i -l "$genome_abs" > "$info_abs"
done

############################
# 16. Add ancestor row names to the final consensus MAF
############################
python "$script_dir/add_ance2mafwithN.py" "$output_maf_abs" "$ancestor_maf_abs" "$anc_name"

############################
# 17. Check chromosome coverage for each aligned taxon against the ancestor-augmented MAF
############################
for tax in "${taxa_args[@]}"; do
    info_abs="$info_tmp_dir/${tax}.info.txt"
    python "$script_dir/check_if_maf_contain_all_taxachrom.py" "$info_abs" "$ancestor_maf_abs" "$tax"
done

############################
# 18. Convert ancestor-augmented MAF to HAL
############################
if (( ${#other_taxa[@]} == 0 )); then
    echo "Error: no non-reference taxa found in the final MAF; cannot build targetGenomes or ancestor model."
    rm -rf "$work_dir"
    exit 1
fi

target_genomes_csv=$(IFS=,; echo "${other_taxa[*]}")
"$patched_maf2hal" --refGenome "$anc_name" --targetGenomes "$target_genomes_csv" "$ancestor_maf_abs" "$ancestor_hal_abs"

############################
# 19. Build the temporary model file for phyloFit/ancestorsML
############################
other_taxa_space="${other_taxa[*]}"
TreeInModel=$(awk -v t="$other_taxa_space" 'BEGIN{
    n=split(t,a," ")
    s=a[1]":1"
    for(i=2;i<=n;i++){
        edge=(i<n?":0!":":0")
        s="(" s "," a[i]":1)" edge
    }
    print s
}')

cp "$model_file" "$model_use_abs"
printf 'TREE: (%s:1,%s);\n' "$reference" "$TreeInModel" >> "$model_use_abs"

############################
# 20. Infer ancestral states and write ancestor FASTA
############################
"$patched_phyloFit" --init-model "$model_use_abs" --subst-mod JC69 --out-root "$jc69_root_abs" --msa-format MAF "$ancestor_maf_abs"
"$patched_ancestorsML" --printWrites "$ancestor_hal_abs" "$anc_name" "${jc69_root_abs}.mod" > "$ancestor_tsv_abs"
env -u LD_LIBRARY_PATH "$patched_halWriteNucleotides" "$ancestor_hal_abs" "$ancestor_tsv_abs" > "$work_dir/halWriteNucleotides.log"
"$patched_hal2fasta" "$ancestor_hal_abs" "$anc_name" > "$ancestor_fasta_abs"

############################
# 21. Delete intermediate temporary workspace
############################
rm -rf "$work_dir"

echo "Pipeline finished."
echo "Final MAF output: $output_maf_abs"
echo "Final FASTA output: $output_fasta_abs"
echo "Intermediate ancestor-augmented MAF (retained for record; not a final output): $ancestor_maf_abs"
echo "Ancestor HAL: $ancestor_hal_abs"
echo "Temporary model file: $model_use_abs"
echo "phyloFit root: $jc69_root_abs"
echo "Ancestor TSV: $ancestor_tsv_abs"
echo "Ancestor FASTA: $ancestor_fasta_abs"
