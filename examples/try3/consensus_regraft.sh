#!/bin/bash
set -euo pipefail

RUN_PATH=/picb/phylcomb/taoqixin/cactus-bin-v2.9.2/trys/try3
patched_halAppendSubtree=/picb/phylcomb/taoqixin/software/cactus-bin-v2.9.2/GALGA_APAVI_CATAU_LEPDI/code0328/tool_used/bin/halAppendSubtree
patched_maf2hal=/picb/phylcomb/taoqixin/software/cactus-bin-v2.9.2/GALGA_APAVI_CATAU_LEPDI/code0328/tool_used/bin/maf2hal

for exe in "$patched_halAppendSubtree" "$patched_maf2hal"; do
    if [[ ! -f "$exe" ]]; then
        echo "Error: bundled executable not found: $exe"
        exit 1
    fi
    chmod +x "$exe"
    if [[ ! -x "$exe" ]]; then
        echo "Error: bundled executable is still not executable: $exe"
        exit 1
    fi
done

root_alltaxa_hal="${RUN_PATH}/ANC0/ANC0_consensus_ANC0.allTaxa.hal"
root_consensus_maf="${RUN_PATH}/ANC0/ANC0_consensus_ANC0.maf"
target_genomes_csv=$(awk '$1=="s"{split($2,a,"."); if(!(a[1] in seen)){seen[a[1]]=1; order[++n]=a[1]}} END{for(i=1;i<=n;i++) printf "%s%s", order[i], (i<n?",":"")}' "$root_consensus_maf")
"$patched_maf2hal" --refGenome ANC0 --targetGenomes "$target_genomes_csv" "$root_consensus_maf" "$root_alltaxa_hal"
"$patched_halAppendSubtree" ${RUN_PATH}/ANC0/ANC0_consensus_ANC0.allTaxa.hal ${RUN_PATH}/ANC0/Anc1/Anc1_aln1.hal Anc1 Anc1 --merge
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/ANC0/jobstorehal2maf_ANC0_allTaxa ${RUN_PATH}/ANC0/ANC0_consensus_ANC0.allTaxa.hal ${RUN_PATH}/ANC0/ANC0_consensus_ANC0.allTaxa.maf
taxa_args=(GALGA APAVI CATAU COLST LEPDI TYTAL)
python /picb/phylcomb/taoqixin/software/cactus-bin-v2.9.2/GALGA_APAVI_CATAU_LEPDI/code0328/maf_to_concat_fasta.py "${taxa_args[@]}" < ${RUN_PATH}/ANC0/ANC0_consensus_ANC0.allTaxa.maf > ${RUN_PATH}/ANC0/ANC0_consensus_ANC0.allTaxa.fasta
rm -f "$root_consensus_maf" "$root_alltaxa_hal"
mkdir -p ${RUN_PATH}/final_output
mv ${RUN_PATH}/ANC0/ANC0_consensus_ANC0.allTaxa.maf ${RUN_PATH}/final_output/ANC0_consensus_ANC0.allTaxa.maf
mv ${RUN_PATH}/ANC0/ANC0_consensus_ANC0.allTaxa.fasta ${RUN_PATH}/final_output/ANC0_consensus_ANC0.allTaxa.fasta
mv ${RUN_PATH}/ANC0/ANC0.fa ${RUN_PATH}/final_output/ANC0.fa
