#!/bin/bash
set -euo pipefail

RUN_PATH=/picb/phylcomb/taoqixin/cactus-bin-v2.9.2/trys/try1
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

root_alltaxa_hal="${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.allTaxa.hal"
root_consensus_maf="${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.maf"
target_genomes_csv=$(awk '$1=="s"{split($2,a,"."); if(!(a[1] in seen)){seen[a[1]]=1; order[++n]=a[1]}} END{for(i=1;i<=n;i++) printf "%s%s", order[i], (i<n?",":"")}' "$root_consensus_maf")
"$patched_maf2hal" --refGenome Anc_root --targetGenomes "$target_genomes_csv" "$root_consensus_maf" "$root_alltaxa_hal"
"$patched_halAppendSubtree" ${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.allTaxa.hal ${RUN_PATH}/Anc_root/Anc_in1/Anc_in1_aln1.hal Anc_in1 Anc_in1 --merge
"$patched_halAppendSubtree" ${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.allTaxa.hal ${RUN_PATH}/Anc_root/Anc_in3/Anc_in3_aln1.hal Anc_in3 Anc_in3 --merge
"$patched_halAppendSubtree" ${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.allTaxa.hal ${RUN_PATH}/Anc_root/Anc_in7/Anc_in7_consensus_Anc_in7.hal Anc_in7 Anc_in7 --merge
"$patched_halAppendSubtree" ${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.allTaxa.hal ${RUN_PATH}/Anc_root/Anc_in7/Anc_in4/Anc_in4_aln1.hal Anc_in4 Anc_in4 --merge
"$patched_halAppendSubtree" ${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.allTaxa.hal ${RUN_PATH}/Anc_root/Anc_in7/Anc_in6/Anc_in6_aln1.hal Anc_in6 Anc_in6 --merge
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/Anc_root/jobstorehal2maf_Anc_root_allTaxa ${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.allTaxa.hal ${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.allTaxa.maf
taxa_args=(GALGA ACACH APAVI CARCR CATAU COLST CORBR FALPE GEOFO HALLE LEPDI MELUN NESNO PICPU TAEGU TYTAL)
python /picb/phylcomb/taoqixin/software/cactus-bin-v2.9.2/GALGA_APAVI_CATAU_LEPDI/code0328/maf_to_concat_fasta.py "${taxa_args[@]}" < ${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.allTaxa.maf > ${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.allTaxa.fasta
rm -f "$root_consensus_maf" "$root_alltaxa_hal"
mkdir -p ${RUN_PATH}/final_output
mv ${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.allTaxa.maf ${RUN_PATH}/final_output/Anc_root_consensus_Anc_root.allTaxa.maf
mv ${RUN_PATH}/Anc_root/Anc_root_consensus_Anc_root.allTaxa.fasta ${RUN_PATH}/final_output/Anc_root_consensus_Anc_root.allTaxa.fasta
mv ${RUN_PATH}/Anc_root/Anc_root.fa ${RUN_PATH}/final_output/Anc_root.fa
