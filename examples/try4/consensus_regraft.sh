#!/bin/bash
set -euo pipefail

RUN_PATH=/picb/phylcomb/taoqixin/cactus-bin-v2.9.2/trys/try4
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

root_alltaxa_hal="${RUN_PATH}/Anc_root/Anc_root_aln1.allTaxa.hal"
cp ${RUN_PATH}/Anc_root/Anc_root_aln1.hal "$root_alltaxa_hal"
"$patched_halAppendSubtree" ${RUN_PATH}/Anc_root/Anc_root_aln1.allTaxa.hal ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_consensus_Anc_in2.hal Anc_in2 Anc_in2 --merge
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/Anc_root/jobstorehal2maf_Anc_root_allTaxa ${RUN_PATH}/Anc_root/Anc_root_aln1.allTaxa.hal ${RUN_PATH}/Anc_root/Anc_root_aln1.allTaxa.maf
taxa_args=(GALGA APAVI CATAU COLST HALLE LEPDI PICPU TYTAL)
python /picb/phylcomb/taoqixin/software/cactus-bin-v2.9.2/GALGA_APAVI_CATAU_LEPDI/code0328/maf_to_concat_fasta.py "${taxa_args[@]}" < ${RUN_PATH}/Anc_root/Anc_root_aln1.allTaxa.maf > ${RUN_PATH}/Anc_root/Anc_root_aln1.allTaxa.fasta
rm -f "$root_alltaxa_hal"
mkdir -p ${RUN_PATH}/final_output
mv ${RUN_PATH}/Anc_root/Anc_root_aln1.allTaxa.maf ${RUN_PATH}/final_output/Anc_root_aln1.allTaxa.maf
mv ${RUN_PATH}/Anc_root/Anc_root_aln1.allTaxa.fasta ${RUN_PATH}/final_output/Anc_root_aln1.allTaxa.fasta
mv ${RUN_PATH}/Anc_root/Anc_root.fa ${RUN_PATH}/final_output/Anc_root.fa
