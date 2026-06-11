#!/bin/bash
set -euo pipefail

RUN_PATH=/picb/phylcomb/taoqixin/cactus-bin-v2.9.2/evolverSimControl_100000long_43taxa_processLostChrom_usetaxadirectly/0530reroot/try2
patched_halAppendSubtree=/picb/phylcomb/taoqixin/cactus-bin-v2.9.2/code0530reroot/tool_used/bin/halAppendSubtree

for exe in "$patched_halAppendSubtree"; do
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

root_alltaxa_hal="${RUN_PATH}/Root/FinalResultAlignment.hal"
# Copy the root-level primary HAL to a stable final-working name.
cp ${RUN_PATH}/Root/consensus.hal "$root_alltaxa_hal"
"$patched_halAppendSubtree" ${RUN_PATH}/Root/FinalResultAlignment.hal ${RUN_PATH}/Root/interior9/consensus.hal interior9 interior9 --merge
"$patched_halAppendSubtree" ${RUN_PATH}/Root/FinalResultAlignment.hal ${RUN_PATH}/Root/interior9/interior10/consensus.hal interior10 interior10 --merge
"$patched_halAppendSubtree" ${RUN_PATH}/Root/FinalResultAlignment.hal ${RUN_PATH}/Root/interior9/interior10/Anc_in2/consensus.hal Anc_in2 Anc_in2 --merge
"$patched_halAppendSubtree" ${RUN_PATH}/Root/FinalResultAlignment.hal ${RUN_PATH}/Root/interior9/interior10/Anc_in2/interior26/interior26_aln1.hal interior26 interior26 --merge
"$patched_halAppendSubtree" ${RUN_PATH}/Root/FinalResultAlignment.hal ${RUN_PATH}/Root/interior9/interior10/Anc_in2/interior25/interior25_aln1.hal interior25 interior25 --merge
"$patched_halAppendSubtree" ${RUN_PATH}/Root/FinalResultAlignment.hal ${RUN_PATH}/Root/interior9/interior10/interior12/interior12_aln1.hal interior12 interior12 --merge
"$patched_halAppendSubtree" ${RUN_PATH}/Root/FinalResultAlignment.hal ${RUN_PATH}/Root/interior9/interior10/interior15/consensus.hal interior15 interior15 --merge
"$patched_halAppendSubtree" ${RUN_PATH}/Root/FinalResultAlignment.hal ${RUN_PATH}/Root/interior9/interior2/interior2_aln1.hal interior2 interior2 --merge
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome TAEGU --noAncestors ${RUN_PATH}/Root/jobstorehal2maf_Root_allTaxa ${RUN_PATH}/Root/FinalResultAlignment.hal ${RUN_PATH}/Root/FinalResultAlignment.maf
taxa_args=(TAEGU ACACH APAVI BUCRH CARCR CATAU COLLI COLST CORBR FALPE GALGA GEOFO HALLE LEPDI MANVI MELUN MERNU MESUN NESNO OPHHO PHORU PICPU PODCR PTEGU TYTAL)
python /picb/phylcomb/taoqixin/cactus-bin-v2.9.2/code0530reroot/maf_to_concat_fasta.py "${taxa_args[@]}" < ${RUN_PATH}/Root/FinalResultAlignment.maf > ${RUN_PATH}/Root/FinalResultAlignment.fasta
mkdir -p ${RUN_PATH}/final_output
mv "$root_alltaxa_hal" ${RUN_PATH}/final_output/FinalResultAlignment.hal
mv ${RUN_PATH}/Root/FinalResultAlignment.maf ${RUN_PATH}/final_output/FinalResultAlignment.maf
mv ${RUN_PATH}/Root/FinalResultAlignment.fasta ${RUN_PATH}/final_output/FinalResultAlignment.fasta
mv ${RUN_PATH}/Root/Root.fa ${RUN_PATH}/final_output/Root.fa
