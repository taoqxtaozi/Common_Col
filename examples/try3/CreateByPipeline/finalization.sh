#!/bin/bash
set -euo pipefail

RUN_PATH=/picb/phylcomb/taoqixin/cactus-bin-v2.9.2/evolverSimControl_100000long_43taxa_processLostChrom_usetaxadirectly/0530reroot/binary
root_alltaxa_hal="${RUN_PATH}/Root/FinalResultAlignment.hal"
# Copy the root-level primary HAL to a stable final-working name.
cp ${RUN_PATH}/Root/Root_aln1.hal "$root_alltaxa_hal"
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/Root/jobstorehal2maf_Root_allTaxa ${RUN_PATH}/Root/FinalResultAlignment.hal ${RUN_PATH}/Root/FinalResultAlignment.maf
taxa_args=(GALGA CATAU CHLUN COLLI MESUN)
python /picb/phylcomb/taoqixin/cactus-bin-v2.9.2/code0530reroot/maf_to_concat_fasta.py "${taxa_args[@]}" < ${RUN_PATH}/Root/FinalResultAlignment.maf > ${RUN_PATH}/Root/FinalResultAlignment.fasta
mkdir -p ${RUN_PATH}/final_output
mv "$root_alltaxa_hal" ${RUN_PATH}/final_output/FinalResultAlignment.hal
mv ${RUN_PATH}/Root/FinalResultAlignment.maf ${RUN_PATH}/final_output/FinalResultAlignment.maf
mv ${RUN_PATH}/Root/FinalResultAlignment.fasta ${RUN_PATH}/final_output/FinalResultAlignment.fasta
mv ${RUN_PATH}/Root/Root.fa ${RUN_PATH}/final_output/Root.fa
