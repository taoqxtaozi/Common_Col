#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -J tqx-galga_LEPDI
#SBATCH --mem 60G
#SBATCH -p hpc
#SBATCH -o ./log10
#SBATCH -e ./error10
#SBATCH -A phylcomb

## Command0
RUN_PATH=/picb/phylcomb/taoqixin/cactus-bin-v2.9.2/trys/try2

## Run the following commands:
## Commands in the same block can be run in parallel if they are marked so.
## You may need to adjust Cactus resource parameters yourself, e.g. --defaultMemory --maxCores --maxMemory.

## Level 1: Command1~Command3 can be run in parallel
## Command1  [Anc_root; guide-tree-specific alignment]
cactus ${RUN_PATH}/Anc_root/jobstorecluster_Anc_root_aln1 ${RUN_PATH}/Anc_root/Anc_root_aln1.txt ${RUN_PATH}/Anc_root/Anc_root_aln1.hal
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/Anc_root/jobstorehal2maf_Anc_root_aln1 ${RUN_PATH}/Anc_root/Anc_root_aln1.hal ${RUN_PATH}/Anc_root/Anc_root_aln1.maf

## Command2  [Anc_root; guide-tree-specific alignment]
cactus ${RUN_PATH}/Anc_root/jobstorecluster_Anc_root_aln2 ${RUN_PATH}/Anc_root/Anc_root_aln2.txt ${RUN_PATH}/Anc_root/Anc_root_aln2.hal
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/Anc_root/jobstorehal2maf_Anc_root_aln2 ${RUN_PATH}/Anc_root/Anc_root_aln2.hal ${RUN_PATH}/Anc_root/Anc_root_aln2.maf

## Command3  [Anc_root; guide-tree-specific alignment]
cactus ${RUN_PATH}/Anc_root/jobstorecluster_Anc_root_aln3 ${RUN_PATH}/Anc_root/Anc_root_aln3.txt ${RUN_PATH}/Anc_root/Anc_root_aln3.hal
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/Anc_root/jobstorehal2maf_Anc_root_aln3 ${RUN_PATH}/Anc_root/Anc_root_aln3.hal ${RUN_PATH}/Anc_root/Anc_root_aln3.maf

## After Command1~Command3 are done, run Command4
## Command4  [Anc_root; extract consensus and infer ancestor]
bash /picb/phylcomb/taoqixin/software/cactus-bin-v2.9.2/GALGA_APAVI_CATAU_LEPDI/code0328/RunPipelineUseThis.sh \
  --input ${RUN_PATH}/Anc_root/Anc_root_aln1.maf ${RUN_PATH}/Anc_root/Anc_root_aln2.maf ${RUN_PATH}/Anc_root/Anc_root_aln3.maf \
  --reference GALGA \
  --pre ${RUN_PATH}/Anc_root/Anc_root_consensus \
  --threads 15 \
  --common_workers 5 \
  --genome_path ${RUN_PATH}/Anc_root/Anc_root.genome.txt \
  --Anc_name Anc_root \
  --ModelFile /picb/phylcomb/taoqixin/software/cactus-bin-v2.9.2/GALGA_APAVI_CATAU_LEPDI/code0328/tryMLstartree.mod

## The planning step is done. The final root-level workspace is under:
## ${RUN_PATH}/Anc_root