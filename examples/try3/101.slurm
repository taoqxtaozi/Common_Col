#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH -J tqx-galga_LEPDI
#SBATCH --mem 60G
#SBATCH -p hpc
#SBATCH -o ./log10
#SBATCH -e ./error10
#SBATCH -A phylcomb

## Copy and run all lines below that do not start with ##
## Command0
RUN_PATH=/picb/phylcomb/taoqixin/cactus-bin-v2.9.2/trys/try3

## Run the following commands:
## Commands in the same block can be run in parallel if they are marked so.
## You may need to adjust Cactus resource parameters yourself, e.g. --defaultMemory --maxCores --maxMemory.

## Command1  [Anc1; fixed single-run]
cactus --root Anc1 ${RUN_PATH}/ANC0/Anc1/jobstorecluster_Anc1_aln1 ${RUN_PATH}/ANC0/Anc1/Anc1_aln1.txt ${RUN_PATH}/ANC0/Anc1/Anc1_aln1.hal
hal2fasta ${RUN_PATH}/ANC0/Anc1/Anc1_aln1.hal Anc1 > ${RUN_PATH}/ANC0/Anc1/Anc1.fa

## Level 2: Command2~Command5 can be run in parallel
## Command2  [ANC0; guide-tree-specific alignment]
cactus ${RUN_PATH}/ANC0/jobstorecluster_ANC0_aln1 ${RUN_PATH}/ANC0/ANC0_aln1.txt ${RUN_PATH}/ANC0/ANC0_aln1.hal
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/ANC0/jobstorehal2maf_ANC0_aln1 ${RUN_PATH}/ANC0/ANC0_aln1.hal ${RUN_PATH}/ANC0/ANC0_aln1.maf

## Command3  [ANC0; guide-tree-specific alignment]
cactus ${RUN_PATH}/ANC0/jobstorecluster_ANC0_aln2 ${RUN_PATH}/ANC0/ANC0_aln2.txt ${RUN_PATH}/ANC0/ANC0_aln2.hal
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/ANC0/jobstorehal2maf_ANC0_aln2 ${RUN_PATH}/ANC0/ANC0_aln2.hal ${RUN_PATH}/ANC0/ANC0_aln2.maf

## Command4  [ANC0; guide-tree-specific alignment]
cactus ${RUN_PATH}/ANC0/jobstorecluster_ANC0_aln3 ${RUN_PATH}/ANC0/ANC0_aln3.txt ${RUN_PATH}/ANC0/ANC0_aln3.hal
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/ANC0/jobstorehal2maf_ANC0_aln3 ${RUN_PATH}/ANC0/ANC0_aln3.hal ${RUN_PATH}/ANC0/ANC0_aln3.maf

## Command5  [ANC0; guide-tree-specific alignment]
cactus ${RUN_PATH}/ANC0/jobstorecluster_ANC0_aln4 ${RUN_PATH}/ANC0/ANC0_aln4.txt ${RUN_PATH}/ANC0/ANC0_aln4.hal
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/ANC0/jobstorehal2maf_ANC0_aln4 ${RUN_PATH}/ANC0/ANC0_aln4.hal ${RUN_PATH}/ANC0/ANC0_aln4.maf

## After Command2~Command5 are done, run Command6
## Command6  [ANC0; extract consensus and infer ancestor]
bash /picb/phylcomb/taoqixin/software/cactus-bin-v2.9.2/GALGA_APAVI_CATAU_LEPDI/code0328/RunPipelineUseThis.sh \
  --input ${RUN_PATH}/ANC0/ANC0_aln1.maf ${RUN_PATH}/ANC0/ANC0_aln2.maf ${RUN_PATH}/ANC0/ANC0_aln3.maf ${RUN_PATH}/ANC0/ANC0_aln4.maf \
  --reference GALGA \
  --pre ${RUN_PATH}/ANC0/ANC0_consensus \
  --threads 15 \
  --common_workers 5 \
  --genome_path ${RUN_PATH}/ANC0/ANC0.genome.txt \
  --Anc_name ANC0 \
  --ModelFile /picb/phylcomb/taoqixin/software/cactus-bin-v2.9.2/GALGA_APAVI_CATAU_LEPDI/code0328/tryMLstartree.mod

## After Command6 is done, run descendant-HAL regrafting:
bash ${RUN_PATH}/consensus_regraft.sh