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
RUN_PATH=/picb/phylcomb/taoqixin/cactus-bin-v2.9.2/trys/try4

## Run the following commands:
## Commands in the same block can be run in parallel if they are marked so.
## You may need to adjust Cactus resource parameters yourself, e.g. --defaultMemory --maxCores --maxMemory.

## Level 1: Command1~Command4 can be run in parallel
## Command1  [Anc_in2; guide-tree-specific alignment]
cactus ${RUN_PATH}/Anc_root/Anc_in2/jobstorecluster_Anc_in2_aln1 ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln1.txt ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln1.hal
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/Anc_root/Anc_in2/jobstorehal2maf_Anc_in2_aln1 ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln1.hal ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln1.maf

## Command2  [Anc_in2; guide-tree-specific alignment]
cactus ${RUN_PATH}/Anc_root/Anc_in2/jobstorecluster_Anc_in2_aln2 ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln2.txt ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln2.hal
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/Anc_root/Anc_in2/jobstorehal2maf_Anc_in2_aln2 ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln2.hal ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln2.maf

## Command3  [Anc_in2; guide-tree-specific alignment]
cactus ${RUN_PATH}/Anc_root/Anc_in2/jobstorecluster_Anc_in2_aln3 ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln3.txt ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln3.hal
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/Anc_root/Anc_in2/jobstorehal2maf_Anc_in2_aln3 ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln3.hal ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln3.maf

## Command4  [Anc_in2; guide-tree-specific alignment]
cactus ${RUN_PATH}/Anc_root/Anc_in2/jobstorecluster_Anc_in2_aln4 ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln4.txt ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln4.hal
cactus-hal2maf --dupeMode single --chunkSize 500000 --refGenome GALGA --noAncestors ${RUN_PATH}/Anc_root/Anc_in2/jobstorehal2maf_Anc_in2_aln4 ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln4.hal ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln4.maf

## After Command1~Command4 are done, run Command5
## Command5  [Anc_in2; extract consensus and infer ancestor]
bash /picb/phylcomb/taoqixin/software/cactus-bin-v2.9.2/GALGA_APAVI_CATAU_LEPDI/code0328/RunPipelineUseThis.sh \
  --input ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln1.maf ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln2.maf ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln3.maf ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_aln4.maf \
  --reference GALGA \
  --pre ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2_consensus \
  --threads 15 \
  --common_workers 5 \
  --genome_path ${RUN_PATH}/Anc_root/Anc_in2/Anc_in2.genome.txt \
  --Anc_name Anc_in2 \
  --ModelFile /picb/phylcomb/taoqixin/software/cactus-bin-v2.9.2/GALGA_APAVI_CATAU_LEPDI/code0328/tryMLstartree.mod

## Command6  [Anc_root; fixed single-run]
cactus --root Anc_root ${RUN_PATH}/Anc_root/jobstorecluster_Anc_root_aln1 ${RUN_PATH}/Anc_root/Anc_root_aln1.txt ${RUN_PATH}/Anc_root/Anc_root_aln1.hal
hal2fasta ${RUN_PATH}/Anc_root/Anc_root_aln1.hal Anc_root > ${RUN_PATH}/Anc_root/Anc_root.fa

## After Command6 is done, run descendant-HAL regrafting:
bash ${RUN_PATH}/consensus_regraft.sh

## The planning step is done. The final root-level workspace is under:
## ${RUN_PATH}/Anc_root

############################################################
# Final output summary
# RUN_PATH = /picb/phylcomb/taoqixin/cactus-bin-v2.9.2/trys/try4
# Final consensus MAF   : /picb/phylcomb/taoqixin/cactus-bin-v2.9.2/trys/try4/final_output/Anc_root_aln1.allTaxa.maf
# Final consensus FASTA : /picb/phylcomb/taoqixin/cactus-bin-v2.9.2/trys/try4/final_output/Anc_root_aln1.allTaxa.fasta
# Final inferred genome FASTA for Anc_root: /picb/phylcomb/taoqixin/cactus-bin-v2.9.2/trys/try4/final_output/Anc_root.fa
############################################################
