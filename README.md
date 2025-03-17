# Common-col extract common column alignments based on coordiantes of several MAF(multiple alignment format) files

Need python >= 3.6
Need to install packages: biopython interval3 collections multiprocessing
```
pip install biopython interval3 collections multiprocessing
```
#### Some prerequisites
If the maf alignments are generated from cactus```(https://github.com/ComparativeGenomicsToolkit/cactus)```, you can directly run step2.

1. Make sure every blocks in all maf files you input have a reference, make it be the first row in every block. If you need to adjust the order of taxa(say, A,B,C,D and A is the reference) in blocks, you can use:
   
   ```
   conda install bioconda::phast
   maf_parse -O A,B,C,D input1.maf > input1_order.maf
   ```
2. Make sure the input maf files have only one same chromosome of reference, like if the input maf1, maf2, maf3, the start of every their blocks can only be A.chr1, or A chr2..., cannot be both A.chr1 and A.chr2. You can use the file ```seperate_maf.py``` to seperate the maf files.

3. Adjust the blocks to make the reference strand be positive and order the blocks by coordinates of reference. You can use ```https://github.com/dentearl/mafTools```
   ```
   conda create -n py2 python=2.7
   conda activate py2
   conda install genomedk::maftools
   mafStrander --maf alignment_Achr1.maf --seq A.chr1 --strand + > alignment_Achr1positive.maf
   mafSorter --maf alignment_Achr1positive --seq A.chr1 > alignment_Achr1positive_sorted.maf
   conda deactivate
   ```
#### Ready to use

Now the files can be used to extracted. 
```
length=`grep -m1 "^s " maf1 | awk '{print $6}'`
python /mnt/c/e/jarvis/assembly/command/linux/try_common_col_withDup_multi_alns_newversion.py --input maf1 maf2 maf3 ... --output common.maf --global_num 5 --global_start 0 --global_end ${length}
```
If the chromosome length of the reference is to long (say, > 10925261), it would be better to seperate the input maf files equally into several parts, divided by ```k```, which can be set by yourself.
You can use ```seperate_maffile_by_index.py```, the whole code is
```
input_maf_files=(maf1 maf2 maf3...)
k=2185052
length=`grep -m1 "^s " maf1 | awk '{print $6}'`
num=$((${length}/${k}))

for i in $(seq 1 ${num}); do start=$((${i}*${k}-${k})); end=$((${i}*${k})); for f in ${input_maf_files[@]}; do python seperate_maffile_by_index.py ${f} ${start} ${end} part${i}_${f}; done; python try_common_col_withDup_multi_alns_newversion.py --input part${i}_maf1 part${i}_maf2 part${i}_maf3 --output part${i}_common.maf --global_num 5 --global_start ${start} --global_end ${end};done

for f in ${input_maf_files[@]}; do python seperate_maffile_by_index.py ${f} $((${num}*${k})) ${length} part$((${num}+1))_${f}; done

python try_common_col_withDup_multi_alns_newversion.py --input part$((${num}+1))_maf1 part$((${num}+1))_maf2 part$((${num}+1))_maf3 --output part$((${num}+1))_common.maf --global_num 5 --global_start $((${num}*${k})) --global_end ${length}

for i in $((seq 1 $((${num}+1)))); do cat part${i}_common.maf >> common.maf
```

