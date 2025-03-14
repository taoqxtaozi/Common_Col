# Common-col extract common column alignments based on coordiantes of several MAF(multiple alignment format) files

Need python >= 3.6
Need to install packages: biopython interval3 collections multiprocessing
```
pip install biopython interval3 collections multiprocessing
```
#### Some prerequisites
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
4. 
