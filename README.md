# Common-col
Extract common column alignments based on coordinates from multiple MAF (Multiple Alignment Format) files.

## Requirements
**Python version**: Python 3.6 or later

**Required packages**:
Install the necessary dependencies using:
```
pip install biopython interval3 collections multiprocessing
```

## Build alignment with Progressive alignment
### Step 1: Decide how many guide-tree-based alignments used to extract the consensus
We build our own script ```generate_random_guidetrees_2models_2modes_usethis_fast2_finalver.py``` to generate fully resolved binary guide trees. In the default mode, we set RF-distance at least 1 and triplet distance at least 2/3, and the species evolve along "yule" mode. In the parameter setting, users can set RF-distance with ```--rf-threshold```, set triplet distance with ```--t-threshold```, and evolution model can be set with ```--model``` choosing from 'yule' or 'uniform'. And need to set the outgroup with ```--outgroup``` and input the ingroup taxa with ```--taxa```.

For the number of guide tree used, in default we use four, for making a balance between alignment quality and resource usage. And with enough computation resources, we suggest that, if ```n```, the number of the ingroup taxa, is no more than ten, we suggest to use ```n``` guide trees, and if it is more than 10, we suggest to use ten. And here the ingroup taxa, we mean the number of unresolved polytomies, that is, if there are groups that the users wish to fix in all guide trees, they are seen as one "taxon".

The example of how to use is: 
If the ingroup taxa are A, B, C, D, E, and F; the outgroup is G
If use the default mode：
```
python generate_random_guidetrees_2models_2modes_usethis_fast2_finalver.py --taxa A B C D E F --outgroup G
```
If you want to relax the constraints among the trees, like make RF distance no less than 0.8, and make triplet distance no less than 0.6; and let taxa evolve along uniform model, also you want to fix (A,B), and (C,D) as two cherries in all guide trees, and want to get 6 guide trees, you can set:
```
python generate_random_guidetrees_2models_2modes_usethis_fast2_finalver.py --taxa "(A,B)" "(C,D)" E F --outgroup G --num_trees 6 --model 'uniform' --rf-threshold 0.8 --t-threshold 0.6
```
Here in default we use fully resolved binary trees as guide trees, if you decide to use ```X``` guide trees, but you want to include one star tree, when generating guide trees with our script, you can set ```--num_trees``` as ```X-1```

### Step 2: Align taxa with Progressive Cactus
Like now we get N guide trees, we need to add that into the alignment set files, ```aln1.txt```, ```aln2.txt```, ```aln3.txt```, ..., ```alnN.txt```, respectively, which is needed in Progressive cactus.


## Prerequisites
If the maf alignments are generated using [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus), you can just see **step2**.

### Step 1: Ensure a consistent reference in all MAF files
Each block in all input MAF files must have a **reference sequence**, which should be positioned as the **first row** in every block.
If you need to adjust the order of taxa(e.g., A,B,C,D with A as the reference) in blocks, use the following command:   
```
conda install bioconda::phast
maf_parse -O A,B,C,D input1.maf > input1_order.maf
```
### Step 2: Ensure that all MAF files contain alignments for only one reference chromosome
All input MAF files (e.g., ```maf1```, ```maf2```, ```maf3```) should contain alignments for only **one chromosome** of the reference.
For example, if the reference sequence in ```maf1```, ```maf2```, and ```maf3``` is ```A.chr1```, all alignment blocks should start with ```A.chr1```.
It is **not allowed** to mix multiple reference chromosomes (e.g., both ```A.chr1``` and ```A.chr2```).
To separate MAF files by chromosome, use the provided script:
```
python separate_maf.py maf output_directory
```
### Step 3: Adjust block orientation and order
To ensure that all blocks are sorted by reference coordinates and that the reference strand is positive, you can use [mafTools](https://github.com/dentearl/mafTools)
```
conda create -n py2 python=2.7
conda activate py2
conda install genomedk::maftools
mafStrander --maf alignment_Achr1.maf --seq A.chr1 --strand + > alignment_Achr1positive.maf
mafSorter --maf alignment_Achr1positive --seq A.chr1 > alignment_Achr1positive_sorted.maf
conda deactivate
```
## Usage

Once the input MAF files are properly formatted, you can run the script to extract common alignments:
```
length=`grep -m1 "^s " maf1 | awk '{print $6}'`
python try_common_col_withDup_multi_alns_newversion.py \
   --input maf1 maf2 maf3 ... \
   --output common.maf --global_num 5 \
   --global_start 0 --global_end ${length}
```
### Handling large reference chromosomes
If the reference chromosome ```length``` is to long (e.g., **> 10,925,261 bp**), it is recommended to **split** the input MAF files into smaller parts.
You can divide the input equally using a custom **partition size (```k```)**, which you can set manually.
Run the following command to split and process the MAF files:
```
bash split_and_process --input maf1 maf2 maf3... -k 2185052 --global_num 5
```
The reults will be saved in ```common.maf```

## Notes
- Ensure that all MAF files follow the correct formatting before running the script.
- If you encounter performance issues, consider **parallelizing** the processing or **optimizing memory usage**.
- For more details, refer to the documentation of the tools mentioned above.
