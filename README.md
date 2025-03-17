# Common-col
Extract common column alignments based on coordinates from multiple MAF (Multiple Alignment Format) files.

## Requirements
**Python version**: Python 3.6 or later

**Required packages**:
Install the necessary dependencies using:
```
pip install biopython interval3 collections multiprocessing
```
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
