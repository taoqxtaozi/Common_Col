#!/usr/bin/env /home/taoqixin/miniconda3/bin/python
# usage: python xxx.py xxx.maf num1 num2 ./folder
from Bio import AlignIO
import sys

f1 = sys.argv[1]
index = sys.argv[2]
k1 = sys.argv[3]
k2 = sys.argv[4]
out_fold = sys.argv[5]

#f1 = "C://F//Download//software//cactus-bin-v2.1.1//GALGA_CHLUN_COLLI_MESUN//aln3//GALGA.1.maf"
lst1 = AlignIO.parse(f1, "maf")

#k1 = 10925261
#k2 = 10925261 * 2

fout = out_fold + "/" + index + "_"  + k1 + "_" + k2 + ".maf"
lst = []
k1_1 = int(k1)
k2_1 = int(k2)
for aln in lst1:
    if aln[0].annotations["start"] >= k1_1:
        if aln[0].annotations["start"] + aln[0].annotations["size"] <= k2_1:
            lst.append(aln)
        elif aln[0].annotations["start"] <= k2_1:
            lst.append(aln)
        else:
            break
    elif aln[0].annotations["start"] < k1_1:
        if aln[0].annotations["start"] + aln[0].annotations["size"] >= k1_1:
            lst.append(aln)
        else:
            continue
AlignIO.write(lst, fout, "maf")
    
