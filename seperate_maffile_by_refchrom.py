#!/usr/bin/env python
# how_to_use python .py orignal.maf outputfolder

import sys
from Bio import AlignIO

f1 = sys.argv[1]
pathout = sys.argv[2]

galga_chrom = []
for multiple_alignmnet in AlignIO.parse(f1,"maf"):
    if len(galga_chrom) == 0:
        galga_chrom.append(multiple_alignmnet)
    elif multiple_alignmnet[0].id == galga_chrom[-1][0].id:
        galga_chrom.append(multiple_alignmnet)
    else:
        AlignIO.write(galga_chrom, pathout + '/' + galga_chrom[-1][0].id+ ".maf", "maf")
        galga_chrom = []
        galga_chrom.append(multiple_alignmnet)
AlignIO.write(galga_chrom, pathout + '/' + galga_chrom[-1][0].id+ ".maf", "maf")

