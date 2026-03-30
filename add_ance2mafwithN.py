# -*- coding: utf-8 -*-
from collections import defaultdict
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.AlignIO import MultipleSeqAlignment
import sys
from Bio.AlignIO.MafIO import MafWriter

# CARCRcommon3.maf CARCRcommon3_plusance.maf ance_CARCRcommon3
f1maf = sys.argv[1]
fnewmaf = sys.argv[2]
ance = sys.argv[3]

srcSize = 0
start = 0

for block in AlignIO.parse(f1maf, "maf"):
    srcSize += len(block[0].seq)

with open(f1maf, "r") as f, open(fnewmaf, "w") as out_f:
    writer = MafWriter(out_f)
    writer.write_header()
    for block in AlignIO.parse(f, "maf"):
        new_block=[]
        seq_length = len(block[0].seq)
        new_record = SeqRecord(
            Seq('N' * seq_length),
            id='.'.join([ance, "node"]),
            description="",
        )
        new_record.annotations = {
            "start": start,
            "size": seq_length,
            "strand": 1,
            "srcSize": srcSize
        }
        new_block.append(new_record)
        for record in block:
            new_block.append(record)
        writer.write_alignment(MultipleSeqAlignment(new_block))
        start += seq_length
