# -*- coding: utf-8 -*-
from Bio import AlignIO
from Bio.AlignIO import MafIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import sys 

#usage python .py interior2.info.txt CARCRcommon4_plusAnceInterior9.maf interior2
fasta_info_file = sys.argv[1]
my_maf = sys.argv[2]
id_to_find = sys.argv[3]

fasta_info = {}
with open(fasta_info_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            seq_id, length = parts
            seq_id = "".join([id_to_find, ".", seq_id])
            fasta_info[seq_id] = int(length)


seen_ids = set()
with open(my_maf, "r") as f:
    for block in AlignIO.parse(f, "maf"):
        for record in block:
            if record.id.startswith(id_to_find):
                seen_ids.add(record.id)

missing_ids = set(fasta_info.keys()) - seen_ids

with open(my_maf, "a") as f_out:
    for missing_id in sorted(missing_ids):
        length = fasta_info[missing_id]

        record = SeqRecord(
            Seq("-"),
            id=missing_id,
            description="",
            annotations={
                "start": 0,
                "size": 0,
                "strand": 1,
                "srcSize": length
            }
        )

        block = MultipleSeqAlignment([record])
        AlignIO.write(block, f_out, "maf")
