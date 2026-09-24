# -*- coding: utf-8 -*-
from Bio import AlignIO
from Bio import SeqIO
from Bio.AlignIO import MafIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import argparse

#usage python .py interior2.info.txt CARCRcommon4_plusAnceInterior9.maf interior2
parser = argparse.ArgumentParser()
parser.add_argument("fasta_info_file")
parser.add_argument("my_maf")
parser.add_argument("id_to_find")
parser.add_argument("--ifAnc", type=int, choices=(0, 1), default=0)
parser.add_argument("--genome", help="FASTA file for the ancestor when --ifAnc 1")
args = parser.parse_args()

if args.ifAnc == 1 and args.genome is None:
    parser.error("--genome is required when --ifAnc 1")
if args.ifAnc == 0 and args.genome is not None:
    parser.error("--genome can only be used when --ifAnc 1")

fasta_info_file = args.fasta_info_file
my_maf = args.my_maf
id_to_find = args.id_to_find

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

ancestor_sequences = {}
if args.ifAnc == 1:
    try:
        with open(args.genome, "r") as genome_file:
            for fasta_record in SeqIO.parse(genome_file, "fasta"):
                ancestor_id = "".join([id_to_find, ".", fasta_record.id])
                if ancestor_id not in missing_ids:
                    continue
                if ancestor_id in ancestor_sequences:
                    parser.error(f"duplicate sequence '{fasta_record.id}' in {args.genome}")

                sequence = str(fasta_record.seq)
                expected_length = fasta_info[ancestor_id]
                if len(sequence) != expected_length:
                    parser.error(
                        f"length mismatch for '{fasta_record.id}': "
                        f"FASTA has {len(sequence)} bases, .info.txt has {expected_length}"
                    )
                if "-" in sequence:
                    parser.error(f"ancestor FASTA sequence '{fasta_record.id}' contains gaps")
                ancestor_sequences[ancestor_id] = sequence
    except OSError as error:
        parser.error(f"cannot read ancestor FASTA '{args.genome}': {error}")

    absent_ids = missing_ids - set(ancestor_sequences.keys())
    if absent_ids:
        parser.error(
            "ancestor FASTA is missing sequence(s): " + ", ".join(sorted(absent_ids))
        )

with open(my_maf, "a") as f_out:
    for missing_id in sorted(missing_ids):
        length = fasta_info[missing_id]

        if args.ifAnc == 1:
            sequence = ancestor_sequences[missing_id]
            size = length
        else:
            sequence = "-"
            size = 0

        record = SeqRecord(
            Seq(sequence),
            id=missing_id,
            description="",
            annotations={
                "start": 0,
                "size": size,
                "strand": 1,
                "srcSize": length
            }
        )

        block = MultipleSeqAlignment([record])
        AlignIO.write(block, f_out, "maf")
