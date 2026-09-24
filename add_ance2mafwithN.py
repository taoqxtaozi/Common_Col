# -*- coding: utf-8 -*-

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.AlignIO import MultipleSeqAlignment
from Bio.AlignIO.MafIO import MafWriter
from contextlib import nullcontext
import os
import sys

# Usage:
#   python add_ance2mafwithN.py input.maf output.maf AncName
#   python add_ance2mafwithN.py input.maf output.maf AncName SeqName
#   python add_ance2mafwithN.py input.maf output.maf AncName SeqName RefGenome.fa
#   python add_ance2mafwithN.py input.maf output.maf AncName SeqName IngroupAncestor
#   python add_ance2mafwithN.py input.maf output.maf AncName SeqName IngroupAncestor --inference-bed multirow.bed
#
# Examples:
#   python add_ance2mafwithN.py consensus.maf consensus_A1.maf A1
#       -> ancestor src name: A1.node
#
#   python add_ance2mafwithN.py consensus.maf consensus_A1.maf A1 chr1
#       -> ancestor src name: A1.chr1
#
#   python add_ance2mafwithN.py consensus.maf consensus_A1.maf A1 chr1 path/ref.fa
#       -> ancestor src name: A1.chr1
#       -> ancestor srcSize is the length of chr1 in ref.fa
#       -> regions not covered by consensus blocks are filled with chr1 reference sequence
#
#   python add_ance2mafwithN.py input.maf output.maf Anc0 chr1 Anc_Input
#       -> blocks containing only Anc_Input.chr1 copy that sequence into Anc0.chr1
#       -> blocks also containing other rows initialize Anc0.chr1 with N
#       -> the ingroup ancestor must cover the chromosome without coordinate gaps
#       -> optional BED limits ancestorsML to blocks with other rows

if len(sys.argv) not in (4, 5, 6, 8) or (len(sys.argv) == 8 and sys.argv[6] != "--inference-bed"):
    sys.stderr.write(
        "Usage:\n"
        "  python add_ance2mafwithN.py input.maf output.maf AncName\n"
        "  python add_ance2mafwithN.py input.maf output.maf AncName SeqName\n"
        "  python add_ance2mafwithN.py input.maf output.maf AncName SeqName RefGenome.fa\n"
        "  python add_ance2mafwithN.py input.maf output.maf AncName SeqName IngroupAncestor\n"
        "  python add_ance2mafwithN.py input.maf output.maf AncName SeqName IngroupAncestor --inference-bed multirow.bed\n"
    )
    sys.exit(1)

f1maf = sys.argv[1]
fnewmaf = sys.argv[2]
ance = sys.argv[3]

if len(sys.argv) >= 5:
    seq_name = sys.argv[4]
else:
    seq_name = "node"

ref_fasta = None
ingroup_ancestor = None
inference_bed = sys.argv[7] if len(sys.argv) == 8 else None
if len(sys.argv) >= 6:
    last_arg = sys.argv[5]
    # The existing five-argument call passes a FASTA path.  A bare ancestor
    # name selects the new mode; recognize existing extensionless FASTA paths.
    if (os.path.isfile(last_arg) or os.path.sep in last_arg or
            last_arg.lower().endswith((".fa", ".fasta", ".fna"))):
        ref_fasta = last_arg
    else:
        ingroup_ancestor = last_arg

fill_missing_with_ref = ref_fasta is not None

ancestor_src = f"{ance}.{seq_name}"
ingroup_src = f"{ingroup_ancestor}.{seq_name}" if ingroup_ancestor else None

if inference_bed is not None and ingroup_src is None:
    sys.stderr.write("ERROR: --inference-bed requires an ingroup ancestor name\n")
    sys.exit(1)

if ingroup_src == ancestor_src:
    sys.stderr.write("ERROR: the new ancestor and ingroup ancestor must have different names\n")
    sys.exit(1)


def get_maf_seq_name(src):
    if "." in src:
        return src.split(".", 1)[1]
    return src


def find_reference_record(block, seq_name):
    # In the current pipeline, the reference row should be the first row
    # of each chromosome-specific consensus MAF block.
    if len(block) > 0 and get_maf_seq_name(block[0].id) == seq_name:
        return block[0]

    matched_records = []
    for record in block:
        if get_maf_seq_name(record.id) == seq_name:
            matched_records.append(record)

    if len(matched_records) == 1:
        return matched_records[0]

    if len(matched_records) == 0:
        sys.stderr.write(
            f"ERROR: no reference-like row with sequence name '{seq_name}' "
            f"was found in one MAF block of {f1maf}\n"
        )
    else:
        sys.stderr.write(
            f"ERROR: multiple rows with sequence name '{seq_name}' were found "
            f"in one MAF block of {f1maf}; cannot determine the reference row\n"
        )

    sys.exit(1)


def read_reference_sequence(ref_fasta, seq_name):
    for record in SeqIO.parse(ref_fasta, "fasta"):
        if record.id == seq_name or record.name == seq_name or record.description.split()[0] == seq_name:
            return str(record.seq).upper()

    sys.stderr.write(
        f"ERROR: sequence '{seq_name}' was not found in reference FASTA: {ref_fasta}\n"
    )
    sys.exit(1)


def make_ancestor_record(seq, start, size, src_size):
    new_record = SeqRecord(
        Seq(seq),
        id=ancestor_src,
        description="",
    )

    new_record.annotations = {
        "start": start,
        "size": size,
        "strand": 1,
        "srcSize": src_size,
    }

    return new_record


def write_ancestor_only_block(writer, seq, start, src_size):
    if len(seq) == 0:
        return

    new_record = make_ancestor_record(
        seq=seq,
        start=start,
        size=len(seq),
        src_size=src_size,
    )

    writer.write_alignment(MultipleSeqAlignment([new_record]))


def get_ingroup_record(block):
    records = [record for record in block if record.id == ingroup_src]
    if len(records) != 1:
        sys.stderr.write(
            f"ERROR: expected exactly one '{ingroup_src}' row per non-empty "
            f"block in {f1maf}; found {len(records)}\n"
        )
        sys.exit(1)
    if any(record.id == ancestor_src for record in block):
        sys.stderr.write(f"ERROR: '{ancestor_src}' already exists in {f1maf}\n")
        sys.exit(1)
    return records[0]


def get_new_ancestor_sequence(block, seq_length):
    if ingroup_src is not None and len(block) == 1:
        sequence = str(block[0].seq)
        return sequence, len(sequence) - sequence.count("-")
    return "N" * seq_length, seq_length


ref_seq = None
if fill_missing_with_ref:
    ref_seq = read_reference_sequence(ref_fasta, seq_name)
    srcSize = len(ref_seq)
else:
    srcSize = 0

start = 0
non_empty_blocks = 0
empty_blocks = 0
ingroup_end = 0
ingroup_src_size = None

# First pass: calculate total ancestor srcSize.
# Skip empty MAF blocks.
for block in AlignIO.parse(f1maf, "maf"):
    if len(block) == 0:
        empty_blocks += 1
        continue

    seq_length = len(block[0].seq)
    if seq_length == 0:
        empty_blocks += 1
        continue

    if ingroup_src is not None:
        ingroup_record = get_ingroup_record(block)
        annotations = ingroup_record.annotations
        ingroup_start = int(annotations["start"])
        ingroup_size = int(annotations["size"])
        record_src_size = int(annotations["srcSize"])
        if ingroup_src_size is None:
            ingroup_src_size = record_src_size
        actual_size = len(ingroup_record.seq) - str(ingroup_record.seq).count("-")
        if (annotations["strand"] != 1 or ingroup_start != ingroup_end or
                ingroup_size != actual_size or record_src_size != ingroup_src_size):
            sys.stderr.write(
                f"ERROR: '{ingroup_src}' does not continuously cover {seq_name} "
                f"in {f1maf}: expected_start={ingroup_end}, "
                f"row_start={ingroup_start}, annotation_size={ingroup_size}, "
                f"non_gap_size={actual_size}\n"
            )
            sys.exit(1)
        ingroup_end += ingroup_size

    if fill_missing_with_ref:
        ref_record = find_reference_record(block, seq_name)
        ref_start = int(ref_record.annotations["start"])
        ref_size = int(ref_record.annotations["size"])
        ref_end = ref_start + ref_size

        if ref_start < start:
            sys.stderr.write(
                f"ERROR: MAF blocks are not sorted by reference coordinate, "
                f"or they overlap on {seq_name}: previous_end={start}, "
                f"current_start={ref_start}\n"
            )
            sys.exit(1)

        if ref_end > srcSize:
            sys.stderr.write(
                f"ERROR: MAF block exceeds reference sequence length for {seq_name}: "
                f"block_end={ref_end}, reference_length={srcSize}\n"
            )
            sys.exit(1)

        start = ref_end
    else:
        _, ancestor_size = get_new_ancestor_sequence(block, seq_length)
        srcSize += ancestor_size
        start += ancestor_size

    non_empty_blocks += 1

if non_empty_blocks == 0:
    sys.stderr.write(f"ERROR: no non-empty alignment blocks found in {f1maf}\n")
    sys.exit(1)

if ingroup_src is not None and ingroup_end != ingroup_src_size:
    sys.stderr.write(
        f"ERROR: '{ingroup_src}' ends at {ingroup_end}, but its srcSize is "
        f"{ingroup_src_size}; ingroup blocks do not cover the whole chromosome\n"
    )
    sys.exit(1)

if empty_blocks > 0:
    sys.stderr.write(
        f"WARNING: skipped {empty_blocks} empty MAF blocks in {f1maf}\n"
    )

# Second pass: add ancestor row to each non-empty block.
with open(f1maf, "r") as f, open(fnewmaf, "w") as out_f, (
    open(inference_bed, "w") if inference_bed else nullcontext()
) as bed_out:
    writer = MafWriter(out_f)
    writer.write_header()

    if fill_missing_with_ref:
        prev_ref_end = 0

        for block in AlignIO.parse(f, "maf"):
            if len(block) == 0:
                continue

            seq_length = len(block[0].seq)
            if seq_length == 0:
                continue

            ref_record = find_reference_record(block, seq_name)
            ref_start = int(ref_record.annotations["start"])
            ref_size = int(ref_record.annotations["size"])
            ref_end = ref_start + ref_size

            if ref_start < prev_ref_end:
                sys.stderr.write(
                    f"ERROR: MAF blocks are not sorted by reference coordinate, "
                    f"or they overlap on {seq_name}: previous_end={prev_ref_end}, "
                    f"current_start={ref_start}\n"
                )
                sys.exit(1)

            # Fill the reference interval before this consensus block.
            if ref_start > prev_ref_end:
                fill_seq = ref_seq[prev_ref_end:ref_start]
                write_ancestor_only_block(
                    writer=writer,
                    seq=fill_seq,
                    start=prev_ref_end,
                    src_size=srcSize,
                )

            ref_block_seq = str(ref_record.seq)
            ancestor_block_seq = "".join(
                "-" if base == "-" else "N"
                for base in ref_block_seq
            )
            ancestor_block_size = sum(1 for base in ancestor_block_seq if base != "-")

            if ancestor_block_size != ref_size:
                sys.stderr.write(
                    f"ERROR: reference row size does not match non-gap length "
                    f"in block on {seq_name}: annotation_size={ref_size}, "
                    f"non_gap_length={ancestor_block_size}, start={ref_start}\n"
                )
                sys.exit(1)

            new_block = []

            new_record = make_ancestor_record(
                seq=ancestor_block_seq,
                start=ref_start,
                size=ancestor_block_size,
                src_size=srcSize,
            )

            new_block.append(new_record)

            for record in block:
                new_block.append(record)

            writer.write_alignment(MultipleSeqAlignment(new_block))
            prev_ref_end = ref_end

        # Fill the reference interval after the last consensus block.
        if prev_ref_end < srcSize:
            fill_seq = ref_seq[prev_ref_end:srcSize]
            write_ancestor_only_block(
                writer=writer,
                seq=fill_seq,
                start=prev_ref_end,
                src_size=srcSize,
            )

        if prev_ref_end > srcSize:
            sys.stderr.write(
                f"ERROR: internal length mismatch: previous_end={prev_ref_end}, "
                f"srcSize={srcSize}\n"
            )
            sys.exit(1)

    else:
        start = 0
        inference_start = None
        inference_end = None

        for block in AlignIO.parse(f, "maf"):
            if len(block) == 0:
                continue

            seq_length = len(block[0].seq)
            if seq_length == 0:
                continue

            ancestor_sequence, ancestor_size = get_new_ancestor_sequence(
                block, seq_length
            )
            if bed_out is not None:
                if len(block) > 1:
                    if inference_start is None:
                        inference_start = start
                    inference_end = start + ancestor_size
                elif inference_start is not None:
                    bed_out.write(f"{seq_name}\t{inference_start}\t{inference_end}\n")
                    inference_start = None
                    inference_end = None
            new_block = []
            new_record = make_ancestor_record(
                seq=ancestor_sequence,
                start=start,
                size=ancestor_size,
                src_size=srcSize,
            )

            new_block.append(new_record)

            for record in block:
                new_block.append(record)

            writer.write_alignment(MultipleSeqAlignment(new_block))
            start += ancestor_size

        if bed_out is not None and inference_start is not None:
            bed_out.write(f"{seq_name}\t{inference_start}\t{inference_end}\n")

        if start != srcSize:
            sys.stderr.write(
                f"ERROR: internal length mismatch: start={start}, srcSize={srcSize}\n"
            )
            sys.exit(1)
