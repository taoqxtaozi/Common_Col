# -*- coding: utf-8 -*-
"""
Replace an ancestor row in a MAF file with the sequence from a one-sequence FASTA file.

Usage:
    python replace_ancestor_seq_in_maf.py new.maf Anc Anc.fa newest.maf
    python replace_ancestor_seq_in_maf.py new.maf Anc Anc.fa newest.maf SeqName

Arguments:
    new.maf      Input MAF file containing rows named Anc.node or Anc.SeqName in some blocks
    Anc          Ancestor name.
                 With no SeqName argument, the script will replace rows named Anc.node.
                 With SeqName, the script will replace rows named Anc.SeqName.
    Anc.fa       FASTA file containing one inferred ancestor sequence
    newest.maf   Output MAF file
    SeqName      Optional sequence/chromosome name, e.g. chr20.
                 If provided, target row is Anc.SeqName instead of Anc.node.

Behavior:
    Blocks containing the target ancestor row:
        replace the target sequence using the inferred ancestor FASTA.

    Blocks not containing the target ancestor row:
        write the block unchanged and do not consume ancestor FASTA sequence.

    Blocks containing more than one target ancestor row:
        stop with an error.
"""

import sys
from pathlib import Path
from typing import List, Optional, Tuple

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO.MafIO import MafWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def die(message: str) -> None:
    sys.stderr.write(f"ERROR: {message}\n")
    sys.exit(1)


def warn(message: str) -> None:
    sys.stderr.write(f"WARNING: {message}\n")


def read_single_fasta_sequence(fasta_path: str) -> str:
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) == 0:
        die(f"No FASTA record found in {fasta_path}")
    if len(records) > 1:
        die(f"{fasta_path} contains more than one FASTA record; expected exactly one.")
    return str(records[0].seq).replace("\n", "").replace("\r", "").upper()


def find_target_record_index(block: MultipleSeqAlignment, target_id: str) -> Optional[int]:
    hits: List[int] = [i for i, rec in enumerate(block) if rec.id == target_id]

    if len(hits) == 0:
        return None

    if len(hits) > 1:
        die(f"Found more than one target row '{target_id}' in one MAF block.")

    return hits[0]


def get_record_annotation_int(record: SeqRecord, key: str, target_id: str) -> int:
    if key not in record.annotations:
        die(f"Target row '{target_id}' is missing required MAF annotation '{key}'.")

    try:
        return int(record.annotations[key])
    except Exception:
        die(
            f"Target row '{target_id}' has non-integer MAF annotation "
            f"{key}={record.annotations[key]!r}."
        )


def non_gap_length(seq: str) -> int:
    return sum(1 for base in seq if base != "-")


def build_replacement_with_gap_pattern(
    ancestor_seq: str,
    record: SeqRecord,
    target_id: str,
) -> str:
    row_seq = str(record.seq)
    start = get_record_annotation_int(record, "start", target_id)
    size = get_record_annotation_int(record, "size", target_id)
    src_size = get_record_annotation_int(record, "srcSize", target_id)

    if start < 0 or size < 0 or src_size < 0:
        die(
            f"Target row '{target_id}' has invalid negative coordinate: "
            f"start={start}, size={size}, srcSize={src_size}."
        )

    end = start + size
    if end > src_size:
        die(
            f"Target row '{target_id}' exceeds its srcSize: "
            f"start={start}, size={size}, srcSize={src_size}."
        )

    if src_size != len(ancestor_seq):
        die(
            f"Length mismatch for {target_id}: FASTA length = {len(ancestor_seq)}, "
            f"but MAF srcSize = {src_size}."
        )

    row_non_gap_len = non_gap_length(row_seq)
    if row_non_gap_len != size:
        die(
            f"MAF row size does not match non-gap sequence length for {target_id}: "
            f"annotation size = {size}, non-gap length = {row_non_gap_len}, "
            f"start = {start}."
        )

    replacement_bases = ancestor_seq[start:end]

    cursor = 0
    out_chars = []
    for base in row_seq:
        if base == "-":
            out_chars.append("-")
        else:
            out_chars.append(replacement_bases[cursor])
            cursor += 1

    if cursor != size:
        die(
            f"Internal replacement cursor error for {target_id}: "
            f"used {cursor} bases, expected {size}."
        )

    return "".join(out_chars)


def get_target_row_summary(maf_path: str, target_id: str) -> Tuple[int, int, int]:
    n_blocks = 0
    n_target_blocks = 0
    required_non_gap_bases = 0

    with open(maf_path, "r") as handle:
        for block in AlignIO.parse(handle, "maf"):
            n_blocks += 1
            idx = find_target_record_index(block, target_id)

            if idx is None:
                continue

            record = block[idx]
            size = get_record_annotation_int(record, "size", target_id)
            row_non_gap_len = non_gap_length(str(record.seq))

            if row_non_gap_len != size:
                start = get_record_annotation_int(record, "start", target_id)
                die(
                    f"MAF row size does not match non-gap sequence length for {target_id}: "
                    f"annotation size = {size}, non-gap length = {row_non_gap_len}, "
                    f"start = {start}."
                )

            required_non_gap_bases += size
            n_target_blocks += 1

    if n_blocks == 0:
        die(f"No alignment blocks found in {maf_path}")

    if n_target_blocks == 0:
        die(f"No block containing target row '{target_id}' was found in {maf_path}")

    return required_non_gap_bases, n_blocks, n_target_blocks


def copy_record_with_new_sequence(record: SeqRecord, new_seq: str) -> SeqRecord:
    new_record = SeqRecord(
        Seq(new_seq),
        id=record.id,
        name=record.name,
        description=record.description,
    )

    # Preserve MAF annotations such as start, size, strand, and srcSize.
    new_record.annotations = dict(record.annotations)

    return new_record


def replace_ancestor_sequence_in_maf(
    input_maf: str,
    ancestor_name: str,
    ancestor_fasta: str,
    output_maf: str,
    seq_name: str = "node",
) -> None:
    target_id = f"{ancestor_name}.{seq_name}"
    ancestor_seq = read_single_fasta_sequence(ancestor_fasta)

    required_non_gap_bases, n_blocks, n_target_blocks = get_target_row_summary(input_maf, target_id)

    if required_non_gap_bases > len(ancestor_seq):
        die(
            f"Length mismatch for {target_id}: target MAF rows require at least "
            f"{required_non_gap_bases} non-gap bases, but FASTA length is "
            f"{len(ancestor_seq)}."
        )

    n_missing_blocks = n_blocks - n_target_blocks
    if n_missing_blocks > 0:
        warn(
            f"{n_missing_blocks} of {n_blocks} MAF blocks do not contain '{target_id}' "
            f"and will be written unchanged."
        )

    out_path = Path(output_maf)
    if out_path.parent and str(out_path.parent) != ".":
        out_path.parent.mkdir(parents=True, exist_ok=True)

    replaced_blocks = 0
    skipped_blocks = 0

    with open(input_maf, "r") as in_handle, open(output_maf, "w") as out_handle:
        writer = MafWriter(out_handle)
        writer.write_header()

        for block in AlignIO.parse(in_handle, "maf"):
            idx = find_target_record_index(block, target_id)

            if idx is None:
                # No ancestor row in this block; keep the block unchanged.
                writer.write_alignment(block)
                skipped_blocks += 1
                continue

            replacement = build_replacement_with_gap_pattern(
                ancestor_seq=ancestor_seq,
                record=block[idx],
                target_id=target_id,
            )

            new_records = []
            for i, record in enumerate(block):
                if i == idx:
                    new_records.append(copy_record_with_new_sequence(record, replacement))
                else:
                    new_records.append(record)

            writer.write_alignment(MultipleSeqAlignment(new_records))
            replaced_blocks += 1

    sys.stderr.write(
        f"Done. Replaced '{target_id}' in {replaced_blocks} blocks; "
        f"skipped {skipped_blocks} blocks without '{target_id}'.\n"
    )


def main() -> None:
    if len(sys.argv) not in (5, 6):
        sys.stderr.write(
            "Usage:\n"
            "    python replace_ancestor_seq_in_maf.py new.maf Anc Anc.fa newest.maf\n"
            "    python replace_ancestor_seq_in_maf.py new.maf Anc Anc.fa newest.maf SeqName\n"
        )
        sys.exit(1)

    input_maf = sys.argv[1]
    ancestor_name = sys.argv[2]
    ancestor_fasta = sys.argv[3]
    output_maf = sys.argv[4]

    if len(sys.argv) == 6:
        seq_name = sys.argv[5]
    else:
        seq_name = "node"

    replace_ancestor_sequence_in_maf(
        input_maf=input_maf,
        ancestor_name=ancestor_name,
        ancestor_fasta=ancestor_fasta,
        output_maf=output_maf,
        seq_name=seq_name,
    )


if __name__ == "__main__":
    main()
