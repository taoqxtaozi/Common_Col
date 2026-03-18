#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO import MafIO


def is_reference(record_id: str, reference: str) -> bool:
    """
    Check whether a MAF record matches the reference.
    Supports two matching modes:
    1. Exact match, e.g. reference='GALGA.1', record_id='GALGA.1'
    2. Species-prefix match, e.g. reference='GALGA', record_id='GALGA.1'
    """
    return record_id == reference or record_id.split(".", 1)[0] == reference


def reorder_block(alignment: MultipleSeqAlignment, reference: str):
    """
    If the reference is present in the block, return a new MultipleSeqAlignment:
    - the reference record is moved to the first position
    - all other records keep their original relative order
    If the reference is absent, return None.
    """
    ref_record = None
    other_records = []

    for record in alignment:
        if ref_record is None and is_reference(record.id, reference):
            ref_record = record
        else:
            other_records.append(record)

    if ref_record is None:
        return None

    new_alignment = MultipleSeqAlignment(
        [ref_record] + other_records,
        annotations=getattr(alignment, "annotations", {}),
        column_annotations=getattr(alignment, "column_annotations", {}),
    )
    return new_alignment


def main():
    parser = argparse.ArgumentParser(
        description="Reorder each MAF block so that the reference sequence is first; discard blocks without the reference."
    )
    parser.add_argument("--input", required=True, help="Input MAF file")
    parser.add_argument("--output", required=True, help="Output MAF file")
    parser.add_argument("--reference", required=True, help="Reference name, e.g. GALGA or GALGA.1")
    args = parser.parse_args()

    kept_blocks = 0
    dropped_blocks = 0

    with open(args.input, "r") as infile, open(args.output, "w") as outfile:
        writer = MafIO.MafWriter(outfile)
        writer.write_header()

        for alignment in AlignIO.parse(infile, "maf"):
            new_alignment = reorder_block(alignment, args.reference)

            if new_alignment is None:
                dropped_blocks += 1
                continue

            writer.write_alignment(new_alignment)
            kept_blocks += 1

    print(
        f"Done. kept_blocks={kept_blocks}, dropped_blocks={dropped_blocks}",
        file=sys.stderr
    )


if __name__ == "__main__":
    main()