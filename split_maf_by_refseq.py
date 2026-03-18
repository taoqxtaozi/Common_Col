#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import sys
from typing import Dict, Tuple, Optional

from Bio import AlignIO
from Bio.AlignIO import MafIO


def is_reference(record_id: str, reference: str) -> bool:
    """
    Check whether a MAF record matches the reference.
    Supports two matching modes:
    1. Exact match, e.g. reference='GALGA.1', record_id='GALGA.1'
    2. Species-prefix match, e.g. reference='GALGA', record_id='GALGA.1'
    """
    return record_id == reference or record_id.split(".", 1)[0] == reference


def sanitize_filename(name: str) -> str:
    """
    Convert a sequence id into a filesystem-safe filename.
    """
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def get_ref_record(alignment, reference: Optional[str]):
    """
    Return the reference record from a block.
    If reference is provided, search for the matching record.
    If reference is None, use the first record in the block.
    """
    if len(alignment) == 0:
        return None

    if reference is None:
        return alignment[0]

    for record in alignment:
        if is_reference(record.id, reference):
            return record

    return None


def make_output_name(ref_id: str, strip_prefix: bool) -> str:
    """
    Generate the output filename stem for a reference sequence id.

    Examples:
    - ref_id='GALGA.1'
        strip_prefix=False -> 'GALGA.1'
        strip_prefix=True  -> '1'
    - ref_id='GALGA.chr1'
        strip_prefix=False -> 'GALGA.chr1'
        strip_prefix=True  -> 'chr1'
    """
    if strip_prefix and "." in ref_id:
        return ref_id.split(".", 1)[1]
    return ref_id


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Split a MAF file into one output file per reference sequence id. "
            "Blocks without the reference are discarded."
        )
    )
    parser.add_argument("--input", required=True, help="Input MAF file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument(
        "--reference",
        default=None,
        help=(
            "Reference name, e.g. GALGA or GALGA.1. "
            "If omitted, the first record of each block is used."
        ),
    )
    parser.add_argument(
        "--strip_prefix",
        action="store_true",
        help=(
            "If set, use only the part after the first dot in the reference id "
            "as the output filename stem. "
            "Example: GALGA.1 -> 1, GALGA.chr1 -> chr1"
        ),
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # key: output filepath
    # value: (file_handle, MafWriter)
    writers: Dict[str, Tuple[object, MafIO.MafWriter]] = {}

    kept_blocks = 0
    dropped_blocks = 0

    try:
        with open(args.input, "r") as infile:
            for alignment in AlignIO.parse(infile, "maf"):
                ref_record = get_ref_record(alignment, args.reference)

                # Drop blocks without the reference
                if ref_record is None:
                    dropped_blocks += 1
                    continue

                out_stem = make_output_name(ref_record.id, args.strip_prefix)
                out_name = sanitize_filename(out_stem) + ".maf"
                out_path = os.path.join(args.output_dir, out_name)

                # Lazily open each output file and write header once
                if out_path not in writers:
                    fh = open(out_path, "w")
                    writer = MafIO.MafWriter(fh)
                    writer.write_header()
                    writers[out_path] = (fh, writer)

                writers[out_path][1].write_alignment(alignment)
                kept_blocks += 1

    finally:
        for fh, _ in writers.values():
            fh.close()

    print(
        f"Done. kept_blocks={kept_blocks}, dropped_blocks={dropped_blocks}, "
        f"output_files={len(writers)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()