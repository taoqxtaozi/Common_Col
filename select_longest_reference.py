#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
select_longest_reference.py

Select the genome with the greatest FASTA sequence length from a two-column
genome-path file.

Input format:
    LABEL    FASTA_PATH

Relative FASTA paths are resolved relative to the genome-path file directory.

Supported FASTA:
    - plain FASTA
    - gzip-compressed FASTA ending in .gz

By default stdout contains ONLY the selected label so the script can be used as:

    --refGenome "$(python select_longest_reference.py Y.genome.txt)"

Use --verbose to print genome lengths and the selected reference to stderr.

Genome length is the total number of non-whitespace sequence characters:
    - FASTA header lines beginning with ">" are ignored
    - N/n are counted
    - lowercase soft-masked bases are counted normally

If multiple genomes have exactly the same length, the lexicographically
smallest label is selected to make the result deterministic.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def read_genome_paths(genome_file: Path) -> List[Tuple[str, Path]]:
    genome_file = genome_file.expanduser().resolve()

    if not genome_file.is_file():
        raise FileNotFoundError(f"Genome-path file not found: {genome_file}")

    base_dir = genome_file.parent
    entries: List[Tuple[str, Path]] = []
    seen_labels: Dict[str, int] = {}

    with genome_file.open("r", encoding="utf-8") as fh:
        for lineno, raw_line in enumerate(fh, 1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split(maxsplit=1)
            if len(parts) != 2:
                raise ValueError(
                    f"Invalid line {lineno} in {genome_file}: "
                    f"expected 'LABEL FASTA_PATH', got {raw_line.rstrip()!r}"
                )

            label, path_text = parts

            if label in seen_labels:
                raise ValueError(
                    f"Duplicate label {label!r} in {genome_file}: "
                    f"first seen on line {seen_labels[label]}, repeated on line {lineno}"
                )
            seen_labels[label] = lineno

            fasta_path = Path(path_text).expanduser()
            if not fasta_path.is_absolute():
                fasta_path = (base_dir / fasta_path).resolve()
            else:
                fasta_path = fasta_path.resolve()

            entries.append((label, fasta_path))

    if not entries:
        raise ValueError(f"No genome entries found in {genome_file}")

    return entries


def open_fasta_binary(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def fasta_total_length(path: Path) -> int:
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    total = 0
    saw_header = False
    saw_sequence = False

    try:
        with open_fasta_binary(path) as fh:
            for raw_line in fh:
                if raw_line.startswith(b">"):
                    saw_header = True
                    continue

                # Count all non-whitespace bytes on sequence lines.
                seq_len = sum(
                    1 for b in raw_line
                    if b not in (9, 10, 11, 12, 13, 32)
                )

                if seq_len:
                    saw_sequence = True
                    total += seq_len

    except (OSError, EOFError) as exc:
        raise OSError(f"Failed to read FASTA file {path}: {exc}") from exc

    if not saw_header:
        raise ValueError(
            f"File does not appear to be FASTA (no '>' header found): {path}"
        )

    if not saw_sequence:
        raise ValueError(f"FASTA contains no sequence data: {path}")

    return total


def choose_longest_reference(
    entries: List[Tuple[str, Path]],
) -> Tuple[str, int, Dict[str, int]]:
    lengths: Dict[str, int] = {}

    for label, fasta_path in entries:
        lengths[label] = fasta_total_length(fasta_path)

    selected_label = min(
        lengths,
        key=lambda label: (-lengths[label], label),
    )

    return selected_label, lengths[selected_label], lengths


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Select the longest genome from a two-column LABEL FASTA_PATH file. "
            "stdout contains only the selected label."
        )
    )
    ap.add_argument(
        "genome_file",
        type=Path,
        help="Two-column file containing LABEL and FASTA_PATH",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help=(
            "Print all calculated genome lengths and the selected reference "
            "to stderr. stdout still contains only the selected label."
        ),
    )
    return ap


def main() -> None:
    ap = build_arg_parser()
    args = ap.parse_args()

    try:
        entries = read_genome_paths(args.genome_file)
        selected_label, selected_length, lengths = choose_longest_reference(entries)

        if args.verbose:
            print("# Genome lengths:", file=sys.stderr)
            width = max(len(label) for label, _ in entries)
            for label, _ in entries:
                print(f"{label:<{width}}  {lengths[label]}", file=sys.stderr)

            print(
                f"# Selected reference: {selected_label} ({selected_length} bp)",
                file=sys.stderr,
            )

        # stdout must contain only the selected label.
        print(selected_label)

    except Exception as exc:
        print(f"select_longest_reference.py: error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
