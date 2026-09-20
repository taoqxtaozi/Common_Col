# -*- coding: utf-8 -*-

"""
Split one chromosome-level MAF into fixed-length coordinate chunks
by scanning the input MAF only once.

Usage:
    python seperate_maffile_by_index.py input.maf sep_length output_folder

Example:
    python seperate_maffile_by_index.py GALGA.1.maf 218505 ./GALGA.1

If the reference chromosome srcSize is 200994015, outputs will be:

    1_0_218505.maf
    2_218505_437010.maf
    ...
    920_200806095_200994015.maf

The first s-line of every alignment block is assumed to be the reference.

IMPORTANT:
The overlap rule intentionally reproduces the behavior of the old script:
a block touching or crossing a chunk boundary is written to both
adjacent chunks.
"""

import os
import sys
from collections import OrderedDict


def usage():
    print(
        "Usage: python seperate_maffile_by_index.py "
        "input.maf sep_length output_folder",
        file=sys.stderr,
    )
    sys.exit(1)


if len(sys.argv) != 4:
    usage()


input_maf = sys.argv[1]

try:
    sep_length = int(sys.argv[2])
except ValueError:
    print("Error: sep_length must be an integer.", file=sys.stderr)
    sys.exit(1)

output_folder = sys.argv[3]

if sep_length <= 0:
    print("Error: sep_length must be > 0.", file=sys.stderr)
    sys.exit(1)

if not os.path.isfile(input_maf):
    print("Error: input MAF not found: {}".format(input_maf), file=sys.stderr)
    sys.exit(1)

os.makedirs(output_folder, exist_ok=True)


def is_alignment_start(line):
    """
    Return True if this line starts a new MAF alignment block.
    """
    stripped = line.lstrip()
    return stripped == "a\n" or stripped == "a" or stripped.startswith("a ")


def iter_blocks(handle, first_a_line):
    """
    Yield one complete MAF alignment block at a time.

    first_a_line is the first 'a' line already read from the file.
    """
    block = [first_a_line]

    for line in handle:
        if is_alignment_start(line):
            yield block
            block = [line]
        else:
            block.append(line)

    if block:
        yield block


def get_reference_coordinates(block):
    """
    Read the first s-line in a MAF block.

    Returns:
        start, size, srcSize, srcName
    """
    for line in block:
        stripped = line.lstrip()

        if not stripped.startswith("s "):
            continue

        fields = stripped.split()

        if len(fields) < 7:
            raise ValueError(
                "Malformed MAF s-line:\n{}".format(line.rstrip())
            )

        src_name = fields[1]
        start = int(fields[2])
        size = int(fields[3])
        src_size = int(fields[5])

        return start, size, src_size, src_name

    raise ValueError("No s-line found in MAF alignment block.")


def normalize_block_text(block):
    """
    Ensure that every written alignment block ends with one blank line.
    """
    text = "".join(block).rstrip("\n")
    return text + "\n\n"


###############################################################################
# 1. Read header and locate the first alignment block
###############################################################################

with open(input_maf, "r", encoding="utf-8") as fin:

    header_lines = []
    first_a_line = None

    for line in fin:
        if is_alignment_start(line):
            first_a_line = line
            break
        header_lines.append(line)

    if first_a_line is None:
        print(
            "Error: no alignment block was found in {}".format(input_maf),
            file=sys.stderr,
        )
        sys.exit(1)

    block_iterator = iter_blocks(fin, first_a_line)

    try:
        first_block = next(block_iterator)
    except StopIteration:
        print(
            "Error: no alignment block was found in {}".format(input_maf),
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        first_start, first_size, chromosome_length, reference_name = \
            get_reference_coordinates(first_block)
    except Exception as exc:
        print(
            "Error while reading first reference block: {}".format(exc),
            file=sys.stderr,
        )
        sys.exit(1)

    if chromosome_length <= 0:
        print(
            "Error: invalid reference srcSize: {}".format(chromosome_length),
            file=sys.stderr,
        )
        sys.exit(1)

    ###########################################################################
    # 2. Determine all chunk intervals
    ###########################################################################

    num_chunks = (chromosome_length + sep_length - 1) // sep_length

    chunk_paths = {}

    for index in range(1, num_chunks + 1):

        start = (index - 1) * sep_length
        end = index * sep_length

        if end > chromosome_length:
            end = chromosome_length

        filename = "{}_{}_{}.maf".format(index, start, end)
        path = os.path.join(output_folder, filename)

        chunk_paths[index] = path

    print("Input MAF: {}".format(input_maf))
    print("Reference: {}".format(reference_name))
    print("Reference srcSize: {}".format(chromosome_length))
    print("Chunk length: {}".format(sep_length))
    print("Number of chunks: {}".format(num_chunks))
    print("Output folder: {}".format(output_folder))

    ###########################################################################
    # 3. Create all output files and write the original MAF header
    #
    # This is important because downstream commands expect every chunk file
    # to exist, even if that chunk contains no alignment blocks.
    ###########################################################################

    header_text = "".join(header_lines)

    if header_text and not header_text.endswith("\n"):
        header_text += "\n"

    if header_text and not header_text.endswith("\n\n"):
        header_text += "\n"

    for index in range(1, num_chunks + 1):
        with open(
            chunk_paths[index],
            "w",
            encoding="utf-8",
        ) as fout:
            fout.write(header_text)

    ###########################################################################
    # 4. Keep only a limited number of output files open simultaneously.
    #
    # A chromosome may have >900 chunks. Keeping all files open at once could
    # exceed the operating system file-descriptor limit.
    ###########################################################################

    MAX_OPEN_FILES = 64
    open_files = OrderedDict()

    def get_output_handle(index):

        if index in open_files:
            handle = open_files.pop(index)
            open_files[index] = handle
            return handle

        handle = open(
            chunk_paths[index],
            "a",
            encoding="utf-8",
            buffering=1024 * 1024,
        )

        open_files[index] = handle

        if len(open_files) > MAX_OPEN_FILES:
            old_index, old_handle = open_files.popitem(last=False)
            old_handle.close()

        return handle

    ###########################################################################
    # 5. Process each MAF block exactly once.
    ###########################################################################

    block_count = 0
    written_block_count = 0

    def process_block(block):

        nonlocal_dummy = None

        try:
            block_start, block_size, block_src_size, block_ref = \
                get_reference_coordinates(block)
        except Exception as exc:
            raise RuntimeError(
                "Failed to parse MAF block: {}".format(exc)
            )

        # The chromosome size should be identical in all reference rows.
        if block_src_size != chromosome_length:
            raise RuntimeError(
                "Reference srcSize changed within the same chromosome MAF: "
                "{} != {}".format(
                    block_src_size,
                    chromosome_length,
                )
            )

        block_end = block_start + block_size

        #######################################################################
        # Reproduce the OLD script's overlap condition exactly:
        #
        #     block_start <= chunk_end
        # and
        #     block_end >= chunk_start
        #
        # Therefore, a block that touches a boundary is present in both
        # neighboring chunks.
        #######################################################################

        # Earliest chunk whose end >= block_start.
        #
        # For example:
        #   block_start = 218505
        #
        # belongs to chunk 1 as well as chunk 2, matching the old program.
        first_chunk = (block_start + sep_length - 1) // sep_length

        if first_chunk < 1:
            first_chunk = 1

        # Latest chunk whose start <= block_end.
        last_chunk = block_end // sep_length + 1

        if last_chunk > num_chunks:
            last_chunk = num_chunks

        if first_chunk > num_chunks:
            return 0

        if last_chunk < 1:
            return 0

        text = normalize_block_text(block)

        n_written = 0

        for chunk_index in range(first_chunk, last_chunk + 1):

            handle = get_output_handle(chunk_index)
            handle.write(text)

            n_written += 1

        return n_written

    ###########################################################################
    # First block
    ###########################################################################

    block_count += 1
    written_block_count += process_block(first_block)

    ###########################################################################
    # Remaining blocks
    ###########################################################################

    for block in block_iterator:

        block_count += 1
        written_block_count += process_block(block)

    ###########################################################################
    # Close cached output files
    ###########################################################################

    for handle in open_files.values():
        handle.close()

    open_files.clear()


print("Finished.")
print("Alignment blocks read: {}".format(block_count))
print("Block writes: {}".format(written_block_count))
print(
    "Generated {} chunk MAF files.".format(num_chunks)
)
