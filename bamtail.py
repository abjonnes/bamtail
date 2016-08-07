#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import io
import os
import re
import struct
import sys


"""Perform a ``tail``-like operation on BAM files.

This tool is used to return the position of an alignment near the end of BAM
file(s). This is useful for monitoring progress during generation of large BAM
files where the resulting file is sorted by position.

The BAM specification is referenced heavily throughout. It can be found at:
https://samtools.github.io/hts-specs/SAMv1.pdf

"""

__version__ = '1.0'


def parse_arguments():
    """Parse the command line arguments.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description='tail for BAMs')
    parser.add_argument(
        'filenames',
        help='BAMs on which to perform the tail operation',
        nargs='*',
        metavar='FILE'
    )
    parser.add_argument(
        '--version',
        '-v',
        help='print the version',
        action='store_true'
    )
    return parser.parse_args()


def decompress_block(byte_stream):
    """Interpret a byte string as a BGZF block and return its decompressed
    contents.

    BGZF format is gunzip compatible, so the gzip module can perform the
    decompression.

    Args:
        byte_stream (bytes): Complete BGZF block to be decompressed.

    Returns:
        Bytes representing the decompressed block.
    """
    byte_stream = io.BytesIO(byte_stream)
    return gzip.GzipFile(fileobj=byte_stream).read()


def get_final_bam_block(filename, size, tailsize=124000):
    """Return the decompressed data from the final BGZF block in the BAM file.

    Data is read from the end of the file, and start positions of BGZF blocks
    are identified by a sequence of 8-bit integers as outlined in the BAM spec.
    The final complete BGZF block is identified as the data between the
    penultimate and final occurrences of this sequence. In the case that there
    are no incomplete blocks at the end of the file (such as when the file is
    complete), the data returned is from the penultimate rather than final
    BGZF block.

    Args:
        filename: Path of BAM.
        size (int): Total filesize of BAM.
        tailsize (int, optional): Number of bytes to look at from the end of
            the BAM file, defaulting to 124000. This should be large enough
            to read an entire unbroken BGZF block in addition to any incomplete
            block at the end of the file.

    Returns:
        Bytes representing the ungzipped final complete block of data in the
            BAM file.

    Raises:
        RuntimeError: Fewer than two BGZF block beginnings were discovered,
            indicating that not enough data was read.
    """
    with open(filename, 'rb') as bamfile:
        bamfile.seek(size-tailsize)
        tail = bamfile.read(tailsize)
    search_bytes = struct.pack('<4B', 31, 139, 8, 4)  # from BAM spec
    block_starts = [m.start() for m in re.finditer(search_bytes, tail)]
    if len(block_starts) < 2:
        raise RuntimeError('not enough data')
    start, end = block_starts[-2:]
    return decompress_block(tail[start:end])


def final_alignment_position(bam_data):
    """Find the reference index and position of the final alignment in the BGZF
    block.

    Assume that alignment blocks do not span BGZF blocks (this is a strong
    assumption (!), but it seems to be true for every BAM file I've
    encountered). Iterate through all the alignment blocks until all data is
    read and return the values for the final alignment.

    Args:
        bam_data (bytes): Decompressed contents of a BGZF block.

    Returns:
        A tuple of ints, where the first int is the index of the reference
            sequence (from the header) or -1 for an unmapped alignment, and the
            second int is the position.
    """
    while bam_data:
        block_size, ref_id, pos = struct.unpack_from('<3i', bam_data)
        bam_data = bam_data[block_size+4:]  # block_size excludes itself
    return ref_id, pos


def get_bam_ref_sequences(filename, headsize=64000):
    """Read the list of reference sequence names from the beginning of the BAM
    file.

    The first BGZF block in the BAM file is read and decompressed, and all the
    reference sequence names are read. Also ensures that the file is actually
    a BAM file by looking for the BAM magic string outlined in the spec. In
    cases with very many reference sequences such that they do not all fit
    within one BGZF block, this will fail! Consult the BAM spec for more detail
    regarding the format.

    Args:
        filename: Path of BAM to get the reference sequences from.
        headsize (int, optional): number of bytes to read so that the complete
        first BGZF block is read. Defaults to 64000.

    Returns:
        List of reference sequence names.

    Raises:
        RuntimeError: The file appears to not be a BAM file.
    """
    with open(filename, 'rb') as bamfile:
        head = bamfile.read(headsize)

    # read a bunch of fields, where the final field is the size of the block
    # minus 1, then decompress the block
    fields = struct.unpack_from('<4BI2BH2B2H', head)
    block_size = fields[-1] + 1
    bam_data = decompress_block(head[:block_size])

    offset = 0  # keep track how much data we've read
    magic, l_text = struct.unpack_from('<4si', bam_data)
    if magic != b'BAM\x01':  # magic string
        raise RuntimeError('{} appears to not be a BAM file'.format(filename))
    offset += 8 + l_text
    n_ref, = struct.unpack_from('<i', bam_data, offset)
    offset += 4
    ref_sequences = list()
    for _ in range(n_ref):
        l_name, = struct.unpack_from('<i', bam_data, offset)
        offset += 4
        name, = struct.unpack_from('<{:d}s'.format(l_name), bam_data, offset)
        ref_sequences.append(name[:-1])  # sequence names are NUL terminated
        offset += l_name + 4  # skipping l_ref
    return ref_sequences


def process_file(filename):
    """Find the chromosome and position of the final alignment in the final
    BGZF block discovered in the BAM file.

    Args:
        filename: Path of BAM on which to perform tail.

    Returns:
        String of the position of the final alignment of the BAM, e.g.
            'chr1:100', or 'unmapped' if the last alignment is not mapped to a
            reference sequence.
    """
    ref_sequences = get_bam_ref_sequences(filename)
    size = os.path.getsize(filename)
    bam_data = get_final_bam_block(filename, size)
    chr, pos = final_alignment_position(bam_data)
    if chr == -1:
        return 'unmapped'
    else:
        chr = ref_sequences[chr].decode('utf-8')
        pos += 1  # BAM positions are 0-indexed
        return '{}:{}'.format(chr, pos)


def main(filenames):
    """Perform the tail operation on all given BAMs. If more than one BAM is
    given, prefix the output with the filename.

    Args:
        filenames: A sequence of paths to BAMs on which to perform tail.
    """
    for filename in filenames:
        if len(filenames) > 1:
            print('{}: '.format(filename), end='')
        print(process_file(filename))


if __name__ == '__main__':
    args = parse_arguments()
    if args.version:
        print('bamtail version {}'.format(__version__))
    if args.filenames:
        main(args.filenames)

