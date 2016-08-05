import argparse
import gzip
import os
import re
import StringIO
import struct
import sys


def parse_arguments():
    """Parse the command line arguments.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description='tail for BAMs')
    parser.add_argument(
        'filenames',
        help='BAMs on which to perform the tail operation',
        nargs='+',
        metavar='FILE'
    )
    return parser.parse_args()


def gunzip_chunk(byte_stream):
    """Interpret a byte string as a gzipped chunk and return its
    unzipped contents.

    Args:
        byte_stream (bytes): Complete chunk of gzipped data to be ungzipped.

    Returns:
        Bytes representing the ungzipped chunk.
    """
    byte_stream = StringIO.StringIO(byte_stream)
    return gzip.GzipFile(fileobj=byte_stream).read()


def get_bam_chunk(filename, size, tailsize=100000):
    """Find the final (or, in rare cases, penultimate) gzipped chunk in the BAM
    file and return its unzipped contents.

    Args:
        filename: Path of BAM.
        size (int): Total filesize of BAM.
        tailsize (int, optional): Number of bytes to look at from the end of
            the BAM file, defaulting to 100000. This should be large enough
            to read an entire unbroken gzipped chunk (64k?).

    Returns:
        Bytes representing the ungzipped final complete chunk of data in the
            BAM file.
    """
    with open(filename, 'rb') as bamfile:
        bamfile.seek(size-tailsize)
        tail = bamfile.read()
    search_bytes = struct.pack('<4B', 31, 139, 8, 4)
    chunk_starts = [m.start() for m in re.finditer(search_bytes, tail)]
    if len(chunk_starts) < 2:
        raise Exception('not enough data')
    start, end = chunk_starts[-2:]
    return gunzip_chunk(tail[start:end])


def final_record_position(bam_data):
    """Find the reference index and position of the final record in the BAM
    chunk.

    Args:
        bam_data (bytes): Ungzipped contents of a chunk of BAM data.

    Returns:
        A tuple of ints, where the first int is the index of the reference
            sequence (from the header) or -1 for an unmapped record, and the
            second int is the position.
    """
    while True:
        block_size, ref_id, pos = struct.unpack_from('<3i', bam_data)
        bam_data = bam_data[block_size+4:]
        if not bam_data:
            break
    return ref_id, pos



def get_bam_ref_sequences(filename, headsize=100000):
    """Read the list of reference sequence names from the beginning of the BAM
    file.

    Args:
        filename: Path of BAM to get the reference sequences from.
        headsize (int, optional): number of bytes to read so that the complete
        reference sequence data is included. Defaults to 10000.

    Returns:
        List of reference sequence names.
    """
    with open(filename, 'rb') as bamfile:
        head = bamfile.read(headsize)
    fields = struct.unpack_from('<4BI2BH2B2H', head)
    block_size = fields[-1] + 1
    bam_data = gunzip_chunk(head[:block_size])
    magic, l_text = struct.unpack_from('<4si', bam_data)
    if magic != 'BAM\x01':
        raise Exception('not a BAM file')
    bam_data = bam_data[8+l_text:]
    n_ref, = struct.unpack_from('<i', bam_data)
    bam_data = bam_data[4:]
    ref_sequences = list()
    for _ in range(n_ref):
        l_name, = struct.unpack_from('<i', bam_data)
        _, name = struct.unpack_from('<i{:d}s'.format(l_name), bam_data)
        ref_sequences.append(name[:-1])
        bam_data = bam_data[8+l_name:]
    return ref_sequences


def process_file(filename):
    """Find the chromosome and position of the final record in the final
    chunk of the BAM file.

    Args:
        filename: Path of BAM on which to perform tail.

    Returns:
        String of the position of the final record of the BAM, e.g. 'chr1:100',
            or 'unmapped' if the last record is not mapped to a reference
            sequence.
    """
    ref_sequences = get_bam_ref_sequences(filename)
    size = os.path.getsize(filename)
    bam_data = get_bam_chunk(filename, size)
    chr, pos = final_record_position(bam_data)
    if chr == -1:
        return 'unmapped'
    else:
        chr = ref_sequences[chr]
        pos += 1
        return '{}:{}'.format(chr, pos)


def main(filenames):
    """Perform the tail operation on all given BAMs. If more than one BAM is
    given, prefix the output with the filename.

    Args:
        filenames: A sequence of paths to BAMs on which to perform tail.
    """
    for filename in filenames:
        if len(filenames) > 1:
            print '{}: '.format(filename),
        print process_file(filename)


if __name__ == '__main__':
    args = parse_arguments()
    main(args.filenames)
