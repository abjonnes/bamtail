#!/usr/bin/env python

import gzip
import os
import re
import StringIO
import struct
import sys


def gunzip_chunk(byte_stream):
    byte_stream = StringIO.StringIO(byte_stream)
    return gzip.GzipFile(fileobj=byte_stream).read()


def get_bam_chunk(filename, size, tailsize=100000):
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
    while True:
        block_size, ref_id, pos = struct.unpack_from('<3i', bam_data)
        bam_data = bam_data[block_size+4:]
        if not bam_data:
            break
    return ref_id, pos



def get_bam_ref_sequences(filename, headsize=100000):
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
    for filename in filenames:
        if len(filenames) > 1:
            print '{}: '.format(filename),
        print process_file(filename)

if __name__ == '__main__':
    import sys

    filenames = sys.argv[1:]
    main(filenames)
