'''
DNA++ (c) DNA++ 2017

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-name-in-module
import os
import subprocess
import tempfile
from Bio import SeqIO

from pysam import AlignmentFile


def get_reads(dirs, min_length=0, dir_filter='pass'):
    '''Converts fast5 files to fasta.'''
    reads = []

    for directory in dirs:
        for dirpath, _, filenames in os.walk(os.path.abspath(directory)):
            for filename in filenames:
                filename = os.path.join(dirpath, filename)

                if filename[-6:] == '.fastq' and dir_filter in filename:
                    with open(filename, 'rU') as fle:
                        reads.extend([record
                                      for record in SeqIO.parse(fle, 'fastq')
                                      if len(record.seq) > min_length])

    return reads


def align(templ_filename, reads, out='align.sam'):
    '''Aligns sequences in barcoded bins.'''
    # Index template:
    subprocess.call(['bwa', 'index', templ_filename])

    # Align and sort:
    sort(mem(reads, templ_filename), out)


def mem(reads, templ_filename, readtype='pacbio'):
    '''Runs BWA MEM.'''
    reads_file = tempfile.NamedTemporaryFile(delete=False)
    out_file = tempfile.NamedTemporaryFile(delete=False)

    SeqIO.write(reads, reads_file.name, 'fasta')

    with open(out_file.name, 'w') as out:
        subprocess.call(['bwa', 'mem',
                         '-x', readtype,
                         templ_filename, reads_file.name],
                        stdout=out)

    return out_file.name


def sort(in_filename, out_filename):
    '''Custom sorts SAM file.'''
    sam_file = AlignmentFile(in_filename, 'r')
    out_file = AlignmentFile(out_filename, 'wh',
                             template=sam_file,
                             header=sam_file.header)

    for read in sorted([read for read in sam_file],
                       key=lambda x: (-x.query_length,
                                      x.reference_start)):
        out_file.write(read)

    out_file.close()

    return out_filename
