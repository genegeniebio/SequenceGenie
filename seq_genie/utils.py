'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-name-in-module
from Bio import Seq, SeqIO, SeqRecord
from synbiochem.utils import io_utils
from sbc_ngs import utils


def reject_indels(sam_filename_in, templ_filename, sam_filename_out):
    '''Rejects indels.'''
    sam_file = Samfile(sam_filename_in, 'r')
    out_file = Samfile(sam_filename_out, 'wh',
                       template=sam_file,
                       header=sam_file.header)
    templ_seq = get_seq(templ_filename)

    all_reads = 0
    passed_reads = 0

    for read in sam_file:
        all_reads += 1

        if read.cigarstring and str(len(templ_seq)) + 'M' in read.cigarstring:
            out_file.write(read)
            passed_reads += 1

    print('%s: %i/%i passed reject_indels filter' % (sam_filename_in,
                                                     passed_reads,
                                                     all_reads))

    out_file.close()


def replace_indels(sam_filename_in, templ_filename, sam_filename_out):
    '''Replace indels, replacing them with wildtype.'''
    sam_filename_out = io_utils.get_filename(sam_filename_out)
    templ_seq = get_seq(templ_filename)
    records = []

    all_reads = 0

    for read in Samfile(sam_filename_in, 'r'):
        # Perform mapping of nucl indices to remove spurious indels:
        all_reads += 1

        seq = ''.join([read.seq[pair[0]]
                       if pair[0]
                       else templ_seq[pair[1]]
                       for pair in read.aligned_pairs
                       if pair[1] is not None])

        if seq:
            records.append(SeqRecord.SeqRecord(Seq.Seq(seq), read.qname,
                                               '', ''))

    reads_filename = io_utils.get_filename(None)

    with open(reads_filename, 'w') as fle:
        SeqIO.write(records, fle, 'fasta')

    utils.mem(templ_filename, reads_filename,
              out_filename=sam_filename_out,
              gap_open=12)

    print('%s: %i/%i passed replace_indels filter' % (sam_filename_in,
                                                      len(records),
                                                      all_reads))

    return sam_filename_out


def get_seq(filename):
    '''Get sequence from Fasta file.'''
    for record in SeqIO.parse(filename, 'fasta'):
        return record.seq

    return None
