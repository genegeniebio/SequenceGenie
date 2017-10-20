'''
DNA++ (c) DNA++ 2017

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-name-in-module
import sys

from Bio import Seq, SeqIO, SeqRecord
from synbiochem.utils import io_utils
from synbiochem.utils.seq_utils import INV_NUCL_CODES

from seq_genie import protein, utils


def _analyse(sam_files, templ_filename, mut_templ_filename,
             fasta_filename='fasta.fasta'):
    '''Analyse.'''
    fasta_filename = io_utils.get_filename(fasta_filename)
    templ_seq = list(utils.parse(templ_filename))[0].seq
    mut_templ_seq = list(utils.parse(mut_templ_filename))[0].seq

    # Replace N with wildtype:
    with open(fasta_filename, 'w') as fle:
        for sam_file in sam_files:
            SeqIO.write(_strip(sam_file, templ_seq, mut_templ_seq), fle,
                        'fasta')

    sam_files = protein.align(templ_filename, [fasta_filename])


def _strip(sam_file, templ_seq, mut_templ_seq):
    '''Replace N and spurious mutations with wildtype.'''
    for read in sam_file:
        # Perform mapping of nucl indices to remove N:
        if read.aligned_pairs:
            templ_idx = zip(*read.aligned_pairs)[1]
            prefix = templ_seq[:min([val for val in templ_idx if val])]
            suffix = templ_seq[max(templ_idx):]
            seq = ''.join([read.seq[pair[0]]
                           if _valid(pair, read, templ_seq, mut_templ_seq)
                           else templ_seq[pair[1]]
                           for pair in read.aligned_pairs
                           if pair[1] is not None])

            if seq:
                yield SeqRecord.SeqRecord(prefix + seq + suffix,
                                          read.qname, '', '')


def _valid(pair, read, templ_seq, mut_templ_seq):
    '''Determines whether nucleotide pair is valid.'''
    return pair[0] and \
        read.seq[pair[0]] != 'N' and \
        (read.seq[pair[0]] in INV_NUCL_CODES[mut_templ_seq[pair[1]]] or
         read.seq[pair[0]] == templ_seq[pair[1]])


def main(args):
    '''main method.'''
    # Align:
    sam_files = protein.align(args[0], args[2:])
    _analyse(sam_files, args[0], args[1])


if __name__ == '__main__':
    main(sys.argv[1:])
