'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-member
# pylint: disable=no-name-in-module
# pylint: disable=superfluous-parens
import sys

from Bio import SeqIO

from synbiochem.utils import seq_utils


def main(args):
    '''main method.'''
    # Align:
    subjects = list(SeqIO.parse(args[0], 'fasta'))
    queries = list(SeqIO.parse(args[1], 'fasta'))

    print(seq_utils.do_blast(subjects, queries))


if __name__ == '__main__':
    main(sys.argv[1:])
