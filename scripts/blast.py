'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-member
# pylint: disable=no-name-in-module
# pylint: disable=wrong-import-order
import sys

from Bio import SeqIO
from synbiochem.utils import seq_utils

import pandas as pd


def main(args):
    '''main method.'''
    # Align:
    subjects = {seq.id: str(seq.seq) for seq in SeqIO.parse(args[0], 'fasta')}

    # TODO: apply virtual PCR of query sequence
    queries = {seq.id: str(seq.seq) for seq in SeqIO.parse(args[1], 'fasta')}

    pd.DataFrame([[seq_id, len(seq)] for seq_id, seq in subjects.items()],
                 columns=['id', 'len']).to_csv('len.csv', index=False)

    for result in seq_utils.do_blast(subjects, queries):
        for alignment in result.alignments:
            for hsp in alignment.hsps:
                print(result.query + '\t' + alignment.title + '\t' +
                      str(hsp.align_length) + '\t' + str(hsp.gaps) + '\t' +
                      str(result.query_length) + '\t' + str(alignment.length))


if __name__ == '__main__':
    main(sys.argv[1:])
