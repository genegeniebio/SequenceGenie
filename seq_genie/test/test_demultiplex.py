'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=wrong-import-order
import os
import unittest

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import pandas as pd
from seq_genie import demultiplex


class Test(unittest.TestCase):
    '''Class to test utils module.'''

    def test_demuliplex_simple(self):
        '''Test demuliplex method.'''
        barcodes = [('AAAAAAGGGGGG', 'AAAAAAGGGGGG')]

        directory = os.path.dirname(os.path.realpath(__file__))
        filename = os.path.join(directory, 'simple_seqs.txt')

        with open(filename) as fle:
            seqs = [SeqRecord(line.strip()) for line in fle.readlines()]

        barcode_seqs = demultiplex.demultiplex(barcodes, seqs, tolerance=1)

        self.assertEqual(len(barcode_seqs[barcodes[0]]), 2)

    def test_demuliplex_complex(self):
        '''Test demuliplex method.'''
        directory = os.path.dirname(os.path.realpath(__file__))

        barcodes_df = \
            pd.read_csv(os.path.join(directory, 'barcodes.csv'))

        barcodes = [tuple(pair)
                    for pair in
                    barcodes_df[['forward', 'reverse']].values.tolist()]

        with open(os.path.join(directory, 'reads.fasta'), 'rU') as fle:
            seqs = list(SeqIO.parse(fle, 'fasta'))

        barcode_seqs = demultiplex.demultiplex(barcodes, seqs, tolerance=10)

        bc_pair = ('GAGTCTTGTGTCCCAGTTACCAGG', 'CGGGCCCTTCATCTCTCAGCCGAT')
        self.assertEqual(len(barcode_seqs[bc_pair]), 406)


if __name__ == "__main__":
    unittest.main()
