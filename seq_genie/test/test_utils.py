'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
import os
import unittest

from Bio.SeqRecord import SeqRecord

from seq_genie import utils


class Test(unittest.TestCase):
    '''Class to test utils module.'''

    def test_bin_seqs(self):
        '''Test bin_seqs method.'''
        barcodes = [('AAAAAAGGGGGG', 'AAAAAAGGGGGG')]

        directory = os.path.dirname(os.path.realpath(__file__))
        filename = os.path.join(directory, 'seqs.txt')

        with open(filename) as fle:
            seqs = [SeqRecord(line.strip()) for line in fle.readlines()]

        barcode_seqs = utils.bin_seqs(barcodes, seqs, search_len=20)

        self.assertEqual(len(barcode_seqs[barcodes[0]]), 2)


if __name__ == "__main__":
    unittest.main()
