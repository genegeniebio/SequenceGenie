'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-member
# pylint: disable=no-name-in-module
import random
import sys

from Bio.Seq import Seq

import numpy as np
from seq_genie.demultiplex import _check_barcode


def compare(barcode1, seq):
    '''Compare barcodes.'''
    for tol in range(len(barcode1)):
        resp = _check_barcode(barcode1, barcode1, seq, len(seq), tol)

        if resp:
            return tol

    return float('NaN')


def test_all_versus_all(barcodes):
    '''Test all versus all.'''
    for idx1, barcode1 in enumerate(barcodes):
        for idx2 in range(idx1 + 1, len(barcodes)):
            compare(barcode1, barcodes[idx2])
            compare(barcode1, str(Seq(barcodes[idx2]).reverse_complement()))


def main(barcodes):
    '''main method.'''
    tols = []

    for _ in range(10000):
        seq = ''.join([random.choice(['A', 'C', 'G', 'T'])
                       for _ in range(128)])

        tols.append(compare(random.choice(barcodes), seq))

    print(np.min(tols))


if __name__ == '__main__':
    main(sys.argv[1:])
