'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=no-name-in-module
import sys
from Bio import pairwise2
from Bio.Seq import Seq
from synbiochem.utils import seq_utils


def compare(barcode, seq):
    '''Compare barcodes.'''
    for tol in range(len(barcode)):
        resp = _check_barcode(barcode, barcode, seq, len(seq), tol)

        if resp:
            return tol

    return float('NaN')


def _check_barcode(orig, barcode, seq, seq_len, tolerance):
    '''Check barcode.'''
    bc_len = len(barcode)

    for substr in [seq[i:i + bc_len] for i in range(seq_len - bc_len + 1)]:
        alignments = pairwise2.align.globalxs(barcode, substr, -1, -1)
        score = max([aln[2] for aln in alignments])

        if bc_len - score <= tolerance:
            return orig

    return None


def main(args):
    '''main method.'''
    barcodes = args[1:]

    seqs = seq_utils.read_fasta(args[0]).values()

    for seq in seqs:
        for barcode in barcodes:

            for s in [seq, str(Seq(seq).reverse_complement())]:
                seq_len = min(len(barcode) + 48, len(s))
                seq_start = s[:seq_len]

                resp = compare(barcode, seq_start)

                print barcode + '\t' + seq_start + '\t' + str(resp)

        print


if __name__ == '__main__':
    main(sys.argv[1:])
