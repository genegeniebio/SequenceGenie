'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=wrong-import-order
from collections import defaultdict
from itertools import izip_longest
import sys
from Bio import Seq
from synbiochem.utils import thread_utils

import numpy as np


def demultiplex(barcodes, sequences, window_size, search_len=32,
                num_threads=0, batch_size=32):
    '''Bin sequences according to barcodes.'''
    barcode_seqs = defaultdict(list)

    num_seqs = len(sequences)

    max_barcode_len = max([len(barcode)
                           for pair in barcodes
                           for barcode in pair])

    if barcodes:
        barcodes = format_barcodes(barcodes)

        if num_threads:
            thread_pool = thread_utils.ThreadPool(num_threads)

            for idx, seqs in enumerate(_get_batch(sequences, batch_size)):
                _report_barcodes(idx, num_seqs, batch_size, barcode_seqs)

                thread_pool.add_task(bin_seqs, seqs, max_barcode_len,
                                     search_len, barcodes, barcode_seqs,
                                     window_size)

            thread_pool.wait_completion()
        else:
            for idx, seqs in enumerate(_get_batch(sequences, batch_size)):
                _report_barcodes(idx, num_seqs, batch_size, barcode_seqs)

                bin_seqs(seqs, max_barcode_len, search_len, barcodes,
                         barcode_seqs, window_size)
    else:
        barcode_seqs['undefined'].extend(sequences)

    return barcode_seqs


def format_barcodes(barcodes):
    '''Format barcodes to reduce number of get_rev_complement calls.'''
    for_brcds = {}
    rev_brcds = {}

    for pair in barcodes:
        for_brcds[pair] = \
            [pair[0], str(Seq.Seq(pair[1]).reverse_complement())]
        rev_brcds[pair] = \
            [pair[1], str(Seq.Seq(pair[0]).reverse_complement())]

    return for_brcds, rev_brcds


def bin_seqs(seqs, max_barcode_len, search_len, barcodes, barcode_seqs,
             window_size):
    '''Bin a batch of sequences.'''
    for seq in seqs:
        if seq:
            for pairs in barcodes:
                if check_seq(seq, max_barcode_len, search_len, pairs,
                             barcode_seqs, window_size):
                    break


def check_seq(seq, max_barcode_len, search_len, pairs, barcode_seqs,
              window_size):
    '''Check sequence against barcode sequences.'''
    search_len = min(max_barcode_len + search_len, len(seq))
    seq_start = seq.seq[:search_len]
    seq_end = seq.seq[-(search_len):]
    selected_barcodes = [None, None]

    # Check all barcodes:
    for orig, bc_pair in pairs.iteritems():
        check_pair(orig, bc_pair, [seq_start, seq_end],
                   selected_barcodes, window_size)

        if all(selected_barcodes):
            barcode_seqs[tuple(selected_barcodes)].append(seq)
            return True

    return False


def check_pair(orig, pair, seqs,
               selected_barcodes, window_size):
    '''Check sliding window scores.'''
    for idx in range(2):
        if not selected_barcodes[idx]:
            resp = check_barcode(orig[idx], pair[idx], seqs[idx], window_size)

            if resp:
                selected_barcodes[idx] = resp


def check_barcode(orig, barcode, seq, window_size):
    '''Check barcode.'''
    for substr in [barcode[i:i + window_size]
                   for i in xrange(len(barcode) - window_size + 1)]:
        if substr in seq:
            return orig

    return None


def _report_barcodes(idx, num_seqs, batch_size, barcode_seqs):
    '''Report barcodes.'''
    seq_lens = [len(seq) for seq in barcode_seqs.values()]

    if seq_lens:
        vals = ((idx + 1) * batch_size,
                num_seqs,
                sum(seq_lens),
                min(seq_lens),
                max(seq_lens),
                np.mean(seq_lens),
                np.median(seq_lens))

        s = 'Seqs: %d/%d\tMatched: %d\tRange: %d-%d\tMean: %.1f\tMedian: %d'
        print s % vals
    else:
        print 'Seqs: %d/%d' % ((idx + 1) * batch_size,  num_seqs)


def _get_batch(iterable, batch_size=64, fillvalue=None):
    '''Get a batch.'''
    args = [iter(iterable)] * batch_size
    return izip_longest(*args, fillvalue=fillvalue)


def main(args):
    '''main.'''
    bin_seqs(*args)


if __name__ == '__main__':
    main(sys.argv[1:])
