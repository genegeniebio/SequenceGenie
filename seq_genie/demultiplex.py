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

from Bio import Seq
from synbiochem.utils import thread_utils

import numpy as np


def demultiplex(barcodes, sequences, tolerance, search_len=48, num_threads=0,
                batch_size=64):
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
                                     tolerance)

            thread_pool.wait_completion()
        else:
            for idx, seqs in enumerate(_get_batch(sequences, batch_size)):
                _report_barcodes(idx, num_seqs, batch_size, barcode_seqs)

                bin_seqs(seqs, max_barcode_len, search_len, barcodes,
                         barcode_seqs, float(tolerance))
    else:
        barcode_seqs['undefined'].extend(sequences)

    return barcode_seqs


def format_barcodes(barcodes):
    '''Format barcodes to reduce number of get_rev_complement calls.'''
    for_brcds = {}
    rev_brcds = {}

    for pair in barcodes:
        for_brcds[pair] = \
            [list(pair[0]), list(Seq.Seq(pair[1]).reverse_complement())]
        rev_brcds[pair] = \
            [list(pair[1]), list(Seq.Seq(pair[0]).reverse_complement())]

    return for_brcds, rev_brcds


def bin_seqs(seqs, max_barcode_len, search_len, barcodes, barcode_seqs,
             tolerance):
    '''Bin a batch of sequences.'''
    for seq in seqs:
        if seq:
            for pairs in barcodes:
                selected_barcodes = [None, None]

                if check_seq(seq, max_barcode_len, search_len, pairs,
                             barcode_seqs, selected_barcodes, tolerance):
                    break


def check_seq(seq, max_barcode_len, search_len, pairs, barcode_seqs,
              selected_barcodes, tolerance):
    '''Check sequence against barcode sequences.'''
    seq_len = min(max_barcode_len + search_len, len(seq))
    seq_start = list(seq.seq[:seq_len])
    seq_end = list(seq.seq[-(seq_len):])

    # Check all barcodes:
    for orig, bc_pair in pairs.iteritems():
        check_pair(orig, bc_pair, [seq_start, seq_end], seq_len,
                   selected_barcodes, tolerance)

        if selected_barcodes[0] and selected_barcodes[1]:
            barcode_seqs[tuple(selected_barcodes)].append(seq)
            return True

    return False


def check_pair(orig, pair, seqs, seq_len, selected_barcodes, tolerance):
    '''Check similarity scores.'''
    for idx in range(2):
        if not selected_barcodes[idx]:
            resp = check_barcode(orig[idx], pair[idx], seqs[idx], seq_len,
                                 tolerance)

            if resp:
                selected_barcodes[idx] = resp


def check_barcode(orig, barcode, seq, seq_len, tolerance):
    '''Check barcode.'''
    bc_len = len(barcode)

    for substr in [seq[i:i + bc_len] for i in range(seq_len - bc_len + 1)]:
        diff = 0

        for idx, s in enumerate(substr):
            if s != barcode[idx]:
                diff += 1

                if diff > tolerance:
                    break

        if diff <= tolerance:
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
