'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=superfluous-parens
# pylint: disable=too-many-arguments
# pylint: disable=wrong-import-order
from collections import defaultdict, OrderedDict

from Bio import Seq

import multiprocessing as mp
import numpy as np
from seq_genie import reads


def demultiplex(barcodes, in_dir, min_length, max_read_files, tolerance,
                num_threads, search_len=48):
    '''Bin sequences according to barcodes.'''
    max_barcode_len = max([len(barcode)
                           for pair in barcodes
                           for barcode in pair])

    pool = mp.Pool(processes=num_threads)

    results = [pool.apply_async(_bin_seqs, args=(fle,
                                                 min_length,
                                                 max_barcode_len,
                                                 search_len,
                                                 _format_barcodes(barcodes),
                                                 tolerance,
                                                 idx,
                                                 max_read_files))
               for idx, fle in enumerate(reads.get_filenames(in_dir,
                                                             max_read_files))]

    return _consolodate_bc_seqs([res.get()
                                 for res in results])


def _format_barcodes(barcodes):
    '''Format barcodes to reduce number of get_rev_complement calls.'''
    for_brcds = OrderedDict()
    rev_brcds = OrderedDict()

    for pair in barcodes:
        for_brcds[pair] = \
            [list(pair[0]), list(Seq.Seq(pair[1]).reverse_complement())]
        rev_brcds[pair] = \
            [list(pair[1]), list(Seq.Seq(pair[0]).reverse_complement())]

    return for_brcds, rev_brcds


def _bin_seqs(reads_filename, min_length, max_barcode_len, search_len,
              barcodes, tolerance, idx, max_read_files):
    '''Bin a batch of sequences.'''
    barcode_seqs = defaultdict(list)

    seqs = reads.get_reads(reads_filename, min_length)

    for seq in seqs:
        if seq:
            for pairs in barcodes:
                selected_barcodes = [None, None]

                if _check_seq(seq, max_barcode_len, search_len, pairs,
                              barcode_seqs, selected_barcodes, tolerance):
                    break

    _report_barcodes(idx, max_read_files, len(seqs), barcode_seqs)

    return barcode_seqs


def _check_seq(seq, max_barcode_len, search_len, pairs, barcode_seqs,
               selected_barcodes, tolerance):
    '''Check sequence against barcode sequences.'''
    seq_len = min(max_barcode_len + search_len, len(seq))
    seq_start = list(seq.seq[:seq_len])
    seq_end = list(seq.seq[-(seq_len):])

    # Check all barcodes:
    for orig, bc_pair in pairs.iteritems():
        _check_pair(orig, bc_pair, [seq_start, seq_end], seq_len,
                    selected_barcodes, tolerance)

        if selected_barcodes[0] and selected_barcodes[1]:
            barcode_seqs[tuple(selected_barcodes)].append(seq)
            return True

    return False


def _check_pair(orig, pair, seqs, seq_len, selected_barcodes, tolerance):
    '''Check similarity scores.'''
    for idx in range(2):
        if not selected_barcodes[idx]:
            resp = _check_barcode(orig[idx], pair[idx], seqs[idx], seq_len,
                                  tolerance)

            if resp:
                selected_barcodes[idx] = resp


def _check_barcode(orig, barcode, seq, seq_len, tolerance):
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


def _consolodate_bc_seqs(all_barcode_seqs):
    '''Consolidate all barcode_seqs.'''
    barcode_seqs = defaultdict(list)

    for dct in all_barcode_seqs:
        for k, v in dct.items():
            barcode_seqs[k].extend(v)

    return barcode_seqs


def _report_barcodes(idx, max_read_files, num_seqs, barcode_seqs):
    '''Report barcodes.'''
    seq_lens = [len(seq) for seq in barcode_seqs.values()]

    if seq_lens:
        vals = ((idx + 1),
                max_read_files,
                sum(seq_lens),
                num_seqs,
                min(seq_lens),
                max(seq_lens),
                np.mean(seq_lens),
                np.median(seq_lens))

        s = 'Seqs: %d/%d\tMatched: %d/%d\tRange: %d-%d\tMean: %.1f\tMedian: %d'
        print(s % vals)
    else:
        print('Seqs: %d' % (idx + 1))
