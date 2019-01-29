'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=superfluous-parens
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=wrong-import-order
from collections import OrderedDict
import os.path
from threading import Thread

from Bio import Seq, SeqIO

import multiprocessing as mp
from seq_genie import reads


class ReadThread(Thread):
    '''Thread-safe class to write demultiplexed reads to Fasta.'''

    def __init__(self, queue, parent_dir):
        self.__queue = queue
        self.__parent_dir = parent_dir
        self.__closed = False
        self.__files = {}
        Thread.__init__(self)

    def run(self):
        '''Run.'''
        while not self.__closed:
            task = self.__queue.get()
            self.__write(task)
            self.__queue.task_done()

    def get_filenames(self):
        '''Get filenames.'''
        return {barcodes: fle.name
                for barcodes, fle in self.__files.iteritems()}

    def close(self):
        '''Close.'''
        for fle in self.__files.values():
            fle.close()

        self.__closed = True

    def __write(self, task):
        barcodes = task[0]

        if barcodes not in self.__files:
            dir_name = os.path.join(self.__parent_dir, '_'.join(barcodes))
            filename = os.path.join(dir_name, 'reads.fasta')
            os.makedirs(dir_name)
            self.__files[barcodes] = open(filename, 'w')

        SeqIO.write(task[1], self.__files[barcodes], 'fasta')


def demultiplex(barcodes, in_dir, min_length, max_read_files, out_dir,
                tolerance, num_threads, search_len=48):
    '''Bin sequences according to barcodes.'''
    max_barcode_len = max([len(barcode)
                           for pair in barcodes
                           for barcode in pair])

    pool = mp.Pool(processes=num_threads)
    write_queue = mp.Manager().Queue()
    read_thread = ReadThread(write_queue, out_dir)
    read_thread.start()

    results = [pool.apply_async(_bin_seqs, args=(fle,
                                                 min_length,
                                                 max_barcode_len,
                                                 search_len,
                                                 _format_barcodes(barcodes),
                                                 tolerance,
                                                 idx,
                                                 max_read_files,
                                                 write_queue))
               for idx, fle in enumerate(reads.get_filenames(in_dir,
                                                             max_read_files))]
    for res in results:
        res.get()

    read_thread.close()
    return read_thread.get_filenames()


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
              barcodes, tolerance, idx, max_read_files, write_queue):
    '''Bin a batch of sequences.'''
    barcode_seqs = 0

    seqs = reads.get_reads(reads_filename, min_length)

    for seq in seqs:
        if seq:
            for pairs in barcodes:
                selected_barcodes = [None, None]

                if _check_seq(seq, max_barcode_len, search_len, pairs,
                              selected_barcodes, tolerance, write_queue):
                    barcode_seqs += 1
                    break

    _report_barcodes(idx, max_read_files, len(seqs), barcode_seqs)

    return None


def _check_seq(seq, max_barcode_len, search_len, pairs, selected_barcodes,
               tolerance, write_queue):
    '''Check sequence against barcode sequences.'''
    seq_len = min(max_barcode_len + search_len, len(seq))
    seq_start = list(seq.seq[:seq_len])
    seq_end = list(seq.seq[-(seq_len):])

    # Check all barcodes:
    for orig, bc_pair in pairs.iteritems():
        _check_pair(orig, bc_pair, [seq_start, seq_end], seq_len,
                    selected_barcodes, tolerance)

        if selected_barcodes[0] and selected_barcodes[1]:
            write_queue.put([tuple(selected_barcodes), seq])
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


def _report_barcodes(idx, max_read_files, num_seqs, barcode_seqs):
    '''Report barcodes.'''
    if barcode_seqs:
        print('Seqs: %d/%d\tMatched: %d/%d' % ((idx + 1),
                                               max_read_files,
                                               barcode_seqs,
                                               num_seqs))
    else:
        print('Seqs: %d' % (idx + 1))
