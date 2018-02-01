'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes
from collections import defaultdict
import math
import multiprocessing
import sys
import time

from Bio import pairwise2, SeqIO, Seq
from synbiochem.utils import mut_utils, seq_utils

import numpy as np


_NUCL_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


class Aligner(object):
    '''Aligner class.'''

    def __init__(self, wt_filename, seqs_filename, mut_strs=None,
                 match=1, mismatch=-1,
                 gapopen=-4, gapextend=-1, wt_prob=0.1):
        self.__match = match
        self.__mismatch = mismatch
        self.__gapopen = gapopen
        self.__gapextend = gapextend
        self.__wt_prob = wt_prob

        # Read template sequence:
        with open(wt_filename, 'rU') as fle:
            for record in SeqIO.parse(fle, 'fasta'):
                self.__wt_seq = str(record.seq)

        # Read sequences:
        with open(seqs_filename, 'rU') as fle:
            self.__seqs = [[record.id, str(record.seq)]
                           for record in SeqIO.parse(fle, 'fasta')]

        self.__nucl_probs = _get_nucl_probs()
        self.__pos_spec_probs = self.__get_pos_spec_probs(mut_strs)

    def align(self):
        '''Aligns sequences.'''
        # Multi processor:
        start = time.time()
        queue = multiprocessing.Queue()
        nprocs = multiprocessing.cpu_count()
        procs = []

        for chunk in _get_chunks(self.__seqs, nprocs):
            proc = multiprocessing.Process(target=self.__align,
                                           args=(chunk, queue))
            proc.Daemon = True
            procs.append(proc)
            proc.start()

        results = {}

        for _ in range(nprocs):
            results.update(queue.get())

        print len(results)
        print time.time() - start

    def __align(self, seqs, queue):
        '''Aligns sequences.'''
        alignments = {}

        for seq_id, seq in seqs:
            best_aln = self.__get_align(seq)
            aln = self.__get_align(Seq.Seq(seq).reverse_complement())

            if aln[2] > best_aln[2]:
                best_aln = aln

            alignments[seq_id] = best_aln

        queue.put(alignments)

    def __get_pos_spec_probs(self, mut_strs):
        '''Gets position specific probabilities.'''
        pos_spec_probs = []
        mut_probs = self.__get_mut_probs(mut_strs)

        for idx, res in enumerate(self.__wt_seq):
            if mut_probs[idx]:
                prob_mat = [self.__nucl_probs[res]]
                prob_mat.extend(mut_probs[idx])
                weights = [self.__wt_prob] + \
                    [(1 - self.__wt_prob) / (len(prob_mat) - 1)] * \
                    (len(prob_mat) - 1)
                pos_spec_prob = list(np.dot(np.transpose(prob_mat), weights))
            else:
                pos_spec_prob = self.__nucl_probs[res]

            pos_spec_probs.append(pos_spec_prob)

        return pos_spec_probs

    def __get_mut_probs(self, mut_strs):
        '''Gets mutations.'''
        mut_probs = defaultdict(list)

        for mut in _parse_mut_strs(mut_strs):
            mut_probs[mut[1]].append(self.__nucl_probs[mut[2]])

        return mut_probs

    def __get_align(self, seq):
        '''Gets alignment from seq.'''
        matcher = Matcher(len(seq), self.__pos_spec_probs)

        for aln in pairwise2.align.globalcs(self.__wt_seq,
                                            seq,
                                            matcher.match_fn,
                                            self.__gapopen,
                                            self.__gapextend,
                                            one_alignment_only=True,
                                            penalize_end_gaps=False):
            return aln

        return ('', '', -1, -1)


def _get_chunks(lst, num):
    '''Split list into approximately equal chunks.'''
    chunks = []

    avg = len(lst) / float(num)
    last = 0.0

    while last < len(lst):
        chunks.append(lst[int(last):int(last + avg)])
        last += avg

    return chunks


def _get_nucl_probs():
    '''Gets nucleotide probabilities.'''
    probs = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])

    for ambig, nucls in seq_utils.INV_NUCL_CODES.iteritems():
        for nucl in nucls:
            probs[ambig][_NUCL_IDX[nucl]] = 1.0 / len(nucls)

    return probs


def _parse_mut_strs(mut_strs):
    '''Parses mutation strings.'''
    mutations = []

    if mut_strs is not None:
        mutations = [mut_utils.parse_mut_str(mut) for mut in mut_strs]
        mutations = [item for sublist in mutations for item in sublist]

    return mutations


def _process_aln(seq_id, aln):
    '''Process alignment.'''
    alignment = []
    start = float('NaN')

    for idx, pos in enumerate(zip(*(aln[0], aln[1]))):
        alignment.append('|' if pos[0] == pos[1] else ' ')

        if math.isnan(start) and pos[0] != '-':
            start = idx

    print seq_id
    print aln[0][start:]
    print ''.join(alignment)[start:]
    print aln[1][start:]
    print


class Matcher(object):
    '''Matcher class implements match_fn.'''

    def __init__(self, seq_len, pos_spec_probs):
        self.__seq_len = seq_len
        self.__pos_spec_probs = pos_spec_probs
        self.__count = 0

    def match_fn(self, _, y):
        '''match function.'''
        wt_pos = self.__count / self.__seq_len
        prob = self.__pos_spec_probs[wt_pos][_NUCL_IDX[y]]
        self.__count = self.__count + 1
        return prob


def main(args):
    '''main method.'''
    mutations = args[2:] if len(args) > 2 else None
    Aligner(args[0], args[1], mutations).align()


if __name__ == '__main__':
    main(sys.argv[1:])
