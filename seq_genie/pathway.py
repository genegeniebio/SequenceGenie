'''
SequenceGenie (c) GeneGenie Bioinformatics Ltd. 2017

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
from __future__ import division

import os
import random
import sys
import uuid

from Bio import SeqIO
import pysam
from synbiochem.utils import ice_utils, seq_utils, thread_utils

import pandas as pd
from seq_genie import utils


class PathwayAligner(object):
    '''Class to align NGS data to pathways.'''

    def __init__(self, barcodes_filename, reads_filename,
                 for_primer, rev_primer,
                 ice_url, ice_username, ice_password, ice_ids,
                 barcode_ice_filename=None,
                 min_length=128, max_seqs=128000):

        self.__primers = [for_primer, rev_primer]

        # Initialise project directory:
        self.__dir_name = str(uuid.uuid4())
        os.makedirs(self.__dir_name)

        # Demultiplex barcoded reads:
        reads = utils.get_reads(reads_filename, min_length)
        self.__barcode_reads = utils.bin_seqs(_get_barcodes(barcodes_filename),
                                              random.sample(reads,
                                                            min(len(reads),
                                                                max_seqs)))

        # Get pathway sequences from ICE:
        self.__ice_files = _get_ice_files(ice_url, ice_username, ice_password,
                                          ice_ids, self.__dir_name)

        # Get barcode-ICE dict:
        self.__barcode_ice = _get_barcode_ice(barcode_ice_filename)

        # Initialise dataframes:
        self.__score_df = pd.DataFrame(columns=self.__ice_files.keys(),
                                       index=self.__barcode_reads.keys())

        self.__mismatches_df = pd.DataFrame(columns=self.__ice_files.keys(),
                                            index=self.__barcode_reads.keys())

    def score_alignments(self, num_threads=8):
        '''Score alignments.'''
        for templ_filename in self.__ice_files.values():
            utils.index(templ_filename)

        thread_pool = thread_utils.ThreadPool(num_threads)

        for barcode, reads in self.__barcode_reads.iteritems():
            thread_pool.add_task(_score_alignment,
                                 self.__dir_name,
                                 barcode,
                                 reads,
                                 self.__ice_files,
                                 self.__primers[0],
                                 self.__primers[1],
                                 self.__barcode_ice,
                                 self.__score_df,
                                 self.__mismatches_df)

        thread_pool.wait_completion()

    def get_results(self):
        '''Get results.'''
        return self.__score_df, self.__mismatches_df


def _get_barcodes(filename):
    '''Get barcode sequences.'''
    barcodes_df = pd.read_csv(filename)
    return {row['well']: row['barcode'] for _, row in barcodes_df.iterrows()}


def _get_ice_files(url, username, password, ice_ids, dir_name):
    '''Get ICE sequences.'''
    ice_client = ice_utils.ICEClient(url, username, password)

    return {ice_id:
            seq_utils.write_fasta({ice_id:
                                   ice_client.get_ice_entry(ice_id).get_seq()},
                                  os.path.join(dir_name, ice_id + '.fasta'))
            for ice_id in ice_ids}


def _get_barcode_ice(barcode_ice_filename):
    '''Get barcode ice dict.'''
    barcode_ice = pd.read_csv(barcode_ice_filename,
                              dtype={'barcode': str, 'ice_id': str}) \
        if barcode_ice_filename else None

    return barcode_ice.set_index('barcode')['ice_id'].to_dict()


def _score_alignment(dir_name, barcode, reads, ice_files,
                     forward_primer, reverse_primer,
                     barcode_ice,
                     score_df, mismatches_df):
    '''Score an alignment.'''
    reads_filename = os.path.join(dir_name, barcode + '.fasta')
    SeqIO.write(reads, reads_filename, 'fasta')

    for ice_id, templ_filename in ice_files.iteritems():
        if not barcode_ice or barcode_ice[barcode] == ice_id:
            _score_barcode_ice(templ_filename, dir_name, barcode, ice_id,
                               reads_filename, forward_primer, reverse_primer,
                               score_df, mismatches_df)

    score_df.to_csv('score.csv')
    mismatches_df.to_csv('mismatches.csv')


def _score_barcode_ice(templ_filename, dir_name, barcode, ice_id,
                       reads_filename, forward_primer, reverse_primer,
                       score_df, mismatches_df):
    '''Score barcode ice pair.'''
    templ_seq = utils.get_reads(templ_filename)[0].seq

    sam_filename = os.path.join(dir_name,
                                barcode + '_' + ice_id + '.sam')

    cons_sam_filename = os.path.join(
        dir_name, barcode + '_' + ice_id + '_cons.sam')

    # Align:
    utils.align(templ_filename, reads_filename, sam_filename)

    # Generate then align consensus:
    fasta_filename = utils.get_consensus(sam_filename, templ_filename,
                                         forward_primer,
                                         reverse_primer)

    utils.align(templ_filename, fasta_filename, cons_sam_filename)

    matches, mismatches = utils.get_mismatches(cons_sam_filename,
                                               templ_seq)

    score_df[ice_id][barcode] = matches
    mismatches_df[ice_id][barcode] = mismatches


def _score(cons_sam_filename):
    '''Scores consensus alignment.'''
    cons_sam_file = pysam.AlignmentFile(cons_sam_filename, 'r')

    score = sum([len(read.positions)
                 for read in cons_sam_file.fetch()
                 if read.reference_length])

    cons_sam_file.close()

    return score


def main(args):
    '''main method.'''
    aligner = PathwayAligner(args[0], args[1], args[2], args[3],
                             args[4], args[5], args[6],
                             args[10:], args[7],
                             int(args[8]), int(args[9]))

    aligner.score_alignments()
    score_df, mismatches_df = aligner.get_results()

    score_df.to_csv('score.csv')
    mismatches_df.to_csv('mismatches.csv')


if __name__ == '__main__':
    main(sys.argv[1:])
