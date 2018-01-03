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
import sys
import uuid

from Bio import SeqIO
import pysam
from synbiochem.utils import ice_utils, seq_utils, thread_utils

import pandas as pd
from seq_genie import utils


class PathwayAligner(object):
    '''Class to align NGS data to pathways.'''

    def __init__(self, barcodes_filename, reads_filename, ice_ids_filename,
                 ice_url, ice_username, ice_password,
                 for_primer, rev_primer,
                 barcode_ice_filename=None,
                 dp_filter=100):

        # Initialise project directory:
        self.__dir_name = str(uuid.uuid4())
        os.makedirs(self.__dir_name)

        # Get pathway sequences from ICE:
        self.__ice_files = _get_ice_files(ice_url, ice_username, ice_password,
                                          ice_ids_filename,
                                          for_primer, rev_primer,
                                          self.__dir_name)

        # Get barcode-ICE dict:
        self.__barcode_ice = _get_barcode_ice(barcode_ice_filename) \
            if barcode_ice_filename else None

        # Demultiplex barcoded reads:
        reads = utils.get_reads(reads_filename)
        self.__barcode_reads = utils.bin_seqs(_get_barcodes(barcodes_filename),
                                              reads)

        # Initialise dataframes:
        self.__score_df = pd.DataFrame(columns=self.__ice_files.keys(),
                                       index=self.__barcode_reads.keys())

        self.__score_df.reindex_axis(sorted(self.__score_df.columns), axis=1)

        self.__mismatches_df = pd.DataFrame(columns=self.__ice_files.keys(),
                                            index=self.__barcode_reads.keys())

        self.__mismatches_df.reindex_axis(sorted(self.__mismatches_df.columns),
                                          axis=1)

        self.__dp_filter = dp_filter

    def score_alignments(self, num_threads=8):
        '''Score alignments.'''
        for templ_filename in self.__ice_files.values():
            utils.index(templ_filename)

        thread_pool = thread_utils.ThreadPool(num_threads)

        for barcode, reads in self.__barcode_reads.iteritems():
            reads_filename = os.path.join(self.__dir_name, barcode + '.fasta')
            SeqIO.write(reads, reads_filename, 'fasta')

            thread_pool.add_task(_score_alignment,
                                 self.__dir_name,
                                 barcode,
                                 reads_filename,
                                 self.__ice_files,
                                 self.__barcode_ice,
                                 self.__score_df,
                                 self.__mismatches_df,
                                 self.__dp_filter)

        thread_pool.wait_completion()

    def get_results(self):
        '''Get results.'''
        return self.__score_df, self.__mismatches_df


def _get_barcodes(filename):
    '''Get barcode sequences.'''
    barcodes_df = pd.read_csv(filename)
    return {row['well']: row['barcode'] for _, row in barcodes_df.iterrows()}


def _get_ice_files(url, username, password, ice_ids_filename,
                   for_primer, rev_primer, dir_name):
    '''Get ICE sequences.'''
    ice_client = ice_utils.ICEClient(url, username, password)

    with open(ice_ids_filename, 'rU') as ice_ids_file:
        ice_ids = [line.strip() for line in ice_ids_file]

    seqs = [utils.pcr(ice_client.get_ice_entry(ice_id).get_seq(),
                      for_primer, rev_primer)
            for ice_id in ice_ids]

    return {ice_id:
            seq_utils.write_fasta({ice_id: seq},
                                  os.path.join(dir_name, ice_id + '.fasta'))
            for ice_id, seq in zip(ice_ids, seqs)}


def _get_barcode_ice(barcode_ice_filename):
    '''Get barcode ice dict.'''
    barcode_ice = pd.read_csv(barcode_ice_filename,
                              dtype={'barcode': str, 'ice_id': str}) \
        if barcode_ice_filename else None

    return barcode_ice.set_index('barcode')['ice_id'].to_dict()


def _score_alignment(dir_name, barcode, reads_filename, ice_files, barcode_ice,
                     score_df, mismatches_df, dp_filter):
    '''Score an alignment.'''
    for ice_id, templ_filename in ice_files.iteritems():
        if not barcode_ice or barcode_ice[barcode] == ice_id:
            _score_barcode_ice(templ_filename, dir_name, barcode, ice_id,
                               reads_filename,
                               score_df, mismatches_df, dp_filter)

    score_df.to_csv('score.csv')
    mismatches_df.to_csv('mismatches.csv')


def _score_barcode_ice(templ_pcr_filename, dir_name, barcode, ice_id,
                       reads_filename, score_df, mismatches_df, dp_filter):
    '''Score barcode ice pair.'''
    sam_filename = os.path.join(dir_name, barcode + '_' + ice_id + '.sam')
    bam_filename = os.path.join(dir_name, barcode + '_' + ice_id + '.bam')

    # Align:
    utils.mem(templ_pcr_filename, reads_filename, sam_filename)

    # Convert sam to bam and sort:
    pysam.view(sam_filename, '-o', bam_filename, catch_stdout=False)
    pysam.sort('-o', bam_filename, bam_filename)
    os.remove(sam_filename)

    # Generate variants file:
    vcf_filename = utils.get_vcf(bam_filename, templ_pcr_filename, dp_filter)

    # Analyse variants file:
    matches, mismatches = utils.get_mismatches(vcf_filename)

    score_df[ice_id][barcode] = matches
    mismatches_df[ice_id][barcode] = mismatches


def main(args):
    '''main method.'''
    aligner = PathwayAligner(*args)
    aligner.score_alignments()
    score_df, mismatches_df = aligner.get_results()

    score_df.to_csv('score.csv')
    mismatches_df.to_csv('mismatches.csv')


if __name__ == '__main__':
    main(sys.argv[1:])
