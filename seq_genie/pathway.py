'''
SequenceGenie (c) GeneGenie Bioinformatics Ltd. 2017

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes
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

    def __init__(self, out_dir,
                 barcodes_filename, reads_filename, ice_ids_filename,
                 ice_url, ice_username, ice_password,
                 for_primer, rev_primer, dp_filter=100):

        # Initialise project directory:
        self.__dir_name = os.path.join(out_dir, str(uuid.uuid4()))
        os.makedirs(self.__dir_name)

        # Get pathway sequences from ICE:
        self.__ice_files = _get_ice_files(ice_url, ice_username, ice_password,
                                          ice_ids_filename,
                                          for_primer, rev_primer,
                                          self.__dir_name)

        # Demultiplex barcoded reads:
        reads = utils.get_reads(reads_filename)
        self.__summary_df = pd.read_csv(barcodes_filename)
        self.__barcode_reads = utils.bin_seqs(
            self.__summary_df['barcode'].tolist(), reads)

        self.__summary_df.set_index('barcode', inplace=True)

        # Initialise dataframes:
        columns = sorted(self.__ice_files.keys())
        self.__identity_df = pd.DataFrame(columns=columns,
                                          index=self.__summary_df.index,
                                          dtype='float')

        self.__mutations_df = pd.DataFrame(columns=columns,
                                           index=self.__summary_df.index)

        self.__deletions_df = pd.DataFrame(columns=columns,
                                           index=self.__summary_df.index)

        self.__dp_filter = dp_filter

    def score_alignments(self, num_threads=8):
        '''Score alignments.'''
        for templ_filename, _ in self.__ice_files.values():
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
                                 self.__identity_df,
                                 self.__mutations_df,
                                 self.__deletions_df,
                                 self.__dp_filter)

        thread_pool.wait_completion()

        # Update summary:
        self.__identity_df.fillna(0, inplace=True)
        self.__summary_df['ice_id'] = self.__identity_df.idxmax(axis=1)
        self.__summary_df['identity'] = self.__identity_df.max(axis=1)
        self.__summary_df['mutations'] = \
            self.__mutations_df.lookup(self.__mutations_df.index,
                                       self.__summary_df['ice_id'])
        self.__summary_df['deletions'] = \
            self.__mutations_df.lookup(self.__deletions_df.index,
                                       self.__summary_df['ice_id'])

        self.__summary_df.to_csv(os.path.join(self.__dir_name, 'summary.csv'))

    def get_results(self):
        '''Get results.'''
        return self.__summary_df, self.__identity_df, self.__mutations_df


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
            (seq_utils.write_fasta({ice_id: seq},
                                   os.path.join(dir_name, ice_id + '.fasta')),
             len(seq))
            for ice_id, seq in zip(ice_ids, seqs)}


def _get_barcode_ice(barcode_ice_filename):
    '''Get barcode ice dict.'''
    barcode_ice = pd.read_csv(barcode_ice_filename,
                              dtype={'barcode': str, 'ice_id': str}) \
        if barcode_ice_filename else None

    return barcode_ice.set_index('barcode')['ice_id'].to_dict()


def _score_alignment(dir_name, barcode, reads_filename, ice_files,
                     identity_df, mutations_df, deletions_df, dp_filter):
    '''Score an alignment.'''
    for ice_id, (templ_filename, templ_len) in ice_files.iteritems():
        _score_barcode_ice(templ_filename, templ_len, dir_name, barcode,
                           ice_id, reads_filename,
                           identity_df, mutations_df, deletions_df,
                           dp_filter)

        identity_df.to_csv(os.path.join(dir_name, 'identity.csv'))
        mutations_df.to_csv(os.path.join(dir_name, 'mutations.csv'))


def _score_barcode_ice(templ_pcr_filename, templ_len, dir_name, barcode,
                       ice_id, reads_filename,
                       identity_df, mutations_df, deletions_df,
                       dp_filter):
    '''Score barcode ice pair.'''
    sam_filename = os.path.join(dir_name, barcode + '_' + ice_id + '.sam')
    bam_filename = os.path.join(dir_name, barcode + '_' + ice_id + '.bam')

    # Align:
    utils.mem(templ_pcr_filename, reads_filename, sam_filename)

    # Convert sam to bam and sort:
    pysam.view(sam_filename, '-o', bam_filename, catch_stdout=False)
    pysam.sort('-o', bam_filename, bam_filename)
    os.remove(sam_filename)

    # Generate and analyse variants file:
    num_matches, mutations, deletions = \
        utils.analyse_vcf(utils.get_vcf(bam_filename,
                                        templ_pcr_filename,
                                        dp_filter))

    identity_df[ice_id][barcode] = num_matches / float(templ_len)
    mutations_df[ice_id][barcode] = mutations
    deletions_df[ice_id][barcode] = deletions


def main(args):
    '''main method.'''
    aligner = PathwayAligner(*args)
    aligner.score_alignments()


if __name__ == '__main__':
    main(sys.argv[1:])
