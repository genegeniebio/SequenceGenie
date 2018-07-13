'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes
# pylint: disable=too-many-locals
# pylint: disable=wrong-import-order
from __future__ import division

import multiprocessing
import os
import sys
import uuid

from Bio import SeqIO
import pysam
from synbiochem.utils import ice_utils, seq_utils, thread_utils

import pandas as pd
from seq_genie import utils, vcf_utils


class PathwayAligner(object):
    '''Class to align NGS data to pathways.'''

    def __init__(self, out_dir, in_dir,
                 ice_url, ice_username, ice_password,
                 for_primer, rev_primer, dp_filter=0.25):
        # Initialise project directory:
        self.__dir_name = os.path.join(out_dir, str(uuid.uuid4()))
        os.makedirs(self.__dir_name)

        # Get pathway sequences from ICE:
        self.__ice_files, self.__pcr_offsets, ice_lengths = \
            _get_ice_files(ice_url, ice_username, ice_password,
                           os.path.join(in_dir, 'ice_ids.txt'),
                           for_primer, rev_primer,
                           self.__dir_name)

        # Get reads:
        self.__reads = utils.get_reads(in_dir, min_length=min(ice_lengths))

        print 'Extracted %d filtered reads' % len(self.__reads)

        # Initialise vcf analyser:
        self.__barcodes_df = pd.read_csv(os.path.join(in_dir, 'barcodes.csv'))
        self.__barcodes_df.dropna(inplace=True)

        self.__vcf_analyser = \
            vcf_utils.VcfAnalyser(sorted(self.__ice_files.keys()),
                                  self.__barcodes_df,
                                  dp_filter,
                                  self.__dir_name)

        self.__barcodes = self.__vcf_analyser.get_src_ids()

    def score_alignments(self, score_threshold, num_threads=0):
        '''Score alignments.'''
        for templ_filename, _ in self.__ice_files.values():
            utils.index(templ_filename)

        barcode_reads = utils.bin_seqs(self.__barcodes, self.__reads,
                                       score_threshold=score_threshold,
                                       num_threads=num_threads)

        if num_threads:
            thread_pool = thread_utils.ThreadPool(num_threads)

            for barcodes, reads in barcode_reads.iteritems():
                reads_filename = os.path.join(self.__dir_name,
                                              '_'.join(barcodes) + '.fasta')
                SeqIO.write(reads, reads_filename, 'fasta')

                thread_pool.add_task(_score_alignment,
                                     self.__dir_name,
                                     barcodes,
                                     reads_filename,
                                     self.__get_ice_files(barcodes),
                                     self.__pcr_offsets,
                                     self.__vcf_analyser)

            thread_pool.wait_completion()
        else:
            for barcodes, reads in barcode_reads.iteritems():
                reads_filename = os.path.join(self.__dir_name,
                                              '_'.join(barcodes) + '.fasta')
                SeqIO.write(reads, reads_filename, 'fasta')

                _score_alignment(self.__dir_name,
                                 barcodes,
                                 reads_filename,
                                 self.__get_ice_files(barcodes),
                                 self.__pcr_offsets,
                                 self.__vcf_analyser)

        # Update summary:
        self.__vcf_analyser.write_summary()

    def get_results(self):
        '''Get results.'''
        return self.__summary_df, self.__identity_df, self.__mutations_df, \
            self.__indels_df, self.__deletions_df

    def __get_ice_files(self, barcodes):
        '''Get appropriate ICE files.'''
        ice_id = \
            self._PathwayAligner__barcodes_df.loc[barcodes, 'actual_ice_id']

        if ice_id:
            return {ice_id: self.__ice_files[ice_id]}

        return self.__ice_files


def _get_ice_files(url, username, password, ice_ids_filename,
                   for_primer, rev_primer, dir_name):
    '''Get ICE sequences.'''
    ice_client = ice_utils.ICEClient(url, username, password)

    with open(ice_ids_filename, 'rU') as ice_ids_file:
        ice_ids = [line.strip() for line in ice_ids_file]

    seqs_offsets = [utils.pcr(ice_client.get_ice_entry(ice_id).get_seq(),
                              for_primer, rev_primer)
                    for ice_id in ice_ids]

    seqs, offsets = zip(*seqs_offsets)

    ice_files = {ice_id:
                 (seq_utils.write_fasta({ice_id: seq},
                                        os.path.join(dir_name,
                                                     ice_id + '.fasta')),
                  len(seq))
                 for ice_id, seq in zip(ice_ids, seqs)}

    pcr_offsets = {ice_id: offset for ice_id, offset in zip(ice_ids, offsets)}

    return ice_files, pcr_offsets, [len(seq) for seq in seqs]


def _get_barcode_ice(barcode_ice_filename):
    '''Get barcode ice dict.'''
    barcode_ice = pd.read_csv(barcode_ice_filename,
                              dtype={'barcode': str, 'ice_id': str}) \
        if barcode_ice_filename else None

    return barcode_ice.set_index('barcode')['ice_id'].to_dict()


def _score_alignment(dir_name, barcodes, reads_filename,
                     ice_files, pcr_offsets,
                     vcf_analyser):
    '''Score an alignment.'''
    for ice_id, (templ_filename, _) in ice_files.iteritems():
        _score_barcodes_ice(templ_filename, dir_name, barcodes,
                            ice_id, pcr_offsets[ice_id], reads_filename,
                            vcf_analyser)


def _score_barcodes_ice(templ_pcr_filename, dir_name, barcodes,
                        ice_id, pcr_offset, reads_filename,
                        vcf_analyser):
    '''Score barcodes ice pair.'''
    barcode_ice = '_'.join(list(barcodes) + [ice_id])
    sam_filename = os.path.join(dir_name, barcode_ice + '.sam')
    bam_filename = os.path.join(dir_name, barcode_ice + '.bam')

    # Align:
    utils.mem(templ_pcr_filename, reads_filename, sam_filename)

    # Convert sam to bam and sort:
    pysam.view(sam_filename, '-o', bam_filename, catch_stdout=False)
    pysam.sort('-o', bam_filename, bam_filename)
    os.remove(sam_filename)

    # Generate and analyse variants file:
    vcf_filename = utils.get_vcf(bam_filename, templ_pcr_filename, pcr_offset)
    vcf_analyser.analyse(vcf_filename, ice_id, barcodes)


def main(args):
    '''main method.'''
    try:
        num_threads = int(args[-1])
    except ValueError:
        if args[-1] == 'True':
            num_threads = multiprocessing.cpu_count()
        else:
            num_threads = 0

    print 'Running pathway with %d threads' % num_threads

    aligner = PathwayAligner(*args[:-2])
    aligner.score_alignments(float(args[-2]), num_threads)


if __name__ == '__main__':
    main(sys.argv[1:])
