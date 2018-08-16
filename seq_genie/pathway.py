'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
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
                 for_primer, rev_primer, min_length=0, dp_filter=0.25):
        # Initialise project directory:
        self.__dir_name = os.path.join(out_dir, str(uuid.uuid4()))
        os.makedirs(self.__dir_name)

        # Get pathway sequences from ICE:
        self.__ice_files, self.__pcr_offsets, _ = \
            _get_ice_files(ice_url, ice_username, ice_password,
                           os.path.join(in_dir, 'ice_ids.txt'),
                           for_primer, rev_primer,
                           self.__dir_name)

        # Get reads:
        self.__reads, total_reads = \
            utils.get_reads(in_dir, min_length=min_length)

        print 'Extracted %d/%d (%.1f%%) filtered reads' \
            % (len(self.__reads),
               total_reads,
               len(self.__reads) / total_reads * 100.0)

        # Initialise vcf analyser:
        self.__barcodes_df = \
            pd.read_csv(os.path.join(in_dir, 'barcodes.csv'))

        self.__barcodes_df.fillna('', inplace=True)

        self.__vcf_analyser = \
            vcf_utils.VcfAnalyser(sorted(self.__ice_files.keys()),
                                  self.__barcodes_df,
                                  dp_filter,
                                  self.__dir_name)

        self.__barcodes = self.__vcf_analyser.get_src_ids()

    def score_alignments(self, window_size, num_threads=0):
        '''Score alignments.'''
        for templ_filename, _ in self.__ice_files.values():
            utils.index(templ_filename)

        barcode_reads = utils.bin_seqs(self.__barcodes, self.__reads,
                                       window_size=window_size,
                                       num_threads=num_threads)

        print 'Extracted %d/%d (%.1f%%) barcoded reads' \
            % (len(barcode_reads),
               len(self.__reads),
               len(barcode_reads) / len(self.__reads) * 100.0)

        if num_threads:
            thread_pool = thread_utils.ThreadPool(num_threads)

            for barcodes, reads in barcode_reads.iteritems():
                reads_filename = \
                    os.path.join(utils.get_dir(self.__dir_name, barcodes),
                                 'reads.fasta')

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
                reads_filename = \
                    os.path.join(utils.get_dir(self.__dir_name, barcodes),
                                 'reads.fasta')

                SeqIO.write(reads, reads_filename, 'fasta')

                _score_alignment(self.__dir_name,
                                 barcodes,
                                 reads_filename,
                                 self.__get_ice_files(barcodes),
                                 self.__pcr_offsets,
                                 self.__vcf_analyser)

        # Update summary:
        self.__vcf_analyser.write_summary()

    def __get_ice_files(self, barcodes):
        '''Get appropriate ICE files.'''
        try:
            ice_id = \
                self._PathwayAligner__barcodes_df.loc[barcodes,
                                                      'actual_ice_id']

            if ice_id:
                return {ice_id: self.__ice_files[ice_id]}
        except KeyError:
            print 'Unexpected barcodes: ' + str(barcodes)
            return {}

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

    # Get Genbank files for subsequent data analysis:
    for ice_id in ice_ids:
        gb_filename = os.path.join(dir_name, ice_id + '.gb')
        ice_client.get_genbank(ice_id, gb_filename)

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
    barcode_dir_name = utils.get_dir(dir_name, barcodes, ice_id)
    sam_filename = os.path.join(barcode_dir_name, 'alignment.sam')
    bam_filename = os.path.join(barcode_dir_name, 'alignment.bam')

    # Align:
    utils.mem(templ_pcr_filename, reads_filename, sam_filename)

    # Convert sam to bam and sort:
    pysam.view(sam_filename, '-o', bam_filename, catch_stdout=False)
    pysam.sort('-o', bam_filename, bam_filename)
    os.remove(sam_filename)

    # Generate and analyse variants file:
    vcf_filename = utils.get_vcf(bam_filename, templ_pcr_filename,
                                 pcr_offset=0)

    vcf_analyser.analyse(vcf_filename, ice_id, barcodes)


def main(args):
    '''main method.'''
    import cProfile

    prf = cProfile.Profile()
    prf.enable()

    try:
        num_threads = int(args[-1])
    except ValueError:
        if args[-1] == 'True':
            num_threads = multiprocessing.cpu_count()
        else:
            num_threads = 0

    print 'Running pathway with %d threads' % num_threads

    aligner = PathwayAligner(*args[:-3], min_length=int(args[-3]))
    aligner.score_alignments(int(args[-2]), num_threads)

    prf.disable()

    prf.print_stats(sort='cumtime')


if __name__ == '__main__':
    main(sys.argv[1:])
