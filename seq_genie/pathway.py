'''
SequenceGenie (c) GeneGenie Bioinformatics Ltd. 2017

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=too-many-arguments
from __future__ import division

import os
import random
import sys
import threading
import uuid

from Bio import SeqIO
import pysam

import pandas as pd
from seq_genie import utils
from synbiochem.utils import ice_utils, seq_utils


def identify(barcodes_filename, reads_filename,
             ice_url, ice_username, ice_password, ice_ids,
             min_length=128, max_seqs=128000):
    '''Identify plamids from sequences.'''
    dir_name = str(uuid.uuid4())

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    barcodes = get_barcodes(barcodes_filename)
    reads = utils.get_reads(reads_filename, min_length)

    print len(reads)

    barcode_reads = utils.bin_seqs(barcodes,
                                   random.sample(reads,
                                                 min(len(reads), max_seqs)))

    print [len(val) for val in barcode_reads.values()]
    print sum([len(val) for val in barcode_reads.values()])

    ice_files = get_ice_files(ice_url, ice_username, ice_password, ice_ids,
                              dir_name)

    return score_alignments(ice_files, barcode_reads, dir_name)


def get_barcodes(filename):
    '''Get barcode sequences.'''
    barcodes_df = pd.read_csv(filename)
    return {row['well']: row['barcode'] for _, row in barcodes_df.iterrows()}


def get_ice_files(url, username, password, ice_ids, dir_name):
    '''Get ICE sequences.'''
    ice_client = ice_utils.ICEClient(url, username, password)

    return {ice_id:
            seq_utils.write_fasta({ice_id:
                                   ice_client.get_ice_entry(ice_id).get_seq()},
                                  os.path.join(dir_name, ice_id + '.fasta'))
            for ice_id in ice_ids}


def score_alignments(ice_files, barcode_reads, dir_name, num_threads=8):
    '''Score alignments.'''
    df = pd.DataFrame(columns=ice_files.keys(), index=barcode_reads.keys())

    for templ_filename in ice_files.values():
        utils.index(templ_filename)

    batches = [barcode_reads.keys()[pos:pos + num_threads]
               for pos in xrange(0, len(barcode_reads.keys()), num_threads)]

    for batched_barcode_reads in batches:
        threads = []

        for barcode in batched_barcode_reads:
            t = threading.Thread(target=_score_alignment,
                                 args=(dir_name, barcode,
                                       barcode_reads[barcode],
                                       ice_files, df))
            t.daemon = True
            t.start()
            threads.append(t)

        for t in threads:
            t.join()

    return df


def _score_alignment(dir_name, barcode, reads, ice_files, df):
    '''Score an alignment.'''
    reads_filename = os.path.join(dir_name, barcode + '.fasta')
    SeqIO.write(reads, reads_filename, 'fasta')

    for ice_id, templ_filename in ice_files.iteritems():
        sam_filename = os.path.join(dir_name,
                                    barcode + '_' + ice_id + '.sam')

        cons_sam_filename = os.path.join(
            dir_name, barcode + '_' + ice_id + '_cons.sam')

        # Align:
        utils.align(templ_filename, reads_filename, sam_filename)

        # Generate then align consensus:
        fasta_filename = utils.get_consensus(sam_filename, templ_filename)
        utils.align(templ_filename, fasta_filename, cons_sam_filename)

        df[ice_id][barcode] = _score(cons_sam_filename)


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
    df = identify(args[0], args[1], args[2], args[3], args[4], args[7:],
                  int(args[5]), int(args[6]))

    df.to_csv('out.csv')


if __name__ == '__main__':
    main(sys.argv[1:])
