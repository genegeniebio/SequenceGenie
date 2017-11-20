'''
SequenceGenie (c) GeneGenie Bioinformatics Ltd. 2017

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
from __future__ import division

import os
import sys
import uuid

from Bio import Seq, SeqRecord
import pysam

import pandas as pd
from seq_genie import utils
from synbiochem.utils import ice_utils, seq_utils


def identify(barcodes_filename, reads_filename,
             ice_url, ice_username, ice_password, ice_ids):
    '''Identify plamids from sequences.'''
    dir_name = str(uuid.uuid4())

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    barcodes = get_barcodes(barcodes_filename)
    reads = utils.get_reads(reads_filename)

    barcode_seqs = utils.bin_seqs(barcodes,
                                  {read.id: str(read.seq)
                                   for read in reads[:128]},
                                  evalue=1e-3)

    ice_files = get_ice_files(ice_url, ice_username, ice_password, ice_ids,
                              dir_name)

    return score_alignments(ice_files, barcode_seqs, dir_name)


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


def score_alignments(ice_files, barcode_seqs, dir_name):
    '''Score alignments.'''
    df = pd.DataFrame(columns=ice_files.keys(), index=barcode_seqs.keys())

    for barcode, seqs in barcode_seqs.iteritems():
        reads = [SeqRecord.SeqRecord(Seq.Seq(seq), id=seq_id)
                 for seq_id, seq in seqs.iteritems()]

        scores = []

        for ice_id, filename in ice_files.iteritems():
            sam_filename = os.path.join(dir_name,
                                        barcode + '_' + ice_id + '.sam')
            utils.align(filename, reads, sam_filename)

            sam_file = pysam.AlignmentFile(sam_filename, 'r')

            for read in sam_file.fetch():
                if read.reference_length:
                    scores.append(read.reference_length / read.query_length)

            if scores:
                df[ice_id][barcode] = sum(scores) / len(scores)

            sam_file.close()

    return df


def main(args):
    '''main method.'''
    df = identify(args[0], args[1], args[2], args[3], args[4], args[5:])
    df.to_csv('out.csv')


if __name__ == '__main__':
    main(sys.argv[1:])
