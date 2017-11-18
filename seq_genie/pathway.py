'''
SequenceGenie (c) GeneGenie Bioinformatics Ltd.

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-member
import os
import sys
import uuid

from Bio import Seq, SeqRecord

import pandas as pd
from seq_genie import utils
from synbiochem.utils import ice_utils, seq_utils


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


def main(args):
    '''main method.'''
    dir_name = str(uuid.uuid4())

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    barcodes = get_barcodes(args[0])
    reads = utils.get_reads(args[1])

    barcode_seqs = utils.bin_seqs(barcodes,
                                  {read.id: str(read.seq)
                                   for read in reads[:128]},
                                  evalue=1e-3)

    ice_files = get_ice_files(args[2], args[3], args[4], args[5:], dir_name)

    for key, values in barcode_seqs.iteritems():
        reads = [SeqRecord.SeqRecord(Seq.Seq(seq), id=seq_id)
                 for seq_id, seq in values.iteritems()]

        for ice_id, filename in ice_files.iteritems():
            sam_filename = os.path.join(dir_name, key + '_' + ice_id + '.sam')
            utils.align(filename, reads, sam_filename)

            import pysam
            sam_file = pysam.AlignmentFile(sam_filename, 'r')

            for read in sam_file.fetch():
                print read.query_alignment_length

            sam_file.close()


if __name__ == '__main__':
    main(sys.argv[1:])
