'''
DNA++ (c) DNA++ 2017

All rights reserved.

@author: neilswainston
'''
import os
import sys

from Bio import Seq, SeqIO, SeqRecord

from seq_genie import utils


def main(args):
    '''main method.'''
    args_idx = len(args)

    for args_idx, arg in enumerate(args):
        if '.fasta' in arg:
            break

    min_templ_len = float('inf')

    for templ_filename in args[args_idx:]:
        with open(templ_filename, 'r') as templ:
            max_templ_len = min(min_templ_len, len(SeqIO.read(templ,
                                                              format='fasta')))

    reads = utils.get_reads([args[0]], min_length=max_templ_len)

    barcode_seqs = \
        utils.bin_seqs({str(idx): barcode
                        for idx, barcode in enumerate(args[1:args_idx])},
                       {read.id: str(read.seq)
                        for read in reads},
                       evalue=1e-3)

    print [(key, len(values)) for key, values in barcode_seqs.iteritems()]

    for key, values in barcode_seqs.iteritems():
        reads = [SeqRecord.SeqRecord(Seq.Seq(seq), id=seq_id)
                 for seq_id, seq in values.iteritems()]

        # with open(key + '.fasta', 'w') as out:
        #    SeqIO.write(reads, out, 'fasta')

        for templ_filename in args[args_idx:]:
            _, tail = os.path.split(templ_filename)
            sam_filename = key + '_' + tail + '.sam'
            utils.align(templ_filename, reads, sam_filename)
            # utils.mpileup(sam_filename, templ_filename,
            #               key + '_' + tail + '.txt')


if __name__ == '__main__':
    main(sys.argv[1:])
