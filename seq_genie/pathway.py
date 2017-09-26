'''
DNA++ (c) DNA++ 2017

All rights reserved.

@author: neilswainston
'''
import sys

from synbiochem.utils import ice_utils, seq_utils

from seq_genie import utils
import matplotlib.pyplot as plt


def blast(ice_client, seq):
    '''Does BLAST search against ICE.'''
    results = ice_client.do_blast(seq)

    if results['results']:
        res = results['results'][0]
        print '\t'.join(str(val)
                        for val in [res['entryInfo']['partId'], res['eValue']])


def analyse(reads):
    '''Analyse records.'''
    for read in reads:
        print '\t'.join(str(val)
                        for val in [read.id, len(read.seq), read.seq])

    plt.hist([len(read.seq) for read in reads], edgecolor='black')
    plt.title('Read lengths')
    plt.xlabel('Read length')
    plt.ylabel('Frequency')
    plt.show()


def main(args):
    '''main method.'''
    reads = utils.get_reads([args[0]], int(args[1]))
    ice_client = ice_utils.ICEClient(args[2], args[3], args[4])

    for ice_id in args[5:]:
        id_seqs = {ice_id: ice_client.get_ice_entry(ice_id).get_seq()}
        templ_filename = ice_id + '.fasta'
        seq_utils.write_fasta(id_seqs, templ_filename)
        utils.align(templ_filename, reads, out=ice_id + '.sam')


if __name__ == '__main__':
    main(sys.argv[1:])
