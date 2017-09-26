'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import sys

from synbiochem.utils import seq_utils


def get_seqs(num_barcodes, length, max_repeat_nuc, evalue=1):
    '''Get barcodes.'''
    barcodes = {}

    while len(barcodes) < num_barcodes:
        barcode = seq_utils.get_random_dna(length, max_repeat_nuc)

        if barcodes and not do_blast(barcodes, {'query': barcode}, evalue):
            continue

        barcodes[len(barcodes)] = barcode
        print str(len(barcodes)) + '\t' + barcode

    return barcodes.values()


def do_blast(subjects, queries, evalue):
    '''Runs BLAST, filtering results.'''
    for result in seq_utils.do_blast(subjects, queries, evalue=evalue,
                                     word_size=4):
        if result.alignments:
            return False

    return True


def main(args):
    '''main method.'''
    get_seqs(int(args[0]), int(args[1]), int(args[2]), float(args[3]))


if __name__ == '__main__':
    main(sys.argv[1:])
