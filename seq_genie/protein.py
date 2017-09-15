'''
DNA++ (c) DNA++ 2017

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-name-in-module
import os
import sys

from pysam import AlignmentFile

from seq_genie import utils
import matplotlib.pyplot as plt


def analyse(sam_files, templ_seq):
    '''Analyse data.'''
    mutations = []

    for sam_file in sam_files:
        vals = []
        mutations.append(vals)

        for read in sam_file:
            if read.alen == len(templ_seq):
                vals.append(len(templ_seq) -
                            sum([(nucl_a == nucl_b)
                                 for (nucl_a, nucl_b)
                                 in zip(read.seq, templ_seq)]))

    bins = max([val for vals in mutations for val in vals]) - \
        min([val for vals in mutations for val in vals])

    plt.hist(mutations, bins=bins, edgecolor='black',
             label=[os.path.split(sam_file.filename)[1]
                    for sam_file in sam_files])
    plt.title('Mutation counts')
    plt.xlabel('Number of mutations')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()


def main(args):
    '''main method.'''
    templ_filename = args[0]
    templ_seq = list(utils.parse(templ_filename))[0].seq

    sam_files = []

    for reads_filename in args[1:]:
        # Align raw file:
        sam_filename = reads_filename + '.sam'
        utils.align(templ_filename, utils.parse(reads_filename),
                    out=sam_filename)

        # Strip indels:
        strip_filename = utils.strip_indels(sam_filename, templ_seq)

        # Re-align:
        utils.align(templ_filename, utils.parse(strip_filename),
                    out=sam_filename)

        # Analyse:
        sam_files.append(AlignmentFile(sam_filename, 'r'))

    analyse(sam_files, templ_seq)


if __name__ == '__main__':
    main(sys.argv[1:])
