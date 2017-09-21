'''
DNA++ (c) DNA++ 2017

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-name-in-module
import os
import sys

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from pysam import AlignmentFile

from seq_genie import utils
import matplotlib.pyplot as plt


def analyse(sam_files, templ_seq):
    '''Analyse data.'''
    # analyse_dna_mut(sam_files, templ_seq)
    analyse_aa_mut(sam_files, templ_seq)


def analyse_dna_mut(sam_files, templ_seq):
    '''Analyse DNA mutations.'''
    mutations = []

    for sam_file in sam_files:
        vals = []
        mutations.append(vals)

        for read in sam_file:
            if read.alen == len(templ_seq):
                vals.append(len(templ_seq) -
                            sum([(nucl_a == nucl_b)
                                 for (nucl_a, nucl_b)
                                 in zip(read.seq[read.qstart:], templ_seq)]))

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


def analyse_aa_mut(sam_files, templ_seq):
    '''Analyse amino acid mutations.'''
    template_aa = templ_seq.translate()
    print template_aa

    all_muts = []

    for _, sam_file in enumerate(sam_files):
        muts = [[] for _ in range(len(template_aa))]
        all_muts.append(muts)

        for read in sam_file:
            read_dna = Seq(read.seq[read.qstart:read.qend])
            read_aa = read_dna.translate()

            if len(read_aa) == len(template_aa):
                for (pos, aas) in enumerate(zip(read_aa, template_aa)):
                    if aas[0] != aas[1]:
                        muts[pos].append(aas[0])


def main(args):
    '''main method.'''
    templ_filename = args[0]
    templ_seq = list(utils.parse(templ_filename))[0].seq

    sam_files = []

    for reads_filename in args[1:]:
        # Align raw file:
        sam_filename = reads_filename + '_raw.sam'
        utils.align(templ_filename, utils.parse(reads_filename),
                    out=sam_filename,
                    gap_open=12)

        # Filter indels:
        sam_filt_filename = reads_filename + '_filtered.sam'
        utils.reject_indels(sam_filename, templ_seq,
                            out_filename=sam_filt_filename)

        # Analyse:
        sam_files.append(AlignmentFile(sam_filt_filename, 'r'))

    analyse(sam_files, templ_seq)


if __name__ == '__main__':
    main(sys.argv[1:])
