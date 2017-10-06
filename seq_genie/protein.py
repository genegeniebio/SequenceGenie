'''
DNA++ (c) DNA++ 2017

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=no-name-in-module
# pylint: disable=relative-import
from collections import defaultdict
from operator import itemgetter
import os
import sys

from Bio.Seq import Seq
from pysam import AlignmentFile
from synbiochem.utils import mut_utils

from seq_genie import utils
import matplotlib.pyplot as plt
import numpy as np


def align(templ_filename, reads_dirs):
    '''Align sequence files.'''
    align_files = []
    templ_seq = list(utils.parse(templ_filename))[0].seq

    for reads_dir in reads_dirs:
        # Align raw file:
        sam_filename = reads_dir + '_raw.sam'
        utils.align(templ_filename, utils.get_reads(reads_dir),
                    out=sam_filename,
                    gap_open=12)

        # Filter indels:
        sam_filt_flename = reads_dir + '_filtered.sam'
        utils.reject_indels(sam_filename, templ_seq,
                            out_filename=sam_filt_flename)

        align_files.append(AlignmentFile(sam_filt_flename, 'r'))

    return align_files


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

    muts = [[[] for _ in range(len(template_aa))]
            for _ in range(len(sam_files))]

    seqs_to_bins = defaultdict(list)

    for sam_idx, sam_file in enumerate(sam_files):
        for read in sam_file:
            mut = analyse_aa_mut_read(read, template_aa)

            if mut:
                muts[sam_idx][mut[0]].append(mut[1])
                seqs_to_bins[mut[2]].append(sam_idx)

    return muts, seqs_to_bins, template_aa


def analyse_aa_mut_read(read, template_aa):
    '''Analyse amino acid mutations in a single read.'''
    read_dna = Seq(read.seq[read.qstart:read.qend])
    read_aa = read_dna.translate()

    if len(read_aa) == len(template_aa):
        read_muts = {}

        for (pos, aas) in enumerate(zip(read_aa, template_aa)):
            if aas[0] != aas[1]:
                read_muts[pos] = aas[0]

        if len(read_muts) == 1:
            # pos = read_muts.keys()[0]
            # print read_dna[pos * 3:(pos + 1) * 3]

            return (read_muts.keys()[0], read_muts.values()[0], read_aa)

    return None


def plot(data):
    '''Plots mutant counts.'''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for z, bin_data in enumerate(data):
        xs = np.arange(len(bin_data))
        ys = [len(set(muts)) for muts in bin_data]
        ax.bar(xs, ys, zs=z, zdir='y', alpha=0.8)

    ax.set_xlabel('Residue')
    ax.set_ylabel('Bin number')
    ax.set_zlabel('Number of mutations')

    plt.show()


def main(args):
    '''main method.'''
    templ_filename = args[0]
    templ_seq = list(utils.parse(templ_filename))[0].seq
    templ_aa_seq = templ_seq.translate()

    sam_files = align(templ_filename, args[1:])

    # Analyse:
    muts, seqs_to_bins, template_aa = analyse_aa_mut(sam_files, templ_seq)

    for seq, bins in seqs_to_bins.iteritems():
        seqs_to_bins[seq] = [seq,
                             mut_utils.get_mutations(template_aa, seq),
                             bins,
                             len(bins),
                             np.mean(bins),
                             np.std(bins)]

    with open('seqs_to_bins.txt', 'w') as outfile:
        for vals in sorted(seqs_to_bins.values(), key=itemgetter(1)):
            outfile.write('\t'.join([str(val) for val in vals]) + '\n')

    with open('mutations_by_pos.txt', 'w') as outfile:
        for idx, vals in enumerate(zip(templ_aa_seq, zip(*muts))):
            outfile.write('\t'.join([str(val)
                                     for val in [idx + 1, vals[0]] +
                                     list(vals[1])]) + '\n')

    # plot(muts)


if __name__ == '__main__':
    main(sys.argv[1:])
