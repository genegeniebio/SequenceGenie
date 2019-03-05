'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=no-name-in-module
# pylint: disable=old-style-class
# pylint: disable=superfluous-parens
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
# pylint: disable=unused-import
# pylint: disable=wrong-import-order
from collections import defaultdict, Counter
import os
import sys

from Bio import Seq
from mpl_toolkits.mplot3d import Axes3D
# from pysal.explore.inequality import gini
import pysam
from synbiochem.utils import mut_utils, seq_utils

import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import pandas as pd
from seq_genie import demultiplex, utils


INDELS_IGNORE = 0
INDELS_REJECT = 1
INDELS_REPLACE = 2


def align(templ_filename, barcodes_filename, in_dir, out_dir, min_length=1000,
          max_read_files=1e16, tolerance=6, indels=INDELS_REPLACE,
          search_len=16):
    '''Align sequence files.'''
    barcodes, _ = demultiplex.get_barcodes(barcodes_filename)

    utils.index(templ_filename)

    barcode_reads = demultiplex.demultiplex(barcodes,
                                            in_dir,
                                            min_length,
                                            max_read_files,
                                            out_dir,
                                            tolerance=tolerance,
                                            num_threads=mp.cpu_count(),
                                            search_len=search_len)

    return [_align(templ_filename, barcode_reads[key], indels)
            for key in [tuple(list(barcode) + ['all'])
                        for barcode in barcodes]]


def _align(templ_filename, reads_filename, indels):
    '''Align a single Fasta file.'''
    name, _ = os.path.splitext(reads_filename)

    sam_filename_in = name + '_raw.sam'
    align_filename = sam_filename_in

    utils.mem(templ_filename, reads_filename,
              out_filename=sam_filename_in,
              gap_open=12)

    if indels == INDELS_REJECT:
        sam_filename_out = name + '_indels_reject.sam'
        align_filename = sam_filename_out
        utils.reject_indels(sam_filename_in, templ_filename, sam_filename_out)
    elif indels == INDELS_REPLACE:
        # Filter indels:
        sam_filename_out = name + '_indels_replace.sam'
        align_filename = sam_filename_out
        utils.replace_indels(sam_filename_in, templ_filename, sam_filename_out)

    return align_filename


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


def analyse_aa_mut(sam_files, templ_aa_seq):
    '''Analyse amino acid mutations.'''
    all_muts = [[[] for _ in range(len(templ_aa_seq))]
                for _ in range(len(sam_files))]

    seqs_to_bins = defaultdict(list)

    for sam_idx, sam_filename in enumerate(sam_files):
        sam_file = pysam.AlignmentFile(sam_filename, 'r')

        for read in sam_file:
            read_muts, read_aa = \
                _analyse_aa_mut(read, templ_aa_seq)

            if read_muts is not None:
                for pos, mut in read_muts.items():
                    all_muts[sam_idx][pos].append(mut)

                seqs_to_bins[read_aa].append(sam_idx + 1)

    return all_muts, seqs_to_bins


def _analyse_aa_mut(read, template_aa):
    '''Analyse amino acid mutations in a single read.'''
    read_dna = Seq.Seq(read.seq[read.qstart:read.qend])
    read_aa = read_dna.translate()

    if len(read_aa) == len(template_aa):
        read_muts = {}

        for (pos, aas) in enumerate(list(zip(read_aa, template_aa))):
            if aas[0] != aas[1]:
                read_muts[pos] = aas[0]

        return read_muts, read_aa

    return None, None


def plot_stacked(data, filename='stacked.png'):
    '''Plots mutant counts as stacked bar chart.'''
    mut_counts = [[len(pos) for pos in mut] for mut in data]
    plt_bars = []

    for idx, datum in enumerate(mut_counts):
        bottom = [0 for _ in range(len(datum))] \
            if idx == 0 else mut_counts[idx - 1]

        plt_bars.append(plt.bar(range(len(datum)), datum, bottom=bottom))

    plt.xlabel('Residue')
    plt.ylabel('Count')
    plt.title('Mutants per residue by activity bin')
    plt.legend(plt_bars, ['Bin ' + str(idx + 1)
                          for idx in range(len(plt_bars))])

    plt.savefig(filename)


def plot_heatmap(seq_data, seq_len, filename='heatmap.png'):
    '''Plot heatmap.'''
    df = pd.DataFrame(0.0, columns=seq_utils.AA_CODES.values(),
                      index=range(1, seq_len + 1))

    for seq_datum in seq_data:
        for mut in seq_datum[0]:
            df[mut.get_mut_res()][mut.get_pos()] = seq_datum[3]

    df = df.transpose()
    df = df.sort_index()

    # Plot heatmap
    heatmap = plt.pcolor(df, cmap='Greens')
    plt.title('Mutations per residue number')
    plt.xlabel('Residue number')
    plt.ylabel('Amino acid')
    plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xticks(np.arange(0, len(df.columns), 25))
    plt.colorbar(heatmap)
    plt.savefig(filename)


def plot3d(data, filename='3d.png'):
    '''Plots mutant counts in 3d.'''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for z, bin_data in enumerate(data):
        xs = np.arange(len(bin_data))
        ys = [len(set(muts)) for muts in bin_data]
        ax.bar(xs, ys, zs=z, zdir='y', alpha=0.8)

    ax.set_xlabel('Residue')
    ax.set_ylabel('Bin number')
    ax.set_zlabel('Number of mutations')

    plt.savefig(filename)


def get_gini(muts):
    '''Get Gini value for each position.'''
    from pysal.explore.inequality import gini

    scores = []
    mut_probs = MutProbs()

    for idx, act_bin in enumerate(muts[1]):
        for mut_res in act_bin:
            if mut_res != '*':
                scores.append(
                    (1 - mut_probs.get_mut_prob(muts[0], mut_res)) * 10 ** idx)

    return gini.Gini(scores).g if scores else 0


def main(args):
    '''main method.'''
    templ_filename = args[0]
    templ_seq = utils.get_seq(templ_filename)
    templ_aa_seq = templ_seq.translate()
    out_dir = args[3]

    # Align:
    sam_files = align(templ_filename, args[1], args[2], out_dir)

    # Analyse:
    muts, seqs_to_bins = analyse_aa_mut(sam_files, templ_aa_seq)

    for seq, bins in seqs_to_bins.items():
        bin_counts = {count + 1: 0 for count in range(len(sam_files))}
        bin_counts.update(dict(Counter(bins)))

        seqs_to_bins[seq] = \
            [mut_utils.get_mutations(templ_aa_seq, seq)] + \
            list(bin_counts.values()) + \
            [len(bins),
             np.mean(bins),
             np.std(bins)]

    columns = ['mutations'] + list(bin_counts.keys()) + \
        ['occurances', 'mean_bin', 'stddev_bin']
    seq_to_bins_df = \
        pd.DataFrame(list(seqs_to_bins.values()),
                     columns=columns)

    seq_to_bins_df.sort_values(['mean_bin', 'stddev_bin'],
                               ascending=[False, True],
                               inplace=True)

    seq_to_bins_df.to_csv(os.path.join(out_dir, 'seqs_to_bins.csv'),
                          index=False)

    with open(os.path.join(out_dir, 'mutations_by_pos.txt'), 'w') as outfile:
        for idx, vals in enumerate(zip(templ_aa_seq, zip(*muts))):
            outfile.write('\t'.join([str(val)
                                     for val in [idx + 1, vals[0],
                                                 get_gini(vals)] +
                                     list(vals[1])]) + '\n')


class MutProbs():
    '''Class to represent mutation probabilities.'''

    def __init__(self):
        directory = os.path.dirname(os.path.realpath(__file__))
        filename = os.path.join(directory, 'blosum62.qij')

        with open(filename) as fle:
            row = 0

            for line in fle:
                if line.startswith('#'):
                    continue

                tokens = line.split()

                if tokens[0][0].isalpha():
                    self.__mut_probs = pd.DataFrame(index=tokens,
                                                    columns=tokens)
                else:
                    self.__mut_probs.iloc[row, :len(tokens)] = \
                        [float(val) for val in tokens]

                    row = row + 1

        # Fill top triangle:
        for idx, column in enumerate(self.__mut_probs):
            self.__mut_probs.iloc[idx] = self.__mut_probs[column]

        # Normalise probabilities for each row / column / residue:
        self.__mut_probs = self.__mut_probs / self.__mut_probs.sum()

    def get_mut_prob(self, wt_res, mut_res):
        '''Gets mutation probabilities.'''
        return self.__mut_probs[wt_res][mut_res]


if __name__ == '__main__':
    main(sys.argv[1:])
