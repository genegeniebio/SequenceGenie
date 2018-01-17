'''
SequenceGenie (c) GeneGenie Bioinformatics Ltd. 2017

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-name-in-module
# pylint: disable=too-many-arguments
from collections import defaultdict
import os
from os.path import splitext
import subprocess

from Bio import Seq, SeqIO, SeqRecord
from fuzzywuzzy import fuzz
from pysam import Samfile
from synbiochem.utils import io_utils, thread_utils

import pandas as pd


def get_reads(reads_filename, min_length=0):
    '''Gets reads.'''
    reads = []

    if os.path.isdir(reads_filename):
        for dirpath, _, filenames in os.walk(os.path.abspath(reads_filename)):
            for filename in filenames:
                filename = os.path.join(dirpath, filename)
                _get_reads(filename, min_length, reads)
    else:
        _get_reads(reads_filename, min_length, reads)

    return reads


def bin_seqs(barcodes, sequences, score_threshold=90, search_len=256,
             num_threads=8):
    '''Bin sequences according to barcodes.'''
    barcode_seqs = defaultdict(list)

    max_barcode_len = max([len(barcode) for barcode in barcodes])

    if barcodes:
        thread_pool = thread_utils.ThreadPool(num_threads)

        for seq in sequences:
            thread_pool.add_task(_bin_seq, seq, max_barcode_len, search_len,
                                 score_threshold, barcodes, barcode_seqs)

        thread_pool.wait_completion()
    else:
        barcode_seqs['undefined'].extend(sequences)

    return barcode_seqs


def index(filename):
    '''Index file.'''
    subprocess.call(['bwa', 'index', filename])


def mem(templ_filename, reads_filename, out_filename=None,
        readtype='ont2d', gap_open=6):
    '''Runs BWA MEM.'''
    out_file = io_utils.get_filename(out_filename)

    with open(out_file, 'w') as out:
        subprocess.call(['bwa', 'mem',
                         '-x', readtype,
                         '-O', str(gap_open),
                         templ_filename, reads_filename],
                        stdout=out)

    return out_file


def get_vcf(bam_filename, templ_filename, dp_filter):
    '''Generates a vcf file.'''
    vcf_filename = bam_filename + '.vcf'

    proc1 = subprocess.Popen(['samtools',
                              'mpileup',
                              '-uvf',
                              templ_filename,
                              '-t', 'DP',
                              bam_filename],
                             stdout=subprocess.PIPE)

    subprocess.call(['bcftools',
                     'filter',
                     '-i', 'FMT/DP>' + str(dp_filter),
                     '-o', vcf_filename],
                    stdin=proc1.stdout)

    proc1.stdout.close()

    return vcf_filename


def sort(in_filename, out_filename):
    '''Custom sorts SAM file.'''
    sam_file = Samfile(in_filename, 'r')
    out_file = Samfile(out_filename, 'wh',
                       template=sam_file,
                       header=sam_file.header)

    for read in sorted([read for read in sam_file],
                       key=lambda x: (-x.query_length,
                                      x.reference_start)):
        out_file.write(read)

    out_file.close()

    return out_filename


def pcr(seq, forward_primer, reverse_primer):
    '''Apply in silico PCR.'''
    for_primer_pos = seq.find(forward_primer.upper())

    rev_primer_pos = \
        seq.find(str(Seq.Seq(reverse_primer).reverse_complement().upper()))

    if for_primer_pos > -1 and rev_primer_pos > -1:
        seq = seq[for_primer_pos:] + \
            seq[:rev_primer_pos + len(reverse_primer)]
    elif for_primer_pos > -1:
        seq = seq[for_primer_pos:]
    elif rev_primer_pos > -1:
        seq = seq[:rev_primer_pos + len(reverse_primer)]

    return seq


def analyse_vcf(vcf_filename):
    '''Analyse vcf file, returning number of matches, mutations and indels.'''
    num_matches = 0
    mutations = []
    deletions = []

    df = _vcf_to_df(vcf_filename)

    for pos in range(1, df.index[-1] + 1):
        try:
            row = df.loc[pos]
            alleles = [row['REF']] + row['ALT'].split(',')

            # PL values are ordered according to following:
            idx_allele = {(x**2 + 3 * x) / 2: allele
                          for x, allele in enumerate(alleles)}

            # Extract PL values and order to find most-likely base:
            pl_index = row['FORMAT'].split(':').index('PL')
            pls = row['DATA'].split(':')[pl_index].split(',')

            allele_pls = [[allele, int(pls[idx])]
                          for idx, allele in idx_allele.iteritems()]
            allele_pls.sort(key=lambda x: x[1])

            # Compare most-likely base to reference:
            if row['REF'] != allele_pls[0][0]:
                mutations.append(row['REF'] + str(pos) + allele_pls[0][0])
            else:
                num_matches += 1
        except KeyError:
            deletions.append(pos)

    return num_matches, mutations, deletions


def reject_indels(sam_filename, templ_seq, out_filename=None):
    '''Rejects indels.'''
    out_filename = io_utils.get_filename(out_filename)

    sam_file = Samfile(sam_filename, 'r')
    out_file = Samfile(out_filename, 'wh',
                       template=sam_file,
                       header=sam_file.header)

    for read in sam_file:
        if read.cigarstring and str(len(templ_seq)) + 'M' in read.cigarstring:
            out_file.write(read)

    out_file.close()

    return out_filename


def replace_indels(sam_filename, templ_seq, out_filename=None):
    '''Replace indels, replacing them with wildtype.'''
    out_filename = io_utils.get_filename(out_filename)

    with open(out_filename, 'w') as fle:
        SeqIO.write(_replace_indels(sam_filename, templ_seq), fle, 'fasta')

    return out_filename


def _get_reads(filename, min_length, reads):
    '''Gets reads.'''
    _, ext = splitext(filename)

    try:
        with open(filename, 'rU') as fle:
            reads.extend([record
                          for record in SeqIO.parse(fle, ext[1:])
                          if len(record.seq) > min_length])
    except (IOError, ValueError), err:
        print err


def _bin_seq(seq, max_barcode_len, search_len, score_threshold, barcodes,
             barcode_seqs):
    '''Bin an individual sequence.'''
    trim_seq = str(seq.seq[:max_barcode_len + search_len])
    max_score = score_threshold
    selected_barcode = None

    for barcode in barcodes:
        score = fuzz.partial_ratio(barcode, trim_seq)

        if score > max_score:
            selected_barcode = barcode
            max_score = score

    if selected_barcode:
        barcode_seqs[selected_barcode].append(seq)


def _vcf_to_df(vcf_filename):
    '''Convert vcf to Pandas dataframe.'''
    data = []

    with open(vcf_filename) as vcf:
        for line in vcf:
            if line.startswith('##'):
                pass
            elif line.startswith('#'):
                columns = line[1:].split()[:-1] + ['DATA']
            else:
                data.append(line.split())

    df = pd.DataFrame(columns=columns, data=data)
    df[['POS']] = df[['POS']].apply(pd.to_numeric)
    df.set_index('POS', inplace=True)
    return df


def _replace_indels(sam_filename, templ_seq):
    '''Replace indels, replacing them with wildtype.'''
    sam_file = Samfile(sam_filename, 'r')

    for read in sam_file:
        # Perform mapping of nucl indices to remove spurious indels:
        seq = ''.join([read.seq[pair[0]]
                       if pair[0]
                       else templ_seq[pair[1]]
                       for pair in read.aligned_pairs
                       if pair[1] is not None])

        if seq:
            yield SeqRecord.SeqRecord(Seq.Seq(seq), read.qname, '', '')


# analyse_vcf(
#    '../results/194b82b1-e81d-4780-964a-21902a24eaab/GGTAGGCGAACACTGAGTCCAACT_3956.bam.vcf')
