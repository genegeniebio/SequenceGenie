'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
import os
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy as np
import pandas as pd


def analyse(vcf_filename, target_id, src_id, dp_filter, write_queue):
    '''Analyse a given vcf file.'''
    num_matches, mutations, indels, deletions, templ_len, \
        consensus_seq, depths = _analyse_vcf(vcf_filename, dp_filter)

    consensus_filename = os.path.join(os.path.dirname(vcf_filename),
                                      'consensus.fasta')
    record = SeqRecord(Seq(consensus_seq), id=vcf_filename)
    SeqIO.write([record], consensus_filename, 'fasta')

    write_queue.put(
        ['identity', num_matches / float(templ_len), target_id, src_id])
    write_queue.put(
        ['mutations', mutations, target_id, src_id])
    write_queue.put(
        ['indels', indels, target_id, src_id])
    write_queue.put(
        ['deletions', deletions, target_id, src_id])
    write_queue.put(
        ['depths', max(depths) if depths else 0, target_id, src_id])


def _analyse_vcf(vcf_filename, dp_filter):
    '''Analyse vcf file, returning number of matches, mutations and
    indels.'''
    num_matches = 0
    mutations = []
    indels = []
    deletions = []
    consensus_seq = []
    depths = []

    df, templ_len = _vcf_to_df(vcf_filename)

    for _, row in df.iterrows():
        if 'INDEL' in row and row.INDEL:
            indels.append(row['REF'] + str(row['POS']) + row['ALT'])
        elif (dp_filter > 1 and row['DP'] > dp_filter) \
                or row['DP_PROP'] > dp_filter:
            alleles = [row['REF']] + row['ALT'].split(',')

            # Extract QS values and order to find most-likely base:
            qs = [float(val)
                  for val in dict([term.split('=')
                                   for term in row['INFO'].split(';')])
                  ['QS'].split(',')]

            # Compare most-likely base to reference:
            hi_prob_base = alleles[np.argmax(qs)]
            consensus_seq.append(hi_prob_base)

            if row['REF'] != hi_prob_base:
                mutations.append(row['REF'] + str(row['POS']) +
                                 hi_prob_base + ' ' + str(max(qs)))
            else:
                num_matches += 1

            depths.append(row['DP'])
        else:
            deletions.append(row['POS'])

    return num_matches, mutations, indels, _get_ranges_str(deletions), \
        templ_len, ''.join(consensus_seq), depths


def _vcf_to_df(vcf_filename):
    '''Convert vcf to Pandas dataframe.'''
    data = []
    templ_len = float('NaN')

    with open(vcf_filename) as vcf:
        for line in vcf:
            if line.startswith('##'):
                mtch = re.match(r'(?:.*),length=(.*)>', line)

                if mtch:
                    templ_len = int(mtch.group(1))
                else:
                    pass
            elif line.startswith('#'):
                columns = line[1:].split()[:-1] + ['DATA']
            else:
                data.append(line.split())

    df = _expand_info(pd.DataFrame(columns=columns, data=data))

    df['POS'] = df['POS'].astype(int)

    if 'DP' in df.columns:
        df['DP'] = df['DP'].astype(int)
        df['DP_PROP'] = df['DP'] / df['DP'].max()

    if 'INDEL' in df.columns:
        df[['INDEL']] = df[['INDEL']].fillna(value=False)

    return df, templ_len


def _expand_info(df):
    '''Expand out INFO column from vcf file.'''
    infos = []

    for row in df.itertuples():
        info = [term.split('=') for term in row.INFO.split(';')]

        infos.append({term[0]: (term[1] if len(term) == 2 else True)
                      for term in info})

    return df.join(pd.DataFrame(infos, index=df.index))


def _get_ranges_str(vals):
    '''Convert list of integers to range strings.'''
    return ['-'.join([str(r) for r in rnge])
            if rnge[0] != rnge[1]
            else rnge[0]
            for rnge in _get_ranges(vals)]


def _get_ranges(vals):
    '''Convert list of integer to ranges.'''
    ranges = []

    if vals:
        sorted_vals = sorted(vals)

        i = 0

        for j in range(1, len(sorted_vals)):
            if sorted_vals[j] > 1 + sorted_vals[j - 1]:
                ranges.append((sorted_vals[i], sorted_vals[j - 1]))
                i = j

        ranges.append((sorted_vals[i], sorted_vals[-1]))

    return ranges
