'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-few-public-methods
import itertools
import os
import re
import sys

import numpy as np
import pandas as pd


class VcfAnalyser(object):
    ''' Class to analyse vcf files.'''

    def __init__(self, columns, src_filename, dir_name='.'):

        # Initialise summary:
        self.__summary_df = pd.read_csv(src_filename)
        self.__summary_df.dropna(inplace=True)

        src_id_cols = ['forward', 'reverse']
        self.__src_ids = [tuple(pair)
                          for pair in
                          self.__summary_df[src_id_cols].values.tolist()]

        self.__summary_df.set_index(src_id_cols, inplace=True)

        # Initialise specific Dataframes:
        index = self.__summary_df.index
        self.__mutations_df = pd.DataFrame(columns=columns, index=index)
        self.__indels_df = pd.DataFrame(columns=columns, index=index)
        self.__deletions_df = pd.DataFrame(columns=columns, index=index)
        self.__identity_df = pd.DataFrame(columns=columns, index=index,
                                          dtype='float')

        self.__dir_name = os.path.abspath(dir_name)

        if not os.path.exists(self.__dir_name):
            os.makedirs(self.__dir_name)

    def get_src_ids(self):
        '''Get source ids.'''
        return self.__src_ids

    def analyse(self, vcf_filename, target_id, src_id, dp_filter):
        '''Analyse a given vcf file.'''
        num_matches, mutations, indels, deletions, templ_len = \
            analyse_vcf(vcf_filename, dp_filter)

        _set_value(self.__identity_df, num_matches / float(templ_len),
                   target_id, *src_id)

        _set_value(self.__mutations_df, mutations, target_id, *src_id)
        _set_value(self.__indels_df, indels, target_id, *src_id)
        _set_value(self.__deletions_df, deletions, target_id, *src_id)

        self.write()

    def write(self):
        '''Write output files.'''
        self.__identity_df.to_csv(os.path.join(self.__dir_name,
                                               'identity.csv'))
        self.__mutations_df.to_csv(os.path.join(self.__dir_name,
                                                'mutations.csv'))
        self.__indels_df.to_csv(os.path.join(self.__dir_name, 'indels.csv'))
        self.__deletions_df.to_csv(os.path.join(self.__dir_name,
                                                'deletions.csv'))

    def write_summary(self):
        '''Write summary output file.'''
        self.__identity_df.fillna(0, inplace=True)
        self.__summary_df['ice_id'] = self.__identity_df.idxmax(axis=1)
        self.__summary_df['identity'] = self.__identity_df.max(axis=1)
        self.__summary_df['mutations'] = \
            self.__mutations_df.lookup(self.__mutations_df.index,
                                       self.__summary_df['ice_id'])
        self.__summary_df['indels'] = \
            self.__indels_df.lookup(self.__indels_df.index,
                                    self.__summary_df['ice_id'])
        self.__summary_df['deletions'] = \
            self.__deletions_df.lookup(self.__deletions_df.index,
                                       self.__summary_df['ice_id'])

        # Remove spurious unidentified entries:
        self.__summary_df = \
            self.__summary_df[self.__summary_df['identity'] != 0]

        self.__summary_df.to_csv(os.path.join(self.__dir_name, 'summary.csv'))


def analyse_vcf(vcf_filename, dp_filter):
    '''Analyse vcf file, returning number of matches, mutations and indels.'''
    num_matches = 0
    mutations = []
    indels = []
    deletions = []

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

            if row['REF'] != hi_prob_base:
                mutations.append(row['REF'] + str(row['POS']) + hi_prob_base +
                                 ' ' + str(max(qs)))
            else:
                num_matches += 1
        else:
            deletions.append(row['POS'])

    return num_matches, mutations, indels, _get_ranges_str(deletions), \
        templ_len


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

    for _, b in itertools.groupby(enumerate(vals), lambda (x, y): y - x):
        b = list(b)
        ranges.append((b[0][1], b[-1][1]))

    return ranges


def _set_value(df, val, col_id, *row_ids):
    '''Set a value in a Dataframe.'''
    df[col_id].loc[row_ids[0], row_ids[1]] = val


def main(args):
    '''main method.'''
    analyser = VcfAnalyser(args[4:], args[0], args[2])

    dp_filter = float(args[3])

    for dirpath, _, filenames in os.walk(os.path.abspath(args[1])):
        for filename in filenames:
            if filename[-4:] == '.vcf':
                mtch = re.match(r'([^\.]*)_([^\.]*)\..*', filename)
                analyser.analyse(os.path.join(dirpath, filename),
                                 mtch.group(2),
                                 mtch.group(1).split('_'),
                                 dp_filter)

    analyser.write_summary()


if __name__ == '__main__':
    main(sys.argv[1:])