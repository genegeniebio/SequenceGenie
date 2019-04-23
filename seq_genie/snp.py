'''
SequenceGenie (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
from collections import defaultdict
import os.path
import sys

import pandas as pd
from seq_genie import vcf_utils


def analyse(dir_name):
    '''Analyse.'''
    dfs = _get_dfs(dir_name)

    for key, values in dfs.items():
        for direction, df in values.items():
            if 'INDEL' in df.columns:
                df.drop(df[df['INDEL']].index, inplace=True)

            nucl_counts = df.apply(_get_nucl_count, axis=1)

            df = pd.concat(
                [df,
                 pd.DataFrame(nucl_counts.tolist(),
                              index=df.index,
                              columns=['A', 'C', 'G', 'T'])],
                axis=1)

            df['correct'] = df[['A', 'C', 'G', 'T']].sum(axis=1) == df['DP']

            # assert((df[['A', 'C', 'G', 'T']].sum(axis=1) == df['DP']).all())

            df.set_index('POS', inplace=True)
            df.to_csv('_'.join(list(key) + [direction]) + '.csv')


def _get_dfs(dir_name):
    '''Get DataFrames.'''
    dfs = defaultdict(defaultdict)

    count = 0

    for dirpath, _, filenames in os.walk(os.path.abspath(dir_name)):
        for filename in filenames:
            if filename.endswith('.vcf'):
                df, _ = vcf_utils.vcf_to_df(os.path.join(dirpath, filename))
                name = os.path.basename(
                    os.path.abspath(os.path.join(dirpath, os.pardir)))
                tokens = name.split('_')
                dfs[tuple(tokens[:2])][tokens[2]] = df

                count += 1

                if count == 3:
                    return dfs

    return dfs


def _get_nucl_count(row):
    '''Get nucleotide count.'''
    counts = {nucl: 0 for nucl in 'ACGT'}
    nucls = (row['REF'] + ',' + row['ALT']).split(',')
    proportions = [round(float(prop) * row['DP'])
                   for prop in row['QS'].split(',')]
    counts.update(dict(zip(nucls, proportions)))
    counts.pop('<*>', None)
    return list(counts.values())


def main(args):
    '''main method.'''
    analyse(args[0])


if __name__ == '__main__':
    main(sys.argv[1:])
