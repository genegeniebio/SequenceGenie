'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=wrong-import-order
import os
import unittest

import pandas as pd
from seq_genie import vcf_utils


class Test(unittest.TestCase):
    '''Class to test utils module.'''

    def test_write_summary(self):
        '''Test demuliplex method.'''
        directory = os.path.dirname(os.path.realpath(__file__))

        summary_df = pd.read_csv(os.path.join(directory, 'summary.csv'))
        identity_df = pd.read_csv(os.path.join(directory, 'identity.csv'))
        mutations_df = pd.read_csv(os.path.join(directory, 'mutations.csv'))
        indels_df = pd.read_csv(os.path.join(directory, 'indels.csv'))
        deletions_df = pd.read_csv(os.path.join(directory, 'deletions.csv'))
        depths_df = pd.read_csv(os.path.join(directory, 'depths.csv'))

        vcf_utils.write_summary(summary_df,
                                identity_df,
                                mutations_df,
                                indels_df,
                                deletions_df,
                                depths_df)

        summary_df.to_csv(os.path.join(directory, 'summary_out.csv'))


if __name__ == "__main__":
    unittest.main()
