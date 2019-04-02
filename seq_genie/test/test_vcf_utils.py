'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=wrong-import-order
import os
import unittest

from seq_genie import vcf_utils


class Test(unittest.TestCase):
    '''Class to test utils module.'''

    def test_vcf_to_df(self):
        '''Test vcf_to_df method.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        vcf_filename = os.path.join(directory, 'variants.vcf')
        df, _ = vcf_utils.vcf_to_df(vcf_filename)
        # df_filename = os.path.join(directory, 'variants_out.csv')
        # df.to_csv(df_filename)
        self.assertFalse(df.empty)

    def test_analyse_vcf(self):
        '''Test analyse_vcf method.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        vcf_filename = os.path.join(directory, 'variants.vcf')
        _, _, indels, _, _, _, _ = vcf_utils.analyse_vcf(vcf_filename, 0)

        self.assertAlmostEqual(indels[0][1], 0.563636, 3)


if __name__ == "__main__":
    unittest.main()
