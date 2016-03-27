import unittest
from pandas.util.testing import assert_frame_equal
import pandas as pd
import hmc_calling

class TestHmc_calling(unittest.TestCase):
    def test_mirror_seq_conversion(self):
        columns = ['chrom', 'pos', 'strand', 'meth_count', 'total_count']
        df = pd.DataFrame([
            ['chr1', 10, '+', 0, 10],
            ['chr1', 20, '-', 5, 12],
        ], columns=columns)
        expected_df = pd.DataFrame([
            ['chr1', 11, '-', 10, 10],
            ['chr1', 19, '+', 7, 12],
        ], columns=columns)

        hmc_calling.mirror_seq_conversion(df)
        assert_frame_equal(df, expected_df)

    def test_mirror_seq_conversion_miss_col(self):
        columns = ['chrom', 'pos', 'meth_count', 'total_count']
        df = pd.DataFrame([
            ['chr1', 10, 0, 10],
        ], columns=columns)

        with self.assertRaises(KeyError):
            hmc_calling.mirror_seq_conversion(df)

if __name__=='__main__':
    unittest.main()
