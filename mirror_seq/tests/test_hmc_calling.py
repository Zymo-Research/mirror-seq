import unittest
from pandas.util.testing import assert_frame_equal
import pandas as pd
from mirror_seq import hmc_calling
import os

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

    def test_meth_call(self):
        import pysam

        samfile = pysam.AlignmentFile(os.path.join(self.data_folder, 'test.bam'))
        results = [
            (0, 11, '-', 'z'),
            (0, 12, '-', 'x'),
            (0, 15, '-', 'z'),
            (0, 22, '-', 'h'),
            (0, 25, '-', 'x'),
            (0, 26, '-', 'h'),
            (0, 27, '-', 'H'),
            (0, 30, '-', 'Z'),
            (0, 33, '-', 'Z'),
            (0, 36, '-', 'x'),
            (0, 37, '-', 'h'),
            (0, 41, '-', 'Z'),
            (0, 44, '-', 'h'),
            (0, 45, '-', 'h'),
            (0, 49, '-', 'x'),
            # Quality is not good.
            # (0, 55, '-', 'x'),
            (0, 56, '-', 'h'),
            (0, 59, '-', 'x'),
        ]
        self.assertEquals(
            results,
            list(hmc_calling.meth_call_for_read(samfile.fetch().next()))
        )
    def test_meth_call_overlap(self):
        import pysam

        samfile = pysam.AlignmentFile(os.path.join(self.data_folder, 'test.bam'))
        results = [
            (3, 11, '-', 'z'),
            (3, 12, '-', 'x'),
            (3, 15, '-', 'z'),
            (3, 22, '-', 'h'),
            (3, 25, '-', 'x'),
            (3, 26, '-', 'h'),
            (3, 27, '-', 'H'),
            (3, 30, '-', 'Z'),
            (3, 33, '-', 'Z'),
            (3, 36, '-', 'x'),
            (3, 37, '-', 'h')
        ]

        self.assertEquals(
            results,
            list(hmc_calling.meth_call_for_read(samfile.fetch('Amplicon4').next()))
        )
    def test_meth_call_by_region(self):
        chrom = 'Amplicon1'
        bam_filename = os.path.join(self.data_folder, 'test.bam')

        df = hmc_calling.meth_call_by_region(bam_filename, chrom=chrom)
        df = df[df['meth_code']=='Z'][['chrom', 'pos', 'strand', 'meth_count', 'total_count']]
        df = df.set_index(['chrom', 'pos', 'strand'])

        expected_df = pd.read_csv(os.path.join(self.data_folder, 'test_CpG.txt'), sep='\t')
        expected_df = expected_df[expected_df['chrom']==chrom]
        expected_df = expected_df.set_index(['chrom', 'pos', 'strand'])

        assert_frame_equal(expected_df, df, check_dtype=False)


    def setUp(self):
        import tempfile

        self.data_folder = os.path.join(os.path.dirname(__file__), 'data')
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil

        shutil.rmtree(self.temp_dir)

if __name__=='__main__':
    unittest.main()
