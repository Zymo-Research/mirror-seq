import unittest
from mirror_seq import trimming

class TestTrimming(unittest.TestCase):
    def test_trim_paired_seqs(self):
        seq1 = 'A'*50 + 'CGA'
        qual1 = 'H' * len(seq1)
        seq2 = 'CG' + 'T'*50
        qual2 = 'H' * len(seq2)
        read_len = 100
        expected_result = ['A'*50, 'H'*50, 'T'*50, 'H'*50]

        result = trimming.trim_paired_seqs(seq1, qual1, seq2, qual2, read_len)
        self.assertSequenceEqual(result, expected_result)

    def tests_trim_paired_seqs_single_end(self):
        seq1 = 'A'*50 + 'CGA'
        qual1 = 'H' * len(seq1)
        seq2 = None
        qual2 = None
        read_len = 100
        expected_result = ['A'*50, 'H'*50, None, None]

        result = trimming.trim_paired_seqs(seq1, qual1, seq2, qual2, read_len)
        self.assertSequenceEqual(result, expected_result)

    def test_trim_paired_seqs_with_no_read1(self):
        seq1 = None
        qual1 = None
        seq2 = None
        qual2 = None
        read_len = 100

        with self.assertRaises(Exception):
            trimming.trim_paired_seqs(seq1, qual1, seq2, qual2, read_len)

    def test_trim_paired_seqs_with_wrong_read2(self):
        seq1 = 'A'*50 + 'CGA'
        qual1 = 'H' * len(seq1)
        seq2 = None
        qual2 = 'H'*50
        read_len = 100

        with self.assertRaises(Exception):
            trimming.trim_paired_seqs(seq1, qual1, seq2, qual2, read_len)

    def test_filled_in_paired_end_trimming(self):
        import tempfile

        with tempfile.NamedTemporaryFile() as fw1, tempfile.NamedTemporaryFile() as fw2:
            trimming.filled_in_paired_end_trimming(self.r1_file.name,
                self.r2_file.name, fw1.name, fw2.name, self.read_len)
            fw1.seek(0)
            fw2.seek(0)
            self.assertEqual(self.trimmed_read1, fw1.read())
            self.assertEqual(self.trimmed_read2, fw2.read())

    def test_filled_in_single_end_trimming(self):
        import tempfile

        with tempfile.NamedTemporaryFile() as fw1:
            trimming.filled_in_paired_end_trimming(self.r1_file.name, None,
                fw1.name, None, self.read_len)
            fw1.seek(0)
            self.assertEqual(self.trimmed_read1, fw1.read())

    def test_filled_in_paired_end_trimming_gzipped_file(self):
        import tempfile
        import subprocess
        import gzip

        gzipped_r1_filename = self.r1_file.name+'.gz'
        gzipped_r2_filename = self.r2_file.name+'.gz'
        out_filename1 = 'out1.fastq.gz'
        out_filename2 = 'out2.fastq.gz'
        subprocess.check_call(('gzip', self.r1_file.name))
        subprocess.check_call(('gzip', self.r2_file.name))

        with gzip.open(out_filename1, 'w') as fw1, gzip.open(out_filename2, 'w') as fw2:
            trimming.filled_in_paired_end_trimming(gzipped_r1_filename,
                gzipped_r2_filename, fw1.name, fw2.name, self.read_len)

        with gzip.open(out_filename1) as fw1, gzip.open(out_filename2) as fw2:
            self.assertEqual(self.trimmed_read1, fw1.read())
            self.assertEqual(self.trimmed_read2, fw2.read())

        subprocess.check_call(('gunzip', gzipped_r1_filename))
        subprocess.check_call(('gunzip', gzipped_r2_filename))
        subprocess.check_call(('rm', out_filename1))
        subprocess.check_call(('rm', out_filename2))

    def test_find_read_len(self):
        expected_result = 51
        result = trimming.find_read_len(self.r1_file.name)

        self.assertEqual(expected_result, result)

    # @unittest.skip('')
    def test_run_trim_galore(self):
        import tempfile
        import subprocess

        out_dir = tempfile.mkdtemp()
        try:
            trimming.run_trim_galore(self.r1_file.name, self.r2_file.name, out_dir,
                self.adapter1, self.adapter2)
            trimming.run_trim_galore(self.r1_file.name, None, out_dir,
                self.adapter1, self.adapter2)
        except:
            raise
            self.fail('Trim Galore! failed.')
        finally:
            subprocess.call(('rm', '-r', out_dir))

    def test_main(self):
        import tempfile
        import subprocess

        out_dir = tempfile.mkdtemp()
        no_adapter_trimming = True
        try:
            trimming.main(self.r1_file.name, self.r2_file.name, out_dir, no_adapter_trimming,
                self.read_len, self.adapter1, self.adapter2)
        except:
            self.fail('Main failed.')
        finally:
            subprocess.call(('rm', '-r', out_dir))


    def setUp(self):
        import tempfile

        self.read_len = 51
        self.adapter1 = 'GATCGGAAGA'
        self.adapter2 = 'GAAGAGCACA'
        self.read1 = '''@HWI-C00124:147:C7MYBANXX:6:2102:13377:23947 1:N:0:TGACCC
CGAACCCGCCGCGTCCCCGTCTCGATCGACACCTCCGAGATCGGAAGACGA
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/
'''
        self.read2 = '''@HWI-C00124:147:C7MYBANXX:6:1304:12834:47844 1:N:0:TGACCC
CGACTCCTCATAAAACTTCCCGCCTTCTATCTCCGAGATCGGAAGAGCACA
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFF/
'''
        self.trimmed_read1 = '''@HWI-C00124:147:C7MYBANXX:6:2102:13377:23947 1:N:0:TGACCC
CGAACCCGCCGCGTCCCCGTCTCGATCGACACCTCCGAGATCGGAAGA
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
'''
        self.trimmed_read2 = '''@HWI-C00124:147:C7MYBANXX:6:1304:12834:47844 1:N:0:TGACCC
ACTCCTCATAAAACTTCCCGCCTTCTATCTCCGAGATCGGAAGAGCACA
+
BBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFF/
'''

        self.r1_file = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False)
        self.r1_file.write(self.read1)
        self.r1_file.seek(0)

        self.r2_file = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False)
        self.r2_file.write(self.read2)
        self.r2_file.seek(0)

    def tearDown(self):
        import subprocess

        subprocess.check_call(('rm', self.r1_file.name))
        subprocess.check_call(('rm', self.r2_file.name))

if __name__=='__main__':
    unittest.main()
