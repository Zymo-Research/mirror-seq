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

    @unittest.skip('It takes too long for test.')
    def test_adapter_trimming(self):
        out_dir = './'
        trimming.run_trim_galore(self.r1_file.name, self.r2_file.name, out_dir)

    def test_filled_in_paired_end_trimming(self):
        import tempfile

        expected_read1 = '''@HWI-C00124:147:C7MYBANXX:6:2102:13377:23947 1:N:0:TGACCC
CGAACCCGCCGCGTCCCCGTCTCGATCGACACCTCCGAGATCGGAAGA
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
'''
        expected_read2 = '''@HWI-C00124:147:C7MYBANXX:6:1304:12834:47844 1:N:0:TGACCC
ACTCCTCATAAAACTTCCCGCCTTCTATCTCCGAGATCGGAAGAGCACA
+
BBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFF/
'''
        with tempfile.NamedTemporaryFile() as fw1, tempfile.NamedTemporaryFile() as fw2:
            trimming.filled_in_paired_end_trimming(self.r1_file.name,
                self.r2_file.name, fw1.name, fw2.name, self.read_len)
            fw1.seek(0)
            fw2.seek(0)
            self.assertEqual(expected_read1, fw1.read())
            self.assertEqual(expected_read2, fw2.read())

    def test_find_read_len(self):
        expected_result = 51
        result = trimming.find_read_len(self.r1_file.name)

        self.assertEqual(expected_result, result)

    def setUp(self):
        import tempfile

        self.read_len = 51
        self.r1_file = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False)
        self.r1_file.write('''@HWI-C00124:147:C7MYBANXX:6:2102:13377:23947 1:N:0:TGACCC
CGAACCCGCCGCGTCCCCGTCTCGATCGACACCTCCGAGATCGGAAGACGA
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/
''')
        self.r1_file.seek(0)

        self.r2_file = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False)
        self.r2_file.write('''@HWI-C00124:147:C7MYBANXX:6:1304:12834:47844 1:N:0:TGACCC
CGACTCCTCATAAAACTTCCCGCCTTCTATCTCCGAGATCGGAAGAGCACA
+
BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFF/
''')
        self.r2_file.seek(0)

    def tearDown(self):
        import subprocess

        subprocess.check_call(('rm', self.r1_file.name))
        subprocess.check_call(('rm', self.r2_file.name))

if __name__=='__main__':
    unittest.main()
