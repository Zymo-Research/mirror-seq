import unittest
import trimming

class TestTrimming(unittest.TestCase):
    def test_trim(self):
        seq1 = 'A'*50 + 'CGA'
        qual1 = 'H' * len(seq1)
        seq2 = 'CG' + 'T'*50
        qual2 = 'H' * len(seq2)
        read_len = 100
        expected_result = ['A'*50, 'H'*50, 'T'*50, 'H'*50]

        result = trimming.trim(seq1, qual1, seq2, qual2, read_len)
        self.assertSequenceEqual(result, expected_result)

    def tests_trim_single_end(self):
        seq1 = 'A'*50 + 'CGA'
        qual1 = 'H' * len(seq1)
        seq2 = None
        qual2 = None
        read_len = 100
        expected_result = ['A'*50, 'H'*50, None, None]

        result = trimming.trim(seq1, qual1, seq2, qual2, read_len)
        self.assertSequenceEqual(result, expected_result)

    def test_trim_with_no_read1(self):
        seq1 = None
        qual1 = None
        seq2 = None
        qual2 = None
        read_len = 100

        with self.assertRaises(Exception):
            trimming.trim(seq1, qual1, seq2, qual2, read_len)

    def test_trim_with_wrong_read2(self):
        seq1 = 'A'*50 + 'CGA'
        qual1 = 'H' * len(seq1)
        seq2 = None
        qual2 = 'H'*50
        read_len = 100

        with self.assertRaises(Exception):
            trimming.trim(seq1, qual1, seq2, qual2, read_len)


if __name__=='__main__':
    unittest.main()
