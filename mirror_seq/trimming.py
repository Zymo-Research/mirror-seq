def trim(seq1, qual1, seq2, qual2, read_len):
    if not seq1 or not qual1:
        raise Exception('seq1 and qual1 are both required.')
    if (seq2==None) ^ (qual2==None):
        raise Exception('Cannot only has seq2 or qual2.')

    if len(seq1)<=read_len and seq1.endswith('CGA'):
        seq1 = seq1[:-3]
        qual1 = qual1[:-3]

    if seq2 and qual2:
        seq2 = seq2[2:]
        qual2 = qual2[2:]

    return seq1, qual1, seq2, qual2
