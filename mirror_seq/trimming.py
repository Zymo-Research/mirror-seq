def trim_paired_seqs(seq1, qual1, seq2, qual2, read_len):
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

def filled_in_paired_end_trimming(read1_filename, read2_filename, out_read1_filename,
    out_read2_filename, read_len):
    import pysam
    import os
    import gzip
    import subprocess

    fastq_file1 = pysam.FastxFile(read1_filename)
    if out_read1_filename.endswith('.gz'):
        fw1 = open(out_read1_filename[:-3], 'w')
    else:
        fw1 = open(out_read1_filename, 'w')

    fastq_file2 = None
    fw2 = None
    if read2_filename:
        fastq_file2 = pysam.FastxFile(read2_filename)
        if out_read2_filename.endswith('.gz'):
            fw2 = open(out_read2_filename[:-3], 'w')
        else:
            fw2 = open(out_read2_filename, 'w')

    for i, read1 in enumerate(fastq_file1):
        if i and i%1000000==0:
            print('{} reads processed'.format(i))

        if fastq_file2:
            read2 = fastq_file2.next()
        else:
            read2 = None

        seq1, qual1, seq2, qual2 = trim_paired_seqs(read1.sequence, read1.quality,
            read2.sequence, read2.quality, read_len)

        read_name_str1 = ' '.join([read1.name, read1.comment])
        fw1.write('@{}\n{}\n+\n{}\n'.format(read_name_str1, seq1, qual1))
        if seq2:
            read_name_str2 = ' '.join([read2.name, read2.comment])
            fw2.write('@{}\n{}\n+\n{}\n'.format(read_name_str2, seq2, qual2))

    fw1.close()
    fastq_file1.close()
    if fw2:
        fw2.close()
        fastq_file2.close()

    if out_read1_filename.endswith('.gz'):
        print('Gziping the files')
        subprocess.check_call(('gzip', '-f', fw1.name))
        if out_read2_filename:
            subprocess.check_call(('gzip', '-f', fw2.name))

def run_trim_galore(read1_filename, read2_filename, out_dir, adapter1, adapter2):
    import subprocess

    cmd = [
        'trim_galore',
        '-o', out_dir,
        '-a', adapter1
    ]

    if read2_filename:
        if adapter2:
            cmd += ['-a2', adapter2]
        cmd += [
            '--paired',
            read1_filename,
            read2_filename,
        ]
    else:
        cmd += [read1_filename]
    subprocess.check_output(cmd)

def main(read1_filename, read2_filename, out_dir, no_adapter_trimming, read_len,
    adapter1, adapter2):
    import subprocess
    import os

    is_gzipped = read1_filename.endswith('.gz')
    out_filename_template = os.path.join(out_dir, '{}_trimmed.fastq')
    if is_gzipped:
        prefix1 = os.path.splitext(os.path.splitext(read1_filename)[0])[0]
        prefix2 = os.path.splitext(os.path.splitext(read2_filename)[0])[0]
        out_filename_template += '.gz'
    else:
        prefix1 = os.path.splitext(read1_filename)[0]
        prefix2 = os.path.splitext(read2_filename)[0]
    prefix1 = os.path.basename(prefix1)
    prefix2 = os.path.basename(prefix2)
    out_read1_filename = out_filename_template.format(prefix1)
    out_read2_filename = out_filename_template.format(prefix2)
    # Trim_galore
    if not no_adapter_trimming:
        run_trim_galore(read1_filename, read2_filename, out_dir, adapter1, adapter2)
        if is_gzipped:
            read1_filename = os.path.join(out_dir,
                '{}_val_1.fq.gz'.format(os.path.basename(
                    os.path.splitext(os.path.splitext(read1_filename)[0])[0]
                )))
            read2_filename = os.path.join(out_dir,
                '{}_val_2.fq.gz'.format(os.path.basename(
                    os.path.splitext(os.path.splitext(read2_filename)[0])[0]
                )))
        else:
            read1_filename = os.path.join(out_dir,
                '{}_val_1.fq'.format(os.path.basename(os.path.splitext(read1_filename)[0])))
            read2_filename = os.path.join(out_dir,
                '{}_val_2.fq'.format(os.path.basename(os.path.splitext(read2_filename)[0])))
    # Fill-in trimming.
    filled_in_paired_end_trimming(read1_filename, read2_filename, out_read1_filename,
        out_read2_filename, read_len)
    print('done!')

def find_read_len(filename):
    import os
    import pysam

    fastq_file = pysam.FastxFile(filename)
    read = fastq_file.next()
    fastq_file.close()

    return len(read.sequence)
