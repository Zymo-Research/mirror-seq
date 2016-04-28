import random
import os

ori_read1_filename = 'in322_1_TGACCA_L006_R1_001.fastq'
ori_read2_filename = 'in322_1_TGACCA_L006_R2_001.fastq'
test_read1_filename = 'test_R1.fastq'
test_read2_filename = 'test_R2.fastq'
read_ids_filename = 'read_ids.txt'

with open(read_ids_filename) as fh:
    read_id_set = set(['@'+x.replace('_', ' ') for x in fh])

with open(ori_read1_filename) as fh, open(test_read1_filename, 'w') as fw:
    to_write = False
    for i, line in enumerate(fh):
        if i%4==0:
            to_write = line in read_id_set
        if to_write:
            fw.write(line)

with open(ori_read2_filename) as fh, open(test_read2_filename, 'w') as fw:
    to_write = False
    for i, line in enumerate(fh):
        if i%4==0:
            to_write = line.replace(' 2', ' 1') in read_id_set
        if to_write:
            fw.write(line)

offsets = set()
filesize = os.path.getsize(ori_read1_filename)
with open(ori_read1_filename) as fh:
    while len(offsets)<500:
        ## set a random pos for the file
        fh.seek(random.randint(0, filesize))
        ## read the rest line
        fh.readline()
        ## get the third line of a read
        i = 0
        while fh.readline() not in ("+\n", ""):
            i +=1
            if i>4:
                raise Exception("Probably wrong fastq format.")
        ## read the fourth line
        if fh.readline()!="":
            offset = fh.tell()
            if fh.readline() not in read_id_set:
                offsets.add(offset)

utils.write_reads_from_file_offsets(ori_read1_filename, 'tmp_1.fastq', offsets)

offsets = set()
filesize = os.path.getsize(ori_read1_filename)
with open(ori_read2_filename) as fh:
    while len(offsets)<500:
        ## set a random pos for the file
        fh.seek(random.randint(0, filesize))
        ## read the rest line
        fh.readline()
        ## get the third line of a read
        i = 0
        while fh.readline().replace(' 2', ' 1') not in ("+\n", ""):
            i +=1
            if i>4:
                raise Exception("Probably wrong fastq format.")
        ## read the fourth line
        if fh.readline()!="":
            offset = fh.tell()
            if fh.readline() not in read_id_set:
                offsets.add(offset)

utils.write_reads_from_file_offsets(ori_read1_filename, 'tmp_2.fastq', offsets)
