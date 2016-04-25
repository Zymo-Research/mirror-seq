#!/bin/bash
cd /scripts
git clone git@github.com:Zymo-Research/mirror-seq.git

sudo mkfs.btrfs /dev/xvdb
sudo mount /dev/xvdb /mnt
sudo chmod 777 -R /mnt/
cd /mnt
docker run --rm -v /mnt/:/mnt anigeo/awscli s3 cp s3://zymo-files/mirror-seq/hg19.tar.gz /mnt/
tar xvf /mnt/hg19.tar.gz
docker run -it --rm -v /mnt/:/mnt \
  -v /scripts/mirror-seq/:/scripts/mirror-seq \
  zymoresearch/mirror-seq bash

cd /scripts/mirror-seq
python bin/mirror-seq -1 /mnt/in367_1_s1_R1.fastq.gz -2 /mnt/in367_1_s1_R2.fastq.gz -g /mnt/hg19 -X 1000 --non_directional --bed
