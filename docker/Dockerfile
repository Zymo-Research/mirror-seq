FROM python:2
MAINTAINER Hunter Chung <hchung@zymoresearch.com>

RUN apt-get update && \
apt-get install -y \
    samtools \
    curl \
    zip \
    unzip \
    --no-install-recommends && \
rm -rf /var/lib/apt/lists/* && \
apt-get clean autoclean && \
apt-get autoremove -y

ADD install_bioinfo_packages.sh install_python_packages.sh /tmp/
RUN bash /tmp/install_bioinfo_packages.sh
RUN bash /tmp/install_python_packages.sh

ENV PATH=${PATH}:/genomicTools:/bowtie2-2.2.6:/bismark_v0.14.5 \
  PYTHONPATH=${PYTHONPATH}:/scripts/mirror-seq \
  TERM=xterm

WORKDIR /mnt
