FROM python:2
MAINTAINER Hunter Chung <b89603112@gmail.com>

RUN apt-get update && \
apt-get install -y \
    samtools \
    curl \
    zip \
    unzip \
    pkg-config \
    --no-install-recommends && \
rm -rf /var/lib/apt/lists/* && \
apt-get clean autoclean && \
apt-get autoremove -y

ADD install_bioinfo_packages.sh /tmp/
RUN bash /tmp/install_bioinfo_packages.sh

RUN pip install mirror_seq

ENV PATH=${PATH}:/genomicTools:/bowtie2-2.2.6:/bismark_v0.14.5:/trim_galore_zip \
  TERM=xterm

WORKDIR /workspace
ENTRYPOINT ["mirror-seq"]
