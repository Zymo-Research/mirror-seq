mkdir /genomicTools
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed -O /genomicTools/bedToBigBed
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedSort -O /genomicTools/bedSort
chmod 755 /genomicTools/*

curl -L \
"http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip" > bowtie2-2.2.6-linux-x86_64.zip
unzip bowtie2-2.2.6-linux-x86_64.zip
rm bowtie2-2.2.6-linux-x86_64.zip

curl -L "http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/bismark_v0.14.5.tar.gz" > /bismark_v0.14.5.tar.gz
tar zxvf /bismark_v0.14.5.tar.gz -C /
rm /bismark_v0.14.5.tar.gz
