# Add Trim Galore! 0.3.7
curl -L "http://www.bioinformatics.bbsrc.ac.uk/projects/trim_galore/trim_galore_v0.3.7.zip" > trim_galore_v0.3.7.zip
unzip trim_galore_v0.3.7.zip
rm trim_galore_v0.3.7.zip
# Add bowtie2 2.2.6
curl -L \
"http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip" > bowtie2-2.2.6-linux-x86_64.zip
unzip bowtie2-2.2.6-linux-x86_64.zip
rm bowtie2-2.2.6-linux-x86_64.zip
# Add Bismark 0.14.5
curl -L "http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/bismark_v0.14.5.tar.gz" > /bismark_v0.14.5.tar.gz
tar zxvf /bismark_v0.14.5.tar.gz -C /
rm /bismark_v0.14.5.tar.gz
