# What is it
Mirror-seq is a [hydroxymethylation](http://www.zymoresearch.com/epigenetics/dna-hydroxymethylation) (hmc) assay invented by [Zymo Research](http://zymoresearch.com) in genomes using [bisulfite sequencing](http://www.zymoresearch.com/bisulfite-beginner-guide). This analysis
tool helps biologists to analyze sequencing data. It takes Fastq files from sequencers and generate hydroxymethylation ratio for CpGs.

# Where Should I Start
We provide three level for scientists to use our tool:
* **Newbies** We have a [small test dataset](https://github.com/Zymo-Research/mirror-seq/wiki/Test-Dataset) for you to play with.
* **Experieced** Follow the [Quick Start](https://github.com/Zymo-Research/mirror-seq/wiki/Quick-Start) to try it with your own data.
* **Expert** You have your homebrew bioinformatics software. Just follow the instruction below to install and run the specific parts for Mirror-seq.

# installation
You need to install the following bioinformatics software below and put them in PATH. `Pip` can install the python package Dependencies for you.
```pip install mirror_seq```

# Dependencies
## Python (2.7)
* [NumPy](http://www.numpy.org/): 1.7.0
* [pandas](http://pandas.pydata.org/): 0.18.0
* [numexpr](https://github.com/pydata/numexpr): 2.5.2
* [pysam](http://pysam.readthedocs.org/en/latest/api.html): 0.9.0
* [cutadapt](http://cutadapt.readthedocs.org/en/stable/guide.html): 1.9.1
* [PyTables](http://www.pytables.org/): 3.2.2

## Bioinformatics software
* [Trim Galore!](http://www.bioinformatics.bbsrc.ac.uk/projects/trim_galore/): 0.3.7
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): 2.2.6
* [Bismark](http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/): 0.14.5

# Usage
We provide three commands for more details of each command, please use `--help`:
* mirror-seq: Run the entire workflow.
* mirror-trimming: Only run the first part of filled-in nucleotides trimming.
* mirror-calling: Only run the last part of hydroxymethylation calling.
