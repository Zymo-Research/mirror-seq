<p align="center">
    <img src="https://travis-ci.org/Zymo-Research/mirror-seq.svg?branch=master" height="130">
</p>

# What is it
Mirror-seq is a [hydroxymethylation](http://www.zymoresearch.com/epigenetics/dna-hydroxymethylation) (hmc) assay invented by [Zymo Research](http://zymoresearch.com) in genomes using [bisulfite sequencing](http://www.zymoresearch.com/bisulfite-beginner-guide). This analysis
tool helps biologists to analyze sequencing data. It takes Fastq files from sequencers and generate hydroxymethylation ratio for CpGs.

# Where Should I Start
There are three levels for scientists to use our tool:
* **Newbies** We provide [tuturial](https://github.com/Zymo-Research/mirror-seq/wiki/tutorial) for you to get familiar with it.
* **Experieced** Follow the [Quick Start](https://github.com/Zymo-Research/mirror-seq/wiki/Quick-Start) to try it with your own data.
* **Expert** You have your homebrew bioinformatics software. No problem! Just follow the instruction below to install and run the specific parts for Mirror-seq.

# Installation
You need to install the following bioinformatics software in the dependencies and put them in PATH to run the full workflow. However, if you have your own trimming and alignment software, you can skip it.

```pip install mirror_seq```
* Note: `pip` can install all the dependencies for you.

## Dependencies
### Python (2.7)
* [NumPy](http://www.numpy.org/): 1.7.0
* [pandas](http://pandas.pydata.org/): 0.18.0
* [pysam](http://pysam.readthedocs.org/en/latest/api.html): 0.9.0
* [cutadapt](http://cutadapt.readthedocs.org/en/stable/guide.html): 1.9.1

### Bioinformatics software
* [Trim Galore!](http://www.bioinformatics.bbsrc.ac.uk/projects/trim_galore/): 0.3.7
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): 2.2.6
* [Bismark](http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/): 0.14.5

# Usage
We provide three commands for more details of each command, please use `--help`:

## Trimming
`mirror-trim` trims off Mirror-seq specific filled-in nucleotides and also do adapter trimming and quality trimming.

### Output file
* **< PREFIX >_trimmed.fastq** The trimmed fastq file.

## Hydroxymethylation Calling
`mirror-call` calls hydroxymethylation ratios for CpGs from alignment files.

### Output files
* **< PREFIX>_CpG.csv.gz** Each row represents a CpG. The columns are:
  * **chrom** The chromosome name of this CpG.
  * **pos** The chromosomal position of this CpG.
  * **strand** Either forward strand or reverse strand.
  * **meth_count** Number of reads aligned at the CpG which are hydroxymethylated.
  * **total_count** The total number of reads aligned at the CpG.
* **< PREFIX >_CpG.bed.gz** Browser tracks can be loaded in [USCS Genome Browser](http://genome.ucsc.edu/) or [igv](https://www.broadinstitute.org/igv/) to visualize hydroxymethylation data. This is the standard [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) with 8 fields. The name and score fields need more description.
  * **name** is formatted as < HYDROXYMETHYLATED READ COUNT >/< TOTAL READ COUNT >(< HYDROXYMETHYLATION RATIO >). For example, 0/3(0%) means non of the three reads at the CpG position is hydroxymethylated. The hydroxymethylation ratio is 0%.
  * **score** hydroxymethylation percentage times 1000.

## Entire Workflow
`mirror-seq` command takes fastq files from sequencer and output the hydroxymethylation calling files.
### Output files
The combination of the two commands above.
