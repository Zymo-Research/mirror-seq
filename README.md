# What is it
Mirror-seq is a [hydroxymethylation]() (hmc) assay invented by [Zymo Research](http://zymoresearch.com) in genomes using [bisulfite sequencing](). This analysis
tool helps biologists to analyze sequencing data. It takes Fastq files from sequencers and generate hydroxymethylation ratio for CpGs.

# What's included
We provide two way to do analysis.

If you are new or like to do a quick analysis, the **Qucik Start** below is the easiest 1 step to run the analysis. Just type in the command below. You will find results in your workplace in several hours.

If you are the bioinformatics expert or anyone eager to try different parameters, we also provide the mirror-seq specific component of initial fill-in nucleotides trimming () and the final hydroxymethylation calling () as standalone commands. You can plug in your favor QC and adapter trimming software and alignment software with your homemade parameters. Please fellow the **installation** section below for more details.

# Quick Start
We created a [Docker]() image to solve the dependency problem and you can use either Windows, MacOS, or Linux to run the analysis.

## 0 - Install Docker
Find your OS and follow the installation instructions of [Windows](https://docs.docker.com/windows/step_one/), [MacOS](https://docs.docker.com/mac/step_one/), or [Linux](https://docs.docker.com/linux/step_one/) from Docker's official website.

## 1 - Run Mirror-seq
You need to create a workplace directory (`<YOUR WORKPLACE>`) and put the following files inside:
* Read 1 and Read 2 Fastq files.
* Bisulfite converted genome index created by Bismark. (We provide [human index](). Unzip the file after downloading.)

Run this command in your console.
```
docker run -it --rm -v <YOUR WORKPLACE>:/workplace \
  zymoresearch/mirror-seq \
  -1 <READ 1 FILENAME> -2 <READ 2 FILENAME> \
  -g <GENOME INDEX FOLDER NAME> --bed
```

## Notes:
Although it is super easy to run the analysis tool, there are several things you need to know in order to run it smoothly.
* The alignment part is memory and CPU intensive. [Bismark](http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/), the aligner we used in our tool, suggests at least 5 cores and > 16GB of RAM.
* Usually Fastq files are several GB with compression. In the first trimming part, the tool may need up to 3X large as the original input. Please make sure your workplace has enough storage space.

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
