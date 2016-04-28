# What is it
Mirror-seq is a assay invented by [Zymo Research](http://zymoresearch.com) to detect
[hydroxymethylation]() (hmc) in genomes using [bisulfite sequencing](). This analysis
tool helps biologists to analyze sequencing data. It takes Fastq files from sequencers and generate hydroxymethylation ratio for each CpGs.

# What's included
We provide two way to do analysis. If you are new or like to do a quick analysis, the **Qucik Start** below is the total solution you would like. It set up the environment for you and has a simple workflow - trimming, alignment, hydroxymethylation calling. Just type in one command. You will find results in your workplace in several hours.

If you are bioinformatics expert or anyone eager to try different parameters, we also provide the component of initial fill-in nucleotides trimming () and the final hydroxymethylation calling () parts as standalone program. You can plug in your favor QC and adapter trimming software and alignment software with your homemade parameters. Please fellow the **installation** section below for more details.

# Quick Start
We created a [Docker]() image to solve the dependency problem and scientists can use either Windows, MacOS, or Linux to run the analysis.

## Install Docker
Find your OS and follow the installation instructions of [Windows](https://docs.docker.com/windows/step_one/), [MacOS](https://docs.docker.com/mac/step_one/), and [Linux](https://docs.docker.com/linux/step_one/) from Docker's official website.

## Run Mirror-seq
You need to create a workplace directory (`<YOUR WORKPLACE>`) and put the following files inside:
* Read 1 and Read 2 Fastq files.
* Genome index (We provide [human index](). Unzip the file after downloading.)

```
docker run -it --rm -v <YOUR WORKPLACE>:/workplace \
  zymoresearch/mirror-seq \
  -1 <READ 1 FILENAME> -2 <READ 2 FILENAME> \
  -g <GENOME INDEX FOLDER NAME>
```

## Notes:
Although it is super easy to run the analysis tool, there are several things you need to know in order to run it smoothly.
* The alignment part is memory intensive and CPU intensive. [Bismark](http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/), the aligner we used in our tool, suggests at least 5 cores and > 16GB of RAM.
* Usually Fastq files are several GB even with compression. In the first trimming part, the tool could need up to 3 times large as the original input. Please make sure your workplace has enough storage space.

# installation

# Dependencies
## Python (2.7)
* [NumPy](http://www.numpy.org/): 1.7.0
* [pandas](http://pandas.pydata.org/): 0.18.0
* [numexpr](https://github.com/pydata/numexpr): 2.5.2
* [pysam](http://pysam.readthedocs.org/en/latest/api.html): 0.9.0
* [cutadapt](http://cutadapt.readthedocs.org/en/stable/guide.html): 1.9.1
* [PyTables](http://www.pytables.org/): 3.2.2

## Bioinformatics software
* [bedToBigBed](http://hgdownload.cse.ucsc.edu/admin/exe/)
* [bedSort](http://hgdownload.cse.ucsc.edu/admin/exe/)
* [Trim Galore!](http://www.bioinformatics.bbsrc.ac.uk/projects/trim_galore/): 0.3.7
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): 2.2.6
* [Bismark](http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/): 0.14.5
