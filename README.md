## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Citation](#citation)


## Introduction
A tool for automated assembly and circularization of mitochondrial genomes in trypanosomatids. The input is reads in BAM or FASTQ format, and the output is circularized minicircles in FASTA format.

komics was developed at the Antwerp Institute of Tropical Medicine and the University of Edinburgh. 
Please report any issues to fvandenbroeck AT itg.be


## Installation
KOMICS has the following dependencies, which need to be installed first:
  * [MEGAHIT](http://www.metagenomics.wiki/tools/assembly/megahit)
  * [TRIMMOMATIC](http://www.usadellab.org/cms/?page=trimmomatic)
  * [BLASTn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  * [SMALT](https://www.sanger.ac.uk/science/tools/smalt-0)
  * [VSEARCH](https://github.com/torognes/vsearch)

Note that you can use the environment variable $TRIMMOMATIC_DIR to specify the directory containing the Trimmomatic jar file. If this environment variable is set, then it is used by KOMICS. Otherwise, you will need to specify the directory with --dir.

Once the dependencies are installed, install KOMICS using pip3:
```
pip3 install komics
```


## Usage
```
Usage: komics <command> [options] <required arguments>

To get minimal usage for a command use:
komics command

To get full help for a command use one of:
komics command -h
komics command --help


Available commands:
all:         	 Performs bam2fq, trimfq, assemble, circularize, polish, qualcheck
bam2fq:      	 Writes paired-end reads from a BAM file to FASTQ files
trimfq:      	 Cleans and filters reads for high quality
assemble:    	 Assembles minicircles using high quality reads from trimfq
circularize: 	 Circularizes minicircles from assemble
polish:      	 Reorientate and filter circular minicircles
qualcheck:       Estimates read counts and depths per minicircle and overall
```


## Citation
If you use this software please cite:

__Ecological divergence and hybridization of Neotropical Leishmania parasites__   
Van den Broeck et al. BIORXIV 2019 doi: [10.1101/824912](https://www.biorxiv.org/content/10.1101/824912v1)

