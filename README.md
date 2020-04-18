## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Tutorial](#tutorial)
  * [Reading](#reading)
  * [Citation](#citation)


## Introduction
A tool for automated assembly and circularization of mitochondrial genomes in trypanosomatids. The input is reads in BAM or FASTQ format, and the output is circularized minicircles in FASTA format.

komics is described in detail here:
__Ecological divergence and hybridization of Neotropical Leishmania parasites__
Van den Broeck et al. BIORXIV 2019 doi: [10.1101/824912](https://www.biorxiv.org/content/10.1101/824912v1)

Please report any issues to fvandenbroeckATitg.be


## Installation
KOMICS has the following dependencies, which need to be installed first:
  * [MEGAHIT](http://www.metagenomics.wiki/tools/assembly/megahit)
  * [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  * [VSEARCH](https://github.com/torognes/vsearch)

Once the dependencies are installed, install the latest version of KOMICS using pip3:
```
pip3 install git+https://github.com/FreBio/komics.git
```

If you are running komics on a supercomputer, you may want to have a look (here)[https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/software/python_package_management.html#alternatives-to-conda] link first on how to setup your local environment.


## Usage
```
Usage: komics <command> [options] <required arguments>

To get minimal usage for a command use:
komics command

To get full help for a command use one of:
komics command -h
komics command --help

Available commands:
all:         	 Performs assemble, circularize and polish
assemble:    	 Assembles minicircles using high quality reads from trimfq
circularize: 	 Circularizes minicircles from assemble
polish:      	 Reorientate and filter circular minicircles
```


## Tutorial

### 1. Preparing your data
If you have a BAM file with sequence reads aligned against a nuclear reference genome, you first need to extract the unaligned reads (i.e those reads that likely originate from the mitochondrial genome) from the BAM file using [samtools](http://www.htslib.org), and then convert the BAM file into FASTQ files using [GATK](https://gatk.broadinstitute.org/hc/en-us):
```
samtools view -b -f 4 -o unmapped.reads.bam reads.bam
gatk SamToFastq -I unmapped.reads.bam -F reads1.fq.gz -F2 reads2.fq.gz
```

Once the sequence reads are in FASTQ format, it is recommended to trim sequences for high quality using e.g. [fastp](https://github.com/OpenGene/fastp):
```
fastp -i reads1.fq.gz -I reads2.fq.gz -o reads1.trimmed.fq.gz -O reads2.trimmed.fq.gz -q 30 -u 10 -5 -3 -W 1 -M 30 --cut_right --cut_right_window_size 10 --cut_right_mean_quality 30 -l 100 -b 125
```
You might need to change the setting -l (minimum read length) to a lower value if your reads are shorter.


### 2. Assemble the mitochondrial genome
Use komics all to automate the assembly, circularization and polishing of the mitochondrial minicircles. This can be done using a single command:
```
komics all --kmin 99 --kmax 119 --kstep 10 --minidentity 95 run1 reads1.trimmed.fq.gz reads2.trimmed.fq.gz
```

Minicircle contigs can be found in the file run1.minicircles.fasta. The komics all command will only retain the circularized minicircles, and discard all other contigs. If you wish to keep all minicircle contigs instead, you need to use komics assemble, circularize and polish independently.

The komics all command will generate independent assemblies for each kmer and then merge the minicircles of all assemblies. For instance, for the above settings, the following steps will be performed:
A) Assemble and circularize contigs using a kmer of 99
B) Assemble and circularize contigs using a kmer of 109
C) Assemble and circularize contigs using a kmer of 119
D) Merge circularized minicircles generated in steps A, B and C
E) Reorient minicircles and put CSB1 at the start of the contig
F) Cluster minicircles that are 95% identical

For some datasets it might be better to assemble contigs using a k-mer sweep approach. You can do this using the komics assemble command, which will optimize a single assembly by sweeping through the different kmers (e.g. 99, 109 and 119). After running komics assemble, you can run komics circularize and polish.

The komics all (and assemble) command will also extract maxicircle contigs. Note that the optimal kmer might be different for maxicircles and minicircles. For minicircles, we recommend using high kmer values that are close to the read length (e.g. kmer of 119 for reads that are 125 bp long). For maxicircles, we recommend using lower kmer values (e.g. 29). It is recommended to try different kmers and verify the output files *.maxicircle.fasta and *.minicircles.fasta for each run. For minicircles, you are probably interested in using a kmer strategy that yields the largest number of circularized minicircles. For maxicircles, you are rather interested in the kmer strategy that yields the longest maxicircle contig.


## Reading
This paper includes a detailed outline on how to assemble, circularize and annotate maxicircles based on homology:
__Mitonuclear Genomics Challenges the Theory of Clonality in Trypanosoma Congolense: Reply to Tibayrenc and Ayala__
Van den Broeck et al. Molecular Ecology doi: [10.1111/mec.14809](https://pubmed.ncbi.nlm.nih.gov/30142241/)

This paper includes a detailed outline on how to study mitochondrial genome complexity in whole-genome sequencing datasets:
__Ecological divergence and hybridization of Neotropical Leishmania parasites__
Van den Broeck et al. BIORXIV 2019 doi: [10.1101/824912](https://www.biorxiv.org/content/10.1101/824912v1)

This paper describes the assembly and annotation of the mitochondrial genome in T. brucei, including characterization of the U indel editing patterns:
__Assembly and annotation of the mitochondrial minicircle genome of a differentiation-competent strain of Trypanosoma brucei__
Cooper et al. Nucleic Acids Research 2019 doi: [10.1093/nar/gkz928](https://academic.oup.com/nar/article/47/21/11304/5609525)


## Citation
If you use this software please cite:

__Ecological divergence and hybridization of Neotropical Leishmania parasites__   
Van den Broeck et al. BIORXIV 2019 doi: [10.1101/824912](https://www.biorxiv.org/content/10.1101/824912v1)

