### Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Tutorial](#tutorial)
  * [Reading](#reading)
  * [Citation](#citation)


### Introduction
A tool for automated assembly and circularization of mitochondrial genomes in trypanosomatids. The input is reads in FASTQ format, and the output is maxicircle and circularized minicircles in FASTA format.

komics is described in detail here:
__Ecological divergence and hybridization of Neotropical Leishmania parasites__
Van den Broeck et al. BIORXIV 2019 doi: [10.1101/824912](https://www.biorxiv.org/content/10.1101/824912v1)

Please report any issues to fvandenbroeckATitg.be


### Installation
KOMICS has the following dependencies, which need to be installed first:
  * [MEGAHIT](http://www.metagenomics.wiki/tools/assembly/megahit)
  * [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  * [VSEARCH](https://github.com/torognes/vsearch)

Once the dependencies are installed, install the latest version of KOMICS using pip3:
```
pip3 install git+https://github.com/FreBio/komics.git
```

If you are running komics on a supercomputer, you may want to have a look [here](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/software/python_package_management.html#alternatives-to-conda) on how to setup your local environment and pip install options.


### Usage
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


### Tutorial

#### Prepare your data
If you have a BAM file with sequence reads aligned against a nuclear reference genome, you first need to extract the unaligned reads (i.e those reads that likely originate from the mitochondrial genome) from the BAM file using [samtools](http://www.htslib.org), and then convert the BAM file into FASTQ files using [GATK](https://gatk.broadinstitute.org/hc/en-us):
```
samtools view -b -f 4 -o unmapped.reads.bam reads.bam
gatk SamToFastq -I unmapped.reads.bam -F reads1.fq.gz -F2 reads2.fq.gz
```

Once the sequence reads are in FASTQ format, it is recommended to trim sequences for high quality using e.g. [fastp](https://github.com/OpenGene/fastp):
```
fastp -i reads1.fq.gz -I reads2.fq.gz -o reads1.trimmed.fq.gz -O reads2.trimmed.fq.gz -q 30 -u 10 -5 -3 -W 1 -M 30 --cut_right --cut_right_window_size 10 --cut_right_mean_quality 30 -l 100 -b 125
```
You might need to change the setting -l (minimum read length) to a lower value if your reads are shorter, but we recommend to keep -b (maximum read length) to 125bp because Illumina sequencing quality decreases with increasing number of cycles (i.e. longer reads).


#### Assemble the mitochondrial genome
Use `komics all` to automate the assembly, circularization and polishing of the mitochondrial minicircles. This can be done using a single command, for instance:
```
komics all --kmin 99 --kmax 119 --kstep 10 --minidentity 95 run1 reads1.trimmed.fq.gz reads2.trimmed.fq.gz
```

The resulting circularized minicircle sequences can be found in the file run1.minicircles.fasta.

`komics all` will generate independent assemblies for each kmer (here 99, 109 and 119) and then merge the minicircles of all assemblies. For instance, for the above settings, the following steps will be performed:
1. Assemble and circularize contigs using a kmer of 99
1. Assemble and circularize contigs using a kmer of 109
1. Assemble and circularize contigs using a kmer of 119
1. Merge circularized minicircles generated during the first three steps
1. Reorient circularized minicircles and put the Conserved Sequence Block 1 (CSB1) at the start of the contig
1. Cluster minicircles that are 95% identical

Note that `komics all` will only retain circularized minicircles, and discard all other contigs. If you also wish to keep the non-circularized minicircles, then you need to use `komics assemble`. In addition, for some datasets it might be better to assemble contigs using a k-mer sweep approach. This too can be done with `komics assemble`, which will optimize a single assembly by sweeping through the different kmers (e.g. 99, 109 and 119).
```
komics assemble --kmin 99 --kmax 119 --kstep 10 run1 reads1.trimmed.fq.gz reads2.trimmed.fq.gz
komics circularize run1 tmp.run1.csb3contigs.fasta
komics polish run1 tmp.run1.circularized.fasta
```

Again, the resulting minicircle sequences can be found in the file run1.minicircles.fasta, and will now include both non-circularized and circularized minicircles.

The `komics all` and `komics assemble` commands will also extract the maxicircle contigs using a blast approach. Note that the optimal kmer might be different for maxicircles and minicircles. For minicircles, we recommend using high kmer values that are close to the read length (e.g. kmer of 119 for reads that are 125 bp long). These might work well too for the maxicircle, but usually you need to choose lower kmers (e.g. 39) to obtain longer maxicircle contigs.


### Reading
This paper includes a detailed outline on how to assemble, circularize and annotate maxicircles based on homology:
__Mitonuclear Genomics Challenges the Theory of Clonality in Trypanosoma Congolense: Reply to Tibayrenc and Ayala__
Van den Broeck et al. 2018 Molecular Ecology doi: [10.1111/mec.14809](https://pubmed.ncbi.nlm.nih.gov/30142241/)

This paper includes a detailed outline on how to study mitochondrial genome complexity in whole-genome sequencing datasets:
__Ecological divergence and hybridization of Neotropical Leishmania parasites__
Van den Broeck et al. 2019 BIORXIV doi: [10.1101/824912](https://www.biorxiv.org/content/10.1101/824912v1)

This paper describes the assembly and annotation of the mitochondrial genome in T. brucei, including characterization of the U indel editing patterns:
__Assembly and annotation of the mitochondrial minicircle genome of a differentiation-competent strain of Trypanosoma brucei__
Cooper et al. 2019 Nucleic Acids Research doi: [10.1093/nar/gkz928](https://academic.oup.com/nar/article/47/21/11304/5609525)


### Citation
If you use this software please cite:

__Ecological divergence and hybridization of Neotropical Leishmania parasites__   
Van den Broeck et al. BIORXIV 2019 doi: [10.1101/824912](https://www.biorxiv.org/content/10.1101/824912v1)

