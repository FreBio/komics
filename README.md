### Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Tutorial](#tutorial)
  * [Citation](#citation)
  * [References](#references)


### Introduction
KOMICS is a python3.7 package that aids in the assembly and circularization of mitochondrial genomes in trypanosomatids (Van den Broeck et al. 2019). The input is reads in FASTQ format, and the output is maxicircle and circularized minicircles in FASTA format.

Please report any issues or questions to fvandenbroeckATitg.be


### Installation
KOMICS has the following dependencies, which need to be installed first:
  * [MEGAHIT](http://www.metagenomics.wiki/tools/assembly/megahit)
  * [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  * [VSEARCH](https://github.com/torognes/vsearch)

Once the dependencies are installed, install the latest version of KOMICS using pip3:
```
pip3 install git+https://github.com/FreBio/komics.git
```

If you are running KOMICS on a supercomputer, you may want to have a look [here](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/software/python_package_management.html#alternatives-to-conda) on how to setup your local environment and pip install options.


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
If this is the first time you use KOMICS, you may want to follow the tutorial using the data provided on our github page. These sequence reads were generated using whole genome sequencing of *Leishmania peruviana* isolate LCA04. The following files contain the sequence reads that did not align to the *Leishmania braziliensis* M2904 reference genome (see step 1 below).
```
wget https://github.com/FreBio/komics/blob/master/testdata/LCA04_1.fq.gz
wget https://github.com/FreBio/komics/blob/master/testdata/LCA04_2.fq.gz
```


#### 1. Extract unaligned reads from BAM file
If you have a BAM file with sequence reads aligned against a nuclear reference genome, you first need to extract the unaligned reads (i.e those reads that likely originate from the mitochondrial genome) from the BAM file using [samtools](http://www.htslib.org), and then convert the BAM file into FASTQ files using [GATK](https://gatk.broadinstitute.org/hc/en-us):
```
samtools view -b -f 4 -o unmapped.reads.bam reads.bam
gatk SamToFastq -I unmapped.reads.bam -F reads1.fq.gz -F2 reads2.fq.gz
```


#### 2. Trim reads for high-quality
Once the sequence reads are in FASTQ format, it is recommended to trim sequences for high quality using e.g. [fastp](https://github.com/OpenGene/fastp):
```
fastp -i LCA04_1.fq.gz -I LCA04_1.fq.gz -o LCA04_trim_1.fq.gz -O LCA04_trim_2.fq.gz -q 30 -u 10 -5 -3 -W 1 -M 30 --cut_right --cut_right_window_size 10 --cut_right_mean_quality 30 -l 100 -b 125
```
You might need to change the setting -l (minimum read length) to a lower value if your reads are shorter, but we recommend to keep -b (maximum read length) to 125bp because Illumina sequencing quality decreases with increasing number of cycles (i.e. longer reads).


#### 3. Assemble the mitochondrial genome
Use `komics assemble` to assemble the mitochondrial maxicircle and minicircles using MEGAHIT:
```
komics assemble --threads 4 --kmin 29 --kmax 119 --kstep 20 LCA04_run1 LCA04_trim_1.fq.gz LCA04_trim_2.fq.gz
grep '>' LCA04_run1.maxicircles.fasta
grep -c '>' tmp.LCA04_run1.csb3contigs.fasta
```
The optimal kmer length depends on the complexity of the mitochondrial genome, and it will also be different for maxicircle and minicircles. For **minicircles**, we recommend using high kmer values that are close to the read length (e.g. kmer of 119 for reads that are 125 bp long). For **maxicircles**, we recommend using low kmer values (e.g. 29). In the example above, we follow a kmer sweep strategy, whereby MEGAhit is run using the following k-list: 29,49,69,89,109,119.

**Maxicircle** contigs are identified using BLAST against maxicircle sequences of Leishmania braziliensis, Trypanosoma brucei, T. equiperdum and T. lewisi. The resulting maxicircle contigs can be found in the file LCA04\_run1.maxicircles.fasta including one contig of length 17,202bp that covers the entire coding region. If you are interested in generating a complete circularized assembly of the maxicircle that includes both coding and divergent region, please read Van den Broeck et al. (2018).

**Minicircle** contigs are extracted based on the presence of the Conserved Sequence Block 3 (CSB3), a 12-bp minicircle motif, also known as the universal minicircle sequence, that is highly conserved across all Kinetoplastida species (Ray 1989). By default, KOMICS uses the known CSB-3 motif GGGGTTGGTGTA, GGGGTTGATGTA and their reverse complements to extract contigs of putative minicircle origin. The resulting minicircle contigs can be found in the file tmp.LCA04\_run1.csb3contigs.fasta including a total of 53 minicircle contigs. This seems like a rather low number of minicircles, and is due to the fact that we included low kmer values. Let's do another assembly, this time using only high kmer values:
```
komics assemble --threads 4 --kmin 99 --kmax 119 --kstep 10 LCA04_run2 LCA04_trim_1.fq.gz LCA04_trim_2.fq.gz
grep -c '>' tmp.LCA04_run2.csb3contigs.fasta
```
This yields better results as we have now assembled a total of 85 minicircle contigs.


#### 4. Circularize the mitochondrial minicircles
Once we are happy with the number of assembled minicircle contigs, we need to circularize the minicircle contigs. KOMICS uses BLAST as a strategy to identify a sequence that is in common at the start and the end of a given minicircle contig. MEGABLAST is run on the entire set of minicircle contigs with the low complexity filter turned off and allowing a maximum e-value of 10-5. The BLAST output is processed to retain only hits among the same minicircle contig (avoiding artificial dimers) with 100% identity and a minimum 20bp overlap at the start and end of a given contig. Whenever an overlap is found, the contig is classified as circular and the duplicated sequence at the start of the contig is removed.
```
komics circularize LCA04_run2 tmp.LCA04\_run2.csb3contigs.fasta
grep -c '>' tmp.LCA04_run2.circularized.fasta
grep -c '_circularized' tmp.LCA04_run2.circularized.fasta
```
Of the 85 minicircle contigs, 74 (87%) were successfully circularized.


#### 5. Polish the circularized minicircles
Finally, we want to align all minicircles by **(a)** reorienting each minicircle contig based on the CSB3-mer, **(b)** putting the Conserved Sequence Block 1 (CSB1) at the start of each circularized minicircle contig and **(c)** cluster contigs based on a minimum percent identity (e.g. 97%) using VSEARCH.
```
komics polish --minidentity 97 LCA04_run2 tmp.LCA04_run2.circularized.fasta
grep -c '>' LCA04_run2.minicircles.fasta 
```


#### 5. Remove intermediate files
Once you are happy with the final set of maxicircles and minicircles, you can remove all intermediate files:
```
rm -r tmp.LCA04_run*
```


#### Automated assembly and circularization
Use `komics all` to automate the assembly, circularization and polishing of the mitochondrial minicircles. This can be done using a single command, for instance:
```
komics all --kmin 99 --kmax 119 --kstep 10 --minidentity 95 run3 LCA04_trim_1.fq.gz LCA04_trim_2.fq.gz
```

Note that `komics all` will only retain circularized minicircles, and discard all other contigs. This is an important difference with `komics assemble` that retains both circularized and non-circularized minicircle contigs. 
Another important difference is that `komics all` will generate independent assemblies for each kmer (here 99, 109 and 119) and then merge the minicircles of all assemblies. For instance, for the above settings, the following steps will be performed:
1. Assemble and circularize contigs using a kmer of 99
1. Assemble and circularize contigs using a kmer of 109
1. Assemble and circularize contigs using a kmer of 119
1. Merge circularized minicircles generated during the first three steps
1. Reorient circularized minicircles and put the Conserved Sequence Block 1 (CSB1) at the start of the contig
1. Cluster minicircles that are 95% identical


### Citation
If you use this software please cite:

__Ecological divergence and hybridization of Neotropical Leishmania parasites__
Van den Broeck et al. BIORXIV 2019 doi: [10.1101/824912](https://www.biorxiv.org/content/10.1101/824912v1)


### References
Ray. 1989 [Conserved sequence blocks in kinetoplast minicircles from diverse species of trypanosomes](https://dx.doi.org/10.1128%2Fmcb.9.3.1365) Molecular and Cellular Biology.

Van den Broeck et al. 2018 [Mitonuclear Genomics Challenges the Theory of Clonality in Trypanosoma Congolense: Reply to Tibayrenc and Ayala](https://pubmed.ncbi.nlm.nih.gov/30142241/) Molecular Ecology.

Van den Broeck et al. 2019 [Ecological divergence and hybridization of Neotropical Leishmania parasites](https://www.biorxiv.org/content/10.1101/824912v1) BIORXIV.