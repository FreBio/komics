import sys
import komics
import argparse

def main():
	parser = argparse.ArgumentParser(
	  prog='komics',
		description='Write unmapped reads from a BAM file to FASTQ files. \
		Only when both paired-end reads are unmapped, they are written to a FASTQ file. \
		Quality scores are taken as is from the BAM file, using (Phred+33) ASCII representations. \
    Depends on python module "pysam".',
		usage = 'komics bam2fq <out> <bam>')
	parser.add_argument('out', help='Prefix used for labeling FASTQ files', metavar='out')
	parser.add_argument('bam', help='BAM file of reads aligned to reference genome', metavar='bam')	
	options = parser.parse_args()
	
	kb = komics.bam2fq.Bam2fq(
		options.out,
		options.bam,
	)
	kb.run()
