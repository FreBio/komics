import sys
import komics
import argparse

def main():
	parser = argparse.ArgumentParser(
  	prog='komics',
		description = 'Circularize contigs using BLAST search for overlapping ends',
		usage = 'komics circularize [options] <out> <fasta>')
	parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
	parser.add_argument('--minoverlap', type=int, help='Minimum overlap (bp) between the contig ends [%(default)s]', default=20, metavar='INT')
	parser.add_argument('out', help='Prefix used for labeling FASTA output file and sequence headers', metavar='out')
	parser.add_argument('fasta', help='FASTA file w/ contigs', metavar = 'fasta')
	options = parser.parse_args()
	
	kc = komics.circularize.Circularizer(
		options.out,
		options.fasta,
		threads=options.threads,
		minoverlap=options.minoverlap,
	)
	kc.run()
	