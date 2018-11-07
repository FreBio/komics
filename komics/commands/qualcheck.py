import sys
import komics
import argparse

def main():
	parser = argparse.ArgumentParser(
  	prog='komics',
		description = 'Validates assembly and circularization',
		usage = 'komics validate [options] <out> <fasta> <reads1> <reads2>')
	parser.add_argument('out', help='Prefix used for labeling FASTA output file and sequence headers', metavar='out')
	parser.add_argument('fasta', help='FASTA file w/ contigs', metavar = 'fasta')
	parser.add_argument('reads1', help='FASTQ file w/ first-in-pair reads', metavar='reads1')
	parser.add_argument('reads2', help='FASTQ file w/ second-in-pair reads', metavar='reads2')
	parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
	parser.add_argument('--word', type=int, help='Specifies the word length for smalt index file [%(default)s]', default=8, metavar='INT')
	parser.add_argument('--step', type=int, help='Specifies how many bases are skipped between indexed words for smalt index file [%(default)s]', default=2, metavar='INT')
		
	options = parser.parse_args()
	
	kq = komics.qualcheck.Tests(
		options.out,
		options.fasta,
    options.reads1,
		options.reads2,
		options.threads,
		options.word,
		options.step
	)
	kq.run()
	