import sys
import komics
import argparse

def main():
	parser = argparse.ArgumentParser(
	  prog='komics',
		description='Assembles contigs using MEGAHIT',
		usage = 'komics assemble [options] <out> <reads1> <reads2>')
	parser.add_argument('out', help='Prefix used for labeling FASTQ files', metavar='out')
	parser.add_argument('reads1', help='FASTQ file w/ first-in-pair reads', metavar='reads1')
	parser.add_argument('reads2', help='FASTQ file w/ second-in-pair reads', metavar='reads1')
	parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
	parser.add_argument('--kmin', help='Minimum k-mer (must be odd number) [%(default)s]', default=39, metavar='INT')
	parser.add_argument('--kmax', help='Maximum k-mer (must be odd number) [%(default)s]', default=119, metavar='INT')
	parser.add_argument('--kstep', help='Steps between k-mers (must be even number) [%(default)s]', default=10, metavar='INT')
	parser.add_argument('--length', help='Minimum length (bp) of contigs to be kept [%(default)s]', default=400, metavar='INT')
	parser.add_argument('--csb3', type=str, help='Specificy one or more Conserved Sequence Block 3 (CSB3) sequences, used for identifying minicircles. When nothing is provided, the following two CSB3-mers are used by default: "GGGGTTGGTGT|GGGGTTGATGT".', metavar='STR')
	options = parser.parse_args()
	
	args = vars(parser.parse_args())
	if args["csb3"] is None:
		csb3mer = 'GGGGTTGGTGT|GGGGTTGATGT'
	elif args["csb3"] is not None:
		csb3mer = args["csb3"]
	
	ka = komics.assemble.Megahit(
		options.out,
		options.reads1,
		options.reads2,
		csb3mer,
		options.threads,
		options.kmin,
		options.kmax,
		options.kstep,
		options.length,
	)
	ka.run()
