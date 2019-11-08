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
#	parser.add_argument('--CSB3', help='Additional CSB3 motif [%(default)s]. \
#	                              The conserved CSB3 motifs are used by komics to find minicircles in the pool of assembled contigs. \
#	                              The default CSB3 motifs are GGGGTTGGTGT and GGGGTTAGTGT (and their reverse complements). \
#	                              With the "CSB3" option, the user can specify one additional motif if needed.', default=None, metavar='')
	options = parser.parse_args()
	
	ka = komics.assemble.Megahit(
		options.out,
		options.reads1,
		options.reads2,
		options.threads,
		options.kmin,
		options.kmax,
		options.kstep,
		options.length,
	)
	ka.run()
