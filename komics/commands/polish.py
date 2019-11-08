import komics
import argparse

def main():
    parser = argparse.ArgumentParser(
            prog='komics',
            description = 'Reorientates minicircles putting CSB3 at start',
            usage = 'komics polish [options] <out> <fasta>')
    parser.add_argument('out', help='Prefix used for labeling FASTA output file and sequence headers', metavar='out')
    parser.add_argument('fasta', help='FASTA file w/ contigs', metavar = 'fasta')
    parser.add_argument('--minidentity', type=int, help='Minimum percent identity between minicircles [%(default)s]', default=95, metavar='INT')
    options = parser.parse_args()
    
    kp = komics.polish.Polisher(
		options.out,
		options.fasta,
		options.minidentity,
	)
    kp.run()