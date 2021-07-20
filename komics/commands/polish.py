import re
import komics
import argparse

def main():
    parser = argparse.ArgumentParser(
            prog='komics',
            description = 'Reorientates minicircles putting CSB3 at start',
            usage = 'komics polish [options] <out> <fasta>')
    parser.add_argument('out', help='Prefix used for labeling FASTA output file and sequence headers', metavar='out')
    parser.add_argument('fasta', help='FASTA file with minicircle contigs', metavar = 'fasta')
    parser.add_argument('--minidentity', type=int, help='Cluster minicircle contigs based on a minimum percent identity [%(default)s].', default=95, metavar='INT')
    parser.add_argument('--csb1', type=str, help='Specificy one or more Conserved Sequence Block 1 (CSB1) sequences used for reorienting the minicircles. Regular expressions are allowed, e.g. "GG.CGTTC|GGG[C/G]GTTC". When nothing is provided, then the following two CSB1-mers are used: "GGGCGT[T/G]C".', metavar='STR')
    options = parser.parse_args()
    
    args = vars(parser.parse_args())
    if args["csb1"] is None:
        csb1mer = re.compile(r'GGGCGTTC|GGGCGTGC')
    elif args["csb1"] is not None:
        csb1mer = re.compile(args["csb1"])

    kp = komics.polish.Polisher(
		options.out,
		options.fasta,
		options.minidentity,
		csb1mer,
	)
    kp.run()