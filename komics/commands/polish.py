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
    parser.add_argument('--csb1', type=str, help='Specificy one or more Conserved Sequence Block 1 (CSB1) sequences, which will put at the start of the alignments. When nothing is provided, the following two CSB1-mers are used by default: "GGGCGTTC|GGGCGTGC".', metavar='STR')
    parser.add_argument('--csb3', type=str, help='Specificy one or more Conserved Sequence Block 3 (CSB3) sequences, used for reorienting the minicircles. When nothing is provided, the following two CSB3-mers are used by default: "GGGGTTGGTGT|GGGGTTGATGT".', metavar='STR')
    options = parser.parse_args()
    
    args = vars(parser.parse_args())
    if args["csb1"] is None:
        csb1mer = 'GGGCGTTC|GGGCGTGC'
    elif args["csb1"] is not None:
        csb1mer = args["csb1"]

    if args["csb3"] is None:
        csb3mer = 'GGGGTTGGTGT|GGGGTTGATGT'
    elif args["csb3"] is not None:
        csb3mer = args["csb3"]

    kp = komics.polish.Polisher(
		options.out,
		options.fasta,
		options.minidentity,
		csb1mer,
		csb3mer,
	)
    kp.run()