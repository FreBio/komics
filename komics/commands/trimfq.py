import komics
import argparse

def main():
	parser = argparse.ArgumentParser(
	  prog='komics',
		description='Trim FASTQ files using TRIMMOMATIC',
		usage = 'komics trimfq [options] <out> <reads1> <reads2>')
	parser.add_argument('out', help='Prefix used for labeling FASTQ files', metavar='out')
	parser.add_argument('reads1', help='FASTQ file w/ first-in-pair reads', metavar='reads1')
	parser.add_argument('reads2', help='FASTQ file w/ second-in-pair reads', metavar='reads1')
	parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
	parser.add_argument('--dir', help='Directory containing the trimmomatic jar file and adapters. \
	                                  If "None" (default), the path environment variable $DIR_TRIMMOMATIC is used [%(default)s].', default=None , metavar='')
	parser.add_argument('--maxlength', type=int, help='Maximum read length after trimming [%(default)s]', default=150, metavar='INT')
	parser.add_argument('--minlength', type=int, help='Minimum read length after trimming [%(default)s]', default=100, metavar='INT')
	parser.add_argument('--minqual', type=int, help='Minimum (phred-scaled) base quality [%(default)s]', default=30, metavar='INT')
	parser.add_argument('--windowbp', type=int, help='Length (bp) of sliding windows [%(default)s]', default=10, metavar='INT')
	options = parser.parse_args()
	
	kt = komics.trimfq.Trimmomatic(
		options.out,
		options.reads1,
		options.reads2,
		options.threads,
		options.dir,
		options.maxlength,
		options.minlength,
		options.minqual,
		options.windowbp,		
	)
	kt.trim_reads()
