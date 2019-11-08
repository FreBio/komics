'''
TODO:
  * put in some checks
'''

import argparse
import os
import sys
import komics

class Error (Exception): pass

def main():
    parser = argparse.ArgumentParser(
        description = 'Runs bamfq, trimfq, assemble, circularize, qualcheck',
        usage = 'komics all [options] <out> <bam>')
    parser.add_argument('out', help='Prefix used for labeling files and sequences', metavar='out')
    parser.add_argument('bam', help='BAM file of reads aligned to reference genome', metavar='bam')
    parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    parser.add_argument('--kmin', help='Minimum k-mer (must be odd number) [%(default)s]', default=29, metavar='INT')
    parser.add_argument('--kmax', help='Maximum k-mer (must be odd number) [%(default)s]', default=119, metavar='INT')
    parser.add_argument('--kstep', help='Steps between k-mers (must be even number) [%(default)s]', default=10, metavar='INT')
    parser.add_argument('--length', help='Minimum length (bp) of contigs to be kept [%(default)s]', default=400, metavar='INT')
    parser.add_argument('--CSB3', help='Additional CSB3 motif [%(default)s].', default=None, metavar='')
    parser.add_argument('--minoverlap', type=int, help='Minimum overlap (bp) between the contig ends [%(default)s]', default=20, metavar='INT')
    parser.add_argument('--word', type=int, help='Specifies the word length for smalt index file [%(default)s]', default=8, metavar='INT')
    parser.add_argument('--step', type=int, help='Specifies how many bases are skipped between indexed words for smalt index file [%(default)s]', default=2, metavar='INT')
    parser.add_argument('--minidentity', type=int, help='Minimum percent identity between minicircles [%(default)s]', default=95, metavar='INT')
    parser.add_argument('--cigar', type=int, help='Specifies the read length [%(default)s]', default=100, metavar='INT')
    options = parser.parse_args()


    # ------ Running assemble -------
    ka = komics.assemble.Megahit(
      options.out,
      options.treads1,
      options.treads2,
      options.threads,
      options.kmin,
      options.kmax,
      options.kstep,
      options.length,
      options.CSB3,
    )
    ka.run()

    # ------ Running circularize ----
    contigs=os.path.abspath('tmp.' + options.out + '.csb3contigs.fasta')

    kc = komics.circularize.Circularizer(
      options.out,
      contigs,
      threads=options.threads,
      minoverlap=options.minoverlap,
    )
    kc.run()
    
    # ------ Running polish ----
    contigs=os.path.abspath('tmp.' + options.out + '.circularized.fasta')
    
    kp = komics.polish.Polisher(
      options.out,
      contigs,
      options.minidentity
    )
    kp.run()
    
    # ------ Running qualcheck ------
    circularcontigs=os.path.abspath(options.out + '.minicircles.fasta')
    
    kq = komics.qualcheck.Tests(
      options.out,
      circularcontigs,
      options.mreads1,
      options.mreads2,
      options.threads,
      options.word,
      options.step,
      options.cigar
    )
    kq.run()

    # ------ Cleaning ---------------
    # rm all tmp files