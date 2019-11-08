'''
TODO:
  * put in some checks
'''


import os
import re
import sys
import komics
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from timeit import default_timer as timer

class Error (Exception): pass

def main():
    parser = argparse.ArgumentParser(
        description = 'Runs bamfq, trimfq, assemble, circularize, polish and qualcheck',
        usage = 'komics all [options] <out> <bam>')
    parser.add_argument('out', help='Prefix used for labeling files and sequences', metavar='out')
    parser.add_argument('bam', help='BAM file of reads aligned to reference genome', metavar='bam')
    parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    parser.add_argument('--dir', help='Trimmomatic directory [%(default)s]', default=None , metavar='')
    parser.add_argument('--maxlength', type=int, help='Maximum read length after trimming [%(default)s]', default=150, metavar='INT')
    parser.add_argument('--minlength', type=int, help='Minimum read length after trimming [%(default)s]', default=100, metavar='INT')
    parser.add_argument('--minqual', type=int, help='Minimum (phred-scaled) base quality [%(default)s]', default=30, metavar='INT')
    parser.add_argument('--windowbp', type=int, help='Length (bp) of sliding windows [%(default)s]', default=10, metavar='INT')
    parser.add_argument('--kmin', type=int, help='Minimum k-mer (must be odd number) [%(default)s]', default=89, metavar='INT')
    parser.add_argument('--kmax', type=int, help='Maximum k-mer (must be odd number) [%(default)s]', default=119, metavar='INT')
    parser.add_argument('--kstep', type=int, help='Steps between k-mers (must be even number) [%(default)s]', default=10, metavar='INT')
    parser.add_argument('--length', type=int, help='Minimum length (bp) of contigs to be kept [%(default)s]', default=400, metavar='INT')
    parser.add_argument('--minoverlap', type=int, help='Minimum overlap (bp) between the contig ends [%(default)s]', default=20, metavar='INT')
    parser.add_argument('--word', type=int, help='Specifies the word length for smalt index file [%(default)s]', default=8, metavar='INT')
    parser.add_argument('--step', type=int, help='Specifies how many bases are skipped between indexed words for smalt index file [%(default)s]', default=2, metavar='INT')
    parser.add_argument('--minidentity', type=int, help='Minimum percent identity between minicircles [%(default)s]', default=95, metavar='INT')
    parser.add_argument('--cigar', type=int, help='Specifies the read length [%(default)s]', default=100, metavar='INT')
    options = parser.parse_args()


    start = timer()
#    sys.stderr.write('\nKOMICS initiated at ' + str(start) + '\n')
    
    # ------ Running bamfq ---------
    if not os.path.exists(options.bam):
      sys.stderr.write('\nERROR: BAM file not found: "' + options.bam + '"\n')
      sys.exit(0)
    
    kb = komics.bam2fq.Bam2fq(
      options.out,
      options.bam,
    )
    kb.run()


    # ------ Running trimfq ---------
    nmppd_reads1=os.path.abspath(options.out + '_unmapped_1P.fq.gz')
    nmppd_reads2=os.path.abspath(options.out + '_unmapped_2P.fq.gz')

    kt = komics.trimfq.Trimmomatic(
      options.out,
      nmppd_reads1,
      nmppd_reads2,
      options.threads,
      options.dir,
      options.maxlength,
      options.minlength,
      options.minqual,
      options.windowbp,		
    )
    kt.trim_reads()


    # ------ Running assemble -------
    trimmed_reads1=os.path.abspath(options.out + '_trimmed_1P.fq.gz')
    trimmed_reads2=os.path.abspath(options.out + '_trimmed_2P.fq.gz')
    
    klist=list(range(options.kmin, options.kmax+1, options.kstep))
    files=[]
    for k in klist:
      OUT = str(options.out + '_' + str(k))
      files.append('tmp.' + OUT + '.circularized.fasta')
      
      ka = komics.assemble.Megahit(
        OUT,
        trimmed_reads1,
        trimmed_reads2,
        options.threads,
        int(k),
        int(k),
        int(0),
        options.length,
      )
      ka.run()


      # ------ Running circularize ----
      contigs=os.path.abspath('tmp.' + OUT + '.csb3contigs.fasta')

      kc = komics.circularize.Circularizer(
        OUT,
        contigs,
        minoverlap=options.minoverlap,
      )
      kc.run()
    
    
    # ------ Combining fasta files; only circular sequences are kept ---
    contigs=os.path.abspath('tmp.' + options.out + '.circularized.fasta')
    
    tmp_result=[]
    for f in files:
      for minicircle in SeqIO.parse(f, "fasta"):
        if re.search('circular', minicircle.id):
          tmp_result.append(SeqRecord(minicircle.seq, id=minicircle.id, description=""))    
    SeqIO.write(tmp_result, contigs, "fasta")
    
    
    # ------ Running polish ----
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
      nmppd_reads1,
      nmppd_reads2,
      options.threads,
      options.cigar
    )
    kq.run()

    end = timer()
#    sys.stderr.write('\nKOMICS ended at ' + str(end) + '\n')
    sys.stderr.write('\nKOMICS finished in: ' + str(end - start) + 'seconds\n')

    # ------ Cleaning ---------------
#    for f in glob.glob('*_megahit'):
#      os.rmdir(f)
#    for f in glob.glob('tmp.' + options.out + '*'):
#      os.remove(f)