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
        description = 'Runs assemble, circularize and polish',
        usage = 'komics all [options] <out> <reads1> <reads2>')
    parser.add_argument('out', help='Prefix used for labeling files and sequences', metavar='out')
    parser.add_argument('reads1', help='FASTQ file w/ first-in-pair reads', metavar='reads1')
    parser.add_argument('reads2', help='FASTQ file w/ second-in-pair reads', metavar='reads1')
    parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    parser.add_argument('--kmin', type=int, help='Minimum k-mer (must be odd number) [%(default)s]', default=89, metavar='INT')
    parser.add_argument('--kmax', type=int, help='Maximum k-mer (must be odd number) [%(default)s]', default=119, metavar='INT')
    parser.add_argument('--kstep', type=int, help='Steps between k-mers (must be even number) [%(default)s]', default=10, metavar='INT')
    parser.add_argument('--length', type=int, help='Minimum length (bp) of contigs to be kept [%(default)s]', default=400, metavar='INT')
    parser.add_argument('--minoverlap', type=int, help='Minimum overlap (bp) between the contig ends [%(default)s]', default=20, metavar='INT')
    parser.add_argument('--minidentity', type=int, help='Minimum percent identity between minicircles [%(default)s]', default=95, metavar='INT')
    options = parser.parse_args()


    start = timer()
#    sys.stderr.write('\nKOMICS initiated at ' + str(start) + '\n')
    

    # ------ Running assemble -------
    if not os.path.exists(options.reads1):
      sys.stderr.write('\nERROR: BAM file not found: "' + options.reads1 + '"\n')
      sys.exit(0)
    if not os.path.exists(options.reads2):
      sys.stderr.write('\nERROR: BAM file not found: "' + options.reads2 + '"\n')
      sys.exit(0)
    
    klist=list(range(options.kmin, options.kmax+1, options.kstep))
    files=[]
    for k in klist:
      OUT = str(options.out + '_' + str(k))
      files.append('tmp.' + OUT + '.circularized.fasta')
      
      ka = komics.assemble.Megahit(
        OUT,
        options.reads1,
        options.reads2,
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
    
    
    end = timer()
#    sys.stderr.write('\nKOMICS ended at ' + str(end) + '\n')
    sys.stderr.write('\nKOMICS finished in: ' + str(end - start) + 'seconds\n')

    # ------ Cleaning ---------------
#    for f in glob.glob('*_megahit'):
#      os.rmdir(f)
#    for f in glob.glob('tmp.' + options.out + '*'):
#      os.remove(f)
