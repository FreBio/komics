'''
TODO:
  * consider removing the extension of minicircles, if they start anyway with CSB1?
  * include coverage plots pdf
'''

from __future__ import division

import os
import re
import sys
import pysam
import numpy
import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Error (Exception): pass

class Tests:
  def __init__(self,
    out,
    fasta,
    reads1,
    reads2,
    threads,
    cigar,
    ):
      self.out = out
      self.input_contigs = os.path.abspath(fasta)
      self.reads1 = os.path.abspath(reads1)
      self.reads2 = os.path.abspath(reads2)
      self.threads = threads
      self.reffile = 'tmp.' + self.out + ".extended.fasta"
      self.indexfile = str(self.reffile) + ".k8s2"
      self.contigstats = self.out + ".contigstats.txt"
      self.overalstats = self.out + ".overalstats.txt"
      self.sam = 'tmp.' + self.out + ".sam"
      self.bam = 'tmp.' + self.out + ".bam"
      self.CSB3 = 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC'
      
      # for internal use only
      self.extend = int(cigar)
      self.samflags = [81, 83, 161, 163, 97, 99, 145, 147]
      self.MQ = 20
      self.cigar = str(cigar) + 'M'
      
      if not os.path.exists(self.input_contigs):
        sys.stderr.write('\nERROR: contigs file not found: "' + self.input_contigs + '"\n')
        sys.exit(0)
      if not os.path.exists(self.reads1):
        sys.stderr.write('\nERROR: reads1 file not found: "' + self.reads1 + '"\n')
        sys.exit(0)
      if not os.path.exists(self.reads2):
        sys.stderr.write('\nERROR: reads2 file not found: "' + self.reads2 + '"\n')
        sys.exit(0)


  def _rev_comp(self, seq):
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    rc = "".join([comp[x] for x in seq[::-1]])
    return rc


  def extend_fasta(self):
    tmp_result=[]    
    for minicircle in SeqIO.parse(self.input_contigs, "fasta"):
      if re.search('circular', minicircle.id):
        newseq = minicircle.seq[(len(minicircle.seq)-self.extend):] + minicircle.seq
        tmp_result.append(SeqRecord(newseq, id=minicircle.id, description=""))
      else:
        tmp_result.append(SeqRecord(minicircle.seq, id=minicircle.id, description=""))
    
    SeqIO.write(tmp_result, self.reffile, "fasta")


  def smalt_run(self):
  
    smalt_index_command = [
      "smalt index",
      "-k 8",
      "-s 2", 
      self.indexfile,
      self.reffile
      ]
    sys.stderr.write('Creating smalt index file:\n')
    subprocess.call(' '.join(smalt_index_command), shell = True)

    smalt_map_command = [
      "smalt map",
      "-i", str(1500),
      "-y", str(0.95),
      "-x",
      "-r", str(0),
      "-n", str(self.threads),
      "-o", self.sam,
      self.indexfile,
      self.reads1,
      self.reads2
      ]
    sys.stderr.write('\n\nRunning smalt map:\n')
    subprocess.call(' '.join(smalt_map_command), shell = True)
    pysam.sort("-o", self.bam, self.sam)
    pysam.index(self.bam)


  def read_stats(self):
    N=0
    N_mapped=0
    N_properpair=0
    N_MQ20=0
    N_CSB3=0
    N_CSB3_mapped=0
    N_CSB3_mapped_pm=0
    N_CSB3_mapped_pp=0
    
    samfile=pysam.AlignmentFile(self.bam, "rb")
    outfile=samfile.header['SQ']
    
    for ref in list(range(0,len(outfile))):
      outfile[ref]['Nmapped']=0
      outfile[ref]['NMQ20']=0
      outfile[ref]['Nproperpair']=0
      outfile[ref]['NCSB3']=0
      outfile[ref]['NCSB3pm']=0
      outfile[ref]['NCSB3pp']=0
      

    sys.stderr.write('\n\nEstimating read counts.\n')
    for read in samfile.fetch(until_eof=True):
      N = N+1
      if re.findall(self.CSB3, str(read.seq)):
        N_CSB3 = N_CSB3 + 1
      if not read.is_unmapped:
        N_mapped = N_mapped+1
        outfile[read.reference_id]['Nmapped'] = outfile[read.reference_id]['Nmapped']+1
        if read.mapping_quality >= self.MQ:
          N_MQ20 = N_MQ20+1
          outfile[read.reference_id]['NMQ20'] = outfile[read.reference_id]['NMQ20']+1
        if read.flag in self.samflags:
          N_properpair = N_properpair+1
          outfile[read.reference_id]['Nproperpair'] = outfile[read.reference_id]['Nproperpair']+1
        if re.findall(self.CSB3, str(read.seq)):
          N_CSB3_mapped = N_CSB3_mapped+1
          outfile[read.reference_id]['NCSB3'] = outfile[read.reference_id]['NCSB3']+1
        if re.findall(self.CSB3, str(read.seq)) and read.cigarstring == self.cigar:
          N_CSB3_mapped_pm = N_CSB3_mapped_pm+1
          outfile[read.reference_id]['NCSB3pm'] = outfile[read.reference_id]['NCSB3pm']+1
        if re.findall(self.CSB3, str(read.seq)) and read.flag in self.samflags:
          N_CSB3_mapped_pp = N_CSB3_mapped_pp+1
          outfile[read.reference_id]['NCSB3pp'] = outfile[read.reference_id]['NCSB3pp']+1

    samfile.close()
    samfile=pysam.AlignmentFile(self.bam, "rb")
    
    with open(self.overalstats, 'w') as f:
      f.write("Number of reads: %s\n" % N)
      f.write("Number of mapped reads: %s\n" % N_mapped)
      f.write("Number of reads w/ MQ>20: %s\n" % N_MQ20)
      f.write("Number of proper pairs: %s\n" % N_properpair)
      f.write("Number of CSB3 reads: %s\n" % N_CSB3)
      f.write("Number of mapped CSB3 reads: %s\n" % N_CSB3_mapped)
      f.write("Number of perfectly matched CSB3 reads: %s\n" % N_CSB3_mapped_pm)
      f.write("Number of proper paired CSB3 reads: %s\n" % N_CSB3_mapped_pp)
    
    
    sys.stderr.write('Estimating read depths.\n')
    
    for contig in list(range(0, len(outfile))):
      start_pos = 0
#      break_pos = outfile[contig]['LN']
      Y = list()
#      Z = 0
      for pileup in samfile.pileup(outfile[contig]['SN'], max_depth = 5000000):
        while pileup.pos != start_pos:
          Y.append(0)
          start_pos = start_pos+1
#        if pileup.pos == break_pos:
#          for pread in pileup.pileups:
#            if pread.alignment.cigarstring == self.cigar:
#              Z=Z+1
              #pread.alignment.reference_start, pread.alignment.reference_end
        Y.append(pileup.n)
        start_pos = start_pos+1
      
      try:
        outfile[contig]['meandepth'] = numpy.mean(Y)
        outfile[contig]['mediandepth'] = numpy.median(Y)
        outfile[contig]['mindepth'] = numpy.min(Y)
        outfile[contig]['maxdepth'] = numpy.max(Y)
#        outfile[contig]['breakdepth'] = Z
      except ValueError:
        sys.stderr.write('WARNING: contig %s has zero read depth.\n' % (outfile[contig]['SN']))

    samfile.close()
    
    with open(self.contigstats, 'w') as f:
      f.write("SN, LN, Nmapped, Nmq20, Npp, NCSB3, NCSB3pm, NCSB3pp, meandepth, mediandepth, mindepth, maxdepth\n")
      for item in outfile:
        try:
          f.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n" % ( \
            item['SN'], \
            item['LN'], \
            item['Nmapped'], \
            item['NMQ20'], \
            item['Nproperpair'], \
            item['NCSB3'], \
            item['NCSB3pm'], \
            item['NCSB3pp'], \
            round(item['meandepth']), \
            item['mediandepth'], \
            item['mindepth'], \
            item['maxdepth'])
            )
        except KeyError:
          sys.stderr.write('WARNING: contig %s has zero read depth and was not included in output.\n' % item['SN'])
    
    sys.stderr.write('Overal read counts were written to %s\n' % (os.path.abspath(self.overalstats)))
    sys.stderr.write('Contig read counts and depths were written to %s\n' % (os.path.abspath(self.contigstats)))
  
  def run(self):
    sys.stderr.write('\nEstimating minicircle read depths and counts\n')
    sys.stderr.write('============================================\n')
    self.extend_fasta()
    self.smalt_run()
    self.read_stats()
    sys.stderr.write('\nkomics qualcheck successfully completed.\n\n')
