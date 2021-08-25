'''
TO DO: 
  * Too slow: due to writing to disk. Speed up by writing in batches or adding threads here:
  * adding thread could perhaps be done using the htslib: https://github.com/pysam-developers/pysam/issues/579
  * Include phred-scale check, now only phred 33 acceepted
  * rev comp is not necessary if read is unmapped, but ok
'''

import os
import sys
import gzip
import pysam


class Error (Exception): pass


class Bam2fq:
  def __init__(self, 
	  out,
	  bam, 
	  ):
      self.sample_name = out
      self.bam = os.path.abspath(bam)
      self.fastq1 = os.path.abspath(self.sample_name + '_unmapped_1P.fq.gz')
      self.fastq2 = os.path.abspath(self.sample_name + '_unmapped_2P.fq.gz')
      if not os.path.exists(self.bam):
        sys.stderr.write('\nERROR: BAM file not found: "' + self.bam + '"\n')
        sys.exit(0)

  def write_fq(self, file, seq):
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}  

    if seq.is_reverse:
      sequence_rc = "".join([comp[x] for x in seq.query_sequence[::-1]])
      quality_rc = seq.query_qualities[::-1]
      print((seq.qname, sequence_rc, "".join([chr(x + 33) for x in quality_rc])))
      file.write(str('@' + seq.qname + '\n' + sequence_rc + '\n+\n' + "".join([chr(x + 33) for x in quality_rc]) + '\n').encode())
#      file.write("@%s\n%s\n+\n%s\n".encode() % (seq.qname, sequence_rc, "".join([chr(x + 33) for x in quality_rc])))
    else:
      file.write(str('@' + seq.qname + '\n' + seq.query_sequence + '\n+\n' + "".join([chr(x + 33) for x in seq.query_qualities]) + '\n').encode())
      #print('@' + seq.qname + '\n' + seq.query_sequence + '\n+\n' + "".join([chr(x + 33) for x in seq.query_qualities]) + '\n')
      # print('\n+\n%s\n'.encode() % (seq.qname, seq.query_sequence, "".join([chr(x + 33) for x in seq.query_qualities])))
#      file.write("@%s\n%s\n+\n%s\n".encode() % (seq.qname, seq.query_sequence, "".join([chr(x + 33) for x in seq.query_qualities])))


  def run(self):
    sys.stderr.write('\nWriting unmapped paired-end reads to FASTQ files\n')
    sys.stderr.write('================================================\n')

    samfile = pysam.AlignmentFile(self.bam, "rb")
    fq1=gzip.GzipFile(self.fastq1, "wb")
    fq2=gzip.GzipFile(self.fastq2, "wb")

    number_of_reads=0
    number_of_unmapped_reads=0
    
    for read in samfile.fetch(until_eof=True):
      number_of_reads = number_of_reads+1
      if read.is_unmapped and read.mate_is_unmapped:
        number_of_unmapped_reads = number_of_unmapped_reads+1
        if read.is_read1:
          self.write_fq(fq1, read)
        elif read.is_read2:
          self.write_fq(fq2, read)

      if number_of_reads % 1000000 == 0:
        sys.stderr.write('Processed %i million reads\n' % (number_of_reads/1000000))

    fq1.close()
    fq2.close()
    samfile.close()
    
    sys.stderr.write('\n%i unmapped reads written to the following FASTQ files:\n%s\n%s\n' % (number_of_unmapped_reads, self.fastq1, self.fastq2))
    sys.stderr.write('\nkomics bam2fq successfully completed.\n\n')