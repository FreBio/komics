import re
import os
import sys
import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline


class Error (Exception): pass


class Megahit:
  def __init__(self, 
    out,
    reads1,
    reads2,
    csb3mer,
    threads,
    kmin,
    kmax,
    kstep,
    length,
    ):
      self.out = out
      self.input_contigs = "tmp." + self.out + "_megahit/" + self.out + ".contigs.fa"
      self.csb3_contigs = os.path.abspath('tmp.' + self.out + '.csb3contigs.fasta')
      self.other_contigs = os.path.abspath('tmp.' + self.out + '.othercontigs.fasta')
      self.maxi_xml = os.path.abspath('tmp.' + self.out + '.maxicircle.xml')
      self.input_ref_maxicircles = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'maxicircles.fasta')
      self.output_max_contigs = os.path.abspath(self.out + '.maxicircles.fasta')
      self.reads1 = os.path.abspath(reads1)
      self.reads2 = os.path.abspath(reads2)
      self.CSB3 = csb3mer
      self.threads = threads
      self.kmin = kmin
      self.kmax = kmax
      self.kstep = kstep
      self.length = length
      #self.CSB3 = 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC'

      if not os.path.exists(self.reads1):
        sys.stderr.write('\nERROR: reads1 file not found: "' + self.reads1 + '"\n')
        sys.exit(0)
      if not os.path.exists(self.reads2):
        sys.stderr.write('\nERROR: reads2 file not found: "' + self.reads2 + '"\n')
        sys.exit(0)


  def _rev_comp(self, seq):
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "|": "|"}
    rc = "".join([comp[x] for x in seq[::-1]])
    return rc


  def runmega(self):
    megahit_cmd = [
      "megahit",
      "-1", self.reads1,
      "-2", self.reads2,
      "-t", str(self.threads),
      "-m 4",
      "-o", "tmp." + self.out + "_megahit",
      "--min-contig-len", str(self.length),
      "--out-prefix", self.out,
      "--k-min", str(self.kmin),
      "--k-max", str(self.kmax),
      "--k-step", str(self.kstep)
      ]

    subprocess.call(' '.join(megahit_cmd), shell = True)
    

  def process_assembly(self):
    minicircles=[]
    othercontigs=[]
    Nminicircles = 0
    for contig in SeqIO.parse(self.input_contigs, 'fasta'):
      if re.findall(str(self.CSB3+"|"+self._rev_comp(self.CSB3)), str(contig.seq)):
        Nminicircles=Nminicircles+1
        minicircles.append(SeqRecord(contig.seq, id=contig.id, description=""))
      else:
        othercontigs.append(SeqRecord(contig.seq, id=contig.id, description=""))

    SeqIO.write(minicircles, self.csb3_contigs, "fasta")
    SeqIO.write(othercontigs, self.other_contigs, "fasta")

    sys.stderr.write('\nFound %s minicircle contigs. They were written to:\n%s\n' % (Nminicircles, os.path.abspath(self.csb3_contigs)))
    sys.stderr.write('\nkomics assemble successfully completed.\n\n')
    
    
  def blast_run(self):
      blastcmd = NcbiblastnCommandline(
              task='megablast',
              query=self.input_contigs,
              subject=self.input_ref_maxicircles,
              out=self.maxi_xml,
              dust="no",
              evalue=0.000001,
              outfmt=5)
      sys.stderr.write('Running BLAST with the following command:\n')
      sys.stderr.write(str(blastcmd)+ '\n\n')
      stdout, stderr=blastcmd()
  
    
  def blast_parse(self):
      sys.stderr.write('Parsing BLAST output:\n')
      
      records=[]
      for record in NCBIXML.parse(open(self.maxi_xml)):
          for alignment in record.alignments:
              records.append(re.split("\s+", record.query)[0])
      return records


  def maxicircle(self):
      maxicircles=self.blast_parse()
      maxicirclesout=[]
      
      for contig in SeqIO.parse(self.input_contigs, "fasta"):
          if contig.id in maxicircles:
              newseq=contig.seq
              newid=self.out + '_contig' + contig.id + '_len' + str(len(newseq))
              maxicirclesout.append(SeqRecord(newseq, id=newid, description=""))
      SeqIO.write(maxicirclesout, self.output_max_contigs, "fasta")
      sys.stderr.write('Maxicircle contigs can be found in %s.\n ' % (self.output_max_contigs))


  def run(self):
    sys.stderr.write('\nAssembling contigs\n')
    sys.stderr.write('==================\n')
    self.runmega()
    sys.stderr.write('\nExtracting maxicircle and minicircle contigs\n')
    sys.stderr.write('============================================\n')
    self.blast_run()
    self.maxicircle()
    self.process_assembly()
		