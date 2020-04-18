'''
TO DO:
  * add blast to find maxicircle seqs
  * add histogram plotting showing minicircle lengths
'''
import re
import os
import sys
import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Error (Exception): pass


class Megahit:
  def __init__(self, 
    out,
    reads1,
    reads2,
    threads=1,
    kmin=39,
    kmax=119,
    kstep=10,
    length=400,
    ):
      self.out = out
      self.csb3_contigs = os.path.abspath('tmp.' + self.out + '.csb3contigs.fasta')
      self.other_contigs = os.path.abspath('tmp.' + self.out + '.othercontigs.fasta')
      self.reads1 = os.path.abspath(reads1)
      self.reads2 = os.path.abspath(reads2)
      self.threads = threads
      self.kmin = kmin
      self.kmax = kmax
      self.kstep = kstep
      self.length = length
      self.CSB3 = 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC'

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


#  def _build_CSB3(self, CSB3):
#    if CSB3 is None:
#      return 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC'
#    else:
#      return 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTGATGT|ACATCAACCCC' + '|' + CSB3 + '|' + self._rev_comp(CSB3)


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
    for contig in SeqIO.parse("tmp." + self.out + "_megahit/" + self.out + ".contigs.fa", "fasta"):
      if re.findall(self.CSB3, str(contig.seq)):
        Nminicircles=Nminicircles+1
        minicircles.append(SeqRecord(contig.seq, id=contig.id, description=""))
      else:
        othercontigs.append(SeqRecord(contig.seq, id=contig.id, description=""))

    SeqIO.write(minicircles, self.csb3_contigs, "fasta")
    SeqIO.write(othercontigs, self.other_contigs, "fasta")

    sys.stderr.write('\nFound %s minicircle contigs. They were written to:\n%s\n' % (Nminicircles, os.path.abspath(self.csb3_contigs)))
    sys.stderr.write('\nkomics assembly successfully completed.\n\n')


  def run(self):
    sys.stderr.write('\nAssembling contigs\n')
    sys.stderr.write('==================\n')
    self.runmega()
    self.process_assembly()