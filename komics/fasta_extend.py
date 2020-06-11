#!/usr/bin/env python


#%% MODULES TO IMPORT 

from __future__ import division
import os
import re
import sys
import numpy
import argparse
import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main(output, input, length):
    tmp_result=[]    
    for minicircle in SeqIO.parse(input, "fasta"):
      if re.search('circular', minicircle.id):
        newseq = minicircle.seq[(len(minicircle.seq)-int(length)):] + minicircle.seq
        tmp_result.append(SeqRecord(newseq, id=minicircle.id, description=""))
      else:
        tmp_result.append(SeqRecord(minicircle.seq, id=minicircle.id, description=""))
    
    SeqIO.write(tmp_result, output, "fasta")

#%% MAIN
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
      description='Extends circular minicircles by a given length for proper mapping.', \
      usage = 'fasta_extend.py <output> <input> <length>')
  parser.add_argument('output', help='Name of the output file', metavar='vcf')
  parser.add_argument('input', help='Name of the input file', metavar='vcf')
  parser.add_argument('length', help='Length in BP to be added', metavar='vcf')
  options = parser.parse_args()
  
  main(options.output, options.input, options.length)