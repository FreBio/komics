#!/usr/bin/env python

import re
import os
import sys
import subprocess
import argparse

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(
		description = 'Realignes minicircles to reduces spurious alignments with insertions and deletions. Commonly used for T. cruzi minicircles that contain multiple conserved regions',
        usage = 'realign.py <out> <fasta>')
parser.add_argument('out', help='Prefix used for labeling FASTA output file and sequence headers', metavar='out')
parser.add_argument('fasta', help='FASTA file with minicircle contigs', metavar = 'fasta')
parser.add_argument('--csb1', type=str, help='Specificy one or more Conserved Sequence Block 1 (CSB1) sequences, which will put at the start of the alignments. When nothing is provided, the following two CSB1-mers are used by default: "GGGCGTTC|GGGCGTGC".', metavar='STR')
options = parser.parse_args()
args = vars(parser.parse_args())

if args["csb1"] is None:
	csb1 = 'GGGCGTTC|GGGCGTGC'
elif args["csb1"] is not None:
	csb1 = args["csb1"]


sys.stderr.write('Running VSEARCH.\n')
vsearch_command = [
	"vsearch",
	"--cluster_fast", os.path.abspath(options.fasta),
	"--id", str(float(99)/float(100)), 
	"--uc", os.path.abspath('tmp.'+ options.out + '.fst.clusters.uc')
	]

subprocess.call(' '.join(vsearch_command), shell = True)

sys.stderr.write('\nProcessing VSEARCH output.\n')
insertions = defaultdict(list)
deletions = defaultdict(list)
with open(os.path.abspath('tmp.'+ options.out + '.fst.clusters.uc'), 'rb') as f:
	for line in f:
		line=bytes.decode(line)
		if line.startswith('H'):
			cigar = line.split('\t')[7].rstrip()
			if cigar != "=":
				insertion = re.search(r'\d\d\dI', cigar)
				deletion = re.search(r'\d\d\dD', cigar)
				if insertion and insertion.start() == 0:
					ins_length = insertion.group().replace('I','')
					del_length = deletion.group().replace('D','')
					if (int(del_length) == int(ins_length)):
						insertions[line.split('\t')[8].rstrip()].append(str(insertion.group().replace('I','')))
				if deletion and deletion.start() == 0:
					ins_length = insertion.group().replace('I','')
					del_length = deletion.group().replace('D','')
					if (int(ins_length) == int(del_length)):
						deletions[line.split('\t')[8].rstrip()].append(str(deletion.group().replace('D','')))

sys.stderr.write('Realigning minicircles.\n')
newseqs = []

for minicircle in SeqIO.parse(os.path.abspath(options.fasta), "fasta"):
	if re.search('circular', minicircle.id):
		if minicircle.id in insertions:
			start = len(str(minicircle.seq)) - int(*insertions[minicircle.id])
			seq=Seq(str(minicircle.seq[start:] + minicircle.seq[:start]))
			newseqs.append(SeqRecord(seq, id=minicircle.id, description=""))
		elif minicircle.id in deletions:
			start = int(*deletions[minicircle.id])
			seq=Seq(str(minicircle.seq[start:] + minicircle.seq[:start]))
			newseqs.append(SeqRecord(seq, id=minicircle.id, description=""))
		else:
			newseqs.append(SeqRecord(minicircle.seq, id=minicircle.id, description=""))

sys.stderr.write('Done.\n')

SeqIO.write(newseqs, os.path.abspath(options.out + '.realigned.minicircles.fasta'), "fasta")

