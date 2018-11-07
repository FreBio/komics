'''
TO DO:
> herorientatie moet kunnen op basis van eender welke sequentie. Default is CSB3, maar als user iets anders geeft gaat dat ook, bijvoorbeeld guide RNA"s or repeats, eender welke sequentie!
> eventueel ook met optie om specifiek slechts 1 sequentie te herorienteren, bijvoorbeeld, in het geval dat er meerdere CSB3s zijn
> WARNING: assumes that _circuralized is in the name, if not, it will not circularize the molecule!
> pymummer>=0.9.0 could perhaps help too identify similar minicircles, instead of blast
''' 

import re
import os
import sys
import string
import subprocess

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline


class Error (Exception): pass


class Polisher:

	def __init__(self,
		out,
		fasta,
		minidentity,
		):
			self.out = out
			self.input_contigs = os.path.abspath(fasta)
			self.oriented_contigs = os.path.abspath('tmp.' + self.out + '.oriented.fasta')
			self.minicircles = os.path.abspath(self.out + '.minicircles.fasta')
			self.csb3 = re.compile(r'GGGGTTGGTGT|GGGGTTGATGT')
			self.csb3_comp =	re.compile(r'CCCCAACCACA|CCCCAACTACA')
			self.csb3_rev =	re.compile(r'TGTGGTTGGGG|TGTAGTTGGGG')
			self.csb3_rev_comp = re.compile(r'ACACCAACCCC|ACATCAACCCC')
#			self.output_xml = os.path.abspath('tmp.' + self.out + '.2.xml')
			self.output_vsearch = os.path.abspath('tmp.' + self.out + '.uc')
			self.min_identity = int(minidentity)
			if not os.path.exists(self.input_contigs):
				raise Error('Contigs file not found:' + self.input_contigs)


	def complement(self, sequence):
		complement = {'A':'T','C':'G','G':'C','T':'A', '-':'-'}
		return "".join([complement.get(nt.upper(), '') for nt in sequence])


	def reorientate(self, sequence):
		if re.findall(self.csb3, str(sequence)):
			return sequence
		elif re.findall(self.csb3_rev, str(sequence)):
			return sequence[::-1]
		elif re.findall(self.csb3_comp, str(sequence)):
			return self.complement(sequence[:-1])
		elif re.findall(self.csb3_rev_comp, str(sequence)):
			return self.complement(sequence[::-1])
		elif re.findall(self.csb3, str(sequence[100:]+sequence[:100])):
			newseq=sequence[100:]+sequence[:100]
			return newseq
		elif re.findall(self.csb3_rev, str(sequence[100:]+sequence[:100])):
			newseq=sequence[100:]+sequence[:100]
			return newseq[::-1]
		elif re.findall(self.csb3_comp, str(sequence[100:]+sequence[:100])):
			newseq=sequence[100:]+sequence[:100]
			return self.complement(newseq[:-1])
		elif re.findall(self.csb3_rev_comp, str(sequence[100:]+sequence[:100])):
			newseq=sequence[100:]+sequence[:100]
			return self.complement(newseq[::-1])
		
	
	def put_csb3_start(self, sequence):
		match=re.search(self.csb3, str(sequence))
		chop=match.span()[0]
		return sequence[chop:]+sequence[:chop]


	def orient_run(self):
		newseqs = []
		for minicircle in SeqIO.parse(self.input_contigs, "fasta"):
			if re.search('circular', minicircle.id):
				tmp=Seq(str(self.put_csb3_start(self.reorientate(minicircle.seq))))
				newseqs.append(SeqRecord(tmp, id=minicircle.id, description=""))
			else:
				tmp=Seq(str(self.reorientate(minicircle.seq)))
				newseqs.append(SeqRecord(tmp, id=minicircle.id, description=""))

		SeqIO.write(newseqs, self.oriented_contigs, "fasta")
	
	
	def vsearch_run(self):
		vsearch_command = [
			"vsearch",
			"--cluster_fast", self.oriented_contigs,
			"--id", str(float(self.min_identity)/float(100)), 
			"--uc", self.output_vsearch
			]
		sys.stderr.write('Running VSEARCH.\n')
		subprocess.call(' '.join(vsearch_command), shell = True)


	def process_uc(self):
		sys.stderr.write('\nProcessing VSEARCH output.\n')
		clusters = defaultdict(list)
		contigC1 = defaultdict(str)
		contigH1 = defaultdict(str)
		contigH2 = defaultdict(str)
		with open(self.output_vsearch, 'rb') as f:
			for line in f:
				if line.startswith('C'):
					contigC1=line.split('\t')[8].rstrip()
  				if contigC1 not in clusters[line.split('\t')[1]]:
  					clusters[line.split('\t')[1]].append(contigC1)
  			if line.startswith('H'):
  				contigH1=line.split('\t')[8].rstrip()
  				contigH2=line.split('\t')[9].rstrip()
  				if contigH1 not in clusters[line.split('\t')[1]]:
  					clusters[line.split('\t')[1]].append(contigH1)
  				if contigH2 not in clusters[line.split('\t')[1]]:
  					clusters[line.split('\t')[1]].append(contigH2)
		
		keep = []
		for c in clusters:
			keep.append(str(clusters[c][1]))
		
		polisher_result=[]
		for minicircle in SeqIO.parse(self.oriented_contigs, "fasta"):
			if minicircle.id in keep:
				polisher_result.append(SeqRecord(minicircle.seq, id=minicircle.id, description=""))

		SeqIO.write(polisher_result, self.minicircles, "fasta")
		
		sys.stderr.write('Found %i unique (>%i) minicircles contigs.\n' % (len(keep), self.min_identity))
		sys.stderr.write('Polished minicircles can be found in %s.\n' % (self.minicircles))
		
	
	def run(self):
		sys.stderr.write('\nPolishing\n')
		sys.stderr.write('=========\n')
		self.orient_run()
		self.vsearch_run()
		self.process_uc()
		sys.stderr.write('\nkomics polish successfully completed.\n\n')

'''
	def blast_run(self):
		"""
		All-against-all blast of the minicircle contigs
		Notes:
			* Uses megablast
			* Turns dust filter off by default
			* Outputs in XML format by default. This behaviour cannot be changed to guarantee correct parsing results.
			* Retains only top-hits to save memory by default
		"""
		
		blastcmd = NcbiblastnCommandline(
		  task='megablast',
			query=self.oriented_contigs,
			subject=self.oriented_contigs,
			out=self.output_xml,
			dust="no",
			max_target_seqs=3,
			outfmt=5)
		
		sys.stderr.write('Running BLAST with the following command:\n')
		sys.stderr.write(str(blastcmd)+ '\n')
		stdout, stderr=blastcmd()


	def blast_parse(self):
		"""
		Parses blast output to identify overlapping minicircles
		"""

		sys.stderr.write('Parsing BLAST output.\n')
		
		blast_parse_result={}
		for record in NCBIXML.parse(open(self.output_xml)):
#		for record in NCBIXML.parse(open('multi.2.xml')):
			for alignment in record.alignments:
				if str(record.query) != str(alignment.hit_id):																			# removes self-hits
					for hsp in alignment.hsps:
						lengths=dict()
						perc_identity=100*hsp.identities/hsp.align_length
						if perc_identity >= self.min_identity:																					# retains alignments with minimum percent identity
							if hsp.align_length >= record.query_length:																		# retains only alignments with length equal or larger than the query length
#								print(str(record.query), str(alignment.hit_id))
								if alignment.hit_id not in blast_parse_result:
									lengths[str(record.query)]=int(record.query_length)
									lengths[str(alignment.hit_id)]=int(alignment.length)
									minlength = min(lengths.keys(), key=(lambda k: lengths[k]))
									blast_parse_result[minlength]=lengths[minlength]
#								else:																																				# removes two-way matches
#									if blast_parse_result[alignment.hit_id] > record.query_length:						# if new matched sequence is shorter than a previously existing match
#										del blast_parse_result[alignment.hit_id]
#										blast_parse_result[record.query]=record.query_length

# a table output would be nice, to know which contigs are similar for the different kmer values?



'''