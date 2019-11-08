import os
import sys

from Bio	import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline


class Error (Exception): pass


class Circularizer:

	def __init__(self, 
		out,
		fasta, 
		minoverlap=20,
		):
			self.sample_name = out
			self.input_contigs = os.path.abspath(fasta)
			self.output_xml = os.path.abspath('tmp.' + self.sample_name + '.1.xml')
			self.output_fasta = os.path.abspath('tmp.' + self.sample_name + '.circularized.fasta')
			self.min_overlap = minoverlap
			self.min_identity = int(100)
			if not os.path.exists(self.input_contigs):
				raise Error('Contigs file not found:' + self.input_contigs)


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
			query=self.input_contigs,
			subject=self.input_contigs,
			out=self.output_xml,
			dust="no",
			max_target_seqs=3,
			outfmt=5)
		
		sys.stderr.write('Running BLAST with the following command:\n')
		sys.stderr.write(str(blastcmd)+ '\n\n')
		stdout, stderr=blastcmd()


	def blast_parse(self):
		"""
		Parses blast output to identify overlapping ends in each minicircle contig
        Currently, nothing is further reported on that
		"""

		sys.stderr.write('Parsing BLAST output:\n')
		
		blast_parse_result={}
		for record in NCBIXML.parse(open(self.output_xml)):
			for alignment in record.alignments:
				if str(record.query) == str(alignment.hit_def):													# retains self-hits only
					for hsp in alignment.hsps:
						if not (hsp.sbjct_start == hsp.query_start and hsp.sbjct_end == hsp.query_end):			# removes self-matches
							if hsp.align_length >= self.min_overlap:												# retains alignments with minimum length
								perc_identity=100*hsp.identities/hsp.align_length
								if perc_identity >= self.min_identity:											# retains alignments with minimum percent identity
									if hsp.query_start == 1 and hsp.sbjct_end == record.query_length:
										#sys.stderr.write('Contig %s is likely circular: %i-%i matches %i-%i with %i perc. identity\n' 
										#									% (record.query.split('_')[1], hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end, perc_identity))
										blast_parse_result[record.query]=hsp.identities
										
		sys.stderr.write('Found %i circular minicircles.\n' % (len(blast_parse_result)))
		return blast_parse_result


	def circularize(self):
		"""
		Chops one of the overlapping ends as identified with blast
		"""

		circularity=self.blast_parse()
		circularize_result=[]
		
		for minicircle in SeqIO.parse(self.input_contigs, "fasta"):
			if minicircle.id in circularity:
				newseq=minicircle.seq[circularity[minicircle.id]:]
				newid=self.sample_name + '_contig' + minicircle.id.split('_')[1] + '_len' + str(len(newseq)) + '_circularized'
				circularize_result.append(SeqRecord(newseq, id=newid, description=""))
			else:
				newid=self.sample_name + '_contig' + minicircle.id.split('_')[1] + '_len' + str(len(minicircle.seq))
				circularize_result.append(SeqRecord(minicircle.seq, id=newid, description=""))
		
		SeqIO.write(circularize_result, self.output_fasta, "fasta")
		sys.stderr.write('Circularized minicircles can be found in %s.\n ' % (self.output_fasta))


	def run(self):
		sys.stderr.write('\nCircularizing minicircle contigs\n')
		sys.stderr.write('================================\n')
		self.blast_run()
		self.circularize()
		sys.stderr.write('\nkomics circularize successfully completed.\n\n')
