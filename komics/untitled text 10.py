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
			query='tmp.LC1412C3.oriented.fasta',
			subject='tmp.LC1412C3.oriented.fasta',
			out='tmp.LC1412C3_2.xml',
			dust="no",
			max_target_seqs=3,
			outfmt=5)

    
    hits=[]
		blast_parse_result=defaultdict(list)
		for record in NCBIXML.parse(open('tmp.LC1412C3_2.xml')):
			for alignment in record.alignments:
				if str(record.query) != str(alignment.hit_id):																			# removes self-hits
					for hsp in alignment.hsps:
						perc_identity=100*float(hsp.identities)/float(hsp.align_length)
						if perc_identity >= 95:																					# retains alignments with minimum percent identity
							if hsp.align_length >= record.query_length:
							  if record.query not in blast_parse_result.keys() and record.query not in str(blast_parse_result.values()):
							    blast_parse_result[str(record.query)].append(str(alignment.hit_id))
							  elif record.query in blast_parse_result.keys() and alignment.hit_id not in blast_parse_result[str(record.query)]:
							    blast_parse_result[str(record.query)].append(str(alignment.hit_id))

							  oneway=['H', hsp.align_length, round(perc_identity, 1), str(record.query), str(alignment.hit_id)]
							  twoway=['H', hsp.align_length, round(perc_identity, 1), str(alignment.hit_id), str(record.query)]
							  if oneway and twoway not in hits:
							    hits.append(oneway)
							  print oneway
							  print twoway
						  
							  alignment.hit_id not in blast_parse_result[str(record.query)]:
and alignment.hit_id 
('LC1412C3_99_contig862_len740_circularized', 'LC1412C3_119_contig223_len740_circularized')
('LC1412C3_99_contig862_len740_circularized', 'LC1412C3_109_contig1228_len740_circularized')

('LC1412C3_109_contig1228_len740_circularized', 'LC1412C3_119_contig223_len740_circularized')

('LC1412C3_109_contig1228_len740_circularized', 'LC1412C3_99_contig862_len740_circularized')

('LC1412C3_119_contig223_len740_circularized', 'LC1412C3_109_contig1228_len740_circularized')
('LC1412C3_119_contig223_len740_circularized', 'LC1412C3_99_contig862_len740_circularized')

							    
							    blast_parse_result[str(record.query)].append(str(alignment.hit_id))
							  elif record.query in blast_parse_result.keys() and alignment.hit_id not in blast_parse_result[str(record.query)]:
							    blast_parse_result[str(record.query)].append(str(alignment.hit_id))
							  
							  
							    # only unique hits, not two-way hits
							    # 
							    str('H', hsp.align_length, round(perc_identity, 1), str(record.query), str(alignment.hit_id))
							    							  
							    							  record.query_length
	blast_parse_result[str(record.query)].append(str(alignment.hit_id))						  
							  if record.query not in blast_parse_result:
							    blast_parse_result[record.query].append(alignment.hit_id)
							    
							  print(str(record.query), str(alignment.hit_id))
							  if alignment.hit_id not in blast_parse_result:
							    lengths[str(record.query)]=int(record.query_length)
							    lengths[str(alignment.hit_id)]=int(alignment.length)
							    minlength = min(lengths.keys(), key=(lambda k: lengths[k]))
							    blast_parse_result[minlength]=lengths[minlength]
#								else:																																				# removes two-way matches
#									if blast_parse_result[alignment.hit_id] > record.query_length:						# if new matched sequence is shorter than a previously existing match
#										del blast_parse_result[alignment.hit_id]
#										blast_parse_result[record.query]=record.query_length

'LC1412C3_119_contig218_len748_circularized': ['LC1412C3_109_contig346_len748_circularized', 'LC1412C3_99_contig1210_len748_circularized']
u'LC1412C3_119_contig218_len748_circularized': [u'LC1412C3_109_contig346_len748_circularized', u'LC1412C3_99_contig1210_len748_circularized']

('LC1412C3_99_contig862_len740_circularized', 'LC1412C3_119_contig223_len740_circularized')
('LC1412C3_99_contig862_len740_circularized', 'LC1412C3_109_contig1228_len740_circularized')
('LC1412C3_109_contig1228_len740_circularized', 'LC1412C3_119_contig223_len740_circularized')

('LC1412C3_109_contig1228_len740_circularized', 'LC1412C3_99_contig862_len740_circularized')
('LC1412C3_119_contig223_len740_circularized', 'LC1412C3_109_contig1228_len740_circularized')
('LC1412C3_119_contig223_len740_circularized', 'LC1412C3_99_contig862_len740_circularized')


# a table output would be nice, to know which contigs are similar for the different kmer values?