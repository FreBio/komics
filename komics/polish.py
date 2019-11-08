import re
import os
import sys
import subprocess

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
            self.csb1 = re.compile(r'GGGCGTTC|GGGCGTGC')
            self.csb3 = re.compile(r'GGGGTTGGTGT|GGGGTTGATGT')
            self.csb3_comp = re.compile(r'CCCCAACCACA|CCCCAACTACA')
            self.csb3_rev = re.compile(r'TGTGGTTGGGG|TGTAGTTGGGG')
            self.csb3_rev_comp = re.compile(r'ACACCAACCCC|ACATCAACCCC')
            self.output_vsearch = os.path.abspath('tmp.' + self.out + '.uc')
            self.min_identity = int(minidentity)
            if not os.path.exists(self.input_contigs):
                raise Error('Contigs file not found:' + self.input_contigs)


    def complement(self, sequence):
        complement = {'A':'T','C':'G','G':'C','T':'A','-':'-'}
        return "".join([complement.get(nt.upper(), '') for nt in sequence])


    def reorientate(self, sequence, id):
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
        else:
            sys.stderr.write("WARNING: no CSB3-mer found for " + str(id) + ". Skipping this one. \n")


    def put_csb_start(self, sequence):
        matches=tuple(re.finditer(self.csb1, str(sequence)))
        if len(matches) == 1:
            try:
                match=re.search(self.csb1, str(sequence))
                chop=match.span()[0]
                return sequence[chop:]+sequence[:chop]
            except AttributeError:
                match = re.search(self.csb1, str(sequence[50:]+sequence[:50]))
                chop=match.span()[0]
                return sequence[chop:]+sequence[:chop]
        if len(matches) > 1:
            fst_csb1 = matches[0].span()[0]
            snd_csb1 = matches[1].span()[0]
            try:
                match_csb3=re.search(self.csb3, str(sequence))
                chop_csb3=match_csb3.span()[0]
            except AttributeError:
                match_csb3 = re.search(self.csb3, str(sequence[50:]+sequence[:50]))
                chop_csb3=match_csb3.span()[0]
            if chop_csb3 > snd_csb1:
                return sequence[snd_csb1:]+sequence[:snd_csb1]
            elif chop_csb3 < fst_csb1:
                return sequence[snd_csb1:]+sequence[:snd_csb1]
            else:
                return sequence[fst_csb1:]+sequence[:fst_csb1]


    def orient_run(self):
        newseqs = []
        for minicircle in SeqIO.parse(self.input_contigs, "fasta"):
            if re.search('circular', minicircle.id):
                tmp=Seq(str(self.put_csb_start(self.reorientate(minicircle.seq, minicircle.id))))
                newseqs.append(SeqRecord(tmp, id=minicircle.id, description=""))
            else:
                tmp=Seq(str(self.reorientate(minicircle.seq, minicircle.id)))
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
        with open(self.output_vsearch, 'rb') as f:
            for line in f:
                line=bytes.decode(line)
                if line.startswith('C'):
                    contigC1=line.split('\t')[8].rstrip()
                    numberC1=line.split('\t')[1].rstrip()
                    clusters[numberC1].append(contigC1)

        keep = []
        for c in clusters:
            keep.append(str(clusters[c][0]))
        
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
