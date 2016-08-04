#Proteomics Data
	#Input: data file (I guess)
	#Tasks: create a data structure from the data file(s)
	#Return: data structure

import peptide

class ProteomicsDataLoader():

	def __init__(self):
		#Creates a dictionary which will store the fasta sequences from file
		self.proteins = {}
		self.proteinNames = {}

	def LoadFastaFile(self, path):
		#This method reads in a fasta file and adds to the proteins dictionary
		FastaFile = open(path,'r')
		accession = ""
		sequence = ""

		for line in FastaFile:
			line = line.strip()
			if line[0] == '>':
				line = line.split('|')
				#Find the accession number
				i = 1
				for l in line:
					if "ref" in l:
						break
					i += 1
				accession = line[i] 
				#Get rid of the version number
				accession = accession[:-2]
				self.proteinNames[accession] = line[4]
			else:
				#After a header file, add the new lines of the sequence to that
				# accession number
				if not accession in self.proteins:
					self.proteins[accession] = line
				else:
					self.proteins[accession] += line

		FastaFile.close()
		return self.proteins

	def LoadPeptideFile(self, path, accessions):
		PeptideFile = open(path,'r')
		header = PeptideFile.readline()
		Peptides = {}

		for line in PeptideFile:
			p = self.LoadPeptide(line)
			if p.pos_start == -1:
				continue
			if p.accession in accessions:
				if p.accession in Peptides:
					Peptides[p.accession].append(p)
				else:
					Peptides[p.accession] = [p]

		PeptideFile.close()
		return Peptides

	def LoadPeptide(self, line):
		#Get the data
		line = line.strip().split()
		p = peptide.Peptide()
		p.sequence = self.ParseSequence(line[0])
		p.accession = self.ParseAccession(line[1])

		#If the accession number is in the full protein sequences we loaded from the  
		# fasta file, then find where the peptide lies within the protein
		if p.accession in self.proteins:
			protein = self.proteins[p.accession]
			p.pos_start = protein.find(p.sequence)
			p.pos_end = p.pos_start + len(p.sequence) - 1

		return p

	def ParseSequence(self, SeqString):
		#Clip ends of sequence
		sequence = SeqString[2:-2]
		#Remove unwanted (non-alphanumeric) characters
		temp = ""
		for ch in sequence:
			if ch.isalnum():
				temp += ch
		return temp

	def ParseAccession(self, IDstring):
		#Separate the accession number from other IDs
		ID = IDstring.split('|')
		accession = ID[1]
		#Remove version number
		return accession[:-2]