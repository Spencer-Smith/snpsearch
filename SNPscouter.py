#This is the main script which manages the other python class files. It takes as arguments the
# paths to the database, peptide file, fasta file, and a list of genes for elisa tests. Using
# this information, it outputs any relevant SNP information

import sys
import getopt
import os
import proteomicsdataloader
import databaseinteraction
import SNPviewer
from datetime import datetime

class SNPscouter:

	def __init__(self):
		self.dbpath = ""
		self.peppath = ""
		self.faspath = ""
		self.elipath = ""
		self.output = "RENAME.txt"

	def Main(self):
		#Prepare database for queries
		dbi = databaseinteraction.DatabaseInteraction(self.dbpath)
		sv = SNPviewer.SNPviewer(self.output)
		globalMAF = .1
		variance = .1

		if self.elipath == "":
			query = dbi.BaseQuery()
			data = dbi.QueryDatabase(query)
			sv.PrintHeader()
			sv.WriteAllData(data, globalMAF, variance)
		else:
			#Figure out what accession number is for each gene in the ELISA file 
			self.LoadElisa(self.elipath)
			accessions = dbi.GetAccession(self.genes)

			if self.faspath == "":
				sv.PrintHeader()
				for acc in accessions:
					query = dbi.BaseQuery() + dbi.AccessionIs(acc)
					data = dbi.QueryDatabase(query)
					sv.WriteGeneData(data, globalMAF, variance)

		if self.faspath != "":
			#Load and prepare peptide and fasta files
			pdl = proteomicsdataloader.ProteomicsDataLoader()
			proteins = pdl.LoadFastaFile(self.faspath)
			overlap = {}
			for pro in accessions:
				if pro in proteins:
					overlap[pro] = 1
			
			if self.peppath == "":
				sv.PrintHeader()
				for pro in overlap:
					query = dbi.BaseQuery() + dbi.AccessionIs(pro)
					data = dbi.QueryDatabase(query)
					sv.WriteGeneData(data, globalMAF, variance)

		if self.peppath != "":
			#Loads the peptides from the peptide file, excluding the  
			# ones which are not in the accessions dictionary (also 
			# ignores repeats)
			peptides = pdl.LoadPeptideFile(self.peppath, overlap)
			#After everything is loaded, we run this loop over the queries
			for peptide in peptides:
				pep = peptides[peptide]
				query = dbi.BaseQuery() + dbi.AccessionIs(pep.accession_number) 
				query += dbi.PeptideLimits(pep.pos_start, pep.pos_end)
				#Gets query data for a single peptide
				data = dbi.QueryDatabase(query)
				if len(data) == 0:
					continue
				data = dbi.LoadSNPdata(data)
				sv.WritePeptideData(pep, data, globalMAF, variance)
		sv.CloseHandles()

	def LoadElisa(self, path):
		#Method creates a list of all genes which appear in the file
		Elisa = open(path, 'r')
		self.genes = {}

		for line in Elisa:
			gene = line.strip()
			self.genes[gene] = 1

		Elisa.close()

	def ParseCommandLine(self, Arguments):
		(Options, Args) = getopt.getopt(Arguments, "d:p:f:e:o:")
		OptionsSeen = {}
		error = False
		#Set our internal variables to the passed in values
		#Error if the file paths aren't found
		for (Option, Value) in Options:
			OptionsSeen[Option] = 1
			if Option == "-d":
				# -d database file
				if not os.path.exists(Value):
					print("  ERROR: couldn't find database '%s'\n\n"%Value)
					error = True
				self.dbpath = Value
			if Option == "-p":
				# -p peptide file
				if not os.path.exists(Value):
					print("  ERROR: couldn't find peptide '%s'\n\n"%Value)
					error = True
				self.peppath = Value	
			if Option == "-f":
				# -f fasta file
				if not os.path.exists(Value):
					print("  ERROR: couldn't find fasta '%s'\n\n"%Value)
					error = True
				self.faspath = Value
			if Option == "-e":
				# -f elisa file
				if not os.path.exists(Value):
					print("  ERROR: couldn't find elisa '%s'\n\n"%Value)
					error = True
				self.elipath = Value
			if Option == "-o":
				# -o output file
				self.output = Value
		# Error out, if we didn't see required options:
		if not "-d" in OptionsSeen:
			print("  Please specify database path (-d)")
			error = True
		if "-f" in OptionsSeen and not "-e" in OptionsSeen:
			print("  ELISA file (-e) required where fasta file (-f) used")
			error = True			
		if "-p" in OptionsSeen and not "-f" in OptionsSeen:
			print("  Fasta file (-f) required where peptide file (-p) used")
			error = True

		#Exit if there's an error anywhere
		if error:
			sys.exit(1)			

if __name__ == "__main__":
	startTime = datetime.now()
	pogo = SNPscouter()
	pogo.ParseCommandLine(sys.argv[1:])
	pogo.Main()
	print("Elapsed time: " + str(datetime.now() - startTime))