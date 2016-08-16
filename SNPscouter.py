#This class uses the other classes to read the input, prepare a container,
# query the database based on the contents of the container, and then sort
# and print the results of the queries.

import sys
import getopt
import os
import elisaloader
import proteomicsdataloader
import querybuilder
import listviewer
import SNPdata
from datetime import datetime

class SNPscouter:

	def __init__(self):
		self.dbpath = "dbfiles\SNP.db"
		self.peppath = ""
		self.faspath = ""
		self.elipath = ""
		self.output = "unnamed_list.txt"
		self.data = {}

	def Main(self):
		#Prepare database for queries
		self.database = querybuilder.QueryBuilder(self.dbpath)

		#Fill a container
		self.PrepareContainer()

		#Query the database with the contents of container
		self.MakeQueries()

		#Sort and print the results of the queries
		self.SortAndPrint()

	def PrepareContainer(self):
		#CONTAINER will be filled with whatever data we're given
		self.container = {}
		self.eloader = elisaloader.eLoader()
		self.pdataloader = proteomicsdataloader.ProteomicsDataLoader()

		if self.elipath != "":
			#Load in the list of genes from the ELISA file
			genes = self.eloader.Load(self.elipath)
			print("%s genes loaded from ELISA file."%len(genes))
			#If all we have is an ELISA file, the genes list is
				# our container
			self.container = genes
		
		#Given a fasta file, the process is a bit more complex
		if self.faspath != "":
			accessions = {}
			#Load in accession numbers from the fasta file
			proteomicsAccessions = self.pdataloader.LoadFastaFile(self.faspath)
			print("%s accessions loaded from fasta file."%len(proteomicsAccessions))

			#If we also have an ELISA file, we'll need to do a bit more work
			if self.elipath != "":
				#Get a list of accession numbers from the genes in the ELISA file
				ELISAaccessions = self.eloader.GetAccessions(self.database)
				print("%s accessions found from ELISA genes."%len(ELISAaccessions))
				for acc in ELISAaccessions:
					if acc in proteomicsAccessions:
						accessions[acc] = 1
				print("%s accessions overlapping between ELISA and fasta accessions."%len(accessions))
			#If not, we'll only use the accessions from the fasta file
			else:
				for acc in proteomicsAccessions:
					accessions[acc] = 1

			#Either way, accessions is our new container of interest
			self.container = accessions
		
		#If we have a peptide file, we'll need to load it in
		if self.peppath != "":
			#Loads the peptides from the peptide file if one is present,  
			# including only those which are in the accessions dictionary
			peptides = self.pdataloader.LoadPeptideFile(self.peppath, accessions)
			print("%s peptides found in peptide file."%len(peptides))

			#If we've gotten this far, peptides becomes our container of interest
			self.container = peptides

	def MakeQueries(self):
		#Query for whatever information (genes, accessions, peptides) has been input
		gene = True
		ptide = False

		for item in self.container:
			if self.container[item] != 1:
				ptide = True
				gene = False
				break
			#Accessions always have the form of '(N/X)P_####'
			if len(item) > 2 and item[2] == '_':
				gene = False
				break

		#Generates queries
		queries = self.database.SmartQuery(self.container, gene)

		print("Querying the database and storing results...")		
		count = 0
		sTime = datetime.now()
		for query in queries:
			results = self.database.QueryDatabase(query)
			self.StoreResults(results, ptide)
			count += 500
			soFar = (datetime.now() - sTime).total_seconds()
			if count > len(self.container):
				count = len(self.container) 
			self.Progress(count, len(self.container), soFar)
		print()

	def Progress(self, count, total, soFar):
		progress = count/total
		progress = round(progress, 4)
		tRemaining = (soFar / progress) - soFar
		mRemaining = str(int(tRemaining / 60))
		sRemaining = str(int(tRemaining % 60))
		tString = ((mRemaining + "m " + sRemaining + "s") if int(mRemaining) > 0 
				else (sRemaining + "s"))
		pString = str(round(progress * 100, 2))
		string = "\tProgress: " + pString + "%   Est. Time Remaining: " + tString + "    "
		print(string, end = "\r")
		
	def StoreResults(self, results, ptide):
		#Divide the data into SNPs
		rows = []
		SNP = ""
		prot_acc = ""
		for row in results:
			if not(row[0] == SNP and row[2] == prot_acc):
				#Load the SNPs into SNPdata structures and store it
				self.AddToData(rows, SNP, prot_acc, ptide)
				rows = []
				SNP = row[0]
				prot_acc = row[2]
			if row[0] == SNP and row[2] == prot_acc:
				rows.append(row)
		#Store the last data in the loop
		self.AddToData(rows, SNP, prot_acc, ptide)

	def AddToData(self, rows, SNP, prot_acc, ptide):
		#Check if it's an indel
		indel = False
		for r in rows:
			if r[4] == 45:
				indel = True
		#If there's only one line that's not enough, unless it's an indel
		if len(rows) >= 2 or indel:
			d = SNPdata.SNPdata(rows)
			index = str(SNP) + "." + prot_acc
			temp = SNPdata.SNPdata(rows)
			#Bail on any data with a variance below 5%
			if temp.variant.variance < 0.05:
				return
			if ptide:
				index = prot_acc
				peptides = self.container[index]
				newpep = []
				for pep in peptides:
					if temp.position <= pep.pos_start or temp.position >= pep.pos_end:
						newpep.append(pep)
				if len(newpep) < 1:
					return
				self.container[index] = newpep
			self.data[index] = temp

	def SortAndPrint(self):
		#This method goes through the data container, printing and then deleting
		# the entry with the highest MAF until all have been printed
		print("Sorting and printing...")
		viewer = listviewer.ListViewer(self.output)
		viewer.SetProteinNames(self.pdataloader.proteinNames)

		if self.peppath == "":
			viewer.PrintHeader("\tAllele\tResidue")
		else:
			viewer.PrintHeader("\tPeptide")
		count = 0
		#Sort the data
		while len(self.data) > 0:
			highkey = self.FindHighest()	
			#If no entry is found, end the loop
			if highkey == '':
				break
			#Print the found dictionary entry
			if self.peppath != "":
				for pep in self.container[highkey]:
					viewer.PrintPeptideData(self.data[highkey], pep)
			else:
				viewer.PrintData(self.data[highkey])
			#Remove the entry from the dictionary
			del self.data[highkey]
			count += 1

		print("%s results output"%count)
		viewer.CloseHandle()

	def FindHighest(self):
		highkey = ""
		highval = 0.0
			#Find the highest MAF still existing in the dictionary
		for key in self.data:
			varfreq = self.data[key].variant.freq
			if varfreq > highval:
				highkey = key
				highval = varfreq
			#On a tie, look for the higher variance between population MAF
			if varfreq == highval:
				current = self.data[highkey].variant
				new = self.data[key].variant
				if new.variance > current.variance:
					highkey = key
		return highkey

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
		if "-p" in OptionsSeen and not "-f" in OptionsSeen:
			print("  ERROR: Fasta file (-f) required where peptide file (-p) used")
			error = True

		#Exit if there's an error anywhere
		if error:
			sys.exit(1)			

if __name__ == "__main__":
	starttime = datetime.now()
	migo = SNPscouter()
	migo.ParseCommandLine(sys.argv[1:])
	migo.Main()
	print("Elapsed time: " + str(datetime.now() - starttime))