#A class to load the parsed, tab-delimited dbSNP files into a database

import database
import os
import sys

class DatabaseLoader:

	def __init__(self):
		self.Directory = ""
		self.DBpath = "SNP.db"
		self.SetUpDatabase()

	def SetUpDatbase(self):
		self.database = database.Database()
		self.database.SetDatabasePath(self.Directory + "\\" + self.DBpath)
		self.database.Connect()

	def Main(self):
		#Here we explicitly name each of our tables and pass that name into
		# the MakeTable method
		self.MakeTable("Allele")
		self.MakeTable("AlleleFreqBySsPop")
		self.MakeTable("SNPAlleleFreq")
		self.MakeTable("SNPContigLocusId")
		self.MakeTable("SNPSubSNPLink")

		self.database.Disconnect()

	def MakeTable(self, TableName):
		path = self.Directory + "\\" + self.TableName + ".Parsed.bcp"
		Handle = open(path, 'r')
		header = Handle.readline()

		#Break up the header to get the name for each field in the table
		fields = {}
		dice = header.strip().split()
		for d in dice:
			fields[d] = "VARCHAR(255)"

		self.database.CreateTable(TableName, fields)

		#Each line is parsed for its info to be inserted into the table
		for line in Handle:
			linebits = line.strip().split()
			counter = 0
			#The fields dictionary, previously used to declare the fields in the
			# CreateTable method, will now hold the info for each line
			for key in fields:
				fields[key] = linebits[0]
				counter += 1
			self.database.InsertRow(TableName, fields)

		Handle.close()

	def ParseCommandLine(self, Arguments):
		(Options, Args) = getopt.getopt(Arguments, "d:")
		OptionsSeen = {}
		for (Option, Value) in Options:
			OptionsSeen[Option] = 1
			if Option == "-d":
				# -d directory
				if not os.path.exists(Value):
					print("** Error: couldn't find directory '%s'\n\n"%Value)
					print(UsageInfo)
					sys.exit(1)
				self.Directory = Value
		# Error out, if we didn't see required options:
		if not OptionsSeen.has_key("-d"):
			print("** Please specify input dir (-d)")
			print(UsageInfo)
			sys.exit(1)

if __name__ == "__main__":
	amigo = DatabaseLoader()
	amigo.ParseCommandLine(sys.argv[1:])
	amigo.Main()
