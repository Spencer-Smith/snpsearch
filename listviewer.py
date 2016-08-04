#Creates a list view for the output of ListMaker

class ListViewer:

	def __init__(self, output):
		self.out = open(output, 'w')

	def SetProteinNames(self, pNames):
		self.protNames = pNames

	def PrintData(self, data):
		#HEY THERE, DATA
		if data.function != "cds-indel":
			self.out.write(self.MakeAString(data, data.reference, "ref"))
		self.out.write(self.MakeAString(data, data.variant, "var"))

	def PrintPeptideData(self, data, peptide):
		alls = []
		if data.function != "cds-indel":
			ref = data.reference
			alls.append(ref)
		var = data.variant
		alls.append(var)

		string = ""
		for al in alls: 
			seq = peptide.sequence
			pos = data.position - peptide.pos_start
			seq = seq[:pos] + '(' + al.residue + ')' + seq[pos+1:]
			string += self.MakeAString(data, al, seq)

		self.out.write(string)	

	def MakeAString(self, data, al, s):
		#Builds the string
		string = data.snp_id + "\t" + data.gene + "\t" + data.accession + "\t" + str(data.position + 1)
		string += "\t" + data.function + "\t" + s
		if len(s) <= 3:
			string += "\t" + al.residue 
		string += self.MakeLine ((al.freq, al.EASfreq, al.EURfreq, al.AFRfreq, al.AMRfreq,
			al.SASfreq, al.variance))
		if len(self.protNames) > 0:
			string += "\t" + self.protNames[data.accession] + "\n"
		else:
			string += "\n"
		return string

	def MakeLine(self, values):
		line = ""
		i = 0
		while i < len(values):
			line += "\t" + str(round(values[i],4))
			i += 1
		return line

	def PrintHeader(self, opt):
		header = "SNP ID\tGene\tAccession\tAA Pos\tFunction" + opt + "\tGlobal\tEAS\tEUR\tAFR\tAMR\tSAS\tVariance"
		if len(self.protNames) > 0:
			header += "\tProtein Description\n"
		else:
			header += "\n"
		self.out.write(header)

	def CloseHandle(self):
		self.out.close()
