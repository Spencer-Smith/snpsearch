#Data structure to hold information from the database

class SNPdata():

	class Allele():

		def __init__(self):
			self.residue = ""
			self.freq = 0
			self.EASfreq = 0
			self.EURfreq = 0
			self.AFRfreq = 0
			self.AMRfreq = 0
			self.SASfreq = 0
			self.variance = 0

	def __init__(self, data):
		#The data argument is a list of the information for each entry in the database
		# for a single protein at a single position (or within a peptide?). The information 
		# is contained in a list with the following information at the specified indices:
			# [0] = snpcontiglocusid.snp_id
			# [1] = snpcontiglocusid.gene
			# [2] = snpcontiglocusid.prot_acc
			# [3] = snpcontiglocusid.aa_pos
			# [4] = snpcontiglocusid.function
			# [5] = snpcontiglocusid.residue
			# [6] = snpallelefreq.freq
			# [7] = allelefreqbysspop.eas
			# [8] = allelefreqbysspop.eur
			# [9] = allelefreqbysspop.afr
			#[10] = allelefreqbysspop.amr
			#[11] = allelefreqbysspop.sas
		self.snp_id = str(data[0][0])
		self.gene = data[0][1]
		self.accession = data[0][2]
		self.position = data[0][3]
		self.function = ""
		self.AssignFunction(data)
		if self.function != "cds-indel":
			self.reference = self.Allele()
		self.variant = self.Allele()

		self.ProcessAlleleData(data)

	def AssignFunction(self, data):
		#Assign the function variable
		functionCodes = {41:"stop gain", 42:"missense", 43:"stop loss",
			44:"frameshift", 45:"cds-indel"}
		
		for allele in data:
			if allele[4] in functionCodes:
				self.function = functionCodes[allele[4]]

	def ProcessAlleleData(self, data):
		#For information about this data, see the method above
		for row in data:
			ally = self.Allele()
			ally.residue = row[5]
			ally.freq = row[6]
			ally.EASfreq = row[7]
			ally.EURfreq = row[8]
			ally.AFRfreq = row[9]
			ally.AMRfreq = row[10]
			ally.SASfreq = row[11]
			ally.variance = self.FindVariance(ally)
			if row[6] > .5:
				self.reference = ally
			else:
				self.variant = ally

	def FindVariance(self, allele):
		high = max(allele.EASfreq,allele.EURfreq,allele.AFRfreq,allele.AMRfreq,allele.SASfreq)
		low = min(allele.EASfreq,allele.EURfreq,allele.AFRfreq,allele.AMRfreq,allele.SASfreq)

		return high - low