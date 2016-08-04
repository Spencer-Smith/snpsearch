#The purpose of this class is to load data from an ELISA file and store it
# in a container. It also has a function to convert the gene list to an
# accession list.

class eLoader:

	def Load(self, input_path):
		self.MakeGeneList(input_path)
		return self.genes		

	def MakeGeneList(self, path):
		#Create a list of all genes which appear in the file
		Elisa = open(path, 'r')
		self.genes = {}

		for line in Elisa:
			gene = line.strip()
			self.genes[gene] = 1

		Elisa.close()

	def GetAccessions(self, interactor):
		#Switch our list of genes for a lsit of accessions. This allows
		# us to compare the input ELISA with input	
		print("Finding accession numbers correspondent to gene names...")
		accessions = interactor.FindAccession(self.genes)
		return accessions
