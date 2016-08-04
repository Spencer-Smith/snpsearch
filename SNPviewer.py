#Class to show information returned from queries to SNP databse
import SNPdata

class SNPviewer:

	def __init__(self, outpath):
		#Get rid of file ending if there is one
		pos = outpath.find('.')
		if pos == -1:
			pos = len(outpath)
		outpath = outpath[:pos]

		self.ghandle = open(outpath + ".global.txt", 'w')
		self.vhandle = open(outpath + ".variance.txt", 'w')
		self.jhandle = open(outpath + ".jon.txt",'w')

	def WritePeptideData(self, peptide, data, globalMAF, variance):
		#Basic (but organized) outputting of data from the queries.
		
		#Make the reference allele more accessible
		ref = data.reference

		#Set some print flags
		gprint = ref.freq < (1 - globalMAF)
		mx = max(ref.EASfreq,ref.EURfreq,ref.AFRfreq,ref.AMRfreq,ref.SASfreq)
		mn = min(ref.EASfreq,ref.EURfreq,ref.AFRfreq,ref.AMRfreq,ref.SASfreq)
		vprint = mx - mn > variance
		jprint = False

		#We want to make a string which indicates with parenthesis the amino acid that the SNP affects
		seq = peptide.peptide
		pos = data.position - peptide.pos_start
		seq = seq[:pos] + '(' + seq[pos] + ')' + seq[pos+1:]

		#Write the "header" data for the peptide
		string = "Peptide:  " + str(peptide.pos_start) + " " + seq + " " + str(peptide.pos_end) + "\n" + "  Gene: "
		string += data.gene + "  Accession: " + data.accession + " at position: " + str(data.position) + "\n"
		
		#Put all the data for variants in a list)
		alleleline = "  Reference: " + ref.residue + " (" + ref.allele + ")"
		totals = [ref.freq]
		eas = [ref.EASfreq]
		eur = [ref.EURfreq]
		afr = [ref.AFRfreq]
		amr = [ref.AMRfreq]
		sas = [ref.SASfreq] 
		for var in data.variants:
			alleleline += "  Variant: " + var.residue + " (" + var.allele + ")" 
			totals.append(var.freq)
			eas.append(var.EASfreq)
			eur.append(var.EURfreq)
			afr.append(var.AFRfreq)
			amr.append(var.AMRfreq)
			sas.append(var.SASfreq)
			vm = max(var.EASfreq,var.EURfreq,var.AFRfreq,var.AMRfreq,var.SASfreq)
			if vm > .1:
				jprint = True

		#Ouput all the data from the lists.
		string += alleleline + "\n"
		string += "\tTotal: " + self.MakeLine(totals) + "\n"
		string += "\t  EAS: " + self.MakeLine(eas) + "\n"
		string += "\t  EUR: " + self.MakeLine(eur) + "\n"
		string += "\t  AFR: " + self.MakeLine(afr) + "\n"
		string += "\t  AMR: " + self.MakeLine(amr) + "\n"
		string += "\t  SAS: " + self.MakeLine(sas) + "\n\n"

		if gprint:
			self.ghandle.write(string)
		if vprint:
			self.vhandle.write(string)
		if jprint:
			self.jhandle.write(string)

	def WriteGeneData(self,data, globalMAF, variance):
		if len(data) == 0:
			return
		#Data is as follows:
			#[0] = snpcontiglocusid.gene
			#[1] = snpcontiglocusid.prot_acc
			#[2] = snpcontiglocusid.aa_pos
			#[3] = snpcontiglocusid.function
			#[4] = snpcontiglocusid.residue
			#[5] = snpallelefreq.freq
			#[6] = allelefreqbysspop.eas
			#[7] = allelefreqbysspop.eur
			#[8] = allelefreqbysspop.afr
			#[9] = allelefreqbysspop.amr
			#[10] = allelefreqbysspop.sas
		string = ""
		cur_pos = 0
		strings = []
		gprint = vprint = jprint = False

		for row in data:
			string = ""
			#Print the buffer and reset flags when we've moved on to a new SNP
			if row[2] != cur_pos:
				cur_pos = row[2]
				self.PrintBuffer(gprint, vprint, jprint, strings)
				#Reset print flags
				gprint = vprint = jprint = False
				#Reset string buffer
				strings = []

			#Validate this gene for each output criteria
			if row[6] > globalMAF and row[6] < .5:
				gprint = True
			if max(row[7:]) - min(row[7:]) > variance:
				vprint = True
			if max(row[7:]) > .1 and max(row[7:]) < .5:
				jprint = True

			#Build string:
				#Grab information we need: gene, position, residue
			string += row[0] + "\t" + row[1] + "\t" + str(row[2]) + "\t"
			if row[5] > .5:
				string += "ref\t"
			else:
				string += "var\t"
			string += row[4] + "\t"
			
			#Round and add in the frequencies
				#[5] is overall freq
				#[6] is EAS freq
				#[7] is EUR freq
				#[8] is AFR freq
				#[9] is AMR freq
				#[10] is SAS freq
			i = 5
			while i < 11:
				string += str(round(row[i],4)) + "\t"
				i += 1
			string = string[:-1] + "\n"

			#Collect all the strings for the gene
			strings.append(string)

	def WriteAllData(self,data, globalMAF, variance):
		curgene = ""
		i = 0
		genedata = []
		for line in data:
			if line[0] == curgene:
				genedata.append(line)
			else:
				if len(genedata) > 0:
					self.WriteGeneData(genedata,globalMAF,variance)
					genedata = []
				curgene = line[0]
				genedata.append(line)

	def PrintBuffer(self,gprint,vprint,jprint,strings):
		for s in strings:
			if gprint:
				self.ghandle.write(s)
			if vprint:
				self.vhandle.write(s)
			if jprint:
				self.jhandle.write(s)

	def PrintHeader(self):
		header = "Gene\tAccession\tAA Pos\tAllele\tResidue\tGlobal\tEAS\tEUR\tAFR\tAMR\tSAS\n"
		self.ghandle.write(header)
		self.vhandle.write(header)
		self.jhandle.write(header)

	def MakeLine(self, values):
		line = ""
		for v in values:
			line += "  " + str(round(v,4))
		return line

	def CloseHandles(self):
		self.ghandle.close()
		self.vhandle.close()
		self.jhandle.close()