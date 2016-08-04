#This class has the test cases containing unit tests for listmaker.

import unittest
import listmaker
import SNPdata
import elisaloader
import querybuilder
import peptide
import proteomicsdataloader

class ListyHelpersTestCase(unittest.TestCase):

	def setUp(self):
		self.querybuilder = querybuilder.QueryBuilder("dbfiles\SNP.db")
		self.pdl = proteomicsdataloader.ProteomicsDataLoader()

	#SNPdata
	def test_assign_function(self):
		"""Can SNPdata correctly assign the function?"""
		expectedSNP = SNPdata.SNPdata([[0 for i in range(12)],[0 for j in range (12)]])
		expectedSNP.function = "missense"
		testData = [[1,"GENE","NP_01",4,8,"R",0.9,0.9,0.9,0.9,0.9,0.9],
			[1,"GENE","NP_01",4,42,"N",0.1,0.1,0.1,0.1,0.1,0.1]]
		testSNP = SNPdata.SNPdata(testData)

		self.assertEqual(expectedSNP.function, testSNP.function)

	def test_process_alleles(self):
		"""Can SNPdata correctly determine which allele is the reference and variant,
			 based on the alleles' frequencies?"""
		testData = [[1,"GENE","NP_01",4,8,"R",0.9,0.9,0.9,0.9,0.9,0.9],
			[1,"GENE","NP_01",4,42,"N",0.1,0.1,0.1,0.1,0.1,0.1]]
		testSNP = SNPdata.SNPdata(testData)

		self.assertTrue(testSNP.reference.freq == 0.9)
		self.assertTrue(testSNP.variant.freq == 0.1)

	#eLoader
	def test_gene_list(self):
		"""Does the MakeGeneList method of eLoader properly load gene names?"""
		expectedgenes = {}
		expectedgenes["APOE"] = 1
		expectedgenes["SELE"] = 1
		expectedgenes["SPAG1"] = 1
		eload = elisaloader.eLoader()
		testgenes = eload.Load("TEST_files\TEST_elisa.txt")
		self.assertDictEqual(expectedgenes, testgenes)

	#QueryBuilder
	def test_find_accessions_single(self):
		"""Can QueryBuilder find the protein accession of a gene?"""
		expectedlist = {"NP_000441":1}
		testlist = self.querybuilder.FindAccession({"SELE":1})
		self.assertDictEqual(expectedlist, testlist)

	def test_find_accessions_multi(self):
		"""Can QueryBuilder find multiple protein accessions of a gene?"""
		expectedlist = {"NP_000032":1, "NP_001289617":1, "NP_001289618":1,
			"NP_001289619":1, "NP_001289620":1}
		testlist = self.querybuilder.FindAccession({"APOE":1})
		self.assertDictEqual(expectedlist, testlist)

	def test_smart_query_gene(self):
		"""Can SmartQuery correctly build a gene query?"""
		testquery = self.querybuilder.SmartQuery({"APOE":1}, True)
		self.assertTrue("SELECT * FROM fulltable WHERE (gene = 'APOE' )" in testquery)

	def test_smart_query_short_gene(self):
		"""Does SmartQuery build a gene query given a small gene?"""
		testquery = self.querybuilder.SmartQuery({"TN":1}, True)
		self.assertTrue("SELECT * FROM fulltable WHERE (gene = 'TN' )" in testquery)

	def test_smart_query_XP_accession(self):
		"""Does SmartQuery build an accession query given an XP accession?"""
		testquery = self.querybuilder.SmartQuery({"XP_01":1}, False)
		self.assertTrue("SELECT * FROM fulltable WHERE (prot_acc = 'XP_01' )" in testquery)

	def test_smart_query_NP_accession(self):
		"""Does SmartQuery build an accession query given an NP accession?"""
		testquery = self.querybuilder.SmartQuery({"NP_01":1}, False)
		self.assertTrue("SELECT * FROM fulltable WHERE (prot_acc = 'NP_01' )" in testquery)

	def test_smart_query_peptide(self):
		"""Can SmartQuery properly build a peptide query?"""
		testpeptide = peptide.Peptide()
		testpeptide.sequence = "TMNT"
		testpeptide.accession = "NP_01"
		testpeptide.pos_start = 2
		testpeptide.pos_end = 5
		testquery = self.querybuilder.SmartQuery({"NP_01":1}, False)
		expectedsubstring = "SELECT * FROM fulltable WHERE (prot_acc = 'NP_01' )"
		self.assertTrue(expectedsubstring in testquery)

	#ProteomicsDataLoader
	def test_parse_sequence(self):
		"""Can ParseSequence() parse a sequence?"""
		testseq = "I.LOVE.U"
		expectedparse = "LOVE"
		self.assertEqual(expectedparse, self.pdl.ParseSequence(testseq))

	def test_parse_accession(self):
		"""Can ParseAccession() locate and parse out the accession number?"""
		testline = "words|NP_01.1|morewords|evenmorewords"
		expectedparse = "NP_01"
		self.assertEqual(expectedparse, self.pdl.ParseAccession(testline))

	def test_load_fasta(self):
		"""Does LoadFastaFile() store protein data correctly?"""
		expectedprot = {}
		sequence = "MIASQFLSALTLVLLIKESGAWSYNTSTEAMTYDEASAYCQQRYTHLVAIQNKEEIEYLN"
		sequence += "SILSYSPSYYWIGIRKVNNVWVWVGTQKPLTEEAKNWAPGEPNNRQKDEDCVEIYIKREK"
		sequence += "DVGMWNDERCSKKKLALCYTAACTNTSCSGHGECVETINNYTCKCDPGFSGLKCEQIVNC"
		sequence += "TALESPEHGSLVCSHPLGNFSYNSSCSISCDRGYLPSSMETMQCMSSGEWSAPIPACNVV"
		sequence += "ECDAVTNPANGFVECFQNPGSFPWNTTCTFDCEEGFELMGAQSLQCTSSGNWDNEKPTCK"
		sequence += "AVTCRAVRQPQNGSVRCSHSPAGEFTFKSSCNFTCEEGFMLQGPAQVECTTQGQWTQQIP"
		sequence += "VCEAFQCTALSNPERGYMNCLPSASGSFRYGSSCEFSCEQGFVLKGSKRLQCGPTGEWDN"
		sequence += "EKPTCEAVRCDAVHQPPKGLVRCAHSPIGEFTYKSSCAFSCEEGFELHGSTQLECTSQGQ"
		sequence += "WTEEVPSCQVVKCSSLAVPGKINMSCSGEPVFGTVCKFACPEGWTLNGSAARTCGATGHW"
		sequence += "SGLLPTCEAPTESNIPLVAGLSAAGLSLLTLAPFLLWLRKCLRKAKKFVPASSCQSLESD"
		sequence += "GSYQKPSYIL"
		expectedprot["NP_000441"] = sequence
		sequence = "MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQT"
		sequence += "LSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVAEETRARLSKELQAAQA"
		sequence += "RLGADMEDVCGRLVQYRGEVQAMLGQSTEELRVRLASHLRKLRKRLLRDADDLQKRLAVY"
		sequence += "QAGAREGAERGLSAIRERLGPLVEQGRVRAATVGSLAGQPLQERAQAWGERLRARMEEMG"
		sequence += "SRTRDRLDEVKEQVAEVRAKLEEQAQQIRLQAEAFQARLKSWFEPLVEDMQRQWAGLVEK"
		sequence += "VQAAVGTSAAPVPSDNH"
		expectedprot["NP_000032"] = sequence

		path = "TEST_files\TEST.fasta"
		testprot = self.pdl.LoadFastaFile(path)
		self.assertDictEqual(expectedprot,testprot)

	def test_load_peptide(self):
		"""Can LoadPeptide() store peptide data from file to data structure?"""
		line = "R.PEPTIDE.C\tref|NP_01.1|gi|111|"
		expectedpep = peptide.Peptide()
		expectedpep.sequence = "PEPTIDE"
		expectedpep.accession = "NP_01"
		expectedpep.pos_start = 1
		expectedpep.pos_end = 7
		self.pdl.proteins["NP_01"] = "APEPTIDESEQUENCE"
		testpep = self.pdl.LoadPeptide(line)
		self.assertPeptideEqual(expectedpep, testpep)		

	def assertPeptideEqual(self, pep1, pep2):
		self.assertEqual(pep1.sequence, pep2.sequence, msg = "Sequences not equal")
		self.assertEqual(pep1.accession, pep2.accession, msg = "Accession numbers not equal")
		self.assertEqual(pep1.pos_start, pep2.pos_start, 
			msg = "Start positions not equal")
		self.assertEqual(pep1.pos_end, pep2.pos_end, msg = "End positions not equal")

	def test_load_peptide_file(self):
		"""Does LoadPeptideFile() store peptide data correctly?"""
		self.pdl.LoadFastaFile("TEST_files\TEST.fasta")

		pep = peptide.Peptide()
		pep.sequence = "MEDVCGRLVQYRGE"
		pep.accession = "NP_000032"
		pep.pos_start = 125
		pep.pos_end = 138
		expectedpep = {}
		expectedpep["NP_000032"] = [pep]

		path = "TEST_files\TEST_peptides.txt"
		testpep = self.pdl.LoadPeptideFile(path, self.pdl.proteins)

		for item in expectedpep:
			if item in testpep:
				count = 0
				for val in expectedpep[item]:
					self.assertPeptideEqual(val, testpep[item][count])
					count += 1
			else:
				self.assertTrue(item in testpep, msg= "Sequence in expected dictionary not found")

class ListMakerTestCase(unittest.TestCase):

	def setUp(self):
		self.listy = listmaker.Listy()
		self.listy.output = "TEST_FILES\TEST_list.txt"
		self.listy.dbpath = "dbfiles\SNP.db"

	def test_listy_elisa(self):
		"""Does Listy work as expected when given only an ELISA file?"""
		self.listy.elipath = "TEST_FILES\TEST_elisa.txt"
		self.listy.Main()

		expectedoutput = "TEST_FILES\expected_listy_elisa.txt"
		self.assertFilesEqual(self.listy.output, expectedoutput)

	def test_listy_fasta(self):
		"""Does Listy work as expected when given only a fasta file?"""
		self.listy.faspath = "TEST_FILES\TEST.fasta"
		self.listy.Main()

		expectedoutput = "TEST_FILES\expected_listy_fasta.txt"
		self.assertFilesEqual(self.listy.output, expectedoutput)

	def test_listy_elisa_and_fasta(self):
		"""Does Listy work as expected when given ELISA and fasta files?"""
		self.listy.elipath = "TEST_FILES\TEST_elisa.txt"
		self.listy.faspath = "TEST_FILES\TEST.fasta"
		self.listy.Main()

		expectedoutput = "TEST_FILES\expected_listy_elifas.txt"
		self.assertFilesEqual(self.listy.output, expectedoutput)

	def test_listy_peptides(self):
		"""Does Listy work as expected when given ELISA, fasta, and peptide files?"""
		self.listy.elipath = "TEST_FILES\TEST_elisa.txt"
		self.listy.faspath = "TEST_FILES\TEST.fasta"
		self.listy.peppath = "TEST_files\TEST_peptides.txt"
		self.listy.Main()

		expectedoutput = "TEST_FILES\expected_listy_peptides.txt"
		self.assertFilesEqual(self.listy.output, expectedoutput)

	def assertFilesEqual(self, file1, file2):
		in1 = open(file1,'r')
		in2 = open(file2,'r')

		lines1 = {}
		lines2 = {}

		for line in in1:
			lines1[line] = 1

		for line in in2:
			lines2[line] = 1

		in1.close()
		in2.close()

		self.assertDictEqual(lines1, lines2)

if __name__ == "__main__":
	unittest.main()