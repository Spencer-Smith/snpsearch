#Prerequisites
*Uses Python 3.5.1

#OVERVIEW
The purpose of listmaker is to find SNPs (small nucleotide polymorphisms) within
the context given to the scouter by input files. Listmaker connects to a SQLite 
database made up of tables from dbSNP which contain data for non-synonomous SNPs
where at least one population has a minor allele frequency of at least 1%. This 
data is almost entirely derived from the 1000 Genmome Project and tracks the
frequency of alleles within 5 populations: East Asia (EAS), Europe (EUR), Africa
(AFR), Americas (AMR), and South Asia (SAS). Unlike SNPscouter, listmaker then
takes the output from database queries and sorts it by highest minor allele 
frequency. It also calculates and prints variance within the populations.

#EXECUTION
Again, the information that is queried from the database and returned is based
on the context of the information given. That context is decided by the user at 
runtime based on the input files given to listmaker. Thus, the SNPs returned
can be filtered to more specific groups. These are in the input/filter options: 

1. Filter to only SNPs in genes associated with ELISAs. Input a .txt file with a
 list of genes:
	> listmaker.py -e elisa_filepath

1. Filter to only SNPs in proteins about which we have protein data. Input the
 a fasta file:
	> listmaker.py -f fasta_filepath

2. Filter SNPs which are in both a list of ELISA and a fasta file:
	> listmaker.py -e elisa_filepath -f fasta_filepath

3. The full scope of listmaker is to return SNPs which lie within in peptides 
 from genes for which an ELISA exists. A fasta file is also used to place each
 peptides within the context of its protein:
	> listmaker.py -e elisa_filepath -f fasta_filepath -p peptide_filepath 

Optionally, each command can be run with the option "-d", which will use the
path following the option as the database (by default, will use the database
located in "dbfiles\SNP.db"); and option "-o" which will specify the output file
(by default, output goes to "unnamed_list.txt").


#OUTPUT
Once the input has been prepared, queried, and sorted, it is output. Only SNP   
data having a variance among populations of at least 5% is written to file.