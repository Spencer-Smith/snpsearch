#This class is the go-between from listmaker to the database. It builds queries using data
# from listmaker and then runs those queries, returning the results back to listmaker.

import database
import SNPdata

class QueryBuilder:

	def __init__(self, path):
		self.db = database.Database()
		self.db.Connect(path)

	def QueryDatabase(self, query):
		cursor = self.db.QueryDatabase(query)
		result = cursor.fetchall()
		return result

	def SmartQuery(self, container, gene):
		#Figures out what kind of information has been given and queries based on that
		queries = []
		query = ""
		count = 0
		for item in container:
			if count % 500 == 0:
				query = query[:-3] + ")"
				if len(query) > 1:
					queries.append(query)
				query = self.BaseQuery() + " WHERE ("
			query += ("gene" if gene else "prot_acc") + " = '" + item + "' OR "
			count += 1

		query = query[:-3] + ")"
		if len(query) > 30:
				queries.append(query)
			
		return queries

	def BaseQuery(self):
		#Builds a string query which will be used to query the database. I had this before in the 
		# database class, but it seems out of place because it only interacts with the database in a 
		# super specific way. Probably should end up in some intermediate class between the database 
		# and the controller which handles database interactions for this program. Just want the 
		# database class to be a bit more reusable
		query = "SELECT * FROM fulltable"
		return query

	def PeptideQuery(self, accession_number, pos_min, pos_max):
		#Add a 'WHERE' clause for a peptide and its accession
		query = self.AccessionQuery(accession_number)
		query += " AND aa_pos BETWEEN " + str(pos_min) + " AND " + str(pos_max)
		return query

	def FindAccession(self, genes):
		#Build queries for 500 genes at a time
		count = 0
		query = ""
		queries = []
		for gene in genes:
			if count % 500 == 0:
				query = query[:-3] + ")"
				if len(query) > 1:
					queries.append(query)
				query = "SELECT prot_acc FROM fulltable WHERE ("
			query += "gene" + " = '" + gene + "' OR "
			count += 1
		
		#Get that last query in the container
		query = query[:-3] + ")"
		if len(query) > 36:
			queries.append(query)

		#Do the queries and add up the accessions
		accessions = {}
		for query in queries:
			result = self.QueryDatabase(query)
			for datum in result:
				accessions[datum[0]] = 1
			
		return accessions