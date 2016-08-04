import sqlite3
import SNPdata
import os

class Database:

	def __init__(self):
		#I know I don't have to initialize variables, but it helps me keep track of them
		self.connection = ""
		self.dbpath = ""
		self.connected = False

	def GetDatabasePath(self):
		return self.dbpath

	def IsConnected(self):
		return self.connected
		
	def Connect(self, path):
		#Establishes the connection to the database, checking that appropriate flags have
		# been set. Returns "False" if the path variable does not existor if a connection 
		# already exists
		if not os.path.exists(path):
			return False
		if self.connected:
			return False
		self.dbpath = path
		self.connection = sqlite3.connect(self.dbpath)
		self.connected = True
		return True
		
	def Disconnect(self):
		#Disconnects from the database, first checking that it's connected to a database. 
		# Returns "False" if not connected to a database
		if self.connected:
			self.connection.commit()
			self.connection.close()
			self.connected = False
			self.dbpath = ""
			return True
		else:
			return False

	def QueryDatabase(self, Query):
		#Queries the database and returns the result
		if self.connected:
			cursor = self.connection.execute(Query)
			return cursor

	def CreateTable(self, TableName, TableFields):
		command = "CREATE TABLE " + TableName + "\n("

		for key in TableFields:
			command += key + "\t" + TableFields[key] + "NOT NULL" + "\n"

		command = command[:-1]
		command += ");"

		self.connection.execute(command)
		self.connection.commit()

	def InsertRow(self, TableName, Row):
		command = "INSERT INTO " + TableName + " ("
		for key in Row:
			command += key + ","
		command[-1] = ")"

		command += " \ VALUES ("

		for key in Row:
			command += Row[key] + ","
		command[-1] = ")"

		self.connection.execute(command)