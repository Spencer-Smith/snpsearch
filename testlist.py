#This script tests the output of the listmaker class

import sys

inf = open(sys.argv[1],'r')

inf.readline() #Get rid of header

val1 = 0
val2 = 0
oneoff = True
contradictions = 0
tolerance = 0.001

for line in inf:
	Bits = line.strip().split("\t")
	if Bits[4] == "cds-indel":
		continue
	if oneoff:
		val1 = float(Bits[7])
	else:
		val2 = float(Bits[7])
		tot = val1 + val2
		tot = round(tot,4)
		if abs(tot - 1) > tolerance:
			s = "\t"
			print(s.join((Bits[0],Bits[1],Bits[7],str(tot))))
			contradictions += 1
	oneoff = not oneoff

print("%s contradictions" %contradictions)
