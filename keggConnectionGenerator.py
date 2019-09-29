import pandas as pd # Pandas
import numpy as np # Numpy

 # Read in Data
keggs = open("WntSignallingConnections.txt", "r")
lines = keggs.readlines()

geneDict = {}
for line in lines:
	lineArr = line.rstrip().split(": ")
	geneDict[lineArr[0].upper()] = lineArr[1].split(",")

relations = open("Canonical.txt","r")
rLines = relations.readlines()

allConnsName = open("CanonicalConnectionsByName.txt", "a+")
allConnsNumber = open("CanonicalConnectionsByNumber.txt", "a+")

for rl in rLines:
	rel = rl.rstrip().split(",")
	group1 = geneDict[rel[0].upper()]
	group2 = geneDict[rel[1].upper()]
	for g in group1:
		gArr = g.split("(") 	# TUPLE: Name, Number)
		for g2 in group2: 
			gArr2 = g2.split("(")
			allConnsNumber.write(gArr[0] + "," + gArr2[0] + "\n")
			allConnsName.write(gArr[1] + "," + gArr2[1] + "\n")
		
	
	










