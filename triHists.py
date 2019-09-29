import math
import os
import re
import pandas as pd # Pandas
import numpy as np # Numpy
import matplotlib.pyplot as plt # Matplotlibrary
import seaborn as sns # Seaborn Library
import random
from scipy.stats.stats import pearsonr
from scipy.stats.stats import ks_2samp

# Functions ------------------------------------------------------------------------------------
def intersection(lst1, lst2): 
    # Use of hybrid method 
    temp = set(lst2) 
    lst3 = [value for value in lst1 if value in temp] 
    return lst3 

def readData(path):
	print("Reading " + path)
        # Read in Data
        GDS1665 = open(path, "r")
        lines = GDS1665.readlines()
        colNames = lines[0]
        data = lines[2:]

        # Make Dictionary Gene Name --> Experiment Result
        geneDict = {}

        for line in data:
                lineArr = line.split("\t")
                try:
			geneDict[lineArr[1]] = np.array(lineArr[3:]).astype(np.float)
		except: 
			print("skip") 
	return geneDict
# END Functions ---------------------------------------------------------------------------------

# Parse WNT Data (Non-specific) -----------------------------------------------------------------

wntGeneData = open("wntGeneNames.txt", "r")
wntLines = wntGeneData.readlines()

# Make list of WNT related genes---------------------------------------------------------------

wntGenes = list()

for hsa in wntLines:
	genesInHsa = hsa.upper().split(":")[2].strip().split(", ")
        wntGenes.extend(genesInHsa)

# Parse WNT DATA (Specific)----------------------------------------------------------------------

wntPairData = open("allWntConnectionsByName.txt", "r")
wntPairLines = wntPairData.readlines()

# Pairwise WNT Genes

wntPairs = list()
for pair in wntPairLines:
	splitPair = pair.upper().strip().split(",")
	wntPairs.insert(0, (splitPair[0], splitPair[1]))

# END Parse--------------------------------------------------------------------------------------

path = "../../disk1/data/human/geneExpression/"
sigFile = open("datasetSignificance.txt","a+") 

for file in os.listdir(path)[897:]:
	if(file[-8:] != ".avg.pcl"):
		continue

	geneDict = readData(path + file)
	geneVals = geneDict.values()

	# Random Pearson Correlations-----------------------------------------------------------

	pearsons = list()
	
	works = True

	for i in range(0,100000):
		try:
			in1 = random.randint(0,len(geneVals)-1)
			in2 = random.randint(0,len(geneVals)-1)
			corr, pVal = pearsonr(geneVals[in1], geneVals[in2])
			if(corr == 1.0):
				corr = 0.99999
			pearsons.append(math.atanh(corr))
		except:
			
			print(file + ": " + str(corr) + ", " + str(geneVals[in1]) + ", " + str(geneVals[in2]))
			works = False
			break

	if(not works):
		continue

	# END Random Pearson Correlations--------------------------------------------------------
	
	# Non Specific Wnt Pearson Correlations--------------------------------------------------

	pearsonsWnt = list()
	
	for wnt1 in wntGenes:
		if wnt1 in geneDict:
			for wnt2 in wntGenes:
				if wnt2 in geneDict and wnt1 != wnt2:
					corr, pVal = pearsonr(geneDict[wnt1], geneDict[wnt2])
		                        if(corr == 1.0):
		                                corr = 0.99999
					pearsonsWnt.append(math.atanh(corr))

	# END Non Specific Wnt Pearson Correlations-----------------------------------------------
	
	# Specific Wnt Pearson Correlations-------------------------------------------------------
	
	pearsonsPairedWnt = list()
	for pair in wntPairs:
		if pair[0] in geneDict and pair[1] in geneDict:
			corr, pVal = pearsonr(geneDict[pair[0]], geneDict[pair[1]])
			if corr == 1.0:
				corr = 0.99999
			pearsonsPairedWnt.append(math.atanh(corr))

	# END Specific Wnt Pearson Correalations--------------------------------------------------

	# Determine if this is statistically significant------------------------------------------

	ks, pValNon = ks_2samp(pearsons, pearsonsWnt)
	ks, pValSpec = ks_2samp(pearsons, pearsonsPairedWnt)
	sigFile.write(file + " " + str(pValNon) + "," + str(pValSpec) + "\n")
	print(file + " " + str(pValNon) + "," + str(pValSpec) + "\n")

	# BUILD HISTOGRAM-------------------------------------------------------------------------

	if pValNon < (0.05 / 769.0) or pValSpec < (0.05 / 769.0):
		print("Number of Pairwise Wnt Data Points: " + str(len(pearsonsPairedWnt)))
		try:
			# the histogram of the random data
			plt.hist(pearsons, 75, density=True, range = [-2,2], facecolor='b', alpha=0.25)			# Blue Random Background Data
			plt.hist(pearsonsWnt, 75, density=True, range = [-2,2], facecolor='g', alpha=0.25)		# Green Non specific Wnt Data
			plt.hist(pearsonsPairedWnt, 75, density = True, range = [-2,2], facecolor='r', alpha=0.25)	# Red specific Wnt Data
		
			plt.xlabel('Correlation')
			plt.ylabel('Density of Correlations')
			plt.title(file)
			plt.grid(True)
			plt.savefig("triHists/" + file[:-8] + ".png")
			plt.clf()	
		except:
			print(file + ": Histogram not built")			

	# BUILD HISTOGRAM---------------------------------------------------------------------------
sigFile.close()
