import sys
import numpy as np # Numpy
from scipy.stats.stats import pearsonr

if len(sys.argv) < 2:
	print("Must input dataset argument")
	sys.quit()

def readData(path):
        print("Reading " + path)
        # Read in Data
        GDS = open(path, "r")
        lines = GDS.readlines()
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

geneDict = readData(sys.argv[1])

allWntConnsJSON = open("miserables.json", "a+")
namesFile = open("canonicalNames.txt", "r")
connsFile = open("CanonicalConnectionsByName.txt", "r")

names = namesFile.readlines()
conns = connsFile.readlines()

# BUILD NODES LIST ============================================================================
allWntConnsJSON.write("{\n\t\"nodes\": [\n\t\t")

# Make list of WNT related genes------------------------------------

wntGenes = {}

for hsa in names:
	wntGenes[hsa.strip()] = 0

allConnsLst = list()
first = True
for pair in conns:
	tup = pair.strip().split(",")
	if tup[0] in geneDict and tup[1] in geneDict:
                corr, pVal = pearsonr(geneDict[tup[0]], geneDict[tup[1]])
		wntGenes[tup[0]] += 1
		wntGenes[tup[1]] += 1
	else:
		continue
	if first:
		first = False
		allConnsLst.append("{\"source\": \"" + tup[0].strip() + "\", \"target\": \"" + tup[1].strip() + "\", \"value\": " + str(abs(corr)) + "}")
	else:
		allConnsLst.append(",\n\t\t{\"source\": \"" + tup[0].strip() + "\", \"target\": \"" + tup[1].strip() + "\", \"value\": " + str(abs(corr)) + "}")


# Add all WNT related genes to Nodes List in JSON-------------------

first = True
for (gene, num) in wntGenes.items():
        if first:
                first = False
                allWntConnsJSON.write("{\"id\": \"" + gene + "\", \"group\": " + str(num) + "}")
        else:
                allWntConnsJSON.write(",\n\t\t{\"id\": \"" + gene + "\", \"group\": " + str(num) + "}")

# BUILD LINKS LIST ===========================================================================
allWntConnsJSON.write("\n\t],\n\t\"links\": [\n\t\t")

for line in allConnsLst:
	allWntConnsJSON.write(line)


allWntConnsJSON.write("\n\t]\n}")







