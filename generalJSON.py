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
namesFile = open("WntSignallingConnections.txt", "r")
connsFile = open("CanonicalGenConnectionsByName.txt", "r")
groupRelsFile = open("Canonical.txt", "r")

names = namesFile.readlines()
conns = connsFile.readlines()
groupRels = groupRelsFile.readlines()

# BUILD RELS DICT ============================================================================

relsDict = {}
for line in groupRels:
	lineArr = line.upper().strip().split(",")
	relsDict[(lineArr[0], lineArr[1])] = list()

# BUILD NODES LIST ============================================================================
allWntConnsJSON.write("{\n\t\"nodes\": [\n\t\t")

# Make list of WNT related genes------------------------------------

wntGenes = {}

for hsa in names:
        geneGroup = hsa.upper().split(":")[0].strip()
	if not (geneGroup == "WNT5" or geneGroup == "WNT11"):
		wntGenes[geneGroup] = 0

first = True
for pair in conns:
	tup = pair.upper().strip().split(",")
	(fst, fstGroup) = tup[0].split(":")
	(snd, sndGroup) = tup[1].split(":")
	if fst in geneDict and snd in geneDict:
                cor, pVal = (pearsonr(geneDict[fst], geneDict[snd]))
		corr = abs(cor)
		if (fstGroup, sndGroup) in relsDict:
			relsDict[(fstGroup, sndGroup)].append(corr)
			wntGenes[fstGroup] += 1
                        wntGenes[sndGroup] += 1	
		else: 
			if (sndGroup, fstGroup) in relsDict:
				relsDisct[(sndGroup, fstGroup)].append(corr)
				wntGenes[fstGroup] += 1
	                        wntGenes[sndGroup] += 1
	else:
		continue


# Add all WNT related gene groups to Nodes List in JSON-------------------
for (gene, num) in wntGenes.items():
        if first:
                first = False
                allWntConnsJSON.write("{\"id\": \"" + gene + "\", \"group\": " + str(num) + "}")
        else:
                allWntConnsJSON.write(",\n\t\t{\"id\": \"" + gene + "\", \"group\": " + str(num) + "}")

# BUILD LINKS LIST ===========================================================================
allWntConnsJSON.write("\n\t],\n\t\"links\": [\n\t\t")

first = True
for ((fst, snd), corrs) in relsDict.items():
	avgCorr = 0.0
	if len(corrs) > 0:
		avgCorr = sum(corrs) / len(corrs)
	if first and avgCorr != 0:
		first = False
		allWntConnsJSON.write("{\"source\": \"" + fst + "\", \"target\": \"" + snd + "\", \"value\": " + str(abs(avgCorr)) + "}")
	else:
		if avgCorr != 0:
			allWntConnsJSON.write(",\n\t\t{\"source\": \"" + fst + "\", \"target\": \"" + snd + "\", \"value\": " + str(abs(avgCorr)) + "}")

allWntConnsJSON.write("\n\t]\n}")







