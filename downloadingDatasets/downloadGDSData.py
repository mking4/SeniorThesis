import os

path = "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/"
data = open("GDSSets.txt", "r")
gdsDirectory = "GDS"
for line in data:
	l = line.strip()
	if len(l) == 7:
		gdsDirectory = line[0:4]		
	cdComm = "cd ~/../../disk1/data/human/geneExpression; "
	curlComm = "curl " + path + gdsDirectory + "nnn/" + l + "/soft/" + l + ".soft.gz > " + l + ".soft.gz; "
	#print(curlComm)
	gunZipComm = "gunzip " + l + ".soft.gz; "
	convertComm = "ruby convertSoft2Pcl.rb " + l + ".soft " + l + ".pcl"
	os.system(cdComm + curlComm + gunZipComm + convertComm)

