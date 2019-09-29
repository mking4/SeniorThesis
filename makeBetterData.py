import os

path = "../../disk1/data/human/geneExpression"
for file in os.listdir(path)[142:]:
        if(file[:-4] + ".knn.pcl" in os.listdir(path) or file[-4:] != ".pcl" or file[:-4] + ".avg.pcl" in os.listdir(path) or file[-7:] == "avg.pcl"):
                continue
	os.system("cd " + path + "; java -Xmx1g -jar MeanGenesThatAgree.jar " + file + "> " + file[:-4] + ".avg.pcl") 
