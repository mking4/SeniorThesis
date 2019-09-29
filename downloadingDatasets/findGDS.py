import json

neededData = open("GDSSets.txt", "a")

with open("humanDataSets.json", "r") as humanData:
	parsed = json.load(humanData)

for dataset in parsed:
	if dataset["etype"]["slug"] == "co-expression": 
		neededData.write(dataset["identifier"]+"\n")

