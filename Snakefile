SAMPLES = config['samples']

def generateOutputFiles():
	"""
	Generates output file paths for all listed input files.
	"""
	outputFiles = []
	for s in SAMPLES:
		if len(SAMPLES[s]['reads'].split(", ")) == 2: 	# Check for paired end data
			ending = ".PE.fasta"
		else:									# Else is single end data
			ending = ".SE.fasta"
		outputFiles.append("output/" + s.split("/")[-1] + ending)
	return(outputFiles)

def getPathOfSingleEndReads():
	"""
	Generates a list with all single end reads listed. Function is needed, so the rule SE_preprocessing knows what files to look for.
	"""
	filePaths = []
	for s in SAMPLES:
		if len(SAMPLES[s]['reads'].split(", ")) == 1:
			filePaths.append(SAMPLES[s])
	return(filePaths)

def getPathOfPairedEndReads():
	"""
	Generates a list with all paired end reads listed. Function is needed, so the rule PE_preprocessing knows what files to look for.
	"""
	filePaths = []
	for s in SAMPLES:
		paths = SAMPLES[s]['reads'].split(", ")
		if len(paths) == 2:
			filePaths.append(paths)
	return(filePaths)

rule all:
	input: generateOutputFiles()

include: "rules/single_end.smk"
include: "rules/paired_end.smk"
