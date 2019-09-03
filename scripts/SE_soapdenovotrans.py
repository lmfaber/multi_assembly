import sys
import os
import json

# Get sample name from Snakefile
fastpJson = snakemake.input["fastpJson"]
reads = snakemake.input["reads"]
sample = str(snakemake.params["configSample"])
baseFolder = snakemake.output[0]
threads = str(snakemake.threads)

os.system("mkdir -p " + baseFolder)

def calculate_kmers_for_soapdenovotrans(fastpFile):
	"""
	Generate automatically odd kmers for soapdenovotrans. Kmer length is roughly half and two thirds of the read length. Returns kmers as list.
	Arguments:
		fastpFile {str} -- fastp json file
	"""
	with open(fastpFile, 'r') as jsonFile:
		data = json.load(jsonFile)
		readLength = data['summary']['after_filtering']['read1_mean_length']

	halfLength = readLength//2
	if halfLength%2 == 0:
		halfLength += 1

	twoThird = readLength * 2//3
	if twoThird % 2 == 0:
		twoThird += 1
	return([halfLength, twoThird])

kmers = calculate_kmers_for_soapdenovotrans(fastpJson)

# Get max read length from the json file of fastp
with open(fastpJson, 'r') as jsonFile:
    data = json.load(jsonFile)
    maxReadLength = str(data['summary']['after_filtering']['read1_mean_length'])

# Create and write the config file
outputConfig = baseFolder + "/config.txt"
with open(outputConfig,)
configFile = open(outputConfig, "w+")
configFile.write(f'# Maximal read length\n\
max_rd_len={maxReadLength}\n\
[LIB]\n\
# Maximal read length in this lib\n\
rd_len_cutof={maxReadLength}\n\
# Average insert size\n\
avg_ins={insertSize}\n\
# If sequence needs to be reversed\n\
reverse_seq={orientation}\n\
# In which part(s) the reads are used\n\
asm_flags=3\n\
# Minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)\n\
map_len=32\n\
# Fastq file for single reads\n\
q={reads}')
configFile.close()

# convert to strings
kmers = [str(k) for k in kmers]

for k in kmers:
    os.system("mkdir -p " + baseFolder + "/" + k)
    #mkdir -p {params.folder}/{params.K}
    os.system("SOAPdenovo-Trans-127mer all -s " + baseFolder + "/config.txt -o " + baseFolder + "/" + k + "/" + k + " -p " + threads + " -K " + k)

# Merge assemblies using cd-hit-estimate and create final output
allKmerContigs = [baseFolder + "/" + k + "/" + k + ".contig" for k in kmers]
temporaryFile = baseFolder + "/temp.fasta"
cdHitEstOutput = baseFolder + "/soapdenovotrans.fasta"
os.system("cat " + " ".join(allKmerContigs) + " > " + temporaryFile)
os.system("cd-hit-est -i " + temporaryFile + " -o " + cdHitEstOutput + " -d 30 -c 1 -T " + threads)

# Remove temporaryFile for a new instance
os.system("rm " + temporaryFile)

## Change header because detonate has a problem with that
writer = open(temporaryFile, "w")

with open(cdHitEstOutput, "r") as soap:
	for line in soap.readlines():
		if line.startswith(">"):
			line = line.replace(" ", "_")
			writer.write(line)
		else:
			writer.write(line)
writer.close()
# Remove the first unprocessed output from cd-hit-est
os.system("rm " + cdHitEstOutput)

# Copy the processed temprary file into the correct output
os.system("mv " + temporaryFile  + " " + cdHitEstOutput)

# Temporary file is not needed anymore so it can be deleted
os.system("rm " + temporaryFile)
