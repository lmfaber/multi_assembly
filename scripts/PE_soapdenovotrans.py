import sys
import os
import json
# Get sample name from Snakefile
fastpJson = snakemake.input["fastpJson"]
R1 = snakemake.input["R1"]
R2 = snakemake.input["R2"]
sample = str(snakemake.params["configSample"])
baseFolder = snakemake.params['output_folder']
threads = str(snakemake.threads)
orientation = snakemake.params['orientation']

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
    insertSize = str(data['insert_size']['peak'])

# Create and write the config file
outputConfig = baseFolder + "/config.txt"
with open(outputConfig, 'w+') as configFile:
	configFile.write(f'# Maximal read length\n\
	max_rd_len={maxReadLength}\n\
	[LIB]\n\
	# Maximal read length in this lib\n\
	rd_len_cutof={maxReadLength}\n\
	# Average insert size\n\
	avg_ins={insertSize}\n\
	# If sequence needs to be reversed\n\
	# This option takes value 0 or 1. It tells the assembler if the read sequences need to be complementarily reversed. Illumima GA produces two types of paired-end libraries: a) forward-reverse, generated from fragmented DNA ends with typical insert size less than 500 bp; b) reverse-forward, generated from circularizing libraries with typical insert size greater than 2 Kb. The parameter "reverse_seq" should be set to indicate this: 0, forward-reverse; 1, reverse-forward.\n\
	reverse_seq={orientation}\n\
	# This indicator decides in which part(s) the reads are used. It takes value 1(only contig assembly), 2 (only scaffold assembly), 3(both contig and scaffold assembly).\n\
	asm_flags=3\n\
	# map_len: This takes effect in the "map" step and is the mininum alignment length between a read and a contig required for a reliable read location. The minimum length for paired-end reads and mate-pair reads is 32 and 35 respectively.\n\
	map_len=32\n\
	# Fastq file for paired end reads in two files\n\
	q1={R1}\n\
	q2={R2}')


# convert to strings
kmers = [str(k) for k in kmers]

for k in kmers:
    os.system(f'mkdir -p {baseFolder}/{k}')
    os.system(f'SOAPdenovo-Trans-127mer all -s {baseFolder}/config.txt -o {baseFolder}/{k}/{k} -p {threads} -K {k}')

# Merge assemblies using cd-hit-estimate and create final output
allKmerContigs = [f'{baseFolder}/{k}/{k}.contig' for k in kmers]
temporaryFile = f'{baseFolder}/temp.fasta'
cdHitEstOutput = f'{baseFolder}/soapdenovotrans.fasta'
contigsString = ' '.join(allKmerContigs)
os.system(f'cat {contigsString} > {temporaryFile}')
os.system(f'cd-hit-est -i {temporaryFile} -o {cdHitEstOutput}  -d 30 -c 1 -T {threads}')

# Remove temporaryFile for a new instance
os.system(f'rm {temporaryFile}')

## Change header because detonate has a problem with that

with open(cdHitEstOutput, "r") as soap, open(temporaryFile, 'w') as writer:
	for line in soap.readlines():
		if line.startswith(">"):
			line = line.replace(" ", "_")
			writer.write(line)
		else:
			writer.write(line)


# Remove the first unprocessed output from cd-hit-est
os.system("rm " + cdHitEstOutput)

# Copy the processed temprary file into the correct output
os.system("mv " + temporaryFile  + " " + cdHitEstOutput)

# Temporary file is not needed anymore so it can be deleted
#os.system("rm " + temporaryFile)
