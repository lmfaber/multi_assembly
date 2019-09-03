import sys
import os
import json
# Get infos from snakefile
fastpJson = snakemake.input["fastpJson"]
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

# convert to strings
kmers = [str(k) for k in kmers]



print(f'KMER SIZES: {kmers}')
reds = snakemake.input['reads']
print(f'PAIRED OR SINGLE END{len(reds)}')

if len(snakemake.input['reads']) == 2:
	print('PAIRED END DATA')
	###################
	# paired end path #
	###################

	# Get max read length from the json file of fastp
	with open(fastpJson, 'r') as jsonFile:
		data = json.load(jsonFile)
		maxReadLength = str(data['summary']['after_filtering']['read1_mean_length'])
		insertSize = str(data['insert_size']['peak'])
	
	
	R1 = snakemake.input["reads"][0]
	R2 = snakemake.input["reads"][1]
	
	# Create and write the config file
	outputConfig = baseFolder + "/config.txt"
	with open(outputConfig, 'w+') as configFile:
		configFile.write(f'# Maximal read length\n'
		f'max_rd_len={maxReadLength}\n'
		f'[LIB]\n'
		f'# Maximal read length in this lib\n'
		f'rd_len_cutof={maxReadLength}\n'
		f'# Average insert size\n'
		f'avg_ins={insertSize}\n'
		f'# If sequence needs to be reversed\n'
		f'# This option takes value 0 or 1. It tells the assembler if the read sequences need to be complementarily reversed. Illumima GA produces two types of paired-end libraries: a) forward-reverse, generated from fragmented DNA ends with typical insert size less than 500 bp; b) reverse-forward, generated from circularizing libraries with typical insert size greater than 2 Kb. The parameter "reverse_seq" should be set to indicate this: 0, forward-reverse; 1, reverse-forward.\n'
		f'reverse_seq={orientation}\n'
		f'# This indicator decides in which part(s) the reads are used. It takes value 1(only contig assembly), 2 (only scaffold assembly), 3(both contig and scaffold assembly).\n'
		f'asm_flags=3\n'
		f'# map_len: This takes effect in the "map" step and is the mininum alignment length between a read and a contig required for a reliable read location. The minimum length for paired-end reads and mate-pair reads is 32 and 35 respectively.\n'
		f'map_len=32\n'
		f'# Fastq file for paired end reads in two files\n'
		f'q1={R1}\n'
		f'q2={R2}')
else:
	print('SINGLE END DATA')
	###################
	# Single end path #
	###################
	
	# Get max read length from the json file of fastp
	with open(fastpJson, 'r') as jsonFile:
		data = json.load(jsonFile)
		maxReadLength = str(data['summary']['after_filtering']['read1_mean_length'])

	R1 = snakemake.input["reads"][0]
	# Create and write the config file
	outputConfig = baseFolder + "/config.txt"
	with open(outputConfig, 'w+') as configFile:
		configFile.write(f'# Maximal read length\n'
		f'max_rd_len={maxReadLength}\n'
		f'[LIB]\n'
		f'# Maximal read length in this lib\n'
		f'rd_len_cutof={maxReadLength}\n'
		f'# If sequence needs to be reversed\n'
		f'reverse_seq={orientation}\n'
		f'# In which part(s) the reads are used\n'
		f'asm_flags=3\n'
		f'# Minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)\n'
		f'map_len=32\n'
		f'# Fastq file for single reads\n'
		f'q={R1}')



for k in kmers:
	kmer_folder = f'{baseFolder}/{k}'
	if not os.path.exists(kmer_folder):
		os.mkdir(kmer_folder)
	os.system(f'SOAPdenovo-Trans-127mer all -s {baseFolder}/config.txt -o {kmer_folder}/{k} -p {threads} -K {k}')


# Merge assemblies using cd-hit-estimate and create final output
allKmerContigs = [f'{baseFolder}/{k}/{k}.contig' for k in kmers]
temporaryFile = f'{baseFolder}/temp.fasta'
cdHitEstOutput = f'{baseFolder}/soapdenovotrans.fasta'
contigsString = ' '.join(allKmerContigs)
os.system(f'cat {contigsString} > {temporaryFile}')
os.system(f'cd-hit-est -i {temporaryFile} -o {cdHitEstOutput} -d 30 -c 1 -T {threads}')

# Remove temporaryFile for a new instance
os.remove(temporaryFile)

## Change header because detonate has a problem with that

with open(cdHitEstOutput, "r") as soap, open(temporaryFile, 'w') as writer:
	for line in soap:
		if line.startswith(">"):
			line = line.replace(" ", "_")
			writer.write(line)
		else:
			writer.write(line)


# Remove the first unprocessed output from cd-hit-est
os.remove(cdHitEstOutput)

# Copy the processed temprary file into the correct output
os.rename(src = temporaryFile, dst = cdHitEstOutput)
