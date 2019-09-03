############################
##### single end reads #####
############################

def orientation_rnaspades(wildcards):
	orientation = config['samples'][wildcards.sample]['orientation']
	if orientation == 'F':
		return '--pe1-fr'
	elif orientation == 'R':
		return '--pe1-rf'
	else:
		print('No valid orientation.')

def orientation_soapdenovotrans(wildcards):
	orientation = config['samples'][wildcards.sample]['orientation']
	if orientation == 'F':
		return '0'
	elif orientation == 'R':
		return '1'
	else:
		print('No valid orientation.')

rule SE_header_correction:
	input: getPathOfPairedEndReads()
	output: 
	 R1 = temp('output/{sample}/SE_{sample}_1.fastq')
	script: '../scripts/parser.py'

rule SE_fastp:
	input: rules.SE_header_correction.output
	output:
		R1 = 'output/{sample}/fastp/{sample}.fastq',
		json = 'output/{sample}/fastp/{sample}_SE.json',
		html = 'output/{sample}/fastp/{sample}_SE.html'
	conda: '../envs/fastp.yml'
	shell: 'fastp -i {input} -o {output.R1} --json {output.json} --html {output.html}'

rule SE_rnaspades:
	input:
		R1 = rules.SE_fastp.output.R1,
		fastpJson = rules.SE_fastp.output.json
	output: 'output/{sample}/SE_assemblies/rnaspades/transcripts.fasta'
	params: 
		orientation = orientation_rnaspades,
		output_folder = 'output/{sample}/SE_assemblies/rnaspades'
	threads: config['threads']
	conda: '../envs/rnaspades.yml'
	benchmark: 'output/{sample}/benchmarks/rnaspades.csv'
	shell: 'rnaspades.py --threads {threads} -s {input.R1} -o {params.output_folder}'


rule SE_soapdenovotrans:
	input:
		fastpJson = rules.SE_fastp.output.json,
		reads = [rules.SE_fastp.output.R1]
	output: 'output/{sample}/SE_assemblies/soapdenovotrans/soapdenovotrans.fasta'
	threads: config['threads']
	params:
		configSample = '{sample}',
		output_folder = 'output/{sample}/SE_assemblies/soapdenovotrans',
		orientation = orientation_soapdenovotrans
	conda: '../envs/soapdenovotrans.yml'
	benchmark: 'output/{sample}/benchmarks/soapdenovotrans.csv'
	script: '../scripts/soapdenovotrans.py'

rule SE_trinity:
	input:
		R1 = rules.SE_fastp.output.R1
	output: 'output/{sample}/SE_assemblies/trinity/Trinity.fasta'
	params: 
		orientation = lambda wildcards: config['samples'][wildcards.sample]['orientation'],
		output_folder = 'output/{sample}/SE_assemblies/trinity',
		memory = config['memory']
	threads: config['threads']
	conda: '../envs/trinity.yml'
	benchmark: 'output/{sample}/benchmarks/trinity.csv'
	shell: 'Trinity --seqType fq --max_memory {params.memory} --single {input.R1} --CPU {threads} --SS_lib_type {params.orientation} --output {params.output_folder}'

rule SE_rename_contigs:
	input:
		rnaspades = rules.SE_rnaspades.output,
		soapdenovotrans = rules.SE_soapdenovotrans.output,
		trinity = rules.SE_trinity.output
	output:
		rnaspades = 'output/{sample}/SE_assemblies/final/rnaspades_solo.fasta',
		soapdenovotrans = 'output/{sample}/SE_assemblies/final/soapdenovotrans.fasta',
		trinity = 'output/{sample}/SE_assemblies/final/trinity.fasta'
	script: '../scripts/contig_unique_header.py'

rule SE_merge:
	input:
		rnaspades = rules.SE_rename_contigs.output.rnaspades,
		soapdenovotrans = rules.SE_rename_contigs.output.soapdenovotrans,
		trinity = rules.SE_rename_contigs.output.trinity
	output: 'output/{sample}.SE.fasta'
	params: dir = directory('output/{sample}/SE_assemblies/final')
	shell: 'cat {params.dir}/*.fasta > {output}'
