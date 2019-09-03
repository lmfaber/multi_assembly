############################
##### paired end reads #####
############################

def orientation_rnaspades(wildcards):
	orientation = config['samples'][wildcards.sample]['orientation']
	if orientation == 'FR':
		return '--pe1-fr'
	elif orientation == 'RF':
		return '--pe1-rf'
	else:
		print('No valid orientation.')

def orientation_soapdenovotrans(wildcards):
	orientation = config['samples'][wildcards.sample]['orientation']
	if orientation == 'FR':
		return '0'
	elif orientation == 'RF':
		return '1'
	else:
		print('No valid orientation.')

rule PE_header_correction:
	input: getPathOfPairedEndReads()
	output:
		R1 = temp('output/{sample}/PE_{sample}_1.fastq'),
		R2 = temp('output/{sample}/PE_{sample}_2.fastq')
	script: '../scripts/parser.py'

rule PE_fastp:
	input: 
		R1 = rules.PE_header_correction.output.R1,
		R2 = rules.PE_header_correction.output.R2
	output:
		R1 = 'output/{sample}/fastp/{sample}_1.fastq',
		R2 = 'output/{sample}/fastp/{sample}_2.fastq',
		json = 'output/{sample}/fastp/{sample}_PE.json',
		html = 'output/{sample}/fastp/{sample}_PE.html'
	conda: '../envs/fastp.yml'
	benchmark: 'output/{sample}/benchmarks/fastp.csv'
	shell: 'fastp --in1 {input.R1} --in2 {input.R2} -o {output.R1} -O {output.R2} --json {output.json} --html {output.html} -c --low_complexity_filter --detect_adapter_for_pe'
		
rule PE_rnaspades:
	input:
		R1 = rules.PE_fastp.output.R1,
		R2 = rules.PE_fastp.output.R2
	output: 'output/{sample}/PE_assemblies/rnaspades/transcripts.fasta'
	threads: config['threads']
	params: 
		orientation = orientation_rnaspades,
		output_folder = 'output/{sample}/PE_assemblies/rnaspades'
	conda: '../envs/rnaspades.yml'
	benchmark: 'output/{sample}/benchmarks/rnaspades.csv'
	shell: 'rnaspades.py --threads {threads} -1 {input.R1} -2 {input.R2} {params.orientation} -o {params.output_folder}'

rule PE_soapdenovotrans:
	input:
		reads = [rules.PE_fastp.output.R1, rules.PE_fastp.output.R2],
		fastpJson = rules.PE_fastp.output.json
	output: 'output/{sample}/PE_assemblies/soapdenovotrans/soapdenovotrans.fasta'
	threads: config['threads']
	params:
		configSample = '{sample}',
		output_folder = 'output/{sample}/PE_assemblies/soapdenovotrans',
		orientation = orientation_soapdenovotrans
	conda: '../envs/soapdenovotrans.yml'
	benchmark: 'output/{sample}/benchmarks/soapdenovotrans.csv'
	script:	'../scripts/soapdenovotrans.py'


rule PE_trinity:
	input:
		R1 = rules.PE_fastp.output.R1,
		R2 = rules.PE_fastp.output.R2
	output: 'output/{sample}/PE_assemblies/trinity/Trinity.fasta'
	threads: config['threads']
	params:
		orientation = lambda wildcards: config['samples'][wildcards.sample]['orientation'],
		memory = config['memory'],
		output_folder = 'output/{sample}/PE_assemblies/trinity'
	conda: '../envs/trinity.yml'
	benchmark: 'output/{sample}/benchmarks/trinity.csv'
	shell: 'Trinity --seqType fq --max_memory {params.memory} --left {input.R1} --right {input.R2} --CPU {threads} --SS_lib_type {params.orientation} --output {params.output_folder}'

rule PE_rename_contigs:
	input:
		rnaspades = rules.PE_rnaspades.output,
		soapdenovotrans = rules.PE_soapdenovotrans.output,
		trinity = rules.PE_trinity.output
	output:
		rnaspades = 'output/{sample}/PE_assemblies/final/rnaspades_solo.fasta',
		soapdenovotrans = 'output/{sample}/PE_assemblies/final/soapdenovotrans.fasta',
		trinity = 'output/{sample}/PE_assemblies/final/trinity.fasta'
	script: '../scripts/contig_unique_header.py'

rule PE_merge:
	input:
		rnaspades = rules.PE_rename_contigs.output.rnaspades,
		soapdenovotrans = rules.PE_rename_contigs.output.soapdenovotrans,
		trinity = rules.PE_rename_contigs.output.trinity
	output: 'output/{sample}.PE.fasta'
	params: dir = directory('output/{sample}/PE_assemblies/final')
	shell: 'cat {params.dir}/*.fasta > {output}'
