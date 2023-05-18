import pandas as pd
import json
import os

"""
MuTect2 Variant Calling Workflow

This code implements a Snakemake workflow for running the MuTect2 variant calling tool
on paired tumor-normal samples. The workflow consists of several rules that define
the steps of the variant calling pipeline.

The workflow follows these main steps:

1. Loading pairing information from a JSON file.
2. Retrieving tumor samples based on the pairings.
3. Defining target files for the workflow.
4. Running MuTect2 on each sample.
5. Generating pileup summaries for inferring contamination.
6. Calculating contamination using the pileup summaries.
7. Correcting read orientation biases for accurate variant calling.
8. Refining variant calls from MuTect2 using customizable thresholds.
"""

configfile: "env/mutect2_paired.yml"

# Load pairing information from JSON file
with open(config['sample_dict'][config['curr_bait']]) as j:
	pairings = json.load(j)
normals = pd.read_table(config['sample_file'][config['curr_bait']][config['curr_file']], sep='\t')['file'].tolist()

# Retrieve tumor samples based on the pairings
tumors = [pairings[i] for i in normals]

curr_indir = config['indir'][config['curr_cohort']]
curr_outdir = config['outdir'][config['curr_cohort']]
curr_ti = config['target_int'][config['curr_bait']]

# Rule defining the target files for the workflow
rule all:
	input:
		expand(curr_outdir+"/{normal}.vcf.gz", normal=normals),
		expand([curr_outdir+"/{normal}.normal.gps.table", curr_outdir+"/{normal}.tumor.gps.table"], normal=normals),
		expand(curr_outdir+"/{normal}.ct.table", normal=normals),
		expand(curr_outdir+"/{normal}.seg.table", normal=normals),
		expand(curr_outdir+"/{normal}.ro_mod.tar.gz", normal=normals),
		expand(curr_outdir+"/{normal}.fv.vcf.gz", normal=normals)

def normalBAM(wildcards): # Function to retrieve normal BAM files to be processed in the matched workflow
	return(curr_indir+"/{normal}.bam")

# Rule for running MuTect2
rule gatk_mt2:
	input: 
		n_bam = normalBAM
	output:
		vcf = curr_outdir+"/{normal}.vcf.gz",
		f1r2 = curr_outdir+"/{normal}.f1r2.tar.gz",
		n_name = temp(curr_outdir+'/{normal}.normal.name'), # Temporary files with sample names
		t_name = temp(curr_outdir+'/{normal}.tumor.name')
	params:
		gatk_dir = config['gatk'],
		tmp_dir = config['gatk_tmp'],
		ref_gen = config['ref_genome'],
		target_int = curr_ti,
		germ_ref = config['germ_ref'],
		pon = config['pon'],
		mem = config['memory'],
		t_bam = lambda wildcards: curr_indir+'/'+pairings[wildcards.normal]+'.bam',
		logname = curr_outdir+"/{normal}.mutect2.log"
	message:
		"""
		Running Mutect2:
		Normal BAM: {input.n_bam}
		Tumor BAM: {params.t_bam}
		"""
	shell:
		"""
		# Extract sample names from BAM files
		{params.gatk_dir}/gatk GetSampleName --verbosity ERROR -I {input.n_bam} -O {output.n_name}
		{params.gatk_dir}/gatk GetSampleName --verbosity ERROR -I {params.t_bam} -O {output.t_name}

		# Run MuTect2
		/usr/bin/time -v {params.gatk_dir}/gatk \
		--java-options \"-Xmx{params.mem} -Djava.io.tmpdir={params.tmp_dir}\" \
		Mutect2 --verbosity INFO -R {params.ref_gen} -L {params.target_int} \
		--germline-resource {params.germ_ref} \
		--panel-of-normals {params.pon} \
		--genotype-germline-sites \
		--interval-padding 75 \
		-I {input.n_bam} -normal $(cat {output.n_name}) \
		-I {params.t_bam} -tumor $(cat {output.t_name}) \
		--f1r2-tar-gz {output.f1r2} \
		-O {output.vcf} &> {params.logname}
		"""

# Rule for generating pileup summaries to infer contamination
rule gatk_gps:
	input:
		n_bam = normalBAM,
		vcf_check = curr_outdir+'/{normal}.vcf.gz'
	output:
		n_gps = curr_outdir+"/{normal}.normal.gps.table",
		t_gps = curr_outdir+"/{normal}.tumor.gps.table"
	params:
		gatk_dir = config['gatk'],
		tmp_dir = config['gatk_tmp'],
		target_int = curr_ti,
		germ_bi = config['germ_bi'],
		mem = config['memory'],
		t_bam = lambda wildcards: curr_indir+'/'+pairings[wildcards.normal]+'.bam',
		logname = curr_outdir+"/{normal}.gps.log"
	shell:
		"""
		# Generate pileup summaries for normal and tumor samples
		{params.gatk_dir}/gatk --java-options \"-Xmx{params.mem} -Djava.io.tmpdir={params.tmp_dir}\" GetPileupSummaries \
		--variant {params.germ_bi} --intervals {params.target_int} --interval-padding 75 \
		--input {input.n_bam} --output {output.n_gps} &> {params.logname}

		{params.gatk_dir}/gatk --java-options \"-Xmx{params.mem} -Djava.io.tmpdir={params.tmp_dir}\" GetPileupSummaries \
		--variant {params.germ_bi} --intervals {params.target_int} --interval-padding 75 \
		--input {params.t_bam} --output {output.t_gps} &>> {params.logname}
		"""

# Rule for calculating contamination
rule gatk_cc:
	input:
		n_gps = rules.gatk_gps.output.n_gps,
		t_gps = rules.gatk_gps.output.t_gps
	output:
		ct = curr_outdir+"/{normal}.ct.table",
		seg = curr_outdir+"/{normal}.seg.table"
	params:
		gatk_dir = config['gatk'],
		tmp_dir = config['gatk_tmp'],
		mem = config['memory'],
		logname = curr_outdir+"/{normal}.ct.log"
	shell:
		"""
		# Calculate contamination using GetPileupSummaries outputs
		{params.gatk_dir}/gatk --java-options \"-Xmx{params.mem} -Djava.io.tmpdir={params.tmp_dir}\" CalculateContamination \
			--input {input.t_gps} --matched {input.n_gps} \
			--output {output.ct} --tumor-segmentation {output.seg} &> {params.logname}
		"""

# Rule for correcting read orientation biases for accurate variant calling
rule gatk_lro:
	input:
		f1r2 = rules.gatk_mt2.output.f1r2
	output:
		ro_mod = curr_outdir+"/{normal}.ro_mod.tar.gz"
	params:
		gatk_dir = config['gatk'],
		tmp_dir = config['gatk_tmp'],
		mem = config['memory'],
		logname = curr_outdir+"/{normal}.lro.log"
	shell:
		"""
			# Learn read orientation model
			{params.gatk_dir}/gatk --java-options \"-Xmx{params.mem} -Djava.io.tmpdir={params.tmp_dir}\" LearnReadOrientationModel \
			-I {input.f1r2} -O {output.ro_mod} &> {params.logname}
		"""

# Rule for refining variant calls from MuTect2 using customizable thresholds
rule gatk_fv:
	input:
		vcf = rules.gatk_mt2.output.vcf,
		cc = rules.gatk_cc.output.ct,
		seg = rules.gatk_cc.output.seg,
		ro_mod = rules.gatk_lro.output.ro_mod
	output:
		vcf = curr_outdir+"/{normal}.fv.vcf.gz"
	params:
		gatk_dir = config['gatk'],
		tmp_dir = config['gatk_tmp'],
		ref_gen = config['ref_genome'],
		mem = config['memory']
	shell:
		"""
		{params.gatk_dir}/gatk --java-options \"-Xmx{params.mem} -Djava.io.tmpdir={params.tmp_dir}\" FilterMutectCalls \
			-V {input.vcf} -R {params.ref_gen} --tumor-segmentation {input.seg} --contamination-table {input.cc} \
			--ob-priors {input.ro_mod} -O {output.vcf}
		"""