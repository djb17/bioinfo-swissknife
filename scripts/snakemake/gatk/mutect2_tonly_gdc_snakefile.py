import pandas as pd
import os

"""
Tumor-only somatic calling using GATK's MuTect2.
This workflow attempts to perform variant calling on individual chromosomes.
[Note]: This was a test on the lab machine to see how much faster it would be.
        Should re-write it such that it doesn't use 'xargs -n1 -P' to run multiple instances.
"""

configfile: 'env/mutect2_tonly.yml'

tumors = pd.read_table(config['sample_file'][config['curr_bait']][config['curr_file']], sep='\t')['file'].tolist()
chrs = pd.read_csv(config['chr_list'], header=None).iloc[:,0].tolist()
curr_indir = config['indir'][config['curr_cohort']]
curr_outdir = config['outdir'][config['curr_cohort']]

rule all:
	input:
		expand(curr_outdir+'/{tumor}.bait_bias_summary_metrics.txt', tumor=tumors),
		expand(curr_outdir+'/{tumor}.error_summary_metrics.txt', tumor=tumors),
		expand(curr_outdir+'/{tumor}.pre_adapter_detail_metrics.txt', tumor=tumors),
		expand(curr_outdir+'/{tumor}.pre_adapter_summary_metrics.txt', tumor=tumors),
		expand(curr_outdir+"/{tumor}.gdc.gps.table", tumor=tumors),
		expand(curr_outdir+"/{tumor}.ct.table", tumor=tumors),
		expand(curr_outdir+"/{tumor}.{chr}.gdc.vcf.gz", tumor=tumors, chr=chrs),
		expand(curr_outdir+"/{tumor}.merged.gdc.vcf.gz", tumor=tumors),
		expand(curr_outdir+"/{tumor}.merged.sorted.gdc.vcf.gz", tumor=tumors),
		expand(curr_outdir+"/{tumor}.filtered.gdc.vcf.gz", tumor=tumors),
		expand(curr_outdir+'/{tumor}.raw_somatic.gdc.vcf.gz', tumor=tumors),
		expand(curr_outdir+'/{tumor}.vep_somatic.gdc.vcf.gz', tumor=tumors)

# 1. Generate OXOG metrics:
rule gatk_csam:
	input: 
		t_bam = curr_indir+'/{tumor}.bam'
	output:
		bbias = curr_outdir+'/{tumor}.bait_bias_summary_metrics.txt',
		ersum = curr_outdir+'/{tumor}.error_summary_metrics.txt',
		adapd = curr_outdir+'/{tumor}.pre_adapter_detail_metrics.txt',
		adaps = curr_outdir+'/{tumor}.pre_adapter_summary_metrics.txt'
	params:
		gatk_dir = config['gatk'][config['curr_gatk']],
		tmp_dir = config['gatk_tmp'],
		ref_gen = config['ref_genome'],
		target_int = config['target_int'][config['curr_bait']],
		germ_ref = config['germ_ref'],
		pon = config['tcga_pon'],
		mem = config['memory'],
		sam_name = curr_outdir+"/{tumor}",
		logname = curr_outdir+"/{tumor}.gdc.csam.log"
	threads: 5
	message:
		"Running CollectSequencingArtifactMetrics {input.t_bam} with {threads} threads..."
	shell:
		"""
		{params.gatk_dir}/gatk --java-options \"-Xmx{params.mem} -Djava.io.tmpdir={params.tmp_dir}\" \
			CollectSequencingArtifactMetrics -R {params.ref_gen} -I {input.t_bam} -O {params.sam_name} \
			--FILE_EXTENSION .txt --VALIDATION_STRINGENCY LENIENT &> {params.logname}
		"""

# 2. Generate pileup summaries on tumor sample:
rule gatk_gps:
	input:
		step1_check = rules.gatk_csam.output.bbias,
		t_bam = curr_indir+'/{tumor}.bam'
	output:
		t_gps = curr_outdir+"/{tumor}.gdc.gps.table"
	params:
		gatk_dir = config['gatk'][config['curr_gatk']],
		tmp_dir = config['gatk_tmp'],
		ref_gen = config['ref_genome'],
		target_int = config['target_int'][config['curr_bait']],
		germ_bi = config['germ_bi'],
		mem = config['memory'],
		logname = curr_outdir+"/{tumor}.gdc.gps.log"
	threads: 5
	shell:
		"""
		{params.gatk_dir}/gatk --java-options \"-Xmx{params.mem} -Djava.io.tmpdir={params.tmp_dir}\" GetPileupSummaries \
		-V {params.germ_bi} --intervals {params.target_int} -R {params.ref_gen} \
		-I {input.t_bam} --output {output.t_gps} &> {params.logname}
		"""

# 3. Calculate contamination on tumor sample
rule gatk_cc:
	input:
		t_gps = rules.gatk_gps.output.t_gps
	output:
		ct = curr_outdir+"/{tumor}.ct.table"
	params:
		gatk_dir = config['gatk'][config['curr_gatk']],
		tmp_dir = config['gatk_tmp'],
		mem = config['memory'],
		logname = curr_outdir+"/{tumor}.gdc.ct.log"
	threads: 5
	shell:
		"""
		{params.gatk_dir}/gatk --java-options \"-Xmx{params.mem} -Djava.io.tmpdir={params.tmp_dir}\" CalculateContamination \
			-I {input.t_gps} -output {output.ct}  &> {params.logname}
		"""

# 4. Find tumor sample name from BAM
rule gatk_mt2: # MuTect2 tumor-only mode
	input: 
		step3_check = rules.gatk_cc.output.ct,
		t_bam = curr_indir+'/{tumor}.bam'
	output:
		vcf = expand(curr_outdir+"/{{tumor}}.{ext}.gdc.vcf.gz", ext=chrs),
		name = curr_outdir+"/{tumor}.name"
	params:
		gatk_dir = config['gatk'][config['curr_gatk']],
		tmp_dir = config['gatk_tmp'],
		ref_gen = config['ref_genome'],
		target_int = config['target_int'][config['curr_bait']],
		target_bed = config['target_bed'][config['curr_bait']],
		chr_list = config['chr_list'],
		germ_ref = config['germ_ref'],
		pon = config['tcga_pon'],
		mem = config['memory'],
		vcf_chr = curr_outdir+"/{tumor}",
		logname_chr = curr_outdir+"/mt2_log/{tumor}.mt2_tonly",
		logname = curr_outdir+"/{tumor}.mt2_tonly.log"
	threads: 5
	message:
		"Running Mutect2 for {input.t_bam} with {threads} threads"
	shell:
		"""
		# Changed for testing: Mutect2 -R {params.ref_gen} -L $0 -tumor $(cat {output.name}) --native-pair-hmm-threads 1
		{params.gatk_dir}/gatk --java-options \"-Xmx{params.mem} -Djava.io.tmpdir={params.tmp_dir}\" \
		GetSampleName -I {input.t_bam} -O {output.name}
		
		cat {params.chr_list} | xargs -n1 -P{threads} bash -c '\
		/usr/bin/time -v {params.gatk_dir}/gatk --java-options \"-Xmx{params.mem} -Djava.io.tmpdir={params.tmp_dir}\" \
		Mutect2 -R {params.ref_gen} -L {params.target_bed}_$0.bed -tumor $(cat {output.name}) --native-pair-hmm-threads 1\
		--germline-resource {params.germ_ref} --panel-of-normals {params.pon} \
		--af-of-alleles-not-in-resource 2.5e-06 --genotype-germline-sites \
		-I {input.t_bam} -O {params.vcf_chr}.$0.gdc.vcf.gz &> {params.logname_chr}.$0.log'
		"""

# 5. Merge vcfs divided by chromosomes
rule gatk_mvcf: 
	input:
		step6_check = rules.gatk_mt2.output.vcf
	output:
		merged_vcf = curr_outdir+"/{tumor}.merged.gdc.vcf.gz"
	params:
		gatk_dir = config['gatk'][config['curr_gatk']],
		vcfs = ' '.join(['-I ' + i for i in expand(curr_outdir+'/{{tumor}}.{ext}.gdc.vcf.gz ', ext=chrs)])
	threads: 5
	shell:
		"""
		{params.gatk_dir}/gatk MergeVcfs {params.vcfs} -O {output.merged_vcf}
		"""

# 6. Sort VCF with Picard
rule gatk_sort:
	input:
		merged_vcf = rules.gatk_mvcf.output.merged_vcf
	output:
		sorted_vcf = curr_outdir+"/{tumor}.merged.sorted.gdc.vcf.gz"
	params:
		gatk_dir = config['gatk'][config['curr_gatk']],
		ref_dict = config['ref_genome_dict']
	threads: 5
	shell:
		"""
		{params.gatk_dir}/gatk SortVcf --SEQUENCE_DICTIONARY {params.ref_dict} --OUTPUT={output.sorted_vcf} \
			--INPUT {input.merged_vcf} --CREATE_INDEX true
		"""

# 7. Filter variant calls from MuTect
rule gatk_filter:
	input:
		sorted_vcf = rules.gatk_sort.output.sorted_vcf,
		ct = rules.gatk_cc.output.ct
	output:
		filtered_vcf = curr_outdir+"/{tumor}.filtered.gdc.vcf.gz"
	params:
		gatk_dir = config['gatk'][config['curr_gatk']],
		target_int = config['target_int'][config['curr_bait']]
	threads: 5
	shell:
		"""
		{params.gatk_dir}/gatk FilterMutectCalls -L {params.target_int} \
			-V {input.sorted_vcf} --contamination-table {input.ct} -O {output.filtered_vcf}
		"""

# 8. Filter variants by orientation bias
rule gatk_obias:
	input:
		filtered_vcf = rules.gatk_filter.output.filtered_vcf,
		adapd = rules.gatk_csam.output.adapd
	output:
		raw_vcf = curr_outdir+'/{tumor}.raw_somatic.gdc.vcf.gz'
	params:
		gatk_dir = config['gatk'][config['curr_gatk']],
		target_int = config['target_int'][config['curr_bait']],
		ref_gen = config['ref_genome']
	threads: 5
	shell:
		"""
		{params.gatk_dir}/gatk FilterByOrientationBias -L {params.target_int} -R {params.ref_gen} \
			-P {input.adapd} -V {input.filtered_vcf} -O {output.raw_vcf} -AM G/T -AM C/T
		"""