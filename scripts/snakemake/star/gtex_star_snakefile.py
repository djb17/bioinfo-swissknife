import pandas as pd
import os

configfile: "env/star.yml"

sample_list = config["sample_dict"][config["curr_batch"]]
sample_index = config["curr_file"]

gtex_fqdir = config["gtex_fqdir"]
star_outdir = config["star_outdir"]

samples = pd.read_table(sample_list[sample_index])["file"]

rule all:
	input: # expected outputs from the workflow
		expand(star_outdir+"/{sample}/SJ.out.tab", sample=samples),
		expand(star_outdir+"/{sample}/aligned/Aligned.sortedByCoord.out.bam", sample=samples)

rule star_pass1:
	input:
		fq1 = gtex_fqdir+"/{sample}.R1.fastq",
		fq2 = gtex_fqdir+"/{sample}.R2.fastq"
	output:
		sj = star_outdir+"/{sample}/SJ.out.tab",
	params:
		star_genomedir = config["star_genomedir"],
		star_sampledir = lambda wildcards: star_outdir+"/"+wildcards.sample+"/", # specify different output path for each sample
		star_tempdir = lambda wildcards: star_outdir+"/"+wildcards.sample+"/tmp/" # specify different tempdir path for each sample

	shell:
		"""
		echo {params.star_sampledir} && echo {output.sj} && \
		STAR --genomeDir {params.star_genomedir} \
			--runThreadN 8 --readFilesCommand cat \
			--readFilesIn {input.fq2},{input.fq1} \
			--outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 \
			--alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 \
			--genomeLoad NoSharedMemory --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 \
			--sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMtype None --outSAMmode None \
			--outTmpDir {params.star_tempdir} \
			--outFileNamePrefix {params.star_sampledir}
		"""
rule star_pass2:
	input:
		fq1 = gtex_fqdir+"/{sample}.R1.fastq",
		fq2 = gtex_fqdir+"/{sample}.R2.fastq",
		sj = rules.star_pass1.output.sj
	output:
		bam = star_outdir+"/{sample}/aligned/Aligned.sortedByCoord.out.bam",
		sa = temp(star_outdir+"/{sample}/intermediate/SA"), # SA, SAindex, Genome files are unnecessary once alignment is finished
		si = temp(star_outdir+"/{sample}/intermediate/SAindex"),
		gn = temp(star_outdir+"/{sample}/intermediate/Genome")
	params:
		intermediate = directory(star_outdir+"/{sample}/intermediate/"),
		aligned = directory(star_outdir+"/{sample}/aligned/"),
		fasta = config["star_fasta"]
	shell:
		"""
		STAR --runMode genomeGenerate --genomeDir {params.intermediate} \
			--genomeFastaFiles {params.fasta} \
			--sjdbOverhang 100 --runThreadN 8 --sjdbFileChrStartEnd {input.sj} && \
		STAR --genomeDir {params.intermediate} \
			--readFilesIn {input.fq2},{input.fq1} \
			--runThreadN 8 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 \
			--alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 \
			--genomeLoad NoSharedMemory --limitBAMsortRAM 0 --readFilesCommand cat --outFilterMatchNminOverLread 0.33 \
			--outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif \
			--outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
			--outSAMheaderHD @HD VN:1.4 --outFileNamePrefix {params.aligned}
		"""
