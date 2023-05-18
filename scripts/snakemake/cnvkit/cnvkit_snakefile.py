import pandas as pd
import os, glob

"""
CNVkit Workflow for Calling Copy Number Variations (CNVs)
This code represents a workflow that utilizes CNVkit, a tool for detecting copy number variations (CNVs) from sequencing data.
The workflow is configured based on the 'cnvkit.yml' file.

The code reads sample information from the provided sample file and sets the input and output directories accordingly.

Depending on the specified function (e.g., 'coverage', 'reference', 'fix', 'segment', or 'call'), different rules are defined to perform specific steps in the CNV calling process.
 - 'coverage' function computes coverage from BAM files, generating target coverage and antitarget coverage files.
 - 'reference' function creates a reference file for subsequent CNV analysis, utilizing normalized coverage files obtained from the 'coverage' step.
 - 'fix' function applies corrections to target coverage and antitarget coverage files using the created reference file.
 - 'segment' function segments the corrected coverage data.
"""

configfile: 'env/cnvkit.yml'

samples = pd.read_table(config['sample_file'][config['curr_cohort']][config['curr_file']], sep='\t')['file'].tolist()

curr_indir = config['indir'][config['curr_cohort']]
curr_outdir = config['outdir'][config['curr_cohort']]

if config['function'] == 'coverage':
    rule all:
        input:
            expand(curr_outdir+'/{sample}.targetcoverage.cnn', sample=samples),
            expand(curr_outdir+'/{sample}.antitargetcoverage.cnn', sample=samples)

    rule cnvkit_cov: # compute coverage from BAM files
        input:
            bam = curr_indir+'/{sample}.bam'
        output:
            tcov = curr_outdir+'/{sample}.targetcoverage.cnn',
            acov = curr_outdir+'/{sample}.antitargetcoverage.cnn'
        params:
            target = config['target'][config['curr_bait']],
            antitarget = config['antitarget'][config['curr_bait']],
            tlog = curr_outdir+'/logs/{sample}.targetcoverage.log',
            alog = curr_outdir+'/logs/{sample}.antitargetcoverage.log'
        shell:
            """
            /usr/bin/time -v cnvkit.py coverage {input.bam} {params.target} -p 2 -o {output.tcov} &> {params.tlog} && \
            /usr/bin/time -v cnvkit.py coverage {input.bam} {params.antitarget} -p 2 -o {output.acov} &> {params.alog}
            """

if config['function'] == 'reference':
    def get_gcov():
        return(glob.glob(curr_outdir+'/g*.cnn'))

    rule all:
        input:
            curr_outdir+'/'+config['curr_refname'][config['curr_cohort']]+'.cnn'
    rule cnvkit_ref:
        input:
            norm_cov = expand(get_gcov())
        output:
            cnv_ref = curr_outdir+'/'+config['curr_refname'][config['curr_cohort']]+'.cnn'
        params:
            ref_genome = config['ref_genome']
        shell:
            """
            /usr/bin/time -v cnvkit.py reference {input.norm_cov} -f {params.ref_genome} -o {output.cnv_ref}
            """

if config['function'] == 'fix':
    rule all:
        input:
            expand(curr_outdir+'/{sample}.cnr', sample=samples)
    rule cnvkit_fix:
        input:
            tcov = curr_outdir+'/{sample}.targetcoverage.cnn',
            acov = curr_outdir+'/{sample}.antitargetcoverage.cnn'
        output:
            cnr = curr_outdir+'/{sample}.cnr'
        params:
            cnv_ref = curr_outdir+'/'+config['curr_refname'][config['curr_cohort']]+'.cnn',
            log = curr_outdir+'/{sample}.cnr.log'
        shell:
            """
            /usr/bin/time -v cnvkit.py fix {input.tcov} {input.acov} {params.cnv_ref} -o {output.cnr} &> {params.log}
            """

if config['function'] == 'segment':
    rule all:
        input:
            expand(curr_outdir+'/{sample}.cns', sample=samples)
    rule cnvkit_seg:
        input:
            cnr = curr_outdir+'/{sample}.cnr'
        output:
            cns = curr_outdir+'/{sample}.cns'
        params:
            log = curr_outdir+'/{sample}.cns.log'
        shell:
            """
            /usr/bin/time -v cnvkit.py segment -m cbs --drop-low-coverage -p {threads} \
                {input.cnr}  -o {output.cns} &> {params.log}
            """

if config['function'] == 'call':
    metadata = pd.read_csv(config['metadata'][config['curr_cohort']], sep='\t', low_memory=False)
    rule all:
        input:
            expand(curr_outdir+'/{sample}.call.cns', sample=samples)
    rule cnvkit_call:
        input:
            cns = curr_outdir+'/{sample}.cns'
        output:
            cns = curr_outdir+'/{sample}.call.cns'
        params:
            log = curr_outdir+'/{sample}.call.cns.log',
            purity = lambda wildcards: metadata[metadata.Sampleid == wildcards.sample.split('.')[0]].Purity.item(),
            ploidy_noround = lambda wildcards: metadata[metadata.Sampleid == wildcards.sample.split('.')[0]].Ploidy.item(),
        shell:
            # """
            """
            # with ploidy & purity
            echo -e "purity: " {params.purity} "ploidy: " {params.ploidy_noround}","{params.ploidy_round}
            /usr/bin/time -v cnvkit.py call -x female --drop-low-coverage \
                {input.cns} -o {output.cns} --purity {params.purity} --ploidy {params.ploidy_round} &> {params.log}
            """