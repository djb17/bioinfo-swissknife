#!/usr/bin/env bash

AFFY_ANNOT=/data/home/data/affy_annotationFiles
AFFY_LIB=$AFFY_ANNOT/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles
REFGENOME=/data/home/data/referenceGenome/hg37.fasta.gz

if [[ $# != 2 ]]; then
	echo -e "\nUsage: $(basename $0) [cel-list] [out-directory]\n";
	exit 1;
fi

cel_list=$1; out_dir=$2;

if [[ ! -e $cel_list ]]; then
	echo -e "\nFile $cel_list does not exist.\n";
	exit 1;
elif [[ ! -d $out_dir ]]; then
	echo -e "\nDirectory $out_dir does not exist. Creating a new directory.\n";
	mkdir $out_dir;
fi

echo -e "\nAnnotation/Library & reference genome paths:"
echo -e "\nAFFY_ANNOT\t$AFFY_ANNOT\nAFFY_LIB\t$AFFY_LIB\nREFGENOME\t$REFGENOME\n"

# Call genotypes using Analysis Power Tools
apt-probeset-genotype \
--out-dir $out_dir \
--analysis-files-path $AFFY_ANNOT \
--read-models-brlmmp /data/home/data/affy_annotationFiles/GenomeWideSNP_6.generic_prior.txt \
--xml-file $AFFY_ANNOT/GenomeWideSNP_6.apt-probeset-genotype.AxiomGT1.xml \
-c $AFFY_LIB/GenomeWideSNP_6.cdf \
--chrX-probes $AFFY_LIB/GenomeWideSNP_6.chrXprobes \
--chrY-probes $AFFY_LIB/GenomeWideSNP_6.chrYprobes \
--special-snps $AFFY_LIB/GenomeWideSNP_6.specialSNPs \
--cel-files $cel_list \
--summaries \
--write-models &&

# Create VCF from the genotype calls
mkdir ./bcftools-sort.XXXXXX && 
bcftools +affy2vcf --no-version -Ou --fasta-ref $REFGENOME --annot $AFFY_ANNOT/GenomeWideSNP_6.na35.annot.csv \
--report $out_dir/AxiomGT1.report.txt \
--calls $out_dir/AxiomGT1.calls.txt \
--confidences $out_dir/AxiomGT1.confidences.txt \
--snp-posteriors $out_dir/AxiomGT1.snp-posteriors.txt \
--summary $out_dir/AxiomGT1.summary.txt | \
bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
bcftools +fixref --no-version -Ou -e 'REF="N" || ALT="N"' -- -f $REFGENOME -m swap --flip-baf | \
bcftools norm --no-version -Ob -o $out_dir/$(basename $cel_list .txt).bcf -c x -f $REFGENOME &&
rm -rf ./bcftools-sort.XXXXXX