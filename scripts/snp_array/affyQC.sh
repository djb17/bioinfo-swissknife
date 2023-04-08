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

timestamp=$(date +"%Y_%m_%d_%T")

apt-geno-qc \
	--cdf-file "$AFFY_LIB"/GenomeWideSNP_6.cdf \
	--qca-file "$AFFY_LIB"/GenomeWideSNP_6.r2.qca \
	--qcc-file "$AFFY_LIB"/GenomeWideSNP_6.r2.qcc \
	--cel-files "$cel_list" \
	--out-file $out_dir/qc_$timestamp.txt