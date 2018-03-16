#!/bin/bash

now=$(date +%Y%m%d);
log=$(date +%Y%m%d%H%M);
if [ ! -d "results/kSNP/$now" ];then
	mkdir results/kSNP/$now;
fi
date > results/kSNP/$now/$log.log;

sp=("Ba" "Bb" "Bhi" "Bho" "Bpp");

for i in "${sp[@]}"
do
	echo $i;
	fasta=/home/yrh8/Documents/Bordetella_species/data/genomes/fasta/$i\_*fasta;

	echo -n >  results/kSNP/$now/$i\_fasta.list;

	for f in $fasta
	do
		n=$(basename $f _Illumina_final.fasta);

		echo -e $f"\t"$n >> results/kSNP/$now/$i\_fasta.list;

	done

	/home/yrh8/Tools/kSNP3/kSNP3.0/kSNP3 -k 23 -core -CPU 4 \
	-in results/kSNP/$now/$i\_fasta.list \
	-outdir results/kSNP/$now/$i\_output


done
