#!/bin/bash

#Description: Just a simple way to sync genome fasta and genbank files from an input list.

dir_fasta="/media/Genomics_Data/Assemblies and Reports/03.Illumina_polished/";
dir_gbk="/media/Genomics_Data/Annotations/Genbank/";
dir_faa="/media/Genomics_Data/Annotations/FASTA_proteins/"

while read line
do

	fasta=$(find "$dir_fasta" -type f -name "*final.fasta" | grep $line);
	gbk=$(find "$dir_gbk" -type f -name "*gbk" | grep $line);
	faa=$(find "$dir_faa" -type f -name "*faa" | grep $line);

	echo $line #"\n\t"$fasta"\n\t"$gbk;
	if [ -f "$fasta" ]; then
		echo -e "\t" $fasta;
		rsync "$fasta" ./data/genomes/fasta/ ;
	fi

	if [ -f "$gbk" ]; then
		echo -e "\t" $gbk;
		rsync "$gbk" ./data/genomes/genbank/;
	fi

	if [ -f "$faa" ]; then
		echo -e "\t" $faa;
		rsync "$faa" ./data/proteins/;
	fi


done < "${1:-/dev/stdin}"
