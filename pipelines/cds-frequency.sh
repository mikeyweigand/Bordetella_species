#!/bin/bash

if [[ "$1" == "" || "$1" == "-h" ]] ; then
   echo "
   This script will run cdhit-est on a directory full of gene sequence multi-fasta files and output files to './results/cd-hit/YYYYMMDD/'
	 If necessary, it will first extract CDS sequences from the Genbank files using gbk2genes.pl (Scriptbox).
	 Expects cdhit-est in PATH.

   Usage: ./pipelines/cds-frequency.sh [folder] [genes]

   folder       Path to the folder containing the genbank files (*.gbk)
   genes        Path to the directory containing the gene files (*.fna)


   " >&2 ;
   exit 1 ;
fi ;



#prepare gene fasta files with simplified headers
if [ ! -d "$2" ];then
	mkdir $2;
fi

for i in $(ls $1/*gbk);
do
	echo $i;
	o=$(basename $i gbk)fna;

	if [ ! -f $2/$o ];then
		gbk2genes.pl -in $i -s > $2/$o;
	fi

done

#cluster genes by sequence identity using cd-hit
now=$(date +%Y%m%d);
log=$(date +%Y%m%d%H%M);
if [ ! -d "results/cd-hit/$now" ];then
	mkdir results/cd-hit/$now;
fi
date > results/cd-hit/$now/$log.log;
which cdhit-est >> results/cd-hit/$now/$log.log;
cdhit-est | head -1 >> results/cd-hit/$now/$log.log;

for i in $(ls $2/*fna);
do
	echo $i;
	o=$(basename $i fna)out;
	s=$(seq -s"," 300);

	if [ ! -f results/cd-hit/$now/$o ];then
		cdhit-est -i $i -d 0 \
		-o results/cd-hit/$now/$o \
		-c 0.95 -n 10 -l 11 -r 0 -G 1 -g 1 -b 20 \
		-s 0.9 -aL 0.90 -aS 0 -T 4 -M 2000 \
		>> results/cd-hit/$now/$log.log;

		clstr_sort_by < results/cd-hit/$now/$o.clstr > results/cd-hit/$now/$o.clstr.sorted

		./src/cdhit_plot_len1.pl results/cd-hit/$now/$o.clstr $s \
		1-99999 > results/cd-hit/$now/$o.clstr.hist;

	fi

done
