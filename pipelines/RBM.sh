#!/bin/bash
#prepare protein fasta files with simplified headers
for i in $(ls data/genomes/genbank/*gbk);
do
	echo $i;
	o=$(basename $i gbk)faa;

	if [ ! -f data/proteins/$o ];then
		gbk2genes.pl -in $i -p -s > data/proteins/$o;
	fi

done

for i in $(ls data/NCBI/*gbk);
do
	echo $i;
	o=$(basename $i gbk)faa;

	if [ ! -f data/NCBI/proteins/$o ];then
		gbk2genes.pl -in $i -p -s > data/NCBI/proteins/$o;
	fi

done
