#!/bin/bash
#prepare gene fasta files with simplified headers
for i in $(ls data/genomes/genbank/*gbk);
do
	echo $i;
	o=$(basename $i gbk)fna;

	if [ ! -f data/genes/$o ];then
		gbk2genes.pl -in $i -s > data/genes/$o;
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

for i in $(ls data/genes/*fna);
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

#H567=read.table("H567_CP012334.out.clstr.hist",sep="\t",skip=2,nrows=299)
#H567[is.na(H567)] <- 0
#ggplot(H567, aes(x=V1,y=V3)) + theme_bw()+geom_bar(stat="identity")+coord_cartesian(ylim=c(0,15),xlim=c(0,70)) +xlab("copies")+ylab("num genes")+ggtitle("Bsp_H567")
