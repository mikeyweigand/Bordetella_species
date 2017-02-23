#!/bin/bash

#set values
k=15;

#count kmers using jellyfish
now=$(date +%Y%m%d);
log=$(date +%Y%m%d%H%M);
if [ ! -d "results/jellyfish/$now" ];then
	mkdir results/jellyfish/$now;
fi
date > results/jellyfish/$now/$log.log;
echo ~/Tools/jellyfish-2.2.6/bin/jellyfish >> results/jellyfish/$now/$log.log;
~/Tools/jellyfish-2.2.6/bin/jellyfish --version >> results/jellyfish/$now/$log.log;
echo "k="$k | tee -a results/jellyfish/$now/$log.log;

for i in $(ls data/genomes/fasta/*fasta);
do
	echo $i | tee -a results/jellyfish/$now/$log.log;
	o=$(basename $i _Illumina_final.fasta);

	~/Tools/jellyfish-2.2.6/bin/jellyfish bc -m $k -s 5M -t 16 \
	-o results/jellyfish/$now/$o\-$k\mer.bc $i;
	~/Tools/jellyfish-2.2.6/bin/jellyfish count -m $k -s 5M -t 16 \
	--bc results/jellyfish/$now/$o\-$k\mer.bc $i \
	-o results/jellyfish/$now/$o\-$k\mer.jf;

	~/Tools/jellyfish-2.2.6/bin/jellyfish histo results/jellyfish/$now/$o\-$k\mer.jf \
	-o results/jellyfish/$now/$o\-$k\mer.hist;

	for j in $(~/Tools/jellyfish-2.2.6/bin/jellyfish dump -c -t results/jellyfish/$now/$o\-$k\mer.jf | cut -f2); \
	do echo -e $o"\t"$j; done > results/jellyfish/$now/$o\-$k\mer.dump;

done

#ggplot
#cat Bb_I328-15mer.dump Bhi_F582-15mer.dump Bho_F613-15mer.dump Bp_J021-15mer.dump Bpp_H904-15mer.dump Bsp_H567-15mer.dump Bt_F581-15mer.dump > dump3;
#dump03 = read.table("dump3");
#qplot( V2, data = dump03, geom = "density" , fill=V1, alpha=.3, ylim=c(0,0.02)) + theme_bw( ) + scale_x_continuous(breaks=seq(0,150,10))
#qplot( V2, data = dump03, geom = "density" , fill=V1, alpha=.2) + scale_x_continuous(breaks=seq(0,150,10),limits=c(10,150)) + coord_cartesian(ylim=c(0,0.4));
