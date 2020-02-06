#!/bin/bash

if [[ "$1" == "" || "$1" == "-h" ]] ; then
   echo "
   This script will run jellyfish on a directory full of genome sequence fasta files and output files to './results/jellyfish/YYYYMMDD/'

   -> Expects 'jellyfish' is in \$PATH.
   -> Final summary table depends on 'Table.merge.pl' (http://enve-omics.ce.gatech.edu/enveomics/. Must be in your PATH)

   Usage: ./pipelines/kmer-frequency.sh [k] [folder]

   k            Your favorite kmer length in bp (eg. 15)
   folder       Path to the directory containing the genome files (*.fasta)


   " >&2 ;
   exit 1 ;
fi ;

#where is your jellyfish?
#myjelly="~/Tools/jellyfish-2.2.6/bin/jellyfish"
myjelly="jellyfish";

#set values
k=15;
k=$1;

#count kmers using jellyfish
now=$(date +%Y%m%d);
log=$(date +%Y%m%d%H%M);
if [ ! -d "results/jellyfish/$now" ];then
	mkdir results/jellyfish/$now;
fi
date > results/jellyfish/$now/$log.log;
echo $myjelly >> results/jellyfish/$now/$log.log;
$myjelly --version >> results/jellyfish/$now/$log.log;
echo "k="$k | tee -a results/jellyfish/$now/$log.log;

for i in $(ls $2/*fasta);
do
	echo $i | tee -a results/jellyfish/$now/$log.log;
	o=$(basename $i .fasta);

	$myjelly bc -m $k -s 5M -t 16 \
	-o results/jellyfish/$now/$o\-$k\mer.bc $i;
	$myjelly count -m $k -s 5M -t 16 \
	--bc results/jellyfish/$now/$o\-$k\mer.bc $i \
	-o results/jellyfish/$now/$o\-$k\mer.jf;

	$myjelly histo results/jellyfish/$now/$o\-$k\mer.jf \
	-o results/jellyfish/$now/$o\-$k\mer.hist;

	for j in $($myjelly dump -c -t results/jellyfish/$now/$o\-$k\mer.jf | cut -f2); \
	do echo -e $o"\t"$j; done > results/jellyfish/$now/$o\-$k\mer.dump;

done

if [ ! -d "results/jellyfish/$now/hist" ];then
	mkdir results/jellyfish/$now/hist;
fi

for i in $(ls results/jellyfish/$now/*.hist | cut -f1 -d"_" | uniq);
do
	n=$(basename $i);
	Table.merge.pl -i " " results/jellyfish/$now/$n\_*.hist | sort -n > results/jellyfish/$now/hist/$n\-$k\mer.hist;

done
