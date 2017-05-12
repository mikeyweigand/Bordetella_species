#!/bin/bash
##prepare protein fasta files with simplified headers
for i in $(ls data/genomes/genbank/*gbk);
do
	o=$(basename $i gbk)faa;

	if [ ! -f data/proteins/$o ];then
		echo $i;
		gbk2genes.pl -in $i -p -s > data/proteins/$o;
	fi

done

for i in $(ls data/NCBI/*gbk);
do
	o=$(basename $i gbk)faa;

	if [ ! -f data/NCBI/proteins/$o ];then
		echo $i;
		gbk2genes.pl -in $i -p -s > data/NCBI/proteins/$o;
	fi

done

##make files of single-copy proteins based on self-blast.
ident="90";
qcov="70";

function singleator {
	i=$1
	echo $i;
	ident=$2;
	qcov=$3;
	o=$(basename $i faa)single-id$ident-len$qcov.faa;
	t=$(basename $i faa)bls;
	s=$(basename $i faa)singles-id$ident-len$qcov.txt;
	tmp=$(basename $i faa)tmp;

	if [ ! -f data/proteins/single-id$ident-len$qcov/$o ];then

		if [ ! -f data/proteins/self-bls/$t ];then
			makeblastdb -in $i -dbtype prot \
			-logfile data/proteins/self-bls/blastdb.log;
			echo "running blastp."
			blastp -query $i -db $i \
			-qcov_hsp_perc $qcov -num_threads 4 \
			-outfmt '6 qseqid sseqid pident qlen slen length' \
			-out data/proteins/self-bls/$t;

			rm $i.p*;

		fi

	#count blast hits for each query
	./src/bls.hit-counter.pl -in data/proteins/self-bls/$t \
	-pident $ident -minlen $qcov -only 1 \
	-out data/proteins/self-bls/$s;

	#new fasta of only single-copy genes
	cut -f1 data/proteins/self-bls/$s > data/proteins/self-bls/$tmp;
	FastA.filter.pl data/proteins/self-bls/$tmp $i > data/proteins/single-id$ident-len$qcov/$o;

	rm data/proteins/self-bls/$tmp;

	fi

}
export -f singleator

fasta=( data/proteins/*faa );
parallel -j 4 -k singleator {} $ident $qcov ::: ${fasta[@]};


##pairwise ortholog matching (RBM)
singles=( data/proteins/single-id$ident-len$qcov/*faa );

today=$(date +%Y%m%d);
if [ ! -d results/RBM/$today ];then
	mkdir results/RBM/$today;
fi

function runRBM {
	outname=$(basename $1 | grep -Po "^[\w\.]+single")_v_$(basename $2 | grep -Po "^[\w\.]+single");
	q=$( echo print $3/100 | perl);
	~/enveomics/Scripts/rbm.rb -q -i 40 -f $q -1 $1 -2 $2 > results/RBM/$4/$outname.rbm ;
}
export -f runRBM

for ((i=0; i<${#singles[@]}-1; i++));
do
        echo -e "\t" ${singles[$i]} "("$(($i + 1)) "of" $((${#singles[@]} - 1))")"...;
        #parallel -k echo ${singles[$i]} {} ::: ${singles[@]:$i+1};
				parallel -j 8 -k runRBM ${singles[$i]} {} $qcov $today ::: ${singles[@]:$i+1};
done
