#!/bin/bash
##prepare protein fasta files with simplified headers
echo -e "\nExtracting protein sequences...";

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

echo -e "\nExtracting gene sequences...";

for i in $(ls data/genomes/genbank/*gbk);
do
	o=$(basename $i gbk)fna;
	if [ ! -f data/genes/$o ];then
		echo $i;
		gbk2genes.pl -in $i -s > data/genes/$o;
	fi
done


##make files of single-copy proteins based on self-blast.
ident="90";
qcov="70";
echo -e "\nPreparing single-copy fasta files...";
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
echo -e "\nRunning pairwise RBM..."
singles=( data/proteins/single-id$ident-len$qcov/*faa );

today=$(date +%Y%m%d);
#today='20170512';
if [ ! -d results/RBM/$today ];then
	mkdir results/RBM/$today;
fi


function runRBM {
	outname=$(basename $1 | grep -Po "^[\w\.]+single")_v_$(basename $2 | grep -Po "^[\w\.]+single");
	q=$( echo print $3/100 | perl);
	if [ ! -f results/RBM/$4/$outname.rbm ];then
		~/enveomics/Scripts/rbm.rb -t 2 -q -i 40 -f $q -1 $1 -2 $2 > results/RBM/$4/$outname.rbm ;
	fi
}
export -f runRBM

for ((i=0; i<${#singles[@]}-1; i++));
do
        echo -e "\t" ${singles[$i]} "("$(($i + 1)) "of" $((${#singles[@]} - 1))")"...;
        #parallel -k echo ${singles[$i]} {} ::: ${singles[@]:$i+1};
				parallel -j 8 -k runRBM ${singles[$i]} {} $qcov $today ::: ${singles[@]:$i+1};
done

##cluster matched orthologs via MCL
echo -e "\nClustering RBMs into orthologs groups...\n";
if [ ! -d results/RBM/$today-ogs ];then
	mkdir results/RBM/$today-ogs;
fi

~/enveomics/Scripts/ogs.mcl.rb -t 12 -d results/RBM/$today/ \
-f '(\S+)\.single_v_(\S+)\.single\.rbm' -o results/RBM/$today-ogs/$today-ogs.txt;

##estimate descriptive stats
~/enveomics/Scripts/ogs.stats.rb -o results/RBM/$today-ogs/$today-ogs.txt \
-t results/RBM/$today-ogs/$today-stats.txt;

##extract core genome/proteome
echo -e "\nExtracting og sequences for alignment...\n";
if [ ! -d results/RBM/$today-core ];then
	mkdir results/RBM/$today-core;
fi
ulimit -n 5120;
~/enveomics/Scripts/ogs.extract.rb -i results/RBM/$today-ogs/$today-ogs.txt \
-o results/RBM/$today-core -c 1 -d 1 -p \
-s data/proteins/%s.faa;



##align ogs, concatenate, and calc phylogeny
echo -e "\nAligning core og sequences with muscle...\n";
if [ ! -d results/RBM/$today-aln ];then
	mkdir results/RBM/$today-aln;
fi

function runMUSCLE {
	aln=$( echo $1 | sed 's/-core/-aln/g' | sed 's/fa$/aln/g');
	#echo -e "\t"$aln;
	muscle3.8.31 -in $1 -out $aln -quiet;
}
export -f runMUSCLE

coreogs=( results/RBM/$today-core/*.fa );
parallel -j 12 -k runMUSCLE {} ::: ${coreogs[@]};

#concatenate
~/enveomics/Scripts/Aln.cat.rb -i '-' results/RBM/$today-aln/*aln > results/RBM/$today-ogs/$today-core.aln;
~/enveomics/Scripts/Aln.cat.rb -i '-' -I results/RBM/$today-aln/*aln > results/RBM/$today-ogs/$today-core.variable.aln;


#ML phylogeny
raxmlHPC-PTHREADS -f a -m PROTGAMMAGTR -x 12345 -p 12345 -T 12 -N 100 --silent \
-n $today-core.variable\
-s results/RBM/$today-ogs/$today-core.variable.aln \
-w ~/Documents/Bordetella_species/results/RBM/$today-ogs/

BLOCK
#bord = read.table("20170522-ogs-bin.var.txt", sep="\t",row.names=1,header=T)
#bordx = as.matrix(bord)
#bordcc=seq(1,ncol(bordx),1)
#heatmap.2(bordx, dendrogram="col",srtCol=45,  adjCol = c(1,1),trace="none",colsep=bordcc,cexCol=0.85,col=c("#F7F7F7","#000000","#CA0020"),breaks=c(-1,-0.5,0.5,1))

#ccolors = c("#984ea3","#337eb8","#984ea3","#ff7f00","#4daf4a","#ff7f00","#f5f5f5","#984ea3","#f5f5f5","#337eb8","#984ea3","#ff7f00","#337eb8","#337eb8","#ff7f00","#337eb8","#4daf4a","#984ea3","#984ea3","#e41a1c","#f5f5f5","#337eb8","#ff7f00","#984ea3","#ff7f00","#ff7f00","#e41a1c","#f5f5f5","#984ea3","#4daf4a","#e41a1c","#984ea3","#4daf4a","#984ea3","#4daf4a","#f5f5f5","#984ea3","#984ea3","#e41a1c","#337eb8","#4daf4a","#337eb8","#ff7f00","#4daf4a","#984ea3","#f5f5f5","#4daf4a","#4daf4a","#f5f5f5","#e41a1c","#337eb8","#984ea3")

#pdf(file="20170522-ogs-bin.var01.pdf", width=8, height=8, pointsize=10)
#heatmap.2(bordx, dendrogram="col",srtCol=45,  adjCol = c(1,1),trace="none",colsep=bordcc,cexCol=0.85,col=c("#F7F7F7","#000000","#CA0020"),breaks=c(-1,-0.5,0.5,1), ColSideColors=ccolors)
#dev.off()
