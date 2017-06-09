#!/bin/bash

# function runCLEAN {
# 	o=$(basename $1).clean;
# 	if [ ! -f results/mauve/Bp-collinear-clean/$o ];then
# 		../Pertussis_n257/Current/mauve.backbone-clean.pl -in $1 \
# 		-genomes data/Bp-list4mauve.txt \
# 		-query ../Pertussis/IS-all.fasta \
# 		-out results/mauve/Bp-collinear-clean/$o;
# 	fi
# }
# export -f runCLEAN
#
# bbone=(results/mauve/Bp-collinear/*backbone );
# parallel -j 12 -k runCLEAN {} ::: ${bbone[@]};
#parallel -j 8 -k echo {} ::: ${bbone[@]};

function runCHECK {
	o=$(basename $1).sum;
	c=$(basename $1).table;
	if [ ! -f results/mauve/Bp-collinear-20170708/gap1500-summary/$o ];then
		../EPS_genome_epi/src/mauve.collinear-check.pl -in $1 -gap 1500 \
		-out results/mauve/Bp-collinear-20170708/gap1500-summary/$o \
		-coords results/mauve/Bp-collinear-20170708/gap1500-table/$c;
	fi
}
export -f runCHECK



bbclean=(results/mauve/Bp-collinear-clean/*backbone.clean );
parallel -j 12 -k runCHECK {} ::: ${bbclean[@]};



#awk -v OFS='\t' '$3="0"' 20170605-cat.sum > ./20170605-cat-zerogap.sum

#for i in $(ls ./results/mauve/Bp-collinear-clean-NR-zerogap/*clean); do o=$(basename $i).inverts; echo $i; ./src/mauve.backbone-inverts.pl -in $i -out results/mauve/Bp-collinear-clean-NR-zerogap-inverts/$o; done
#grep -c "^-" ./Bp-collinear-clean-NR-zerogap-inverts/* | grep :1$ > singles-20170606-zerogap.txt
#for i in $(cat singles-20170606-zerogap.txt); do n=$(basename $i .xmfa.backbone.clean.inverts:1); j=$(echo $n | sed 's/-/,/g'); echo $j; done > singles-zerogap-list-20170606.csv


#> zerogap = read.csv("singles-zerogap-list-20170606.csv",header=F)
#> zerogap.net = network(zerogap,directed=F,matrix.type='edgelist')
#> ggnet2(zerogap.net, label=T,label.size=3.5,size="degree",layout.par = list(cell.jitter = 0.75))
