#!/bin/bash

function runCLEAN {
	o=$(basename $1).clean;
	if [ ! -f results/mauve/20180221-Bp/02.colinear-clean-20180425/$o ];then
		~/Documents/Pertussis_n257/Current/mauve.backbone-clean.pl -in $1 \
		-genomes data/list-Bpertussis-IScollapsed-20180221.list \
		-query data/IS-rRNA-MFS.fasta \
		-out results/mauve/20180221-Bp/02.colinear-clean-20180425/$o;
	fi
}
export -f runCLEAN
#
function wtfECHO {
	echo $1;
}
export -f wtfECHO

#bbone=( find results/mauve/20180221-Bp/01.colinear-alignments/*backbone );
#parallel -j 12 -k runCLEAN {} ::: ${bbone[@]};
#parallel -j 8 -k wtfECHO {} ::: ${bbone[@]};
#find results/mauve/20180221-Bp/01.colinear-alignments/ -name "*.backbone" | sort | parallel -j 8 -k wtfECHO {};

find results/mauve/20180221-Bp/01.colinear-alignments/ -name "*.backbone" | sort | parallel -j 12 -k runCLEAN {};

function runCHECK {
	o=$(basename $1).sum;
	c=$(basename $1).table;
	if [ ! -f results/mauve/20180221-Bp/03.colinear-check-gap1500-20180425/$o ];then
		~/Documents/EPS_genome_epi/src/mauve.collinear-check.pl -in $1 -gap 1500 \
		-out results/mauve/20180221-Bp/03.colinear-check-gap1500-20180425/$o \
		-coords results/mauve/20180221-Bp/03.colinear-coords-gap1500-20180425/$c;
	fi
}
export -f runCHECK
#
# #bbclean=(results/mauve/20170619-Bp-collinear/Backbone-clean/*backbone.clean );
# #parallel -j 12 -k runCHECK {} ::: ${bbclean[@]};
# #find results/mauve/20180221-Bp/02.colinear-clean/ -name "*.backbone.clean" | sort | parallel -j 8 -k wtfECHO {};

find results/mauve/20180221-Bp/02.colinear-clean-20180425/ -name "*.backbone.clean" | sort | parallel -j 12 -k runCHECK {};



# function runMATCH {
# 	g0=$(echo $1 | cut -d "," -f1 | sed 's/Bp_//g');
# 	g1=$(echo $1 | cut -d "," -f2 | sed 's/Bp_//g');
# 	bb=Bp_$g0\_Illumina_final_IScollapsed-Bp_$g1\_Illumina_final_IScollapsed.xmfa.backbone;
# 	out=$g0\_$g1\-rRNAmatches.txt;
#
# 	./src/mauve.pairwiseISmatcher.pl -name0 $g0 -name1 $g1 \
# 	-bbone results/mauve/20170619-Bp-collinear/$bb \
# 	-out results/mauve/20170619-Bp-collinear/ISmatch-singles-20170626/$out \
# 	-seq0 data/genomes/Bpertussis-fasta-IScollapsed/Bp_$g0\_Illumina_final_IScollapsed.fasta \
# 	-seq1 data/genomes/Bpertussis-fasta-IScollapsed/Bp_$g1\_Illumina_final_IScollapsed.fasta \
# 	-query ~/Documents/Pertussis/E976-rRNA.fasta;
#
# }
# export -f runMATCH
#
# singles=$(head -40 results/mauve/20170619-Bp-collinear/Collinear-mcl/20170620-collinear-gap1500.NR.sum.1invert-1gap-sym.csv);
# parallel -j 12 -k runMATCH {} ::: ${singles[@]};

# function runMATCH {
# 	g0=$(echo $1 | cut -d "," -f1 | sed 's/Bp_//g');
# 	g1=$(echo $1 | cut -d "," -f2 | sed 's/Bp_//g');
# 	bb=Bp_$g0\_Illumina_final_IScollapsed-Bp_$g1\_Illumina_final_IScollapsed.xmfa.backbone.clean.table;
# 	out=$g0\_$g1\-rRNAmatches-boundary.txt;
#
# 	./src/mauve.boundary-matcher.pl -name0 $g0 -name1 $g1 \
# 	-bbone results/mauve/20170619-Bp-collinear/Collinear-table-gap1500/$bb \
# 	-out results/mauve/20170619-Bp-collinear/ISmatch-singles-20170626/$out \
# 	-seq0 data/genomes/Bpertussis-fasta-IScollapsed/Bp_$g0\_Illumina_final_IScollapsed.fasta \
# 	-seq1 data/genomes/Bpertussis-fasta-IScollapsed/Bp_$g1\_Illumina_final_IScollapsed.fasta \
# 	-query ~/Documents/Pertussis/E976-rRNA.fasta -flank 5500;
# 	#-query data/genomes/IS-all.fasta -boundary;
#
# }
# export -f runMATCH
#
# singles=$(head -40 results/mauve/20170619-Bp-collinear/Collinear-mcl/20170620-collinear-gap1500.NR.sum.1invert-1gap-sym.csv);
# parallel -j 12 -k runMATCH {} ::: ${singles[@]};

# function invertSize {
# 	g0=$(echo $1 | cut -d "," -f1 | sed 's/Bp_//g');
# 	g1=$(echo $1 | cut -d "," -f2 | sed 's/Bp_//g');
# 	bb=Bp_$g0\_Illumina_final_IScollapsed-Bp_$g1\_Illumina_final_IScollapsed.xmfa.backbone.clean.table;
# 	out=$g0\_$g1\-invertsize.txt;
# 	g1len=$(FastA.length.pl ~/Documents/Bordetella_species/data/genomes/Bpertussis-fasta-IScollapsed/Bp_$g1\_Illumina_final_IScollapsed.fasta | cut -f2);
# 	g1ter=$(echo "aattcgcataatgtatattatgtaaagt" | blastn -outfmt '6 sstart' -subject ~/Documents/Bordetella_species/data/genomes/Bpertussis-fasta-IScollapsed/Bp_$g1\_Illumina_final_IScollapsed.fasta);
# 	./src/mauve.invert-symmetric.pl -in results/mauve/20170619-Bp-collinear/Collinear-table-gap1500/$bb \
# 	-t $g1ter -size -seq1 $g1len > results/mauve/20170619-Bp-collinear/Collinear-summary-gap1500-NR-singles-symmetry2/$out;
#
#
# }
# export -f invertSize
#
# singles=$(head -40 results/mauve/20170619-Bp-collinear/Collinear-mcl/20170620-collinear-gap1500.NR.sum.1invert-1gap-sym.csv);
# parallel -j 12 -k invertSize {} ::: ${singles[@]};


# function invertAll {
# 	#g0=$(echo $1 | cut -d "," -f1 | sed 's/Bp_//g');
# 	#g1=$(echo $1 | cut -d "," -f2 | sed 's/Bp_//g');
# 	o=$(basename $1 .sum);
# 	g0=${o:3:4};
# 	g1=${o:38:4};
# 	g0single=$(grep "$g0" results/mauve/20170619-Bp-collinear/Collinear-mcl/20170620-collinear-gap1500.NR.sum.1invert-1gap-sym.csv);
# 	g1single=$(grep "$g1" results/mauve/20170619-Bp-collinear/Collinear-mcl/20170620-collinear-gap1500.NR.sum.1invert-1gap-sym.csv);
# 	if [ -z "$g0single" ] || [ -z "$g1single" ]; then  #unless both are part of single invert/del network
# 		echo $o $g0 $g1;
# 		bb=Bp_$g0\_Illumina_final_IScollapsed-Bp_$g1\_Illumina_final_IScollapsed.xmfa.backbone.clean.table;
# 		out=$g0\_$g1\-invertALL.txt;
# 		g1len=$(FastA.length.pl ~/Documents/Bordetella_species/data/genomes/Bpertussis-fasta-IScollapsed/Bp_$g1\_Illumina_final_IScollapsed.fasta | cut -f2);
# 		g1ter=$(echo "aattcgcataatgtatattatgtaaagt" | blastn -outfmt '6 sstart' -subject ~/Documents/Bordetella_species/data/genomes/Bpertussis-fasta-IScollapsed/Bp_$g1\_Illumina_final_IScollapsed.fasta);
#
# 	 ./src/mauve.invert-symmetric.pl -in results/mauve/20170619-Bp-collinear/Collinear-table-gap1500/$bb \
# 	 -t $g1ter -all -seq1 $g1len > results/mauve/20170619-Bp-collinear/Collinear-table-gap1500-NR-all/$out;
#
# 	fi
#
# }
# export -f invertAll
#
# singles=(results/mauve/20170619-Bp-collinear/Collinear-summary-gap1500-NR/*backbone.clean.sum);
# parallel -j 12 -k invertAll {} ::: ${singles[@]};


# cat results/mauve/20170619-Bp-collinear/Collinear-summary-gap1500-NR/*.sum | cut --output-delimiter="," -f1,2,3,4 | grep ",0,1$" | sed 's/_Illumina_final_IScollapsed//g' | sed 's/Bp_//g'



#> invertgaps0620=read.csv("./20170620-collinear-gap1500.NR.sum.1invert-1gap.csv", header=F)
#> invertgaps0620.net = network(invertgaps0620, directed=F, matrix.type='edgelist', ignore.eval=F, names.eval="weights")
#> set.edge.attribute(invertgaps0620.net, "lty", ifelse(invertgaps0620.net %e% "weights" > 1, 2, 1))
#> ggnet2(invertgaps0620.net, label=T, size="degree", label.size=3, edge.lty="lty")

#> symm0620=read.csv("./20170620-collinear-gap1500.NR.sum.1invert-1gap-sym.csv", header=F)
#> symm0620.net = network(symm0620, directed=F, matrix.type='edgelist', ignore.eval=F, names.eval=c("weights","symm"))
#> set.edge.attribute(symm0620.net, "lty", ifelse(symm0620.net %e% "weights" > 1, 2, 1))
#> set.edge.attribute(symm0620.net, "color", ifelse(symm0620.net %e% "symm" > 1, "red", "black"))
#> ggnet2(symm0620.net, label=T, size="degree", label.size=2, edge.lty="lty", edge.color = "color", layout.par = list(cell.jitter = 0.5))

#> ggplot(invsize, aes(V2,V1)) + geom_boxplot() + geom_jitter(width=0.2)

#awk -v OFS='\t' '$3="0"' 20170605-cat.sum > ./20170605-cat-zerogap.sum



#> zerogap = read.csv("singles-zerogap-list-20170606.csv",header=F)
#> zerogap.net = network(zerogap,directed=F,matrix.type='edgelist')
#> ggnet2(zerogap.net, label=T,label.size=3.5,size="degree",layout.par = list(cell.jitter = 0.75))
