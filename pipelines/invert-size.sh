#!/bin/bash

function invertAll {
	g0=$1;
	g1=$2;

	bb=$(find ./results/mauve/20180221-Bp/03.colinear-coords-gap1500/$g0\-$g1\.xmfa.backbone.clean.table);
	f1=$(find ./data/genomes/Bpertussis-fasta-IScollapsed-20180221/$g1\.fasta);
	out=$g0\-$g1\.invertALL.txt;

	g1len=$(FastA.length.pl $f1 | cut -f2);
	g1ter=$(echo "aattcgcataatgtatattatgtaaagt" | blastn -outfmt '6 sstart' -subject $f1);

	# echo -e "\t"$bb"\n\t"$f1"\n\t"$out;
	# echo -e "\t"$g1len"\t"$g1ter"\n";

	echo -ne $g0"\t"$g1"\t" > results/mauve/20180221-Bp/05.colinear-invertAll-gap1500/$out;

	./src/mauve.invert-symmetric.pl -in $bb -t $g1ter -size -seq1 $g1len >> results/mauve/20180221-Bp/05.colinear-invertAll-gap1500/$out;



}
export -f invertAll


while read line
do
	g0=$(echo $line | cut -d " " -f1);
	g1=$(echo $line | cut -d " " -f2);
	echo -e $g0"\t"$g1;

	invertAll $g0 $g1;

done < "${1:-/dev/stdin}"
