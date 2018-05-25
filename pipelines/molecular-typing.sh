#!/bin/bash

if [ ! -d "results/molecular-typing" ]; then
  mkdir results/molecular-typing;
fi

now=$(date +%Y%m%d);

echo -n > results/molecular-typing/$now-fimH.txt
echo -n > results/molecular-typing/$now-prn.txt
echo -n > results/molecular-typing/$now-ptxP.txt
echo -n > results/molecular-typing/$now-ptxA.txt
echo -n > results/molecular-typing/$now-ptxB.txt

for genome in $(ls data/genomes/Bpertussis-fasta/*fasta);
do

	name=$(basename $genome .fasta);
	pdl=${name:3:4};
	echo -n $pdl;

	echo -ne "\tprn";
	/home/yrh8/Documents/EPS_genome_epi/src/typing.prn.pl \
	-in $genome \
	-wt /home/yrh8/Documents/EPS_genome_epi/data/allele-references/prn-alleles-PDL.fasta \
	-ref /home/yrh8/Documents/EPS_genome_epi/data/allele-references/prn-mut-20160112.fasta \
	-name $pdl \
	-dir results/molecular-typing/blastout-prn \
	>> results/molecular-typing/$now"-prn.txt";

	echo -ne "\tfimH";
	/home/yrh8/Documents/EPS_genome_epi/src/typing.fimH.pl \
	-in $genome \
	-ref /home/yrh8/Documents/EPS_genome_epi/data/allele-references/fimH-alleles-20160314.fasta \
	-name $pdl \
	-dir results/molecular-typing/blastout-fimH \
	>> results/molecular-typing/$now"-fimH.txt";

	echo -ne "\tptxP";
	/home/yrh8/Documents/EPS_genome_epi/src/typing.ptxP.pl \
	-in $genome \
	-ref /home/yrh8/Documents/EPS_genome_epi/data/allele-references/ptxP-alleles-20160705.fasta \
	-name $pdl \
	-dir results/molecular-typing/blastout-ptxP \
	>> results/molecular-typing/$now"-ptxP.txt";

	echo -ne "\tptxA";
	/home/yrh8/Documents/EPS_genome_epi/src/typing.ptxA.pl \
	-in $genome \
	-ref /home/yrh8/Documents/EPS_genome_epi/data/allele-references/ptxA-alleles-20160316.fasta \
	-name $pdl \
	-dir results/molecular-typing/blastout-ptxA \
	>> results/molecular-typing/$now"-ptxA.txt";

	echo -ne "\tptxB";
	/home/yrh8/Documents/EPS_genome_epi/src/typing.ptxB.pl \
	-in $genome \
	-ref /home/yrh8/Documents/EPS_genome_epi/data/allele-references/ptxB-alleles-20160316.fasta \
	-name $pdl \
	-dir results/molecular-typing/blastout-ptxB \
	>> results/molecular-typing/$now"-ptxB.txt";

	echo

done

#combine typing results into table
Table.merge.pl -s -e UNK -h CDC_ID \
results/molecular-typing/$now"-prn.txt" \
results/molecular-typing/$now"-fimH.txt" \
results/molecular-typing/$now"-ptxP.txt" \
results/molecular-typing/$now"-ptxA.txt" \
results/molecular-typing/$now"-ptxB.txt" \
| sed 's/'$now'-//g' > results/molecular-typing/$now"-combined.txt";

head -1 results/molecular-typing/$now"-combined.txt" > results/molecular-typing/$now"-combined.sorted.txt";

tail -n+2 results/molecular-typing/$now"-combined.txt" | sort >> results/molecular-typing/$now"-combined.sorted.txt";
