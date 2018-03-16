#!/usr/bin/perl -w
use strict;
use Getopt::Long;


&GetOptions(	'in=s' => \my$infile,		#
		'g=s' => \my$gaps,		#
		'list=s' => \my$list,
		'i=s' => \my$inverts);		#
($infile) or &HELP_MESSAGE;
(($gaps or $gaps == 0) and ($inverts or $inverts ==0)) or &HELP_MESSAGE;

my%genomes = ();
if($list){
	open LIST, "$list";
	while(my$l = <LIST>){
		chomp$l;
		$genomes{ $l } = 1;
	}
}


open IN, "$infile";
while(my$i = <IN>){
	chomp$i;
	my@si=split("\t",$i);

	if(($si[2] == $gaps) && ($si[3] == $inverts)){

		if($list){
			if( exists($genomes{$si[0]}) && exists($genomes{$si[1]}) ){
	 			print $i."\n";
			}

		}else{
			print $i."\n";
		}
	}
}



sub HELP_MESSAGE { die "
.Description:
   Filters concatenated output of 'mauve.collinear-check.pl' based on defined gap and/or inversion counts.

.Usage: $0 -in [in.txt] -g [int] -i [int] > out.txt

   [mandatory]
	 -in <in.txt>	Output from 'mauve.collinear-check.pl'.
	 			(Assumes format: GenomeA GenomeB [gaps] [inversions])
	 -g <int>	Number of desired gaps.
	 -i <int>	Number of desired inversions.

   [optional]
	 -list <list.txt>	Only consider genomes from this list.

   [dependencies]


" }
