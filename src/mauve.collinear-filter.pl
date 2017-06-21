#!/usr/bin/perl -w
use strict;
use Getopt::Long;


&GetOptions(	'in=s' => \my$infile,		#
		'g=s' => \my$gaps,		#
		'i=s' => \my$inverts);		#
($infile) or &HELP_MESSAGE;
(($gaps or $gaps == 0) and ($inverts or $inverts ==0)) or &HELP_MESSAGE;

open IN, "$infile";
while(my$i = <IN>){
	chomp$i;
	my@si=split("\t",$i);

	if(($si[2] == $gaps) && ($si[3] == $inverts)){
		print $i."\n";
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


   [dependencies]


" }
