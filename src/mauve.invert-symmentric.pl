#!/usr/bin/perl -w
use strict;
use Getopt::Long;


&GetOptions(	'in=s' => \my$infile,		#
		't=s' => \my$term);		#
($infile and $term) or &HELP_MESSAGE;

open IN, "$infile";
while(my$i = <IN>){
	unless($i =~ /seq/){
		chomp$i;
		my@si=split("\t",$i);
		if(($si[0] < 1) && ($si[1] < 1)){
			if(($si[2] < $term) && ($si[3] > $term)){

				print $i."\tSymmetric\n";
			}else{
				print $i."\tAsymmetric\n";
			}
		}
	}
}


sub HELP_MESSAGE { die "
.Description:
   Characterize single inversions in pairwise mauve backbone files as 'symmetric' or 'asymmetric'.

.Usage: $0 -in [in.txt] -t [int] > out.txt

   [mandatory]
	 -in <in.txt>	Mauve backbone file, probably the simplified coordinate table from from 'mauve.collinear-check.pl'.
	 -t <int>	Coordinate of replication terminus.

   [optional]


   [dependencies]


" }
