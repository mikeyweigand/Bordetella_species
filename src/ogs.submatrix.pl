#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;


&GetOptions(	'in=s' => \my$inogs,		#
		'list=s' => \my$list,
		'out=s' => \my$out);		#
($inogs and $list and $out) or &HELP_MESSAGE;

my@header=();
my%genomes=();
open LIST, "$list";
while(my$l = <LIST>){
	chomp$l;
	$genomes{ $l } = 1;
}
close LIST;

open OGS, "$inogs";
open OUT, ">$out";
while(my$i = <OGS>){
	chomp$i;
	my@si=split("\t",$i);
	if( scalar(@header) < 1){
		for(my$c = 0; $c < @si; $c++){
			if(exists( $genomes{ $si[$c] })){
				push(@header,$c);
				delete $genomes{ $si[$c]};
			}
		}
		print OUT join("\t",@si[@header])."\n";
		print "Genomes matched: ".scalar(@header)."\n";
		print "Genomes missing: ".scalar(keys%genomes)."\n";
		if(scalar(keys%genomes) > 0){
			print "\t\t\t".join("\n\t\t\t",(sort(keys%genomes)))."\n";
		}
	}else{
		print OUT join("\t",@si[@header])."\n";
	}
}
close OUT;


################################
sub HELP_MESSAGE { die "
.Description:
   Takes ogs matrix and makes new sub-matrix from list of genomes.

.Usage: $0 -in [ogs] -out [out] -list [list]

   [mandatory]
	 -in	<ogs>	Tab-separated matrix of ogs, like from ogs.mcl.rb.
	 -out	<out>	New sub-matrix of ogs in list.
	 -list	<list>	List of genomes/names to match in header of ogs matrix.

   [optional]
	 -r	Reverse list, returns sub-matrix of genomes NOT in list.


   [dependencies]

" }
