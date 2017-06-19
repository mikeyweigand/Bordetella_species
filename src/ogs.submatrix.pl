#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;


&GetOptions(	'in=s' => \my$inogs,		#
		'list=s' => \my$list,
		'ids=s' => \my$prefix,
		'min=s' => \my$min,
		'ind' => \my$ind,
		'max=s' => \my$max,
		'out=s' => \my$out);		#
($inogs and $out) or &HELP_MESSAGE;
($list or $prefix or ($min and $max)) or &HELP_MESSAGE;
#if($list and ($min or $max)){ &HELP_MESSAGE };
if(($min) and (!$max)){ &HELP_MESSAGE};
if(($max) and (!$min)){ &HELP_MESSAGE};

if($prefix){
	print "\nWARNING: parameter '-ids' overrides all other options. A column of og identifiers have been added but no filtering was performed.\n\n";
}

my@header=();
my%genomes=();
my$j=0;

if($list){
	open LIST, "$list";
	while(my$l = <LIST>){
		chomp$l;
		$genomes{ $l } = 1;
	}
	close LIST;
}

open OGS, "$inogs";
open OUT, ">$out";
while(my$i = <OGS>){
	chomp$i;
	my@si=split("\t",$i);
	#first line
	if( scalar(@header) < 1){
		if($prefix){
			@header = @si;
			print OUT "\t". join("\t",@si)."\n";
		}elsif($list){
			shift@si;
			for(my$c = 0; $c < @si; $c++){
				if(exists( $genomes{ $si[$c] })){
					push(@header,$c);
					delete $genomes{ $si[$c]};
				}
			}
			if($ind){ print OUT "\t" };
			if(($min or $min == 0) and ($max or $max ==0)){
				print OUT join("\t",@si)."\n";
			}else{
				print OUT join("\t",@si[@header])."\n";
			}
			print "Genomes matched: ".scalar(@header)."\n";
			print "Genomes missing: ".scalar(keys%genomes)."\n";

		}else{
			@header = @si;
			print OUT join("\t",@si)."\n";
		}

		if($list and scalar(keys%genomes) > 0){
			print "\t\t\t".join("\n\t\t\t",(sort(keys%genomes)))."\n";
		}
	#all subsequent lines
	}else{
		if($prefix){
			my$p = &countogs( \@si );
			my$n=sprintf("%04d", $j);
			if($p > 1){
				print OUT join("_",($prefix,$n,$p))."\t";
				print OUT join("\t",@si)."\n";
			}
			$j++;
			next;

		}elsif($list){
			my$pan1='';
			if($ind){ $pan1 = shift@si };
			if(($min or $min == 0) and ($max or $max ==0)){
				my@check = @si[@header];
				my$num = &countogs( \@check );
				if(($num >= $min) && ($num <= $max)){
					if($ind){ print OUT $pan1."\t" };
					print OUT join("\t",@si)."\n";
				}

			}else{
				if($ind){ print OUT $pan1."\t" };
				print OUT join("\t",@si[@header])."\n";
			}

		}elsif(($min or $min == 0) and ($max or $max ==0)){
			my$pan='';
			if($ind){ $pan = shift@si };
			my$present = &countogs( \@si );
			if(($present >= $min) && ($present <= $max)){
				if($ind){ print OUT $pan."\t" };
				print OUT join("\t",@si)."\n";
			}
		}
	}

}
close OUT;


################################
sub countogs {
	my@in = @{$_[0]};
	my$zero = 0;
	foreach my$i (@in){
		if($i =~ /^\d+$/){
			if($i == 0){ $zero++};
		}elsif($i eq '-'){
			$zero++;
		}
	}
	return( scalar@in - $zero );
}


sub HELP_MESSAGE { die "
.Description:
   Takes ogs matrix and makes new sub-matrix based on: [1] list of genomes, or [2] frequency min and/or max.
   If -list and -min/max defined, then new sub-matrix filters ogs by frequency in list of genomes.

.Usage: $0 -in [ogs] -out [out] -list [list] -min [int] -max [int]

   [mandatory]
	 -in	<ogs>	Tab-separated matrix of ogs, like from ogs.mcl.rb.
	 -out	<out>	New sub-matrix of ogs in list.

   [optional]
	 -list	<list>	List of genomes/names to match in header of ogs matrix.
	 -r		Reverse list, returns sub-matrix of genomes NOT in list. (NOT FUNCTIONAL)
	 -min	<int>	Only output ogs present in at _least_ this many genomes.
	 -max	<int>	Only output ogs present in at _most_ this many genomes.
	 -ind		Input ogs matrix has indices in first column (Default: 'off').
	 -ids	<pre>	Add first column with og identifiers using given prefix, OVERRIDES ALL OTHER OPTIONS (eg 'PAN').


   [dependencies]

" }
