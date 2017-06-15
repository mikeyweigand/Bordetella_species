#!/usr/bin/perl -w
use strict;
use Getopt::Long;

#converts og matrix into binary presencen/absence. **optional add first column as PAN/og number**
&GetOptions(	'in=s' => \my$inogs,		#
		'ids=s' => \my$prefix,
		'ind' => \my$ind,
		'out=s' => \my$out);		#
($inogs and $out) or &HELP_MESSAGE;

open IN, "$inogs";
open OUT, ">$out";

my$j=0;
while(my$i = <IN>){
	if($j==0){
		if($prefix){ print OUT "\t"};
		print OUT $i;
		$j=1;
	}else{
		chomp$i;
		my@si=split("\t",$i);
		my@binary=();

		my$p=0;
		my$pan='';
		if($ind){ $pan = shift@si };
		for (my$g=0; $g < @si; $g++){
			if( length($si[$g]) < 2 && $si[$g] =~ /\-/){
				push(@binary,"0");
			}else{
				push(@binary,"1");
				$p++;
			}
		}

		my$n=sprintf("%04d", $j);
		if($p > 1){ # && ($p < @si)){
			if($prefix){ print OUT join("_",($prefix,$n,$p))."\t"};
			if($ind){ print OUT $pan."\t"};
			print OUT join("\t",@binary)."\n";
		}
		$j++;
	}
}

################################
sub HELP_MESSAGE { die "
.Description:
   Converts og matrix into binary presencen/absence matrix.

.Usage: $0 -in [ogs] -out [out] -ids [prefix]

   [mandatory]
	 -in	<ogs>	Tab-separated matrix of ogs, like from ogs.mcl.rb.
	 -out	<out>	New binary matrix.

   [optional]
	 -ids	<pre>	Add first column with og identifiers using given prefix (eg 'PAN').
	 -ind		Input ogs matrix already has indices in first column (Default: 'off').


   [dependencies]

" }
