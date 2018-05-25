#!/usr/bin/perl -w
use strict;
use Getopt::Long;


&GetOptions(	'inverts=s' => \my$infile,		#
		'sizes=s' => \my$sizes,
		'gaps=s' => \my$gaps,
		'missing=s' => \my$missing,
		'simple' => \my$simple,
		'rename=s' => \my$rename);		#
($infile) or &HELP_MESSAGE;

if($missing){
	&MISSING( $missing, $infile, $gaps );
	exit;
}

#store list or new names
my%names = ();
if($rename){
	open NAMES, "$rename";
	while(my$re = <NAMES>){
		chomp$re;
		my@sre = split("\t",$re);
		$names{ $sre[0] } = $sre[1];
	}
	close NAMES;
}

#process invert table and prepare csv
open INVERTS, "$infile";
while(my$i = <INVERTS>){
	chomp$i;
	my@sinv = split("\t",$i);

	if($rename){
		$sinv[0] = $names{ $sinv[0] };
		$sinv[1] = $names{ $sinv[1] };
	}

	print join(",",@sinv[0,1]);
	if($i =~ m/Symm/){
		print ",1,1";
	}elsif($i =~ m/Asymm/){
		print ",1,2";
	}
	print "\n";
}
close INVERTS;

#if gaps given, add those to output
if($gaps){
	open GAPS, "$gaps";
	while(my$g = <GAPS>){
		chomp$g;
		my@sg = split("\t",$g);

		if($rename){
			$sg[0] = $names{ $sg[0] };
			$sg[1] = $names{ $sg[1] };
		}

		print join(",",@sg[0,1]).",2,1\n";
	}
	close GAPS;
}

#################
sub MISSING {
	my%names=();
	open INVERTS, "$_[1]";
	while(my$i = <INVERTS>){
		chomp$i;
		my@sinv = split("\t",$i);
		$names{ $sinv[0] } = 1;
		$names{ $sinv[1] } = 1;
	}
	close INVERTS;
	if($_[2]){
		open GAPS, "$_[2]";
		while(my$g = <GAPS>){
			chomp$g;
			my@sg = split("\t",$g);
			$names{ $sg[0] } = 1;
			$names{ $sg[1] } = 1;
		}
		close GAPS;
	}

	open LIST, "$_[0]";
	while(my$l = <LIST>){
		chomp$l;
		my@sl = split("\t",$l);
		unless( exists( $names{ $sl[0] } ) ){
			print $l."\n";
		}
	}
	close LIST;
}

sub HELP_MESSAGE { die "
.Description:
   Prepare comma-separated edge list for network analysis of single inverts and gaps.

.Usage: $0 -inverts [in.txt] > out.csv

   [mandatory]
	 -inverts	<in.txt>	Tab-separated file of pairwise invert info, probably concatenated outputs from 'mauve.invert-symmetric.pl'.

   [optional]
	 -simple			Input table lacks symmetry data, maybe the output from 'mauve.collinear-filter.pl'.
	 -gaps		<gaps.txt>	Tab-separated file of pairwise gap info, probably the output from 'mauve.collinear-filter.pl'.
	 -rename	<names.tsv>	Tab-separated file to rename nodes (expects: 'old' 'new').
	 -sizes		<sizes.tsv>	Tab-separated file of node sizes, probably cluster sizes from 'mauve.collinear-mcl.pl' (WARNING: must match _new_ name if provided with -rename).
	 -missing	<list.txt>	Skip everything and only return node names from list NOT present in the inverts and/or gaps file(s).

   [dependencies]


" }
