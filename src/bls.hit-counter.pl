#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Statistics::R;

#set default cut-offs:
my$minident = 90;
my$minlen = 85;

&GetOptions(	'in=s' => \my$inbls,		#
		'pident=s' => $minident,
		'minlen=s' => $minlen,
		'only=s' => \my$only,
		'index=s' => \my$indices,
		'out=s' => \my$out);		#
($inbls and $out) or &HELP_MESSAGE;

my%counts=();

open BLS, "$inbls";
open OUT, ">$out";
while( my$i = <BLS>){
	chomp$i;
	my@s = split("\t",$i);

	if( ($s[2] >= $minident) && (($s[5]/$s[3])*100 >= $minlen) ){
		if( exists( $counts{$s[0]} )){
			$counts{$s[0]}++;
		}else{
			$counts{$s[0]} = 1;
		}
	}
}

foreach my$k (sort keys%counts){
	if($only){
		if($counts{$k} == $only){
			print OUT $k ."\t". $counts{$k}."\n";
		}
	}else{
		print OUT $k ."\t". $counts{$k}."\n";
	}
}

################################
sub HELP_MESSAGE { die "
.Description:
   Counts the number of hits for each query in a tab-separated blast output file.

.Usage: $0 -in [in.bls] -out [out.txt]

   [mandatory]
	 -in	<in.bls>	Tab separated blast output file.
	 			(expects: outfmt '6 qseqid sseqid pident qlen slen length')
	 -out	<out.txt>	Output table of query match counts.

   [optional]
	 -pident	<Real, 0-100>	Minimum percent identity for matches (Default = 90.0).
	 -minlen	<Real, 0-100>	Minimum query coverage for matches (Default = 85.0).
	 -only		<integer>	Only output queries with exactly this number of matches.
	 -index		<list>		[not available] Column indices for qseqid, sseqid, pident, qlen, and length if different from expected (Eg. '1,2,3,4,6').
	 -noself			[not available] Ignore self matches.


   [dependencies]

" }
