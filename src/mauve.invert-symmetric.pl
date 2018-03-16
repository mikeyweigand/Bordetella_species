#!/usr/bin/perl -w
use strict;
use Getopt::Long;


&GetOptions(	'in=s' => \my$infile,		#
		'size' => \my$size,
		'all' => \my$all,
		'seq1=s' => \my$total,
		't=s' => \my$term);		#
($infile and $term) or &HELP_MESSAGE;
if($size){ $total or &HELP_MESSAGE };
if($all){ $total or &HELP_MESSAGE };
if($all and $size){ &HELP_MESSAGE };


open IN, "$infile";
while(my$i = <IN>){
	unless($i =~ /seq/){
		chomp$i;
		my@si=split("\t",$i);
		if(($si[0] < 1) && ($si[1] < 1)){
			print $i;

			if($all){
				unless($si[3] == $total){
					print "\t-\t-\t";
					my$d1 = abs($si[3] - $term); #distance from term
					my$d2 = $total - $si[3];	#distance from end
					my@dsort = sort{$a <=> $b}($d2,$si[3]);
					if( $dsort[0] < $d1 ){
						print $dsort[0]."\t".$dsort[0]."\tori\n";
					}else{
						print $d1."\t".$d1."\tterm\n";
					}
				}


			}else{

				if($size){
					print "\t". &INVsize($si[2], $si[3], $total);
				}

				if(($si[2] < $term) && ($si[3] > $term)){
					print "\tSymmetric";
					if($size){
						my$leftside=0;
						my$rightside=0;
						if($si[2] < ($term - $si[2])){
							$leftside = $total - $si[3];
							print "\t".$si[2]."\t".$leftside."\tori";
						}else{
							$rightside = $term - $si[2];
							$leftside =  $si[3] - $term;
							print "\t".$rightside."\t". $leftside."\tterm";
						}
					}
					print "\n";

				}else{
					print "\tAsymmetric";
					if($size){ print "\t-\t-\t-"};
					print "\n";
				}
			}
		}
	}
}
#####################
sub INVsize {
	my$size = ($_[1] - $_[0]);
	if( $size > ($_[2] / 2) ){
		$size = ( $_[2] - $_[1] + $_[0]); #ori
	}
	return( $size );
}

sub HELP_MESSAGE { die "
.Description:
   Characterize single inversions in pairwise mauve backbone files as 'symmetric' or 'asymmetric'.

.Usage: $0 -in [in.txt] -t [int] > out.txt

   [mandatory]
	 -in	<in.txt>	Mauve backbone file, probably the simplified coordinate table from from 'mauve.collinear-check.pl'.
	 -t	<int>		Coordinate of replication terminus in Seq1.

   [optional]
	 -size		Also output inversion size; left, right, and total (requires -seq1, conflicts with -all).
	 -seq1	<int>	Seq1 total length.
	 -all		Don't predict symmetry, just output distance from ori/term for all block boundaries (requires -seq1, conflicts with -size).

   [dependencies]


" }
