#!/usr/bin/perl -w
use strict;
use Getopt::Long;

&GetOptions(	'mcl=s' => \my$mcl,		#
		'd=s' => \my$dist,
		'names=s' => \my$names,
		'filter=s' => \my$filter,
		'h' => \my$h,
		'help' => \my$help,
		'sum' => \my$summary,
		#'t=s' => \my$term,
		);		#

($mcl and $dist and $names) or &HELP_MESSAGE;
if($h or $help){ &HELP_MESSAGE };

my%clusters = ();
my%filters = ();

if($filter){
	open FILTER, "$filter";
	while(my$f = <FILTER>){
		chomp$f;
		$filters{ $f } = 1;
	}
	close FILTER;
}

open MCL, "$mcl";
while(my$m = <MCL>){
	unless($m =~ /^Single/){
		chomp$m;
		my@sm=split("\t",$m);
		my$clustr = shift@sm;
		$clustr =~ s/Cluster-//g;
		my$clust = "Cluster-" . sprintf("%02d", $clustr);
		shift@sm;

		foreach my$iso (@sm){
			my($pdlid) = $iso =~ /(?<=Bp_)(\w+)(?=_Illumina)/;
			#print join("\t",($clust,$iso,$pdlid))."\n";

			if($filter){
				if( exists( $filters{$pdlid} )){
					$clusters{$pdlid} = $clust;
				}else{
					#print join("\t",($pdlid,$clust)). "\n";
				}
			}else{
				$clusters{$pdlid} = $clust;
			}

		}
	}
}
close MCL;

open DIST, "$dist";
while(my$din = <DIST>){
	chomp$din;
	my@sdin = split("\t",$din);
	if( exists($clusters{$sdin[0]}) && exists($clusters{$sdin[1]}) ){
		if( $clusters{$sdin[0]} eq $clusters{$sdin[1]} ){
			print $din."\t". $clusters{$sdin[0]}. "\n";
		}
	}

}




sub HELP_MESSAGE { die "
.Description:
   Determine pairwise distances among isolates clustered in 'colinear' groups.

.Usage: $0 -mcl [mcl.tsv] -d [dist.tsv] -names [names.tsv] > out.txt

   [mandatory]
	 -mcl	<mcl.tsv>	Input of clusters, probably the output from 'mauve.collinear-mcl.pl'.
	 -d	<dist.tsv>	3-column table of pariwise distances between all isolates.
	 -names	<names.tsv>	2-column table matching isolate names between input files.

   [optional]
	 -sum			Only output summary values per cluster, not individual values for each isolate pair.
	 -filter <list.txt>	Only output pairwise distances among isolates on this list.


   [dependencies]


" }
