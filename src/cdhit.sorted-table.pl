#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

my$min=1;

&GetOptions(	'in=s' => \my$in,		#
		'min=s' => \$min,
		'q' => \my$q,
		'out=s' => \my$out);		#
($in and $out) or &HELP_MESSAGE;



open IN, "$in";
open OUT, ">$out";
my$cluster='';
my%info=();
my$id='';
while(my$i = <IN>){
	chomp$i;
	if($i =~ /^>Cluster/){
		($id) = $i =~ /(\d+)$/;
		#print $i."\t".$id."\n";
		@{$info{ $id }} = (0,0,0);
	}else{
		#$info[0]++;
		@{$info{ $id }}[0]+=1;
		if($i =~ /\*$/){
			my($replen) = $i =~ /(\d+)(?=nt)/;
			my($repid) = $i =~ /(?<=>)(\w+)(?=\.)/;
			@{$info{ $id }}[1] = $replen;
			@{$info{ $id }}[2] = $repid;
			#print join("\t",($i,$replen,$repid))."\n";
		}
	}
}
my$counter=0;
foreach my$o (sort{$a <=> $b}(keys(%info))){
	if( @{$info{ $o }}[0] >= $min){
		print OUT "Cluster-".sprintf("%03d", $o)."\t";
		print OUT join("\t",@{$info{ $o }})."\n";
		$counter++;
	}
}

unless($q){
	print "\n\tFound $counter/". keys(%info)." clusters with >= $min genes.\n\n";
}





sub HELP_MESSAGE { die "
.Description:
   Makes a summary table from a sorted cdhit cluster file.

.Usage: $0 -in [in] -out [out]

   [mandatory]
	 -in	<in>	Input cluster file. MUST be sorted.
	 -out	<out>	Output summary table.

   [optional]
	 -min	<int>	Minimum cluster size (Default = 1).
	 -q		Run quietly.

   [dependencies]

" }
