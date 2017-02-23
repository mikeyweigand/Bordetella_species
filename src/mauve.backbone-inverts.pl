#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

my%coords=();

while(my$in = <>){
	unless($in =~ /^seq/){
		chomp$in;
		my@sin = split("\t",$in);
		unless($sin[0] == 0 || $sin[2] == 0){
			$coords{ abs($sin[0]) }[0] = $sin[0];
			if($sin[0] > 0){
				$coords{ abs($sin[0]) }[1] = 0; #positive
			}else{
				$coords{ abs($sin[0]) }[1] = 1;	#negative
			}
		}
	}
}
my$n=1;
my$sign=0;
my@sorted = sort{$a <=> $b}(keys(%coords));
foreach my$s (@sorted){
	#print $s."\t".$coords{$s}[0]."\t".$coords{$s}[1]."\n";
	unless( $coords{$s}[1] == $sign ){
		print $n."\n";
		$sign = $coords{$s}[1];
		$n = $coords{$s}[0];
	}
}
print $n."\n";
