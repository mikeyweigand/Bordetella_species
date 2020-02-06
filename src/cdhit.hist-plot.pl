#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Statistics::R;


&GetOptions(	'in=s' => \my$inlist,		#
		'density' => \my$densplot,
		'out=s' => \my$out);		#
($inlist and $out) or &HELP_MESSAGE;

my$ind=1;
my@xcoord=();
my@ycoord=();
my@cx=();
my@genomes=();
my@ticks=();
my@dens1=();
my@dens2=();

open IN, "$inlist";
while(my$i = <IN>){
	chomp$i;
	my$name = basename($i);
	$name =~ s/\.out[\w\W]+$//g;
	push(@genomes, $name);
	push(@ticks, $ind);
	#print join("\t", ($i,$name,$ind))."\n";

	open HIST, "$i";
	while(my$hist = <HIST>){
		chomp$hist;
		if($hist =~ /^\d/){
			my@sh = split("\t",$hist);

			if($sh[1] > 0){
				for(my$d=0; $d<$sh[2]; $d++){
					push(@dens1, $name);
					push(@dens2, $sh[0]);
				}

				if($sh[0] > 1){
					#print join("\t",($name,$ind,$sh[0],$sh[2]))."\n";
					push(@xcoord, $ind);
					push(@ycoord, $sh[0]);
					push(@cx, $sh[2]);
				}
			}
		}
	}
	$ind++;
}

#set largest value to 10, scale other values accordingly
my@cxsorted = sort{$a <=> $b}@cx;
#print $cxsorted[-1]."\n";
my$m = 10 / $cxsorted[-1];
#print $m * $cxsorted[-1]."\n";
foreach my$x (@cx) { $x = $x * $m; };
for(my$ww = 0; $ww < scalar(@dens1); $ww++){
	#print $dens1[$ww]."\t".$dens2[$ww]."\n";
}

my$R = Statistics::R->new();
if($densplot){
	$R->run( qq' library(ggplot2) ' );
	$R->set('d1', \@dens1 );
	$R->set('d2', \@dens2 );
	$R->run( qq' dd <- data.frame(d1,d2) ');
	#my$dfhead = $R->get( qq' colnames(dd[1]) ');
	#my$ddclass = $R->get( qq' class(dd) ');

	#print $dfhead."\t".$ddclass."\n";
	$R->run( qq' png("$out",width=1400,height=800,pointsize=14) ');
	$R->run( qq' qplot( d2, data = dd, geom = "density" , fill=d1, alpha=.2,adjust=0.05) + scale_x_continuous(breaks=seq(0,200,5)) + coord_cartesian(ylim=c(0,0.15),xlim=c(3,70)) ');
	$R->run( qq' dev.off() ' );
}else{
	$R->set('ind', \@xcoord );
	$R->set('counts', \@ycoord );
	$R->set('diam', \@cx );
	$R->set('genomes', \@genomes );
	$R->set('ttt', \@ticks );
	$R->run( qq' png("$out",width=1000,height=800,pointsize=14) ');
	$R->run( qq' plot(ind, counts, type="p",pch=21,bg=rgb(0,1,0,0.3), cex=diam, xlab="Genome index", ylab="Gene occurance (>1)",main = "CDHIT-est", xaxt ="n") ');
	$R->run( qq' axis(side=1, at=ttt, tck=-0.02) ' );
	$R->run( qq' dev.off() ' );
}



sub HELP_MESSAGE { die "
.Description:
   Takes set of cdhit cluster distribution files from cdhit.plot_len.pl and draws a simple plot.

.Usage: $0 -in [in.txt] -out [out.png]

   [mandatory]
	 -in	<in.txt>	List of paths to cluster distribution files.
	 -out	<out.png>	Output image file of plot.

   [optional]
	 -density		Make density plot using 'ggplot2'.


   [dependencies]
	 R (Statistics::R, ggplot2)

" }
