#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Statistics::R;
use File::Basename;

&GetOptions(	'fasta=s' => \my$infile,		#
		'oric=s' => \my$oric,
		'term=s' => \my$term,
		'log' => \my$log,
		);		#
($infile and $oric and $term) or &HELP_MESSAGE;

my$seqlen = qx( FastA.length.pl $infile);
chomp$seqlen;
my@seq = split("\t",$seqlen);
my$c_term = &FIND_TERM( $term, $infile );
my$c_ori = &FIND_ORIC( $oric, $infile );
my$offset = &OFFSET( $seq[1], $c_ori);

my@replichores=();
$replichores[0] = $c_term - $offset;
$replichores[1] = $seq[1] - $c_term + $offset;
my@sorted =  sort {$a <=> $b} @replichores;
my$out = ($sorted[1]/$sorted[0]);

print $seq[0]."\t";
if($log){
	print sprintf("%.4f",log($out));
}else{
	print sprintf("%.4f", $out);
}
print "\n";


######################################

sub OFFSET {
	my($seqlen, $ori) = @_;
	my$offset=0;
	unless($ori == 1){
		if( ($seqlen - $ori) < $ori ){ #left of 1
			$offset = ($seqlen - $ori) * -1;
		}else{
			$offset = $ori;
		}
	}
	return($offset);
}


sub FIND_TERM {
	my$term = qx( blastn -task 'blastn-short' -outfmt '6 sstart send' -query $_[0] -subject $_[1] -qcov_hsp_perc 75 | head -1 );
	chomp$term;
	my@coords = split("\t",$term);
	my$center = ($coords[0] + $coords[1])/2;
	return($center);
}

sub FIND_ORIC {
	my$ori = qx( blastn -outfmt '6 sstart send' -query $_[0] -subject $_[1] -qcov_hsp_perc 75 );
	chomp$ori;
	my@coords = split("\t",$ori);
	my$center = ($coords[0] + $coords[1])/2;
	return($center);
}

sub HELP_MESSAGE { die "
.Description:
   Calculates replichore size ratio (Left/Right), normalized to always be >= 1.

.Usage: $0 -fasta [in.fasta] -oric [oric.fasta] -term [term.fasta]

   [mandatory]
	 -fasta	<in.fasta>	Input genome (single contig only).
	 -oric	<ori.fasta>	Fasta file of known OriC to find coordinates.
	 -term	<term.fasta>	Fasta file of known term or dif seq to find coordinates.

   [optional]
	 -log			Return natural log of ratio.

   [dependencies]


" }
