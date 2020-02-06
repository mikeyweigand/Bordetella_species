#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Statistics::R;
use File::Basename;

&GetOptions(	'fasta=s' => \my$infile,		#
		'query=s' => \my$query,
		'inverts=s' => \my$inverts,
		'out=s' => \my$out,
		'oric=s' => \my$oric,
		'temp=s' => \my$temp,
		);		#
($infile and $query and $inverts and $out and $temp) or &HELP_MESSAGE;

my$srcdir = dirname($0);
unless( -x "$srcdir/inverts.pred.R" ){
	print "\nERROR: Cannot find '$srcdir/inverts.pred.R'\n";
	&HELP_MESSAGE;
}
unless( -x "$srcdir/FastA.length.pl" ){
	print "\nERROR: Cannot find '$srcdir/FastA.length.pl'\n";
	&HELP_MESSAGE;
}

my$seqlen = qx( FastA.length.pl $infile | cut -f2);
chomp$seqlen;
my$term = &FIND_TERM( $infile );

my$ori=1;
if($oric){
	$ori = &FIND_ORIC( $oric, $infile );
}
#print join("\t",($term,$ori,$seqlen))."\n";

#find blastn hits for query in subject and find shortest dist to ori/term
my$bls = &BLASTN( $query, $infile );
my@blastout = split("\n",$bls);

my%bls_hits = ();
my@centers = ();
foreach my$s1 (@blastout){
	my@bls2 = split("\t",$s1);
	my$center = ($bls2[0] + $bls2[1]) / 2;

	if($bls2[0] < $bls2[1]){ push(@bls2, 1);
	}else{ push(@bls2,-1) };

	my@inv = &INVDIS($term, $seqlen, $center, $ori);
	@{$bls_hits{ $center }} = (@bls2,@inv);
	push(@centers,$center);
	#print join("__",(@bls2,$center))."\t".join("__", @inv)."\n";
}

# temp file for RScript
open TEMP, ">$temp";
foreach my$pos (@centers){
	my@short = @{$bls_hits{$pos}};
	splice(@short,2,1); #removes qseqid
	#print TEMP $pos."\t". join("\t",(  @{$bls_hits{$pos}} ))."\n";
	print TEMP $pos."\t". join("\t",( @short ))."\n";
}
close TEMP;

# linear regression of observed inverts and predict new boundaries
my$temp2 = $temp;
$temp2 =~ s/\.txt$/\-pred.txt/;
#system( qq( Rscript /home/yrh8/Documents/Bordetella_species/src/inverts.pred.R -r $inverts -f $temp -o $temp2 ) );
system( qq( Rscript $srcdir/inverts.pred.R -r $inverts -f $temp -o $temp2 ) );

# search new boundary predictions for matching blastn hits
open PRED, "$temp2";
open OUT, ">$out";
my%matched = ();
my$count = 0;
while(my$p = <PRED>){
	chomp$p;
	my@sp = split("\t",$p);
	my@window = &WINDOW( \@sp, $term, $seqlen, $ori );
	my@pred=();
	@pred = &MATCHER( \@sp, \@centers, \@window, \%bls_hits);
	#print $sp[0]."\t".scalar(@pred)."\n";

	foreach my$j (@pred){
		unless( exists( $matched{$sp[0]}{$j}) ){

			if($bls_hits{$j}[2] eq $bls_hits{$sp[0]}[2]){	#must have matching qseqids
				print OUT $bls_hits{$j}[2]."\t";
				print OUT join("\t",( $sp[0],$sp[1],$sp[2],$sp[3]))."\t";
				print OUT join("\t",( $j, $bls_hits{$j}[0],$bls_hits{$j}[1],$bls_hits{$j}[3] ))."\t";
				print OUT join("\t",( $bls_hits{$sp[0]}[-2], $bls_hits{$j}[-2] ))."\n";

				$matched{$sp[0]}{$j} = 1;
				$matched{$j}{$sp[0]} = 1;
				$count++;
			}
		}
	}
}
print "\n\tFound $count possible inversions between BLASTn hits of '".basename($query)."' in '".basename($infile)."'\n\n";

######################################

sub INVDIS {
	my($term, $seqlen, $pos, $ori) = @_;
	my@out = (0,0);
	my$offset = &OFFSET( $seqlen, $ori);

	if($term > $pos){
		if( ($term - $pos) > ($pos - $offset) ){
			@out = (($pos - $offset), "RtOri");
		}else{
			@out = (($term - $pos), "RtTerm");
		}
	}else{
		if( ($pos - $term) > ($seqlen + $offset - $pos) ){
			@out = (($seqlen + $offset - $pos), "LfOri");
		}else{
			@out = (($pos - $term), "LfTerm");
		}
	}
	return( @out );
}

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

sub BLASTN {
	my$bls = qx( blastn -query $_[0] -subject $_[1] -outfmt '6 sstart send qseqid' -qcov_hsp_perc 80 | sort -n );
	chomp$bls;
	return($bls);
}

sub WINDOW {
	my@sp = @{$_[0]};
	my$term = $_[1];
	my$seqlen = $_[2];
	my$ori = $_[3];
	my$offset = &OFFSET( $seqlen, $ori );
	my@coords = ('1',$seqlen);

	if($sp[5] eq 'RtOri'){
		if($sp[7] < 0){
			$coords[1] = $seqlen; # + $offset;
		}else{
			$coords[1] = $seqlen + $offset - $sp[7];
		}
		$coords[0] = $seqlen + $offset - $sp[8];

	}elsif($sp[5] eq 'RtTerm'){
		if($sp[7] < 0){
			$coords[0] = $term;
		}else{
			$coords[0] = $term + $sp[7];
		}
		$coords[1] = $term + $sp[8];

	}elsif($sp[5] eq 'LfTerm'){
		if($sp[7] < 0){
			$coords[1] = $term;
		}else{
			$coords[1] = $term - $sp[7];
		}
		$coords[0] = $term - $sp[8];

	}else{ 	#LfOri
		if($sp[7] < 0){
			$coords[0] = 1;
		}else{
			$coords[0] = $sp[7] + $offset;
		}
		$coords[1] = $sp[8] + $offset;
	}

	return( @coords );

}

sub MATCHER {
	my@sp = @{$_[0]};
	my@centers = @{$_[1]};
	my@window = @{$_[2]};
	my%bls_hits = %{$_[3]};
	my$dir=$sp[3];

	my@found=();
	foreach my$c (@centers){
		#my$qseqid = $bls_hits{ $c }[2];
		if($window[0] < $c && $c < $window[1]){
			unless( $dir == $bls_hits{ $c }[3]){	#strand
				push(@found,$c);
			}
		}
	}
	return ( @found );
}

sub FIND_TERM {
	my$term = qx( echo "aattcgcataatgtatattatgtaaagt" | blastn -task 'blastn-short' -outfmt '6 sstart send' -subject $_[0] -qcov_hsp_perc 75 );
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
   Predict all possible symmetric inversions in given genome based on table of observed inversions.

   NOTE: the Bordetella dif sequence is hard-coded into this script and may not be suitable for other organism: 'aattcgcataatgtatattatgtaaagt'

.Usage: $0 -fasta [in.fasta] -query [query.fasta] -inverts [inv.txt] -out [out.txt] -temp [temp.txt]

   [mandatory]
	 -fasta		<in.fasta>	Input genome (single contig only).
	 -query		<q.fasta>	Query sequence for predicting inversion boundaries.
	 -inverts	<inv.txt>	Table of observed inversions for modeling, probably output from 'mauve.invert-symmetric.pl'.
	 -out		<out.txt>	Output table of predicted inversion boundaries.
	 -temp		<tmp.txt	Temp file to store coordiates for RScript.

   [optional]
	 -oric		<ori.fasta>	Fasta file of known OriC to find coordinates, otherwise use 1.

   [dependencies]
	 R			(Statistics::R)
	 inverts.pred.R		(Must be in the same directory as this script)
	 FastA.length.pl	(http://enve-omics.ce.gatech.edu/enveomics/. Must be in your \$PATH)
	 blastn			(Must be in your \$PATH)


" }
