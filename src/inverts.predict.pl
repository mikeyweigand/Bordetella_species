#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Statistics::R;


&GetOptions(	'fasta=s' => \my$infile,		#
		'query=s' => \my$query,
		'inverts=s' => \my$inverts,
		'out=s' => \my$out,
		'temp=s' => \my$temp,
		);		#
($infile and $query and $inverts and $out and $temp) or &HELP_MESSAGE;

my$seqlen = qx( FastA.length.pl $infile | cut -f2);
chomp$seqlen;
my$term = &FIND_TERM( $infile );
#print $term."\t".$seqlen."\n";

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
	#print join("--",(@bls2,$center))."\n";

	my@inv = &INVDIS($term, $seqlen, $center);
	@{$bls_hits{ $center }} = (@bls2,@inv);
	push(@centers,$center);
}

#print @centers."\n";

# temp file for RScript
open TEMP, ">$temp";
foreach my$pos (@centers){
	print TEMP $pos."\t". join("\t",(  @{$bls_hits{$pos}} ))."\n";
}
close TEMP;

# linear regression of observed inverts and predict new boundaries
my$temp2 = $temp;
$temp2 =~ s/\.txt$/\-pred.txt/;
#print $temp2."\n";

system( qq( Rscript /home/yrh8/Documents/Bordetella_species/src/inverts.pred.R -r $inverts -f $temp -o $temp2 ) );


# search new boundary predictions for matching blastn hits
open PRED, "$temp2";
my%matched = ();
while(my$p = <PRED>){
	chomp$p;
	my@sp = split("\t",$p);

	my@window = &WINDOW( \@sp, $term, $seqlen );
	my@pred=();
	@pred = &MATCHER( \@sp, \@centers, \@window, \%bls_hits);

	#print join("\t",(@sp,@window))."\t".@pred."\n";
	foreach my$j (@pred){

		unless( exists( $matched{$sp[0]}{$j}) ){ #} || exists( $matched{$j}{$sp[0]} ) ){
			print join("\t",( $sp[0],$sp[1],$sp[2],$sp[3]))."\t";
			print join("\t",( $j, $bls_hits{$j}[0],$bls_hits{$j}[1],$bls_hits{$j}[2] ))."\n";

			$matched{$sp[0]}{$j} = 1;
			$matched{$j}{$sp[0]} = 1;

		}
	}
}


######################################

sub INVDIS {
	my($term, $seqlen, $pos) = @_;
	my@out = (0,0);
	if($term > $pos){
		if( ($term - $pos) > $pos ){
			@out = ($pos, "RtOri");
		}else{
			@out = (($term - $pos), "RtTerm");
		}
	}else{
		if( ($pos - $term) > ($seqlen - $pos) ){
			@out = (($seqlen - $pos), "LfOri");
		}else{
			@out = (($pos - $term), "LfTerm");
		}
	}
	return( @out );
}

sub BLASTN {
	my$bls = qx( blastn -query $_[0] -subject $_[1] -outfmt '6 sstart send' -qcov_hsp_perc 80 | sort -n );
	chomp$bls;
	return($bls);
}

sub WINDOW {
	my@sp = @{$_[0]};
	my$term = $_[1];
	my$seqlen = $_[2];
	my@coords = ('1',$_[2]);

	if($sp[5] eq 'RtOri'){
		if($sp[7] < 0){
			$coords[1] = $seqlen;
		}else{
			$coords[1] = $seqlen - $sp[7];
		}
		$coords[0] = $seqlen - $sp[8];

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
			$coords[0] = $sp[7];
		}
		$coords[1] = $sp[8];
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
		if($window[0] < $c && $c < $window[1]){
			unless( $dir == $bls_hits{ $c }[2]){
				push(@found,$c);
			}
		}
	}
	return ( @found );
}

sub FIND_TERM {
	my$term = qx( echo "aattcgcataatgtatattatgtaaagt" | blastn -outfmt '6 sstart' -subject $_[0] );
	chomp$term;
	return($term);
}

sub HELP_MESSAGE { die "
.Description:
   Predict all possible symmetric inversions in given genome based on table of observed inversions.

.Usage: $0 -in [in.txt] -t [int] > out.txt

   [mandatory]
	 -fasta	<in.fasta>	Input genome (single contig only).
	 -query	<q.fasta>		Query sequence for predicting inversion boundaries.
	 -inverts	<inv.txt>	Table of observed inversions for modeling, probably output from 'mauve.invert-symmetric.pl'.
	 -out	<out.txt>		Output table of predicted inversion boundaries.
	 -temp	<tmp.txt	Temp file to store coordiates for RScript.

   [optional]


   [dependencies]
	 inverts.pred.R

" }