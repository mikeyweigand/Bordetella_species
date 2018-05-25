#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

&GetOptions(	'bbone=s' => \my$bbone,		# progressiveMauve backbone file
		'q' => \my$q,
		'h' => \my$h,
		'seq0=s' => \my$seq0,
		'seq1=s' => \my$seq1,
		'name0=s' => \my$name0,
		'flank=s' => \my$bflank,
		'name1=s' => \my$name1,
		'query=s' => \my$query,
		'out=s' => \my$out);		# output matrix of matched IS-elements

($bbone and $seq0 and $seq1 and $out and $query and $name1 and $name0 and $bflank) or &HELP_MESSAGE;
if($h){ &HELP_MESSAGE};

my$genome0 = basename($seq0, ".fasta");
my$genome1 = basename($seq1, ".fasta");
my$len0 = &FASTA_LEN( $seq0 );
my$len1 = &FASTA_LEN( $seq1 );
my%matches = ();
open BBONE, "$bbone";
open OUT, ">$out";
while(my$bb = <BBONE>){
	if($bb =~ /^seq/){

	}else{
		#print "\n".$bb;
		chomp$bb;
		my@sbb = split("\t",$bb);

		my$rtseq0bls='';
		my$rtseq1bls='';
		my$lfseq0bls='';
		my$lfseq1bls='';
		#my$bflank = 1200;


		if(($sbb[0] < 0 ) || ($sbb[2] < 0)){
			$lfseq0bls = &doBLASTn( $query, $seq0, (abs($sbb[0]) - $bflank), (abs($sbb[0]) + $bflank) );
			$lfseq1bls = &doBLASTn( $query, $seq1, (abs($sbb[2]) - $bflank), (abs($sbb[2]) + $bflank) );
			$rtseq0bls = &doBLASTn( $query, $seq0, (abs($sbb[1]) - $bflank), (abs($sbb[1]) + $bflank) );
			$rtseq1bls = &doBLASTn( $query, $seq1, (abs($sbb[3]) - $bflank), (abs($sbb[3]) + $bflank) );

		}elsif($sbb[0] == 0){
			$lfseq1bls = &doBLASTn( $query, $seq1, (abs($sbb[2]) - $bflank), (abs($sbb[2]) + $bflank) );
			$rtseq1bls = &doBLASTn( $query, $seq1, (abs($sbb[3]) - $bflank), (abs($sbb[3]) + $bflank) );

		}elsif($sbb[2] == 0){
			$lfseq0bls = &doBLASTn( $query, $seq0, (abs($sbb[0]) - $bflank), (abs($sbb[0]) + $bflank) );
			$rtseq0bls = &doBLASTn( $query, $seq0, (abs($sbb[1]) - $bflank), (abs($sbb[1]) + $bflank) );

		}

		if(($sbb[0] < 0 ) || ($sbb[2] < 0) || ($sbb[0] == 0) || ($sbb[2] ==0)){
			my@lf0bls = split("\n",$lfseq0bls);
			my@lf1bls = split("\n",$lfseq1bls);
			my@lmatch = &BOUNDRY_MATCH( \@lf0bls, \@lf1bls );
			my@rt0bls = split("\n",$rtseq0bls);
			my@rt1bls = split("\n",$rtseq1bls);
			my@rmatch = &BOUNDRY_MATCH( \@rt0bls, \@rt1bls );

			print OUT join("\t",($bb,
					($name0."_".$lmatch[0]),
					($name0."_".$rmatch[0]),
					($name1."_".$lmatch[1]),
					($name1."_".$rmatch[1])		))."\n";

			unless($q){
				print join("\t",($bb,
					($name0."_".$lmatch[0]),
					($name0."_".$rmatch[0]),
					($name1."_".$lmatch[1]),
					($name1."_".$rmatch[1])		))."\n";
			}

		}

	}
}

#############################
sub BOUNDRY_MATCH {
	my@bls0 = @{$_[0]};
	my@bls1 = @{$_[1]};
	my%h0 = ();
	my%h1 = ();
	my@match = ();
	my@sh=();
	if(scalar@bls0 == 0){
		@match = ("-");
	}elsif(scalar@bls0 > 1){
		@match = ("multi");
	}else{
		foreach my$hit0 (@bls0){
			@sh = split("\t",$hit0);
			# if( exists$h0{$sh[0]} ){
			# 	@match = (("multi-".$sh[0]));
			# }else{
			# 	@{ $h0{$sh[0]} } = ($sh[-2], $sh[-1]);
			# }
			@match = (($sh[0]."_".$sh[-2]."-".$sh[-1]));
		}
	}
	if(scalar@bls1 == 0){
		push(@match, "-");
	}elsif(scalar@bls1 > 1){
		push(@match, "multi");
	}else{
		foreach my$hit1 (@bls1){
			@sh = split("\t",$hit1);
		# if( exists$h1{$sh[0]} ){
		# 	@match = (("multi-".$sh[0]));
		# }else{
		# 	@{ $h1{$sh[0]} } = ($sh[-2], $sh[-1]);
		# 	if( exists($h0{$sh[0]}) ){
		# 		push(@match,($sh[0]."_".$h0{$sh[0]}[0]."-".$h0{$sh[0]}[1], $sh[0]."_".$sh[-2]."-".$sh[-1] ));
		# 	}
		# }
			push(@match, ($sh[0]."_".$sh[-2]."-".$sh[-1]));
		}
	}
	return( @match );
}
sub doBLASTn {
	my($query, $subject, $start, $stop) = @_;
	my$blaster = qx( blastn -query $query -subject $subject -subject_loc $start-$stop -outfmt '6 qseqid sseqid pident sstart send' -qcov_hsp_perc 99 -perc_identity 90 );
	chomp$blaster;
	return( $blaster );
}

sub FASTA_LEN {
	open FA, "$_[0]";
	my$len = 0;
	while(<FA>){
		unless(m/^>(\S+)\s?/){
         s/[^A-Za-z]//g;
         $len+= length $_;
      }
	}
	close FA;
	return( $len );
}

sub HELP_MESSAGE { die "
.Description:
   Matches boundary-flanking elements between 2 genomes aligned with progressiveMauve.

.Usage: $0 -bbone in.bbone -out matrix.txt -seq0 genome0.fasta -seq1 genome1.fasta -query is.fasta -flank [int]

   [mandatory]
   -bbone	<in.bbone>	Input backbone file from progressiveMauve.
   -query	<is.fasta>	Query sequence file of IS-element(s).
   -out		<matrix.txt>	New matrix of matched IS-element insertions.
   -seq0	<seq0.fasta>	Fasta file matching 'seq0' in backbone file.
   -seq1	<seq1.fasta>	Fasta file matching 'seq1' in backbone file.
   -name0	<name>		Prefix name for 'seq0' matches (eg. strain name; 'A123').
   -name1	<name>		Prefix name for 'seq1' matches.
   -flank	<int>		Number of flanking bp to search.

   [optional]
   -q				Run quietly.
   -h				This helpful message.

" }
