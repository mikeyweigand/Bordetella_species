#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Sort::Key::Natural qw(natsort);
use List::Util qw(max);

&GetOptions(	'bbone=s' => \my$bbone,		# progressiveMauve backbone file
		'q' => \my$q,
		'h' => \my$h,
		'seq0=s' => \my$seq0,
		'seq1=s' => \my$seq1,
		'name0=s' => \my$name0,
		'name1=s' => \my$name1,
		'query=s' => \my$query,
		'min=s' => \my$minblock,
		'sum=s' => \my$summ,		# output optional summary table of IS-element counts in each genome
		'out=s' => \my$out);		# output matrix of matched IS-elements

($bbone and $seq0 and $seq1 and $out and $query and $name1 and $name0) or &HELP_MESSAGE;
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
		unless($q){print $genome0."--vs--".$genome1."\n"};
	}else{
		#print "\n".$bb;
		chomp$bb;
		my@sbb = split("\t",$bb);

		unless(($sbb[0] == 0 ) || ($sbb[2] == 0)){	#ignores alignment gaps
			#blastn search within backbone blocks
			my$seq0bls = &doBLASTn( $query, $seq0, abs($sbb[0]), abs($sbb[1]) );
			my$seq1bls = &doBLASTn( $query, $seq1, abs($sbb[2]), abs($sbb[3]) );
			my@sseq0bls = split("\n",$seq0bls);
			my@sseq1bls = split("\n",$seq1bls);

			#store blastn hits from seq0
			my%seq0=();
			foreach my$hit0 (@sseq0bls){
				my@bls0 = split("\t",$hit0);
				my@rel0 = &REL_COORDS( $sbb[0], $sbb[1], $bls0[-2], $bls0[-1] );
				#print join("\t",($sbb[0],$sbb[1]))."\t||\t". $hit0."\t||\t".join("-",@rel0)."\n";
				@{ $seq0{ $bls0[0] }{$rel0[0] } } = ($rel0[1],$rel0[2],$bls0[-2],$bls0[-1]);
			}

			#look for matches to blastn hits in seq1
			foreach my$hit1 (@sseq1bls){
				my@bls1 = split("\t",$hit1);
				my@rel1 = &REL_COORDS( $sbb[2], $sbb[3], $bls1[-2], $bls1[-1] );
				my$match = &FIND_MATCH( \@rel1, \%seq0, $bls1[0]);
				unless($match == 0){
					#print $bls1[0]."\t". join("\t",@rel1)."\t||\t".join("\t",($match,$seq0{$bls1[0]}{$match}[0],$seq0{$bls1[0]}{$match}[1])). "\n";
					my$m0 = $bls1[0]."_".join("-",($seq0{$bls1[0]}{$match}[2],$seq0{$bls1[0]}{$match}[3]));
					my$m1 = $bls1[0]."_".join("-",( $bls1[-2], $bls1[-1]));
					$matches{ $m0 } = $m1;
					delete $seq0{$bls1[0]}{$match};
				}else{
					unless($q){print "\tWARNING: no match for ".$name1."_".$bls1[0]."_".join("-",($bls1[-2],$bls1[-1]))."\t".join("-",@rel1)."\n" };
				}
			}
			#what's left in the hash of seq0 blastn hits?
			foreach my$miss (keys%seq0){
				foreach my$missing (sort{$a <=> $b}(keys$seq0{$miss})){
					unless($q){
						print "\tWARNING: no match for $name0\_$miss\_".$seq0{$miss}{$missing}[-2]."-".$seq0{$miss}{$missing}[-1]."\t";
						print join("-",($missing,$seq0{$miss}{$missing}[0],$seq0{$miss}{$missing}[1] ))."\n";
				 }
				}
			}

			#also check near boundaries (+/- 1000bp)
			unless(($sbb[0] < 1 ) || ($sbb[2] < 1)){
				my$rtseq0bls='';
				my$rtseq1bls='';
				my$lfseq0bls='';
				my$lfseq1bls='';
				my$bflank = 1200;
				unless(($sbb[1] == $len0) || ($sbb[3] == $len1)){
					$rtseq0bls = &doBLASTn( $query, $seq0, (abs($sbb[1]) - $bflank), (abs($sbb[1]) + $bflank) );
					$rtseq1bls = &doBLASTn( $query, $seq1, (abs($sbb[3]) - $bflank), (abs($sbb[3]) + $bflank) );
				}
				unless(($sbb[0] == 1)){
					$lfseq0bls = &doBLASTn( $query, $seq0, (abs($sbb[0]) - $bflank), (abs($sbb[0]) + $bflank) );
					$lfseq1bls = &doBLASTn( $query, $seq1, (abs($sbb[2]) - $bflank), (abs($sbb[2]) + $bflank) );
				}
				#check for right boundary matches
				my@bmatch=();
				if((length$rtseq0bls > 1) && (length$rtseq1bls > 1)){
					my@rt0bls = split("\n",$rtseq0bls);
					my@rt1bls = split("\n",$rtseq1bls);
					@bmatch = &BOUNDRY_MATCH( \@rt0bls, \@rt1bls );
					for(my$bm = 0; $bm < @bmatch-1; $bm+=2){
						if($bmatch[$bm] eq "multi"){
							print "\tWARNING: Multiple right boundary hits:\t".join("\t",($name0,$sbb[1],$bmatch[$bm+1]))."\n";
						}elsif($bmatch[$bm] eq "none"){
							print "\tWARNING: No right boundary matches found:\t".join("\t",($name0,$sbb[1],$name1,$sbb[3]))."\n";
						}elsif( exists$matches{ $bmatch[$bm] } ){
							print "\tWARNING: Already found right boundary match:\t".join("\t",($name0,$bmatch[$bm],$matches{$bmatch[$bm]}))."\n";
						}else{
							$matches{ $bmatch[$bm] } = $bmatch[$bm+1];
							print "\tFOUND: Right boundary match:\t".join("\t",($name0,$sbb[1],$bmatch[$bm],"||",$name1,$sbb[3],$bmatch[$bm+1]))."\n";
						}
					}
				}
					#check for left boundary matches
				if((length$lfseq0bls > 1) && (length$lfseq1bls > 1)){
					my@lf0bls = split("\n",$lfseq0bls);
					my@lf1bls = split("\n",$lfseq1bls);
					@bmatch = &BOUNDRY_MATCH( \@lf0bls, \@lf1bls );
					#print join("\n",@lf0bls)."\n".join("\n",@lf1bls)."\n";
					for(my$lbm = 0; $lbm < @bmatch-1; $lbm+=2){
						if($bmatch[$lbm] eq "multi"){
							print "\tWARNING: Multiple left boundary hits:\t".join("\t",($name0,$sbb[0],$bmatch[$lbm+1]))."\n";
						}elsif($bmatch[$lbm] eq "none"){
							print "\tWARNING: No left boundary matches found:\t".join("\t",($name0,$sbb[0],$name1,$sbb[2]))."\n";
						}elsif( exists$matches{ $bmatch[$lbm] } ){
							print "\tWARNING: Already found left boundary match:\t".join("\t",($name0,$bmatch[$lbm],$matches{$bmatch[$lbm]}))."\n";
						}else{
							$matches{ $bmatch[$lbm] } = $bmatch[$lbm+1];
							print "\tFOUND: Left boundary match:\t".join("\t",($name0,$sbb[0],$bmatch[$lbm],$name1,$sbb[2],$bmatch[$lbm+1]))."\n";
						}
					}
				}
					# print length$lfseq0bls."\t".length$lfseq1bls."\n";
					# print length$lfseq0bls."\t".length$lfseq1bls."\n";
			}


			# 		if(($rtseq0bls !~ /\n/) && ($rtseq1bls !~ /\n/)){
			# 			my@r0 = split("\t",$rtseq0bls);
			# 			my@r1 = split("\t",$rtseq1bls);
			# 			my$b0 = $r0[0]."_".join("-",($r0[-2],$r0[-1])); #."\t";
			# 			my$b1 = $r1[0]."_".join("-",($r1[-2],$r1[-1])); #."\n";
			# 			unless(exists($matches{ $b0 })){
			# 				$matches{ $b0 } = $b1;
			# 			}else{
			# 				unless($q){print "\tWARNING: Duplicate matches: ".join("\t",($b0,$matches{$b0},$b1))."\n"};
			# 			}
			# 		}else{
			# 			unless($q){
			# 				print "\tWARNING: Multiple matches:\t$name0\t$sbb[1]\t||\t". $rtseq0bls."\n";
			# 				print "\tWARNING: Multiple matches:\t$name1\t$sbb[3]\t||\t". $rtseq1bls."\n";
			# 			}
			# 		}
			# 	}else{
			# 		unless(($sbb[1] == $len0) || ($sbb[3] == $len1)){
			# 			unless($q){
			# 				print "\tWARNING: No RIGHT boundary match:\t$name0\t$sbb[1]\t||\t". $rtseq0bls."\n";
			# 				print "\t\t\t\t\t\t$name1\t$sbb[3]\t||\t". $rtseq1bls."\n";
			# 			}
			# 		}
			# 	}
			# 	if((length$lfseq0bls > 1) && (length$lfseq1bls > 1)){
			# 		if(($lfseq0bls !~ /\n/) && ($lfseq1bls !~ /\n/)){
			# 			my@r0 = split("\t",$lfseq0bls);
			# 			my@r1 = split("\t",$lfseq1bls);
			# 			my$b0 = $r0[0]."_".join("-",($r0[-2],$r0[-1])); #."\t";
			# 			my$b1 = $r1[0]."_".join("-",($r1[-2],$r1[-1])); #."\n";
			# 			unless(exists($matches{ $b0 })){
			# 				$matches{ $b0 } = $b1;
			# 			}else{
			# 				unless($q){print "\tWARNING: Duplicate matches: ".join("\t",($b0,$matches{$b0},$b1))."\n"};
			# 			}
			# 		}else{
			# 			unless($q){
			# 				print "\tWARNING: Multiple boundary matches:\t$name0\t||\t". $lfseq0bls."\n";
			# 				print "\tWARNING: Multiple boundary matches:\t$name1\t||\t". $lfseq1bls."\n";
			# 			}
			# 		}
			# 	}else{
			# 		unless(($sbb[0] == 1) || ($sbb[2] ==1)){
			# 			unless($q){
			# 				print "\tWARNING: No LEFT boundary match:\t$name0\t$sbb[0]\t||\t". $lfseq0bls."\n";
			# 				print "\t\t\t\t\t\t$name1\t$sbb[2]\t||\t". $lfseq1bls."\n";
			# 			}
			# 		}
			# 	}
			# }

		}else{
			#ignore gaps = strain-specific insertions
		}
	}
}
#output matches
unless($q){ print "\tIdentified matches = ". scalar(keys%matches)."\n\n" };
foreach my$o (sort(keys%matches)){
	print OUT $name0."_".$o."\t".$name1."_".$matches{$o}."\n";
}
close BBONE;
close OUT;



#########################################
sub BOUNDRY_MATCH {
	my@bls0 = @{$_[0]};
	my@bls1 = @{$_[1]};
	my%h0 = ();
	my%h1 = ();
	my@match = ();
	my@sh=();
	foreach my$hit0 (@bls0){
		@sh = split("\t",$hit0);
		if( exists$h0{$sh[0]} ){
			@match = ("multi", $sh[0]);
		}else{
			@{ $h0{$sh[0]} } = ($sh[-2], $sh[-1]);
		}
	}
	foreach my$hit1 (@bls1){
		@sh = split("\t",$hit1);
		if( exists$h1{$sh[0]} ){
			@match = ("multi", $sh[0]);
		}else{
			@{ $h1{$sh[0]} } = ($sh[-2], $sh[-1]);
			if( exists($h0{$sh[0]}) ){
				push(@match,($sh[0]."_".$h0{$sh[0]}[0]."-".$h0{$sh[0]}[1], $sh[0]."_".$sh[-2]."-".$sh[-1] ));
			}
		}
	}
	if(scalar@match < 1){
		@match = ("none", 0);
	}
	return( @match );
}

sub doBLASTn {
	my($query, $subject, $start, $stop) = @_;
	my$blaster = qx( blastn -query $query -subject $subject -subject_loc $start-$stop -outfmt '6 qseqid sseqid pident sstart send' -qcov_hsp_perc 90 -perc_identity 90 );
	chomp$blaster;
	return( $blaster );
}

sub REL_COORDS {
	my($block1,$block2,$is1,$is2) = @_;
	my$mid = (abs($block2 + $block1) / 2);
	my@rel=();
	if($block1 > 0){
		if((($block2 - $block1) > 1000000) && ($is1 >= $mid) && ($is2 >= $mid)){
			@rel = &REL_NEG( $block1, $block2, $is1, $is2);
		}else{
			@rel=(($is1-$block1+1),($is2-$block1+1));
		}
		if($is1 < $is2){
			push(@rel, 'FWD');
		}else{
			push(@rel, 'REV')
		}

	}else{
		if((abs($block2 - $block1) > 1000000) && ($is1 <= $mid) && ($is2 <= $mid)){
			@rel = &REL_NEG( $block1, $block2, $is1, $is2);
		}else{
			@rel=((abs($block2)-$is1+1),(abs($block2)-$is2+1));
		}
		if($is1 < $is2){
			push(@rel, 'REV');
		}else{
			push(@rel, 'FWD');
		}
	}
	return @rel;
}

sub REL_NEG {
	my($block1,$block2,$is1,$is2) = @_;
	my@rel=();
	my$mid = (($block2 + $block1) / 2);
	if($block1 > 0){
		@rel=( (-1 * ($block2-$is1)) , (-1 * ($block2-$is2)) );
	}else{
		@rel=( (-1 * ($is1-abs($block1)+1)) , (-1 * ($is2-abs($block1)+1)) );
	}
	return @rel;
}

sub FIND_MATCH {
	my@rel = @{$_[0]};
	my%seq0 = %{$_[1]};
	my$q = $_[2];

	my$flex = 100; #Adjust this to depend on relative coord, larger coord = larger flex, to account for increasing number of indels further from block ends.
	if(abs($rel[0]) <= 50000){
		$flex = 20;
	}elsif(abs($rel[0]) <= 100000){
		$flex = 50;
	}elsif(abs($rel[0]) <= 200000){
		$flex = 75;
	}

	my$return = 0;
	my@seq0keys = keys$seq0{$q};
	foreach my$s0 (@seq0keys) {
		my@seq0hit = @{ $seq0{$q}{ $s0 } };
		if((abs($s0 - $rel[0]) <= $flex ) && (abs($seq0hit[0] - $rel[1]) <= $flex )){
			$return = $s0;
		}
	}
	return( $return );
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
   Matches IS-element insertions between 2 genomes aligned with progressiveMauve.

.Usage: $0 -bbone in.bbone -out matrix.txt -seq0 genome0.fasta -seq1 genome1.fasta -query is.fasta

   [mandatory]
   -bbone	<in.bbone>	Input backbone file from progressiveMauve.
   -query	<is.fasta>	Query sequence file of IS-element(s).
   -out		<matrix.txt>	New matrix of matched IS-element insertions.
   -seq0	<seq0.fasta>	Fasta file matching 'seq0' in backbone file.
   -seq1	<seq1.fasta>	Fasta file matching 'seq1' in backbone file.
   -name0	<name>		Prefix name for 'seq0' matches (eg. strain name; 'A123').
   -name1	<name>		Prefix name for 'seq1' matches.

   [optional]
   -q				Run quietly.
   -h				This helpful message.
   -min		<int>		Minimum block size to consider. (NOT AVAILABLE)

" }
