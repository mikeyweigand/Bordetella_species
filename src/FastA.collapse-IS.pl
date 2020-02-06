#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Bio::SeqIO;

my$flank = 0;
&GetOptions(	'in=s' => \my$fasta,		#
		'list=s' => \my$list,
		'query=s' => \my$query,
		'q' => \my$q,
		'flank=s' => \$flank,
		'out=s' => \my$out);		#
($fasta and $query and $out) or &HELP_MESSAGE;

#[1] blastn search
my$bls = qx( blastn -query $query -subject $fasta -outfmt '6 qseqid sseqid length pident sstart send' -qcov_hsp_perc 90 -perc_identity 90 );
chomp$bls;
my@hits = split("\n",$bls);
my%iselements=();

#[2] store blastn hits in hash
foreach my$ise (@hits){
	my@IShit = split("\t",$ise);
	if($IShit[-2] < $IShit[-1]){
		@{ $iselements{ $IShit[-2]} } = ($IShit[-1],1); #forward hit
	}else{
		@{ $iselements{ $IShit[-1]} } = ($IShit[-2],-1);	#reverse hit
	}
}

#[3] Look for neighboring matches
my%dups = ();
my@sorted = sort{$a <=> $b}keys%iselements;
for(my$c = 0; $c < @sorted-1; $c++){
	my@dup = &IS_DUP( $sorted[$c], $sorted[$c+1], \%iselements );
	if(scalar@dup > 1){
		#@{ $dups{$dup[0]} } = ($dup[1],$dup[2]);
		@{ $dups{ $dup[0]} } = ($dup[1],$dup[2],1,$sorted[$c], $iselements{$sorted[$c]}[0] );

	}
}

#[4] Consolidate duplications (triplications, etc.)
my@sorted2 = sort{$b <=> $a}keys%dups;
for(my$c2 = 0; $c2 < @sorted2-1; $c2++){
	my@dup2 = &IS_DUP( $sorted2[$c2+1], $sorted2[$c2], \%dups );
	if(scalar@dup2 > 1){
		#@{ $dups{ $sorted2[$c2+1] }} = ($dup2[1],$dup2[2]);
		$dups{ $sorted2[$c2+1] }[0] = $dup2[1];
		$dups{ $sorted2[$c2+1] }[1] = $dup2[2];
		$dups{ $sorted2[$c2+1] }[2] += 1;
		delete $dups{ $sorted2[$c2] };
	}
}

#[5] Delete neighboring IS sequences.
my$seqI = Bio::SeqIO->new(-file=> "$fasta", -format=> 'fasta');
my$original = '';
my$header = '';
my$count=0;
while(my$seq = $seqI->next_seq){
  if($count == 1){ &ERROR };
  $original .= $seq->seq();
	$header = $seq->display_id();
  $count++;
}
my$seqIN = Bio::Seq->new(-seq => $original, -display_id=>"original", -format=>'fasta');
my$seqOUT = '';
my$del=0;
my$deltotal=0;
if(scalar(keys%dups) == 0){
	$seqOUT = $seqIN;
}else{
	my@sortdups = (sort{$a <=> $b}keys%dups);
	my$end = $seqIN->length;
	my$c1=0;
	my$c2=0;
	for (my$j=0; $j < @sortdups; $j++){
		if($j==0){
			$c1 = $sortdups[$j] - 1 ;
			$seqOUT .= $seqIN->subseq(1,$c1);
			#print "1--".$c1."\t";
			$dups{$sortdups[$j]}[5] = ($dups{$sortdups[$j]}[3] - $deltotal);
			$dups{$sortdups[$j]}[6] = ($dups{$sortdups[$j]}[4] - $deltotal);
			$del =($dups{$sortdups[$j]}[0] - $sortdups[$j] + 1 + $flank );
			#$dups{$sortdups[$j]}[7] = $deltotal;
			$deltotal += $del;
			#print $del."\t".$deltotal."\n";

		}else{
			$c1 = $dups{$sortdups[$j-1]}[0] + 1 + $flank;
			$c2 = $sortdups[$j] - 1;
			$seqOUT .= $seqIN->subseq($c1,$c2);
			#print "$c1--$c2\t";
			$dups{$sortdups[$j]}[5] = ($dups{$sortdups[$j]}[3] - $deltotal);
			$dups{$sortdups[$j]}[6] = ($dups{$sortdups[$j]}[4] - $deltotal);
			#$dups{$sortdups[$j]}[7] = $deltotal;
			$del =($dups{$sortdups[$j]}[0] - $sortdups[$j] + 1 + $flank );
			$deltotal += $del;
			#print $del."\t".$deltotal."\n";


		}
	}
	$c1 = $dups{$sortdups[-1]}[0] + 1 + $flank;
	$seqOUT .= $seqIN->subseq($c1,$end);
	#print "$c1--$end\n";
}

#[6] Output new sequence.
open OUT, ">$out";
print OUT ">".$header."-IScollapsed\n";
for( my$o=0; $o <= length($seqOUT); $o+=60){
        print OUT substr($seqOUT,$o,60) . "\n";
}
close OUT;

unless($q){
	print "\n$fasta:\n";
	print "\tQuery: $query\n";
	print "\tCollapsed loci with neighboring IS-elements: ".scalar(keys%dups)."\n";
	print "\tOriginal length = ".$seqIN->length."\n";
	print "\tNew length = ".length($seqOUT)."\n\n";
}

#[7] Output coords of merged IS-elements
if($list){
	open LIST, ">$list";
	foreach my$log (sort{$a <=> $b}keys%dups){
		print LIST join(",",(@{ $dups{$log} }[3,4,5,6,1,2] ) )."\n";
	}
	close LIST;
}

#############
sub IS_DUP {
	my$coord = $_[0];
	my$next = $_[1];
	my%hash = %{ $_[-1] };
	my@out=("nope");
	if($hash{$coord}[1] == $hash{$next}[1]){
		if( abs( $next - $hash{$coord}[0] ) <= 10 ){
			@out = ($next, $hash{$next}[0], $hash{$next}[1]);
		}
	}
	return( @out );
}

sub ERROR {
  print "\nERROR: Looks like your input file has >1 sequence. Sorry, this script can only handle single-contig genomes.\n\n";
  &HELP_MESSAGE;
}

sub HELP_MESSAGE { die "
.Description:
   Collapses neighboring (duplicated, triplicated, etc.) IS-elements into single insertions.

.Usage: $0 -in [fasta] -out [fasta] -query [fasta] -list [list]

   [mandatory]
	 -in	<fasta>	Single-contig genome sequence in fasta format.
	 -out	<fasta>	New genome sequnece fasta with neighboring ISEs collapsed.
	 -query <fasta>	Query sequence file for blastn.

   [optional]
	 -list	<list>	Output csv of collapsed ISE coords.
	 -flank	<int>	Number of 3' flanking bp to also remove from neighboring insertions (eg. 6).
	 -q		Run quietly.

   [dependencies]
	 BioPerl	(Bio::SeqIO)
	 blastn		(Must be in your \$PATH)
" }
