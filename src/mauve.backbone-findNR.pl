#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

&GetOptions(	'indir=s' => \my$indir,		#
		'list=s' => \my$list,
		'outdir=s' => \my$outdir);		#
($indir and $outdir and $list) or &HELP_MESSAGE;

my%hash = ();

open LIST, "$list";
while(my$l = <LIST>){
	chomp$l;
	$hash{ $l } = 1;
}
close LIST;

my$infiles = qx( find $indir -name "*backbone*" );
my@files = split("\n",$infiles);
#print @files . "\n";

foreach my$f (@files){
	my($ext) = $f =~ /(\.xmfa[\.\w]+$)/;
	my$name = basename($f, $ext); # (".xmfa.backbone",".xmfa.backbone.clean"));
	my@genomes=split("-",$name);
	if((exists($hash{ $genomes[0] })) && (exists($hash{ $genomes[1] })) ){
		#print $f."\t". $name."\n";
		system("cp $f $outdir/ ")
	}

}



sub HELP_MESSAGE { die "
.Description:
   Find and copy mauve backbone files for pairwise alignment of genomes in list.

.Usage: $0 -indir [dir] -outdir [dir] -list [list.txt]

   [mandatory]
	 -indir		Directory of all backbone files.
	 -outdir	Directory to copy backbone files.
	 -list		List of genomes.

   [optional]


   [dependencies]

" }
