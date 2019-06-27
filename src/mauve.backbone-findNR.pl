#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

&GetOptions(	'indir=s' => \my$indir,		#
		'list=s' => \my$list,
		'ext=s' => \my$ext,
		'outdir=s' => \my$outdir);		#
($indir and $outdir and $list and $ext) or &HELP_MESSAGE;

my%hash = ();

open LIST, "$list";
while(my$l = <LIST>){
	chomp$l;
	$hash{ $l } = 1;
}
close LIST;

my$infiles = qx( find $indir -name "*$ext*" );
my@files = split("\n",$infiles);
print @files . "\n";

foreach my$f (@files){
	#my($ext) = $f =~ /(\.xmfa[\.\w]+$)/;
	my$name = basename($f, $ext); # (".xmfa.backbone",".xmfa.backbone.clean"));
	my@genomes=split("-",$name);
	if((exists($hash{ $genomes[0] })) && (exists($hash{ $genomes[1] })) ){
		#print $f."\t". $name."\n";
		system("cp $f $outdir/ ")
	}

}



sub HELP_MESSAGE { die "
.Description:
   Find and copy mauve (or related) files for pairwise alignment of genomes in list from a directory of many alignment output files.

.Usage: $0 -indir [dir] -outdir [dir] -list [list.txt] -ext [str]

   [mandatory]
	 -indir		Source directory of many alignment files.
	 -ext		File extension to match (eg. 'backbone' or 'backbone.clean').
	 -outdir	Destination directory to copy files.
	 -list		List of genomes.

   [optional]


   [dependencies]

" }
