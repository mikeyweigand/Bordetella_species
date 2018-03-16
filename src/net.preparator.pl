#!/usr/bin/perl -w
use strict;
use Getopt::Long;


&GetOptions(	'inverts=s' => \my$infile,		#
		'sizes=s' => \my$sizes,
		'gaps=s' => \my$gaps,
		'simple' => \my$simple,
		'rename=s' => \my$rename,
		't=s' => \my$term);		#
($infile and $term) or &HELP_MESSAGE;







sub HELP_MESSAGE { die "
.Description:
   Prepare comma-separated edge list for network analysis of single inverts and gaps.

.Usage: $0 -inverts [in.txt] > out.csv

   [mandatory]
	 -inverts	<in.txt>	Tab-separated file of pairwise invert info, probably concatenated outputs from 'mauve.invert-symmetric.pl'.

   [optional]
	 -simple		Input table lacks symmetry data, maybe the output from 'mauve.collinear-filter.pl'.
	 -gaps	<gaps.txt>	Tab-separated file of pairwise gap info, probably the output from 'mauve.collinear-filter.pl'.
	 -rename	<names.tsv>	Tab-separated file to rename nodes (expects: 'old' 'new').
	 -sizes		<sizes.tsv>	Tab-separated file of node sizes, probably cluster sizes from 'mauve.collinear-mcl.pl' (WARNING: must match _new_ name if provided with -rename).

   [dependencies]


" }
