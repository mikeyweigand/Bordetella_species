# Bordetella Genomics Project Outline
------
Revised: 2017-01-18
------
## Project goals:
1. Produce reference-quality genome assemblies for various Bordetella species.
1. Look for rearrangement in genomes of other Bordetella species.
1. Identify functional constraints to rearrangement absent in *B. pertussis*.

##### Hypothesis:
Repeat content (IS481), not functional gene content, determines frequency/extent of genome rearrangment in Bordetella species. -- But does this explain very high rearrangement in *B. pertussis* specifically? or perhaps reflects (diversifying) selection pressure?

## Project tasks:
1. Assemble genomes of various species
	+ *B. bronchiseptica*: 20
	+ *B. parapertussis*: 16
	+ *B. holmesii*: 13
	+ *B. hinzii*: 9
	+ *B. avium*: 8
	+ *B. trematum*: 2
	+ *B. pertussis*: representative subset (n=12?)
	+ should also include NCBI genomes? (Bb=4, Bpp=2, Bho=2, Ba=1, Bt=1)
1. kSNP phylogeny
	+ all and within species
	+ Compare pairwise SNP distributions within each species.
1. Whole-genome alignment within each species.
	+ progressiveMauve (sw=16, hmm=0.85)
	+ Identify rearrangements.
		+ Yes/no, genes at boundaries, etc.
	+ Evaluate gene loss/gain (as within-species variation?)
1. Functional gene content comparison between species.
	+ Match orthologs and annotate variable functions (protein seqs).
1. Repeat content within genomes.
	+ k-mer frequency distribution (jellyfish)
	+ ortholog frequency/duplication (cdhit-est)
1. Compare methylation motifs between species(?).
