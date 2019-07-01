# Bordetella Genomics Project Outline  
------  
Revised: 2019-07-01  

------
## Project goals:
1. Produce reference-quality genome assemblies for various Bordetella species.
1. Look for rearrangement in genomes of other Bordetella species.
1. Predict all possible symmetric inversions based on IS element content.

##### Hypothesis:
Repeat content (IS481) determines frequency/extent of genome rearrangment in Bordetella species. -- But does this explain very high rearrangement in *B. pertussis* specifically? or perhaps reflects (diversifying) selection pressure?

## Project tasks:
1. Assemble genomes of various species
	+ *B. bronchiseptica*: 12
	+ *B. parapertussis*: 15
	+ *B. holmesii*: 14
	+ *B. hinzii*: 3
	+ *B. avium*: 0
	+ *B. trematum*: 2
	+ *B. pertussis*: All CDC available (n=469)
1. Phylogeny
	+ Genus: mash + mashtree (NJ)
	+ Within species: kSNP3. Compare pairwise SNP distributions within each species   
1. Whole-genome alignment within each species.
	+ progressiveMauve (sw=16, hmm=0.85)
	+ Identify rearrangements.
		+ Yes/no, genes at boundaries, etc.
		+ Circos plots
1. Repeat content within genomes.
	+ k-mer frequency distribution (jellyfish)
	+ ortholog frequency (cdhit-est)
	  + match ISE stability among genomes?
