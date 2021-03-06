## Workflow for rearrangement analysis:  
For characterizing genome rearrangements in complete assemblies of *Bordetella* species as in [Weigand *et al.* 2019](https://msystems.asm.org/content/4/6/e00702-19). In principle, this analysis could be performed on all genomes at once, but for runtime (and sanity) it is recommended that each species under consideration be aligned and clustered separately. Modified from [Weigand *et al.* 2017](https://jb.asm.org/content/199/8/e00806-16.long)  

1. Exhaustive pairwise alignment with __progressiveMauve__ ([Darling *et al.* 2010](http://www.ncbi.nlm.nih.gov/pubmed/20593022)).  
 + __Pertussis_n257/Current/pairwise-mauve-all.pl__  
 + __Pertussis_n257/Current/pairwise-mauve-parallel.sh__  
1. Identify and cluster colinear genomes.  
 + Follow __Pertussis_n257/Current/WORKFLOW-rearrangement.md__
1. Identify a non-redundant subset of genomes which represent each of the clusters and then find all their pairwise alignment files.
 + `cut -f3 [mcl-clusters.txt]`  
 + __mauve.backbone-findNR.pl__
1. Identify single inversion or insertion/deletion alignment pairs.  
 + __Pertussis_n257/Current/mauve.collinear-check.pl__  
 + __mauve.collinear-filter.pl__  
    + Single inversions = `-g 0 -i 1`  
    + Single insertion/deletion = `-g 1 -i 0`  
1. Characterize symmetry of single inversion alignment pairs.
 + Find the coordinate positions of the replication terminus, maybe by __blastn__ alignment of known *dif* sequences.
 + __mauve.invert-symmetric.pl__  
1. Construct network from alignment pairs and analyze.
 + __net.preparator.pl__
 + Follow __networks.R__

 ---
 NOTES and WARNINGS:
 + The expected input for this workflow is a collection of complete (1-contig) bacterial genome assemblies, each saved in a separate fasta file. Intermediate steps use filename parsing "GenomeA-GenomeB" and so the original fasta filenames must contain only [A-Za-z0-9_].
 + To date, these scripts have only been tested with *Bordetella* species and may require modification for use with other species.
