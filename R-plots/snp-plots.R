# This script plots pairwise SNP distances within species or colinear structure clusters. Intended to be run interactively in Rstudio.

setwd("~/Documents/Bordetella_species/results/mauve/20180221-Bp/")

library(ggplot2)
library(scales)
library(reshape)

##### Pairwise SNP distances within each of the 32 colinear clusters observed in B. pertussis #
clust.dist <- read.table("./04.colinear-mcl/20180425-Bp-check-gap1500-mcl-cluster-dist.tsv", sep="\t", header=F)

first15 = sort(unique(clust.dist$V4))[1:15]
first12 = sort(unique(clust.dist$V4))[1:12]
head(clust.dist)
(ggplot( data = clust.dist, aes(x=V4, y=V3))
  + geom_jitter(width = 0.25, size=0.5, color="red")
  + geom_boxplot(outlier.size = 0, outlier.stroke = 0, alpha=0)
  #+ scale_x_discrete(limits=first12)
  + scale_y_continuous(limits=c(0,35))
  + theme_minimal(base_size = 16)
  + theme(axis.text.x=element_text(angle=45, hjust=1), axis.text=element_text(color='black'))
  + labs(x="Colinear cluster", y="Pairwise SNP distance")
)
#ggsave("./99.figures/20180430-cluster-SNPs.pdf", device = 'pdf', width = 8, height = 6, units = 'in', useDingbats=F)
#ggsave("./99.figures/Bordetella_species/20180723-cluster-SNPs.pdf", device = 'pdf', width = 6.25, height = 3, units = 'in', useDingbats=F)

##### Alternative density plots of pairwise SNP distances in the top 12 clusters #
# clust.wtf=clust.dist[clust.dist$V4 %in% first12,]
# (ggplot( data = clust.wtf, aes(x=V3, fill=V4, color=V4, alpha=0.5))
#   + geom_density(trim=T, adjust=1)
#   #+ scale_y_continuous(limits=c(0,1))
#   + scale_x_continuous(limits = c(0,35))
#   + theme_classic()
#   + facet_grid(V4 ~ .)
#   + theme(strip.text.x=element_blank(), panel.background=element_blank(),
#           panel.spacing.y=unit(-2.5, "cm"), panel.grid=element_blank(),
#           axis.ticks.y=element_blank(), axis.text.y=element_blank())
#   + labs(y="Colinear cluster", x="Pairwise SNP distance")
# )


##### Temporal sampling of colinear clusters and singletons in B. pertussis #
clust.years <- read.table("./04.colinear-mcl/20180425-Bp-check-gap1500-mcl-years.tsv", sep="\t",header=F)
clust.years2k <- read.table("./04.colinear-mcl/20180425-Bp-check-gap1500-mcl-years2000.tsv", sep="\t",header=F)
clust.years2kp <- read.table("./04.colinear-mcl/20180425-Bp-check-gap1500-mcl-years2000-pre.tsv", sep="\t",header=F)


head(clust.years)
y.years = seq(2000, 2016, 1)
y2 = seq(1935,2015,5)
first12s = sort(first12, decreasing = T)
first12s = as.character(first12s)
first12s = c("Singleton", first12s)
first12s = as.factor(first12s)
rev=sort(unique(clust.years$V3)) #, decreasing =F)

(ggplot( data = clust.years2kp, aes(x=as.factor(V2), y=V3))
  + geom_jitter(width = 0.1, height=0.35, size=1, color="red")
  + scale_y_discrete(limits=sort(unique(clust.years$V3),decreasing=T) )
  + theme_minimal(base_size = 16)
  + theme(axis.text.x=element_text(angle=45, hjust=1),
          panel.grid.minor.x = element_blank(),
          axis.text=element_text(color='black'))
  + labs(y="Colinear cluster", x="Year of Isolation")

  )
#ggsave("./99.figures/20180723-cluster-years.pdf", device = 'pdf', width = 6.25, height = 4.5, units = 'in', useDingbats=F)
#ggsave("./99.figures/20180430-cluster-years-02.pdf", device = 'pdf', width = 6, height = 3.3, units = 'in', useDingbats=F)

#wtf=clust.years2kp[clust.years2kp$V3 %in% first12s,]


##### Bordetella species comparison of within-species SNP distances #

#pairwise.snps <- read.table("../../kSNP/20180501-SNP-dist/20180501-pairwiseSNPs-down-plusComplexes.tsv",sep="\t",header=F);
pairwise.snps <- read.table("../../kSNP/20180501-SNP-dist/20190826-pairwiseSNPs-down-CORRECTED.tsv",sep="\t",header=F);
pairwise.snps$V1 <- factor(pairwise.snps$V1, levels = c("Bbronchiseptica","Bb_complexI","Bb_complexIV","Bholmesii","Bparapertussis","Bpertussis"), ordered=T)

#pairwise.snps$V1
(ggplot( data = pairwise.snps, aes(x=V1, y=V2))
 + geom_jitter(width = 0.2, size=1.5, color="red")
 + geom_boxplot(outlier.size = 0, outlier.stroke = 0, alpha=0, color='black')
 + theme_minimal(base_size = 16)
 + theme(axis.text.x=element_text(angle=25, hjust=1), axis.text=element_text(color='black'))
 + labs(x="Species", y="Pairwise SNP distance")
 + scale_y_continuous( labels = comma)
)
#ggsave("./99.figures/20180726-pairwiseSNPs-down.pdf", device = 'pdf', width = 6, height = 5, units = 'in', useDingbats=F)


# Subset B. pertussis, B.holmesii, and B. parapertussis only
pairwise.pert <-subset(pairwise.snps,V1 == "Bparapertussis" | V1 == "Bpertussis" | V1 == "Bholmesii")
#nrow(pairwise.pert)
#mean(subset(pairwise.snps, V1 == "Bholmesii")[,2])

(ggplot( data = pairwise.pert, aes(x=V1, y=V2))
  + geom_jitter(width = 0.2, size=1.2, color="red")
  + geom_boxplot(outlier.size = 0, outlier.stroke = 0, alpha=0, color='black')
  + theme_minimal(base_size = 14)
  + labs(x="Species", y="Pairwise SNP distance")
  + scale_y_continuous( labels = comma , limits=c(0,500))
  + theme( panel.grid.major.x = element_blank(),
           panel.grid.minor = element_blank(),
           axis.text.x = element_blank(),
           #axis.text.x=element_text(angle=25, hjust=1),
           axis.text=element_text(color='black'))
)
#ggsave("./99.figures/20190826-pairwiseSNPs-Bp-CORRECTED.pdf", device = 'pdf', width = 3, height = 3.5, units = 'in', useDingbats=F)

# Subset Bb only
pairwise.bb <-subset(pairwise.snps,V1 == "Bbronchiseptica" | V1 == "Bb_complexI" | V1 == "Bb_complexIV")
#nrow(pairwise.bb)

(ggplot( data = pairwise.bb, aes(x=V1, y=V2))
  + geom_jitter(width = 0.2, size=1.2, color="red")
  + geom_boxplot(outlier.size = 0, outlier.stroke = 0, alpha=0, color='black')
  + theme_minimal(base_size = 14)
  + labs(x="Species", y="Pairwise SNP distance")
  + scale_y_continuous( labels = comma )#, limits=c(0,500))
  + theme( panel.grid.major.x = element_blank(),
           panel.grid.minor = element_blank(),
           axis.text.x = element_blank(),
           #axis.text.x=element_text(angle=25, hjust=1),
           axis.text=element_text(color='black'))
)
ggsave("./99.figures/20190823-pairwiseSNPs-Bb-CORRECTED.pdf", device = 'pdf', width = 3, height = 3.5, units = 'in', useDingbats=F)
