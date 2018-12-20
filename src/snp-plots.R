setwd("~/Documents/Bordetella_species/results/mauve/20180221-Bp/")

library(ggplot2)
library(scales)
library(reshape)

#clust.dist <-  read.table("./04.colinear-mcl/test-cluster-dist.tsv", sep="\t",header=F)
clust.dist <- read.table("./04.colinear-mcl/20180425-Bp-check-gap1500-mcl-cluster-dist.tsv", sep="\t", header=F)
#clust.size <-  read.table("./04.colinear-mcl/test-cluster-counts.tsv", sep="\t",header=F)
#clust.dist2k10 <- read.table("./04.colinear-mcl/20180425-Bp-check-gap1500-mcl-cluster-dist2k10.tsv", sep="\t",header=F)

#head(clust.dist)
#head(clust.size)
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
ggsave("/media/cdc_documents/Bordetella/Bordetella_species/20180723-cluster-SNPs.pdf", device = 'pdf', width = 6.25, height = 3, units = 'in', useDingbats=F)


# head(clust.wtf)
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
#   
# )



# head(clust.dist2k10)
# (ggplot( data = clust.dist2k10, aes(x=V4, y=V3))
#   + geom_jitter(width = 0.2, size=0.5, color="red")
#   + geom_boxplot(outlier.size = 0, outlier.stroke = 0, alpha=0)
#   + scale_x_discrete(limits=first12)
#   + scale_y_continuous(limits=c(0,35))
#   + theme_minimal()
#   + theme(axis.text.x=element_text(angle=45, hjust=1))
#   + labs(x="Colinear cluster", y="Pairwise SNP distance")
# )
#ggsave("./99.figures/20180430-cluster-SNPs2k10.pdf", device = 'pdf', width = 8, height = 6, units = 'in', useDingbats=F)




######
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

# (ggplot( data = clust.years, aes(x=V2, y=V3))
#   + geom_jitter(width = 0.2, height=0.2, size=0.5, color="red")
#   + scale_y_discrete(limits=rev)
#   + scale_x_continuous(breaks = y2)
#   + theme(axis.text.x=element_text(angle=45, hjust=1))
#   + labs(y="Colinear cluster", x="Year of Isolation")
# )

(ggplot( data = clust.years2kp, aes(x=as.factor(V2), y=V3))
  + geom_jitter(width = 0.1, height=0.35, size=1, color="red")
  #+ scale_y_discrete(limits=first12s)
  + scale_y_discrete(limits=sort(unique(clust.years$V3),decreasing=T) )
    #+ scale_x_continuous(breaks = y.years)
  + theme_minimal(base_size = 16)
  + theme(axis.text.x=element_text(angle=45, hjust=1), 
          panel.grid.minor.x = element_blank(),
          axis.text=element_text(color='black'))
  + labs(y="Colinear cluster", x="Year of Isolation")

  )
ggsave("/media/cdc_documents/Bordetella/Bordetella_species/20180723-cluster-years.pdf", device = 'pdf', width = 6.25, height = 4.5, units = 'in', useDingbats=F)
#ggsave("/media/cdc_documents/Abstracts/ASM2018/poster-figs/20180430-cluster-years-02.pdf", device = 'pdf', width = 6, height = 3.3, units = 'in', useDingbats=F)

head(clust.dist)
head(wtf)

wtf=clust.years2kp[clust.years2kp$V3 %in% first12s,]
?geom_density
#(ggplot( data = wtf, aes(x=V2, fill=V3, color=V3, alpha=0.3))
#(ggplot( data = clust.years2kp, aes(x=V2, fill=V3, color=V3, alpha=0.3))
#  + geom_density(trim=F, kernel="g", adjust=0.75)
  #+ geom_freqpoly(binwidth=1)
  #+ scale_y_continuous(limits=c(0,1))
#   + scale_x_continuous(limits = c(2000,2016))
#   + theme_classic()
#   + facet_grid(V3 ~ .)
#   + theme(strip.text.x=element_blank(), panel.background=element_blank(),
#           panel.spacing.y=unit(-2.5, "cm"), panel.grid=element_blank(),
#           axis.ticks.y=element_blank(), axis.text.y=element_blank())
#   + labs(y="Colinear cluster", x="Year of Isolation")
# 
# )

##### Bordetella species comparison

pairwise.snps <- read.table("../../kSNP/20180501-SNP-dist/20180501-pairwiseSNPs-down-plusComplexes.tsv",sep="\t",header=F);

pairwise.snps$V1 <- factor(pairwise.snps$V1, levels = c("Bbronchiseptica","Bb_complexI","Bb_complexIV","Bholmesii","Bparapertussis","Bpertussis"), ordered=T)
pairwise.snps$V1
(ggplot( data = pairwise.snps, aes(x=V1, y=V2))
 + geom_jitter(width = 0.2, size=1.5, color="red")
 + geom_boxplot(outlier.size = 0, outlier.stroke = 0, alpha=0)
 + theme_minimal(base_size = 16)
 + theme(axis.text.x=element_text(angle=25, hjust=1), axis.text=element_text(color='black'))
 + labs(x="Species", y="Pairwise SNP distance")
 + scale_y_continuous( labels = comma)
)
ggsave("../../kSNP/20180501-SNP-dist/20180726-pairwiseSNPs-down.pdf", device = 'pdf', width = 6, height = 5, units = 'in', useDingbats=F)


pairwise.pert <-subset(pairwise.snps,V1 == "Bparapertussis" | V1 == "Bpertussis")
#pairwise.pert <- read.table("../../kSNP/20180501-SNP-dist/20180501-pairwiseSNPs-Bp-short.tsv",sep="\t",header=F);
nrow(pairwise.pert)

(ggplot( data = pairwise.pert, aes(x=V1, y=V2))
  + geom_jitter(width = 0.2, size=1.5, color="red")
  + geom_boxplot(outlier.size = 0, outlier.stroke = 0, alpha=0)
  + theme_minimal(base_size = 16)
  + labs(x="Species", y="Pairwise SNP distance")
  + scale_y_continuous( labels = comma )#, limits=c(0,500))
  + theme( panel.grid.major.x = element_blank(),
           axis.text.x=element_text(angle=25, hjust=1), 
           axis.text=element_text(color='black'))
)
ggsave("../../kSNP/20180501-SNP-dist/20180726-pairwiseSNPs-Bp2.pdf", device = 'pdf', width = 4, height = 3.5, units = 'in', useDingbats=F)

