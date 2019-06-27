# This script draws a coorelation plot from Tempest and a nice tree from the BEAST output

setwd("~/Documents/Bordetella_species/results/beast/")
library("ape")
library("ggtree")

beast.tree <- read.beast("./C734-core-ALL_drop75cov_plusYear.HKY-gamma4-strict.mcc.common.tre")
clusters <- read.table("clusters-01.tsv",sep="\t",header=F)
clusters.sub <- subset(clusters, V2 %in% c("Cluster-01","Cluster-02","Cluster-03","Cluster-04","Cluster-05","Cluster-06","Cluster-07","Cluster-08","Cluster-10","Singleton"))
colnames(clusters.sub) <- c("Isolate","Cluster")

beast <- (ggtree(beast.tree, right=T, ladderize=T, mrsd="2016-01-01",size=.35)
  + theme_tree2()
  + scale_x_continuous(breaks=c(1980, 1985,1990, 1995, 2000, 2005, 2010, 2015), minor_breaks=seq(1980, 2016, 1), limits = c(1980,2016))
  + theme(panel.grid.major   = element_line(color="black", size=.2),
          panel.grid.minor   = element_line(color="grey", size=.2, linetype=2),
          legend.text=element_text(size=7),
          legend.title=element_text(size=10),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) 
  + labs(x="Year")
)

beast <- beast %<+% clusters.sub + geom_tippoint(aes(color=Cluster),size=0.75)

(beast + theme(legend.position="right")
  + geom_text2(aes(subset=(node %in% c(730,728,439,637)), label=(round((2016-height),digits=0))), hjust=1.1, vjust=1.7, size=3, color="red")
  + geom_point2(aes(subset=(node %in% c(730,728,439,637))),color="red",shape="*",size=7)
)
#ggsave("./C734-core-ALL_drop75cov_plusYear.HKY-gamma4-strict.mcc.common.ggtree01-20190423.pdf", width=6, height=7, units="in",useDingbats=F,device='pdf')

# extract a few internal node dates of interest
MRCA(beast.tree, tip=c('H348_2010','I735_2013')) #cluster-01 vs cluster-07
MRCA(beast.tree, tip=c('I113_2012','H348_2010')) #cluster-07 vs cluster-04
MRCA(beast.tree, tip=c('I113_2012','C734_2000')) #fimH1 vs fimH2
MRCA(beast.tree, tip=c('F687_2008','H559_2010')) #cluster-4 vs cluster-5


############################
# Plots the tempest coorelation plot
tempest <- read.table("./C734-core-ALL_drop75cov_plusYear.rtt.txt",sep="\t",header=T,
                      colClasses = c("character", "numeric","numeric","numeric") )

(ggplot( data = tempest, aes(x=date,y=distance))
  + geom_jitter(fill="light grey", alpha=0.5)
  + geom_smooth(method=lm, formula = y ~ x, fullrange=T, se=F)
  + geom_vline(xintercept=1983, size=0,color="white") 
  + geom_vline(xintercept=2020, size=0,color="white")
  + theme_classic(base_size = 14)
  + coord_cartesian(xlim=c(1980,2020), ylim=c(0,0.025))
  + scale_y_continuous(expand = c(0,0))
  + labs(x="Year", y="Root-to-tip divergence" )
  + theme(axis.text = element_text(color='black'))
)

#ggsave("./C734-core-ALL_drop75cov_plusYear.rtt.pdf", device = 'pdf', width = 3, height = 3, units = 'in', useDingbats=F)

