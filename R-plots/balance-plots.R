# This script plots the density of observed replichore balance or single-inversion balance. Intended to be run interactively in Rstudio.

setwd("~/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/balance/")
library(ggplot2)
library(scales)

##### Replichore balance observed in 107 unique B. pertussis chromosome structures
log.df = read.table("./Bp-check-gap1500-mcl.NR.balance-log.txt",header=F, sep="\t")
(ggplot( data = log.df, aes(x=V2))
  + geom_density(color="black",fill="grey",alpha=0.4, size=0.75)
  #  + geom_histogram(binwidth = 0.04,color="blue",fill="blue",alpha=0.4)
  + theme_classic(base_size = 14)
  + coord_cartesian(ylim=c(0,6),xlim=c(0,1.3))
  + scale_x_continuous(expand = c(0,0), breaks=c(0,0.25,0.5,0.75,1,1.25))
  + scale_y_continuous(expand = c(0,0))
  + labs(x="Replichore balance |ln(ratio)|", y="Density" )
  + theme(axis.text = element_text(color='black'))
)

#ggsave("./20180803-Replichores.pdf", device = 'pdf', width = 3, height = 3, units = 'in', useDingbats=F)


##### Single-inversion balance observed in B. pertussis and predicted with the linear model
refs = read.table("../invert-predict/Bsp-combined-20180719.inverts.txt", header=F,sep="\t")
refs.symmetric = subset(refs, V11 == "ori" | V11 == "term")
ratio = abs(log(as.numeric(as.character(refs.symmetric$V9))/as.numeric(as.character(refs.symmetric$V10))))
refs.ratio = cbind(refs.symmetric, ratio)

#model prediction
J549.pred = read.table("../invert-predict/J549-ISall-inverts-20180719.txt", header=F, sep="\t")
J549.ratio = abs(log(as.numeric(as.character(J549.pred$V10))/as.numeric(as.character(J549.pred$V11))))
J549.pred = cbind(J549.pred, J549.ratio)

(ggplot()
  + geom_histogram( data = refs.ratio, mapping=aes(x=ratio, fill=V11,color=V11),
                    binwidth = 0.2,alpha=0.5,position='dodge')
  + theme_classic(base_size = 14)
  + coord_cartesian(ylim=c(0,12),xlim=c(0,3.55))
  + scale_x_continuous( breaks=seq(0,4.5,0.5))
  + scale_y_continuous(expand = c(0,0))
  + labs(x="Breakpoint balance |ln(ratio)|", y="Count or density" )
  + theme(axis.text = element_text(color='black'),
          legend.position = c(0.75,0.5))
  + geom_density( data = refs.ratio, mapping=aes(x=ratio, y=..scaled..*10), color='red', size=0.75, linetype=1 )
  + geom_density( data = J549.pred, mapping=aes(x=J549.ratio, y=..scaled..*10), color='black', size=0.75, linetype=3)
)

#ggsave("./20180822-symmetric-inverts.pdf", device = 'pdf', width = 3, height = 3, units = 'in', useDingbats=F)


##### Replichore balance vs colinear cluster size
bal = read.table("./cluster-size-v-balance.csv", header=F,sep=",")
top5=head(bal[order(-bal$V2),], n=6)

(ggplot() #data = bal, aes(x=V3,y=V2))
  + geom_point(data = subset(bal, V2 > 1), aes(x=V3,y=V2),
               shape = 19, color="black", alpha=0.5,size=2)
  + geom_point(data = subset(bal, V2 == 1), aes(x=V3,y=V2),
               shape = 19, color="red", alpha=0.5,size=2)
  + theme_classic(base_size = 10)
  + theme(axis.text = element_text(color='black'))
  + scale_x_continuous(limits=c(0,1.25),expand = c(0,0))
  + scale_y_continuous(limits=c(0,100),expand = c(0,0))
  + labs(y="Cluster size", x="Replichore balance |ln(ratio)|" )
  + geom_text( data = top5, aes(x=V3, y=V2, label=V1),
               vjust=2.5, size=3, angle = 0, hjust = 0.01)
)

#ggsave("./20180814-cluster-v-balance.pdf", device = 'pdf', width = 4, height = 4, units = 'in', useDingbats=F)

bho = read.table("/home/yrh8/Documents/Bordetella_species/results/mauve/Bsp-combined-20190903/Bho_cluster-size-v-balance.csv", header=F,sep=",")
topho=head(bho[order(-bho$V2),], n=4)

(ggplot() #data = bal, aes(x=V3,y=V2))
  + geom_point(data = subset(bho, V2 > 1), aes(x=V3,y=V2),
               shape = 19, color="black", alpha=0.5,size=2)
  + geom_point(data = subset(bho, V2 == 1), aes(x=V3,y=V2),
               shape = 19, color="red", alpha=0.5,size=2)
  + theme_classic(base_size = 10)
  + theme(axis.text = element_text(color='black'))
  + scale_x_continuous(limits=c(0,0.8),expand = c(0,0))
  + scale_y_continuous(limits=c(0,50),expand = c(0,0))
  + labs(y="Cluster size", x="Replichore balance |ln(ratio)|" )
  + geom_text( data = topho, aes(x=V3, y=V2, label=V1),
               vjust=2.5, size=3, angle = 0, hjust = 0.01)
)
#ggsave("/home/yrh8/Documents/Bordetella_species/results/mauve/Bsp-combined-20190903/20190904-Bho-cluster-v-balance.pdf", device = 'pdf', width = 4, height = 4, units = 'in', useDingbats=F)


bpp = read.table("/home/yrh8/Documents/Bordetella_species/results/mauve/Bsp-combined-20190903/Bpp_cluster-size-v-balance.csv", header=F,sep=",")
toppp=head(bpp[order(-bpp$V2),], n=2)

(ggplot() #data = bal, aes(x=V3,y=V2))
  + geom_point(data = subset(bpp, V2 > 1), aes(x=V3,y=V2),
               shape = 19, color="black", alpha=0.5,size=2)
  + geom_point(data = subset(bpp, V2 == 1), aes(x=V3,y=V2),
               shape = 19, color="red", alpha=0.5,size=2)
  + theme_classic(base_size = 10)
  + theme(axis.text = element_text(color='black'))
  + scale_x_continuous(limits=c(0,0.8),expand = c(0,0))
  + scale_y_continuous(limits=c(0,50),expand = c(0,0))
  + labs(y="Cluster size", x="Replichore balance |ln(ratio)|" )
  + geom_text( data = toppp, aes(x=V3, y=V2, label=V1),
               vjust=2.5, size=3, angle = 0, hjust = 0.01)
)
#ggsave("/home/yrh8/Documents/Bordetella_species/results/mauve/Bsp-combined-20190903/20190904-Bpp-cluster-v-balance.pdf", device = 'pdf', width = 4, height = 4, units = 'in', useDingbats=F)
