setwd("~/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/balance/")
library(ggplot2)
library(scales)

# balance.df = read.table("~/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/balance/Bp-check-gap1500-mcl.NR.balance.txt",header=F, sep="\t")
# head(balance.df)

log.df = read.table("~/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/balance/Bp-check-gap1500-mcl.NR.balance-log.txt",header=F, sep="\t")
head(log.df)
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

#ggsave("~/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/balance/20180803-Replichores.pdf", device = 'pdf', width = 3, height = 3, units = 'in', useDingbats=F)


##############

refs = read.table("/home/yrh8/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/invert-predict/Bsp-combined-20180719.inverts.txt", header=F,sep="\t")
head(refs)
refs.symmetric = subset(refs, V11 == "ori" | V11 == "term")
ratio = abs(log(as.numeric(as.character(refs.symmetric$V9))/as.numeric(as.character(refs.symmetric$V10))))
refs.ratio = cbind(refs.symmetric, ratio)
head(refs.ratio)

J549.pred = read.table("/home/yrh8/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/invert-predict/J549-ISall-inverts-20180719.txt", header=F, sep="\t")
J549.ratio = abs(log(as.numeric(as.character(J549.pred$V10))/as.numeric(as.character(J549.pred$V11))))
head(J549.pred)
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

#ggsave("~/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/balance/20180822-symmetric-inverts.pdf", device = 'pdf', width = 3, height = 3, units = 'in', useDingbats=F)


###################

bal = read.table("~/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/balance/cluster-size-v-balance.csv", header=F,sep=",")
head(bal)

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

#ggsave("~/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/balance/20180814-cluster-v-balance.pdf", device = 'pdf', width = 4, height = 4, units = 'in', useDingbats=F)


