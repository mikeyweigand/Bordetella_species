# This script calculates a rare faction curve from the colinear structure clusters observed in years 2000 - 2010.

setwd("/home/yrh8/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/rarefaction")
library(vegan)

Bp.y2k10 = read.table("20180425-Bp-check-gap1500-mcl-rarefaction-years2k10.txt", sep="\t", header=T)
row.names(Bp.y2k10) = c("Structures")
rarecurve(Bp.y2k10, step=2, label=F, ylab="Unique Structures")

rmax = min(rowSums(Bp.y2k10))
N <- seq(2, rmax, by=2)
S <- rarefy(Bp.y2k10, N, se = TRUE)

rc = data.frame(rbind(N, S[1,], S[2,]))
rownames(rc) <- c("N","C","SE")

(ggplot()
  + labs(x="Sample size", y="Unique structures")
  + theme_classic(base_size = 16)
  + theme(axis.text = element_text(color='black'))
  + coord_cartesian(xlim=c(0,400),ylim=c(0,80))
  + geom_ribbon(data= data.frame(t(rc)),aes(ymin=(C - SE*2), ymax=(C + SE*2), x=N),alpha=0.25 )
  + geom_line(data = data.frame(t(rc)),aes(x=N,y=C),size=2)
)

#ggsave("./20180425-Bp-check-gap1500-mcl-rarefaction-years2k10.pdf", device = 'pdf', width = 6, height = 4, units = 'in', useDingbats=F)
