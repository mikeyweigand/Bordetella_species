setwd("~/Documents/Bordetella_species/results/mauve/20180221-Bp/")

library(GGally)
library(network)
library(sna)
library(ggplot2)

edges0316 <- read.csv("./04.colinear-mcl/20180221-Bp-check-gap1500-cat.NR.network-20180316.csv", header = F)
net0316 <- network(edges0316, directed=F, matrix.type='edgelist', ignore.eval=F, names.eval=c("weights","symm"))
nodesize0316 <- read.csv("./04.colinear-mcl/20180221-Bp-check-gap1500-cat.NR.nodesize-20180316.csv", header=F,row.names=1)

set.edge.attribute(net0316, "lty", ifelse(net0316 %e% "weights" > 1, 2, 1))
set.edge.attribute(net0316, "color", ifelse(net0316 %e% "symm" > 1, "red", "black"))
set.vertex.attribute(net0316, "cluster", which(network.vertex.names(net0316) %in% nodesize0316$V1))
ggnet2(net0316, label=T, size="cluster", size.min = 0, label.size=3, edge.lty="lty", edge.color = "color", legend.size = 0)

length(network.vertex.names(net0316))

which(network.vertex.names(net0316) %in% nodesize0316[,1])
nodesize0316["Cluster-1_H540",]
grep(network.vertex.names(net0316), nodesize0316$V1)