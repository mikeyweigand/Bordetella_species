# This script calculates the single inversion and deletion/insertion network and related summary statistics as well as draw the relvant plots.

setwd("~/Documents/Bordetella_species/results/mauve/20180221-Bp/")

detach("package:igraph", unload=T)
detach("package:circlize", unload=T)
library(reshape)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(scales)


## import network edge list and calculate/plot network
edges0426 <- read.csv("./04.colinear-mcl/20180425-Bp-check-gap1500-cat.NR.network.csv", header = F)
net0426 <- network(edges0426, directed=F, matrix.type='edgelist', ignore.eval=F, names.eval=c("weights","symm"))
nodesize0426 <- read.csv("./04.colinear-mcl/20180425-Bp-check-gap1500-cat.NR.nodesize.csv", 
                         header=F,row.names=1)

set.edge.attribute(net0426, "lty", ifelse(net0426 %e% "weights" > 1, 2, 1))
set.edge.attribute(net0426, "color", ifelse(net0426 %e% "symm" > 1, "red", "black"))
set.vertex.attribute(net0426, "cluster", 
                     apply(data.frame(network.vertex.names(net0426)), 1, function(x) nodesize0426[x,])  )
set.vertex.attribute(net0426, "singles",
                     ifelse(grepl("Singleton",network.vertex.names(net0426)), "Singleton", "Cluster") )

netgraph2=ggnet2(net0426, label=T, size="cluster", size.min = 0, label.size=3, max_size = 20,
       mode = "fruchtermanreingold", layout.par = list(niter = 2000 ,ncell =1500 ),# cell.jitter = 0.25),
       color = ifelse(net0426 %v% "singles" == "Singleton", "blue", "gray"),
       edge.lty="lty", edge.color = "color")
netgraph2

#ggsave("./99.figures/20180430-netgraph2.pdf", device = 'pdf', width = 8, height = 6, units = 'in', useDingbats=F)

# netgraph2.dat <- netgraph2$data
# netgraph2.top5=head(netgraph2.dat[order(-netgraph2.dat$size),], n=8)
# (netgraph2
#   + geom_text(data=netgraph2.top5, aes(x=x,y=y,label = label)) #substr(label,9,1)))
# )
#ggsave("./99.figures/20180430-netgraph2-lables.pdf", device = 'pdf', width = 8, height = 6, units = 'in', useDingbats=F)



###### Calculate node centrality and plot relative to cluster size

df = data.frame(cbind( network.vertex.names(net0426), 
                       evcent(net0426, gmode="graph",ignore.eval = T),
                       degree(net0426, gmode="graph", ignore.eval = T),
                       get.vertex.attribute(net0426, attrname = "cluster") ) )

colnames(df) <- c("Cluster","Centrality","Degree","Size")
df$Size <- as.numeric(as.character(df$Size))
df$Centrality <- as.numeric(as.character(df$Centrality))
df$Degree <- as.numeric(as.character(df$Degree))
#write.table(x = df, file="./04.colinear-mcl/20180425-Bp-check-gap1500-cat.NR.network-output.tsv", sep="\t",quote = F)
top5=head(df[order(-df$Degree),], n=7)
(ggplot( data = df, aes(x=Degree, y=Size))
  + geom_jitter(shape = 19, color="black", alpha=0.5, width=0.2,height=0.2, size=2.5)
  + theme_classic(base_size = 10)
  + theme(axis.text = element_text(color='black'))
  + scale_y_continuous(limits=c(0,100),expand = c(0,0))
  + scale_x_continuous(limits=c(0,15),expand = c(0,0))
  + labs(y="Cluster size", x="Degree centrality" )
  + geom_text( data = top5, aes(y=Size, x=Degree, label=Cluster),
              vjust=-1.2, size=3, angle = 0, hjust = 1) #-0.01)
)
#ggsave("./99.figures/20180906-centrality.pdf", device = 'pdf', width = 4, height = 4, units = 'in', useDingbats=F)



get.vertex.attribute(net0426, attrname = "cluster")
degree(net0426, gmode="graph", ignore.eval = T)
centralization(net0426, degree, mode="graph")
cugtest(dat = net0426, gcor, gmode="graph" )

###### Plot inversion sizes
#symm vs asymm
sym0430 <- read.table("./04.colinear-mcl/20180425-Bp-check-gap1500-cat.NR.invertALL.txt", sep="\t", header = F, stringsAsFactors = F)
?read.table
(ggplot( data=sym0430, aes(x=V8, y=(V7/1000000), color=V11) )
  + geom_boxplot(outlier.size = 0, outlier.stroke = 0)  
  + geom_jitter(width = 0.15, size=2.5, alpha=0.6) #, aes(fill=V11))
  + labs(x="Inversion type", y="Inversion size (Mbp)" )
  + theme_classic()
  + scale_y_continuous( labels = comma, breaks = pretty(sym0430$V7/1000000, n = 10))
  + theme( legend.position = c(0.1,0.8), 
           legend.background=element_rect(fill="white", size=0.5, color="black"),
           panel.grid.minor.y = element_blank(),
           panel.grid.major.x = element_blank())
)
#ggsave("./99.figures/Inversion-size.pdf", device = 'pdf', width = 6, height = 6, units = 'in', useDingbats=F)
