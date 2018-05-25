setwd("~/Documents/Bordetella_species/results/mauve/20180221-Bp/")

library(reshape)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(igraph)
library(scales)


edges0316 <- read.csv("./04.colinear-mcl/20180221-Bp-check-gap1500-cat.NR.network-20180316.csv", header = F)
net0316 <- network(edges0316, directed=F, matrix.type='edgelist', ignore.eval=F, names.eval=c("weights","symm"))
nodesize0316 <- read.csv("./04.colinear-mcl/20180221-Bp-check-gap1500-cat.NR.nodesize-20180316.csv", 
                         header=F,row.names=1)

set.edge.attribute(net0316, "lty", ifelse(net0316 %e% "weights" > 1, 2, 1))
set.edge.attribute(net0316, "color", ifelse(net0316 %e% "symm" > 1, "red", "black"))
set.vertex.attribute(net0316, "cluster", 
                     apply(data.frame(network.vertex.names(net0316)), 1, function(x) nodesize0316[x,])  )
set.vertex.attribute(net0316, "singles",
                     ifelse(grepl("Singleton",network.vertex.names(net0316)), "Singleton", "Cluster") )

ggnet2(net0316, label=T, size="cluster", size.min = 0, label.size=3, max_size = 10,
       mode = "fruchtermanreingold", #layout.par = list(cell.jitter = 0.5),
       #mode = "kamadakawai", #layout.par = list(evsel = "size"),
       color = ifelse(net0316 %v% "singles" == "Singleton", "blue", "gray"),
       edge.lty="lty", edge.color = "color")

## updated 20180426
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

netgraph2=ggnet2(net0426, label=F, size="cluster", size.min = 0, label.size=3, max_size = 20,
#netgraph2=ggnet2(net0426, label=F, size=0, size.min = 0, label.size=3, max_size = 20,
       mode = "fruchtermanreingold", layout.par = list(niter = 2000 ,ncell =1500 ),# cell.jitter = 0.25),
       #mode = "target", #layout.par = list(evsel = "size"),
       color = ifelse(net0426 %v% "singles" == "Singleton", "blue", "gray"),
       #shape = ifelse(net0426 %v% "singles" == "Singleton", 15, 19),
       edge.lty="lty", edge.color = "color")
netgraph2
ggsave("./99.figures/20180430-netgraph2.pdf", device = 'pdf', width = 8, height = 6, units = 'in', useDingbats=F)
net0426

netgraph
netgraph2.dat <- netgraph2$data
netgraph2.top5=head(netgraph2.dat[order(-netgraph2.dat$size),], n=8)
(netgraph2
  + geom_text(data=netgraph2.top5, aes(x=x,y=y,label = label)) #substr(label,9,1)))
)
ggsave("./99.figures/20180430-netgraph2-lables.pdf", device = 'pdf', width = 8, height = 6, units = 'in', useDingbats=F)






?gplot.layout
?edge_betweenness
###############################################
df = data.frame(cbind( network.vertex.names(net0316), 
       degree(net0316, gmode="graph", ignore.eval = T),
       get.vertex.attribute(net0316, attrname = "cluster") ) )
colnames(df) <- c("Cluser","Centrality","Size")
df[sort(-as.numeric(df$Size)),]
plot(get.vertex.attribute(net0316, attrname = "cluster"),  
     degree(net0316, gmode="graph", ignore.eval = T),
     xlab = "Cluster size", ylab = "Centrality" )

centralization(net0316, degree, mode="graph")
cugtest(dat = net0316, gcor, gmode="graph" )

## updated 20180426
?degree
?betweenness
?closeness
?centralization
df = data.frame(cbind( network.vertex.names(net0426), 
                       degree(net0426, gmode="graph", ignore.eval = T),
                       get.vertex.attribute(net0426, attrname = "cluster") ) )
centralization(net0426,g=10,closeness)
df = data.frame(cbind( network.vertex.names(net0426), 
                       evcent(net0426, gmode="graph",ignore.eval = T),
                       degree(net0426, gmode="graph", ignore.eval = T),
                       get.vertex.attribute(net0426, attrname = "cluster") ) )

colnames(df) <- c("Cluster","Centrality","Degree","Size")
df$Size <- as.numeric(as.character(df$Size))
df$Centrality <- as.numeric(as.character(df$Centrality))
df$Degree <- as.numeric(as.character(df$Degree))


top5=head(df[order(-df$Centrality),], n=4)
top4 = head(top5, n=4)

(ggplot( data = df, aes(x=Degree, y=Centrality))
  + geom_point(shape = 19, color="black")
  + theme_minimal()
#  + scale_x_continuous(limits=c(0,100))
  #+ scale_y_continuous(limits=c(0,15))
#  + labs(x="Cluster size", y="Centrality" )
#  + geom_point( data = top4, aes(x=Size, y=Centrality), shape = 1, color ="red", size=6)
  + geom_text( data = top5, aes(x=Degree, y=Centrality, label=Cluster),
               vjust=-1.2, size=3, angle = 0, hjust = 1) #-0.01)
#  + geom_segment(data = top5, aes(x = top5[3,3], y = top5[3,2], xend=top5[4,3], yend=top5[4,2]), lty=1, size=0.2)
#  + geom_segment(data = top5, aes(x = top5[1,3], y = top5[1,2], xend=top5[2,3], yend=top5[2,2]), lty=1, size=0.2)
#  + geom_segment(data = top5, aes(x = top5[3,3], y = top5[3,2], xend=top5[1,3], yend=top5[1,2]), lty=1, size=0.2)
#  + geom_segment(data = top5, aes(x = top5[4,3], y = top5[4,2], xend=top5[2,3], yend=top5[2,2]), lty=1, size=0.2)
  )
ggsave("./99.figures/20180430-centrality.pdf", device = 'pdf', width = 6, height = 6, units = 'in', useDingbats=F)



get.vertex.attribute(net0426, attrname = "cluster")
degree(net0426, gmode="graph", ignore.eval = T)
centralization(net0426, degree, mode="graph")
cugtest(dat = net0426, gcor, gmode="graph" )

###########################################
#symm vs asymm
sym0430 <- read.table("./04.colinear-mcl/20180425-Bp-check-gap1500-cat.NR.invertALL.txt", sep="\t", header = F, stringsAsFactors = F)
?read.table
(ggplot( data=sym0430, aes(x=V8, y=(V7/1000000), color=V11) )
  + geom_boxplot(outlier.size = 0, outlier.stroke = 0)  
  + geom_jitter(width = 0.15, size=2.5, alpha=0.6) #, aes(fill=V11))
  + labs(x="Inversion type", y="Inversion size (Mbp)" )
  + theme_classic()
#  + theme_minimal()
#  + theme_gray()
  + scale_y_continuous( labels = comma, breaks = pretty(sym0430$V7/1000000, n = 10))
  + theme( legend.position = c(0.1,0.8), 
           legend.background=element_rect(fill="white", size=0.5, color="black"),
           panel.grid.minor.y = element_blank(),
           panel.grid.major.x = element_blank())
)
ggsave("./99.figures/Inversion-size.pdf", device = 'pdf', width = 6, height = 6, units = 'in', useDingbats=F)



sym0430.filt = subset(sym0430, V11 == "ori" | V11 == "term")
head(sym0430.ratio)
sym.ratio = (as.numeric(as.character(sym0430.filt$V9))/as.numeric(as.character(sym0430.filt$V10)))
head(sym.ratio)
nrow(sym0430.filt)
sym0430.ratio = cbind(sym0430.filt, sym.ratio)
break1 = c(0, 250,500,750,1000,1250,1500) 
sym0430.filt2 = subset(sym0430.filt, sym.ratio < 3 & sym.ratio > 0.3 )
not.out = sym0430.ratio[!sym0430.ratio$sym.ratio %in% boxplot.stats(sym0430.ratio$sym.ratio)$out,]
nrow(not.out)

nrow(sym0430.filt2)
#sym0430.filt2 = subset(sym0430.filt, as.numeric( 3 > as.character(V9))/as.numeric(as.character(V10)) > 0.3 )

sym0430.fit = lm(as.numeric(as.character(V10)) ~ as.numeric(as.character(V9)), data = sym0430.filt2)
summary(sym0430.fit)
?geom_smooth

flip = cbind(as.numeric(as.character(sym0430.filt$V10)),as.numeric(as.character(sym0430.filt$V9)))
colnames(flip) = c('V9','V10')
head(flip)



just.coords = sym0430.filt[,c('V9','V10')]
wtf.all = rbind(just.coords,flip)
nrow()
(ggplot() # data=sym0430.filt, aes(x=as.numeric(as.character(V9))/1000, y=as.numeric(as.character(V10))/1000))#, color=V11) )
  #+ geom_point(size=2.0,alpha=0.6, aes(color = V11))
  + geom_point(data=sym0430.filt, aes(x=as.numeric(as.character(V9))/1000, 
                                      y=as.numeric(as.character(V10))/1000, 
                                      color=V11),
               size=2.0, alpha=0.6)
  + geom_smooth(method=lm, data=sym0430.filt,aes(x=as.numeric(as.character(V9))/1000, 
                                                 y=as.numeric(as.character(V10))/1000)
                )
  + geom_smooth(method=lm, data=wtf.df3, aes(x=as.numeric(as.character(X1))/1000, 
                                      y=as.numeric(as.character(X2))/1000),fullrange=F) 
  + geom_point(data=wtf.df3, mapping = aes(x=as.numeric(as.character(X1))/1000, 
                                 y=as.numeric(as.character(X2))/1000))
#  + theme_classic()
  + theme_minimal()
  + theme( legend.position = c(0.79,0.3), 
           legend.background=element_rect(fill="white", size=0.5, color="black"))
  + scale_x_continuous(labels = comma, breaks = break1)
  + scale_y_continuous(labels = comma, breaks = break1)
  + coord_cartesian(xlim=c(0,1500),ylim=c(0,1500))
  + labs(x="Right replichore size (kbp)", y="Left replichore size (kbp)" )
  + geom_abline(aes(intercept=0,slope=1),lty=2, size=0.4) 
  #+ geom_abline(aes(intercept=parameters[2]/1000, slope=parameters[1]))
  + geom_line(data = new.df, mapping = aes(x=X1/1000, y=lwr/1000),color="red",linetype='dashed')
  + geom_line(data = new.df, mapping = aes(x=X1/1000, y=upr/1000),color="red",linetype='dashed')
  #+ geom_linerange(data = new.df, mapping = aes(x=X1/1000,ymin=lwr/1000, ymax=upr/1000),color="red")
  
)

ggsave("./99.figures/Replichore-symm.pdf", device = 'pdf', width = 6, height = 6, units = 'in', useDingbats=F)


###################
wtf.df <- data.frame(cbind(as.numeric(as.character(sym0430.filt2$V9)),as.numeric(as.character(sym0430.filt2$V10))) )
class(wtf.df3)
head(wtf.df3)
nrow(wtf.df3)
wtf.df2 = cbind(wtf.df$X2,wtf.df$X1)
colnames(wtf.df2) = c('X1','X2')
wtf.df3 = rbind(wtf.df,wtf.df2)

linearRSS = function(parameters, data){
  predicted = parameters[1] * wtf.df3[,'X1'] + parameters[2]
  sum((predicted - wtf.df3[,'X2'])^2)
}
parameters = optim(par=c(1,1), fn = linearRSS, data = wtf.df3)$par
optized = optim(par=c(1,1), fn = linearRSS, data = wtf.df3)
optized
print(parameters)

sym.lm = lm( X2 ~ X1, data=wtf.df3)
summary(sym.lm)
newtest = data.frame(X1 = wtf.df$X1)
head(newtest)
newpred = predict(sym.lm, newdata = newtest, interval="pred",level=0.95)
new.df = cbind(newtest, newpred)
head(new.df)
head(newpred)

xval = seq(0,1500000, 50000)
yval = xval * parameters[1] + parameters[2]

plot(x=wtf.df$X1, y=wtf.df$X2)
lines(xval,yval)

(ggplot( data = wtf.df, aes(x=X1, y=X2))
  + geom_point()
  + coord_cartesian(xlim=c(0,1500000),ylim=c(0,1500000))
  + geom_abline(aes(intercept=parameters[2], slope=parameters[1]))
  + geom_abline(aes(intercept=0,slope=1),lty=2, size=0.4) 
)

head(sym0430.ratio)
(ggplot( data=sym0430.ratio, aes(x=V8, y=sym.ratio)) #, color=V11) )
  + geom_boxplot( outlier.size = 4, outlier.stroke = 0, outlier.color = "red")  
  + geom_jitter(width = 0.15, size=2.5, alpha=0.6) #, aes(fill=V11))
  #+ labs(x="Inversion type", y="Inversion size (Mbp)" )
  + theme_classic()
  #+ coord_cartesian(ylim=c(0,6))
  #+ ylim(0,3)
)

#############
J549.IS481 = read.table("./04.colinear-mcl/invert-predict/J549-IScollapsed.txt", sep="\t",header=F)
head(J549.IS481)
IS481.cent = (J549.IS481$V1 + J549.IS481$V2) / 2
head(IS481.cent)
to.pred = data.frame(X1 = IS481.cent)

J549.pred = predict(sym.lm, newdata = to.pred, interval="pred",level=0.95)
J549.ISpred = cbind(J549.IS481,IS481.cent,J549.pred)
head(J549.ISpred)

write.table(J549.ISpred,file="./04.colinear-mcl/invert-predict/J549-ISpredicted.txt",quote=F,sep="\t",row.names=F)
