setwd("~/Documents/Bordetella_species/results/jellyfish/")
library(ggplot2)
library(scales)
library(reshape)

### B. bronch ###
# hist.Bb = read.table("./20180222/hist/Bb-15mer.hist", header=T, sep="\t")
# hist.Bb2 <- cbind(hist.Bb, rowMeans(hist.Bb[-1]))
# colnames(hist.Bb2)[dim(hist.Bb2)[2]] = "Mean"
# hist.Bb2.melted <- melt(hist.Bb2, id = "Tag")
# head(hist.Bb2.melted)
# (ggplot( data = hist.Bb2.melted, aes(x=Tag, y=value, color = variable))
#   + geom_point(size=0.5)
#   + geom_line()
#   + coord_cartesian(ylim=c(0,1000),xlim=c(0,80))
# )


### B. pertussis ###
# hist.Bp = read.table("./20180222/hist/Bp-15mer.hist", header=T, sep="\t")
# hist.Bp2 <- cbind(hist.Bp, rowMeans(hist.Bp[-1]))
# colnames(hist.Bp2)[dim(hist.Bp2)[2]] = "Mean"
# hist.Bp2.melted <- melt(hist.Bp2, id = "Tag")
# 
# (ggplot( data = hist.Bp2.melted, aes(x=Tag, y=value, color = variable))
#   + geom_point(size=0.5)
#   + geom_line()
#   + coord_cartesian(ylim=c(0,750),xlim=c(110,140))
# )
# (ggplot( hist.Bp2, aes(x=Tag, y=Bp_H559.15mer))
#   + geom_line()
#   + geom_point(size = 0.5)
#   + scale_y_continuous(label=comma,name = "Average Number of 15-mers")
#   + scale_x_continuous(name = "Frequency")
#   + coord_cartesian(ylim=c(0,1000),xlim=c(100,150))
# )

### B. parapertussis ###
# hist.Bpp = read.table("./20180222/hist/Bpp-15mer.hist", header=T, sep="\t")
# hist.Bpp2 <- cbind(hist.Bpp, rowMeans(hist.Bpp[-1]))
# colnames(hist.Bpp2)[dim(hist.Bpp2)[2]] = "Mean"
# hist.Bpp2.melted <- melt(hist.Bpp2, id = "Tag")
# 
# (ggplot( data = hist.Bpp2.melted, aes(x=Tag, y=value, color = variable))
#   + geom_point(size=0.5)
#   + geom_line()
#   + coord_cartesian(ylim=c(0,1500),xlim=c(0,80))
# )

### B. holmesii ###
# hist.Bho = read.table("./20180222/hist/Bho-15mer.hist", header=T, sep="\t")
# hist.Bho2 <- cbind(hist.Bho, rowMeans(hist.Bho[-1]))
# colnames(hist.Bho2)[dim(hist.Bho2)[2]] = "Mean"
# hist.Bho2.melted <- melt(hist.Bho2, id = "Tag")
# 
# (ggplot( data = hist.Bho2.melted, aes(x=Tag, y=value, color = variable))
#   + geom_point(size=0.5)
#   + geom_line()
#   + coord_cartesian(ylim=c(0,1500),xlim=c(0,50))
# )

### B. hinzii ###
# hist.Bhi = read.table("./20180222/hist/Bhi-15mer.hist", header=T, sep="\t")
# hist.Bhi2 <- cbind(hist.Bhi, rowMeans(hist.Bhi[-1]))
# colnames(hist.Bhi2)[dim(hist.Bhi2)[2]] = "Mean"
# hist.Bhi2.melted <- melt(hist.Bhi2, id = "Tag")
# 
# (ggplot( data = hist.Bhi2.melted, aes(x=Tag, y=value, color = variable))
#   + geom_point(size=0.5)
#   + geom_line()
#   + coord_cartesian(ylim=c(0,1000),xlim=c(0,30))
# )

### B. avium ###
# hist.Ba = read.table("./20180222/hist/Ba-15mer.hist", header=T, sep="\t")
# hist.Ba2 <- cbind(hist.Ba, rowMeans(hist.Ba[-1]))
# colnames(hist.Ba2)[dim(hist.Ba2)[2]] = "Mean"
# hist.Ba2.melted <- melt(hist.Ba2, id = "Tag")
# 
# (ggplot( data = hist.Ba2.melted, aes(x=Tag, y=value, color = variable))
#   + geom_point(size=0.5)
#   + geom_line()
#   + coord_cartesian(ylim=c(0,700),xlim=c(0,60))
# )

### B. trematum ###
# hist.Bt = read.table("./20180222/hist/Bt-15mer.hist", header=T, sep="\t")
# hist.Bt2 <- cbind(hist.Bt, rowMeans(hist.Bt[-1]))
# colnames(hist.Bt2)[dim(hist.Bt2)[2]] = "Mean"
# hist.Bt2.melted <- melt(hist.Bt2, id = "Tag")
# 
# (ggplot( data = hist.Bt2.melted, aes(x=Tag, y=value, color = variable))
#   + geom_point(size=0.5)
#   + geom_line()
#   + coord_cartesian(ylim=c(0,100),xlim=c(0,200))
# )

### B. sp ###
# hist.Bsp = read.table("./20180222/hist/Bsp-15mer.hist", header=T, sep="\t")
# hist.Bsp2 <- cbind(hist.Bsp, rowMeans(hist.Bsp[-1]))
# colnames(hist.Bsp2)[dim(hist.Bsp2)[2]] = "Mean"
# hist.Bsp2.melted <- melt(hist.Bsp2, id = "Tag")
# 
# (ggplot( data = hist.Bsp2.melted, aes(x=Tag, y=value, color = variable))
#   + geom_point(size=0.5)
#   + geom_line()
#   + coord_cartesian(ylim=c(0,1000),xlim=c(0,60))
# )
# 
# (ggplot( hist.Bhi2, aes(x=Tag, y=Bhi_H720.15mer))
#   + geom_line()
#   + geom_point(size = 0.5)
#   + scale_y_continuous(label=comma,name = "Average Number of 15-mers")
#   + scale_x_continuous(name = "Frequency")
#   + coord_cartesian(ylim=c(0,1500),xlim=c(0,50))
# )

### Representative subset ###
#hist.Breps = read.table("./20180222/hist/Reps-15mer.hist", header=T, sep="\t")
hist.Breps = read.table("./20180222/hist/Reps-15mer.20181129.hist", header=T, sep="\t")
hist.Breps.melted <- melt(hist.Breps, id = "Tag")
head(hist.Breps.melted)
(ggplot( data = subset(hist.Breps.melted, Tag > 1 ), aes(x=Tag, y=value, color = variable))
  #+ geom_col(position="dodge") #size=0.25)
  + geom_line()
  + coord_cartesian(ylim=c(0,1550),xlim=c(0,27))
  + theme_classic(base_size = 14)
  + scale_x_continuous(expand = c(0,0))
  + scale_y_continuous(expand = c(0,0))
  + theme( legend.position = c(0.5,0.5),
           legend.background=element_rect(fill="white", size=0.5, color="black"),
           legend.text=element_text(size=6),
           legend.title=element_text(size=6),
           legend.direction = 'horizontal',
           axis.text = element_text(color='black')
  )
  + labs(x="Copy number", y="15-mer count" )
)
ggsave("./20180222/99.figures/20181129-Reps-15mer.pdf", device = 'pdf', width = 6, height = 3, units = 'in', useDingbats=F)

hist.bp.bho = subset(hist.Breps.melted, subset = variable %in% c('Bho_F615.15mer','Bp_H627.15mer'))
head(hist.bp.bho)
(ggplot( data = subset(hist.bp.bho, Tag > 1 ), aes(x=Tag, y=value, color = variable))
  #+ geom_col(position="dodge") #size=0.25)
  + geom_line()
  + coord_cartesian(ylim=c(0,1550),xlim=c(0,155))
  + theme_classic(base_size = 14)
  + scale_x_continuous(expand = c(0,0))
  + scale_y_continuous(expand = c(0,0))
  + theme( legend.position = c(0.75,0.6),
           legend.background=element_rect(fill="white", size=0.5, color="black"),
           legend.text=element_text(size=7),
           axis.text = element_text(color='black')
  )
  + labs(x="Copy number", y="15-mer count" )
)
ggsave("./20180222/99.figures/20180803-Reps-15mer-Bp-Bho.pdf", device = 'pdf', width = 3, height = 3, units = 'in', useDingbats=F)



# (ggplot( hist.Breps, aes(x=Tag, y=Bpetrii_DSM12804.15mer))
#   + geom_line()
#   + geom_point(size = 0.5)
#   + scale_y_continuous(label=comma,name = "Average Number of 15-mers")
#   + scale_x_continuous(name = "Frequency")
#   + coord_cartesian(ylim=c(0,1500),xlim=c(0,50))
# )
