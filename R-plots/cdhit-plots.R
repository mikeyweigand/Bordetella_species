# This script takes the outputs from cdhit and draws bar plots of gene copy number abundance. Intended to be run interactively in Rstudio.

setwd("~/Documents/Bordetella_species/results/cd-hit/")
library(ggplot2)
library(scales)
library(reshape)

### Representative subset ###
cdhit.Breps = read.table("./20180222/table/Reps.clster.hist.txt", header=T, sep="\t")
cdhit.Breps = read.table("./20180222/table/Reps.clster.hist.20181130.txt", header=T, sep="\t")

cdhit.Breps.melted <- melt(cdhit.Breps, id = "Tag")
#head(cdhit.Breps)
head(cdhit.Breps.melted)
subset(cdhit.Breps.melted, value != 0 )
(ggplot( data = subset(subset(cdhit.Breps.melted, value != 0), Tag > 1), aes(x=Tag, y=value, color = variable, fill=variable))
  #+ geom_bar(position="dodge") #size=0.25)
  + geom_col(position="dodge")
  + coord_cartesian(ylim=c(0,22),xlim=c(0,26))
  + theme_classic(base_size = 14)
  + scale_x_continuous(expand = c(0,0))
  + scale_y_continuous(expand = c(0,0))
  + theme( legend.position = c(0.5,0.6),
           legend.background=element_rect(fill="white", size=0.5, color="black"),
           legend.text=element_text(size=6),
           legend.title=element_text(size=6),
           legend.direction = 'horizontal',
           axis.text = element_text(color='black')
  )
  + labs(x="Copy number", y="Gene count" )
)
#ggsave("./20180222/99.figures/20181130-Reps-clster-hist.pdf", device = 'pdf', width = 6, height = 3, units = 'in', useDingbats=F)
#ggsave("./20180222/99.figures/Reps-clster-hist02.pdf", device = 'pdf', width = 8, height = 6, units = 'in', useDingbats=F)

bp.bho = subset(cdhit.Breps.melted, subset = variable %in% c('F615_CP018899','H627_CP010962'))
head(bp.bho)
(ggplot( data = subset(subset(bp.bho, value != 0), Tag > 1), aes(x=Tag, y=value, color = variable, fill=variable))
  #+ geom_bar(position="dodge") #size=0.25)
  + geom_col(position="dodge")
  + coord_cartesian(ylim=c(0,2.2),xlim=c(26,255))
  + theme_classic(base_size = 14)
  #+ scale_color_brewer(palette = "Set1")
  + scale_x_continuous(expand = c(0,0))
  + scale_y_continuous(expand = c(0,0))
  + theme( legend.position = c(0.7,0.6),
           legend.background=element_rect(fill="white", size=0.5, color="black"),
           axis.text = element_text(color='black')
  )
  + labs(x="Copy number", y="Gene count" )
)
#ggsave("./20180222/99.figures/20180803-Reps-clster-hist-Bp-Bho.pdf", device = 'pdf', width = 3, height = 3, units = 'in', useDingbats=F)
