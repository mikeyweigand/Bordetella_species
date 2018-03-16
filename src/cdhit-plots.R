setwd("~/Documents/Bordetella_species/results/cd-hit/")
library(ggplot2)
library(scales)
library(reshape)

### Representative subset ###
cdhit.Breps = read.table("./20180222/table/Reps.clster.hist.txt", header=T, sep="\t")
cdhit.Breps.melted <- melt(cdhit.Breps, id = "Tag")

(ggplot( data = cdhit.Breps.melted, aes(x=Tag, y=value, color = variable))
  #+ geom_point(size=0.25)
  + geom_line()
  + coord_cartesian(ylim=c(0,10),xlim=c(0,300))
)
