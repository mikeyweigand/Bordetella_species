setwd("~/Documents/Bordetella_species/results/cd-hit/")
library(ggplot2)
library(scales)
library(reshape)

### Representative subset ###
cdhit.Breps = read.table("./20180222/table/Reps.clster.hist.txt", header=T, sep="\t")
cdhit.Breps.melted <- melt(cdhit.Breps, id = "Tag")
head(cdhit.Breps)
head(cdhit.Breps.melted)
subset(cdhit.Breps.melted, value != 0 )
(ggplot( data = subset(cdhit.Breps.melted, value != 0 ), aes(x=Tag, y=value, color = variable, fill=variable))
  #+ geom_point(size=0.25)
  + geom_col(position="dodge")
  + coord_cartesian(ylim=c(0,10),xlim=c(0,250))
  + theme_minimal()
  + theme( legend.position = c(0.7,0.6),
           legend.background=element_rect(fill="white", size=0.5, color="black")
           )
  + labs(x="Copy number", y="Gene count" )
)
ggsave("./20180222/99.figures/Reps-clster-hist.pdf", device = 'pdf', width = 8, height = 6, units = 'in', useDingbats=F)
ggsave("/media/cdc_documents/Abstracts/ASM2018/poster-figs/Reps-clster-hist02.pdf", device = 'pdf', width = 8, height = 6, units = 'in', useDingbats=F)

