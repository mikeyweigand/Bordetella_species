#!/usr/bin/R Rscript
library("optparse")
library(ggplot2)
library(scales)

option_list = list(
  make_option(c("-r", "--ref"), type="character", default=NULL, 
              help="set of observed single inversions", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="relative positions of IS481", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="pred-out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$ref) | is.null(opt$file)){
  print_help(opt_parser)
  stop("Required input file(s) missing.", call.=FALSE)
}

## For testing:
#opt$ref = ("/home/yrh8/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/20180425-Bp-check-gap1500-cat.NR.invertALL.txt")
#opt$file = ("/home/yrh8/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/invert-predict/tmp-centers.txt")

### Set of observed single inversions for model building ###
refs.df <- read.table(opt$ref,header=F,sep="\t")
refs.symmetric = subset(refs.df, V11 == "ori" | V11 == "term")
ratio = (as.numeric(as.character(refs.symmetric$V9))/as.numeric(as.character(refs.symmetric$V10)))
refs.ratio = cbind(refs.symmetric, ratio)

# remove outliers based on Left:Right ratio
refs.passed = refs.ratio[!refs.ratio$ratio %in% boxplot.stats(refs.ratio$ratio)$out,]

# duplicate and flip all pairs to balance ratios
just.lens = refs.passed[,c('V9','V10')]
just.lens[] <- lapply(just.lens, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
flip.lens = just.lens[,c('V10','V9')]
colnames(flip.lens) = c('V9','V10')
comb.lens = rbind(just.lens,flip.lens)

# linear regression
sym.lm = lm( V10 ~ V9, data=comb.lens)
#summary(sym.lm)

### Predict boundaries from new coordinates ###
IS481.df <- read.table(opt$file, header=F, sep="\t")
IS481.test = data.frame(V9 = IS481.df$V5)

IS481.pred = predict(sym.lm, newdata = IS481.test, interval="pred",level=0.95)
IS481.df = cbind(IS481.df, format(IS481.pred,digits=5))
#head(IS481.df)

#outfile <- sub(".txt","-out.txt", opt$file)
write.table(IS481.df, file=opt$out,quote=F,sep="\t",row.names=F,col.names=F)

### plotting ###
#calculate prediction interval boundaries

# new.test = data.frame(V9 = just.lens$V9)
# new.pred = predict(sym.lm, newdata = new.test, interval="pred",level=0.95)
# new.df = cbind(new.test, new.pred)
# break1 = c(0, 250,500,750,1000,1250,1500)
# 
# (ggplot()
#   + geom_point(data=refs.symmetric, aes(x=as.numeric(as.character(V9))/1000,
#                                       y=as.numeric(as.character(V10))/1000,
#                                       color=V11), size=2.2, alpha=0.6)
#   + geom_smooth(method=lm, data=comb.lens,aes(x=V9/1000, y=V10/1000))
#   + geom_point(data=comb.lens, mapping = aes(x=V9/1000,y=V10/1000), size=0.5)
#   #+ theme_classic()
#   + theme_minimal()
#   + theme( legend.position = c(0.79,0.3),
#            legend.background=element_rect(fill="white", size=0.75, color="black"))
#   + scale_x_continuous(labels = comma, breaks = break1)
#   + scale_y_continuous(labels = comma, breaks = break1)
#   + coord_cartesian(xlim=c(0,1500),ylim=c(0,1500))
#   + labs(x="Right replichore size (kbp)", y="Left replichore size (kbp)" )
#   + geom_abline(aes(intercept=0,slope=1),lty=2, size=0.4)
#   + geom_line(data = new.df, mapping = aes(x=V9/1000, y=lwr/1000),color="red",linetype='dashed')
#   + geom_line(data = new.df, mapping = aes(x=V9/1000, y=upr/1000),color="red",linetype='dashed')
# )


