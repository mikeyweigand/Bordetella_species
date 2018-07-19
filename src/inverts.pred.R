#!/usr/bin/R Rscript
library(optparse)
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
#opt$ref = ("/home/yrh8/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/20180614-Bp-check-gap1500-cat.NR.invertALL.txt")
#opt$ref = ("/home/yrh8/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/invert-predict/Bsp-combined-20180719.inverts.txt")
#opt$file = ("/home/yrh8/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/invert-predict/tmp-centers.txt")

### Set of observed single inversions for model building ###
refs.df <- read.table(opt$ref,header=F,sep="\t")
refs.symmetric = subset(refs.df, V11 == "ori" | V11 == "term")
ratio = log(as.numeric(as.character(refs.symmetric$V9))/as.numeric(as.character(refs.symmetric$V10)))
refs.ratio = cbind(refs.symmetric, ratio)

# duplicate ratios and invert #
refs.invert = cbind(refs.ratio[1:8],refs.ratio[c("V10","V9","V11","V12")], lapply(refs.ratio[c("ratio")], function(x) x*-1))
colnames(refs.invert)[9] = "V9"
colnames(refs.invert)[10] = "V10"
refs.double = rbind.data.frame(refs.ratio, refs.invert)
rownames(refs.double) <- 1:nrow(refs.double)


# remove outliers based on Left:Right ratio
refs.passed = refs.double[!refs.double$ratio %in% boxplot.stats(refs.double$ratio)$out,]

double.lens = refs.passed[,c('V9','V10')]
double.lens[] <- lapply(double.lens, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})

# linear regression
sym.lm2 = lm( V10 ~ V9 - 1, data=double.lens)
#summary(sym.lm2)


### Predict boundaries from new coordinates ###
IS481.df <- read.table(opt$file, header=F, sep="\t")
IS481.test = data.frame(V9 = IS481.df$V5)
IS481.pred = predict(sym.lm2, newdata = IS481.test, interval="pred",level=0.95)
IS481.df = cbind(IS481.df, format(IS481.pred,digits=5))

#outfile <- sub(".txt","-out.txt", opt$file)
write.table(IS481.df, file=opt$out,quote=F,sep="\t",row.names=F,col.names=F)



### plotting ###
#calculate prediction interval boundaries

# new.test = data.frame(V9 = double.lens$V9)
# new.pred = predict(sym.lm2, newdata = new.test, interval="pred",level=0.95)
# new.df = cbind(new.test, new.pred)
# head(new.df)
# 
# wtf.matches = ("/home/yrh8/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/invert-predict/J549-IS481-inverts-20180719.txt")
# wtf.df <- read.table(wtf.matches,header=F,sep="\t")
# sapply(wtf.df, class)
# head(wtf.df)
# break1 = c(0, 250,500,750,1000,1250,1500)
# 
# (ggplot()
#   + geom_point(data=wtf.df, mapping=aes(x=V10/1000,y=V11/1000), size =1, color="gray", alpha=0.7,stroke=0)
#   + geom_smooth(method=lm, data=double.lens,aes(x=V9/1000, y=V10/1000), formula = y ~ x - 1, fullrange=T)
#   + geom_point(data=subset(refs.symmetric, V12 == 'Bp'),
#                aes(x=(as.numeric(as.character(V9))/1000),
#                    y=(as.numeric(as.character(V10))/1000),
#                    color=V11),
#                size=2.4, alpha=0.6, shape=19,stroke=0)
#   + geom_point(data=subset(refs.symmetric, V12 == 'Bpp'),
#                aes(x=(as.numeric(as.character(V9))/1000),
#                    y=(as.numeric(as.character(V10))/1000),
#                    color=V11),
#                size=2.4, alpha=0.6, shape=15,stroke=0)
#   + geom_point(data=subset(refs.symmetric, V12 == 'Bho'),
#                aes(x=(as.numeric(as.character(V9))/1000),
#                    y=(as.numeric(as.character(V10))/1000),
#                    color=V11),
#                size=2.4, alpha=0.6, shape=17,stroke=0)
#   + geom_point(data=subset(refs.passed, V12 == 'Bp'),
#                aes(x=(as.numeric(as.character(V9))/1000),
#                    y=(as.numeric(as.character(V10))/1000),color=V11),
#                size=2.2, shape=1 ,stroke=1)
#   + geom_point(data=subset(refs.passed, V12 == 'Bpp'),
#                aes(x=(as.numeric(as.character(V9))/1000),
#                    y=(as.numeric(as.character(V10))/1000),color=V11),
#                size=2.2, shape=0 ,stroke=1)
#   + geom_point(data=subset(refs.passed, V12 == 'Bho'),
#                aes(x=(as.numeric(as.character(V9))/1000),
#                    y=(as.numeric(as.character(V10))/1000),color=V11),
#                size=2.2, shape=2 ,stroke=1)
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
#   #+ coord_fixed(ratio=1)
# 
# )
# 
# 
# ggsave("/home/yrh8/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/invert-predict/Bsp-combined-model-20180614.pdf", device = 'pdf', width = 6, height = 6, units = 'in', useDingbats=F)

# (ggplot( data = refs.ratio, aes(y=V7/1000,x=V11,color=V11)) #inverstion size
#   + geom_jitter(width = 0.2, size=0.75) #, color="red")
#   + geom_boxplot(outlier.size = 2, outlier.stroke = 0, alpha=0)
#   #+ geom_density(alpha=0.2, adjust=0.25)
#   #+ geom_histogram(binwidth=100, alpha=0.4, position='dodge')
#   + theme_minimal()
# )
# 
# (ggplot( data = refs.double, aes(x=ratio,color=V11, fill=V11))
#   + geom_density(alpha=0.1)
#   #+ geom_histogram(binwidth = 0.2, alpha =0.4)
#   + scale_x_continuous(breaks = round(seq(min(refs.double$ratio), max(refs.double$ratio), by = 0.4),1))
#   + theme_minimal()
#   + xlim(0,3)
#   + labs( x = "Log ratio")
# )

