# This script plots the observed and predicted inversions using circos. Intended to be run interactively in Rstudio.

#setwd("~/Documents/Bordetella_species/results/mauve/20180221-Bp/04.colinear-mcl/invert-predict/")
setwd("~/Documents/Bordetella_species/results/mauve/Bsp-combined-20190903/")
detach("package:sna", unload=T)
library(circlize)

## import all observed single invert coordinates
BP.obs <- read.table("./20180614-Bp-check-gap1500-cat.NR.invertALL.txt",
                     sep="\t",header=F)
BP.symm = subset( BP.obs, V8 == 'Symmetric')
BP.asym = subset( BP.obs, V8 == 'Asymmetric')

BP.s1 = data.frame("Bp_J549", BP.symm$V5, BP.symm$V5, BP.symm$V7)
BP.s2 = data.frame("Bp_J549", BP.symm$V6, BP.symm$V6, BP.symm$V7)
colnames(BP.s1) = c("chr","start","end","value1")
colnames(BP.s2) = c("chr","start","end","value1")
BP.a1 = data.frame("Bp_J549", BP.asym$V5, BP.asym$V5, BP.asym$V7)
BP.a2 = data.frame("Bp_J549", BP.asym$V6, BP.asym$V6, BP.asym$V7)
colnames(BP.a1) = c("chr","start","end","value1")
colnames(BP.a2) = c("chr","start","end","value1")
head(BP.s1)

## import positions of ISE in example B. pertussis J549
#J549.pred <- read.table("./J549-ISall-inverts-20180719.txt",
J549.pred <- read.table("./predict/Bp_J549_IScollapsed-ISall.20190911.txt",
                        sep="\t",header=F,
                        colClasses = c("character", "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
J549.1 = data.frame("Bp_J549",as.numeric(J549.pred$V3), J549.pred$V4, J549.pred$V5)
colnames(J549.1) = c("chr","start","end","value1")
J549.2 = data.frame("Bp_J549",as.numeric(J549.pred$V7), J549.pred$V8, J549.pred$V9)
colnames(J549.2) = c("chr","start","end","value1")

J549.bls <- read.table("../20180221-Bp/04.colinear-mcl/invert-predict/J549-ISall-20180719.bed",
                        sep="\t",header=F,
                        colClasses = c("character", "numeric","numeric"))
J549.ISE = data.frame("Bp_J549", J549.bls$V2,J549.bls$V3, 0)
colnames(J549.ISE) = c("chr","start","end","value1")
J549.ISE$value1 = with(J549.ISE, ifelse(start < end, 1, -1))
J549.ISEbed = J549.ISE
J549.ISEbed$start = with(J549.ISE, ifelse(start < end, start, end))
J549.ISEbed$end = with(J549.ISE, ifelse(start > end, start, end))
head(J549.ISE)
nrow(J549.ISEbed)

## import all observed breakpoint coordinates between J549 and strains of varied ptxP backgrounds
J549.ptxp1 <- read.table("../20180221-Bp/04.colinear-mcl/J549-breakpoints/J549-breakpoints.ptxp12.sorted.uniq.merged1000.txt",
                         sep="\t",header=F,
                         colClasses = c("numeric"))
J549.bks1 = data.frame("Bp_J549",J549.ptxp1$V1, J549.ptxp1$V1+1,1)
colnames(J549.bks1) = c("chr","start","end","value1")
head(J549.bks1)

J549.ptxp3 <- read.table("../20180221-Bp/04.colinear-mcl/J549-breakpoints/J549-breakpoints.ptxp3.sorted.uniq.merged1000.txt",
                         sep="\t",header=F,
                         colClasses = c("numeric"))
J549.bks3 = data.frame("Bp_J549",J549.ptxp3$V1, J549.ptxp3$V1+1,1)
colnames(J549.bks3) = c("chr","start","end","value1")

J549.bkpts1k <- read.table("../20180221-Bp/04.colinear-mcl/J549-breakpoints/J549-breakpoints.sorted.uniq.merged1000.txt",
                            sep="\t",header=F,
                            colClasses = c("numeric"))
J549.breaks1k = data.frame("Bp_J549",J549.bkpts1k$V1, J549.bkpts1k$V1+1,1)
colnames(J549.breaks1k) = c("chr","start","end","value1")


bp = data.frame(
  name  = c("Bp_J549"),
  start = c(1),
  end   = c(4106572))

oric = data.frame(
  chr = "Bp_J549",
  start = c(2049460,4104849),
  end = c(2049460,4106572),
  value = c(0,0)
)

## Draw circular plot of J549 with observed breakpoints, IS elements, and single inversions.
#pdf(file = "../J549-breakpoints/J549-circlize-20180823.pdf",
pdf(file = "./circlize/J549-circlize-20190911.pdf",
    width = 6, height = 6, useDingbats = F)
circos.par("track.height"=0.08, start.degree=90, gap.degree=0)
circos.genomicInitialize(bp, major.by = 250000, sector.names = NA)
circos.genomicTrack(J549.ISEbed, # ylim=c(-1,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, baseline = 0, type = "h", col="blue")
                    })
circos.genomicTrack(J549.bks1,  ylim=c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h", col="orange")
                    })
circos.genomicTrack(J549.bks3,  ylim=c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h", col="green")
                    })
circos.genomicLink(J549.1,J549.2, col = "grey", border=NA) #predicted
circos.genomicLink(BP.s1,BP.s2, col = "black", border=NA) #observed, symmetric
circos.genomicLink(BP.a1,BP.a2, col = "red", border=NA) #observed, asymmetric
circos.clear()
dev.off()


##### Draw plot for B. holmesii C690

bho = data.frame(
  name  = c("Bho_C690"),
  start = c(1),
  end   = c(3698268))

C690.bls <- read.table("../20180221-Bp/04.colinear-mcl/invert-predict/C690-ISall-20180719.bed",
                       sep="\t",header=F,
                       colClasses = c("character", "numeric","numeric"))
C690.ISE = data.frame("Bho_C690", C690.bls$V2,C690.bls$V3, 0)
colnames(C690.ISE) = c("chr","start","end","value1")
C690.ISE$value1 = with(C690.ISE, ifelse(start < end, 1, -1))
C690.ISEbed = C690.ISE
C690.ISEbed$start = with(C690.ISE, ifelse(start < end, start, end))
C690.ISEbed$end = with(C690.ISE, ifelse(start > end, start, end))
head(C690.ISEbed)
nrow(C690.ISEbed)

#C690.pred <- read.table("./C690-ISall-inverts-20180719.txt",
C690.pred <- read.table("./predict/Bho_C690-ISall.20190911.txt",
                        sep="\t",header=F,
                        colClasses = c("character", "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
C690.1 = data.frame("Bho_C690",as.numeric(C690.pred$V3), C690.pred$V4, C690.pred$V5)
colnames(C690.1) = c("chr","start","end","value1")
C690.2 = data.frame("Bho_C690",as.numeric(C690.pred$V7), C690.pred$V8, C690.pred$V9)
colnames(C690.2) = c("chr","start","end","value1")
head(C690.1)

#C690.inverts <- read.table("/home/yrh8/Documents/Bordetella_species/results/mauve/20180501-Bsp/Bho-pairwise/C690-inverts.txt",
C690.inverts <- read.table("./Bho-collinear-20190821.02-clean.mcl.NR.1invert.sym.txt",
                           sep="\t",header=F,
                           colClasses = c("character","character","numeric","numeric","numeric","numeric","numeric","character","numeric","numeric","character"))

Bho.s1 = data.frame("Bho_C690", C690.inverts$V5, C690.inverts$V5, 1)
Bho.s2 = data.frame("Bho_C690", C690.inverts$V6, C690.inverts$V6, 1)
colnames(Bho.s1) = c("chr","start","end","value1")
colnames(Bho.s2) = c("chr","start","end","value1")
head(Bho.s2)


#pdf(file = "../J549-breakpoints/C690-circlize-20180823.pdf",
pdf(file = "./circlize/C690-circlize-20190911.pdf",
    width = 3, height = 3, useDingbats = F)
circos.par("track.height"=0.08, start.degree=90, gap.degree=0)
circos.genomicInitialize(bho, major.by = 500000, sector.names = NA)
circos.genomicTrack(C690.ISEbed, # ylim=c(-1,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, baseline = 0, type = "h", col="blue")
                    })
circos.genomicLink(C690.1,C690.2, col = "grey", border=NA) #predicted
circos.genomicLink(Bho.s1,Bho.s2, col = "black", border=NA) #observed, symmetric
circos.clear()
dev.off()


##### Draw plot for B. parapertussis B271

bpp = data.frame(
  name  = c("Bpp_B271"),
  start = c(1),
  end   = c(4773897))

B271.bls <- read.table("../20180221-Bp/04.colinear-mcl/invert-predict/B271-ISall-20180719.bed",
                       sep="\t",header=F,
                       colClasses = c("character", "numeric","numeric"))
B271.ISE = data.frame("Bpp_B271", B271.bls$V2,B271.bls$V3, 0)
colnames(B271.ISE) = c("chr","start","end","value1")
B271.ISE$value1 = with(B271.ISE, ifelse(start < end, 1, -1))
B271.ISEbed = B271.ISE
B271.ISEbed$start = with(B271.ISE, ifelse(start < end, start, end))
B271.ISEbed$end = with(B271.ISE, ifelse(start > end, start, end))
head(B271.ISEbed)
nrow(B271.ISEbed)

#B271.pred <- read.table("./B271-ISall-inverts-20180719.txt",
B271.pred <- read.table("./predict/Bpp_B271-ISall.20190911.txt",
                        sep="\t",header=F,
                        colClasses = c("character", "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
B271.1 = data.frame("Bpp_B271",as.numeric(B271.pred$V3), B271.pred$V4, B271.pred$V5)
colnames(B271.1) = c("chr","start","end","value1")
B271.2 = data.frame("Bpp_B271",as.numeric(B271.pred$V7), B271.pred$V8, B271.pred$V9)
colnames(B271.2) = c("chr","start","end","value1")
head(B271.1)

#B271.inverts <- read.table("/home/yrh8/Documents/Bordetella_species/results/mauve/20180501-Bsp/Bpp-pairwise/A005-E738.invert.txt",
B271.inverts <- read.table("./Bpp-collinear-20190822.02-clean.mcl.NR.1invert.sym.txt",
                          sep="\t",header=F,
                           colClasses = c("character","character","numeric","numeric","numeric","numeric","numeric","character","numeric","numeric","character"))

Bpp.s1 = data.frame("Bpp_B271", B271.inverts$V5, B271.inverts$V5, 1)
Bpp.s2 = data.frame("Bpp_B271", B271.inverts$V6, B271.inverts$V6, 1)
colnames(Bpp.s1) = c("chr","start","end","value1")
colnames(Bpp.s2) = c("chr","start","end","value1")
head(Bpp.s2)

#pdf(file = "../J549-breakpoints/B271-circlize-20180823.pdf",
pdf(file = "./circlize/B271-circlize-20190911.pdf",
    width = 3, height = 3, useDingbats = F)
circos.par("track.height"=0.08, start.degree=90, gap.degree=0)
circos.genomicInitialize(bpp, major.by = 500000, sector.names = NA)
circos.genomicTrack(B271.ISEbed, # ylim=c(-1,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, baseline = 0, type = "h", col="blue")
                    })
circos.genomicLink(B271.1,B271.2, col = "grey", border=NA) #predicted
circos.genomicLink(Bpp.s1,Bpp.s2, col = "black", border=NA) #observed, symmetric
circos.clear()
dev.off()
