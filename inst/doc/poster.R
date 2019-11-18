### R code from vignette source 'poster.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style-Sweave
###################################################
options(prompt = "R> ", continue = "+  ", width = 60, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: digest
###################################################
library(protViz)
irt.peptide <- as.character(
protViz::iRTpeptides$peptide)
irt.pim <- parentIonMass(irt.peptide)

op <- par(mfrow = c(1,2), 
pch = 16, 
col = rgb(0.5, 0.5, 0.5, alpha = 0.5))
hist(irt.pim, xlab="peptide mass [in Da]")

irt.ssrc <- sapply(irt.peptide, ssrc)
plot(irt.pim ~ irt.ssrc,
cex = 2,
main = 'In-silico LC-MS map')
par(op)


###################################################
### code chunk number 3: poster.Rnw:149-152
###################################################
defaultIon
peptides <- c('HTLNQIDSVK')
fi <- fragmentIon(peptides)


###################################################
### code chunk number 4: insilico
###################################################
par <- par(mfrow = c(1, 1))
pim <- parentIonMass(peptides)
for (i in 1:length(peptides)){
    plot(0,0,
        xlab='m/Z',
        ylab='',
        xlim=range(c(fi[i][[1]]$b,fi[i][[1]]$y)),
        ylim=c(0, 1),
        type='n',
        axes=FALSE,
        sub=paste( pim[i], "Da"));
    box()
    axis(1,fi[i][[1]]$b,round(fi[i][[1]]$b,2))
    pepSeq <- strsplit(peptides[i],"")
    axis(3,fi[i][[1]]$b,pepSeq[[1]])

    abline(v=fi[i][[1]]$b, col='red',lwd=2) 
    abline(v=fi[i][[1]]$c, col='orange') 
    abline(v=fi[i][[1]]$y, col='blue',lwd=2)
    abline(v=fi[i][[1]]$z, col='cyan')
}
par(op)


###################################################
### code chunk number 5: poster.Rnw:187-192
###################################################
fi.HTLNQIDSVK.1 <- 
fragmentIon('HTLNQIDSVK')[[1]]
Hydrogen <- 1.007825
fi.HTLNQIDSVK.2 <-
(fi.HTLNQIDSVK.1 + Hydrogen) / 2


###################################################
### code chunk number 6: xtable2
###################################################
df <- as.data.frame(cbind(fi.HTLNQIDSVK.1, fi.HTLNQIDSVK.2))
names(df) <- c(paste(names(fi.HTLNQIDSVK.1),1,sep=''),
  paste(names(fi.HTLNQIDSVK.2),2,sep=''))
library(xtable)
print.xtable(xtable(df,
	caption = "Singly and doubly charged fragment ions of the HTLNQIDSVK tryptic peptide of the SwissProt P12763 FETUA BOVIN Alpha-2-HS-glycoprotein protein are listed.",
	label = "Table:xtable2"), 
  include.rownames = FALSE, 
  table.placement = "H",
  scalebox = 0.8)


###################################################
### code chunk number 7: poster.Rnw:221-236
###################################################
spec <- list(scans=1138,
   title = "178: (rt=22.3807) [20080816_23_fetuin_160.RAW]",
   rtinseconds = 1342.8402,
   charge = 2,
   mZ = c(195.139940, 221.211970, 239.251780, 290.221750, 
316.300770, 333.300050, 352.258420, 448.384360, 466.348830, 
496.207570, 509.565910, 538.458310, 547.253380, 556.173940, 
560.358050, 569.122080, 594.435500, 689.536940, 707.624790, 
803.509240, 804.528220, 822.528020, 891.631250, 909.544400, 
916.631600, 973.702160, 990.594520, 999.430580, 1008.583600, 
1017.692500, 1027.605900),
   intensity=c(931.8, 322.5, 5045, 733.9, 588.8, 9186, 604.6,
1593, 531.8, 520.4, 976.4, 410.5, 2756, 2279, 5819, 2.679e+04,
1267, 1542, 979.2, 9577, 3283, 9441, 1520, 1310, 1.8e+04,
587.5, 2685, 671.7, 3734, 8266, 3309))


###################################################
### code chunk number 8: poster.Rnw:239-248
###################################################
peptideSequence <- 'HTLNQIDSVK'
str(spec, nchar.max = 25, vec.len = 2)

fi <- fragmentIon(peptideSequence)
n <- nchar(peptideSequence)
by.mZ <- c(fi[[1]]$b, fi[[1]]$y)
idx <- findNN(by.mZ, spec$mZ)
mZ.error <- abs(spec$mZ[idx]-by.mZ)
which(mZ.error < 0.3)


###################################################
### code chunk number 9: peakplot
###################################################
p <- peakplot('HTLNQIDSVK', spec)


###################################################
### code chunk number 10: LFQtrellis
###################################################
library(lattice)
data(fetuinLFQ)

cv <- 1-1:7/10
t<-trellis.par.get("strip.background")
t$col<-(rgb(cv,cv,cv))
trellis.par.set("strip.background",t)

print(xyplot(abundance ~ conc | prot * method,
    groups = prot,
    xlab = "Fetuin concentration spiked into experiment [fmol]",
    ylab = "Abundance",
    aspect = 1,
    data = fetuinLFQ$t3pq[fetuinLFQ$t3pq$prot
        %in% c('Fetuin', 'P15891', 'P32324', 'P34730'),],
    panel = function(x, y, subscripts, groups) {
        if (groups[subscripts][1] == "Fetuin")  {
            panel.fill(col="#ffcccc")
        }
        panel.grid(h=-1,v=-1)
        panel.xyplot(x, y)
        panel.loess(x,y, span=1)
        if (groups[subscripts][1] == "Fetuin")  {
            panel.text(min(fetuinLFQ$t3pq$conc),
                max(fetuinLFQ$t3pq$abundance),
                paste("R-squared:", 
                round(summary(lm(x~y))$r.squared,2)),
                cex=0.75,
                pos=4)
        }
    }
))


###################################################
### code chunk number 11: poster.Rnw:328-330
###################################################
data(pgLFQfeature)
data(pgLFQprot)


###################################################
### code chunk number 12: featureDensityPlot
###################################################
par(mfrow = c(1,1)); 
data(pgLFQfeature)
data(pgLFQprot)
protViz:::.featureDensityPlot(
asinh(
pgLFQfeature$"Normalized abundance"), 
nbins=25)


###################################################
### code chunk number 13: image1
###################################################
op <- par(mfrow=c(1,1),
mar = c(18,18,4,1),
cex=0.5)

samples <-
names(pgLFQfeature$"Normalized abundance")

image(cor(
  asinh(
    pgLFQfeature$"Normalized abundance")),
col = gray(seq(0,1,length=20)),
asp = 1,
main = 'pgLFQfeature correlation',
axes=FALSE)

axis(1, 
at=seq(from = 0, to = 1, 
length.out=length(samples)), 
labels=samples, las=2)

axis(2,
at=seq(from = 0, to = 1, 
length.out=length(samples)), 
labels=samples, las=2)
par(op)


###################################################
### code chunk number 14: image2
###################################################
op<-par(mfrow=c(1,1),mar=c(18,18,4,1), cex=0.5)
image(cor(asinh(pgLFQprot$"Normalized abundance")),
    main = 'pgLFQprot correlation',
    asp = 1,
    axes = FALSE,
    col = gray(seq(0,1,length=20)))
axis(1,at=seq(from=0, to=1, 
    length.out=length(samples)), labels=samples, las=2)
axis(2,at=seq(from=0, to=1, 
    length.out=length(samples)), labels=samples, las=2)
par(op)


###################################################
### code chunk number 15: ANOVA
###################################################
par(mfrow=c(1,4),mar=c(6,3,4,1))
ANOVA <- pgLFQaov(
  pgLFQprot$"Normalized abundance", 
    groups=as.factor(pgLFQprot$grouping), 
    names=pgLFQprot$output$Accession,
    idx=c(15,16,196,107),
    plot=TRUE)


###################################################
### code chunk number 16: iTRAQqqnorm
###################################################
data(iTRAQ)
par(mfrow = c(2,4),
    mar = c(6,4,3,0.5));
for (i in 3:10){
    qqnorm(asinh(iTRAQ[,i]), 
    asp = 1,
        main=names(iTRAQ)[i])
    qqline(asinh(iTRAQ[,i]), col='grey')
}


###################################################
### code chunk number 17: iTRAQboxplot
###################################################
b <- boxplot(asinh(iTRAQ[, c(3:10)]),
main='boxplot iTRAQ')


###################################################
### code chunk number 18: iTRAQboxplot2
###################################################
data(iTRAQ)
group1Protein<-numeric()
group2Protein<-numeric()

for (i in c(3,4,5,6))
    group1Protein<-cbind(group1Protein,
        asinh(tapply(iTRAQ[,i], paste(iTRAQ$prot), sum, na.rm=TRUE)))
         
for (i in 7:10)
    group2Protein<-cbind(group2Protein,
        asinh(tapply(iTRAQ[,i], paste(iTRAQ$prot), sum, na.rm=TRUE)))
                  
                  
par(mfrow = c(1,4), mar = c(6,3,4,1))
for (i in 1:nrow(group1Protein)){
    boxplot.color="#ffcccc"
    tt.p_value <- t.test(as.numeric(group1Protein[i,]), 
        as.numeric(group2Protein[i,]))$p.value       

    if (tt.p_value < 0.05)
        boxplot.color='lightgreen'

    b <- boxplot(as.numeric(group1Protein[i,]), 
        as.numeric(group2Protein[i,]),
        main=row.names(group1Protein)[i],
        sub=paste("t-Test: p-value =", round(tt.p_value,2)),
        col=boxplot.color,
        axes=F)
    axis(1, 1:2, c('group_1','group_2')); axis(2); box()

    points(rep(1,b$n[1]), as.numeric(group1Protein[i,]), col='blue')
    points(rep(2,b$n[2]), as.numeric(group2Protein[i,]), col='blue')
}


