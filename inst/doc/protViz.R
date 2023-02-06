### R code from vignette source 'protViz.Rnw'

###################################################
### code chunk number 1: protViz.Rnw:44-45
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: protViz.Rnw:88-91
###################################################
library(protViz)
fname <- system.file("extdata", name='P12763.fasta', package = "protViz")
F <- Fasta$new(fname)


###################################################
### code chunk number 3: protViz.Rnw:95-96
###################################################
substr(F$getSequences(), 1, 60)


###################################################
### code chunk number 4: protViz.Rnw:99-100
###################################################
(fetuin <- F$getTrypticPeptides())


###################################################
### code chunk number 5: protViz.Rnw:111-112
###################################################
mass <- protViz::parentIonMass(fetuin)


###################################################
### code chunk number 6: protViz.Rnw:116-117
###################################################
hydrophobicity <- protViz::ssrc(fetuin)


###################################################
### code chunk number 7: xtable1
###################################################
library(xtable)
print(xtable(data.frame(peptide = names(hydrophobicity), 
  mass = parentIonMass(names(hydrophobicity)), hydrophobicity=hydrophobicity),
  caption="parent ion mass and hydrophobicity values of the tryptic digested protein \texttt{P12763}.",  label="Table:hydrophobicity"), include.rownames=FALSE, scalebox="0.75")


###################################################
### code chunk number 8: protViz.Rnw:133-139
###################################################
op <- par(mfrow = c(1, 1))
plot(hydrophobicity ~ mass, 
  log = 'xy', pch = 16, col = '#88888888', cex = 2,
  main = "sp|P12763|FETUA_BOVIN Alpha-2-HS-glycoprotein",
  sub = 'tryptic peptides')
text(mass, hydrophobicity, fetuin, pos=2, cex=0.5, col = '#CCCCCC88')


###################################################
### code chunk number 9: protViz.Rnw:147-148
###################################################
defaultIon


###################################################
### code chunk number 10: protViz.Rnw:151-176
###################################################
## plot in-silico fragment ions of
peptides <- c('HTLNQIDSVK', 'ALGGEDVR', 'TPIVGQPSIPGGPVR')

pims <- peptides |> protViz::parentIonMass()
fis <- peptides |> protViz::fragmentIon()

par(mfrow = c(3, 1)); 
rv <- mapply(FUN = function(fi, pim, peptide){
    plot(0,0,
        xlab='m/Z', ylab='',
        xlim = range(c(fi$b, fi$y)),
        ylim = c(0,1),
        type = 'n', axes = FALSE,
        sub=paste(pim, "Da"));

    axis(1, fi$b,round(fi$b, 2))

    pepSeq <- strsplit(peptide, "")
    axis(3, fi$b, pepSeq[[1]])

    abline(v = fi$b, col='red', lwd=2) 
    abline(v = fi$y, col='blue',lwd=2)
    abline(v = fi$c, col='orange') 
    abline(v = fi$z, col='cyan')
  }, fis, pims, peptides)


###################################################
### code chunk number 11: protViz.Rnw:180-183
###################################################
Hydrogen<-1.007825
(fi.HTLNQIDSVK.1 <- fragmentIon('HTLNQIDSVK'))[[1]]
(fi.HTLNQIDSVK.2 <-(fi.HTLNQIDSVK.1[[1]] + Hydrogen) / 2)


###################################################
### code chunk number 12: protViz.Rnw:196-235
###################################################
    peptideSequence <- 'HTLNQIDSVK'
    spec <- list(scans=1138,
        title="178: (rt=22.3807) [20080816_23_fetuin_160.RAW]",
        rtinseconds=1342.8402,
        charge=2,
        mZ=c(195.139940, 221.211970, 239.251780, 290.221750, 
    316.300770, 333.300050, 352.258420, 448.384360, 466.348830, 
    496.207570, 509.565910, 538.458310, 547.253380, 556.173940, 
    560.358050, 569.122080, 594.435500, 689.536940, 707.624790, 
    803.509240, 804.528220, 822.528020, 891.631250, 909.544400, 
    916.631600, 973.702160, 990.594520, 999.430580, 1008.583600, 
    1017.692500, 1027.605900),
        intensity=c(931.8, 322.5, 5045, 733.9, 588.8, 9186, 604.6,
    1593, 531.8, 520.4, 976.4, 410.5, 2756, 2279, 5819, 2.679e+05,
    1267, 1542, 979.2, 9577, 3283, 9441, 1520, 1310, 1.8e+04,
    587.5, 2685, 671.7, 3734, 8266, 3309))

    fi <- protViz::fragmentIon(peptideSequence)
    n <- nchar(peptideSequence)

    by.mZ <- c(fi[[1]]$b, fi[[1]]$y)
    by.label <- c(paste("b",1:n,sep=''), paste("y",n:1,sep=''))

    # should be a R-core function as findInterval!
    idx <- protViz::findNN(by.mZ, spec$mZ) 

    mZ.error <- abs(spec$mZ[idx]-by.mZ)

    plot(mZ.error[mZ.error.idx<-order(mZ.error)],
        main="Error Plot",
        pch='o',
        cex=0.5,
        sub='The error cut-off is 0.6Da (grey line).',
        log='y')
    abline(h=0.6,col='grey')
    text(1:length(by.label), 
        mZ.error[mZ.error.idx],  
        by.label[mZ.error.idx],
        cex=0.75,pos=3) 


###################################################
### code chunk number 13: protViz.Rnw:242-267
###################################################
library(protViz)

ptm.0 <- cbind(AA="-", 
    mono=0.0, avg=0.0, desc="unmodified", unimodAccID=NA)

ptm.616 <- cbind(AA='S',
    mono=-27.010899, avg=NA, desc="Substituition", 
    unimodAccID=616)

ptm.651 <- cbind(AA='N',
    mono=27.010899, avg=NA, desc="Substituition", 
    unimodAccID=651)


m <- as.data.frame(rbind(ptm.0, ptm.616, ptm.651))

genMod(c('TAFDEAIAELDTLNEESYK','TAFDEAIAELDTLSEESYK'), m$AA)

fi <- protViz::fragmentIon(c('TAFDEAIAELDTLSEESYK', 
    'TAFDEAIAELDTLNEESYK', 'TAFDEAIAELDTLSEESYK', 
    'TAFDEAIAELDTLNEESYK'), 
        modified=c('0000000000000200000', 
        '0000000000000100000', '0000000000000000000', 
        '0000000000000000000'), 
    modification=m$mono)


###################################################
### code chunk number 14: protViz.Rnw:275-280
###################################################
data(msms)
op <- par(mfrow = c(2, 1))
protViz::peakplot("TAFDEAIAELDTLNEESYK", msms[[1]])
protViz::peakplot("TAFDEAIAELDTLSEESYK", msms[[2]])
par(op)


###################################################
### code chunk number 15: peptideSearch
###################################################
.peptideFragmentIonSpectrumMatch <- function (x, 
                           peptideSet, 
                           framentIonMassToleranceDa = 0.1) 
{
  ## Here we ignore the peptide mass
  # peptideMassTolerancePPM = 5
  # query.mass <- ((x$pepmass[1] * x$charge)) - (1.007825 * (x$charge - 1))
  # eps <- query.mass * peptideMassTolerancePPM * 1e-06
  # pimIdx <- protViz::parentIonMass(peptideSequence)
  # lower <- protViz::findNN(query.mass - eps, pimIdx)
  # upper <- protViz::findNN(query.mass + eps, pimIdx)
  
  
  rv <- lapply(peptideSet, FUN = protViz::psm, spec = x, plot = FALSE) |>
    lapply(FUN = function(p) {
      ## determine peaks considered as hits
      idx <- abs(p$mZ.Da.error) < framentIonMassToleranceDa
      intensityRatio <- sum(x$intensity[idx]) / sum(x$intensity)
      
      ## derive objectives for a good match
      data.frame(nHits=sum(idx), intensityRatio = intensityRatio)
    }) |>
    Reduce(f=rbind)
  
  
  idx.tophit <- which(rv$nHits == max(rv$nHits))[1]
  
  data.frame(peptideMatch = peptideSet[idx.tophit],
             nHits = max(rv$nHits),
             nPeaks = length(x$mZ),
             intensityRatio = rv$intensityRatio[idx.tophit]
  )
}


###################################################
### code chunk number 16: protViz.Rnw:323-324
###################################################
peptideSet <- c("ELIVSK", 'TAFDEAIAELDTLNEESYK','TAFDEAIAELDTLSEESYK')


###################################################
### code chunk number 17: protViz.Rnw:328-334
###################################################
mZ <- protViz::fragmentIon("TAFDEAIAELDTLNEESYK")[[1]] |>
  unlist() |> sort()

intensity <- mZ |> length() |> sample()

msms.insilico <- list(mZ = mZ, intensity = intensity)


###################################################
### code chunk number 18: protViz.Rnw:338-342
###################################################
peptideSet.rev <- peptideSet |>
  sapply(FUN = function(x){
    strsplit(x, "")[[1]] |> rev() |> paste0(collapse = "")
  })


###################################################
### code chunk number 19: protViz.Rnw:347-352
###################################################
lapply(list(msms[[1]], msms[[2]], msms.insilico),
       FUN = .peptideFragmentIonSpectrumMatch,
       peptideSet = c(peptideSet, peptideSet.rev),
       framentIonMassToleranceDa = 0.05) |>
  Reduce(f=rbind)


###################################################
### code chunk number 20: protViz.Rnw:372-404
###################################################
library(lattice)
data(fetuinLFQ)

cv<-1-1:7/10
t<-trellis.par.get("strip.background")
t$col<-(rgb(cv,cv,cv))
trellis.par.set("strip.background",t)

print(xyplot(abundance~conc|prot*method,
    groups=prot,
    xlab="Fetuin concentration spiked into experiment [fmol]",
    ylab="Abundance",
    aspect=1,
    data=fetuinLFQ$t3pq[fetuinLFQ$t3pq$prot 
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
### code chunk number 21: protViz.Rnw:416-441
###################################################
data(pgLFQfeature)
data(pgLFQprot)

featureDensityPlot<-function(data, n=ncol(data), nbins=30){
    my.col<-rainbow(n);
    mids<-numeric()
    density<-numeric()
    for (i in 1:n) { 
        h<-hist(data[,i],nbins, plot=FALSE)
        mids<-c(mids, h$mids)
        density<-c(density, h$density)
    }
    plot(mids,density, type='n')
    for (i in 1:n) { 
        h<-hist(data[,i],nbins, plot=FALSE)
        lines(h$mids,h$density, col=my.col[i])
    }
    legend("topleft", names(data), cex=0.5,
        text.col=my.col
    )
}

par(mfrow=c(1,1)); 
featureDensityPlot(asinh(pgLFQfeature$"Normalized abundance"),
    nbins=25)


###################################################
### code chunk number 22: protViz.Rnw:447-461
###################################################
op<-par(mfrow=c(1,1),mar=c(18,18,4,1),cex=0.5)
samples<-names(pgLFQfeature$"Normalized abundance")
image(cor(asinh(pgLFQfeature$"Normalized abundance")),
    col=gray(seq(0,1,length=20)),
    main='pgLFQfeature correlation',
    axes=FALSE)

axis(1,at=seq(from=0, to=1, 
    length.out=length(samples)), 
    labels=samples, las=2)

axis(2,at=seq(from=0, to=1, 
    length.out=length(samples)), labels=samples, las=2)
par(op)


###################################################
### code chunk number 23: protViz.Rnw:467-477
###################################################
op<-par(mfrow=c(1,1),mar=c(18,18,4,1),cex=0.5)
image(cor(asinh(pgLFQprot$"Normalized abundance")),
    main='pgLFQprot correlation',
    axes=FALSE,
    col=gray(seq(0,1,length=20)))
axis(1,at=seq(from=0, to=1, 
    length.out=length(samples)), labels=samples, las=2)
axis(2,at=seq(from=0, to=1, 
    length.out=length(samples)), labels=samples, las=2)
par(op)


###################################################
### code chunk number 24: protViz.Rnw:485-491
###################################################
par(mfrow=c(2,2),mar=c(6,3,4,1))
ANOVA<-pgLFQaov(pgLFQprot$"Normalized abundance", 
    groups=as.factor(pgLFQprot$grouping), 
    names=pgLFQprot$output$Accession,
    idx=c(15,16,196,107),
    plot=TRUE)


###################################################
### code chunk number 25: protViz.Rnw:505-514
###################################################
data(iTRAQ)
x<-rnorm(100)
par(mfrow=c(3,3),mar=c(6,4,3,0.5));
for (i in 3:10){
    qqnorm(asinh(iTRAQ[,i]), 
        main=names(iTRAQ)[i])
    qqline(asinh(iTRAQ[,i]), col='grey')
}
b<-boxplot(asinh(iTRAQ[,c(3:10)]), main='boxplot iTRAQ')


###################################################
### code chunk number 26: protViz.Rnw:521-554
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
                  
                  
par(mfrow=c(2,3),mar=c(6,3,4,1))
for (i in 1:nrow(group1Protein)){
    boxplot.color="#ffcccc"
    tt.p_value <- t.test(as.numeric(group1Protein[i,]), 
        as.numeric(group2Protein[i,]))$p.value       

    if (tt.p_value < 0.05)
        boxplot.color='lightgreen'

    b <- boxplot(as.numeric(group1Protein[i,]), 
        as.numeric(group2Protein[i,]),
        main=row.names(group1Protein)[i],
        sub=paste("t.test: p-value =", round(tt.p_value,2)),
        col=boxplot.color,
        axes=FALSE)
    axis(1, 1:2, c('group_1','group_2')); axis(2); box()

    points(rep(1,b$n[1]), as.numeric(group1Protein[i,]), col='blue')
    points(rep(2,b$n[2]), as.numeric(group2Protein[i,]), col='blue')
}


###################################################
### code chunk number 27: protViz.Rnw:564-571
###################################################
data(iTRAQ)
q <- iTRAQ2GroupAnalysis(data=iTRAQ, 
    group1=c(3,4,5,6), 
    group2=7:10, 
    INDEX=paste(iTRAQ$prot,iTRAQ$peptide), 
    plot=FALSE)
q[1:10,]


###################################################
### code chunk number 28: protViz.Rnw:582-584
###################################################
data(pressureProfile)
ppp(pressureProfile)


###################################################
### code chunk number 29: protViz.Rnw:594-610
###################################################
pp.data<-pps(pressureProfile, time=seq(25,40,by=5))
print(xyplot(Pc ~ as.factor(file) | paste("time =", 
    as.character(time), "minutes"),
    panel = function(x, y){
        m<-sum(y)/length(y)
        m5<-(max(y)-min(y))*0.05
        panel.abline(h=c(m-m5,m,m+m5), 
            col=rep("#ffcccc",3),lwd=c(1,2,1))
        panel.grid(h=-1, v=0)
        panel.xyplot(x, y)
    },
    ylab='Pc [psi]',
    layout=c(1,4),
    sub='The three red lines indicate the average plus min 5%.',
    scales = list(x = list(rot = 45)),
    data=pp.data))


###################################################
### code chunk number 30: protViz.Rnw:616-621
###################################################
pp.data<-pps(pressureProfile, time=seq(0,140,length=128))
print(levelplot(Pc ~ time * as.factor(file),
    main='Pc(psi)',
    data=pp.data,
    col.regions=rainbow(100)[1:80]))


###################################################
### code chunk number 31: sessioninfo
###################################################
toLatex(sessionInfo())


