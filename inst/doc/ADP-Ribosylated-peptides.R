### R code from vignette source 'ADP-Ribosylated-peptides.Rnw'

###################################################
### code chunk number 1: setup
###################################################



###################################################
### code chunk number 2: loaddata
###################################################
library(protViz)


###################################################
### code chunk number 3: ADP-Ribosylated-peptides.Rnw:93-127 (eval = FALSE)
###################################################
## 
## rawUrl <- paste0("http://ftp.pride.ebi.ac.uk",
##   "/pride/data/archive/2021/05/PXD017013/20171220_15_Muscle_HCD35.raw")
## 
## f <- basename(rawUrl)
## download.file(rawUrl, f )
## 
## scans <- c(9210, 13738, 14908, 7590, 10718) 
## 
## ## read spectra
## ## remove peaks with no intensity
## ADPR.ms2 <- rawrr::readSpectrum(f, scans) |>
##   lapply(function(x){
##     idx <- x$intensity > 0
##     list(mZ=x$mZ[idx], intensity=x$intensity[idx], scan=x$scan)
##   })
## 
## ## peak assignments
## ADPR.annotation <-
##   readr::read_delim("/Users/cp//Downloads/2020-05-27_InputToLabelSpectra.tsv",
##                 "\t", escape_double = FALSE, trim_ws = TRUE)
## 
## ## subsetting
## ADPR.annotation <-
##   ADPR.annotation[,c('scanNr', 'PepSeq', 'mz', 'LabelLow', 'color')] |>
##   as.data.frame()
## 
## ## some render metadata 
## ADPR.lim <- readr::read_delim("/Users/cp/Downloads/lim.txt",
##                   ",", escape_double = FALSE, trim_ws = TRUE) |>
##   as.data.frame()
## 
## save(ADPR.annotation, ADPR.ms2, ADPR.lim,
##      file="/tmp/ADPR.RData", compression_level = 9, compress = TRUE)


###################################################
### code chunk number 4: ADP-Ribosylated-peptides.Rnw:133-153
###################################################
## Heuristic to determine a useful y-axis range.
## While we deal with profile data we have to
## find the most intense peak within a mass window.
.findLocalMaxIntensity <-
  function(q, mZ, intensity, stepsize = 20, eps = 0.005){
  n <- length(mZ)
  idx <- protViz::findNN(q, mZ) |>
    vapply(function(i){
    i.max <- i
    
    for (j in seq(i - stepsize, i + stepsize)){
      if(0 < j & j <= n)
        if (intensity[j] > intensity[i.max])
          i.max <- j
    }
    i.max
  }, FUN.VALUE = 1)

  intensity[idx]
}


###################################################
### code chunk number 5: ADP-Ribosylated-peptides.Rnw:157-194
###################################################
## Adapted protViz::peakplot plot function
.peakplot <-
  function(x, mZ, intensity, lim, ...){
    p.i <- .findLocalMaxIntensity(x$mz, mZ, intensity)
    sn <- unique(x$scanNr)
    cutoff <- max(p.i) * lim$rintensity / 100
    
    plot(intensity ~ mZ,
         type = 'h',
         xlab = 'm/z',
         ylab = 'Relative Intensity [%]',
         col = 'lightgrey',
         xlim = c(lim$xmin, lim$xmax),
         ylim = c(0,  cutoff),
         axes = FALSE);
    
    legend("topright", "",  title= unique(x$PepSeq), bty='n',cex=2)
    legend("right", sprintf("% 10.3f   %s", x$mz,x$LabelLow),
           title= "Fragment Ions", bty='n',cex=0.75)
    
    axis(2, seq(0, max(intensity), length=11), round(seq(0, 100, length = 11)))
    
    points(x$mz, p.i,  col=x$color, type='h', lwd=2)
    points(x$mz, p.i,  col=x$color, pch=16,cex=0.5)
    
    select <- p.i < 0.75 * max(intensity)
    
    text(x$mz, p.i + 0.0125 * cutoff,
         x$LabelLow, adj = c(0,0), cex=1.0, srt=90, , col=x$color)
    idx <- p.i > cutoff
    
    axis(1)
    axis(3, x$mz[idx],
         paste(x$LabelLow[idx], "(", round(100 * p.i[idx] / max(p.i)), "%)", sep=''),
         cex=0.3)
    box()
  }


###################################################
### code chunk number 6: ADP-Ribosylated-peptides.Rnw:201-209
###################################################
scan <- 9210
idx <- which(vapply(ADPR.ms2, function(x)x$scan, 1) == scan) 
lim <- ADPR.lim[ADPR.lim$scan==scan,]
.peakplot(
    x = ADPR.annotation[ADPR.annotation$scanNr == scan,],
    mZ = ADPR.ms2[[idx]]$mZ,
    intensity = ADPR.ms2[[idx]]$intensity,
    lim)


###################################################
### code chunk number 7: ADP-Ribosylated-peptides.Rnw:217-225
###################################################
scan <- 13738
idx <- which(vapply(ADPR.ms2, function(x)x$scan, 1) == scan) 
lim <- ADPR.lim[ADPR.lim$scan==scan,]
.peakplot(
    x = ADPR.annotation[ADPR.annotation$scanNr == scan,],
    mZ = ADPR.ms2[[idx]]$mZ,
    intensity = ADPR.ms2[[idx]]$intensity,
    lim = lim)


###################################################
### code chunk number 8: ADP-Ribosylated-peptides.Rnw:234-242
###################################################
scan <- 14908
idx <- which(vapply(ADPR.ms2, function(x)x$scan, 1) == scan) 
lim <- ADPR.lim[ADPR.lim$scan==scan,]
.peakplot(
    x = ADPR.annotation[ADPR.annotation$scanNr == scan,],
    mZ = ADPR.ms2[[idx]]$mZ,
    intensity = ADPR.ms2[[idx]]$intensity,
    lim)


###################################################
### code chunk number 9: ADP-Ribosylated-peptides.Rnw:250-258
###################################################
scan <- 7590
idx <- which(vapply(ADPR.ms2, function(x)x$scan, 1) == scan) 
lim <- ADPR.lim[ADPR.lim$scan==scan,]
.peakplot(
    x = ADPR.annotation[ADPR.annotation$scanNr == scan,],
    mZ = ADPR.ms2[[idx]]$mZ,
    intensity = ADPR.ms2[[idx]]$intensity,
    lim)


###################################################
### code chunk number 10: ADP-Ribosylated-peptides.Rnw:266-274
###################################################
scan <- 10718
idx <- which(vapply(ADPR.ms2, function(x)x$scan, 1) == scan) 
lim <- ADPR.lim[ADPR.lim$scan==scan,]
.peakplot(
    x = ADPR.annotation[ADPR.annotation$scanNr == scan,],
    mZ = ADPR.ms2[[idx]]$mZ,
    intensity = ADPR.ms2[[idx]]$intensity,
    lim)


