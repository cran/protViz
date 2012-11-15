#R

peakplot.putlabel<-function(MASS, INTENSITY, LABEL, 
    l.col="green", 
    delta=0, 
    yMin, 
    maxIntensity=max(INTENSITY)) {

	noise<-seq(0, 2*delta, length=length(MASS))

	if(delta==0){
		segments(min(MASS), 0, max(MASS), 0, col="green", cex=2, lwd=3)
	}


	#segments(min(s$MASS), yMin, max(s$MASS), yMin, col=l.col, cex=2, lwd=3)
	#segments(min(s$MASS), yMin+delta, max(s$MASS), yMin+delta, col=l.col, cex=2, lwd=1)

	#charge<-as.character(s$CHARGE)

	segments(MASS, 
		INTENSITY+0.03*maxIntensity,
		MASS, 
		1.1*maxIntensity,
        lty=2,
		#yMin+noise,
		col=l.col, 
		pch=5, 
		cex=0.5,
		lwd=0.5)

	segments(MASS, INTENSITY, MASS, rep(0, length(MASS)), lwd=1.5, col=l.col)

#	text(MASS, yMin+noise,
#		paste(LABEL), 
#		pos=3, 
#		cex=0.66,
#		col=l.col)


	text(MASS, yMin+noise, round(MASS,2), 
		cex=0.40,
		pos=1,
		#offset=0.0,
		col="red",
		srt=90)
}

peakplot.label<-function(spec, match, yMax){

    bin.n<-10
    lab.n<-3
    bin<-seq(min(spec$mZ), max(spec$mZ), length=(bin.n+1))
    bin.min<-bin[seq(1,length(bin)-1)];
    bin.max<-bin[seq(2,length(bin))];
    bin.range<-1:lab.n


    # filtering the ions
    LABEL.abc<-(abs(match$mZ.Da.error) < 0.6) & (regexpr("[abc].*", match$label) > 0)
    LABEL.xyz<-(abs(match$mZ.Da.error) < 0.6) & (regexpr("[xyz].*", match$label) > 0)

    points(spec$mZ[match$idx[LABEL.abc]], spec$intensity[match$idx[LABEL.abc]], col="black", cex=1.0, pch=22)
    points(spec$mZ[match$idx[LABEL.xyz]], spec$intensity[match$idx[LABEL.xyz]], col="blue", cex=1.0, pch=22)


    for (i in c(1:length(bin)-1)) {
        bin.label <- bin.min[i] <= spec$mZ[match$idx] & spec$mZ[match$idx] < bin.max[i]

        bin.label.abc <- bin.label & LABEL.abc 
        bin.label.xyz <- bin.label & LABEL.xyz 


        if (sum(spec$intensity[match$idx[bin.label.abc | bin.label.xyz]]) > 0){
            d<-(yMax - max(spec$intensity[match$idx[bin.label.abc | bin.label.xyz]], na.rm=TRUE)) / 2

            if (sum(spec$intensity[match$idx[bin.label.abc]]) > 0)
            peakplot.putlabel(MASS=spec$mZ[match$idx[bin.label.abc]], 
                INTENSITY=spec$intensity[match$idx[bin.label.abc]], 
                LABEL=match$label[bin.label.abc],
                l.col="black",
                yMin=1.1 * d + max(spec$intensity[match$idx[bin.label.abc | bin.label.xyz]]), 
                delta=d, 
                maxIntensity=max(spec$intensity))

            if (sum(spec$intensity[match$idx[bin.label.xyz]]) > 0)
            peakplot.putlabel(MASS=spec$mZ[match$idx[bin.label.xyz]], 
                INTENSITY=spec$intensity[match$idx[bin.label.xyz]], 
                LABEL=match$label[bin.label.xyz],
                l.col="blue",
                yMin=d + max(spec$intensity[match$idx[bin.label.xyz | bin.label.xyz]]), 
                delta=d, 
                maxIntensity=max(spec$intensity))
        }

    }

}

peakplot<-function(peptideSequence, spec, FUN=defaultIons){ 

    m<-psm(peptideSequence, spec, FUN, plot=FALSE)

    max.intensity<-max(spec$intensity, na.rm=TRUE)
    yMax <- 1.0 * max.intensity

    op<-par(mar=c(5,5,5,5), cex=0.75)
    plot(spec$mZ,spec$intensity,
        xlab='m/z',
        ylab='Intensity',
        type='h',
        ylim=c(0,yMax),
        main=peptideSequence,
        axes='F',
        xlim=c(min(m$fragmentIons), max(m$fragmentIons))
    ) 

    n<-nchar(peptideSequence)
    #axis(3,m$fragmentIons$b, substring(peptideSequence, 1:n, 1:n))
    LABEL.abc<-(abs(m$mZ.Da.error) < 0.6) & (regexpr("[abc].*", m$label) > 0)
    LABEL.xyz<-(abs(m$mZ.Da.error) < 0.6) & (regexpr("[xyz].*", m$label) > 0)

    axis(1,spec$mZ[m$idx[LABEL.abc]], m$label[LABEL.abc])
    axis(2)
    axis(3,spec$mZ[m$idx[LABEL.xyz]], m$label[LABEL.xyz])


    peakplot.label(spec=spec, match=m, yMax=yMax)
    axis(4,seq(0,yMax,length=6), seq(0,100,length=6))
    box()
    par(op)

    return(m)
}                      
 

peakplot.pie<-function(spec, match){ 

    LABEL.abc<-abs(match$mZ.Da.error < 0.6) & (regexpr("[abc].*", match$label) > 0)
    LABEL.xyz<-abs(match$mZ.Da.error < 0.6) & (regexpr("[xyz].*", match$label) > 0)


    i.abc<-spec$intensity[match$idx[LABEL.abc]]
    i.xyz<-spec$intensity[match$idx[LABEL.xyz]]

    l.abc<-match$label[LABEL.abc]
    l.xyz<-match$label[LABEL.xyz]

    i.rest<-sum(spec$intensity)-sum(i.abc)-sum(i.xyz)

    pie(c(i.abc,i.xyz,i.rest), c(l.abc, l.xyz, "rest"), col=c(rep("blue",length(i.abc)), rep("grey",length(i.abc)), "white"))
}                      

#peakplot('HTLNQIDSVK', spec)
# source("/home/cp/__svncheckouts/R/protViz/R/peakplot.R")
