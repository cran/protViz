#R

# TODO 
# compute score by sum of error div. by number of hits

psm<-function(sequence, spec, FUN=defaultIons,
    plot=TRUE, 
    fi=fragmentIons(sequence, FUN=FUN)[[1]],
    fragmentIonError=0.6) { 

    n<-nchar(sequence)

    pim<-fi$y[nrow(fi)]

    # consider only b and y ions
    # by.mZ<-c(fi$b, fi$y)
    # by.label<-c(paste("b",1:n,sep=''), paste("y",1:n,sep=''))

    by.mZ<-numeric()
    by.label<-character()
    fi.names<-names(fi)

    for (i in 1:ncol(fi)){
        by.mZ <- c(by.mZ, fi[,i])
        by.label <- c(by.label, paste(fi.names[i],1:n,sep=''))
    }


    # this is an alternative to the NN c funktion
    # ansi-c bsearch  does not work since it is searching for an 
    # exact match
    #idx<-findInterval(by.mZ, spec$mZ)
    #mZ.error<-spec$mZ[idx+1]-by.mZ

    out <- .C("findNN",
        nbyion=as.integer(length(by.mZ)),
        nmZ=as.integer(length(spec$mZ)),
        byion=as.double(by.mZ),
        mZ=as.double(spec$mZ),
        NN=as.integer(rep(-1, length(by.mZ))))


    mZ.error<-spec$mZ[out$NN+1] - by.mZ

    if (plot == TRUE){
        plot(mZ.error[mZ.error.idx<-order(mZ.error)],
            main=paste("Error of", sequence, "(parent ion mass =", round(pim,2) ,"Da)"),
            ylim=c(-5*fragmentIonError, 5*fragmentIonError),
            pch='o',
            cex=0.5,
            sub=paste('The error cut-off is', 
                fragmentIonError, 'Da (grey line).'),
            )

        abline(h=fragmentIonError,col='grey')
        abline(h=-fragmentIonError,col='grey')
        abline(h=0,col='grey',lwd=2)

        text(1:length(by.label), 
            mZ.error[mZ.error.idx],  
            by.label[mZ.error.idx],
            cex=0.75,pos=3) 

        hits=(abs(mZ.error) < fragmentIonError)
        nHits<-sum(hits)

        sumMZerror=round(sum(abs(mZ.error[hits])),2)

        avgMZerror=round(sumMZerror / nHits, 2)
        cover=round(nHits/(nrow(fi)*ncol(fi)),2)

        legend("topleft", paste(c('nHits','sumMZerror','avgMZerror','cover'),
            as.character(c(nHits, sumMZerror, avgMZerror, cover)),sep='=')) 
    }


    return (list(mZ.Da.error=mZ.error, 
        mZ.ppm.error=1E+6*mZ.error/by.mZ,
        idx=out$NN+1,
        label=by.label, 
        score=-1, 
        sequence=sequence,
        fragmentIons=fi))
}

