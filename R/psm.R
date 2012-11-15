#R

# TODO 
# compute score by sum of error div. by number of hits

psm<-function(sequence, spec, FUN=defaultIons,
    plot=TRUE, 
    fragmentIonError=0.6) { 

    n<-nchar(sequence)
    pim<-parentIonMass(sequence)

    fi<-fragmentIons(sequence, FUN=FUN)

    by.mZ<-c(fi[[1]]$b, fi[[1]]$y)
    by.label<-c(paste("b",1:n,sep=''), paste("y",1:n,sep=''))

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
            main=paste("Error of", sequence, "(patent ion mass =", pim ,"Da)"),
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
    }

    return (list(mZ.Da.error=mZ.error, 
        mZ.ppm.error=1E+6*mZ.error/by.mZ,
        idx=out$NN+1,
        label=by.label, 
        score=-1, 
        sequence=sequence,
        fragmentIons=fi[[1]]))
}
