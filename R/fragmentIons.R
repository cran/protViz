#R

defaultIons<-function(fi){
    Hydrogen <- 1.007825
    Oxygen <- 15.994915
    Nitrogen <- 14.003074

    y_0 <- fi$y - Oxygen - Hydrogen - Hydrogen
    c<- fi$b + (Nitrogen + (3 * Hydrogen))
    z<- fi$y - (Nitrogen + (3 * Hydrogen))
    
    return(cbind(y_0,c,z))
}


fragmentIons<-function(sequence, FUN=defaultIons) {
    if (!is.character(sequence)) {
        stop ("argument x must be a character")
    }else{
                                                    
        FUN<-match.fun(FUN)

        R<-list()
        pim<-parentIonMass(sequence)

        Hydrogen <- 1.007825
        Carbon <- 12.000000
        Nitrogen <- 14.003074
        Oxygen <- 15.994915
        Electron <- 0.000549
        C_term <- 17.002740
        N_term <- 1.007825



        for (i in 1:length(sequence)){
            pepseq<-sequence[i]
            pepseq.pim<-pim[i]

            out <- .C("computeFragmentIons",
                n=as.integer(nchar(pepseq)),
                pepSeq=as.character(pepseq),
                pim=as.double(pepseq.pim),
                b=as.double(rep(0.0,nchar(pepseq))),
                y=as.double(rep(0.0,nchar(pepseq))))


            fi<-as.data.frame(cbind(b=out$b, y=out$y))
            fi<-cbind(fi,as.data.frame(ff<-FUN(fi)))

            R[[length(R)+1]] <- fi
        }
    }
    return(R)
}
