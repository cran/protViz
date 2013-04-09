#R

defaultIons<-function(fi){
    Hydrogen <- 1.007825
    Oxygen <- 15.994915
    Nitrogen <- 14.003074

    #yo <- fi$y - Oxygen - Hydrogen - Hydrogen
    c<- fi$b + (Nitrogen + (3 * Hydrogen))
    z<- fi$y - (Nitrogen + (3 * Hydrogen))
    
    return(cbind(c,z))
}


fragmentIons<-function(sequence, FUN=defaultIons, modified=numeric(), modification=numeric()) {
    if (!is.character(sequence)) {

        R<-list()

        input.n <- length(sequence)

        out <- .C("_computeFragmentIons", n=input.n, 
            W_=as.double(sequence), 
            b_=as.double(rep(0,input.n)), 
            y_=as.double(rep(0,input.n)))

        fi<-as.data.frame(cbind(b=out$b, y=out$y))

        fi<-cbind(fi,as.data.frame(ff<-FUN(fi)))

        R[[1]] <- fi

    }else if (length(modification) > 1) {
                                                    
        FUN<-match.fun(FUN)

        R<-list()
        pim<-parentIonMass(sequence)

        C_term <- 17.002740
        N_term <- 1.007825
        Oxygen <- 15.994915
        Carbon <- 12.000000
        Hydrogen <- 1.007825
        Nitrogen <- 14.003074
        Electron <- 0.000549

        for (i in 1:length(sequence)){
            input.sequence<-sequence[i]
            input.n<-nchar(input.sequence)
            input.modified <- as.integer(strsplit(modified[i], '')[[1]])
            input.pim<-pim[i]+(sum(as.double(as.character(modification[input.modified+1]))))

            if (input.n != length(input.modified))
                stop (paste("unvalid argument",i,"- number of AA and modification config differ! stop."))

            out <- .C("computeFragmentIonsModification",
                n=as.integer(input.n),
                pepSeq=as.character(input.sequence),
                pim=as.double(input.pim),
                b=as.double(rep(0.0, input.n)),
                y=as.double(rep(0.0, input.n)),
                modified=input.modified,
                modification=as.double(as.character(modification)))

            fi<-as.data.frame(cbind(b=out$b, y=out$y))
            fi<-cbind(fi,as.data.frame(ff<-FUN(fi)))

            R[[length(R)+1]] <- fi
        }
    } else{
        FUN<-match.fun(FUN)

        R<-list()
        pim<-parentIonMass(sequence)

        C_term <- 17.002740
        N_term <- 1.007825
        Oxygen <- 15.994915
        Carbon <- 12.000000
        Hydrogen <- 1.007825
        Nitrogen <- 14.003074
        Electron <- 0.000549

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
