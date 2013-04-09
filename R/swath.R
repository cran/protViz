#R

# <cp@fgcz.ethz.ch>
# Wed Nov 28 15:12:28 CET 2012

genSwathIonLib <- function(data, mascotIonScoreCutOFF=20, proteinIDPattern='', outputFile="genSwathIonLib.txt"){

    rr<-numeric()
    n<-length(data)
    idx<-1:n

    for (i in idx){
        #i<-39760
        # TODO: 
        #m<-psm(data[[i]]$peptideSequence, data[[i]], FUN=function(x){})
        #group_id peptide_sequence q1 q3 decoy prec_z frg_type frg_nr frg_z intensity_predicted irt


        if ((!is.na(data[[i]]$mascotScore)) & data[[i]]$mascotScore > mascotIonScoreCutOFF & regexpr(proteinIDPattern, data[[i]]$proteinInformation) > 1){
        fi<-fragmentIons(data[[i]]$peptideSequence)[[1]]
        findNN.idx<-findNN(c(fi$b,fi$y), data[[i]]$mZ)

        m<-2*nchar(data[[i]]$peptideSequence)
        mZ.error<-data[[i]]$mZ[findNN.idx] - c(fi$b,fi$y)

    # prepare table for output
        group_id <- rep(paste(data[[i]]$peptideSequence,".", data[[i]]$charge,";", data[[i]]$pepmass, sep=''), m)
        peptide_sequence <- rep(data[[i]]$peptideSequence, m)
        q1<-rep(data[[i]]$pepmass, m)
        q3<-data[[i]]$mZ[findNN.idx]
        decoy<-rep(0, m)
        prec_z<-rep(data[[i]]$charge, m)
        frg_type<-c(rep('b',length(fi$b)), rep('y',length(fi$y)))
        frg_nr<-rep(1:length(fi$b),2)
        frg_z<-rep(1, m)
        intensity_predicted<-100*round(data[[i]]$intensity[findNN.idx]/max(data[[i]]$intensity[findNN.idx]),2)
        irt<-rep(data[[i]]$rtinseconds/60, m)

        r<-as.data.frame(cbind(group_id,peptide_sequence,q1,q3,decoy,prec_z,frg_type,frg_nr,frg_z,intensity_predicted,irt))


        rr<-rbind(rr,r[mZ.error < 0.1,])
        }
    }
    write.table(rr, file=outputFile, 
        row.names = FALSE,
        col.names = TRUE, 
        quote = FALSE, sep = "\t")
    return(rr)
}

