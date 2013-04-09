\name{genMod}

\alias{genMod}

\title{Generates all possible modification of a given peptide seqence and modification pattern}

\description{

    This function can be used to screen precursor masses for meaningful mass shifts which could correspond to post translational modifications.
    We suggest to only use a maxiumum of a handful of mass shifts.
}

\usage{

    genMod(peptideSequences, modificationPattern, nModification=2)

}

\arguments{

  \item{peptideSequences}{xxx xxx xx xxxx}

  \item{modificationPattern}{xxx xxx xxx}

  \item{nModification}{number of maximal modifications.}

}

\details{

  t.b.d. xxx xxxxx 
}

\author{

   Hubert Rehrauer and Christian Panse 2012

}

\examples{
    ptm.0<-cbind(AA="-",
        mono=0.0, avg=0.0, desc="unmodified", unimodAccID=NA)

    ptm.616<-cbind(AA='S',
        mono=-27.010899, avg=NA, desc="Substituition",
            unimodAccID=616)

    ptm.651<-cbind(AA='N',
        mono=27.010899, avg=NA, desc="Substituition",
            unimodAccID=651)


    m<-as.data.frame(rbind(ptm.0, ptm.616, ptm.651))

    genMod(c('TAFDEAIAELDTLNEESYK','TAFDEAIAELDTLSEESYK'), m$AA)

}
