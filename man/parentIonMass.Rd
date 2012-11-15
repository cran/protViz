\name{parentIonMass}

\alias{parentIonMass}
\alias{pim}
\alias{pepmass}

\title{Compute Parent Ion Mass of a Peptide Sequence}

\description{
    It computes the parent ion mass of a vector of peptide sequences. The weights of the amino acids and the C-Term are hard coded in the C function.
}

\usage{

    parentIonMass(sequence)

}

\arguments{
  \item{sequence}{peptide sequence encoded as character sequence using the 20 amino acid letters.}
}

\details{

    The \code{parentIonMass} function reqires one argument \code{secuence}.

}

\author{

   Christian Panse 2006

}
\examples{
        
        fetuin<-list( sequence=c('MK', 'SFVLLFCLAQLWGCHSIPLDPVAGYK',
        'EPACDDPDTEQAALAAVDYINK',
        'HLPR', 'GYK', 'HTLNQIDSVK', 'VWPR',
        'RPTGEVYDIEIDTLETTCHVLDPTPLANCSVR',
        'QQTQHAVEGDCDIHVLK', 'QDGQFSVLFTK',
        'CDSSPDSAEDVR', 'K', 'LCPDCPLLAPLNDSR',
        'VVHAVEVALATFNAESNGSYLQLVEISR',
        'AQFVPLPVSVSVEFAVAATDCIAK',
        'EVVDPTK', 'CNLLAEK', 'QYGFCK',
        'GSVIQK', 'ALGGEDVR',
        'VTCTLFQTQPVIPQPQPDGAEAEAPSAVPDAAGPTPSAAGPPVASVVVGPSVVAVPLPLHR',
        'AHYDLR', 'HTFSGVASVESSSGEAFHVGK',
        'TPIVGQPSIPGGPVR', 'LCPGR', 'IR', 'YFK', 'I'),
        description=
            "FETUA_BOVIN Alpha-2-HS-glycoprotein ")

        plot(peptideMass<-parentIonMass(sequence=fetuin$sequence), 
            sub=list(fetuin$description, cex=0.75), 
            main="tryptic digested peptides")
}

