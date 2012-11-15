\name{fragmentIons}

\alias{fragmentIons}
\alias{bymatrix}
\alias{ionseries}
\alias{defaultIons}

\title{Compute b and y Fragment Ions of a Peptide Sequence}

\description{

    The function computes the fragment ions of a peptide sequence.

    If multiple peptide sequences are given it returns a list of fragment ion table.
    
}

\usage{

    fragmentIons(sequence, FUN)

}

\arguments{

  \item{sequence}{peptide sequence encoded as character sequence using the 20 amino acid letters.}

  \item{FUN}{the function to be applied to compute further ions. 
  If no function is assigned \code{fragmentIons} will use \code{defaultIons}.}

}

\details{

    The fragment ions of a peptide can be computed following the rules proposed in PMID:6525415.
    Beside the \code{b} and \code{y} ions the \code{FUN} argument of \code{fragmentIons} defines
    which ions are computed.
    the default ions beeing computed are defined in the function \code{defaultIons}.
    The are no limits for defining other forms of fragment ions for ETD (c and z ions) CID (b and y ions).

}

\references{

      \url{http://www.ncbi.nlm.nih.gov/pubmed/6525415}

}

\author{

   Christian Panse, Bertran Gerrits 2006; first used in \code{peakplot}.

}

\examples{
        myIons<-function(fi){
            Hydrogen <- 1.007825
            Oxygen <- 15.994915
            Nitrogen <- 14.003074

            y_0 <- fi$y - Oxygen - Hydrogen - Hydrogen
            c<- fi$b + (Nitrogen + (3 * Hydrogen))
            z<- fi$y - (Nitrogen + (3 * Hydrogen))
            return(cbind(y_0,c,z))
        }
        
        peptides<-c('HTLNQIDSVK', 'ALGGEDVR', 'TPIVGQPSIPGGPVR')

        pim<-parentIonMass(peptides)
        fi<-fragmentIons(peptides)

        par(mfrow=c(3,1)); 
        for (i in 1:length(peptides)){
            plot(0,0,
                xlab='m/Z',
                ylab='',
                xlim=range(c(fi[i][[1]]$b,fi[i][[1]]$y)),
                ylim=c(0,1),
                type='n',
                axes=FALSE,
                sub=paste( pim[i], "Da"));
            box()
            axis(1,fi[i][[1]]$b,round(fi[i][[1]]$b,2))

            pepSeq<-strsplit(peptides[i],"")
            axis(3,fi[i][[1]]$b,pepSeq[[1]])

            abline(v=fi[i][[1]]$b, col='red',lwd=2) 
            abline(v=fi[i][[1]]$c, col='orange') 
            abline(v=fi[i][[1]]$y, col='blue',lwd=2)
            abline(v=fi[i][[1]]$z, col='cyan')
        }
}
