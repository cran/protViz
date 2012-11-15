\name{peakplot}

\alias{peakplot}
\alias{peakplot.label}
\alias{peakplot.putlabel}
\alias{peakplot.pie}

\title{Labelling of Peptide Fragment Mass Spectra}

\description{

     \code{peakplot} labels a mass spectrum  from a peptide sequence  
     assignment by the \code{matching} function with the appropriate 
     fragment ion labels.

}

\usage{

    peakplot(peptideSequence, spec, FUN=defaultIons) 

}

\arguments{

  \item{peptideSequence}{peptide sequence encoded as character sequence using the 20 amino acid letters.}
  \item{spec}{ a tandem mass spec data structure which is a list having two eq sized vectors calles mZ and intensity. mZ values are sorted.}
  \item{FUN}{ the function to be applied to compute further ions. If no function is assigned \code{fragmentIons} will
    use defaultIons.}

}

\details{
\code{peakplot} computes the in-silico fragment ion and matches them again the MS2.
If the mass error between the in-silico peak and the meassured one is below a 0.6 Da the
match is considered has hit.

The major objective of the labelling is avoidance of overlapping
labels for which purpose PeakPlot applies filtering. A label is only drawn if
the corresponding ion count of the m/z peak is higher than a given
threshold. Experience from several hundred annotations shows that the 0.9
percentile is a good cut off value. To most efficiently use the limited screen
and printout space and to ensure that labels representing important local
peaks are also considered for drawing, 
\code{peakplot} divides the display space
into a number of bins depending on the peptide sequence length along the
m/z axis. From these bins the top n labels are ordered according to
abundance. For the visual clustering effect the abc and xyz ions are drawn
on different y-axis levels using different colours. Ion types considered for
labelling is dependent on the instrument setting applied during the initial
search. 
}

\author{

   Bertran Gerrits, Christian Panse 2006-2012; first time used in 
    \url{http://www.ncbi.nlm.nih.gov/pubmed/17114649}.

}


\examples{
    spec<-list(scans=1138,
        title="178: (rt=22.3807) [20080816_23_fetuin_160.RAW]",
        rtinseconds=1342.8402,
        charge=2,
        mZ=c(195.139940, 221.211970, 239.251780, 290.221750, 
    316.300770, 333.300050, 352.258420, 448.384360, 466.348830, 
    496.207570, 509.565910, 538.458310, 547.253380, 556.173940, 
    560.358050, 569.122080, 594.435500, 689.536940, 707.624790, 
    803.509240, 804.528220, 822.528020, 891.631250, 909.544400, 
    916.631600, 973.702160, 990.594520, 999.430580, 1008.583600, 
    1017.692500, 1027.605900),
        intensity=c(931.8, 322.5, 5045, 733.9, 588.8, 9186, 604.6,
    1593, 531.8, 520.4, 976.4, 410.5, 2756, 2279, 5819, 2.679e+05,
    1267, 1542, 979.2, 9577, 3283, 9441, 1520, 1310, 1.8e+04,
    587.5, 2685, 671.7, 3734, 8266, 3309))

    peakplot('HTLNQIDSVK', spec)

}

