\docType{data}

\name{HexNAc}

\alias{HexNAc}

\title{HexNAc - Analysis of N-HexNAc glycopeptides by LC-MS/MS, using HCD and ETD
fragmentation techniques}

\format{A data set consists of eleven tandem mass spectra.}

\description{

A sample containing enriched glycopeptides from yeast was
analysed by LC-MS/MS on a LTQ-Orbitrap Velos instrument using sequential
HCD and ETD fragmentation techniques.
The dataset contains only 11 tandem mass spectra extracted from this
experiment. Three pairs of HCD/ETD spectra correspond to peptides carrying
N-HexNAc modification identified with high confidence (6 spectra in total).
One pair of HCD/ETD spectra corresponds to peptides carrying N-HexNAc
modification, identified with very low confidence (2 spectra in total). The
remaining 3 spectra are from unmodified peptides.

}

\author{

  Paolo Nanni, Christian Panse, 2013

}

\examples{

    data(HexNAc)
    HexNAc[[1]]

    plot(HexNAc[[1]]$mZ, HexNAc[[1]]$intensity, type='h')

}

\keyword{datasets}
