\name{swath}

\alias{genSwathIonLib}

\title{SWATH exporter}

\description{
    This function generates an ion library for SWATH analysis. It reads mascot dat-files and can be used with an Ion score filter.
}

\usage{

    genSwathIonLib(data, 
        mascotIonScoreCutOFF=20, 
        proteinIDPattern='', 
        outputFile="genSwathIonLib.txt")

}

\arguments{

  \item{data}{data set containing mZ and peptide sequence.}

  \item{mascotIonScoreCutOFF}{a value for filtering the specs.}

  \item{proteinIDPattern}{a filter for protein.}

  \item{outputFile}{the output file name.}

}

\details{

    The \code{parentIonMass} function reqires one argument \code{secuence}.

}

\author{

   Christian Panse 2012

}

