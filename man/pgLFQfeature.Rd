\docType{data}
\name{pgLFQfeature}

\alias{pgLFQfeature}

\title{pgLFQfeature - A data set with a featuremap export ProgenesisLCMS}

\format{A data object from ProgenesisImporter.R}

\description{
This data set contains a fraction of the most abundant features from 24 LC-MSMS runs measured on 
an Orbitrap-Velos. The data contains digested proteins from human HeLa cells infected with Shigella bacteria grown on a time course.
It is structured in a way, that 6 biological replicates from 4 conditions are measured (Not_infected, Infected_1hr, Infected_2hr, Infected_3hr). 
The 24 LC-MSMS runs are aligned with each other and a so called mastermap (feature map) is generated where "normalized volumes" on MS1 are extracted for all features in the respective LC-MSMS runs. Some of the features are further annotated with peptide sequences and protein accessions (using the Mascot search algorithm and an ion score cutoff of 25).



More information on the commercial software can be found here:
\url{http://www.nonlinear.com/products/progenesis/lc-ms/overview/}.

}

\author{

   Christian Panse, Jonas Grossmann 2012

}

\examples{
data(pgLFQfeature)
op<-par(mfrow=c(1,1),mar=c(18,18,4,1),cex=0.5)
samples<-names(pgLFQfeature$"Normalized abundance")
image(cor(asinh(pgLFQfeature$"Normalized abundance")),
    col=gray(seq(0,1,length=20)),
    main='pgLFQfeature correlation',
    axes=FALSE)

axis(1,at=seq(from=0, to=1, 
    length.out=length(samples)), 
    labels=samples, las=2)

axis(2,at=seq(from=0, to=1, 
    length.out=length(samples)), labels=samples, las=2)
par(op)

}

\references{

      \url{http://www.ncbi.nlm.nih.gov/pubmed/20576481}

}

\keyword{datasets}
