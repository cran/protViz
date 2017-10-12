\name{as.psmSet}
\alias{as.psmSet}
\alias{is.psmSet}
\alias{as.psmSet.mascot}
\alias{as.psmSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  psmSet - a set of peptide spectrum matches
}
\description{
defines a class of peptide spectrum matches.
}
\usage{

is.psmSet(object)

 \method{as.psmSet}{mascot}(object, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  an \code{psmSet} S3 class object
}
  \item{\dots}{
...
}
}
\details{
the \code{is.psmSet} method checks if the as input given objects fullfills the 
proterties to be an \code{psmSet} object.

\code{as.psmSet} transformas an object into a \code{psmSet} object.
}
\value{
while \code{is.psmSet} returns \code{TRUE} or \code{FALSE}, 
\code{as.psmSet} returns an instance of \code{psmSet} or \code{NULL}.
}
\author{
Christian Panse
}
