\name{pressureProfilePlot}

\alias{pressureProfilePlot}

\title{plotting pressure profile data}

\description{
    plots the pressure profiles data on a scatter plot \code{Pc} versus \code{time} grouped by time range. 
}

\usage{

    pressureProfilePlot(data)

}

\arguments{

  \item{data}{A data set like \code{pressureProfile}. The data must have an attribute \code{file}, \code{time}, and \code{Pc}.}

}

\details{

    t.b.d.

}

\author{

    Bernd Roschitzki, Christian Trachsel, Christian Panse 2012

}

\examples{

    data(pressureProfile)
    pressureProfilePlot(pressureProfile)

}
