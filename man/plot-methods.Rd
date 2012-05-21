\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY-method}
\alias{plot,aws-method}
\alias{plot,awssegment-method}
\alias{plot,ICIsmooth-method}
\alias{plot,kernsm-method}
\title{Methods for Function `plot' from package 'graphics' in Package `aws'}
\description{
Visualization of objects of class \code{"aws"}, \code{"awsswgment"}, 
\code{"kernsm"} and \code{"ICIsmooth"}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{
Generic function: see \code{\link[graphics]{plot}}. 
}

\item{\code{signature(x = "aws")}}{
Visualization of objects of class \code{"aws"}
}

\item{\code{signature(x = "awssegment")}}{
Visualization of objects of class \code{"awssegment"}

}

\item{\code{signature(x = "ICIsmooth")}}{
Visualization of objects of class \code{"ICIsmooth"}
}

\item{\code{signature(x = "kernsm")}}{
Visualization of objects of class \code{"kernsm"}
}
}}
\author{ 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\linkS4class{aws}}, 
  \code{\linkS4class{awssegment}}, 
  \code{\linkS4class{ICIsmooth}}
  \code{\linkS4class{kernsm}}
}
\keyword{methods}
\keyword{ hplot }
