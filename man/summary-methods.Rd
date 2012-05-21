\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,ANY-method}
\alias{summary,aws-method}
\alias{summary,awssegment-method}
\alias{summary,ICIsmooth-method}
\alias{summary,kernsm-method}
\title{Methods for Function `summary' from package 'base' in Package `aws' }
\description{
The method provides summary information for objects of class \code{"aws"}. }
\arguments{
  \item{object}{ Object of class \code{"dti"}, \code{"dtiData"}, \code{"dtiTensor"}, \code{"dwiMixtensor"}, \code{"dtiIndices"},  \code{"dwiQball"} or \code{"dwiFiber"}. }
  \item{\dots}{ Additional arguments in \code{\dots} are passed to function \code{quantile}, e.g. argument \code{probs} may be specified here.}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ANY")}}{
Generic function: see \code{\link[base]{summary}}.
}

\item{\code{signature(object = "aws")}}{
The function provides summary information for objects of class \code{"aws"}
}

\item{\code{signature(object = "awssegment")}}{
The function provides summary information for objects of class \code{"awssegment"}
}

\item{\code{signature(object = "ICIsmooth")}}{
The function provides summary information for objects of class \code{"ICIsmooth"}

}

\item{\code{signature(object = "kernsm")}}{
The function provides summary information for objects of class \code{"kernsm"}

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
\keyword{ utilities }
