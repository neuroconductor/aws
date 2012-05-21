\name{print-methods}
\docType{methods}
\alias{print-methods}
\alias{print,ANY-method}
\alias{print,aws-method}
\alias{print,awssegment-method}
\alias{print,ICIsmooth-method}
\alias{print,kernsm-method}
\title{Methods for Function `print' from package 'base' in Package `aws'}
\description{
  The function provides information on data dimensions, creation of the object and 
  existing slot-names for objects of class \code{"aws"}, \code{"awssegment"},
 \code{"ICIsmooth"} and  \code{"kernsm"}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{
Generic function: see \code{\link[base]{print}}.
}

\item{\code{signature(x = "aws")}}{
Provide information on data dimensions, creation of the object and 
  existing slot-names for objects of class \code{"aws"}
}

\item{\code{signature(x = "awssegment")}}{
Provide information on data dimensions, creation of the object and 
  existing slot-names for objects of class \code{"awssegment"}
}

\item{\code{signature(x = "ICIsmooth")}}{
Provide information on data dimensions, creation of the object and 
  existing slot-names for objects of class \code{"ICIsmooth"}
}

\item{\code{signature(x = "kernsm")}}{
Provide information on data dimensions, creation of the object and 
  existing slot-names for objects of class \code{"kernsm"}
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
