\name{show-methods}
\docType{methods}
\alias{show-methods}
\alias{show,ANY-method}
\alias{show,aws-method}
\alias{show,awssegment-method}
\alias{show,ICIsmooth-method}
\alias{show,kernsm-method}
\title{Methods for Function `show' in Package `aws'}
\description{
The function provides information on data dimensions, data source and 
existing slot-names for objects of class \code{"aws"}, \code{"awssegment"},
 \code{"ICIsmooth"} and \code{"kernsm"} in package \pkg{aws}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ANY")}}{
Generic function. 
}

\item{\code{signature(object = "aws")}}{
Provide information on data dimensions, data source and existing slot-names for objects of class \code{"dti"}
                    and classes that extent \code{"aws"}.}

\item{\code{signature(object = "awssegment")}}{
Provide information on data dimensions, data source and existing slot-names for objects of class \code{"dti"}
                    and classes that extent \code{"awssegment"}.}


\item{\code{signature(object = "ICIsmooth")}}{
Provide information on data dimensions, data source and existing slot-names for objects of class \code{"dti"}
                    and classes that extent \code{"ICIsmooth"}.}

\item{\code{signature(object = "kernsm")}}{
Provide information on data dimensions, data source and existing slot-names for objects of class \code{"dti"}
                    and classes that extent \code{"kernsm"}.}

}}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\linkS4class{aws}}, 
  \code{\linkS4class{awssegment}}, 
  \code{\linkS4class{ICIsmooth}}
  \code{\linkS4class{kernsm}}
}
\keyword{methods}
\keyword{ utiities }
