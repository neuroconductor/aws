\name{extract-methods}
\docType{methods}
\alias{extract-methods}
\alias{extract,ANY-method}
\alias{extract,aws-method}
\alias{extract,awssegment-method}
\alias{extract,ICIsmooth-method}
\alias{extract,kernsm-method}
\title{ Methods for Function \code{extract} in Package \pkg{aws}}
\description{
  The method extract and/or compute specified statistics from object of class 
\code{"aws"}, \code{"awssegment"}, \code{ICIsmooth} and \code{"kernsm"}. }
\usage{
  \S4method{extract}{aws}(x, what="y")
  \S4method{extract}{awssegment}(x, what="y")
  \S4method{extract}{ICIsmooth}(x, what="y")
  \S4method{extract}{kernsm}(x, what="y")
}
\arguments{
\item{x}{ object }
\item{what}{Statistics to extract, defaults to \code{what="y"} corresponding to the 
          original data (response variable). Alternatives are
          \code{what="yhat"} for the smoothed response, \code{what="vhat"} for the 
          estimated variance of the smoothed response, \code{what="sigma2"} for the 
          estimated error variance of the original data, \code{what="vred"} for the variance
          reduction achieved and in case of \code{signature(x = "ICIsmooth")} 
          \code{what="hbest"}
          for the selected bandwidth. A vector of any of these choices may be provided.}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{
Returns a message that method \code{extract}  is not defined.
}

\item{\code{signature(x = "aws")}}{
Returns a list with components containing the requested statistics. Component names correspond to
\code{tolower(what)}  
}

\item{\code{signature(x = "awssegment")}}{
Returns a list with components containing the requested statistics. Component names correspond to
\code{tolower(what)}  
}

\item{\code{signature(x = "ICIsmooth")}}{
Returns a list with components containing the requested statistics. Component names correspond to
\code{tolower(what)}.  
}

\item{\code{signature(x = "kernsm")}}{
Returns a list with components containing the requested statistics. Component names correspond to
\code{tolower(what)}.  
}
}}
\keyword{methods}
\keyword{manip}
