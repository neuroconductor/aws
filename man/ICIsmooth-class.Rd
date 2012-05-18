\name{ICIsmooth-class}
\Rdversion{1.1}
\docType{class}
\alias{ICIsmooth-class}

\title{Class \code{"ICIsmooth"}}
\description{The \code{"ICIsmooth"} class is 
used for objects obtained by functions \code{ICIsmooth} and \code{ICIcombined}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ICIsmooth", ...)}
or by functions \code{ICIsmooth} and \code{ICIcombined}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}, usually empty. }
    \item{\code{y}:}{Object of class \code{"array"} containing the original (response) data }
    \item{\code{dy}:}{Object of class \code{"numeric"} dimension attribute of \code{y} }
    \item{\code{x}:}{Object of class \code{"numeric"} if provided the design points}
    \item{\code{hmax}:}{Object of class \code{"numeric"} maximum bandwidth }
    \item{\code{hinc}:}{Object of class \code{"numeric"} initial bandwidth }
    \item{\code{thresh}:}{Object of class \code{"numeric"} threshold used for bandwidth selection }
    \item{\code{kern}:}{Object of class \code{"character"} kernel, can be any of
     \code{c("Gaussian","Uniform","Triangle","Epanechnicov","Biweight","Triweight")}.
Defaults to \code{kern="Gaussian"}. }
    \item{\code{m}:}{Object of class \code{"integer"} vector of length \code{length(dy)} 
       determining the order of derivatives specified for the coordinate directios.}
    \item{\code{nsector}:}{Object of class \code{"integer"} number of sectors to use. }
    \item{\code{sector}:}{Object of class \code{"integer"} sector used. }
    \item{\code{symmetric}:}{Object of class \code{"logical"} sectors are symmetric with respect to the origin. }
    \item{\code{yhat}:}{Object of class \code{"array"} smoothed response variable }
    \item{\code{vhat}:}{Object of class \code{"array"} estimated variance of smoothed response variable }
    \item{\code{hbest}:}{Object of class \code{"array"} selected bandwidth(s)) }
    \item{\code{sigma}:}{Object of class \code{"numeric"} estimated standard deviation of errors in \code{y} }
    \item{\code{call}:}{Object of class  \code{"call"} that created the object.  }
  }
}
\section{Methods}{
  \describe{
    \item{extract}{\code{signature(x = "ICIsmooth")}: ... }
    \item{risk}{\code{signature(y = "ICIsmooth")}: ... }
    \item{plot}{ Method for Function `plot' in Package `aws'. }
    \item{show}{ Method for Function `show' in Package `aws'. }
    \item{print}{ Method for Function `print' in Package `aws'. }
    \item{summary}{ Method for Function `summary' in Package `aws'. }
	 }
}
\references{
V. Katkovnik, K. Egiazarian and J. Astola, 
\emph{Local Approximation Techniques in Signal And Image Processing},
 SPIE Society of Photo-Optical Instrumentation Engin., 2006, PM157
}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{
\code{\link{ICIsmooth}}, \code{\link{ICIcombined}}, \code{\link{kernsm}}, \code{\link{aws}}
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("ICIsmooth")
}
\keyword{classes}
