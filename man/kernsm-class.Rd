\name{kernsm-class}
\Rdversion{1.1}
\docType{class}
\alias{kernsm-class}

\title{Class \code{"kernsm"}}
\description{This class refers to objects created by function kernsm. These objects
contain 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("kernsm", ...). they are usually
created by a call to function{kernsm}}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}, usually empty. }
    \item{\code{y}:}{Object of class \code{"array"} containing the response in
          nonparametric regression. The design is assumed to be a 1D, 2D or 3D grid,
          with dimensionality determined by \code{dim(y)}. }
    \item{\code{dy}:}{Object of class \code{"numeric"} containing dim(y). }
    \item{\code{x}:}{Object of class \code{"numeric"} currently not used.  }
    \item{\code{h}:}{Object of class \code{"numeric"} containing the bandwidth employed. }
    \item{\code{kern}:}{Object of class \code{"character"} determining the kernel that was used,
        can be one of \code{c("Gaussian","Uniform","Triangle","Epanechnikov","Biweight","Triweight")} }
    \item{\code{m}:}{Object of class \code{"integer"} with length \code{length(dy)} determining the
        order of derivatives in the corresponding coordinate directions. If \code{m[i6>0]} a
        dirivative kernel derived from \code{kern} has been used for the 
        corresponding coordinate direction. }
    \item{\code{nsector}:}{Object of class \code{"integer"}. If \code{nsector>1} positive weights 
          are restricted to a segment of a circle (1D or 2D only). The segment is given by
          \code{sector}.}
    \item{\code{sector}:}{Object of class \code{"integer"} containing the number of the segment 
          used in case of \code{nsector>1} }
    \item{\code{symmetric}:}{Object of class \code{"logical"} determines if the sector is mirrored at the origin. }
    \item{\code{yhat}:}{Object of class \code{"array"} with same size and dimension as \code{y} providing
             the convolution of \code{y} with the chosen kernel. }
    \item{\code{vred}:}{Object of class \code{"array"} Variance reduction achieved by convolution
          assuming independence. }
    \item{\code{call}:}{Object of class \code{"function"}, call that created the object.}
  }
}
\section{Methods}{
  \describe{
    \item{extract}{\code{signature(x = "aws")}: ... }
    \item{risk}{\code{signature(y = "aws")}: ... }
    \item{plot}{ Method for Function `plot' in Package `aws'. }
    \item{show}{ Method for Function `show' in Package `aws'. }
    \item{print}{ Method for Function `print' in Package `aws'. }
    \item{summary}{ Method for Function `summary' in Package `aws'. }
	 }
}
\author{
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{kernsm}}, 
\code{\link{ICIsmooth}}, 
\code{\link{ICIcombined}}, 
\code{\linkS4class{ICIsmooth}} 
}
\examples{
showClass("kernsm")
}
\keyword{classes}
