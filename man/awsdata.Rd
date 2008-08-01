\name{awsdata}
\alias{awsdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract information from an object of class aws}
\description{
Extract data and estimates  from an object of class aws
}
\usage{
awsdata(awsobj, what)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{awsobj}{ an object of class \code{aws} }
  \item{what}{ can be "data" (extracts observed response),
"theta" (estimated parameters),
"est" (estimated regression function), 
"var" (approx. variance of estimated regression function), 
"sd" (approx. standard deviation of estimated regression function),
"sigma2" (error variance),
"mae" (mean absolute error for each iteration step, if available),
"ni" (number of observations per bin),
"mask" (logical indicator for bins where the regression function is estimated).
 }
}
\details{
The returned object is formatted as an array if appropriate. 
The returned object may be NULL if the information is not available.
}
\value{
an vector or array containing the specified information.
}
\references{ 
Joerg Polzehl, Vladimir Spokoiny,  Adaptive Weights Smoothing with applications to image restoration,  J. R. Stat. Soc. Ser. B Stat. Methodol.  62 ,  (2000)  , pp. 335--354

Joerg Polzehl, Vladimir Spokoiny,  Propagation-separation approach for local likelihood estimation,  Probab. Theory Related Fields  135 (3),  (2006)  , pp. 335--362.

Joerg Polzehl, Vladimir Spokoiny, in V. Chen, C.; Härdle, W. & Unwin, A. (ed.) Handbook of Data Visualization Structural adaptive smoothing by propagation-separation methods Springer-Verlag, 2008, 471-492

 }
\author{ Joerg Polzehl \email{polzehl@wias-berlin.de}  }
\seealso{\code{link{awsdata}},\code{\link{aws}}, \code{\link{aws.irreg}}}
\examples{
library(aws)
# 1D local constant smoothing
\dontrun{demo(aws_ex1)}
\dontrun{demo(aws_ex2)}
# 2D local constant smoothing
\dontrun{demo(aws_ex3)}
# 1D local polynomial smoothing
\dontrun{demo(lpaws_ex1)}
# 2D local polynomial smoothing
\dontrun{demo(lpaws_ex2)}
# 1D irregular design
\dontrun{demo(irreg_ex1)}
# 2D irregular design 
\dontrun{demo(irreg_ex2)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ smooth }
\keyword{ manip }
