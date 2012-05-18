\name{awssegment-class}
\Rdversion{1.1}
\docType{class}
\alias{awssegment-class}

\title{Class \code{"awssegment"}}
\description{
The \code{"aws"} class is 
used for objects obtained by functions \code{aws.segment}
}
\section{Objects from the Class}{
Objects are created by calls to functions \code{aws.segment}
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}, usually empty. }
    \item{\code{y}:}{Object of class \code{"array"} containing the original (response) data }
    \item{\code{dy}:}{Object of class \code{"numeric"} dimension attribute of \code{y} }
    \item{\code{x}:}{Object of class \code{"numeric"} if provided the design points}
    \item{\code{ni}:}{Object of class \code{"numeric"} sum of weights used in final estimate }
    \item{\code{mask}:}{Object of class \code{"logical"} mask of design points where computations are performed }
    \item{\code{segment}:}{Object of class \code{"array"} segmentation results 
            (3 segments coded by \code{c(-1, 0, 1)})}
    \item{\code{level}:}{Object of class \code{"numeric"} center of segment \code{0} }
    \item{\code{delta}:}{Object of class \code{"numeric"} half width of segment \code{0} }
    \item{\code{theta}:}{Object of class \code{"array"} ~~ }
    \item{\code{theta}:}{Object of class \code{"array"} containes the smoothed object and in case
              of function \code{lpaws} its derivatives up to the specified degree. 
              Dimension is \code{dim(theta)=c(dy,p)}}
    \item{\code{mae}:}{Object of class \code{"numeric"} Mean absolute error with respect to 
        array in argument \code{u} if provided. }
    \item{\code{var}:}{Object of class \code{"numeric"} pointwise variance of 
            \code{theta[...,1]}}
    \item{\code{xmin}:}{Object of class \code{"numeric"} not used }
    \item{\code{xmax}:}{Object of class \code{"numeric"} not used}
    \item{\code{wghts}:}{Object of class \code{"numeric"} weights used in location penalty for 
            different coordinate directions }
    \item{\code{degree}:}{not used }
    \item{\code{hmax}:}{Object of class \code{"numeric"} maximal bandwidth }
    \item{\code{sigma2}:}{Object of class \code{"numeric"} estimated error variance }
    \item{\code{scorr}:}{Object of class \code{"numeric"} estimated spatial correlation }
    \item{\code{family}:}{Object of class \code{"character"} distribution of \code{y},
          can be any of \code{c("Gaussian","Bernoulli","Poisson","Exponential",
           "Volatility","Variance")} }
    \item{\code{shape}:}{Object of class \code{"numeric"} possible shape parameter of distribution of \code{y}}
    \item{\code{lkern}:}{Object of class \code{"integer"} location kernel, can be 
       any of \code{c("Triangle","Quadratic","Cubic","Plateau","Gaussian")}, defauts to 
       \code{"Triangle"}}
    \item{\code{lambda}:}{Object of class \code{"numeric"} scale parameter used in adaptation}
    \item{\code{ladjust}:}{Object of class \code{"numeric"} factor to adjust scale parameter with respect to its
        predetermined default. }
    \item{\code{aws}:}{Object of class \code{"logical"} Adaptation by Propagation-Separation  }
    \item{\code{memory}:}{Object of class \code{"logical"} Adaptation by Stagewise Aggregation }
    \item{\code{homogen}:}{Object of class \code{"logical"} detect regions of homogeneity (used to speed up
              the calculations) currently \code{FALSE}}
    \item{\code{earlystop}:}{Object of class \code{"logical"} currently \code{FALSE} }
    \item{\code{varmodel}:}{Object of class \code{"character"} variance model used  
        currently \code{"Gaussian"} }
    \item{\code{vcoef}:}{Object of class \code{"numeric"} contains \code{NULL} }
    \item{\code{call}:}{Object of class  \code{"call"} that created the object. }
  }
}
\section{Methods}{
  \describe{
    \item{extract}{\code{signature(x = "awssegment")}: ... }
    \item{plot}{\code{signature(x = "awssegment")}: ... }
    \item{print}{\code{signature(x = "awssegment")}: ... }
    \item{risk}{\code{signature(y = "awssegment")}: ... }
    \item{show}{\code{signature(object = "awssegment")}: ... }
    \item{summary}{\code{signature(object = "awssegment")}: ... }
	 }
}
\author{
Joerg Polzehl, \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{aws.segment}}
}
\examples{
showClass("awssegment")
}
\keyword{classes}
