\name{aws-class}
\Rdversion{1.1}
\docType{class}
\alias{aws-class}

\title{Class \code{"aws"}}
\description{
The \code{"aws"} class is
used for objects obtained by functions \code{aws}, \code{lpaws}, \code{aws.irreg} and \code{aws.gaussian}.}
\section{Objects from the Class}{
Objects are created by calls to functions \code{aws}, \code{lpaws}, \code{aws.irreg} and \code{aws.gaussian}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}, usually empty. }
    \item{\code{y}:}{Object of class \code{"array"} containing the original (response) data }
    \item{\code{dy}:}{Object of class \code{"numeric"} dimension attribute of \code{y} }
		\item{\code{nvec}:}{Object of class \code{"integer"} leading dimension of \code{y} in vector valued data.}
    \item{\code{x}:}{Object of class \code{"numeric"} if provided the design points}
    \item{\code{ni}:}{Object of class \code{"numeric"} sum of weights used in final estimate }
    \item{\code{mask}:}{Object of class \code{"logical"} mask of design points where computations are performed }
    \item{\code{theta}:}{Object of class \code{"array"} containes the smoothed object and in case
              of function \code{lpaws} its derivatives up to the specified degree.
              Dimension is \code{dim(theta)=c(dy,p)}}
    \item{\code{hseq}:}{Sequence of bandwidths employed.}
    \item{\code{mae}:}{Object of class \code{"numeric"} Mean absolute error with respect to
        array in argument \code{u} if provided. }
    \item{\code{psnr}:}{Object of class \code{"numeric"} Peak Signal to Noise Ratio (PSNR) with respect to
        array in argument \code{u} if provided. }
    \item{\code{var}:}{Object of class \code{"numeric"} pointwise variance of
            \code{theta[...,1]}}
    \item{\code{xmin}:}{Object of class \code{"numeric"} min of \code{x} in case of irregular design }
    \item{\code{xmax}:}{Object of class \code{"numeric"} max of \code{x} in case of irregular design}
    \item{\code{wghts}:}{Object of class \code{"numeric"} weights used in location penalty for
            different coordinate directions, corresponds to ratios of distances in coordinate directions 2 and 3 to
            and distance in coordinate direction 1.}
    \item{\code{degree}:}{Object of class \code{"integer"} degree of local polynomials used in
         function \code{lpaws} }
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
              the calculations) }
    \item{\code{earlystop}:}{Object of class \code{"logical"} further speedup in function \code{lpaws}
          estimates are fixed if sum of weigths does not increase with iterations. }
    \item{\code{varmodel}:}{Object of class \code{"character"} variance model used in
        function \code{aws.gaussian} }
    \item{\code{vcoef}:}{Object of class \code{"numeric"} estimates variance parameters
        in function \code{aws.gaussian} }
    \item{\code{call}:}{Object of class  \code{"call"} that created the object. }
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
\references{
Joerg Polzehl, Vladimir Spokoiny,  Adaptive Weights Smoothing with applications to image restoration,  J. R. Stat. Soc. Ser. B Stat. Methodol.  62 ,  (2000)  , pp. 335--354

Joerg Polzehl, Vladimir Spokoiny,  Propagation-separation approach for local likelihood estimation,  Probab. Theory Related Fields  135 (3),  (2006)  , pp. 335--362.
}
\author{
Joerg Polzehl, \email{polzehl@wias-berlin.de}
}
\seealso{
\code{\link{aws}}, \code{\link{lpaws}}, \code{\link{aws.irreg}}, \code{\link{aws.gaussian}}
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("aws")
}
\keyword{classes}
