\name{risk-methods}
\docType{methods}
\alias{risk-methods}
\alias{risk,ANY-method}
\alias{risk,array-method}
\alias{risk,aws-method}
\alias{risk,awssegment-method}
\alias{risk,ICIsmooth-method}
\alias{risk,kernsm-method}
\alias{risk,numeric-method}
\title{Compute risks characterizing the quality of smoothing results}
\description{
Methods function \code{risk} in package \pkg{aws}. For an given array \code{u} 
the following statistics are computed :
Root Mean Squared Error \code{RMSE <- sqrt(mean((y-u)^2))},
Signal to Noise Ratio \code{SNR <- 10*log(mean(u^2)/MSE,10)},
Peak Signal to Noise Ratio \code{PSNR <- 10*log(max(u^2)/MSE,10)},
Mean Absolute Error \code{MAE <- mean(abs(y-u))},
Maximal Absolute Error \code{MaxAE <- max(abs(y-u))},
Universal Image Quality Index (UIQI) (Wang and Bovik (2002)).
}
\usage{
  \S4method{risk}{array}(y, u=0)
  \S4method{risk}{aws}(y, u=0)
  \S4method{risk}{awssegment}(y, u=0)
  \S4method{risk}{ICIsmooth}(y, u=0)
  \S4method{risk}{kernsm}(y, u=0)
  \S4method{risk}{numeric}(y, u=0)
}
\arguments{
\item{y}{ object }
\item{u}{ array of dimension \code{dim(y)} or \code{dim(extract(y,what="yhat")$y)}
      or scalar value used in comparisons. }
}
\section{Methods}{
\describe{

\item{\code{signature(y = "ANY")}}{
The method extract and/or compute specified statistics from object of class 
}

\item{\code{signature(y = "array")}}{
Returns a list with components \code{RMSE, SNR, PSNR, MAE, MaxAE, UIQI}
}

\item{\code{signature(y = "aws")}}{
Returns a list with components \code{RMSE, SNR, PSNR, MAE, MaxAE, UIQI}
}

\item{\code{signature(y = "awssegment")}}{
Returns a list with components \code{RMSE, SNR, PSNR, MAE, MaxAE, UIQI}
}

\item{\code{signature(y = "ICIsmooth")}}{
Returns a list with components \code{RMSE, SNR, PSNR, MAE, MaxAE, UIQI}
}

\item{\code{signature(y = "kernsm")}}{
Returns a list with components \code{RMSE, SNR, PSNR, MAE, MaxAE, UIQI}
}

\item{\code{signature(y = "numeric")}}{
Returns a list with components \code{RMSE, SNR, PSNR, MAE, MaxAE, UIQI}
}
}}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}
\references{ 
V. Katkovnik, K. Egiazarian and J. Astola, 
\emph{Local Approximation Techniques in Signal And Image Processing},
 SPIE Society of Photo-Optical Instrumentation Engin., 2006, PM157

Z. Wang  and A. C. Bovik, \emph{A universal image quality index}, 
IEEE Signal Processing Letters,
vol. 9, N3, pp. 81-84, 2002.

}
\keyword{methods}
\keyword{smooth}
