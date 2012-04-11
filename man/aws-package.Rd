\name{aws-package}
\alias{aws-package}
\docType{package}
\title{Adaptive Weights Smoothing}
\description{
The package contains R-functions implementing the Propagation-Separation Approach to adaptive smoothing as described
in J. Polzehl and V. Spokoiny (2006) Propagation-Separation Approach for Local Likelihood Estimation,
Prob. Theory and Rel. Fields 135(3):335-362. and
J. Polzehl and V. Spokoiny (2004) Spatially adaptive regression estimation: Propagation-separation approach,
WIAS-Preprint 998. Additionally it contains an implementation of selected LPA-ICI pointwise adaptive
smoothing algorithms from the book V. Katkovnik, K. Egiazarian and J. Astola (2006). Local Approximation
Techniques in Signal and Image Processing, SPIE Press Monograph Vol. PM 157. 
 }
\details{
\tabular{ll}{
Package: \tab aws\cr
Version: \tab 1.6\cr
Date: \tab 2009-04-07\cr
License: \tab GPL (>=2)\cr
Copyright: \tab 2008 Weierstrass Institute for\cr
           \tab Applied Analysis and Stochastics.\cr
URL: \tab http://www.wias-berlin.de/project-areas/stat/\cr
}

Index:
\preformatted{
aws                     AWS for local constant models on a grid
aws.gaussian            Adaptive weights smoothing for Gaussian data
                        with variance depending on the mean.
aws.irreg               local constant AWS for irregular (1D/2D) design
aws.segment             Segmentation by adaptive weights for Gaussian
                        models.
awsdata                 Extract information from an object of class aws
binning                 Binning in 1D, 2D or 3D
lpaws                   Local polynomial smoothing by AWS
kernsm                  1D, 2D, 3D nonparametric kernel smoothing via fft
ICIsmooth               pointwise adaptive kernel smoothing  
ICIcombined             pointwise adaptive kernel smoothing with fusing
}
}
\author{
Joerg Polzehl <polzehl@wias-berlin.de>

Maintainer: Joerg Polzehl <polzehl@wias-berlin.de>
}
\references{ J. Polzehl and V. Spokoiny (2006) 
             Propagation-Separation Approach for Local Likelihood Estimation,
             \emph{Prob. Theory and Rel. Fields} \bold{135(3)}, 335-362.

             J. Polzehl and V. Spokoiny (2004)
             Spatially adaptive regression estimation: Propagation-separation approach,
             WIAS-Preprint 998.

             V. Katkovnik, K. Egiazarian and J. Astola (2006) 
             \emph{Local Approximation Techniques in Signal and Image Processing},
              SPIE Press Monograph Vol. PM 157
}
\keyword{ package }
