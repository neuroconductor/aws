\name{auxiliary}
\alias{gethani}
\alias{getvofh}
\alias{sofmchi}
\alias{residualVariance}
\alias{residualSpatialCorr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Auxiliary functions (for internal use)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gethani(x, y, lkern, value, wght, eps = 0.01)
getvofh(bw, lkern, wght)
sofmchi(L, to = 50, delta = 0.01)
residualVariance(residuals, mask, resscale = 1, compact = FALSE)
residualSpatialCorr(residuals, mask, lags = c(5, 5, 3), compact = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{lkern}{
%%     ~~Describe \code{lkern} here~~
}
  \item{value}{
%%     ~~Describe \code{value} here~~
}
  \item{wght}{
%%     ~~Describe \code{wght} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
\item{bw}{
%%     ~~Describe \code{eps} here~~
}
\item{L}{
%%     ~~Describe \code{eps} here~~
}
\item{to}{
%%     ~~Describe \code{eps} here~~
}
\item{delta}{
%%     ~~Describe \code{eps} here~~
}
\item{residuals}{
%%     ~~Describe \code{eps} here~~
}
\item{mask}{
%%     ~~Describe \code{eps} here~~
}
\item{resscale}{
%%     ~~Describe \code{eps} here~~
}
\item{compact}{
%%     ~~Describe \code{eps} here~~
}
\item{lags}{
%%     ~~Describe \code{eps} here~~
}






}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, lkern, value, wght, eps = 0.01)
{
    .Fortran(C_gethani, as.double(x), as.double(y), as.integer(lkern),
        as.double(value), as.double(wght), as.double(eps), bw = double(1))$bw
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
