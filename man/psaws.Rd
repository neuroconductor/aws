\name{psaws}
\alias{psaws}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
psaws(y, x = NULL, p = 0, hmax = NULL, model = "Regression", family = "Gaussian", homoskedastic = TRUE, sigma2 = NULL, shape = NULL, lkern = "Triangle", NN = FALSE, u = NULL, graph = FALSE, demo = FALSE, aws = TRUE, control = TRUE, wghts = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{p}{ ~~Describe \code{p} here~~ }
  \item{hmax}{ ~~Describe \code{hmax} here~~ }
  \item{model}{ ~~Describe \code{model} here~~ }
  \item{family}{ ~~Describe \code{family} here~~ }
  \item{homoskedastic}{ ~~Describe \code{homoskedastic} here~~ }
  \item{sigma2}{ ~~Describe \code{sigma2} here~~ }
  \item{shape}{ ~~Describe \code{shape} here~~ }
  \item{lkern}{ ~~Describe \code{lkern} here~~ }
  \item{NN}{ ~~Describe \code{NN} here~~ }
  \item{u}{ ~~Describe \code{u} here~~ }
  \item{graph}{ ~~Describe \code{graph} here~~ }
  \item{demo}{ ~~Describe \code{demo} here~~ }
  \item{aws}{ ~~Describe \code{aws} here~~ }
  \item{control}{ ~~Describe \code{control} here~~ }
  \item{wghts}{ ~~Describe \code{wghts} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, x = NULL, p = 0, hmax = NULL, model = "Regression", 
    family = "Gaussian", homoskedastic = TRUE, sigma2 = NULL, 
    shape = NULL, lkern = "Triangle", NN = FALSE, u = NULL, graph = FALSE, 
    demo = FALSE, aws = TRUE, control = TRUE, wghts = NULL) 
{
    ergs <- "not implemented"
    if (aws) 
        qlambda <- NULL
    else qlambda <- 1
    if (control) 
        heta <- NULL
    else heta <- 1e+10
    if (model == "Regression") {
        if (p == 0) {
            if (homoskedastic) {
                ergs <- caws(y, x, lkern = lkern, family = family, 
                  sigma2 = sigma2, qlambda = qlambda, heta = heta, 
                  shape = shape, hmax = hmax, NN = NN, u = u, 
                  graph = graph, demo = demo, wghts = wghts, 
                  demo = demo)
            }
            else {
                ergs <- chaws(y, x, lkern = lkern, family = family, 
                  sigma2 = sigma2, qlambda = qlambda, heta = heta, 
                  shape = shape, hmax = hmax, NN = NN, u = u, 
                  graph = graph, demo = demo, wghts = wghts, 
                  demo = demo)
            }
        }
        else {
            if (homoskedastic) {
                ergs <- cpaws(y, x, p = p, lkern = lkern, sigma2 = sigma2, 
                  qlambda = qlambda, heta = heta, hmax = hmax, 
                  NN = NN, u = u, graph = graph, demo = demo, 
                  wghts = wghts, demo = demo)
            }
            else {
                ergs <- cphaws(y, x, p = p, lkern = lkern, sigma2 = sigma2, 
                  qlambda = qlambda, heta = heta, hmax = hmax, 
                  NN = NN, u = u, graph = graph, demo = demo, 
                  wghts = wghts, demo = demo)
            }
        }
    }
    ergs
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
