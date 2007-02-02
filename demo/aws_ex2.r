require(aws)
f1 <- function(x){
     xj <- c(20,30,55,65,80)/100
     hj <- c(50,-30,50,-40,30)
     Kern <- function(x) (1-sign(x))/2
     (16.5+apply(Kern(outer(xj,x,"-"))*hj,2,sum))/50
     }
x <- seq(0,1,length=2048)
ifamily <- as.integer(readline("Select distribution: 1 for Bernoulli,  2 for Poisson, 3 for Exponential, 4 for Volatility model"))
if(!(ifamily %in% (1:4))) ifamily <- 2 
factor <- as.numeric(switch(ifamily,0.5,
     readline("Enter mean intensitity:"),
     readline("Enter mean intensitity:"),
     readline("Enter mean standard deviation:")))
if(is.na(factor)) factor <- switch(ifamily,1,10,1,1)
if(factor <= 0) factor <- switch(ifamily,1,10,1,1)
u <- factor*f1(x)
family <- switch(ifamily,"Bernoulli","Poisson","Exponential","Volatility")
plot(x,u,type="l",col=2,lwd=2)
title(paste(switch(ifamily,"Probability","Intensity","Intensity","Standard deviation"),"in",family,"model"))
y <- switch(ifamily,rbinom(x,1,u),rpois(x,u),rexp(x,1/u),rnorm(x,0,u))
plot(x,u,type="l",col=2,ylim=range(u,y),lwd=2)
title(paste(" Data and", switch(ifamily,"Probability","Intensity","Intensity","Standard deviation"),"in",family,"model"))
points(x,y)
if(ifamily==4) lines(x,-u,col=2, lwd=2)
hmax <- readline("Maximal bandwidth:\n Press 'Enter' for hmax=250, otherwise provide value of hmax:")
if(is.na(as.numeric(hmax))) hmax <- 250 else hmax <- as.numeric(hmax)
if(hmax <= 1) hmax <- 250
spmin <- readline("Kernel form:\n Press 'Enter' for spmin=0, otherwise provide value of spmin:")
if(is.na(as.numeric(spmin))) spmin <- 0 else spmin <- pmax(0,pmin(1,as.numeric(spmin)))
cat("Run aws \n")
yhat <- aws(y,hmax=hmax,family=family,graph=TRUE,spmin=spmin,qtau=1)
readline("Press ENTER to show results")
if(ifamily==4) {
y <- abs(y)
yhat$theta <- sqrt( yhat$theta)
}
plot(x,u,type="l",col=2,ylim=range(u,y),lwd=2)
points(x,y)
lines(x,yhat$theta,col=3,lwd=2)
title("Data, fitted and true values")
rm(f1,x,u,y,hmax,yhat,factor,ifamily,family)
