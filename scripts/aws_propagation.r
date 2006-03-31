#
#   propagation condition for aws  (alpha=0.1)   skern = "Triangle"
#
#  univariate Gaussian
set.seed(1)
y <- rnorm(50000)
yhat<-aws(y,qlambda=.98,qtau=1,hmax=1000,testprop=TRUE,lseq=c(1.5),graph=TRUE)    
#  bivariate Gaussian
set.seed(1)
y <- matrix(rnorm(512^2),512,512)  
yhat<-aws(y,qlambda=.97,hmax=10,qtau=1,testprop=TRUE,lseq=c(1.8,1.3,1.2,1.2,1.1,1.1,1.1),graph=TRUE)  
#  3D Gaussian
set.seed(1)
y <- array(rnorm(64^3),c(64,64,64))  
yhat<-aws(y,qlambda=.97,hmax=10,qtau=1,testprop=TRUE,lseq=c(1.9,1.5,1.3,1.3,1.3,1.3,rep(1.1,8)))

#  univariate Poisson
set.seed(1)
y <- rpois(50000,1)
yhat<-aws(y,qlambda=.965,qtau=1,family="Poisson",hmax=1000,testprop=TRUE,lseq=1,u=1)    
#  bivariate Poisson
set.seed(1)
y <- matrix(rpois(512^2,1),512,512)  
yhat<-aws(y,qlambda=.965,hmax=10,qtau=1,family="Poisson",testprop=TRUE,lseq=1,u=1)  
#  3D Poisson
set.seed(1)
y <- array(rpois(64^3,1),c(64,64,64))  
yhat<-aws(y,qlambda=.965,hmax=6,qtau=1,family="Poisson",testprop=TRUE,lseq=1,u=1)

#  univariate Bernoulli
set.seed(1)
y <- rbinom(50000,1,.5)
yhat<-aws(y,qlambda=.97,qtau=1,hmax=1000,family="Bernoulli",testprop=TRUE,lseq=1,u=.5,graph=TRUE)    
#  bivariate Bernoulli
set.seed(1)
y <- matrix(rbinom(512^2,1,.5),512,512)  
yhat<-aws(y,qlambda=.965,hmax=10,qtau=1,family="Bernoulli",testprop=TRUE,lseq=1,u=.5)  
#  3D Bernoulli
set.seed(1)
y <- array(rbinom(64^3,1,.5),c(64,64,64))  
yhat<-aws(y,qlambda=.97,hmax=6,qtau=1,family="Bernoulli",testprop=TRUE,lseq=1,u=.5)

#  univariate Exponential
set.seed(1)
y <- rexp(50000,1)
yhat<-aws(y,qlambda=.965,qtau=1,hmax=1000,family="Exponential",testprop=TRUE,lseq=3.9,u=1,graph=TRUE)    
#  bivariate Exponential
set.seed(1)
y <- matrix(rexp(512^2,1),512,512)  
yhat<-aws(y,qlambda=.965,hmax=10,qtau=1,family="Exponential",testprop=TRUE,lseq=4.8,u=1)  
#  3D Exponential
set.seed(1)
y <- array(rexp(64^3,1),c(64,64,64))  
yhat<-aws(y,qlambda=.975,hmax=6,qtau=1,family="Exponential",testprop=TRUE,lseq=4.7,u=1)

#  univariate Variance
set.seed(1)
y <- rchisq(50000,1)
yhat<-aws(y,qlambda=.965,qtau=1,hmax=1000,family="Variance",testprop=TRUE,lseq=9,u=1,graph=TRUE)    
#  bivariate Variance
set.seed(1)
y <- matrix(rchisq(512^2,1),512,512)  
yhat<-aws(y,qlambda=.965,hmax=10,qtau=1,family="Variance",testprop=TRUE,lseq=11,u=1)  
#  3D Variance
set.seed(1)
y <- array(rchisq(64^3,1),c(64,64,64))  
yhat<-aws(y,qlambda=.975,hmax=6,qtau=1,family="Variance",testprop=TRUE,lseq=11,u=1)

#  univariate Volatility
set.seed(1)
y <- rnorm(50000,0,1)
yhat<-aws(y,qlambda=.965,qtau=1,hmax=1000,family="Volatility",testprop=TRUE,lseq=9,u=1,graph=TRUE)    
#  bivariate Volatility
set.seed(1)
y <- matrix(rnorm(512^2,0,1),512,512)  
yhat<-aws(y,qlambda=.965,hmax=10,qtau=1,family="Volatility",testprop=TRUE,lseq=11,u=1)  
#  3D Volatility
set.seed(1)
y <- array(rnorm(64^3,0,1),c(64,64,64))  
yhat<-aws(y,qlambda=.975,hmax=6,qtau=1,family="Volatility",testprop=TRUE,lseq=11,u=1)


#
#   propagation condition for aws  (alpha=0.1)   skern = "Exponential"
#
set.seed(1)
y <- rnorm(50000)
yhat<-aws(y,qlambda=.992,skern="Exponential",qtau=1,hmax=1000,testprop=TRUE,lseq=c(1.5))    
#  bivariate Gaussian
set.seed(1)
y <- matrix(rnorm(512^2),512,512)  
yhat<-aws(y,qlambda=.992,skern="Exponential",hmax=10,qtau=1,testprop=TRUE,lseq=c(1.8,1.3,1.2,1.2,1.1,1.1,1.1))  
#  3D Gaussian
set.seed(1)
y <- array(rnorm(64^3),c(64,64,64))  
yhat<-aws(y,qlambda=.995,skern="Exponential",hmax=10,qtau=1,testprop=TRUE,lseq=c(1.9,1.5,1.3,1.3,1.3,1.3,rep(1.1,8)))

#  univariate Poisson
set.seed(1)
y <- rpois(50000,1)
yhat<-aws(y,qlambda=.981,skern="Exponential",qtau=1,family="Poisson",hmax=1000,testprop=TRUE,lseq=1,u=1)    
#  bivariate Poisson
set.seed(1)
y <- matrix(rpois(512^2,1),512,512)  
yhat<-aws(y,qlambda=.99,skern="Exponential",hmax=10,qtau=1,family="Poisson",testprop=TRUE,lseq=1,u=1)  
#  3D Poisson
set.seed(1)
y <- array(rpois(64^3,1),c(64,64,64))  
yhat<-aws(y,qlambda=.995,skern="Exponential",hmax=6,qtau=1,family="Poisson",testprop=TRUE,lseq=1,u=1)

#  univariate Bernoulli
set.seed(1)
y <- rbinom(50000,1,.5)
yhat<-aws(y,qlambda=.985,skern="Exponential",qtau=1,hmax=1000,family="Bernoulli",testprop=TRUE,lseq=1,u=.5)    
#  bivariate Bernoulli
set.seed(1)
y <- matrix(rbinom(512^2,1,.5),512,512)  
yhat<-aws(y,qlambda=.99,skern="Exponential",hmax=10,qtau=1,family="Bernoulli",testprop=TRUE,lseq=1,u=.5)  
#  3D Bernoulli
set.seed(1)
y <- array(rbinom(64^3,1,.5),c(64,64,64))  
yhat<-aws(y,qlambda=.995,skern="Exponential",hmax=6,qtau=1,family="Bernoulli",testprop=TRUE,lseq=1,u=.5)

#  univariate Exponential
set.seed(1)
y <- rexp(50000,1)
yhat<-aws(y,qlambda=.985,skern="Exponential",qtau=1,hmax=1000,family="Exponential",testprop=TRUE,lseq=3.9,u=1)    
#  bivariate Exponential
set.seed(1)
y <- matrix(rexp(512^2,1),512,512)  
yhat<-aws(y,qlambda=.985,skern="Exponential",hmax=10,qtau=1,family="Exponential",testprop=TRUE,lseq=4.8,u=1)  
#  3D Exponential
set.seed(1)
y <- array(rexp(64^3,1),c(64,64,64))  
yhat<-aws(y,qlambda=.995,skern="Exponential",hmax=6,qtau=1,family="Exponential",testprop=TRUE,lseq=4.7,u=1)

#  univariate Variance
set.seed(1)
y <- rchisq(50000,1)
yhat<-aws(y,qlambda=.985,skern="Exponential",qtau=1,hmax=1000,family="Variance",testprop=TRUE,lseq=9,u=1)    
#  bivariate Variance
set.seed(1)
y <- matrix(rchisq(512^2,1),512,512)  
yhat<-aws(y,qlambda=.99,skern="Exponential",hmax=10,qtau=1,family="Variance",testprop=TRUE,lseq=11,u=1)  
#  3D Variance
set.seed(1)
y <- array(rchisq(64^3,1),c(64,64,64))  
yhat<-aws(y,qlambda=.995,skern="Exponential",hmax=6,qtau=1,family="Variance",testprop=TRUE,lseq=11,u=1)

#  univariate Volatility
set.seed(1)
y <- rnorm(50000,0,1)
yhat<-aws(y,qlambda=.985,skern="Exponential",qtau=1,hmax=1000,family="Volatility",testprop=TRUE,lseq=9,u=1)    
#  bivariate Volatility
set.seed(1)
y <- matrix(rnorm(512^2,0,1),512,512)  
yhat<-aws(y,qlambda=.99,skern="Exponential",hmax=10,qtau=1,family="Volatility",testprop=TRUE,lseq=11,u=1)  
#  3D Volatility
set.seed(1)
y <- array(rnorm(64^3,0,1),c(64,64,64))  
yhat<-aws(y,qlambda=.995,skern="Exponential",hmax=6,qtau=1,family="Volatility",testprop=TRUE,lseq=11,u=1)

#
#   propagation condition for stagewise aggregation 
#
# based on aggkern="Uniform"   limit number of estimates that are fixed to 2 out of 10000
#
set.seed(1)
y <- rnorm(50000)
yhat<-aws(y,qlambda=1,qtau=.4,heta=1.25,hmax=1000,testprop=TRUE,graph=TRUE)    
#  bivariate Gaussian
set.seed(1)
y <- matrix(rnorm(512^2),512,512)  
yhat<-aws(y,qlambda=1,qtau=.35,heta=1.25,hmax=10,testprop=TRUE,graph=TRUE)    
#  3D Gaussian
set.seed(1)
y <- array(rnorm(64^3),c(64,64,64))  
yhat<-aws(y,qlambda=1,qtau=.35,heta=1.25,hmax=6,testprop=TRUE)    

#  univariate Poisson
set.seed(1)
y <- rpois(50000,1)
yhat<-aws(y,qlambda=1,qtau=.45,family="Poisson",hmax=1000,testprop=TRUE,graph=TRUE,u=1)    
#  bivariate Poisson
set.seed(1)
y <- matrix(rpois(512^2,1),512,512)  
yhat<-aws(y,qlambda=1,qtau=.35,family="Poisson",hmax=10,testprop=TRUE,graph=TRUE,u=1)    
#  3D Poisson
set.seed(1)
y <- array(rpois(64^3,1),c(64,64,64))  
yhat<-aws(y,qlambda=1,qtau=.35,family="Poisson",hmax=6,testprop=TRUE,u=1)    

#  univariate Bernoulli
set.seed(1)
y <- rbinom(50000,1,.5)
yhat<-aws(y,qlambda=1,qtau=.45,heta=1.25,family="Bernoulli",hmax=1000,testprop=TRUE,graph=TRUE,u=.5)    
y <- rbinom(50000,1,.05)
yhat<-aws(y,qlambda=1,qtau=.45,heta=1.25,family="Bernoulli",hmax=1000,testprop=TRUE,graph=TRUE,u=.5)    
#  bivariate Bernoulli
set.seed(1)
y <- matrix(rbinom(512^2,1,.5),512,512)  
yhat<-aws(y,qlambda=1,qtau=.55,heta=1.25,family="Bernoulli",hmax=10,testprop=TRUE,graph=TRUE,u=.5)    
#  3D Bernoulli
set.seed(1)
y <- array(rbinom(64^3,1,.005),c(64,64,64))  
yhat<-aws(y,qlambda=1,qtau=.55,heta=1.25,family="Bernoulli",hmax=6,testprop=TRUE,u=.05)    

#  univariate Exponential
set.seed(1)
y <- rexp(50000,1)
yhat<-aws(y,qlambda=1,qtau=.8,hmax=1000,family="Exponential",testprop=TRUE,graph=TRUE,u=1)    
#  bivariate Exponential
set.seed(1)
y <- matrix(rexp(512^2,1),512,512)  
yhat<-aws(y,qlambda=1,qtau=.8,hmax=1000,family="Exponential",testprop=TRUE,graph=TRUE,u=1)    
#  3D Exponential
set.seed(1)
y <- array(rexp(64^3,1),c(64,64,64))  
yhat<-aws(y,qlambda=.975,skern="Exponential",hmax=6,qtau=1,family="Exponential",testprop=TRUE,lseq=4.7,u=1)

#  univariate Variance
set.seed(1)
y <- rchisq(50000,1)
yhat<-aws(y,qlambda=.965,skern="Exponential",qtau=1,hmax=1000,family="Variance",testprop=TRUE,lseq=9,u=1)    
#  bivariate Variance
set.seed(1)
y <- matrix(rchisq(512^2,1),512,512)  
yhat<-aws(y,qlambda=.965,skern="Exponential",hmax=10,qtau=1,family="Variance",testprop=TRUE,lseq=11,u=1)  
#  3D Variance
set.seed(1)
y <- array(rchisq(64^3,1),c(64,64,64))  
yhat<-aws(y,qlambda=.975,skern="Exponential",hmax=6,qtau=1,family="Variance",testprop=TRUE,lseq=11,u=1)

#  univariate Volatility
set.seed(1)
y <- rnorm(50000,0,1)
yhat<-aws(y,qlambda=.965,skern="Exponential",qtau=1,hmax=1000,family="Volatility",testprop=TRUE,lseq=9,u=1)    
#  bivariate Volatility
set.seed(1)
y <- matrix(rnorm(512^2,0,1),512,512)  
yhat<-aws(y,qlambda=.965,skern="Exponential",hmax=10,qtau=1,family="Volatility",testprop=TRUE,lseq=11,u=1)  
#  3D Volatility
set.seed(1)
y <- array(rnorm(64^3,0,1),c(64,64,64))  
yhat<-aws(y,qlambda=.975,skern="Exponential",hmax=6,qtau=1,family="Volatility",testprop=TRUE,lseq=11,u=1)
