setClass("aws",
         representation(.Data = "list",
                        y = "numeric",
                        dy = "numeric",
                        x = "numeric",
                        ni = "integer",
                        mask = "logical",
                        theta = "numeric",
                        mae = "numeric",
                        var = "numeric",
                        xmin = "numeric",
                        xmax = "numeric",
                        wghts = "numeric",
                        degree = "integer",
                        hmax  = "numeric", # 
                        sigma2  = "numeric", # error variance
                        scorr = "numeric", # spatial correlation
                        family = "character",
                        shape = "numeric",
                        lkern  = "integer",
                        lambda = "numeric",
                        ladjust = "numeric",
                        aws = "logical",
                        memory = "logical",
                        homogen = "logical",
                        earlystop = "logical",
                        varmodel = "character",
                        vcoef = "numeric",
                        call = "function")
         )

awsobj <- function(y,theta,var,hmax,sigma2,lkern,lambda,ladjust,aws,memory,
call,homogen=FALSE,earlystop=FALSE,family="Gaussian",degree=0,
x=numeric(0),ni=as.integer(1),xmin=numeric(0),xmax=numeric(0),
wghts=numeric(0),scorr=0,mae=numeric(0),shape=numeric(0),
varmodel="Constant",vcoef=NULL,mask=logical(0)){
              dy <- if(is.null(dim(y))) length(y) else dim(y)
              invisible(new("aws",
                        y = as.numeric(y),
                        dy = dy,
                        x = as.numeric(x),
                        ni = as.integer(ni),
                        mask = as.logical(mask),
                        theta = as.numeric(theta),
                        mae = as.numeric(mae),
                        var = as.numeric(var),
                        xmin = as.numeric(xmin),
                        xmax = as.numeric(xmax),
                        wghts = as.numeric(wghts),
                        degree = as.integer(degree),
                        hmax  = as.numeric(hmax), # 
                        sigma2  = as.numeric(sigma2), # error variance
                        scorr = as.numeric(scorr), # spatial correlation
                        family = family,
                        shape = as.numeric(shape),
                        lkern  = as.integer(lkern),
                        lambda = as.numeric(lambda),
                        ladjust = as.numeric(ladjust),
                        aws = aws,
                        memory = memory,
                        homogen = homogen,
                        earlystop = earlystop,
                        varmodel = varmodel,
                        vcoef = if(is.null(vcoef)) mean(sigma2) else vcoef,
                        call = call)
            )
}

awsdata <- function(awsobj,what){
what <- tolower(what)
switch(what,"y"=array(awsobj@y,awsobj@dy),
                 "data"=array(awsobj@y,awsobj@dy),
                 "theta"=array(awsobj@theta,c(awsobj@dy,awsobj@degree+1)),
                 "est"=array(awsobj@theta,c(awsobj@dy)),
                 "var"=array(awsobj@var,awsobj@dy),
                 "sd"=array(sqrt(awsobj@var),awsobj@dy),
                 "sigma2"=array(awsobj@sigma2,awsobj@dy),
                 "mae"=awsobj@mae,
                 "ni"=array(awsobj@ni,awsobj@dy),
                 "mask"=array(awsobj@mask,awsobj@dy))
}