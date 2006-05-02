require(aws)
fw0 <- function(){
         xy <- rbind(rep(0:255,256),rep(0:255,rep(256,256)))
         indw <- c(1:12,29:48,73:100,133:168,209:256)
         w0 <- matrix(rep(FALSE,256*256),ncol=256)
         w0[indw,] <- TRUE
         w0[,indw] <- !w0[,indw]
         w0 <- w0-.5
	 z <- (xy[1,]-129)^2+(xy[2,]-129)^2
         w0[z<=10000&z>=4900] <- 0
         w0[abs(xy[1,]-xy[2,])<=20&z<4900] <- 0
	 z <- (xy[1,]-225)^2+2*(xy[2,]-30)^2-(xy[1,]-225)*(xy[2,]-30)
         w0[z<=625] <- 0
         w0[z<=625&xy[2,]>27&xy[2,]<31] <- -.5
         w0[z<=625&xy[1,]>223&xy[1,]<227] <- .5
	 z <- ((xy[2,]-225)^2+2*(xy[1,]-30)^2)+(xy[2,]-225)*(xy[1,]-30)
         w0[z<=625] <- 0
         w0[z<=625&xy[1,]>27&xy[1,]<31] <- -.5
         w0[z<=625&xy[2,]>223&xy[2,]<227] <- .5
         w0[((xy[2,]-225)^2+(xy[1,]-225)^2)+1*(xy[2,]-225)*(xy[1,]-225)<=400] <- 0
         w0[((xy[2,]-30)^2+(xy[1,]-30)^2)<=256] <- 0
	 w0
	 }
w0 <- fw0()
image(w0,col=gray((0:255)/255))
title("Original image")
sigma <- readline("Standard deviation of noise:\n Press 'Enter' for sigma=0.5, otherwise provide value of sigma:")
if(is.na(as.numeric(sigma))) sigma <- 0.5 else sigma <- as.numeric(sigma)
if(sigma <= 0) sigma <- 0.1
y <- w0+rnorm(w0,0,sigma)
image(y,col=gray((0:255)/255))
title("Noisy image")
hmax <- readline("Maximal bandwidth:\n Press 'Enter' for hmax=10, otherwise provide value of hmax:")
if(is.na(as.numeric(hmax))) hmax <- 10 else hmax <- as.numeric(hmax)
if(hmax <= 1) hmax <- 10
risk <- readline("Report risks (N/Y):")
if(risk %in% c("y","Y")) u <-w0 else u <- NULL
cat("Run aws \n")
yhat <- aws(y,hmax=hmax,graph=TRUE,u=u,qtau=1)
rm(fw0,w0,y,sigma,hmax,yhat,u,risk)
