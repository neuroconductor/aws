Varcor<-function(lkern,h,d=1){
#
#   Calculates a correction for the variance estimate obtained by (IQRdiff(y)/1.908)^2
#
#   in case of colored noise that was produced by smoothing with lkern and bandwidth h
#
if(lkern=="Gaussian") h<-h/2.3548
ih<-switch(lkern,Gaussian=trunc(4*h)+1,trunc(h)+1)
dx<-2*ih+1
x<- ((-ih):ih)/h
if(d==2) x<-sqrt(outer(x^2,x^2,"+"))
if(d==3) x<-sqrt(outer(x^2,outer(x^2,x^2,"+"),"+"))
penl<-switch(lkern,Triangle=pmax(0,1-x^2),
                   Uniform=as.numeric(abs(x)<=1),
                   Quadratic=pmax(0,1-x^2)^2,
                   Cubic=pmax(0,1-x^2)^3,
                   Gaussian=dnorm(x))
.5/sum(diff(penl)^2)*sum(abs(diff(penl)))^2
}
Varcor.gauss<-function(h){
#
#   Calculates a correction for the variance estimate obtained by (IQRdiff(y)/1.908)^2
#
#   in case of colored noise that was produced by smoothing with lkern and bandwidth h
#
h<-h/2.3548
ih<-trunc(4*h)+1
dx<-2*ih+1
d<-length(h)
penl <- dnorm(((-ih[1]):ih[1])/h[1])
if(d==2) penl <- outer(penl,dnorm(((-ih[2]):ih[2])/h[2]),"*")
if(d==3) penl <- outer(penl,outer(dnorm(((-ih[2]):ih[2])/h[2]),dnorm(((-ih[3]):ih[3])/h[3]),"*"),"*")
.5/sum(diff(penl)^2)*sum(abs(diff(penl)))^2
}

SpatialCorr<-function(lkern,h,d=1){
#
#   Calculates the correlation of 
#
#   colored noise that was produced by smoothing with lkern and bandwidth h
#
#  !!! the result is not monotone in h for   lkern="Triangle" (all d) and lkern="Uniform" (d>1)
#
if(lkern=="Gaussian") h<-h/2.3548
ih<-switch(lkern,Gaussian=trunc(4*h+1),trunc(h)+1)
dx<-2*ih+1
x<- ((-ih):ih)/h
if(d==2) x<-sqrt(outer(x^2,x^2,"+"))
if(d==3) x<-sqrt(outer(x^2,outer(x^2,x^2,"+"),"+"))
penl<-as.vector(switch(lkern,Triangle=pmax(0,1-x^2),
                   Uniform=as.numeric(abs(x)<=1),
                   Quadratic=pmax(0,1-x^2)^2,
                   Cubic=pmax(0,1-x^2)^3,
                   Gaussian=dnorm(x)))
dim(penl)<-rep(dx,d)
if(d==1) z<-sum(penl[-1]*penl[-dx])/sum(penl^2)
if(d==2) z<-sum(penl[-1,]*penl[-dx,])/sum(penl^2)
if(d==3) z<-sum(penl[-1,,]*penl[-dx,,])/sum(penl^2)
z
}

SpatialCorr.gauss<-function(h){
#
#   Calculates the correlation of 
#
#   colored noise that was produced by smoothing with "gaussian" kernel and bandwidth h
#
#   Result does not depend on d for "Gaussian" kernel !!
#
h<-h/2.3548
ih<-trunc(4*h+1)
dx<-2*ih+1
penl<-dnorm(((-ih):ih)/h)
sum(penl[-1]*penl[-dx])/sum(penl^2)
}

Spatialvar<-function(lkern,lkern0,h,h0,d){
#
#   Calculates the factor of variance reduction obtained at bandwidth h in 
#
#   case of colored noise that was produced by smoothing with lkern0 and bandwidth h0
#
#   Spatialvariance(lkern,h,h0,d)/Spatialvariance(lkern,h,1e-5,d) gives the 
#   a factor for lambda to be used with bandwidth h 
#
if(lkern=="Gaussian") h<-h/2.3548
ih<-switch(lkern,Gaussian=trunc(4*h),trunc(h))
ih<-max(1,ih)
dx<-2*ih+1
x<- ((-ih):ih)/h
if(d==2) x<-sqrt(outer(x^2,x^2,"+"))
if(d==3) x<-sqrt(outer(x^2,outer(x^2,x^2,"+"),"+"))
penl<-switch(lkern,Triangle=pmax(0,1-x^2),
                   Uniform=as.numeric(abs(x)<=1),
                   Quadratic=pmax(0,1-x^2)^2,
                   Cubic=pmax(0,1-x^2)^3,
                   Gaussian=dnorm(x))
dim(penl)<-rep(dx,d)
if(lkern0=="Gaussian") h0<-h0/2.3548
ih<-switch(lkern0,Gaussian=trunc(4*h0),trunc(h0))
ih<-max(1,ih)
dx0<-2*ih+1
x<- ((-ih):ih)/h0
if(d==2) x<-sqrt(outer(x^2,x^2,"+"))
if(d==3) x<-sqrt(outer(x^2,outer(x^2,x^2,"+"),"+"))
penl0<-switch(lkern0,Triangle=pmax(0,1-x^2),
                   Uniform=as.numeric(abs(x)<=1),
                   Quadratic=pmax(0,1-x^2)^2,
                   Cubic=pmax(0,1-x^2)^3,
                   Gaussian=dnorm(x))
dim(penl0)<-rep(dx0,d)
penl0<-penl0/sum(penl0)
dz<-dx+dx0-1
z<-array(0,rep(dz,d))
if(d==1){
for(i1 in 1:dx0) {
ind1<-c(0:(i-1),(dz-dx0+i):dz+1)
ind1<-ind1[ind1<=dz][-1]
z[-ind1]<-z[-ind1]+penl*penl0[i]
}
} else if(d==2){
for(i1 in 1:dx0) for(i2 in 1:dx0){
ind1<-c(0:(i1-1),(dz-dx0+i1):dz+1)
ind1<-ind1[ind1<=dz][-1]
ind2<-c(0:(i2-1),(dz-dx0+i2):dz+1)
ind2<-ind2[ind2<=dz][-1]
z[-ind1,-ind2]<-z[-ind1,-ind2]+penl*penl0[i1,i2]
}
} else if(d==3){
for(i1 in 1:dx0) for(i2 in 1:dx0) for(i3 in 1:dx0){
ind1<-c(0:(i1-1),(dz-dx0+i1):dz+1)
ind1<-ind1[ind1<=dz][-1]
ind2<-c(0:(i2-1),(dz-dx0+i2):dz+1)
ind2<-ind2[ind2<=dz][-1]
ind3<-c(0:(i3-1),(dz-dx0+i3):dz+1)
ind3<-ind3[ind3<=dz][-1]
z[-ind1,-ind2,-ind3]<-z[-ind1,-ind2,-ind3]+penl*penl0[i1,i2,i3]
}
}
sum(z^2)/sum(z)^2
}

Spatialvar.gauss<-function(h,h0,d){
#
#   Calculates the factor of variance reduction obtained for Gaussian Kernel and bandwidth h in 
#
#   case of colored noise that was produced by smoothing with Gaussian kernel and bandwidth h0
#
#   Spatialvariance(lkern,h,h0,d)/Spatialvariance(lkern,h,1e-5,d) gives the 
#   a factor for lambda to be used with bandwidth h 
#
h<-h/2.3548
ih<-trunc(4*h)
ih<-max(1,ih)
dx<-2*ih+1
x<- ((-ih):ih)/h
penl<-dnorm(x)
if(d==2) penl<-outer(penl,penl,"*")
if(d==3) penl<-outer(penl,outer(penl,penl,"*"),"*")
dim(penl)<-rep(dx,d)
h0<-h0/2.3548
if(length(h0)==1) h0<-rep(h0,d)
ih<-trunc(4*h0)
ih<-pmax(1,ih)
dx0<-2*ih+1
x<- ((-ih[1]):ih[1])/h0[1]
penl0<-dnorm(((-ih[1]):ih[1])/h0[1])
if(d==2) penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),dnorm(((-ih[2]):ih[2])/h0[2]),"*")
if(d==3) penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),outer(dnorm(((-ih[2]):ih[2])/h0[2]),dnorm(((-ih[3]):ih[3])/h0[3]),"*"),"*")
dim(penl0)<-dx0
penl0<-penl0/sum(penl0)
dz<-dx+dx0-1
z<-array(0,dz)
if(d==1){
for(i1 in 1:dx0) {
ind1<-c(0:(i-1),(dz-dx0+i):dz+1)
ind1<-ind1[ind1<=dz][-1]
z[-ind1]<-z[-ind1]+penl*penl0[i]
}
} else if(d==2){
for(i1 in 1:dx0[1]) for(i2 in 1:dx0[2]){
ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
ind1<-ind1[ind1<=dz[1]][-1]
ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
ind2<-ind2[ind2<=dz[2]][-1]
z[-ind1,-ind2]<-z[-ind1,-ind2]+penl*penl0[i1,i2]
}
} else if(d==3){
for(i1 in 1:dx0[1]) for(i2 in 1:dx0[1]) for(i3 in 1:dx0[1]){
ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
ind1<-ind1[ind1<=dz[1]][-1]
ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
ind2<-ind2[ind2<=dz[2]][-1]
ind3<-c(0:(i3-1),(dz[3]-dx0[3]+i3):dz[3]+1)
ind3<-ind3[ind3<=dz[3]][-1]
z[-ind1,-ind2,-ind3]<-z[-ind1,-ind2,-ind3]+penl*penl0[i1,i2,i3]
}
}
sum(z^2)/sum(z)^2
}

geth.gauss<-function(corr,step=1.01){
#   get the   bandwidth for lkern corresponding to a given correlation
h<-.1
z<-0
# 
#  keep it simple result does not depend on d
#
while(z<corr){
h<-h*step
z<-SpatialCorr.gauss(h)
}
h
}