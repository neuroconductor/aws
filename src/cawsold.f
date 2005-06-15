CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        one iteration of univariate local constant aws on a grid
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawsuni(y,fix,n,hakt,lambda,theta,bi,bi0,ai,
     1                   model,kern,spmax)
C   
C   y        observed values of regression function
C   n        number of observations
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C
      implicit logical (a-z)
      external kldist,lkern
      real*8 kldist,lkern
      integer n,model,kern
      logical aws,fix(n)
      real*8 y(n),hakt,theta(n),bi(n),bi0(n),ai(n),lambda,spmax
      integer ih,i,j,ja,je
      real*8 thetai,bii,sij,swj,swjy,z,wj,swj0,bii0
      ih=hakt
      aws=lambda.lt.1d40
      DO i=1,n
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         thetai=theta(i)
         bii=bi(i)/lambda
C   scaling of sij outside the loop
         bii0=bi0(i)
         ja=max0(1,i-ih)
         je=min0(n,i+ih)
         swj=0.d0
         swj0=0.d0
         swjy=0.d0
         DO j=ja,je
C  first stochastic term
            z=(i-j)/hakt
            wj=lkern(kern,z*z)
            swj0=swj0+wj
            IF (aws) THEN
               sij=bii*kldist(model,thetai,theta(j),bii0)
               IF (sij.gt.spmax) CYCLE
               wj=wj*dexp(-sij)
            ENDIF
            swj=swj+wj
            swjy = swjy+wj*y(j)
         END DO
         ai(i)=swjy
         bi(i)=swj
         bi0(i)=swj0
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        one iteration of univariate local constant aws on a grid 
C        heteroskedastic case
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chawsuni(y,fix,si2,n,hakt,lambda,theta,bi,bi0,ai,
     1                    kern,spmax)
C   
C   y        observed values of regression function
C   n        number of observations
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel                  
C   spmax    specifies the truncation point of the stochastic kernel
C
      implicit logical (a-z)
      external lkern
      real*8 lkern
      integer n,kern
      logical aws,fix(n)
      real*8 y(n),hakt,theta(n),bi(n),bi0(n),ai(n),lambda,
     1       si2(n),spmax
      integer ih,i,j,ja,je
      real*8 thetai,bii,sij,swj,swjy,z,wj,swj0
      ih=hakt
      aws=lambda.lt.1d40
      DO i=1,n
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         thetai=theta(i)
         bii=bi(i)/lambda
C   scaling of sij outside the loop
         ja=max0(1,i-ih)
         je=min0(n,i+ih)                  
         swj=0.d0
         swj0=0.d0
         swjy=0.d0
         DO j=ja,je
C  first stochastic term
            z=(i-j)/hakt
            wj=lkern(kern,z*z)*si2(j)
            swj0=swj0+wj
            IF (aws) THEN
               z=thetai-theta(j)
               sij=bii*z*z
               IF (sij.gt.spmax) CYCLE
               wj=wj*dexp(-sij)
            END IF
            swj=swj+wj
            swjy=swjy+wj*y(j)
         END DO
         ai(i)=swjy
         bi(i)=swj
         bi0(i)=swj0
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant bivariate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawsbi(y,fix,n1,n2,hakt,lambda,theta,bi,bi0,ai,
     1                  model,kern,spmax,wght)
C   
C   y        observed values of regression function
C   n1,n2    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second dimension (larger values shrink)
C   
      implicit logical (a-z)                  
      external kldist,lkern
      real*8 kldist,lkern
      integer n1,n2,model,kern
      logical aws,fix(n1,n2)
      real*8 y(n1,n2),hakt,theta(n1,n2),bi(n1,n2),ai(n1,n2),
     1       lambda,spmax,wght,bi0(n1,n2)
      integer ih1,ih2,i1,i2,j1,j2,ja1,je1,ja2,je2
      real*8 thetai,bii,sij,swj,swjy,z1,z2,wj,hakt2,swj0,bii0
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO  i1=1,n1
         DO  i2=1,n2
            IF (fix(i1,i2)) CYCLE
C    nothing to do, final estimate is already fixed by control 
            thetai=theta(i1,i2)
            bii=bi(i1,i2)/lambda
C   scaling of sij outside the loop
            bii0=bi0(i1,i2)
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            swj=0.d0
            swj0=0.d0
            swjy=0.d0
            DO j1=ja1,je1
               z1=(i1-j1)
               z1=z1*z1
               ih2=dsqrt(hakt2-z1)/wght
               ja2=max0(1,i2-ih2)
               je2=min0(n2,i2+ih2)
               DO j2=ja2,je2
C  first stochastic term
                  z2=(i2-j2)*wght
                  wj=lkern(kern,(z1+z2*z2)/hakt2)
                  swj0=swj0+wj
                  IF (aws) THEN
                  sij=bii*kldist(model,thetai,theta(j1,j2),bii0)
                  IF (sij.gt.spmax) CYCLE
                  wj=wj*dexp(-sij)
                  ENDIF
                  swj=swj+wj
                  swjy=swjy+wj*y(j1,j2)
               END DO
            END DO
            ai(i1,i2)=swjy
            bi(i1,i2)=swj
            bi0(i1,i2)=swj0
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant bivariate aws (gridded)
C   Heteroscedastic Gaussian case
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chawsbi(y,fix,si2,n1,n2,hakt,lambda,theta,bi,bi0,ai,
     1                   kern,spmax,wght)

C   
C   y        observed values of regression function
C   n1,n2    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second dimension (larger values shrink)
C   
      implicit logical (a-z)
      external lkern
      real*8 lkern
      integer n1,n2,kern
      logical aws,fix(n1,n2)
      real*8 y(n1,n2),hakt,theta(n1,n2),bi(n1,n2),ai(n1,n2),
     1       lambda,spmax,wght,bi0(n1,n2),si2(n1,n2)
      integer ih1,ih2,i1,i2,j1,j2,ja1,je1,ja2,je2
      real*8 thetai,bii,sij,swj,swjy,z1,z2,wj,hakt2,swj0,z
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i1=1,n1
         DO i2=1,n2
            IF (fix(i1,i2)) CYCLE
C    nothing to do, final estimate is already fixed by control 
            thetai=theta(i1,i2)
            bii=bi(i1,i2)/lambda
C   scaling of sij outside the loop
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            swj=0.d0
            swj0=0.d0
            swjy=0.d0
            DO j1=ja1,je1
               z1=(i1-j1)
               z1=z1*z1
               ih2=dsqrt(hakt2-z1)/wght
               ja2=max0(1,i2-ih2)
               je2=min0(n2,i2+ih2)
               DO j2=ja2,je2
C  first stochastic term
                  z2=(i2-j2)*wght
                  wj=lkern(kern,(z1+z2*z2)/hakt2)*si2(j1,j2)
                  swj0=swj0+wj  
                  IF (aws) THEN
                     z=thetai-theta(j1,j2)
                     sij=bii*z*z
                     IF (sij.gt.spmax) CYCLE
                     wj=wj*dexp(-sij)
                  END IF
                  swj=swj+wj
                  swjy = swjy+wj*y(j1,j2)
               END DO
            END DO
            ai(i1,i2)=swjy
            bi(i1,i2)=swj
            bi0(i1,i2)=swj0
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawstri(y,fix,n1,n2,n3,hakt,lambda,theta,bi,bi0,ai,
     1                   model,kern,spmax,wght)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldist,lkern
      real*8 kldist,lkern
      integer n1,n2,n3,model,kern
      logical aws,fix(n1,n2,n3)
      real*8 y(n1,n2,n3),theta(n1,n2,n3),bi(n1,n2,n3),bi0(n1,n2,n3),
     1       ai(n1,n2,n3),lambda,spmax,wght(2),hakt,bii0
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,ja1,je1,ja2,je2,ja3,je3
      real*8 thetai,bii,sij,swj,swj0,swjy,z1,z2,z3,wj,hakt2
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               IF (fix(i1,i2,i3)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               thetai=theta(i1,i2,i3)
               bii=bi(i1,i2,i3)/lambda
C   scaling of sij outside the loop
               bii0=bi0(i1,i2,i3)
               ja1=max0(1,i1-ih1)
               je1=min0(n1,i1+ih1)
               swj=0.d0
               swj0=0.d0
               swjy=0.d0
               DO j1=ja1,je1
                  z1=(i1-j1)
                  z1=z1*z1
                  ih2=dsqrt(hakt2-z1)/wght(1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
                  DO j2=ja2,je2
                     z2=(i2-j2)*wght(1)
                     z2=z2*z2
                     ih3=dsqrt(hakt2-z1-z2)/wght(2)
                     ja3=max0(1,i3-ih3)
                     je3=min0(n3,i3+ih3)
                     DO j3=ja3,je3
C  first stochastic term
                        z3=(i3-j3)*wght(2)
                        wj=lkern(kern,(z1+z2+z3*z3)/hakt2)
                        swj0=swj0+wj
                        IF (aws) THEN
                     sij=bii*kldist(model,thetai,theta(j1,j2,j3),bii0)
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*exp(-sij)
                        END IF
                        swj=swj+wj
                        swjy=swjy+wj*y(j1,j2,j3)
                     END DO
                  END DO
               END DO
               ai(i1,i2,i3)=swjy
               bi(i1,i2,i3)=swj
               bi0(i1,i2,i3)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chawstri(y,fix,si2,n1,n2,n3,hakt,lambda,theta,bi,bi0,
     1                    ai,kern,spmax,wght)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external lkern
      real*8 lkern
      integer n1,n2,n3,kern
      logical aws,fix(n1,n2,n3)
      real*8 y(n1,n2,n3),theta(n1,n2,n3),bi(n1,n2,n3),bi0(n1,n2,n3),
     1       ai(n1,n2,n3),lambda,spmax,si2(n1,n2,n3),
     2       wght(2),hakt
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,ja1,je1,ja2,je2,ja3,je3
      real*8 thetai,bii,sij,swj,swj0,swjy,z1,z2,z3,wj,hakt2,z
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               IF (fix(i1,i2,i3)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               thetai=theta(i1,i2,i3)
               bii=bi(i1,i2,i3)/lambda
C   scaling of sij outside the loop
               ja1=max0(1,i1-ih1)
               je1=min0(n1,i1+ih1)
               swj=0.d0
               swj0=0.d0
               swjy=0.d0
               DO j1=ja1,je1
                  z1=(i1-j1)
                  z1=z1*z1
                  ih2=dsqrt(hakt2-z1)/wght(1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
                  DO j2=ja2,je2
                     z2=(i2-j2)*wght(1)
                     z2=z2*z2
                     ih3=dsqrt(hakt2-z1-z2)/wght(2)
                     ja3=max0(1,i3-ih3)
                     je3=min0(n3,i3+ih3)
                     DO j3=ja3,je3
                        z3=(i3-j3)*wght(2)
                   wj=lkern(kern,(z1+z2+z3*z3)/hakt2)*si2(j1,j2,j3)
                        swj0=swj0+wj
                        IF (aws) THEN
                           z=thetai-theta(j1,j2,j3)
                           sij=bii*z*z
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*dexp(-sij)
                        END IF
                        swj=swj+wj
                        swjy=swjy+wj*y(j1,j2,j3)
                     END DO
                  END DO
               END DO
               ai(i1,i2,i3)=swjy
               bi(i1,i2,i3)=swj
               bi0(i1,i2,i3)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Local constant aws for irregular design
C
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawsmul(y,fix,x,n,dx,hakt,lambda,theta,bi,bi0,ai,
     1                   model,kern,spmax)
C    
C     n          number of design points
C     y          observed values at design points
C     hakt      actual number of nearest neighbors     
C     lambda   lambda or lambda*sigma2 for Gaussian models
C     theta    estimates from last step   (input)
C     bi       \sum  Wi   (output)
C     ai       \sum  Wi Y     (output)
C     model    specifies the probablilistic model for the KL-Distance
C     kern     specifies the location kernel
C     spmax    specifies the truncation point of the stochastic kernel
C     wght     scaling factor for second dimension (larger values shrink)
C     
      implicit logical (a-z)
      external kldist,lmkern
      real*8 kldist,lmkern
      integer n,dx,model,kern
      logical aws,fix(n)
      real*8 y(n),x(dx,n),hakt,theta(n),bi(n),bi0(n),ai(n),
     1       lambda,spmax
      integer i,j
      real*8 thetai,bii,sij,swj,swjy,wj,hakt2,swj0,bii0
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      hakt2=hakt*hakt
      aws=lambda.gt.100
      DO i=1,n
C        loop over design points
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         thetai=theta(i)
         bii=bi(i)/lambda
C   scaling of sij outside the loop
         bii0=bi0(i)
         swj=0.d0
         swj0=0.d0
         swjy=0.d0
         DO j=1,n
            wj=lmkern(kern,dx,x(1,i),x(1,j),hakt2)
            IF (wj.le.0.d0) CYCLE
C    thats the location penalty
C  first stochastic term
            swj0=swj0+wj
            IF (aws) THEN
               sij=bii*kldist(model,thetai,theta(j),bii0)
               IF (sij.gt.spmax) CYCLE
               wj=wj*dexp(-sij)
            END IF
            swj=swj+wj
            swjy = swjy+wj*y(j)
         END DO
         ai(i)=swjy
         bi(i)=swj
         bi0(i)=swj0
      END DO      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawsmuld(y,fix,x,n,nn,dx,hakt,lambda,theta,bi,bi0,ai,
     1                    model,kern,spmax)
C    
C     n          number of design points in TS
C     nn        number of design points in TS + new 
C     y          observed values at design points
C     hakt      actual number of nearest neighbors     
C     lambda   lambda or lambda*sigma2 for Gaussian models
C     theta    estimates from last step   (input)
C     bi       \sum  Wi   (output)
C     ai       \sum  Wi Y     (output) 
C     model    specifies the probablilistic model for the KL-Distance
C     kern     specifies the location kernel
C     spmax    specifies the truncation point of the stochastic kernel
C     wght     scaling factor for second dimension (larger values shrink)
C     
      implicit logical (a-z)
      external kldist,lmkern
      real*8 kldist,lmkern
      integer n,nn,dx,model,kern
      logical aws,fix(nn)
      real*8 y(nn),x(dx,nn),hakt,theta(nn),bi(nn),bi0(nn),ai(nn),
     1       lambda,spmax
      integer i,j
      real*8 thetai,bii,sij,swj,swjy,wj,hakt2,swj0,bii0
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      hakt2=hakt*hakt
      aws=lambda.lt.100
      DO i=1,nn
C        loop over design points
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         thetai=theta(i)
         bii=bi(i)/lambda
C   scaling of sij outside the loop
         bii0=dmax1(bi0(i),1.d0)
         swj=0.d0
         swj0=0.d0
         swjy=0.d0
         DO j=1,n
            wj=lmkern(kern,dx,x(1,i),x(1,j),hakt2)
            IF (wj.le.0.d0) CYCLE
C    thats the location penalty
C  first stochastic term
            swj0=swj0+wj
            IF (aws) THEN
               sij=bii*kldist(model,thetai,theta(j),bii0)
               IF (sij.gt.spmax) CYCLE
               wj=wj*dexp(-sij)
            END IF
            swj=swj+wj
            swjy=swjy+wj*y(j)
         END DO
         ai(i)=swjy
         bi(i)=swj
         bi0(i)=swj0
      END DO      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws 
C   Heteroskedastic Gaussian case
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chawsmul(y,fix,x,si2,n,dx,hakt,lambda,theta,bi,bi0,ai,
     1                    kern,spmax)
C    
C     n          number of design points
C     y          observed values at design points
C     hakt      actual number of nearest neighbors     
C     lambda   lambda or lambda*sigma2 for Gaussian models
C     theta    estimates from last step   (input)
C     bi       \sum  Wi   (output)
C     ai       \sum  Wi Y     (output)
C     model    specifies the probablilistic model for the KL-Distance
C     kern     specifies the location kernel
C     spmax    specifies the truncation point of the stochastic kernel
C     wght     scaling factor for second dimension (larger values shrink)
C     
      implicit logical (a-z)
      external lmkern
      real*8 lmkern
      integer n,dx,kern
      logical aws,fix(n)
      real*8 y(n),x(dx,n),hakt,theta(n),bi(n),bi0(n),ai(n),
     1       lambda,spmax,si2(n)
      integer i,j
      real*8 thetai,bii,sij,swj,swjy,z,wj,hakt2,swj0
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      hakt2=hakt*hakt
      aws=lambda.lt.1.e20
      DO i=1,n
C        loop over design points
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         thetai=theta(i)
         bii=bi(i)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         swj0=0.d0
         swjy=0.d0
         DO j=1,n
	    wj=lmkern(kern,dx,x(1,i),x(1,j),hakt2)*si2(j)
	    IF (wj.le.0.d0) CYCLE
C    thats the location penalty
C  first stochastic term
            swj0=swj0+wj
            IF (aws) THEN
	       z=thetai-theta(j)
               sij=bii*z*z
               IF (sij.gt.spmax) CYCLE
	       wj=wj*dexp(-sij)
	    END IF
            swj=swj+wj
            swjy = swjy+wj*y(j)
         END DO
         ai(i)=swjy
         bi(i)=swj
         bi0(i)=swj0
      END DO      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws (NN)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawsmnn(y,fix,nn,n,inmax,iakt,lambda,theta,bi,bi0,ai,
     1                   model,kern,spmax)
C    
C     y          observed values at design points
C     nn        matrix of nearest neighbors
C     n          number of design points
C     inmax    maximal number of nearest neighbors
C     iakt      actual number of nearest neighbors     
C     lambda   lambda or lambda*sigma2 for Gaussian models
C     theta    estimates from last step   (input)
C     bi       \sum  Wi   (output)
C     ai       \sum  Wi Y     (output)
C     model    specifies the probablilistic model for the KL-Distance
C     kern     specifies the location kernel
C     spmax    specifies the truncation point of the stochastic kernel
C     wght     scaling factor for second dimension (larger values shrink)
C     
      implicit logical (a-z)
      external kldist,lkern
      real*8 kldist,lkern
      integer n,inmax,nn(inmax,n),iakt,model,kern
      logical aws,fix(n)
      real*8 y(n),theta(n),bi(n),bi0(n),ai(n),
     1       lambda,spmax
      integer i,j,jj
      real*8 thetai,bii,sij,swj,swj0,swjy,z,wj,hakt,bii0
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      hakt=iakt+.5
      aws=lambda.lt.1e20
      DO i=1,n
C        loop over design points        
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         thetai=theta(i)
         bii=bi(i)/lambda
C   scaling of sij outside the loop
         bii0=dmax1(bi0(i),1.d0)
         swj=0.d0
         swj0=0.d0
         swjy=0.d0
         DO jj=1,iakt
	    j=nn(jj,i)
            z=(jj-1)/hakt
	    wj=lkern(kern,z*z)
            swj0=swj0+wj
	    IF (aws) THEN
               sij=bii*kldist(model,thetai,theta(j),bii0)
               IF (sij.gt.spmax) CYCLE
	       wj=wj*dexp(-sij)
	    END IF
            swj=swj+wj
            swjy = swjy+wj*y(j)
         END DO
         ai(i)=swjy
         bi0(i)=swj0
         bi(i)=swj
      END DO      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local constant aws (NN)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chawsmnn(y,fix,si2,nn,n,inmax,iakt,lambda,theta,bi,
     1                    bi0,ai,kern,spmax)
C    
C     y          observed values at design points
C     nn        distance matrix 
C     n          number of design points
C     inmax    maximal number of nearest neighbors
C     iakt      actual number of nearest neighbors     
C     lambda   lambda or lambda*sigma2 for Gaussian models
C     theta    estimates from last step   (input)
C     bi       \sum  Wi   (output)
C     ai       \sum  Wi Y     (output)
C     model    specifies the probablilistic model for the KL-Distance
C     kern     specifies the location kernel
C     spmax    specifies the truncation point of the stochastic kernel
C     wght     scaling factor for second dimension (larger values shrink)
C     
      implicit logical (a-z)
      external lkern
      real*8 lkern
      integer n,inmax,nn(inmax,n),iakt,kern
      logical aws,fix(n)
      real*8 y(n),theta(n),bi(n),bi0(n),ai(n),lambda,spmax,si2(n)
      integer i,j,jj
      real*8 thetai,bii,sij,swj,swj0,swjy,wj,hakt,z
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      hakt=iakt+.5
      aws=lambda.lt.1e20
      DO i=1,n
C        loop over design points
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         thetai=theta(i)
         bii=bi(i)/lambda
C   scaling of sij outside the loop
         swj=0.d0
         swj0=0.d0
         swjy=0.d0
         DO jj=1,iakt
	    j=nn(jj,i)
            z=(i-j)/hakt
	    wj=lkern(kern,z*z)*si2(j)
            swj0=swj0+wj
	    IF (aws) THEN
	       z=thetai-theta(j)
               sij=bii*z*z
               IF (sij.gt.spmax) CYCLE
	       wj=wj*dexp(-sij)
	    END IF
            swj=swj+wj
            swjy = swjy+wj*y(j)
         END DO
         ai(i)=swjy
         bi0(i)=swj0
         bi(i)=swj
      END DO      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Multivariate local polynomial ( p=0 or 1 ) aws ( p=0 or 1 )
C
C   Nongridded design and Nearest Neighbor approach
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local polynomial aws 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cpawsmul(n,px,dp1,dp2,x,y,fix,nn,distm,ihakt,theta,
     1                    bi,bi0,ai,lam,h,kern,
     2                    dmat,thij,psix,psiy,spmax)
C    
C     n          number of design points
C     px         dimension of x
C     dp1        number of parameters  (1+p*px)
C     dp2        number of components in bi  (dp1*(dp1+1)/2)
C     x          design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     distm      distances ordered by nearest neighbors
C     ihakt      initial number of nearest neighbors     
C     theta      old estimates from last step (input)
C     bi        \sum \Psi^T Wi^k \Psi        (input/output)
C     bi0        \sum \Psi^T Wi0 \Psi    (input/output)
C     ai         \sum \Psi^T Wi^k Y           (input/output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     dmat       working arrays  dp1 times dp1
C     thij                     vector for parameter differences
C     psix,psiy,xij                working memory 
C     
      implicit logical (a-z)
      integer n,dp1,dp2,px,i,j,k,l,ja,nwij,
     1        ihakt,m,nn(ihakt,n),kern      
      logical aws,fix(n),pos
      real*8 x(px,n),y(n),psix(dp1),psiy(dp1),theta(dp1,n),
     1 bi(dp2,n),ai(dp1,n),lam,thij(dp1),
     2 dmat(dp1,dp1),distm(ihakt,n),lkern,skern,
     3 bi0(dp2,n),wijl,lambda,h,spmax,d,wij,z,xij,s0i,ha,ha2,eps
      external lkern,skern
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      eps=.2
      ha=h
      ha2=ha*ha
      aws=lam.le.1d20
      DO i=1,n
         IF (fix(i)) CYCLE
C        loop over design points
C  disable extension penalty if estimate is very unstable
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmat
C    expand Bi in dmat
         s0i=bi0(1,i)
         m=1
         DO k=1,dp1
            DO l=1,k
               dmat(l,k)=bi(m,i)
               dmat(k,l)=bi(m,i)
               m=m+1
	    END DO
	 END DO
         DO l=1,dp2
            bi0(l,i)=0.d0
         END DO
8999     nwij=0
         DO l=1,dp2
            bi(l,i)=0.d0
            IF (l.gt.dp1) CYCLE
            ai(l,i)=0.d0
         END DO
C      if not enough points with positive weights (counted in nwij)
C      lambda will be increased to reduce panalization by stochastic term
         DO j=1,ihakt
            z=distm(j,i)
            IF (z.ge.ha) CYCLE
            wijl=lkern(kern,z*z/ha2)
	    wij=wijl
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
            ja=nn(j,i)
            DO l=1,dp1
               thij(l)=theta(l,ja)
            END DO
            psix(1)=1
            psiy(1)=y(ja)
            IF (dp1.gt.1) THEN
               DO l=2,dp1
                  xij=(x(l-1,i)-x(l-1,ja))
                  thij(1)=thij(1)+xij*theta(l,ja)
                  psix(l)=-xij
                  psiy(l)=-xij*y(ja)
	       END DO        
	    END IF
C     thij contais theta_{j} in the model centered at xi
C
C     now we have everything to compute  w_{ij}
C
	    IF (aws) THEN
               d=0.0d0
               DO l=1,dp1
                  thij(l)=theta(l,i)-thij(l)
               END DO
C
C   thats the difference between thetai and thetaij
C
               DO l=1,dp1
                  DO k=1,dp1
                     d=d+dmat(k,l)*thij(l)*thij(k)
                  END DO
	       END DO
               wij=wij*skern(d/lam,spmax)
	    END IF
            IF (wij.gt.eps) nwij=nwij+1
C        now compute contributions to bi(i),ai(i)  
            m=1
	    pos=wij.gt.0.d0
            DO k=1,dp1
               IF (pos) ai(k,i)=ai(k,i)+wij*psiy(k)
               DO l=1,k
	          z=psix(k)*psix(l)
                  IF (pos) bi(m,i)=bi(m,i)+wij*z
                  bi0(m,i)=bi0(m,i)+wijl*z
                  m=m+1
               END DO 
            END DO                 
         END DO 
C    this was the j - loop
         IF (nwij.ge.dp1)  CYCLE
         lambda=lambda*1.5
C    increase lambda to weaken stochastic penalization
         goto 8999
      END DO      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in multivariate local polynomial aws (NN)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cpawsmnn(n,px,dp1,dp2,x,y,fix,nn,ihakt,theta,bi,bi0,
     1  ai,lam,h,kern,dmat,thij,psix,psiy,spmax)
C    
C     n          number of design points
C     px         dimension of x
C     dp1        number of parameters  (1+p*px)
C     dp2        number of components in bi  (dp1*(dp1+1)/2)
C     x          design points
C     y          observed values at design points
C     nn         indices of ihinit nearest neighbors
C     ihakt      initial number of nearest neighbors     
C     theta      old estimates from last step (input)
C     bi        \sum \Psi^T Wi^k \Psi        (output)
C     bi0        \sum \Psi^T Wi0 \Psi    (input/output)
C     ai         \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     dmat, dmat      working arrays  dp1 times dp1
C     thij                     vector for parameter differences
C     psix,psiy,xij                working memory 
C     
      implicit logical (a-z)
      integer n,dp1,dp2,px,i,j,k,l,ja,nwij,ij,
     1        ihakt,m,nn(ihakt,n),kern      
      logical aws,fix(n),pos
      real*8 x(px,n),y(n),psix(dp1),psiy(dp1),theta(dp1,n),
     1 bi(dp2,n),ai(dp1,n),lam,lkern,skern,
     2 dmat(dp1,dp1),thij(dp1),spmax,
     3 bi0(dp2,n),pij,wijl,lambda,h,
     4 d,wij,z,xij,s0i,ha,ha2,eps
      external lkern,skern
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      
      eps=.2
      ha=h-.9d0
      ha2=ha*ha
      aws=lam.le.1d20
      DO i=1,n
         IF (fix(i)) CYCLE
C        loop over design points
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmat 
C    expand Bi in dmat
         s0i=bi0(1,i)
         m=1
         DO k=1,dp1
            DO l=1,k
               dmat(l,k)=bi(m,i)
               dmat(k,l)=bi(m,i)
               m=m+1
            END DO
         END DO
         DO l=1,dp2
            bi0(l,i)=0.d0
         END DO
8999     nwij=0
         DO l=1,dp2
            bi(l,i)=0.d0
            IF (l.gt.dp1) CYCLE
            ai(l,i)=0.d0
         END DO
C      if not enough points with positive weights (counted in nwij)
C      lambda will be increased to reduce panalization by stochastic term)
         DO j=1,ihakt
            ja=nn(j,i)
            ij=j-1
            pij=ij*ij/ha2
	    wijl=lkern(kern,pij)
	    wij=wijl
C    thats spatial penalization only, now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
            DO l=1,dp1
               thij(l)=theta(l,ja)
            END DO
            psix(1)=1
            psiy(1)=y(ja)
            IF (dp1.gt.1) THEN
               z=1.d0
               DO l=2,dp1
                  xij=(x(l-1,i)-x(l-1,ja))
                  thij(1)=thij(1)+xij*theta(l,ja)
                  psix(l)=-xij
                  psiy(l)=-xij*y(ja)
               END DO
	    END IF
C     thij contains theta_{j} in the model centered at xi
C
C     now we have everything to compute  w_{ij}
C
	    IF (aws) THEN
               d=0.0d0
               DO l=1,dp1
                  thij(l)=theta(l,i)-thij(l)
               END DO
C
C   thats the difference between thetai and thetaij
C
               DO l=1,dp1
                  DO k=1,dp1
                     d=d+dmat(k,l)*thij(l)*thij(k)
                  END DO    
               END DO 
               wij=wij*skern(d/lam,spmax)
	    END IF
            IF (wij.gt.eps) nwij=nwij+1
C        now compute contributions to bi(i),ai(i)  
            m=1
	    pos=wij.gt.0.d0
            DO k=1,dp1
               IF (pos) ai(k,i)=ai(k,i)+wij*psiy(k)
               DO l=1,k
	          z=psix(k)*psix(l)
                  IF (pos) bi(m,i)=bi(m,i)+wij*z
                  bi0(m,i)=bi0(m,i)+wijl*z
                  m=m+1
               END DO    
            END DO 
	 END DO
C    this was the j - loop
         IF (nwij.ge.dp1) CYCLE
         lambda=lambda*1.25
C    increase lambda and tau to weaken stochastic and influence penalization
         goto 8999
      END DO      
      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (multivariate polynomial aws)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mpawsmul(n,dp1,dp2,ai,bi,theta,dmat)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi    
C     theta      new parameter estimate
C     dmat       working array
C
      integer n,dp1,dp2
      real*8 ai(dp1,n),bi(dp2,n),theta(dp1,n),dmat(dp1,dp1)
      integer i,j,k,l,m,info
      real*8 d
      DO i=1,n
         l=1
         m=1
         DO j=1,dp1
            DO k=1,j
               dmat(k,j)=bi(m,i)
               m=m+1
            END DO
	 END DO
         call invers(dmat,dp1,info)
         IF (info.ne.0) CYCLE 
C     just keep the old estimate
C     now dmat contains inverse of B_i 
C     now calculate theta as B_i^{-1} A_i
         DO j=1,dp1
            d=0.0d0
            DO k=1,dp1
               d=d+dmat(j,k)*ai(k,i)  
            END DO
            theta(j,i)=d
         END DO
      END DO
      RETURN
      END      
