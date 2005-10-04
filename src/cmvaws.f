CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cmvaws(y,fix,sigma2,df,n1,n2,n3,hakt,lambda,theta,
     1                  s2hat,bi,bi2,bi0,ai,si,kern,spmax,wght)
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
      external kldistmv,lkern
      real*8 kldistmv,lkern
      integer n1,n2,n3,kern
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,spmax,wght(2),
     1       bi2(1),hakt,df,sigma2(1),s2hat(1),si(1)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,ja1,je1,ja2,je2,ja3,je3,
     1        iind,jind,jind3,jind2
      real*8 thetai,bii,sij,swj,swj2,swj0,swjy,swjv,z1,z2,z3,wj,hakt2,
     1        bii0,s2i
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
	       iind=i1+(i2-1)*n1+(i3-1)*n1*n2
               IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               thetai=theta(iind)
	       s2i=s2hat(iind)
               bii=bi(iind)/lambda
C   scaling of sij outside the loop
               bii0=bi0(iind)
	       ih3=hakt/wght(2)
               ja3=max0(1,i3-ih3)
               je3=min0(n3,i3+ih3)
               swj=0.d0
	       swj2=0.d0
               swj0=0.d0
               swjy=0.d0
	       swjv=0.d0
               DO j3=ja3,je3
                  z3=(i3-j3)*wght(2)
                  z3=z3*z3
                  ih2=dsqrt(hakt2-z3)/wght(1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
		  jind3=(j3-1)*n1*n2
                  DO j2=ja2,je2
                     z2=(i2-j2)*wght(1)
                     z2=z3+z2*z2
                     ih1=dsqrt(hakt2-z2)
                     ja1=max0(1,i1-ih1)
                     je1=min0(n1,i1+ih1)
		     jind2=jind3+(j2-1)*n1
                     DO j1=ja1,je1
C  first stochastic term
                        jind=j1+jind2
                        z1=(i1-j1)
                        wj=lkern(kern,(z1*z1+z2)/hakt2)
                        swj0=swj0+wj
                        IF (aws) THEN
                           sij=bii*kldistmv(thetai,theta(jind),
     1                                      s2i,s2hat(jind),df)
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*exp(-sij)
                        END IF
                        swj=swj+wj
                        swj2=swj2+wj*wj
                        swjy=swjy+wj*y(jind)
			swjv=swjv+wj*sigma2(jind)
                     END DO
                  END DO
               END DO
               ai(iind)=swjy
	       si(iind)=swjv
               bi(iind)=swj
               bi2(iind)=swj2
               bi0(iind)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute the Kullback-Leibler Distance
C
C          Model=1    Gaussian   
C          Model=2    Bernoulli   
C          Model=3    Poisson   
C          Model=4    Exponential   
C
C     computing dlog(theta) and dlog(1.d0-theta) outside the AWS-loops 
C     will reduces computational costs at the price of readability
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function kldistmv(thi,thj,s2i,s2j,df)
      implicit logical (a-z)
      real*8 thi,thj,s2i,s2j,df,zm,zv
C        Gaussian
         zm=thi-thj
	 zv=s2i/s2j
         kldistmv=zm*zm/s2j+df*(-dlog(zv)-1+zv)
      RETURN
      END
