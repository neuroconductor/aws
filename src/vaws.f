C
C    Copyright (C) 2004 Weierstrass-Institut fuer 
C                       Angewandte Analysis und Stochastik (WIAS)
C
C    Author:  Joerg Polzehl
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
C  USA.
C
C  The following routines are part of the aws package and contain  
C  FORTRAN 77 code needed in R function aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Local constant aws on a grid      
C
C   this is a reimplementation of the original aws procedure
C
C   should be slightly slower for non-Gaussian models (see function )
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   Local constant aws on a grid      
C
C   this is a reimplementation of the original aws procedure
C
C   should be slightly slower for non-Gaussian models (see function kldist)
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        one iteration of univariate local constant aws on a grid
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vawsuni(y,fix,d,n,hakt,lambda,theta,bi,bi0,ai,
     1                   model,kern,spmax,thetai,swjy,vw)
C   
C   y        observed values of regression function
C   d        number of components in y
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
      external vkldist,lkern
      real*8 vkldist,lkern
      integer n,d,model,kern
      logical aws,fix(n)
      real*8 y(d,n),hakt,theta(d,n),bi(n),bi0(n),ai(d,n),lambda,spmax,
     1       thetai(d),swjy(d),vw(d)
      integer ih,i,j,ja,je,k
      real*8 bii,sij,swj,z,wj,swj0,bii0
      ih=hakt
      aws=lambda.lt.1d40
      DO i=1,n
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         DO k=1,d
            thetai(k)=theta(k,i)
            swjy(k)=0.d0
	 END DO
         bii=bi(i)/lambda
C   scaling of sij outside the loop
         bii0=bi0(i)
         ja=max0(1,i-ih)
         je=min0(n,i+ih)
         swj=0.d0
         swj0=0.d0
         DO j=ja,je
C  first stochastic term
            z=(i-j)/hakt
            wj=lkern(kern,z*z)
            swj0=swj0+wj
            IF (aws) THEN
               sij=bii*vkldist(model,d,thetai,theta(1,j),bii0,vw)
               IF (sij.gt.spmax) CYCLE
               wj=wj*dexp(-sij)
            ENDIF
            swj=swj+wj
	    DO k=1,d
               swjy(k) = swjy(k)+wj*y(k,j)
	    END DO
         END DO
	 DO k=1,d
            ai(k,i)=swjy(k)
	 END DO
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
      subroutine vhawsuni(y,fix,si2,d,n,hakt,lambda,theta,bi,bi0,ai,
     1                    kern,spmax,thetai,swjy,vw)
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
      external vkldist,lkern
      real*8 vkldist,lkern
      integer n,kern,d,model
      logical aws,fix(n)
      real*8 y(d,n),hakt,theta(d,n),bi(n),bi0(n),ai(d,n),lambda,
     1       si2(n),spmax,thetai(d),swjy(d),vw(d)
      integer ih,i,j,ja,je,k
      real*8 bii,bii0,sij,swj,z,wj,swj0
      model=1
      ih=hakt
      aws=lambda.lt.1d40
      DO i=1,n
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         DO k=1,d
            thetai(k)=theta(k,i)
            swjy(k)=0.d0
	 END DO
	 bii0=bi0(i)
         bii=bi(i)/lambda
C   scaling of sij outside the loop
         ja=max0(1,i-ih)
         je=min0(n,i+ih)                  
         swj=0.d0
         swj0=0.d0
         DO j=ja,je
C  first stochastic term
            z=(i-j)/hakt
            wj=lkern(kern,z*z)*si2(j)
            swj0=swj0+wj
            IF (aws) THEN
               sij=bii*vkldist(model,d,thetai,theta(1,j),bii0,vw)
               IF (sij.gt.spmax) CYCLE
               wj=wj*dexp(-sij)
            ENDIF
            swj=swj+wj
	    DO k=1,d
               swjy(k) = swjy(k)+wj*y(k,j)
	    END DO
         END DO
	 DO k=1,d
            ai(k,i)=swjy(k)
	 END DO
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
      subroutine vawsbi(y,fix,d,n1,n2,hakt,lambda,theta,bi,bi0,ai,
     1                  model,kern,spmax,wght,thetai,swjy,vw)
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
      external vkldist,lkern
      real*8 vkldist,lkern
      integer d,n1,n2,model,kern
      logical aws,fix(n1,n2)
      real*8 y(d,n1,n2),hakt,theta(d,n1,n2),bi(n1,n2),ai(d,n1,n2),
     1       lambda,spmax,wght,bi0(n1,n2),thetai(d),swjy(d)
      integer ih1,ih2,i1,i2,j1,j2,ja1,je1,ja2,je2,k
      real*8 bii,sij,swj,z1,z2,wj,hakt2,swj0,bii0,vw(d)
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO  i1=1,n1
         DO  i2=1,n2
            IF (fix(i1,i2)) CYCLE
C    nothing to do, final estimate is already fixed by control 
            DO k=1,d
               thetai(k)=theta(k,i1,i2)
               swjy(k)=0.d0
	    END DO
            bii=bi(i1,i2)/lambda
C   scaling of sij outside the loop
            bii0=bi0(i1,i2)
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            swj=0.d0
            swj0=0.d0
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
                      sij=bii*vkldist(model,d,thetai,
     1                                theta(1,j1,j2),bii0,vw)
                  IF (sij.gt.spmax) CYCLE
                  wj=wj*dexp(-sij)
                  ENDIF
                  swj=swj+wj
		  DO k=1,d
                     swjy(k)=swjy(k)+wj*y(k,j1,j2)
		  END DO
               END DO
            END DO
	    DO k=1,d
               ai(k,i1,i2)=swjy(k)
	    END DO
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
      subroutine vhawsbi(y,fix,si2,d,n1,n2,hakt,lambda,theta,bi,bi0,
     1                   ai,kern,spmax,wght,thetai,swjy,vw)

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
      external vkldist,lkern
      real*8 vkldist,lkern
      integer d,n1,n2,kern
      logical aws,fix(n1,n2)
      real*8 y(d,n1,n2),hakt,theta(d,n1,n2),bi(n1,n2),ai(d,n1,n2),
     1       lambda,spmax,wght,bi0(n1,n2),si2(n1,n2),vw(d)
      integer ih1,ih2,i1,i2,j1,j2,ja1,je1,ja2,je2,k,model
      real*8 thetai(d),bii,bii0,sij,swj,swjy(d),z1,z2,wj,hakt2,swj0
      model=1
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i1=1,n1
         DO i2=1,n2
            IF (fix(i1,i2)) CYCLE
C    nothing to do, final estimate is already fixed by control 
            DO k=1,d
               thetai(k)=theta(k,i1,i2)
               swjy(k)=0.d0
	    END DO
            bii=bi(i1,i2)/lambda
            bii0=bi0(i1,i2)
C   scaling of sij outside the loop
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            swj=0.d0
            swj0=0.d0
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
                     sij=bii*vkldist(model,d,thetai,
     1                               theta(1,j1,j2),bii0,vw)
                  IF (sij.gt.spmax) CYCLE
                  wj=wj*dexp(-sij)
                  ENDIF
                  swj=swj+wj
		  DO k=1,d
                     swjy(k)=swjy(k)+wj*y(k,j1,j2)
		  END DO
               END DO
            END DO
	    DO k=1,d
               ai(k,i1,i2)=swjy(k)
	    END DO
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
      subroutine vawstri(y,fix,d,n1,n2,n3,hakt,lambda,theta,bi,bi0,ai,
     1                   model,kern,spmax,wght,thetai,swjy,vw)
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
      external vkldist,lkern
      real*8 vkldist,lkern
      integer d,n1,n2,n3,model,kern
      logical aws,fix(n1,n2,n3)
      real*8 y(d,n1,n2,n3),theta(d,n1,n2,n3),bi(n1,n2,n3),vw(d),
     1       bi0(n1,n2,n3),ai(d,n1,n2,n3),lambda,spmax,wght(2),hakt
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,ja1,je1,ja2,je2,ja3,je3,k
      real*8 thetai(d),bii,sij,swj,swj0,swjy(d),z1,z2,z3,wj,hakt2,bii0
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               IF (fix(i1,i2,i3)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               DO k=1,d
                  thetai(k)=theta(k,i1,i2,i3)
		  swjy(k)=0.d0
	       END DO
               bii=bi(i1,i2,i3)/lambda
C   scaling of sij outside the loop
               bii0=bi0(i1,i2,i3)
               ja1=max0(1,i1-ih1)
               je1=min0(n1,i1+ih1)
               swj=0.d0
               swj0=0.d0
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
                           sij=bii*vkldist(model,d,thetai,
     1                                     theta(1,j1,j2,j3),bii0,vw)
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*exp(-sij)
                        END IF
                        swj=swj+wj
			DO k=1,d
                           swjy(k)=swjy(k)+wj*y(k,j1,j2,j3)
			END DO
                     END DO
                  END DO
               END DO
	       DO k=1,d
                  ai(k,i1,i2,i3)=swjy(k)
	       END DO
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
      subroutine vhawstri(y,fix,si2,d,n1,n2,n3,hakt,lambda,theta,bi,
     1                    bi0,ai,kern,spmax,wght,thetai,swjy,vw)
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
      external vkldist,lkern
      real*8 vkldist,lkern
      integer d,n1,n2,n3,kern,model
      logical aws,fix(n1,n2,n3)
      real*8 y(d,n1,n2,n3),theta(d,n1,n2,n3),bi(n1,n2,n3),
     1       bi0(n1,n2,n3),ai(d,n1,n2,n3),lambda,spmax,si2(n1,n2,n3),
     2       wght(2),hakt,thetai(d),swjy(d),vw(d)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,ja1,je1,ja2,je2,ja3,je3,k
      real*8 bii,bii0,sij,swj,swj0,z1,z2,z3,wj,hakt2
      model=1
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               IF (fix(i1,i2,i3)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               DO k=1,d
                  thetai(k)=theta(k,i1,i2,i3)
		  swjy(k)=0.d0
	       END DO
               bii=bi(i1,i2,i3)/lambda
               bii0=bi0(i1,i2,i3)
C   scaling of sij outside the loop
               ja1=max0(1,i1-ih1)
               je1=min0(n1,i1+ih1)
               swj=0.d0
               swj0=0.d0
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
                           sij=bii*vkldist(model,d,thetai,
     1                                     theta(1,j1,j2,j3),bii0,vw)
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*exp(-sij)
                        END IF
                        swj=swj+wj
			DO k=1,d
                           swjy(k)=swjy(k)+wj*y(k,j1,j2,j3)
			END DO
                     END DO
                  END DO
               END DO
	       DO k=1,d
                  ai(k,i1,i2,i3)=swjy(k)
	       END DO
               bi(i1,i2,i3)=swj
               bi0(i1,i2,i3)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vpawsuni(n,d,dp1,dp2,x,y,fix,theta,bi,bi0,
     1         ai,lam,h,kern,cb,dmat,thij,psix,psiy,bii,spmax,vw)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     x          design points (ordered)
C     y          observed values at design points
C     theta      old estimates from last step (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     si0        \sum \Psi^T Wi0 \Psi [1,]    (input)
C     bi0        \sum \Psi^T Wi0 \Psi    (input/output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     cb         binomial coefficients     
C     dmat,          working arrays  dp1 times dp1
C     thij                     vector for parameter differences
C     psix,psiy                working memory 
C     
      implicit logical (a-z)
      integer n,d,dp1,dp2,i,j,k,l,m,je,ja,nwij,kern 
      logical aws,fix(n)
      real*8 lkern,x(n),y(d,n),psix(dp2),psiy(dp1,d),theta(dp1,d,n),
     1 bi(dp2,n),bii(dp2),ai(dp1,d,n),lam,spmax,dmat(dp1,dp1),
     2 thij(dp1,d),cb(dp1,dp1),bi0(dp2,n),lambda,h,vw(d),
     3 sij,wij,z,xij,ha,ha2,eps,xi,zj,thijl
      external lkern,skern
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      eps=1.e-10
      aws=lam.gt.1.d20
      DO i=1,n
C        loop over design points, initialization
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
	 lambda=lam
         ha=h
         xi=x(i)
C     first search for ja and je
1099     DO j=i,1,-1
            IF (x(j).le.xi-ha) EXIT
            ja=j
         END DO
         DO j=i,n
            IF (x(j).gt.xi+ha) EXIT
            je=j
         END DO
         IF (je-ja-1.le.dp1) THEN
            ha=ha*1.25
C            not enough points in neighborhood to estimate parameters
C            increase ha
            goto 1099
         END IF
         ha2=ha*ha
8999     nwij=0
         IF (aws) THEN
C
C            prepare for computation of stochastic and extension penalty 
C
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmat
C    expand Bi in dmat
            l=1
            DO j=1,dp1
               DO k=1,dp1
                  dmat(j,k)=bi(j+k-1,i)/lambda
	       END DO
            END DO
         END IF
         DO l=1,dp2
            bii(l)=0.d0
            bi0(l,i)=0.d0
            IF (l.gt.dp1) CYCLE
	    DO m=1,d
               ai(l,m,i)=0.d0
	    end do
         END DO
C      if not enough points with positive weights (counted in nwij)
C      lambda will be increased 
C      (to reduce panalization by stochastic and influence term)
C
C             loop over local neighborhood
C
         DO j=ja,je
C
C              get location weights
C
            xij=(x(j)-x(i))
            wij=lkern(kern,xij*xij/ha2)
            z=1.d0
            DO k=1,dp2
               psix(k)=z
               IF (k.le.dp1) THEN
	          DO m=1,d
	             psiy(k,m)=z*y(m,j)
		  END DO
	       END IF
               z=z*xij
            END DO
            DO k=1,dp2
               bi0(k,i)=bi0(k,i)+wij*psix(k)
	    END DO
            IF (aws) THEN
C    now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
               DO l=1,dp1
	          DO m=1,d
                     thij(l,m)=theta(l,m,j)
		  END DO
               END DO
               IF (dp1.gt.1) THEN
                  z=1.d0
                  zj=1.d0
                  DO l=2,dp1
                     z=-z*xij
                     DO k=1,dp1-l+1
		        DO m=1,d
                           thij(k,m)=thij(k,m)+
     1                               cb(k+l-1,k)*z*theta(k+l-1,m,j)
			END DO
                     END DO
                  END DO
               END IF
C     thij contains theta_{j} in the model centered at xi
C
C                 now get sij
C
	       DO l=1,dp1
	          DO m=1,d
                     thij(l,m)=theta(l,m,i)-thij(l,m)
		  END DO
               END DO
C
C   thats the difference between thetai and thetaij
C
               sij=0.0d0
               DO l=1,dp1
		  DO m=1,d
	             thijl=thij(l,m)
		     sij=sij+dmat(l,l)*thijl*thijl*vw(m)
		     IF (l.eq.dp1) CYCLE
                     DO k=l+1,dp1
                        sij=sij+2.d0*dmat(k,l)*thijl*thij(k,m)*vw(m)
                     END DO
		  END DO
               END DO
C
C     now we have everything to compute  w_{ij}
C       
               IF (sij.gt.spmax) CYCLE
               wij=wij*dexp(-sij)
            END IF
            IF (wij.gt.eps) nwij=nwij+1
C        now compute contributions to bi(i),bi0(i),ai(i)  
            DO k=1,dp2
               bii(k)=bii(k)+wij*psix(k)
               IF (k.le.dp1) THEN
	          DO m=1,d
	             ai(k,m,i)=ai(k,m,i)+wij*psiy(k,m)
		  END DO
	       END IF
            END DO    
         END DO    
C    this was the j - loop
         IF (nwij.lt.dp1.and.aws) THEN
            lambda=lambda*1.25
C    increase lambda to weaken stochastic penalization
            goto 8999
	 END IF
	 DO k=1,dp2
	    bi(k,i)=bii(k)
	 END DO
      END DO      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vphawsun(n,d,dp1,dp2,x,y,fix,si2,theta,bi,si0,
     1   bi0,ai,lam,h,kern,cb,dmat,thij,psix,psiy,bii,spmax,vw)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     x          design points (ordered)
C     y          observed values at design points
C     theta      old estimates from last step (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     si0        \sum \Psi^T Wi0 \Psi [1,]    (input)
C     bi0        \sum \Psi^T Wi0 \Psi    (input/output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     cb         binomial coefficients     
C     dmat,        working arrays  dp1 times dp1
C     thij                vector for parameter differences
C     psix,psiy                working memory 
C     
C      implicit logical (a-z)
      integer n,d,dp1,dp2,i,j,k,l,je,ja,nwij,kern,m
      logical aws,fix(n)
      real*8 x(n),y(d,n),psix(dp2),psiy(dp1,d),theta(dp1,d,n),
     1 bi(dp2,n),bii(dp2),ai(dp1,d,n),spmax,si2(n),dmat(dp1,dp1),
     2 thij(dp1,d),cb(dp1,dp1),bi0(dp2,n),lambda,h,si0(n),vw(d),
     4 sij,wij,z,xij,ha,ha2,eps,xi,zj,wijs,thijl,lkern,skern,lam
      external lkern,skern
      eps=1.e-10
      aws=lam.gt.1.d20
      DO i=1,n
C        loop over design points, initialization
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
	 lambda=lam
         ha=h
         xi=x(i)
C     first search for ja and je
1099     DO j=i,1,-1
            IF (x(j).le.xi-ha) EXIT
            ja=j
         END DO
         DO j=i,n
            IF (x(j).gt.xi+ha) EXIT
            je=j
         END DO
         IF (je-ja-1.le.dp1) THEN
            ha=ha*1.25
C            not enough points in neighborhood to estimate parameters
C            increase ha
            goto 1099
         END IF
         ha2=ha*ha
8999     nwij=0
          IF (aws) THEN
C
C            prepare for computation of stochastic and extension penalty 
C
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmat 
C    expand Bi in dmat
           l=1
            DO j=1,dp1
               DO k=1,dp1
                  dmat(j,k)=bi(j+k-1,i)/lambda
	       END DO
            END DO
         END IF
         DO l=1,dp2
            bii(l)=0.d0
            bi0(l,i)=0.d0
            IF (l.gt.dp1) CYCLE
	    DO m=1,d
               ai(l,m,i)=0.d0
	    END DO
         END DO
C      if not enough points with positive weights (counted in nwij)
C      lambda will be increased 
C      (to reduce panalization by stochastic and influence term)
C
C             loop over local neighborhood
C
         DO j=ja,je
C
C              get location weights
C
            xij=(x(j)-x(i))
            wij=lkern(kern,xij*xij/ha2)*si2(j)
	    wijs=1.d0
            z=1.d0
            DO k=1,dp2
               psix(k)=z
               IF (k.le.dp1) THEN
	          DO m=1,d
	             psiy(k,m)=z*y(m,j)
		  END DO
	       END IF
               z=z*xij
            END DO
            DO k=1,dp2
               bi0(k,i)=bi0(k,i)+wij*psix(k)
	    END DO
            IF (aws) THEN
C    now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
               DO l=1,dp1
	          DO m=1,d
                     thij(l,m)=theta(l,m,j)
		  END DO
               END DO
               IF (dp1.gt.1) THEN
                  z=1.d0
                  zj=1.d0
                  DO l=2,dp1
                     z=-z*xij
                     DO k=1,dp1-l+1
		        DO m=1,d
                           thij(k,m)=thij(k,m)+
     1                               cb(k+l-1,k)*z*theta(k+l-1,m,j)
			END DO
                     END DO
                  END DO
               END IF
C     thij contains theta_{j} in the model centered at xi
C
C                 now get sij
C
	       DO l=1,dp1
	          DO m=1,d
                     thij(l,m)=theta(l,m,i)-thij(l,m)
		  END DO
               END DO
C
C   thats the difference between thetai and thetaij
C
               sij=0.0d0
               DO l=1,dp1
		  DO m=1,d
	             thijl=thij(l,m)
		     sij=sij+dmat(l,l)*thijl*thijl*vw(m)
		     IF (l.eq.dp1) CYCLE
                     DO k=l+1,dp1
                        sij=sij+2.d0*dmat(k,l)*thijl*thij(k,m)*vw(m)
                     END DO
		  END DO
               END DO
C
C     now we have everything to compute  w_{ij}
C       
               IF (sij.gt.spmax) CYCLE
               wij=wij*dexp(-sij)
            END IF
            IF (wij.gt.eps) nwij=nwij+1
C        now compute contributions to bi(i),bi0(i),ai(i)  
            DO k=1,dp2
               bii(k)=bii(k)+wij*psix(k)
               IF (k.le.dp1) THEN
	          DO m=1,d
	             ai(k,m,i)=ai(k,m,i)+wij*psiy(k,m)
		  END DO
	       END IF
            END DO    
         END DO    
C    this was the j - loop
         IF (nwij.lt.dp1.and.aws) THEN
            lambda=lambda*1.25
C    increase lambda to weaken stochastic penalization
            goto 8999
	 END IF
	 DO k=1,dp2
	    bi(k,i)=bii(k)
	 END DO
      END DO      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in bivariate local polynomial aws (gridded) 
C
C   p > 0   only  !!!! 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vpawsbi(n1,n2,d,dp1,dp2,y,fix,theta,bi,bi0,ai,
     1    lam,h,kern,dmat,thij,psix,si,si0,siy,ind,wght,spmax,vw)
C
C     n1         number of points in first dimension
C     n2         number of points in second dimension
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     y          observed values at design points
C     theta      estimates from step (k-1)   (input)
C     bi         \sum \Psi^T Wi \Psi  from step (k-1)
C     bi0        \sum \Psi^T Wi0 \Psi  from step (k-1)
C     ai        \sum \Psi^T Wi Y     (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     dmat       working arrays  dp1 times dp1
C     thij               vector for parameter differences
C     psix                     working memory for Psij
C     si         working array  for \sum \Psi^T Wi \Psi
C     si0        working array  for \sum \Psi^T Wi0 \Psi
C     siy        working array  for \sum \Psi^T Wi Y 
C
      implicit logical (a-z)
      integer n1,n2,d,dp1,dp2,i1,i2,j1,j2,k,l,info,je1,ja1,je2,ja2,
     1        j,ih1,ih2,ind(dp1,dp1),m,kern
      logical aws,fix(n1,n2)
      real*8 bi(dp2,n1,n2),bi0(dp2,n1,n2),lkern,ai(dp1,d,n1,n2),wght,
     1       theta(dp1,d,n1,n2),siy(dp1,d),dmat(dp1,dp1),thij(dp1),
     2       spmax,y(d,n1,n2),h,psix(dp2),si(dp2),si0(dp2),sii,wijy,z,
     3       z11,z12,z22,ha2,lambda,wij,s0i,z1,z2,lam,vw(d),z0,thijl
      external lkern
C     
C     in case of dp1==1  lawsbi should be called  (p=0)
C
      ha2=h*h
      ih1=h
      aws=lam.lt.1.d20
      DO i1=1,n1
         DO i2=1,n2
	    IF (fix(i1,i2)) CYCLE
C    nothing to do, estimate in (i1,i2) is already fixed by control
            lambda=lam
C  first fill si and siy with 0's
C  fields are used to sum components of ai, bi and bi0
8999        DO j=1,dp2
               si(j)=0.d0
               si0(j)=0.d0
               IF (j.le.dp1) THEN
	          DO m=1,d
	             siy(j,m)=0.d0
		  END DO
	       END IF
            END DO
            IF (aws) THEN
               DO j=1,dp1
                  DO k=1,dp1
                     m=ind(j,k)
                     dmat(j,k)=bi(m,i1,i2)/lambda
                  END DO
               END DO
               s0i=bi0(1,i1,i2)
               sii=bi(1,i1,i2)
            END IF
C
C     Prepare for loop over j's
C
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            DO j1=ja1,je1
               z1=(i1-j1)
               z11=z1*z1
               ih2=dsqrt(ha2-z11)/wght
               ja2=max0(1,i2-ih2)
               je2=min0(n2,i2+ih2)
               DO j2=ja2,je2
                  z2=(i2-j2)*wght
                  z22=z2*z2
                  z12=z1*z2
C          first compute location part
                  wij=lkern(kern,(z11+z22)/ha2)
                  si0(1)=si0(1)+wij
                  si0(2)=si0(2)-z1*wij
                  si0(3)=si0(3)-z2*wij
                  si0(4)=si0(4)+z11*wij
                  si0(5)=si0(5)+z12*wij
                  si0(6)=si0(6)+z22*wij
                  IF (dp1.gt.3) THEN
                     si0(7)=si0(7)-z11*z1*wij
                     si0(8)=si0(8)-z11*z2*wij
                     si0(9)=si0(9)-z1*z22*wij
                     si0(10)=si0(10)-z2*z22*wij
                     si0(11)=si0(11)+z11*z11*wij
                     si0(12)=si0(12)+z11*z12*wij
                     si0(13)=si0(13)+z11*z22*wij
                     si0(14)=si0(14)+z12*z22*wij
                     si0(15)=si0(15)+z22*z22*wij
		  END IF
C          this is the location penalty, now fill si0
                  IF (aws) THEN 
C          now fill psix 
C           now translate thetaj into model centered in xi
C                    
                     z=0.d0
                     DO m=1,d
                        thij(1)=theta(1,m,j1,j2)
                        thij(2)=theta(2,m,j1,j2)
                        thij(3)=theta(3,m,j1,j2)
                        thij(1)=thij(1)+thij(2)*z1+thij(3)*z2
                        IF (dp1.gt.3) THEN
                           thij(4)=theta(4,m,j1,j2)
                           thij(5)=theta(5,m,j1,j2)
                           thij(6)=theta(6,m,j1,j2)
                           thij(1)=thij(1)+thij(4)*z11+thij(5)*z12+
     1                                     thij(6)*z22
                           thij(2)=thij(2)+thij(5)*z2+2.d0*thij(4)*z1
                           thij(3)=thij(3)+thij(5)*z1+2.d0*thij(6)*z2
		        END IF
C  
C           get difference of thetas
C
                        DO l=1,dp1
                           thij(l)=theta(l,m,i1,i2)-thij(l)
                        END DO
C
C           get stochastic penalty
C
                        z0=0.d0
                        DO l=1,dp1
			   thijl=thij(l)
			   z0=z0+dmat(l,l)*thijl*thijl
			   IF(l.eq.dp1) CYCLE
                           DO k=l+1,dp1
                              z0=z0+2.d0*dmat(k,l)*thijl*thij(k)
                           END DO
                        END DO
			z=z+z0*vw(m)
		     END DO
C
C           now compute weights 
C           stochastic and extension penalty can be added because kernel is exp
C
                     IF (z.gt.spmax) CYCLE
                     wij=wij*dexp(-z)
	          END IF
C           now compute contributions to bi(i),ai(i)  
                  si(1)=si(1)+wij
                  si(2)=si(2)-z1*wij
                  si(3)=si(3)-z2*wij
                  si(4)=si(4)+z11*wij
                  si(5)=si(5)+z12*wij
                  si(6)=si(6)+z22*wij
		  DO m=1,d
                     wijy=y(m,j1,j2)*wij
                     siy(1,m)=siy(1,m)+wijy
                     siy(2,m)=siy(2,m)-z1*wijy
                     siy(3,m)=siy(3,m)-z2*wijy
		  END DO
        	  IF (dp1.le.3) CYCLE
                  si(7)=si(7)-z11*z1*wij
                  si(8)=si(8)-z11*z2*wij
                  si(9)=si(9)-z1*z22*wij
                  si(10)=si(10)-z2*z22*wij
                  si(11)=si(11)+z11*z11*wij
                  si(12)=si(12)+z11*z12*wij
                  si(13)=si(13)+z11*z22*wij
                  si(14)=si(14)+z12*z22*wij
                  si(15)=si(15)+z22*z22*wij
		  DO m=1,d
                     wijy=y(m,j1,j2)*wij
                     siy(4,m)=siy(4,m)+z11*wijy
                     siy(5,m)=siy(5,m)+z12*wijy
		     siy(6,m)=siy(6,m)+z22*wijy
		  END DO
               END DO
            END DO
C        prepare matrix to test for singularity of Bi 
C        this should be changed to SVD at some point
            DO k=1,dp1
               DO j=1,dp1
                  IF (j.gt.k) THEN
                     dmat(j,k)=0.d0
                  ELSE
                     dmat(j,k)=si(ind(j,k))
                  ENDIF
               END DO
            END DO
C
C     compute choleski decomposition
C
            call invers(dmat,dp1,info)
            IF (info.gt.0) THEN
C
C          if singular relax stochastic and extension penalty 
C
               lambda=1.5*lambda
               goto 8999
            END IF
C     
C     now fill ai, bi and bi0
C
            DO j=1,dp1
	       DO m=1,d
                  ai(j,m,i1,i2)=siy(j,m)
	       END DO
            END DO
            DO j=1,dp2
               bi(j,i1,i2)=si(j)
               bi0(j,i1,i2)=si0(j)
            END DO
         END DO     
      END DO     
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in bivariate local polynomial aws (gridded) 
C
C   p > 0   only  !!!! 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vphawsbi(n1,n2,d,dp1,dp2,y,fix,si2,theta,bi,bi0,
     1    ai,lam,h,kern,dmat,thij,psix,si,si0,siy,ind,wght,spmax,vw)
C
C
C     n1         number of points in first dimension
C     n2         number of points in second dimension
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     y          observed values at design points
C     theta      estimates from step (k-1)   (input)
C     bi         \sum \Psi^T Wi \Psi  from step (k-1)
C     bi0        \sum \Psi^T Wi0 \Psi  from step (k-1)
C     ai         \sum \Psi^T Wi Y     (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     dmat       working arrays  dp1 times dp1
C     thij               vector for parameter differences
C     psix                     working memory for Psij
C     si         working array  for \sum \Psi^T Wi \Psi
C     si0        working array  for \sum \Psi^T Wi0 \Psi
C     siy        working array  for \sum \Psi^T Wi Y 
C
      implicit logical (a-z)
      integer n1,n2,d,dp1,dp2,i1,i2,j1,j2,k,l,info,je1,ja1,je2,ja2,
     1        j,ih1,ih2,ind(dp1,dp1),m,kern
      logical aws,fix(n1,n2)
      real*8 bi(dp2,n1,n2),bi0(dp2,n1,n2),lkern,skern,ai(dp1,d,n1,n2),
     1       wght,si2(n1,n2),theta(dp1,d,n1,n2),dmat(dp1,dp1),spmax,h,
     2       siy(dp1,d),thij(dp1),y(d,n1,n2),lam,psix(dp2),si(dp2),
     3       si0(dp2),sii,wijy,z,z11,z12,z22,ha2,lambda,wij,s0i,z1,z2,
     4       vw(d),z0,thijl
      external lkern,skern
C     
C     in case of dp1==1  lawsbi should be called  (p=0)
C
      ha2=h*h
      ih1=h
      aws=lam.lt.1.d20
      DO i1=1,n1
         DO i2=1,n2
	    IF (fix(i1,i2)) CYCLE
C    nothing to do, estimate in (i1,i2) is already fixed by control
            lambda=lam
C  first fill si and siy with 0's
C  fields are used to sum components of ai, bi and bi0
8999        DO j=1,dp2
               si(j)=0.d0
               si0(j)=0.d0
               IF (j.le.dp1) THEN
	          DO m=1,d
	             siy(j,m)=0.d0
		  END DO
	       END IF
            END DO
            IF (aws) THEN
               DO j=1,dp1
                  DO k=1,dp1
                     m=ind(j,k)
                     dmat(j,k)=bi(m,i1,i2)/lambda
                  END DO
               END DO
               s0i=bi0(1,i1,i2)
               sii=bi(1,i1,i2)
            END IF
C
C     Prepare for loop over j's
C
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            DO j1=ja1,je1
               z1=(i1-j1)
               z11=z1*z1
               ih2=dsqrt(ha2-z11)/wght
               ja2=max0(1,i2-ih2)
               je2=min0(n2,i2+ih2)
               DO j2=ja2,je2
                  z2=(i2-j2)*wght
                  z22=z2*z2
                  z12=z1*z2
C          first compute location part
                  wij=lkern(kern,(z11+z22)/ha2)*si2(j1,j2)
                  si0(1)=si0(1)+wij
                  si0(2)=si0(2)-z1*wij
                  si0(3)=si0(3)-z2*wij
                  si0(4)=si0(4)+z11*wij
                  si0(5)=si0(5)+z12*wij
                  si0(6)=si0(6)+z22*wij
                  IF (dp1.gt.3) THEN
                     si0(7)=si0(7)-z11*z1*wij
                     si0(8)=si0(8)-z11*z2*wij
                     si0(9)=si0(9)-z1*z22*wij
                     si0(10)=si0(10)-z2*z22*wij
                     si0(11)=si0(11)+z11*z11*wij
                     si0(12)=si0(12)+z11*z12*wij
                     si0(13)=si0(13)+z11*z22*wij
                     si0(14)=si0(14)+z12*z22*wij
                     si0(15)=si0(15)+z22*z22*wij
		  END IF
C          this is the location penalty, now fill si0
                  IF (aws) THEN 
C          now fill psix 
C           now translate thetaj into model centered in xi
C                    
                     z=0.d0
                     DO m=1,d
                        thij(1)=theta(1,m,j1,j2)
                        thij(2)=theta(2,m,j1,j2)
                        thij(3)=theta(3,m,j1,j2)
                        thij(1)=thij(1)+thij(2)*z1+thij(3)*z2
                        IF (dp1.gt.3) THEN
                           thij(4)=theta(4,m,j1,j2)
                           thij(5)=theta(5,m,j1,j2)
                           thij(6)=theta(6,m,j1,j2)
                           thij(1)=thij(1)+thij(4)*z11+thij(5)*z12+
     1                                     thij(6)*z22
                           thij(2)=thij(2)+thij(5)*z2+2.d0*thij(4)*z1
                           thij(3)=thij(3)+thij(5)*z1+2.d0*thij(6)*z2
		        END IF
C  
C           get difference of thetas
C
                        DO l=1,dp1
                           thij(l)=theta(l,m,i1,i2)-thij(l)
                        END DO
C
C           get stochastic penalty
C
                        z0=0.d0
                        DO l=1,dp1
			   thijl=thij(l)
			   z0=z0+dmat(l,l)*thijl*thijl
			   IF(l.eq.dp1) CYCLE
                           DO k=l+1,dp1
                              z0=z0+2.d0*dmat(k,l)*thijl*thij(k)
                           END DO
                        END DO
			z=z+z0*vw(m)
		     END DO
C
C           now compute weights 
C           stochastic and extension penalty can be added because kernel is exp
C
                     IF (z.gt.spmax) CYCLE
		     if(z.lt.0.d0) call dblepr("z",1,z,1)
                     wij=wij*dexp(-z)
	          END IF
C           now compute contributions to bi(i),ai(i)  
                  si(1)=si(1)+wij
                  si(2)=si(2)-z1*wij
                  si(3)=si(3)-z2*wij
                  si(4)=si(4)+z11*wij
                  si(5)=si(5)+z12*wij
                  si(6)=si(6)+z22*wij
		  DO m=1,d
                     wijy=y(m,j1,j2)*wij
                     siy(1,m)=siy(1,m)+wijy
                     siy(2,m)=siy(2,m)-z1*wijy
                     siy(3,m)=siy(3,m)-z2*wijy
		  END DO
        	  IF (dp1.le.3) CYCLE
                  si(7)=si(7)-z11*z1*wij
                  si(8)=si(8)-z11*z2*wij
                  si(9)=si(9)-z1*z22*wij
                  si(10)=si(10)-z2*z22*wij
                  si(11)=si(11)+z11*z11*wij
                  si(12)=si(12)+z11*z12*wij
                  si(13)=si(13)+z11*z22*wij
                  si(14)=si(14)+z12*z22*wij
                  si(15)=si(15)+z22*z22*wij
		  DO m=1,d
                     wijy=y(m,j1,j2)*wij
                     siy(4,m)=siy(4,m)+z11*wijy
                     siy(5,m)=siy(5,m)+z12*wijy
		     siy(6,m)=siy(6,m)+z22*wijy
		  END DO
               END DO
            END DO
C        prepare matrix to test for singularity of Bi 
C        this should be changed to SVD at some point
            DO k=1,dp1
               DO j=1,dp1
                  IF (j.gt.k) THEN
                     dmat(j,k)=0.d0
                  ELSE
                     dmat(j,k)=si(ind(j,k))
                  ENDIF
               END DO
            END DO
C
C     compute choleski decomposition
C
            call invers(dmat,dp1,info)
            IF (info.gt.0) THEN
C
C          if singular relax stochastic and extension penalty 
C
               lambda=1.5*lambda
               goto 8999
            END IF
C     
C     now fill ai, bi and bi0
C
            DO j=1,dp1
	       DO m=1,d
                  ai(j,m,i1,i2)=siy(j,m)
	       END DO
            END DO
            DO j=1,dp2
               bi(j,i1,i2)=si(j)
               bi0(j,i1,i2)=si0(j)
            END DO
         END DO     
      END DO     
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (univariate polynomial aws)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vmpawsun(n,d,dp1,dp2,ai,bi,theta,dmat)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     dpm        number of components in di  (dp1+1)*dp1/2
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi
C     di         inverse of bi     
C     di0         inverse of bi0     
C     theta      new parameter estimate
C     dmat       working array
C
C      implicit logical(a-z)
      integer n,d,dp1,dp2
      real*8 ai(dp1,d,n),bi(dp2,n),theta(dp1,d,n),dmat(dp1,dp1)
      integer i,j,k,info,m
      real*8 z
      DO i=1,n
         DO j=1,dp1
            DO k=1,dp1
               IF (j.gt.k) then 
                  dmat(j,k)=0.0d0
               ELSE
                  dmat(j,k)=bi(j+k-1,i)
               END IF
            END DO
	 END DO
         call invers(dmat,dp1,info)
         IF (info.ne.0) CYCLE 
C      now dmat contains inverse of B_i 
C      now calculate theta as B_i^{-1} A_i
         DO m=1,d
	    DO j=1,dp1
               z=0.0d0
               DO k=1,dp1
                  z=z+dmat(j,k)*ai(k,m,i)  
               END DO
               theta(j,m,i)=z
            END DO
	 END DO
C     just keep the old estimate if info > 0
      END DO
      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (bivariate case)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vmpawsbi(n,d,dp1,dp2,ai,bi,theta,dmat,ind)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi    
C     theta      new parameter estimate
C     dmat       working arrays
C
C      implicit logical (a-z)
      integer n,d,dp1,dp2
      real*8 ai(dp1,d,n),bi(dp2,n),theta(dp1,d,n),dmat(dp1,dp1)
      integer i,j,k,m,info,ind(dp1,dp1)
      real*8 z
      DO i=1,n
         DO k=1,dp1
            DO j=1,dp1
               IF (j.gt.k) then
                  dmat(j,k)=0.d0
               ELSE
                  dmat(j,k)=bi(ind(j,k),i)
               END IF
            END DO
	 END DO
         call invers(dmat,dp1,info)
         IF (info.gt.0) CYCLE  
C     just keep the old estimate
C     now dmat contains inverse of B_i 
C     now calculate theta as B_i^{-1} A_i
         DO m=1,d
            DO j=1,dp1
               z=0.0d0
               DO k=1,dp1
                  z=z+dmat(j,k)*ai(k,m,i)  
               END DO
               theta(j,m,i)=z
            END DO
         END DO
      END DO
      RETURN
      END      
