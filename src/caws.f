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
C  FORTRAN 77 code needed in R functions aws, vaws, 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    location penalty for multivariate non-gridded aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function lmkern(kern,dx,xi,xj,h2)
      implicit logical (a-z)
      external lkern
      integer kern,dx,i
      real*8 xi(dx),xj(dx),h2,z,zd,lkern
      z=0.d0
      do 1 i=1,dx
         zd=xi(i)-xj(i)
         z=z+zd*zd
         if(z.gt.h2) goto 2
1     continue
      lmkern=lkern(kern,z/h2)
      goto 999
2     lmkern=0.d0
999   return
      end
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
      real*8 function kldist(model,thi,thj,bi0)
      implicit logical (a-z)
      integer model
      real*8 thi,thj,z,thij,bi0,eta,tthi
      IF (model.eq.1) THEN
C        Gaussian
         z=thi-thj
         kldist=z*z
      ELSE IF (model.eq.2) THEN
C        Bernoulli
         kldist=0.d0
         eta=0.5d0/bi0
         thij=(1.d0-eta)*thj+eta*thi
         tthi=(1.d0-thi)
         IF (thi.gt.1.d-10) kldist=kldist+thi*log(thi/thij)
         IF (tthi.gt.1.d-10) kldist=kldist+tthi*log(tthi/(1.d0-thij))
      ELSE IF (model.eq.3) THEN
C        Poisson
         kldist=0.d0
         eta=0.5d0/bi0
         thij=(1.d0-eta)*thj+eta*thi
         IF (thi.gt.1.d-10) kldist=thi*log(thi/thij)-thi+thij
      ELSE IF (model.eq.4) THEN
C        Exponential
         kldist=thi/thj-1.d0-log(thi/thj)
      ELSE IF (model.eq.5) THEN
C        Variance
         kldist=thi/thj-1.d0-log(thi/thj)
      ELSE
C        use Gaussian
         z=thi-thj
         kldist=z*z
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute Location Kernel (Compact support only, based on x^2
C                                   ignores scaling)
C
C          Kern=1     Uniform
C          Kern=2     Epanechnicov
C          Kern=3     Biweight
C          Kern=4     Triweight
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function lkern(kern,xsq)
      implicit logical (a-z)
      integer kern
      real*8 xsq,z
      IF (xsq.ge.1) THEN
         lkern=0.d0
      ELSE IF (kern.eq.1) THEN
         IF(xsq.le.0.5d0) THEN
            lkern=1.d0
         ELSE
            lkern=2.d0*(1.d0-xsq)
         END IF
      ELSE IF (kern.eq.2) THEN
         lkern=1.d0-xsq
      ELSE IF (kern.eq.3) THEN
         z=1.d0-xsq
         lkern=z*z
      ELSE IF (kern.eq.4) THEN
         z=1.d0-xsq
         lkern=z*z*z
      ELSE IF (kern.eq.5) THEN
         lkern=exp(-xsq*8.d0)
      ELSE
C        use Epanechnikov
         lkern=1.d0-xsq
      ENDIF
      RETURN 
      END   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        Compute truncated Exponential Kernel 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function skern(x,xmin,xmax)
      implicit logical (a-z)
      real*8 x,xmin,xmax,spf
      spf=xmax/(xmax-xmin)
      IF (x.le.xmin) THEN
         skern=1.d0
      ELSE IF (x.gt.xmax) THEN
         skern=0.d0
      ELSE
         skern=exp(-spf*(x-xmin))
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine caws(y,fix,n1,n2,n3,hakt,hhom,lambda,theta,bi,bi2,
     1                bi0,ai,model,kern,spmin,lwght,wght)
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
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)

      external kldist,lkern
      real*8 kldist,lkern
      integer n1,n2,n3,model,kern
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,wght(2),
     1       bi2(1),hakt,lwght(1),spmin,spf,hhom(1)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12
      real*8 thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,hakt2,bii0,
     1       hmax2,hhomi,hhommax,w1,w2
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih1=hakt
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=hakt/w2
      ih2=hakt/w1
      ih1=hakt
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      hmax2=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=sqrt(hakt2-z3)/w1
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
               if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ai,bi,bi0,bi2,hhom,n1,n2,n3,hakt2,hmax2,theta,
C$OMP& ih3,lwght,wght,y,fix)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& model,spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,hhomi,hhommax,thetai,bii,bii0,swj,swj2,
C$OMP& swj0,swjy,sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1         
         hhomi=hhom(iind)
         hhomi=hhomi*hhomi
         hhommax=hmax2
         IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         thetai=theta(iind)
         bii=bi(iind)/lambda
C   scaling of sij outside the loop
         bii0=bi0(iind)
         swj=0.d0
         swj2=0.d0
         swj0=0.d0
         swjy=0.d0
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=sqrt(hakt2-z3)/w1
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  wj=lwght(jw1+clw1+1+jwind2)
                  swj0=swj0+wj
                  z1=jw1
                  z1=z2+z1*z1
                  IF (aws.and.z1.ge.hhomi) THEN
                     sij=bii*kldist(model,thetai,theta(jind),bii0)
                     IF (sij.gt.1.d0) THEN
                        hhommax=min(hhommax,z1)
                        CYCLE
                     END IF
                     IF (sij.gt.spmin) THEN
                        wj=wj*(1.d0-spf*(sij-spmin))
                        hhommax=min(hhommax,z1)
                     END IF
                  END IF
                  swj=swj+wj
                  swj2=swj2+wj*wj
                  swjy=swjy+wj*y(jind)
               END DO
            END DO
         END DO
         ai(iind)=swjy
         bi(iind)=swj
         bi2(iind)=swj2
         bi0(iind)=swj0
         hhom(iind)=sqrt(hhommax)
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2,hhom)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C   used in awstestprop only 
C   bi0 contains sum of weights (without invers variances) !!! 
C   no lambda, fix, spmin, hhom, theta
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine caws1(y,n1,n2,n3,hakt,bi,bi2,
     1                bi0,ai,model,kern,lwght,wght)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)

      external kldist,lkern
      real*8 kldist,lkern
      integer n1,n2,n3,model,kern
      real*8 y(1),bi(1),bi0(1),ai(1),wght(2),
     1       bi2(1),hakt,lwght(1)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12
      real*8 swj,swj2,swj0,swjy,z1,z2,z3,wj,hakt2,hmax2,w1,w2
      hakt2=hakt*hakt
      w1=wght(1)
      w2=wght(2)
      ih1=hakt
C
C   first calculate location weights
C
      ih3=hakt/w2
      ih2=hakt/w1
      ih1=hakt
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      hmax2=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=sqrt(hakt2-z3)/w1
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
               if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ai,bi,bi0,bi2,n1,n2,n3,hakt2,hmax2
C$OMP& ,lwght,wght,y)
C$OMP& FIRSTPRIVATE(ih1,ih2,model,dlw1,clw1,dlw2,clw2,dlw3,clw3,w1,w2,
C$OMP& n12,dlw12)
C$OMP& PRIVATE(i1,i2,i3,iind,swj,swj2,swj0,swjy,wj
C$OMP& ,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,j1,jw1,jind,z1)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
C    nothing to do, final estimate is already fixed by control 
         swj=0.d0
         swj2=0.d0
         swj0=0.d0
         swjy=0.d0
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=sqrt(hakt2-z3)/w1
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  wj=lwght(jw1+clw1+1+jwind2)
                  swj0=swj0+wj
                  z1=jw1
                  z1=z2+z1*z1
                  swj=swj+wj
                  swj2=swj2+wj*wj
                  swjy=swjy+wj*y(jind)
               END DO
            END DO
         END DO
         ai(iind)=swjy
         bi(iind)=swj
         bi2(iind)=swj2
         bi0(iind)=swj0
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chaws(y,fix,si2,n1,n2,n3,hakt,lambda,theta,bi,bi2,
     1           bi0,vred,ai,model,kern,spmin,lwght,wght)
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
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      external kldist,lkern
      real*8 kldist,lkern
      integer n1,n2,n3,model,kern
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,wght(2),
     1       bi2(1),hakt,lwght(1),si2(1),vred(1),spmin
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12
      real*8 thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,hakt2,bii0,
     1        sv1,sv2,spf,w1,w2,wjsi2
      w1=wght(1)
      w2=wght(2)
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih1=hakt
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=hakt/w2
      ih2=hakt/w1
      ih1=hakt
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=sqrt(hakt2-z3)/w1
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ai,bi,bi0,bi2,si2,vred,n1,n2,n3,hakt2,hakt,theta
C$OMP& ,lwght,wght,y,fix)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,dlw12,n12
C$OMP& ,model,spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,w1,w2)
C$OMP& PRIVATE(iind,thetai,bii,bii0,swj
C$OMP& ,swj2,swj0,swjy,sij,sv1,sv2,i1,i2,i3,wj
C$OMP& ,j3,jw3,jind3,z3,jwind3
C$OMP& ,j2,jw2,jind2,z2,jwind2
C$OMP& ,j1,jw1,jind,z1,wjsi2)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
         IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         thetai=theta(iind)
         bii=bi(iind)/lambda
C   scaling of sij outside the loop
         bii0=bi0(iind)
         swj=0.d0
         swj2=0.d0
         swj0=0.d0
         swjy=0.d0
         sv1=0.d0
         sv2=0.d0
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=sqrt(hakt2-z3)/w1
            jwind3=(jw3+clw3)*dlw12
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               jwind2=jwind3+(jw2+clw2)*dlw1
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  wj=lwght(jw1+clw1+1+jwind2)
                  swj0=swj0+wj*si2(jind)      
                  IF (aws) THEN
                     sij=bii*kldist(model,thetai,theta(jind),bii0)
                     IF (sij.gt.1.d0) CYCLE
                     IF (sij.gt.spmin) wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
                  sv1=sv1+wj
                  sv2=sv2+wj*wj
                  wjsi2=wj*si2(jind)
                  swj=swj+wjsi2
                  swj2=swj2+wj*wjsi2
                  swjy=swjy+wjsi2*y(jind)
               END DO
            END DO
         END DO
         ai(iind)=swjy
         bi(iind)=swj
         bi2(iind)=swj2
         bi0(iind)=swj0
         vred(iind)=sv2/sv1/sv1
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2,vred)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant non-adaptive three-variate aws (gridded)
C   used in awstestprop only 
C   bi0 contains sum of weights (without invers variances) !!! 
C   no lambda, fix, spmin
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine chaws1(y,si2,n1,n2,n3,hakt,bi,bi2,
     1           bi0,vred,ai,model,kern,lwght,wght)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      external lkern
      real*8 lkern
      integer n1,n2,n3,model,kern
      real*8 y(1),bi(1),bi0(1),ai(1),wght(2),
     1       bi2(1),hakt,lwght(1),si2(1),vred(1)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12
      real*8 swj,swj2,swjy,z1,z2,z3,wj,hakt2,sv1,sv2,w1,w2
      hakt2=hakt*hakt
      w1=wght(1)
      w2=wght(2)
      ih1=hakt
C
C   first calculate location weights
C
      ih3=hakt/w2
      ih2=hakt/w1
      ih1=hakt
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=sqrt(hakt2-z3)/w1
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
             END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ai,bi,bi0,bi2,si2,vred,n1,n2,n3,hakt2,hakt
C$OMP& ,lwght,wght,y)
C$OMP& FIRSTPRIVATE(ih1,ih2,model,dlw1,clw1,dlw2,clw2,dlw3,clw3,w1,w2,
C$OMP& dlw12,n12)
C$OMP& PRIVATE(iind,swj,swj2,swjy,sv1,sv2,i1,i2,i3,wj
C$OMP& ,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2
C$OMP& ,j1,jw1,jind,z1)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
         swj=0.d0
         swj2=0.d0
         swjy=0.d0
         sv1=0.d0
         sv2=0.d0
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=sqrt(hakt2-z3)/w1
            jwind3=(jw3+clw3)*dlw12
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               jwind2=jwind3+(jw2+clw2)*dlw1
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  wj=lwght(jw1+clw1+1+jwind2)
                  sv1=sv1+wj
                  sv2=sv2+wj*wj
                  swj=swj+wj*si2(jind)
                  swj2=swj2+wj*wj*si2(jind)
                  swjy=swjy+wj*si2(jind)*y(jind)
               END DO
            END DO
         END DO
         ai(iind)=swjy
         bi(iind)=swj
         bi2(iind)=swj2
         bi0(iind)=sv1
         vred(iind)=sv2/sv1/sv1
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2,vred)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded) with variance - mean model
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cgaws(y,fix,mask,si2,n1,n2,n3,hakt,hhom,lambda,
     1        theta,bi,bi2,bi0,gi,vred,ai,kern,spmin,lwght,wght)
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
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit logical (a-z)
      external lkern
      real*8 lkern
      integer n1,n2,n3,kern
      logical aws,fix(1),mask(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,wght(2),
     1       bi2(1),hakt,lwght(1),si2(1),vred(1),spmin,hhom(1),gi(1)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12
      real*8 thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,hakt2,bii0,
     1        sv1,sv2,spf,z,hhomi,hhommax,hmax2,w1,w2
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      w1=wght(1)
      w2=wght(2)
      ih1=hakt
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=hakt/w2
      ih2=hakt/w1
      ih1=hakt
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      hmax2=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=sqrt(hakt2-z3)/w1
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
               if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ai,bi,bi0,bi2,si2,hhom,n1,n2,n3,hakt2,hmax2,theta
C$OMP& ,lwght,wght,y,fix,mask,vred,gi)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,dlw12
C$OMP& ,model,spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,w1,w2)
C$OMP& PRIVATE(iind,hhomi,hhommax,thetai,bii,bii0,swj
C$OMP& ,swj2,swj0,swjy,sij,sv1,sv2,i1,i2,i3,wj
C$OMP& ,j3,jw3,jind3,z3,jwind3
C$OMP& ,j2,jw2,jind2,z2,jwind2
C$OMP& ,j1,jw1,jind,z1,z)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
         hhomi=hhom(iind)
         hhomi=hhomi*hhomi
         hhommax=hmax2
         IF (fix(iind).or..not.mask(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         thetai=theta(iind)
         bii=bi(iind)/lambda
C   scaling of sij outside the loop
         bii0=bi0(iind)
         swj=0.d0
         swj2=0.d0
         swj0=0.d0
         swjy=0.d0
         sv1=0.d0
         sv2=0.d0
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=sqrt(hakt2-z3)/w1
            jwind3=(jw3+clw3)*dlw12
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               jwind2=jwind3+(jw2+clw2)*dlw1
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  if(.not.mask(jind)) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  swj0=swj0+wj*si2(jind)
                  z1=-jw1
                  z1=z2+z1*z1
                  IF (aws.and.z1.ge.hhomi) THEN
C
C      gaussian case only
C
                     z=(thetai-theta(jind))
                     sij=bii*z*z
                     IF (sij.gt.1.d0) THEN
                        hhommax=min(hhommax,z1)
                        CYCLE
                     END IF
                     IF (sij.gt.spmin) THEN
                        wj=wj*(1.d0-spf*(sij-spmin))
                        hhommax=min(hhommax,z1)
                     END IF
                  END IF
                  sv1=sv1+wj
                  sv2=sv2+wj*wj
                  swj=swj+wj*si2(jind)
                  swj2=swj2+wj*wj*si2(jind)
                  swjy=swjy+wj*si2(jind)*y(jind)
               END DO
            END DO
         END DO
         ai(iind)=swjy
         bi(iind)=swj
         bi2(iind)=swj2
         bi0(iind)=swj0
         hhom(iind)=sqrt(hhommax)
         gi(iind)=sv1
         vred(iind)=sv2/sv1/sv1
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(ai,bi,bi0,bi2,hhom,gi,vred)
      RETURN
      END
