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
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine corcaws(y,fix,n1,n2,n3,hakt,acor,cw1,cw2,cw3,wind,w,
     1                hhom,lambda,theta,bi,bi2,bi0,ai,model,kern,
     2                spmin,lwght,wght)
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
      implicit none

      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,n3,model,kern,cw1,cw2,cw3,wind(3,*)
      logical aws,fix(*)
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,wght(2),
     1       bi2(*),hakt,lwght(*),spmin,spf,hhom(*),w(*),
     2       acor(cw1,cw2,cw3)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,nw,k,l
      double precision thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,
     1       hakt2,hmax2,hhomi,hhommax,w1,w2,snw
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
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
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
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
      DO iind=1,n1*n2*n3
         nw=0
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
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
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
                     sij=bii*kldist(model,thetai,theta(jind))
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
                  nw=nw+1
                  wind(1,nw)=jw1
                  wind(2,nw)=jw2
                  wind(3,nw)=jw3
                  w(nw)=wj
               END DO
            END DO
         END DO
C now get scale factor for statistical penalty
         snw=swj2
         DO k=1,nw-1
            DO l=k+1,nw
               jw1=abs(wind(1,k)-wind(1,l))+1
               if(jw1.gt.cw1) CYCLE
               jw2=abs(wind(2,k)-wind(2,l))+1
               if(jw2.gt.cw2) CYCLE
               jw3=abs(wind(3,k)-wind(3,l))+1
               if(jw3.gt.cw3) CYCLE
               snw=snw+2d0*acor(jw1,jw2,jw3)*w(k)*w(l)
            END DO
         END DO
         ai(iind)=swjy
         bi(iind)=swj
C         bi2(iind)=swj*swj2/snw
         bi2(iind)=swj*swj/snw
C this takes care of correlation  penalty factor returned in bi2
         bi0(iind)=swj0
         hhom(iind)=sqrt(hhommax)
      END DO
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine corchaws(y,fix,si2,n1,n2,n3,hakt,lambda,theta,bi,bi2,
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
      implicit none
      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,n3,model,kern
      logical aws,fix(*)
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,wght(2),
     1       bi2(*),hakt,lwght(*),si2(*),vred(*),spmin
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12
      double precision thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,
     1        hakt2,sv1,sv2,spf,w1,w2,wjsi2
      w1=wght(1)
      w2=wght(2)
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
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
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
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
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,dlw12,n12,
C$OMP& model,spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,w1,w2)
C$OMP& PRIVATE(iind,thetai,bii,swj
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
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jwind3=(jw3+clw3)*dlw12
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jwind2=jwind3+(jw2+clw2)*dlw1
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  wj=lwght(jw1+clw1+1+jwind2)
                  swj0=swj0+wj*si2(jind)
                  IF (aws) THEN
                     sij=bii*kldist(model,thetai,theta(jind))
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
