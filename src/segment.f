CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded) with variance - mean model
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine segment(y,pos,fix,level,delta,si2,n1,n2,n3,hakt,
     1        lambda,theta,bi,bi2,bi0,gi,gi2,thetan,kern,spmin,lwght,
     2        wght,segm,segmn,beta,thresh,ext,fov,varest)
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
      external lkern,fpchisq
      double precision lkern,fpchisq
      integer n1,n2,n3,kern,segm(*),segmn(*),pos(*),fix(*)
      logical aws
      double precision y(*),theta(*),bi(*),bi0(*),thetan(*),lambda,
     1       wght(2),bi2(*),hakt,lwght(*),si2(*),gi2(*),spmin,gi(*),
     2       level,delta,beta,thresh,ext,varest(*),fov
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,n12,
     2        segmi,iindp,jindp
      double precision bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,hakt2,
     1        bii0,sv1,sv2,spf,z,a,b,thi,s2i,si,ti,cofh,extthr,
     2        wght1,wght2,si2j
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      wght1=wght(1)
      wght2=wght(2)
      ih3=FLOOR(hakt/wght2)
      ih2=FLOOR(hakt/wght1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      n12=n1*n2
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      z2=0.d0
      z3=0.d0
      swj0=0.d0
      extthr=thresh+ext
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*wght2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/wght1)
            jind3=(j3+clw3)*dlw1*dlw2
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*wght1
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
               wj=lkern(kern,(z1*z1+z2)/hakt2)
               swj0=swj0+wj
               lwght(jind)=wj
            END DO
         END DO
      END DO
      a = level-delta
      b = level+delta
      call rchkusr()
      IF(hakt.gt.1.25) THEN
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               iind=i1+(i2-1)*n1+(i3-1)*n12
               iindp=pos(iind)
               if(iindp.eq.0) CYCLE
               thi = theta(iindp)
               s2i = si2(iindp)
            cofh = sqrt(beta*log(varest(iindp)*s2i*fov))
           if(max(a-thi,thi-b)/sqrt(varest(iindp))-cofh.gt.extthr) THEN
      if(segm(iindp).eq.0) segm(iindp)=FLOOR(sign(1.d0,thi-level))
               ELSE
                  ti=max(0.d0,max(a-thi,thi-b))
               END IF
            END DO
         END DO
      END DO
      END IF
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thetan,bi,bi0,bi2,si2,n1,n2,n3,theta,kern,hakt,
C$OMP& lwght,y,gi2,gi,segm,segmn,varest,pos,fix)
C$OMP& FIRSTPRIVATE(lambda,aws,beta,fov,a,b,hakt2,n12,wght1,wght2,
C$OMP& spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,thresh,ih1,ih2)
C$OMP& PRIVATE(iind,bii,bii0,swj,thi,cofh,segmi,iindp,jindp,
C$OMP& swj2,swj0,swjy,si,sij,sv1,sv2,i1,i2,i3,wj,j3,jw3,jind3,z3,
C$OMP& jwind3,j2,jw2,jind2,z2,jwind2,j1,jw1,jind,z1,z,si2j)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         iindp=pos(iind)
         if(iindp.eq.0) CYCLE
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1
         segmi=segm(iindp)
         thi=theta(iindp)
         bii=bi(iindp)/lambda
C   scaling of sij outside the loop
         bii0=bi0(iindp)
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
            z3=jw3*wght2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/wght1)
            jwind3=(jw3+clw3)*dlw1*dlw2
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1+jind3
               z2=jw2*wght1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jwind2=jwind3+(jw2+clw2)*dlw1
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  jindp=pos(jind)
                  if(jindp.eq.0) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  swj0=swj0+wj*si2(jindp)
                  IF (aws) THEN
C
C      gaussian case only
C
                     z=(thi-theta(jindp))
                     sij=bii*z*z
                     IF (sij.gt.1.d0) CYCLE
                     IF (sij.gt.spmin) wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
                  sv1=sv1+wj
                  sv2=sv2+wj*wj
                  si2j=si2(jindp)*wj
                  swj=swj+si2j
                  swj2=swj2+wj*si2j
                  swjy=swjy+si2j*y(jindp)
               END DO
            END DO
         END DO
         thetan(iindp)=swjy/swj
         bi(iindp)=swj
         bi2(iindp)=swj2
         bi0(iindp)=swj0
         si=swj2/swj/swj
         varest(iindp)=si
         cofh = sqrt(beta*log(si*si2(iindp)*fov))
C    both are equivalent for  homogeneous si2
         si=sqrt(si)
         If(fix(iind).eq.0) THEN
         IF((thi-a)/si+cofh.lt.-thresh) THEN
            segmn(iindp)=-1
         ELSE IF ((thi-b)/si-cofh.gt.thresh) THEN
            segmn(iindp)=1
         ELSE
            segmn(iindp)=0
         END IF
      ELSE
         IF(segmi.lt.0) thetan(iindp)=min(thetan(iindp),a)
         IF(segmi.gt.0) thetan(iindp)=max(thetan(iindp),b)
      END IF
         gi(iindp)=sv1
         gi2(iindp)=sv2
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thetan,bi,bi0,bi2,gi,gi2,varest,segmn)
      RETURN
      END
