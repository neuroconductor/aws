CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded) with variance - mean model
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine segment(y,fix,level,delta,si2,n1,n2,n3,hakt,
     1        lambda,theta,bi,bi2,bi0,gi,vred,thetan,kern,spmin,lwght,
     2        wght,pvalue,segm,segmn,beta,thresh,ext,fov,varest)
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
      external lkern,fpchisq
      real*8 lkern,fpchisq
      integer n1,n2,n3,kern,segm(1),segmn(1)
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),thetan(1),lambda,wght(2),
     1       bi2(1),hakt,lwght(1),si2(1),vred(1),spmin,gi(1),
     2       level,delta,pvalue(1),beta,thresh,ext,varest(1),fov
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,n12,
     2        segmi
      real*8 bii,sij,swj,swj2,swj0,swjy,z1,z2,z3,wj,hakt2,bii0,pvi,
     1        sv1,sv2,spf,z,a,b,thi,s2i,si,ti,cofh,extthr,wght1,wght2,
     2        si2j
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih1=hakt
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih3=hakt/wght(2)
      ih2=hakt/wght(1)
      ih1=hakt
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      n12=n1*n2
      wght1=wght(1)
      wght2=wght(2)
      clw1=ih1+1
      clw2=ih2+1
      clw3=ih3+1
      dlw1=ih1+clw1
      dlw2=ih2+clw2
      dlw3=ih3+clw3
      z2=0.d0
      z3=0.d0
      swj0=0.d0
      extthr=thresh+ext
      DO j3=1,dlw3
         if(n3.gt.1) THEN
            z3=(clw3-j3)*wght2
            z3=z3*z3
            ih2=sqrt(hakt2-z3)/wght1
            jind3=(j3-1)*dlw1*dlw2
         ELSE
            jind3=0
         END IF
         DO j2=clw2-ih2,clw2+ih2
            if(n2.gt.1) THEN
               z2=(clw2-j2)*wght1
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               jind2=jind3+(j2-1)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=clw1-ih1,clw1+ih1
C  first stochastic term
               jind=j1+jind2
               z1=clw1-j1
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
               if(fix(iind)) CYCLE
               thi = theta(iind)
               s2i = si2(iind)
            cofh = sqrt(beta*log(varest(iind)*s2i*fov))
           if(max(a-thi,thi-b)/sqrt(varest(iind))-cofh.gt.extthr) THEN
                  fix(iind)=.TRUE.
                  if(segm(iind).eq.0) segm(iind)=sign(1.d0,thi-level)
C we need to assign a value to segment before we can fix the decision
               ELSE
                  fix(iind)=.FALSE.
                  ti=max(0.d0,max(a-thi,thi-b))
                  pvalue(iind)=fpchisq(ti*s2i,1.d0,1,0)
               END IF
            END DO
         END DO
      END DO
      END IF
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thetan,bi,bi0,bi2,si2,n1,n2,n3,theta,kern,hakt,
C$OMP& lwght,y,fix,vred,gi,segm,segmn,pvalue,varest)
C$OMP& FIRSTPRIVATE(lambda,aws,beta,fov,a,b,hakt2,n12,wght1,wght2,
C$OMP& spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,thresh,ih1,ih2)
C$OMP& PRIVATE(iind,bii,bii0,swj,thi,cofh,segmi,pvi,
C$OMP& swj2,swj0,swjy,si,sij,sv1,sv2,i1,i2,i3,wj,j3,jw3,jind3,z3,
C$OMP& jwind3,j2,jw2,jind2,z2,jwind2,j1,jw1,jind,z1,z,si2j)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n1/n2+1         
         segmi=segm(iind)
         pvi=pvalue(iind)
         IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         thi=theta(iind)
         bii=bi(iind)/lambda
C   scaling of sij outside the loop
         bii0=bi0(iind)
         swj=0.d0
         swj2=0.d0
         swj0=0.d0
         swjy=0.d0
         sv1=0.d0
         sv2=0.d0
         DO jw3=1,dlw3
            j3=jw3-clw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jind3=(j3-1)*n12
            z3=(clw3-jw3)*wght2
            z3=z3*z3
            if(n2.gt.1) ih2=sqrt(hakt2-z3)/wght1
            jwind3=(jw3-1)*dlw1*dlw2
            DO jw2=clw2-ih2,clw2+ih2
               j2=jw2-clw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jind2=(j2-1)*n1+jind3
               z2=(clw2-jw2)*wght1
               z2=z3+z2*z2
               ih1=sqrt(hakt2-z2)
               jwind2=jwind3+(jw2-1)*dlw1
               DO jw1=clw1-ih1,clw1+ih1
C  first stochastic term
                  j1=jw1-clw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  wj=lwght(jw1+jwind2)
                  swj0=swj0+wj*si2(jind)
                  IF (aws) THEN
C
C      gaussian case only
C
                     z=(thi-theta(jind))
                     sij=bii*z*z
                     IF(segmi*segm(jind).gt.0)
     1                  sij=max(pvi,pvalue(jind))*sij
                     IF (sij.gt.1.d0) CYCLE
                     IF (sij.gt.spmin) wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
                  sv1=sv1+wj
                  sv2=sv2+wj*wj
                  si2j=si2(jind)*wj
                  swj=swj+si2j
                  swj2=swj2+wj*si2j
                  swjy=swjy+si2j*y(jind)
               END DO
            END DO
         END DO
         thetan(iind)=swjy/swj
         bi(iind)=swj
         bi2(iind)=swj2
         bi0(iind)=swj0
         si=swj2/swj/swj
         varest(iind)=si
         cofh = sqrt(beta*log(si*si2(iind)*fov))
C    both are equivalent for  homogeneous si2
         si=sqrt(si)
         IF((thi-a)/si+cofh.lt.-thresh) THEN
            segmn(iind)=-1
         ELSE IF ((thi-b)/si-cofh.gt.thresh) THEN
            segmn(iind)=1
         ELSE
            segmn(iind)=0
         END IF               
         gi(iind)=sv1
         vred(iind)=sv2/sv1/sv1
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thetan,bi,bi0,bi2,gi,vred,varest,segmn)
      RETURN
      END
