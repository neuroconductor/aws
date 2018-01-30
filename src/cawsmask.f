CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawsmask(y,mask,ni,fix,n1,n2,hakt,lambda,theta,
     1         bi,bi2,bi0,ai,model,kern,spmin,lwght,wght)
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
C   wght     scaling factor for second dimension (larger values shrink)
C
      implicit none
      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,model,kern,ni(*)
      logical aws,fix(*),mask(*)
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,wght,
     1       bi2(*),hakt,lwght(*),spmin,spf
      integer ih1,ih2,i1,i2,j1,j2,jw1,jw2,jwind2,
     1        iind,jind,jind2,clw1,clw2,dlw1,dlw2
      double precision thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,wj,hakt2
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih2=FLOOR(hakt/wght)
      ih1=FLOOR(hakt)
      if(n2.eq.1) ih2=0
      clw1=ih1+1
      clw2=ih2+1
      dlw1=ih1+clw1
      dlw2=ih2+clw2
      z2=0.d0
      DO j2=1,dlw2
         if(n2.gt.1) THEN
            z2=(clw2-j2)*wght
            z2=z2*z2
            ih1=FLOOR(sqrt(hakt2-z2))
            jind2=(j2-1)*dlw1
         ELSE
            jind2=0
         END IF
         DO j1=clw1-ih1,clw1+ih1
C  first stochastic term
            jind=j1+jind2
            z1=clw1-j1
            lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
         END DO
      END DO
      call rchkusr()
      DO i2=1,n2
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ai,bi,bi0,bi2,n1,n2,hakt2,theta
C$OMP& ,lwght,wght,y,fix,mask,ni)
C$OMP& FIRSTPRIVATE(ih1,i2,lambda,aws
C$OMP& ,model,dlw1,clw1,dlw2,clw2)
C$OMP& PRIVATE(iind,thetai,bii,swj
C$OMP& ,swj2,swj0,swjy,sij,i1,wj
C$OMP& ,j2,jw2,jind2,z2,jwind2
C$OMP& ,j1,jw1,jind)
C$OMP DO SCHEDULE(GUIDED)
         DO i1=1,n1
            iind=i1+(i2-1)*n1
            if(.not.mask(iind)) CYCLE
            IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control
            thetai=theta(iind)
            bii=bi(iind)/lambda
C   scaling of sij outside the loop
            swj=0.d0
            swj2=0.d0
            swj0=0.d0
            swjy=0.d0
            DO jw2=1,dlw2
               j2=jw2-clw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=(jw2-1)*dlw1
               jind2=(j2-1)*n1
               z2=(clw2-jw2)*wght
               z2=z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=clw1-ih1,clw1+ih1
C  first stochastic term
                  j1=jw1-clw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  if(ni(jind).eq.0) CYCLE
                  wj=lwght(jw1+jwind2)*ni(jind)
                  swj0=swj0+wj
                  IF (aws) THEN
                    sij=bii*kldist(model,thetai,theta(jind))
                     IF (sij.gt.1.d0) CYCLE
                     wj=wj*(1.d0-sij)
C   if sij <= spmin  this just keeps the location penalty
C    spmin = 0 corresponds to old choice of K_s
C   new kernel is flat in [0,spmin] and then decays exponentially
                  END IF
                  swj=swj+wj
                  swj2=swj2+wj*wj
                  swjy=swjy+wj*y(jind)
               END DO
            END DO
            ai(iind)=swjy
            bi(iind)=swj
            bi2(iind)=swj2
            bi0(iind)=swj0
         END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
         call rchkusr()
      END DO
C$OMP FLUSH(ai,bi,bi0,bi2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cgawsmas(y,mask,ni,fix,si2,n1,n2,hakt,lambda,theta,
     1  bi,bi2,bi0,vred,ai,model,kern,spmin,lwght,wght)
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
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit none
      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,model,kern,ni(*)
      logical aws,fix(*),mask(*)
      double precision y(*),theta(*),bi(*),bi0(*),ai(*),lambda,wght,
     1       bi2(*),hakt,lwght(*),si2(*),vred(*),spmin
      integer ih1,ih2,i1,i2,j1,j2,jw1,jw2,jwind2,
     1        iind,jind,jind2,clw1,clw2,dlw1,dlw2
      double precision thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,wj,
     1        hakt2,sv1,sv2,spf,wj0
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      aws=lambda.lt.1d40
C
C   first calculate location weights
C
      ih2=FLOOR(hakt/wght)
      ih1=FLOOR(hakt)
      if(n2.eq.1) ih2=0
      clw1=ih1+1
      clw2=ih2+1
      dlw1=ih1+clw1
      dlw2=ih2+clw2
      z2=0.d0
         DO j2=1,dlw2
            if(n2.gt.1) THEN
               z2=(clw2-j2)*wght
               z2=z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=(j2-1)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=clw1-ih1,clw1+ih1
C  first stochastic term
               jind=j1+jind2
               z1=clw1-j1
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
         call rchkusr()
         DO i2=1,n2
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(ai,bi,bi0,bi2,si2,n1,n2,hakt2,theta
C$OMP& ,lwght,wght,y,fix,mask,ni,vred)
C$OMP& FIRSTPRIVATE(ih1,i2,lambda,aws
C$OMP& ,model,dlw1,clw1,dlw2,clw2)
C$OMP& PRIVATE(iind,thetai,bii,swj
C$OMP& ,swj2,swj0,swjy,sij,sv1,sv2,i1,wj,wj0
C$OMP& ,j2,jw2,jind2,z2,jwind2
C$OMP& ,j1,jw1,jind)
C$OMP DO SCHEDULE(GUIDED)
             DO i1=1,n1
               iind=i1+(i2-1)*n1
               if(.not.mask(iind)) CYCLE
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
               DO jw2=1,dlw2
                  j2=jw2-clw2+i2
                  if(j2.lt.1.or.j2.gt.n2) CYCLE
                  jind2=(j2-1)*n1
                  z2=(clw2-jw2)*wght
                  z2=z2*z2
                  ih1=FLOOR(sqrt(hakt2-z2))
                  jwind2=(jw2-1)*dlw1
                  DO jw1=clw1-ih1,clw1+ih1
C  first stochastic term
                     j1=jw1-clw1+i1
                     if(j1.lt.1.or.j1.gt.n1) CYCLE
                     jind=j1+jind2
                     if(ni(jind).eq.0) CYCLE
                     wj=lwght(jw1+jwind2)
                     swj0=swj0+wj*si2(jind)
                     IF (aws) THEN
                  sij=bii*kldist(model,thetai,theta(jind))
                        IF (sij.gt.1.d0) CYCLE
                        wj=wj*(1.d0-sij)
                     END IF
                     wj0 = wj
                     wj = wj*ni(jind)
                     sv1=sv1+wj
                     sv2=sv2+wj*wj0
                     wj = wj*si2(jind)
                     swj=swj+wj
                     swj2=swj2+wj*wj0
                     swjy=swjy+wj*y(jind)
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
         call rchkusr()
         END DO
C$OMP FLUSH(ai,bi,bi0,bi2,vred)
      RETURN
      END
      subroutine mask(maskin,maskout,n1,n2,h)
      integer n1,n2,h
      logical maskin(n1,n2),maskout(n1,n2)
      integer i1,i2
      DO i1=1,n1
         j1a=max0(1,i1-h)
         j1e=min0(n1,i1+h)
         DO i2=1,n2
            if(.not.maskin(i1,i2)) CYCLE
            j2a=max0(1,i2-h)
            j2e=min0(n2,i2+h)
            DO j1=j1a,j1e
               DO j2=j2a,j2e
                  maskout(j1,j2)=.TRUE.
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
