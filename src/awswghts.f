CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C   compute weighting scheme
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawsw(n1,n2,n3,hakt,lambda,theta,bi,
     1                model,kern,spmin,lwght,wght)
C   compute weights for all combinations of design points
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
C
      implicit none

      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,n3,model,kern
      double precision theta(*),bi(*),lambda,
     1       wght(n1,n2,n3,n1,n2,n3),hakt,lwght(*),spmin,spf
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12
      double precision thetai,bii,sij,z1,z2,z3,wj,hakt2
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
C
C   first calculate location weights
C
      ih3=FLOOR(hakt)
      ih2=ih3
      ih1=ih3
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
            z3=j3
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3))
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2
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
C$OMP& SHARED(bi,n1,n2,n3,hakt2,theta,
C$OMP& ih3,lwght,wght)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,n12,
C$OMP& model,spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12)
C$OMP& PRIVATE(i1,i2,i3,iind,thetai,bii,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
         thetai=theta(iind)
         bii=bi(iind)/lambda
C   scaling of sij outside the loop
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3))
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  wj=lwght(jw1+clw1+1+jwind2)
                  sij=bii*kldist(model,thetai,theta(jind))
                  IF (sij.gt.1.d0) CYCLE
                  IF (sij.gt.spmin) THEN
                     wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
                  wght(i1,i2,i3,j1,j2,j3) = wj
               END DO
            END DO
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(wght)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws
C   for selected points
C   compute weighting scheme
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawsw1(n1,n2,n3,inx,iny,inz,anz,hakt,lambda,theta,bi,
     1                model,kern,spmin,lwght,wght)
C   compute weights for all combinations of design points
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
C
      implicit none

      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,n3,model,kern,anz,inx(anz),iny(anz),inz(anz)
      double precision theta(*),bi(*),lambda,
     1       wght(n1,n2,n3,anz),hakt,lwght(*),spmin,spf
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12
      double precision thetai,bii,sij,z1,z2,z3,wj,hakt2
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
C
C   first calculate location weights
C
      ih3=FLOOR(hakt)
      ih2=ih3
      ih1=ih3
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
            z3=j3
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3))
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2
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
      DO iind=1,anz
         i1=inx(iind)
         i2=iny(iind)
         i3=inz(iind)
         thetai=theta(i1+(i2-1)*n1+(i3-1)*n1*n2)
         bii=bi(iind)/lambda
C   scaling of sij outside the loop
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3))
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  wj=lwght(jw1+clw1+1+jwind2)
                  sij=bii*kldist(model,thetai,theta(jind))
                  IF (sij.gt.1.d0) CYCLE
                  IF (sij.gt.spmin) THEN
                     wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
                  wght(j1,j2,j3,iind) = wj
               END DO
            END DO
         END DO
      END DO
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Patch based aws compute weighting scheme
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pawswght(n1,n2,n3,i1,i2,i3,hakt,lambda,theta,bi,
     1              model,kern,spmin,lwght,wght,npsize,wi)
C
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (input)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   npsize       patch size
C
      implicit none

      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,n3,i1,i2,i3,model,kern,npsize
      double precision theta(*),bi(*),lambda,wght(2),
     1       wi(*),hakt,lwght(*),spmin,spf
      integer ih1,ih2,ih3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,ip1,ip2,ip3,nph1,nph2,nph3,
     3        ipind,jp1,jp2,jp3,jpind,np1,np2,np3
      double precision sij,z1,z2,z3,wj,
     1       hakt2,w1,w2
       hakt2=hakt*hakt
       spf=1.d0/(1.d0-spmin)
C
C   first calculate location weights
C
       w1=wght(1)
       w2=wght(2)
       ih3=FLOOR(hakt/w2)
       ih2=FLOOR(hakt/w1)
       ih1=FLOOR(hakt)
       np1=2*npsize+1
       np2=min(n2,2*npsize+1)
       np3=min(n3,2*npsize+1)
       if(n3.eq.1) ih3=0
       if(n2.eq.1) ih2=0
       clw1=ih1
       clw2=ih2
       clw3=ih3
       dlw1=ih1+clw1+1
       dlw2=ih2+clw2+1
       dlw3=ih3+clw3+1
       dlw12=dlw1*dlw2
       nph1=(np1-1)/2
       nph2=(np2-1)/2
       nph3=(np3-1)/2
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
                jind=j1+clw1+1+jind2
                z1=j1
                lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
          END DO
      END DO
      call rchkusr()
      iind=i1+(i2-1)*n1+(i3-1)*n12
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
                sij=0.d0
                DO ip1=i1-nph1,i1+nph1
                    if(sij.gt.1.d0) CYCLE
                    if(ip1.le.0.or.ip1.gt.n1) CYCLE
                    jp1=ip1+jw1
                    if(jp1.le.0.or.jp1.gt.n1) CYCLE
                    DO ip2=i2-nph2,i2+nph2
                      if(sij.gt.1.d0) CYCLE
                      if(ip2.le.0.or.ip2.gt.n2) CYCLE
                      jp2=ip2+jw2
                      if(jp2.le.0.or.jp2.gt.n2) CYCLE
                      DO ip3=i3-nph3,i3+nph3
                          if(sij.gt.1.d0) CYCLE
                          if(ip3.le.0.or.ip3.gt.n3) CYCLE
                          jp3=ip3+jw3
                          if(jp3.le.0.or.jp3.gt.n3) CYCLE
                          jpind=jp1+(jp2-1)*n1+(jp3-1)*n12
                          ipind=ip1+(ip2-1)*n1+(ip3-1)*n12
                          sij=max(sij,bi(ipind)/lambda*
     1                kldist(model,theta(ipind),theta(jpind)))
                      END DO
                    END DO
                END DO
                IF(sij.gt.1.d0) CYCLE
                IF (sij.gt.spmin) THEN
                    wj=wj*(1.d0-spf*(sij-spmin))
                END IF
                wi(jind) = wj
            END DO
        END DO
      END DO
      RETURN
      END
C
C
C
