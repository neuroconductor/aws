CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded) with variance - mean model
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cgawsdti(y,mask,si2,n1,n2,n3,hakt,lambda,theta,
     1        bi,gi,gi2,thetan,lwght,wght)
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
      integer n1,n2,n3
      integer mask(n1,n2,n3)
      double precision y(n1,n2,n3),theta(n1,n2,n3),bi(n1,n2,n3),
     1       thetan(n1,n2,n3),lambda,wght(2),hakt,lwght(*),
     2       si2(n1,n2,n3),spmin,gi(n1,n2,n3),gi2(n1,n2,n3)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3
      double precision thetai,bii,sij,swj,swjy,z1,z2,z3,wj,hakt2,
     1        sv1,sv2,spf,z
      hakt2=hakt*hakt
      spmin=0.25d0
      spf=1.d0/(1.d0-spmin)
      ih1=int(hakt)
C
C   first calculate location weights
C
      ih3=int(hakt/wght(2))
      ih2=int(hakt/wght(1))
      ih1=int(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1+1
      clw2=ih2+1
      clw3=ih3+1
      dlw1=ih1+clw1
      dlw2=ih2+clw2
      dlw3=ih3+clw3
      z2=0.d0
      z3=0.d0
      DO j3=1,dlw3
         if(n3.gt.1) THEN
            z3=(clw3-j3)*wght(2)
            z3=z3*z3
            ih2=int(sqrt(hakt2-z3)/wght(1))
            jind3=(j3-1)*dlw1*dlw2
         ELSE
            jind3=0
         END IF
         DO j2=clw2-ih2,clw2+ih2
            if(n2.gt.1) THEN
               z2=(clw2-j2)*wght(1)
               z2=z3+z2*z2
               ih1=int(sqrt(hakt2-z2))
               jind2=jind3+(j2-1)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=clw1-ih1,clw1+ih1
C  first stochastic term
               jind=j1+jind2
               z1=clw1-j1
               lwght(jind)=max(0.d0,1.d0-(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      call rchkusr()
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
               IF (mask(i1,i2,i3).eq.0) CYCLE
C    nothing to do, final estimate is already fixed by control
               thetai=theta(i1,i2,i3)
               bii=bi(i1,i2,i3)/lambda
C   scaling of sij outside the loop
               ih3=int(hakt/wght(2))
               swj=0.d0
               swjy=0.d0
               sv1=0.d0
               sv2=0.d0
               DO jw3=1,dlw3
                  j3=jw3-clw3+i3
                  if(j3.lt.1.or.j3.gt.n3) CYCLE
                  z3=(clw3-jw3)*wght(2)
                  z3=z3*z3
                  if(n2.gt.1) ih2=int(sqrt(hakt2-z3)/wght(1))
                  jwind3=(jw3-1)*dlw1*dlw2
                  DO jw2=clw2-ih2,clw2+ih2
                     j2=jw2-clw2+i2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     z2=(clw2-jw2)*wght(1)
                     z2=z3+z2*z2
                     ih1=int(sqrt(hakt2-z2))
                     jwind2=jwind3+(jw2-1)*dlw1
                     DO jw1=clw1-ih1,clw1+ih1
C  first stochastic term
                        j1=jw1-clw1+i1
                        if(j1.lt.1.or.j1.gt.n1) CYCLE
                        if(mask(j1,j2,j3).eq.0) CYCLE
                        wj=lwght(jw1+jwind2)
                        z1=(clw1-jw1)
                        z1=z2+z1*z1
                        z=(thetai-theta(j1,j2,j3))
                        sij=bii*z*z
                        IF (sij.gt.1.d0) CYCLE
                        IF (sij.gt.spmin) THEN
                           wj=wj*(1.d0-spf*(sij-spmin))
                        END IF
                        sv1=sv1+wj
                        sv2=sv2+wj*wj
                        swj=swj+wj*si2(j1,j2,j3)
                        swjy=swjy+wj*si2(j1,j2,j3)*y(j1,j2,j3)
                     END DO
                  END DO
               END DO
               thetan(i1,i2,i3)=swjy/swj
               bi(i1,i2,i3)=swj
               gi(i1,i2,i3)=sv1
               gi2(i1,i2,i3)=sv2
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END
