CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform 3D smoothing on a grid (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smooth3d(y,si2,pos,wlse,nvox,n1,n2,n3,dv,hakt,
     1                    thn,kern,lwght,wght,swjy)
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
      implicit none
      integer nvox,n1,n2,n3,kern,dv,wlse,mask(n1,n2,n3)
      double precision y(nvox,dv),thn(nvox,dv),wght(2),
     1       si2(nvox),hakt,lwght(*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,
     1        clw1,clw2,clw3,dlw1,dlw2,dlw3,k,n,dlw12,iindp,jindp
      double precision swj,swjy(dv),wj,hakt2
      hakt2=hakt*hakt
C
C   first calculate location weights
C
      ih3=FLOOR(hakt/wght(2))
      ih2=FLOOR(hakt/wght(1))
      ih1=FLOOR(hakt)
      n=n1*n2*n3
      dlw1=min(2*n1-1,2*ih1+1)
      dlw2=min(2*n2-1,2*ih2+1)
      dlw3=min(2*n3-1,2*ih3+1)
      dlw12=dlw1*dlw2
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
C
C    get location weights
C
      call locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      call rchkusr()
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               iindp = pos(i1,i2,i3)
               if(iindp.eq.0) CYCLE
               swj=0.d0
               DO k=1,dv
                  swjy(k)=0.d0
               END DO
               DO jw3=1,dlw3
                  j3=jw3-clw3+i3
                  if(j3.lt.1.or.j3.gt.n3) CYCLE
                  DO jw2=1,dlw2
                     j2=jw2-clw2+i2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     DO jw1=1,dlw1
C  first stochastic term
                        j1=jw1-clw1+i1
                        if(j1.lt.1.or.j1.gt.n1) CYCLE
                        jindp = pos(j1,j2,j3)
                        IF(jindp.eq.0) CYCLE
                        wj=lwght(jw1+(jw2-1)*dlw1+(jw3-1)*dlw12)
                        if(wj.le.0.d0) CYCLE
                        if(wlse.ne.0) THEN
                           wj=wj*si2(jindp)
                        END IF
                        swj=swj+wj
                        DO k=1,dv
                           swjy(k)=swjy(k)+wj*y(jindp,k)
                        END DO
                     END DO
                  END DO
               END DO
               DO k=1,dv
                  thn(iindp,k)=swjy(k)/swj
               END DO
               call rchkusr()
            END DO
         END DO
      END DO
      RETURN
      END

C
C   calculate location weights in lwght
C
      subroutine locwghts(dlw1,dlw2,dlw3,wght,hakt2,kern,lwght)
      implicit none
      integer dlw1,dlw2,dlw3,kern
      double precision wght(2),hakt2,lwght(dlw1,dlw2,dlw3),lkern
      external lkern
      double precision z1,z2,z3
      integer j1,j2,j3,clw1,clw2,clw3,ih1,ih2
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
      DO j3=1,dlw3
         Do j2=1,dlw2
            DO j1=1,dlw1
               lwght(j1,j2,j3)=0.d0
            END DO
         END DO
         z3=(clw3-j3)*wght(2)
         z3=z3*z3
         ih2=FLOOR(sqrt(hakt2-z3)/wght(1))
         DO j2=clw2-ih2,clw2+ih2
            IF(j2.lt.1.or.j2.gt.dlw2) CYCLE
            z2=(clw2-j2)*wght(1)
            z2=z3+z2*z2
            ih1=FLOOR(sqrt(hakt2-z2))
            DO j1=clw1-ih1,clw1+ih1
               IF(j1.lt.1.or.j1.gt.dlw1) CYCLE
               z1=clw1-j1
               lwght(j1,j2,j3)=lkern(kern,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      RETURN
      END
C
C
C
