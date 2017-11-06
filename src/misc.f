      subroutine getvofh(bw,kern,wght,vol)
C
C   wght(1) is voxel extension x / voxel extension y,  i.e. zero in univariate situations
C   wght(2) is voxel extension x / voxel extension z,  i.e. zero in univariate and bivariate situations
      implicit logical(a-z)
      integer kern
      double precision bw,wght(2),vol,sofw
      external sofw
      vol=sofw(bw,kern,wght)
      RETURN
      END
      double precision function sofw(bw,kern,wght)
      implicit logical(a-z)
      integer kern
      double precision bw,wght(2)
      integer j1,j2,j3,dlw1,dlw2,dlw3,clw1,clw2,clw3,ih1,ih2,ih3
      double precision sw,sw2,h2,lkern,z1,z2,z3,z
      external lkern
      h2=bw*bw
C
C   first calculate location weights
C
      ih3=FLOOR(bw*wght(2))
      ih2=FLOOR(bw*wght(1))
      ih1=FLOOR(bw)
      dlw1=2*ih1+1
      dlw2=2*ih2+1
      dlw3=2*ih3+1
      clw1=(dlw1+1)/2
      clw2=(dlw2+1)/2
      clw3=(dlw3+1)/2
      sw=0.d0
      sw2=0.d0
      DO j1=1,dlw1
         z1=(clw1-j1)
         z1=z1*z1
         if(wght(1).gt.0.d0) THEN
            ih2=FLOOR(sqrt(h2-z1)*wght(1))
            DO j2=clw2-ih2,clw2+ih2
               z2=(clw2-j2)/wght(1)
               z2=z1+z2*z2
               if(wght(2).gt.0.d0) THEN
                  ih3=FLOOR(sqrt(h2-z2)*wght(2))
                  DO j3=clw3-ih3,clw3+ih3
                     z3=(clw3-j3)/wght(2)
                     z=lkern(kern,(z3*z3+z2)/h2)
                     sw=sw+z
                     sw2=sw2+z*z
                  END DO
               ELSE
                  z=lkern(kern,z2/h2)
                  sw=sw+z
                  sw2=sw2+z*z
               END IF
            END DO
         ELSE
            z=lkern(kern,z1/h2)
            sw=sw+z
            sw2=sw2+z*z
         END IF
      END DO
      sofw=sw*sw/sw2
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   determine sum of location weights for a given geometry a(3) and given
C   bandwidth
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Algorithmus zur Nullstellenbestimmung einer monotonen Funktion auf(0,\infty)
      subroutine gethani(x,y,kern,value,wght,eps,bw)
      implicit logical(a-z)
      integer kern
      double precision x,y,value,wght(2),eps,bw
      double precision fw1,fw2,fw3,z
      double precision sofw
      external sofw
      if(x.ge.y) RETURN
      fw1=sofw(x,kern,wght)
      fw2=sofw(y,kern,wght)
      DO WHILE(fw1.gt.value)
         x=x*x/y
         fw1=sofw(x,kern,wght)
      END DO
      DO WHILE(fw2.le.value)
         y=y*y/x
         fw2=sofw(y,kern,wght)
      END DO
      DO WHILE(min(fw2/value,value/fw1).gt.1.d0+eps.and.
     1             abs(y-x).gt.1d-6)
         z=x+(value-fw1)/(fw2-fw1)*(y-x)
         fw3=sofw(z,kern,wght)
         if(fw3.le.value) THEN
            x=z
            fw1=fw3
         ENDIF
         if(fw3.ge.value) THEN
            y=z
            fw2=fw3
         ENDIF
         call rchkusr()
      END DO
      if(fw2/value.gt.value/fw1) THEN
          bw=x+(value-fw1)/(fw2-fw1)*(y-x)
      ELSE
          bw=y-(fw2-value)/(fw2-fw1)*(y-x)
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Calculate exceedence probabilities in awstestprop
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exceed(x,n,z,nz,exprob)
      implicit logical (a-z)
      integer n,nz
      double precision x(n),z(nz),exprob(nz)
      integer i,j,k
      double precision sk,zj
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n,nz,x,z,exprob)
C$OMP& PRIVATE(i,j,k,sk,zj)
C$OMP DO SCHEDULE(GUIDED)
      DO j=1,nz
         k=0
         zj=z(j)
         DO i=1,n
            if(x(i).gt.zj) k=k+1
         END DO
         sk = k
         exprob(j)=sk/n
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(exprob)
      Return
      End
