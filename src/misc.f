CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Hypergeometric 1F1 NIST HB 13.2.2, 13.2.39, 13.2(iv)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine hg1f1(a,b,z,n,fz)
      implicit none
      integer n
      double precision a,b,z(n),fz(n)
      integer i
      double precision x,y,d,eps,zi,ezi,ai,gofbai
      double precision gammaf
      external gammaf
      eps=1.d-15
      gofbai=gammaf(b)/gammaf(b-a)
      DO i=1,n
         d = 1.d0
         zi = z(i)
         IF(zi.lt.0) THEN
            ezi=exp(zi/2)
            ai=b-a
            if(zi.lt.-1400) THEN
               fz(i) = exp((-a)*log(-zi))*gofbai+5.6e-3+1.9e-3*b
C   add +5.6e-3+1.9e-3*b to keep the function monotone
               CYCLE
            END IF
         ELSE
            ezi=1.d0
            ai=a
         ENDIF
         x = ezi
         y = ezi
         DO WHILE (abs(y).gt.abs(x)*eps)
            y = -y*(ai+d-1.d0)/(b+d-1.d0)*zi/d
            x = x+y
            d = d+1.d0
         END DO
         fz(i) = ezi*x
      END DO
      RETURN
      END

      subroutine getvofh(bw,kern,wght,vol)
C
C   wght(1) is voxel extension x / voxel extension y,  i.e. zero in univariate situations
C   wght(2) is voxel extension x / voxel extension z,  i.e. zero in univariate and bivariate situations
      implicit none
      integer kern
      double precision bw,wght(2),vol,sofw
      external sofw
      vol=sofw(bw,kern,wght)
      RETURN
      END
      double precision function sofw(bw,kern,wght)
      implicit none
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
      implicit none
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
      implicit none
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Calculate exceedence probabilities in awstestprop
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exceedm(x,n,z,nz,exprob,mask)
      implicit none
      integer n,nz
      double precision x(n),z(nz),exprob(nz)
      integer mask(n)
      integer i,j,k,ni
      double precision sk,zj
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(n,nz,x,z,exprob,mask)
C$OMP& PRIVATE(i,j,k,sk,zj,ni)
C$OMP DO SCHEDULE(GUIDED)
      DO j=1,nz
        k=0
        zj=z(j)
        ni=0
        DO i=1,n
          if(mask(i).ne.0) THEN
             if(x(i).gt.zj) k=k+1
             ni=ni+1
          END IF
        END DO
        sk = k
        exprob(j)=sk/ni
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(exprob)
      Return
      End

      subroutine imcorr(res,mask,n1,n2,n3,nv,scorr,l1,l2,l3)

      implicit none
      integer n1,n2,n3,nv,l1,l2,l3,lag(3)
      double precision scorr(l1,l2,l3),res(nv,n1,n2,n3)
      integer mask(n1,n2,n3)
      integer i1,i2,i3
      Do i1=1,l1
         lag(1)=i1-1
         DO i2=1,l2
            lag(2)=i2-1
            DO i3=1,l3
               lag(3)=i3-1
               call imcorrl(res,mask,n1,n2,n3,nv,scorr(i1,i2,i3),lag)
               call rchkusr()
            END DO
         END DO
      END DO
      return
      end

      subroutine imcorrl(res,mask,n1,n2,n3,nv,scorr,lag)

      implicit none
      integer n1,n2,n3,nv,lag(3)
      double precision scorr,res(nv,n1,n2,n3)
      integer mask(n1,n2,n3)
      double precision z2,y2,resi,resip1,vrm,vrmp1,zk,zcorr,z
      integer i1,i2,i3,i4,l1,l2,l3,k
      zk=nv
      l1=lag(1)
      l2=lag(2)
      l3=lag(3)
      z=0.d0
      k=0
C  correlation in x
      do i1=1,n1-l1
         do i2=1,n2-l2
            do i3=1,n3-l3
         if ((mask(i1,i2,i3)*mask(i1+l1,i2+l2,i3+l3)).eq.0) CYCLE
               z2=0.d0
               y2=0.d0
               zcorr=0.d0
               do i4=1,nv
                  resi=res(i4,i1,i2,i3)
                  resip1=res(i4,i1+l1,i2+l2,i3+l3)
                  z2=z2+resi*resi
                  y2=y2+resip1*resip1
                  zcorr=zcorr+resi*resip1
               enddo
               vrm=z2/zk
               vrmp1=y2/zk
               vrm=vrm*vrmp1
               if(vrm.gt.1e-10) THEN
                  z=z+zcorr/zk/sqrt(vrm)
                  k=k+1
               end if
            enddo
         enddo
      enddo
      if(k.ne.0) THEN
         scorr=z/k
      ELSE
        scorr=0.d0
      END IF
      return
      end

      subroutine ivar(res,resscale,nvoxel,nt,var)
C
C   compute variance estimates from residuals in spm !!! (not the inverse)
C
      implicit none
      integer nvoxel,nt
      double precision resscale,var(nvoxel),res(nt,nvoxel)
      double precision z2,zk,resi,ressc2,z1
      integer i,it
      zk=nt
      ressc2=resscale*resscale

C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(res,resscale,nvoxel,nt,var,zk,ressc2)
C$OMP& PRIVATE(i,it,z1,z2,resi)
C$OMP DO SCHEDULE(GUIDED)
      do i=1,nvoxel
         z2=0.d0
         z1=0.d0
         do it=1,nt
            resi=res(it,i)
            z1=z1+resi
            z2=z2+resi*resi
         end do
         z1 = z1/zk
         z2 = z2/zk
         var(i)=(z2-z1*z1)*ressc2
      end do
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(var)
      return
      end

      subroutine sweepm(res,nvoxel,nt)

      implicit none
      integer nvoxel,nt
      double precision res(nt,nvoxel)
      integer i,k
      double precision z

C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(res,nvoxel,nt)
C$OMP& PRIVATE(i,k,z)
C$OMP DO SCHEDULE(GUIDED)
      Do i=1,nvoxel
         z=0.d0
         DO k=1,nt
            z=z+res(k,i)
         END DO
         z=z/nt
         DO k=1,nt
            res(k,i)=res(k,i)-z
         END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(res)
      return
      end
