CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C          Compute the Kullback-Leibler Distance
C
C          Gaussian, Diagonal covariance matrix
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function kldistd(thi,thj,n,wght,nwght)
      implicit logical (a-z)
      integer n,nwght,i,k,thi(1),thj(1)
      real*8 z,wght(nwght)
      kldistd=0.d0
      i=1
         z=thi(i)-thj(i)
         kldistd=z*z
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
C
C    Estimate variance parameters
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine esigmac(y,n,theta,bi,quant,varcoef,mvar)
      implicit logical (a-z)
      integer n,y(n),theta(n),quant
      real*8 bi(n),varcoef,mvar
      integer i,k
      real*8 z,bii,sumres,sumwght,wght,res
      sumres=0.d0
      sumwght=0.d0
      DO i=1,n
         bii=bi(i)
         if(bii.le.1.d0.or.y(i).ge.quant) CYCLE
         wght=bii-1.d0
         res=(y(i)-theta(i))
         res=res*res*bii/wght
         sumres=sumres+res*wght
         sumwght=sumwght+wght
      END DO
      z=sumres/sumwght
      varcoef=z
      mvar=z
      RETURN
      END
      subroutine esigmal(y,n,theta,bi,quant,varcoef,mvar)
      implicit logical (a-z)
      integer n,y(n),theta(n),quant
      real*8 bi(n),varcoef(2),mvar,res
      integer i,k
      real*8 z,bii,s0,s1,s2,t0,t1,d,wght,thi,mth
      s0=0.d0
      s1=0.d0
      s2=0.d0
      t0=0.d0
      t1=0.d0
      mth=0.d0
      DO i=1,n
         bii=bi(i)
         mth=mth+theta(i)
         if(bii.le.1.d0.or.y(i).ge.quant) CYCLE
         wght=bii-1.d0
         thi=theta(i)
         res=(y(i)-thi)
         res=res*res*bii/wght
         s0=s0+wght
         z=wght*thi
         s1=s1+z
         s2=s2+z*thi
         t0=t0+wght*res
         t1=t1+z*res
      END DO
      d=s2*s0-s1*s1
      varcoef(1)=(s2*t0-s1*t1)/d
      varcoef(2)=(-s1*t0+s0*t1)/d
      mvar=varcoef(1)+varcoef(2)*mth/n
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
C
C    Estimate correlations
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine estcorr1(res,n,scorr)
      implicit logical (a-z)
      integer n
      real*8 res(n),scorr(1)
      integer i,j,k,n,m,l
      real*8 vres,z,z1,z2,resij
      z1=0.d0
      z2=0.d0
      DO i=1,n
         resij=res(i)
         z1=z1+resij
         z2=z2+resij*resij
      END DO
      z2=z2/n
      z1=z1/n
      vres=n/(n-1)*(z2-z1*z1)
      DO i=1,n1
         res(i)=res(i)-z1
      END DO
      z=0.d0
      DO i=1,n1-1
         z=z+res(i)*res(i+1)
      END DO
      RETURN
      END
      subroutine estcorr2(res,n1,n2,scorr)
      implicit logical (a-z)
      integer n1,n2
      real*8 res(n1,n2),scorr(2)
      integer i,j,k,n,m,l
      real*8 vres(4),z,z1,z2,resij
      n=n1*n2
      z1=0.d0
      z2=0.d0
      DO i=1,n1
         DO j=1,n2
            resij=res(i,j)
            z1=z1+resij
            z2=z2+resij*resij
         END DO
      END DO
      z2=z2/n
      z1=z1/n
      vres=n/(n-1)*(z2-z1*z1)
      DO i=1,n1
         DO j=1,n2
            res(i,j)=res(i,j)-z1
         END DO
      END DO
      z=0.d0
      DO i=1,n1-1
         DO j=1,n2
            z=z+res(i,j)*res(i+1,j)
         END DO
      END DO
      scorr(1)=z/n2/(n1-1)/vres
      z=0.d0
      DO i=1,n1
         DO j=1,n2-1
            z=z+res(i,j)*res(i,j+1)
         END DO
      END DO
      scorr(2)=z/n1/(n2-1)/vres
      RETURN
      END
      subroutine estcorr3(res,n1,n2,n3,scorr)
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 res(n1,n2),scorr(3)
      integer i,j,k,n,m,l
      real*8 vres(4),z,z1,z2,resij
      n=n1*n2*n3
      z1=0.d0
      z2=0.d0
      DO i=1,n1
         DO j=1,n2
            DO k=1,n3
               resij=res(i,j,k)
               z1=z1+resij
               z2=z2+resij*resij
            END DO
         END DO
      END DO
      z2=z2/n
      z1=z1/n
      vres=n/(n-1)*(z2-z1*z1)
      DO i=1,n1
         DO j=1,n2
            DO k=1,n3
               res(i,j,k)=res(i,j,k)-z1
            END DO
         END DO
      END DO
      z=0.d0
      DO i=1,n1-1
         DO j=1,n2
            DO k=1,n3
               z=z+res(i,j,k)*res(i+1,j,k)
            END DO
         END DO
      END DO
      scorr(1)=z/n2/n3/(n1-1)/vres
      z=0.d0
      DO i=1,n1
         DO j=1,n2-1
            DO k=1,n3
               z=z+res(i,j,k)*res(i,j+1,k)
            END DO
         END DO
      END DO
      scorr(2)=z/n1/n3/(n2-1)/vres
      z=0.d0
      DO i=1,n1
         DO j=1,n2
            DO k=1,n3-1
               z=z+res(i,j,k)*res(i,j,k+1)
            END DO
         END DO
      END DO
      scorr(3)=z/n1/n2/(n3-1)/vres
      RETURN
      END
