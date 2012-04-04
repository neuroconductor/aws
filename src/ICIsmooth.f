      subroutine lkern1(x,n,h,kern,m,khx)
      implicit logical(a-z)
      integer n,kern,m
      real*8 x(n),h,khx(n)
      integer i
      real*8 mu0, mu2, mu4, xi, xih, skhx
      select case (kern)            
         case (1)
C  Gaussian
            mu0 = 2.506628274631d0
            mu2 = 1.d0
            mu4 = 3.d0
         case (2) 
C  Uniform
            mu0 = 0.5d0
            mu2 = 0.333333333333d0
            mu4 = 0.2d0
         case (3) 
C  Triangular
            mu0 = 1
            mu2 = 0.166666666666d0
            mu4 = 0.066666666666d0
         case (4) 
C  Epanechnikov
            mu0 = 4.d0/3.d0
            mu2 = 0.2d0
            mu4 = 3.d0/3.5d1
         case (5) 
C  Biweight
            mu0 = 1.6d1/1.5d1
            mu2 = 1.d0/7.d0
            mu4 = 1.d0/2.1d1
         case (6) 
C  Triweight
            mu0 = 3.2d1/3.5d1
            mu2 = 1.d0/9.d0
            mu4 = 1.d0/3.3d1
         case default
C  Gaussian
            mu0 = 2.506628274631d0
            mu2 = 1.d0
            mu4 = 3.d0
      end select                    
      skhx = 0.d0
      DO i=1,n
         xi=x(i)
         xih=xi/h
         select case (kern)            
            case (1)
               khx(i)=exp(-xih*xih/2.d0)/mu0
            case (2)
               xih=abs(xih)
               if(xih.le.1.d0) khx(i)=1.d0/mu0           
            case (3)
               xih=abs(xih)
               if(xih.le.1.d0) khx(i)=(1.d0-xih)/mu0         
            case (4)
               if(abs(xih).le.1.d0) THEN
                  xih=xih*xih
                  khx(i)=(1.d0-xih)/mu0
               ENDIF           
            case (5)
               if(abs(xih).le.1.d0) THEN
                  xih=1.d0-xih*xih
                  khx(i)=xih*xih/mu0
               ENDIF           
            case (6)
               if(abs(xih).le.1.d0) THEN
                  xih=1.d0-xih*xih
                  khx(i)=xih*xih*xih/mu0
               ENDIF           
            case default
               khx(i)=exp(-xih*xih/2.d0)/mu0
         end select
         if(m.eq.0) skhx=skhx+khx(i)
         if(m.eq.1) THEN
C  compute first order derivatives 
            khx(i)=mu2*xi/h*khx(i)
            skhx=skhx+khx(i)*xi
         END IF
         IF(m.eq.2) THEN
C  compute second order derivatives
            xi=xi*xi
            khx(i)=(xi-mu2)/(mu4-mu2*mu2)*khx(i)
            skhx=skhx+khx(i)*xi
         END IF
      END DO
      DO i=1,n
         khx(i)=khx(i)/skhx
      END DO
      RETURN
      END
      subroutine sector(x1,n1,x2,n2,nsect,sect,symm,insect)
      implicit logical (a-z)
      logical symm
      integer n1,n2,nsect,sect
      real*8 x1(n1),x2(n2),insect(n1,n2)
      integer i1,i2,isect
      real*8 alpha,xi1,xi2,xnorm,ax
      if(symm) THEN 
         alpha=3.14159265358979d0/nsect
      ELSE
         alpha=6.28318530717959d0/nsect
      END IF
      Do i1=1,n1
         xi1=x1(i1)
         DO i2=1,n2
            xi2=x2(i2)
            xnorm=sqrt(xi1*xi1+xi2*xi2)
            if(xnorm.le.1d-10) THEN
               insect(i1,i2)=1.d0/nsect
            ELSE
               ax=acos(xi1/xnorm)
            if(xi2.lt.0.d0.and.xi1.lt.1.d0) ax=ax+3.14159265358979d0
               isect=ax/alpha
               IF(symm.and.isect.gt.nsect) isect=isect-nsect
               if(isect.eq.sect-1) insect(i1,i2)=1.d0
            END IF
         END DO
      END DO
      RETURN
      END
      subroutine median1d(y,n,yhat)
      implicit logical (a-z)
      integer n
      real*8 y(n),yhat(n)
      integer i
      real*8 ys(3)
      yhat(1)=y(1)
      yhat(n)=y(n)
      DO i=2,n-1
         ys(1)=y(i-1)
         ys(2)=y(i)
         ys(3)=y(i+1)
         call qsort3(ys,1,3)
         yhat(i)=ys(2)
      END DO
      RETURN
      END
      subroutine median2d(y,n1,n2,yhat)
      implicit logical (a-z)
      integer n1,n2
      real*8 y(n1,n2),yhat(n1,n2)
      integer i1,i2
      real*8 ys(9)
      i1 = 1
      DO i2=1,n2
         yhat(i1,i2)=y(i1,i2)
      END DO
      i1 = n1
      DO i2=1,n2
         yhat(i1,i2)=y(i1,i2)
      END DO
      i2 = 1
      DO i1=1,n1
         yhat(i1,i2)=y(i1,i2)
      END DO
      i2 = n2
      DO i1=1,n1
         yhat(i1,i2)=y(i1,i2)
      END DO
      DO i1=2,n1-1
         DO i2=2,n2-1
            ys(1)=y(i1-1,i2-1)
            ys(2)=y(i1,i2-1)
            ys(3)=y(i1+1,i2-1)
            ys(4)=y(i1-1,i2)
            ys(5)=y(i1,i2)
            ys(6)=y(i1+1,i2)
            ys(7)=y(i1-1,i2+1)
            ys(8)=y(i1,i2+1)
            ys(9)=y(i1+1,i2+1)
            call qsort3(ys,1,9)
            yhat(i1,i2)=ys(5)
         END DO
      END DO
      RETURN
      END
      subroutine median3d(y,n1,n2,n3,yhat)
      implicit logical (a-z)
      integer n1,n2,n3
      real*8 y(n1,n2,n3),yhat(n1,n2,n3)
      integer i1,i2,i3
      real*8 ys(27)
      i1 = 1
      DO i2=1,n2
         DO i3=1,n3
            yhat(i1,i2,i3)=y(i1,i2,i3)
         END DO
      END DO
      i1 = n1
      DO i2=1,n2
         DO i3=1,n3
            yhat(i1,i2,i3)=y(i1,i2,i3)
         END DO
      END DO
      i2 = 1
      DO i1=1,n1
         DO i3=1,n3
            yhat(i1,i2,i3)=y(i1,i2,i3)
         END DO
      END DO
      i2 = n2
      DO i1=1,n1
         DO i3=1,n3
            yhat(i1,i2,i3)=y(i1,i2,i3)
         END DO
      END DO
      i3 = 1
      DO i1=1,n1
         DO i2=1,n2
            yhat(i1,i2,i3)=y(i1,i2,i3)
         END DO
      END DO
      i2 = n3
      DO i1=1,n1
         DO i2=1,n2
            yhat(i1,i2,i3)=y(i1,i2,i3)
         END DO
      END DO
      DO i1=2,n1-1
         DO i2=2,n2-1
            DO i3=2,n3-1
               ys(1)=y(i1-1,i2-1,i3-1)
               ys(2)=y(i1,i2-1,i3-1)
               ys(3)=y(i1+1,i2-1,i3-1)
               ys(4)=y(i1-1,i2,i3-1)
               ys(5)=y(i1,i2,i3-1)
               ys(6)=y(i1+1,i2,i3-1)
               ys(7)=y(i1-1,i2+1,i3-1)
               ys(8)=y(i1,i2+1,i3-1)
               ys(9)=y(i1+1,i2+1,i3-1)
               ys(10)=y(i1-1,i2-1,i3)
               ys(11)=y(i1,i2-1,i3)
               ys(12)=y(i1+1,i2-1,i3)
               ys(13)=y(i1-1,i2,i3)
               ys(14)=y(i1,i2,i3)
               ys(15)=y(i1+1,i2,i3)
               ys(16)=y(i1-1,i2+1,i3)
               ys(17)=y(i1,i2+1,i3)
               ys(18)=y(i1+1,i2+1,i3)
               ys(19)=y(i1-1,i2-1,i3+1)
               ys(20)=y(i1,i2-1,i3+1)
               ys(21)=y(i1+1,i2-1,i3+1)
               ys(22)=y(i1-1,i2,i3+1)
               ys(23)=y(i1,i2,i3+1)
               ys(24)=y(i1+1,i2,i3+1)
               ys(25)=y(i1-1,i2+1,i3+1)
               ys(26)=y(i1,i2+1,i3+1)
               ys(27)=y(i1+1,i2+1,i3+1)
               call qsort3(ys,1,27)
               yhat(i1,i2,i3)=ys(14)
            END DO
         END DO
      END DO
      RETURN
      END
      
      