CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C   Kullback-Leibler Distance of two general Gaussian distributions 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C      real*8 function klmnorm(mu1,mu2,s1,s2,d,s2inv)
      real*8 function klmnorm(mu1,mu2,s2,s1,d,s2inv)
      implicit logical (a-z)
      integer d,i,j,info
      real*8 s1(d,d),s2(d,d),s2inv(d,d),mu1(d),mu2(d)
      real*8 z,det2,det1,mud
      call deta(s1,d,s2inv,det1,info)
      if(info.gt.0) THEN
         klmnorm=0.d0
	 RETURN
      END IF
      call invdeta(s2,d,s2inv,det2,info)
      if(info.gt.0) THEN
         klmnorm=0.d0
	 RETURN
      END IF
      z=det2/det1
      z=dlog(z)-d
      do i=1,d
         do j=1,d
            z=z+s1(i,j)*s2inv(i,j)
         END DO
      END DO
      DO i=1,d
         mud=mu1(i)-mu2(i)
	 z=z+mud*mud*s2inv(i,i)
	 IF(i.eq.d) CYCLE
         DO j=i+1,d
	    z=z+2.d0*mud*(mu1(j)-mu2(j))*s2inv(i,j)
         END DO
      END DO
      klmnorm=z
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C   Kullback-Leibler Distance of two general Gaussian distributions 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C      subroutine klmnormn(mu1,mu2,th1,th2,d,ds,n,s1,s2,s2inv,kld)
      subroutine klmnormn(mu1,mu2,th1,th2,d,ds,n,s2,s1,s2inv,kld)
      implicit logical (a-z)
      integer d,ds,n,i,j,l,m,info
      real*8 th1(ds,n),th2(ds,n),s1(d,d),s2(d,d),s2inv(d,d),z,
     1       det2,det1,kld(n),mu1(d,n),mu2(d,n),mud
      DO l=1,n
         m=1
         DO i=1,d
	    DO j=1,i
	       s1(i,j)=th1(m,l)
	       s1(j,i)=th1(m,l)
	       s2(i,j)=th2(m,l)
	       s2(j,i)=th2(m,l)
	       m=m+1
	    END DO
	 END DO
         call deta(s1,d,s2inv,det1,info)
	 IF(info.gt.0) THEN
	    kld(l)=0
C	    call intpr("l1",2,l,1)
	    CYCLE
	 END IF
         call invdeta(s2,d,s2inv,det2,info)
	 IF(info.gt.0) THEN
	    kld(l)=0
C	    call intpr("l2",2,l,1)
	    CYCLE
	 END IF
         z=dlog(det2/det1)-d
         do i=1,d
            do j=1,d
               z=z+s1(i,j)*s2inv(i,j)
            END DO
         END DO
         DO i=1,d
            mud=mu1(i,l)-mu2(i,l)
	    z=z+mud*mud*s2inv(i,i)
	    IF(i.eq.d) CYCLE
            DO j=i+1,d
	       z=z+2.d0*mud*(mu1(j,l)-mu2(j,l))*s2inv(i,j)
            END DO
         END DO
         kld(l)=z
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        one iteration of univariate local constant aws for general Gaussian case
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cvawsnun(y,ys,fix,n,d,ds,hakt,lambda,mu,sigma,bi,bi0,
     1                    ami,asi,kern,spmax,swjs,swjm,work,s1,s2)
C   
C   y        observed values of regression function
C   n        number of observations
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C
      implicit logical (a-z)
      external klmnorm,lkern
      real*8 klmnorm,lkern
      integer n,kern,d,ds
      logical aws,fix(n)
      real*8 y(d,n),ys(ds,n),hakt,sigma(ds,n),mu(d,n),bi(n),bi0(n),
     1       ami(d,n),lambda,spmax,swjs(ds),swjm(d),work(d,d),s1(d,d),
     2       asi(ds,n),s2(d,d)
      integer ih,i,j,k,l,m,ja,je
      real*8 bii,sij,swj,z,wj,swj0,bii0
      ih=hakt
      aws=lambda.lt.1d40
      DO i=1,n
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         bii=bi(i)/lambda
C   scaling of sij outside the loop
         bii0=bi0(i)
         ja=max0(1,i-ih)
         je=min0(n,i+ih)
         swj=0.d0
         swj0=0.d0
	 DO k=1,ds
            swjs(k)=0.d0
	 END DO
	 DO k=1,d
            swjm(k)=0.d0
	 END DO
         DO j=ja,je
C  first stochastic term
            z=(i-j)/hakt
            wj=lkern(kern,z*z)
            swj0=swj0+wj
            IF (aws) THEN
               m=1
               DO k=1,d
	          DO l=1,k
	             s1(l,k)=sigma(m,i)
	             s1(k,l)=sigma(m,i)
	             s2(l,k)=sigma(m,j)
	             s2(k,l)=sigma(m,j)
	             m=m+1
	          END DO
	       END DO
               sij=bii*klmnorm(mu(1,i),mu(1,j),s1,s2,d,work)
               IF (sij.gt.spmax) CYCLE
               wj=wj*dexp(-sij)
            ENDIF
            swj=swj+wj
	    DO k=1,ds
               swjs(k) = swjs(k)+wj*ys(k,j)
	    END DO
	    DO k=1,d
               swjm(k) = swjm(k)+wj*y(k,j)
	    END DO
         END DO
	 DO k=1,ds
           asi(k,i)=swjs(k)
	 END DO
	 DO k=1,d
           ami(k,i)=swjm(k)
	 END DO
         bi(i)=swj
         bi0(i)=swj0
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        one iteration of univariate local constant aws for general Gaussian case
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cvawsnbi(y,ys,fix,n1,n2,d,ds,hakt,lambda,mu,sigma,bi,
     1               bi0,ami,asi,kern,spmax,swjs,swjm,work,s1,s2,wght)
C   
C   y        observed values of regression function
C   n        number of observations
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C
      implicit logical (a-z)
      external klmnorm,lkern
      real*8 klmnorm,lkern
      integer n1,n2,kern,d,ds
      logical aws,fix(n1,n2)
      real*8 y(d,n1,n2),ys(ds,n1,n2),hakt,sigma(ds,n1,n2),mu(d,n1,n2),
     1       bi(n1,n2),bi0(n1,n2),ami(d,n1,n2),lambda,spmax,swjs(ds),
     2       swjm(d),work(d,d),s1(d,d),asi(ds,n1,n2),s2(d,d),wght
      integer ih1,ih2,i1,i2,j1,j2,k,l,m,ja1,je1,ja2,je2
      real*8 bii,sij,swj,z1,z2,wj,swj0,bii0,hakt2
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i1=1,n1
         DO i2=1,n2
            IF (fix(i1,i2)) CYCLE
C    nothing to do, final estimate is already fixed by control 
            bii=bi(i1,i2)/lambda
C   scaling of sij outside the loop
            bii0=bi0(i1,i2)
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            swj=0.d0
            swj0=0.d0
	    DO k=1,ds
               swjs(k)=0.d0
	    END DO
	    DO k=1,d
               swjm(k)=0.d0
	    END DO
            DO j1=ja1,je1
               z1=(i1-j1)
               z1=z1*z1
               ih2=dsqrt(hakt2-z1)/wght
               ja2=max0(1,i2-ih2)
               je2=min0(n2,i2+ih2)
               DO j2=ja2,je2
C  first stochastic term
                  z2=(i2-j2)*wght
                  wj=lkern(kern,z1+z2*z2)
                  swj0=swj0+wj
                  IF (aws) THEN
                     m=1
                     DO k=1,d
	                DO l=1,k
	                   s1(l,k)=sigma(m,i1,i2)
	                   s1(k,l)=sigma(m,i1,i2)
	                   s2(l,k)=sigma(m,j1,j2)
	                   s2(k,l)=sigma(m,j1,j2)
	                   m=m+1
	                END DO
	             END DO
             sij=bii*klmnorm(mu(1,i1,i2),mu(1,j1,j2),s1,s2,d,work)
                     IF (sij.gt.spmax) CYCLE
                     wj=wj*dexp(-sij)
                  ENDIF
                  swj=swj+wj
	          DO k=1,ds
                     swjs(k) = swjs(k)+wj*ys(k,j1,j2)
	          END DO
	          DO k=1,d
                     swjm(k) = swjm(k)+wj*y(k,j1,j2)
	          END DO
               END DO
            END DO
	    DO k=1,ds
               asi(k,i1,i2)=swjs(k)
	    END DO
	    DO k=1,d
               ami(k,i1,i2)=swjm(k)
	    END DO
            bi(i1,i2)=swj
            bi0(i1,i2)=swj0
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        one iteration of univariate local constant aws for general Gaussian case
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cvawsntr(y,ys,fix,n1,n2,n3,d,ds,hakt,lambda,mu,sigma,
     1            bi,bi0,ami,asi,kern,spmax,swjs,swjm,work,s1,s2,wght)
C   
C   y        observed values of regression function
C   n        number of observations
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C
      implicit logical (a-z)
      external klmnorm,lkern
      real*8 klmnorm,lkern
      integer n1,n2,n3,kern,d,ds
      logical aws,fix(n1,n2,n3)
      real*8 y(d,n1,n2,n3),ys(ds,n1,n2,n3),hakt,sigma(ds,n1,n2,n3),
     1       mu(d,n1,n2,n3),bi(n1,n2,n3),bi0(n1,n2,n3),
     2       ami(d,n1,n2,n3),lambda,spmax,swjs(ds),swjm(d),work(d,d),
     3       s1(d,d),asi(ds,n1,n2,n3),s2(d,d),wght(2)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,k,l,m,ja1,je1,ja2,je2,
     1        ja3,je3
      real*8 bii,sij,swj,z1,z2,z3,wj,swj0,bii0,hakt2
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               IF (fix(i1,i2,i3)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               bii=bi(i1,i2,i3)/lambda
C   scaling of sij outside the loop
               bii0=bi0(i1,i2,i3)
               ja1=max0(1,i1-ih1)
               je1=min0(n1,i1+ih1)
               swj=0.d0
               swj0=0.d0
               DO k=1,ds
                  swjs(k)=0.d0
               END DO
               DO k=1,d
                  swjm(k)=0.d0
               END DO
               DO j1=ja1,je1
                  z1=(i1-j1)
                  z1=z1*z1
                  ih2=dsqrt(hakt2-z1)/wght(1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
                  DO j2=ja2,je2
C  first stochastic term
                     z2=(i2-j2)*wght(1)
                     z2=z2*z2
                     ih3=dsqrt(hakt2-z1-z2)/wght(2)
                     ja3=max0(1,i3-ih3)
                     je3=min0(n3,i3+ih3)
                     DO j3=ja3,je3
C  first stochastic term
                        z3=(i3-j3)*wght(2)
                        wj=lkern(kern,z1+z2*z2)
                        swj0=swj0+wj
                        IF (aws) THEN
                           m=1
                           DO k=1,d
	                      DO l=1,k
	                         s1(l,k)=sigma(m,i1,i2,i3)
	                         s1(k,l)=sigma(m,i1,i2,i3)
	                         s2(l,k)=sigma(m,j1,j2,j3)
	                         s2(k,l)=sigma(m,j1,j2,j3)
	                         m=m+1
	                      END DO
	                   END DO
          sij=bii*klmnorm(mu(1,i1,i2,i3),mu(1,j1,j2,j3),s1,s2,d,work)
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*dexp(-sij)
                        ENDIF
                        swj=swj+wj
	                DO k=1,ds
                           swjs(k) = swjs(k)+wj*ys(k,j1,j2,j3)
                        END DO
                        DO k=1,d
                           swjm(k) = swjm(k)+wj*y(k,j1,j2,j3)
                        END DO
                     END DO
                  END DO
               END DO
	       DO k=1,ds
                  asi(k,i1,i2,i3)=swjs(k)
	       END DO
	       DO k=1,d
                  ami(k,i1,i2,i3)=swjm(k)
	       END DO
               bi(i1,i2,i3)=swj
               bi0(i1,i2,i3)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
