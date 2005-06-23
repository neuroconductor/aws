CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C   Compute Inverse and determinant of symmetric positive definite a
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine invdeta(a,n,ainv,det,info)
      implicit logical (a-z)
      integer n,info,i,j,job
      real*8 a(n,n),det,det0(2),ainv(n,n)
      job=11
      DO i=1,n
         DO j=1,n
            ainv(i,j)=a(i,j)
         END DO
      END DO
      call dpofa(ainv,n,n,info)
      if(info.gt.0) THEN
      END IF
      call dpodi(ainv,n,n,det0,job)
      det=det0(1)*10.0**det0(2)
      DO i=1,n
         if(i.eq.n) CYCLE
         DO j=i+1,n
            ainv(j,i)=ainv(i,j)
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C   Compute the determinant of symmetric positive definite a
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine deta(a,n,ainv,det,info)
      implicit logical (a-z)
      integer n,info,i,j,job
      real*8 a(n,n),det,det0(2),ainv(n,n)
      job=10
      DO i=1,n
         DO j=1,n
            ainv(i,j)=a(i,j)
         END DO
      END DO
      call dpofa(ainv,n,n,info)
      if(info.gt.0) THEN
      END IF
      call dpodi(ainv,n,n,det0,job)
      det=det0(1)*10.0**det0(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C   Kullback-Leibler Distance of volatility matrices
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C      real*8 function klvola(s1,s2,d,s2inv)
      real*8 function klvola(s2,s1,d,s2inv)
      implicit logical (a-z)
      integer d,i,j,info
      real*8 s1(d,d),s2(d,d),s2inv(d,d)
      real*8 z,det2,det1
      call deta(s1,d,s2inv,det1,info)
      if(info.gt.0) THEN
         klvola=0.d0
	 RETURN
      END IF
      call invdeta(s2,d,s2inv,det2,info)
      if(info.gt.0) THEN
         klvola=0.d0
	 RETURN
      END IF
      z=det2/det1
      z=dlog(z)-d
      do i=1,d
         do j=1,d
            z=z+s1(i,j)*s2inv(i,j)
         END DO
      END DO
      klvola=z
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C   Kullback-Leibler Distance of volatility matrices
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C      subroutine klvolan(th1,th2,d,ds,n,s1,s2,s2inv,kld)
      subroutine klvolan(th1,th2,d,ds,n,s2,s1,s2inv,kld)
      implicit logical (a-z)
      integer d,ds,n,i,j,l,m,info
      real*8 th1(ds,n),th2(ds,n),s1(d,d),s2(d,d),s2inv(d,d),z,
     1       det2,det1,kld(n)
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
	    CYCLE
	 END IF
         call invdeta(s2,d,s2inv,det2,info)
	 IF(info.gt.0) THEN
	    kld(l)=0
	    CYCLE
	 END IF
         z=dlog(det2/det1)-d
         do i=1,d
            do j=1,d
               z=z+s1(i,j)*s2inv(i,j)
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
C        one iteration of univariate local constant aws on a grid
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cvawsvol(y,fix,n,d,ds,hakt,lambda,theta,bi,bi0,ai,
     1                    kern,spmax,swjy,work,s1,s2)
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
      external klvola,lkern
      real*8 klvola,lkern
      integer n,kern,d,ds
      logical aws,fix(n)
      real*8 y(ds,n),hakt,theta(ds,n),bi(n),bi0(n),ai(ds,n),lambda,
     1       spmax,swjy(ds),work(d,d),s1(d,d),s2(d,d)
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
            swjy(k)=0.d0
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
	             s1(l,k)=theta(m,i)
	             s1(k,l)=theta(m,i)
	             s2(l,k)=theta(m,j)
	             s2(k,l)=theta(m,j)
	             m=m+1
	          END DO
	       END DO
               sij=bii*klvola(s1,s2,d,work)
               IF (sij.gt.spmax) CYCLE
               wj=wj*dexp(-sij)
            ENDIF
            swj=swj+wj
	    DO k=1,ds
               swjy(k) = swjy(k)+wj*y(k,j)
	    END DO
         END DO
	 DO k=1,ds
           ai(k,i)=swjy(k)
	 END DO
         bi(i)=swj
         bi0(i)=swj0
      END DO
      RETURN
      END
