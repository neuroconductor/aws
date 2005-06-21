CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C   Kullback-Leibler Distance of two general Gaussian distributions 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      real*8 function klnorm(mu1,mu2,s1,s2)
      implicit logical (a-z)
      real*8 s1,s2,mu1,mu2,mud
      mud=mu1-mu2
C      klnorm=dlog(s2/s1)-1+(s1+mud*mud)/s2
      klnorm=dlog(s1/s2)-1+(s2+mud*mud)/s1
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C   Kullback-Leibler Distance of two general Gaussian distributions 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine klnormn(mu1,mu2,s1,s2,n,kld)
      implicit logical (a-z)
      integer n,l
      real*8 s1(n),s2(n),kld(n),mu1(n),mu2(n),mud
      DO l=1,n
         mud=mu1(l)-mu2(l)
C         kld(l)=dlog(s2(l)/s1(l))-1+(s1(l)+mud*mud)/s2(l)
         kld(l)=dlog(s1(l)/s2(l))-1+(s2(l)+mud*mud)/s1(l)
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
      subroutine cawsnuni(y,ys,fix,n,hakt,lambda,mu,sigma,bi,bi0,
     1                    ami,asi,kern,spmax)
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
      external klnorm,lkern
      real*8 klnorm,lkern
      integer n,kern
      logical aws,fix(n)
      real*8 y(n),ys(n),hakt,sigma(n),mu(n),bi(n),bi0(n),
     1       ami(n),lambda,spmax,swjs,swjm,s1,asi(n)
      integer ih,i,j,ja,je
      real*8 bii,sij,swj,z,wj,swj0,bii0,mui
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
         swjs=0.d0
         swjm=0.d0
	 s1=sigma(i)
	 mui=mu(i)
         DO j=ja,je
C  first stochastic term
            z=(i-j)/hakt
            wj=lkern(kern,z*z)
            swj0=swj0+wj
            IF (aws) THEN
               sij=bii*klnorm(mui,mu(j),s1,sigma(j))
               IF (sij.gt.spmax) CYCLE
               wj=wj*dexp(-sij)
            ENDIF
            swj=swj+wj
            swjs = swjs+wj*ys(j)
            swjm = swjm+wj*y(j)
         END DO
         asi(i)=swjs
         ami(i)=swjm
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
      subroutine cawsnbi(y,ys,fix,n1,n2,hakt,lambda,mu,sigma,bi,bi0,
     1                    ami,asi,kern,spmax,wght)
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
      external klnorm,lkern
      real*8 klnorm,lkern
      integer n1,n2,kern
      logical aws,fix(n1,n2)
      real*8 y(n1,n2),ys(n1,n2),hakt,sigma(n1,n2),mu(n1,n2),bi(n1,n2),
     1       bi0(n1,n2),ami(n1,n2),lambda,spmax,swjs,swjm,asi(n1,n2)
      integer ih1,ih2,i1,i2,j1,j2,ja1,je1,ja2,je2
      real*8 bii,sij,swj,z1,z2,wj,swj0,bii0,mui,s1,wght,hakt2
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO  i1=1,n1
         DO  i2=1,n2
            IF (fix(i1,i2)) CYCLE
C    nothing to do, final estimate is already fixed by control 
            bii=bi(i1,i2)/lambda
C   scaling of sij outside the loop
            bii0=bi0(i1,i2)
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            swj=0.d0
            swj0=0.d0
            swjs=0.d0
            swjm=0.d0
	    s1=sigma(i1,i2)
	    mui=mu(i1,i2)
            DO j1=ja1,je1
               z1=(i1-j1)
               z1=z1*z1
               ih2=dsqrt(hakt2-z1)/wght
               ja2=max0(1,i2-ih2)
               je2=min0(n2,i2+ih2)
               DO j2=ja2,je2
C  first stochastic term
                  z2=(i2-j2)*wght
                  wj=lkern(kern,(z1+z2*z2)/hakt2)
                  swj0=swj0+wj
                  IF (aws) THEN
                     sij=bii*klnorm(mui,mu(j1,j2),s1,sigma(j1,j2))
                     IF (sij.gt.spmax) CYCLE
                     wj=wj*dexp(-sij)
                  ENDIF
                  swj=swj+wj
                  swjs = swjs+wj*ys(j1,j2)
                  swjm = swjm+wj*y(j1,j2)
               END DO
            END DO
            asi(i1,i2)=swjs
            ami(i1,i2)=swjm
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
      subroutine cawsntri(y,ys,fix,n1,n2,n3,hakt,lambda,mu,sigma,bi,
     1                    bi0,ami,asi,kern,spmax,wght)
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
      external klnorm,lkern
      real*8 klnorm,lkern
      integer n1,n2,n3,kern
      logical aws,fix(n1,n2,n3)
      real*8 y(n1,n2,n3),ys(n1,n2,n3),hakt,sigma(n1,n2,n3),
     1       mu(n1,n2,n3),bi(n1,n2,n3),bi0(n1,n2,n3),ami(n1,n2,n3),
     2       lambda,spmax,swjs,swjm,asi(n1,n2,n3)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,ja1,je1,ja2,je2,ja3,je3
      real*8 bii,sij,swj,z1,z2,z3,wj,swj0,bii0,mui,s1,wght(2),hakt2
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO  i1=1,n1
         DO  i2=1,n2
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
               swjs=0.d0
               swjm=0.d0
	       s1=sigma(i1,i2,i3)
	       mui=mu(i1,i2,i3)
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
                        wj=lkern(kern,(z1+z2+z3*z3)/hakt2)
                        swj0=swj0+wj
                        IF (aws) THEN
                           sij=bii*klnorm(mui,mu(j1,j2,j3),s1,
     1                                        sigma(j1,j2,j3))
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*dexp(-sij)
                        ENDIF
                        swj=swj+wj
                        swjs = swjs+wj*ys(j1,j2,j3)
                        swjm = swjm+wj*y(j1,j2,j3)
                     END DO
                  END DO
               END DO
               asi(i1,i2,i3)=swjs
               ami(i1,i2,i3)=swjm
               bi(i1,i2,i3)=swj
               bi0(i1,i2,i3)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        one iteration of univariate local constant aws for general Gaussian case
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawsnorm(y,ys,fix,n1,n2,n3,hakt,lambda,mu,sigma,bi,
     1                    bi0,ami,asi,kern,spmax,wght)
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
      external klnorm,lkern
      real*8 klnorm,lkern
      integer n1,n2,n3,kern
      logical aws,fix(1)
      real*8 y(1),ys(1),hakt,sigma(1),mu(1),bi(1),bi0(1),ami(1),
     1       lambda,spmax,swjs,swjm,asi(1)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,ja1,je1,ja2,je2,ja3,je3,
     1       iind,jind,jind3,jind2 
      real*8 bii,sij,swj,z1,z2,z3,wj,swj0,bii0,mui,s1,wght(2),hakt2
      hakt2=hakt*hakt
      ih1=hakt
      aws=lambda.lt.1d40
      DO i3=1,n3
         DO i2=1,n2
             DO i1=1,n1
	       iind=i1+(i2-1)*n1+(i3-1)*n1*n2
               IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               bii=bi(iind)/lambda
C   scaling of sij outside the loop
               bii0=bi0(iind)
	       ih3=hakt/wght(2)
               ja3=max0(1,i3-ih3)
               je3=min0(n3,i3+ih3)
               swj=0.d0
               swj0=0.d0
               swjs=0.d0
               swjm=0.d0
	       s1=sigma(iind)
	       mui=mu(iind)
               DO j3=ja3,je3
                  z3=(i3-j3)*wght(2)
                  z3=z3*z3
                  ih2=dsqrt(hakt2-z3)/wght(1)
                  ja2=max0(1,i2-ih2)
                  je2=min0(n2,i2+ih2)
		  jind3=(j3-1)*n1*n2
                  DO j2=ja2,je2
                     z2=(i2-j2)*wght(1)
                     z2=z3+z2*z2
                     ih1=dsqrt(hakt2-z2)
                     ja1=max0(1,i1-ih1)
                     je1=min0(n1,i1+ih1)
		     jind2=jind3+(j2-1)*n1
                     DO j1=ja1,je1
C  first stochastic term
                        jind=j1+jind2
                        z1=(i1-j1)
                        wj=lkern(kern,(z1*z1+z2)/hakt2)
                        swj0=swj0+wj
                        IF (aws) THEN
                           sij=bii*klnorm(mui,mu(jind),s1,
     1                                        sigma(jind))
                           IF (sij.gt.spmax) CYCLE
                           wj=wj*dexp(-sij)
                        ENDIF
                        swj=swj+wj
                        swjs = swjs+wj*ys(jind)
                        swjm = swjm+wj*y(jind)
                     END DO
                  END DO
               END DO
               asi(iind)=swjs
               ami(iind)=swjm
               bi(iind)=swj
               bi0(iind)=swj0
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cpawsnun(n,dp1,dp2,x,y,fix,theta,sigma,bi,bi0,
     1        ai,lam,h,kern,cb,dmat,thij,psix,psiy,bii,spmax)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     x          design points (ordered)
C     y          observed values at design points
C     theta      old estimates from last step (input)
C     bi         \sum \Psi^T Wi^(k-1) \Psi    (input)
C     si0        \sum \Psi^T Wi0 \Psi [1,]    (input)
C     bi0        \sum \Psi^T Wi0 \Psi    (input/output)
C     ai         \sum \Psi^T Wi^(k-1) Y       (input)
C     ai        \sum \Psi^T Wi^k Y           (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     cb         binomial coefficients     
C     dmat       working arrays  dp1 times dp1
C     thij                     vector for parameter differences
C     psix,psiy                working memory 
C     
      implicit logical (a-z)
      integer n,dp1,dp2,i,j,k,l,je,ja,nwij,kern 
      logical aws,fix(n)
      real*8 lkern,x(n),y(n),psix(dp2),psiy(dp1),theta(dp1,n),
     1 bi(dp2,n),ai(dp1,n),lam,spmax,dmat(dp1,dp1),bii(dp2),
     2 thij(dp1),cb(dp1,dp1),bi0(dp2,n),lambda,h,sigma(n),
     3 sij,wij,z,xij,ha,ha2,eps,xi,zj,thijl,swijs,res
      external lkern,skern
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      eps=1.e-10
      aws=lam.gt.1.d20
      DO i=1,n
C        loop over design points, initialization
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
	 lambda=lam
         ha=h
         xi=x(i)
C     first search for ja and je
1099     DO j=i,1,-1
            IF (x(j).le.xi-ha) EXIT
            ja=j
         END DO
         DO j=i,n
            IF (x(j).gt.xi+ha) EXIT
            je=j
         END DO
         IF (je-ja-1.le.dp1) THEN
            ha=ha*1.25
C            not enough points in neighborhood to estimate parameters
C            increase ha
            goto 1099
         END IF
         ha2=ha*ha
8999     nwij=0
         IF (aws) THEN
C
C            prepare for computation of stochastic and extension penalty 
C
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmat
C    expand Bi in dmat
            l=1
            DO j=1,dp1
               DO k=1,dp1
                  dmat(j,k)=bi(j+k-1,i)/lambda
	       END DO
            END DO
         END IF
         DO l=1,dp2
            bii(l)=0.d0
            bi0(l,i)=0.d0
            IF (l.gt.dp1) CYCLE
            ai(l,i)=0.d0
         END DO
	 swijs=0.d0
C      if not enough points with positive weights (counted in nwij)
C      lambda will be increased 
C      (to reduce panalization by stochastic and influence term)
C
C             loop over local neighborhood
C
         DO j=ja,je
C
C              get location weights
C
            xij=(x(j)-x(i))
            wij=lkern(kern,xij*xij/ha2)
            z=1.d0
            DO k=1,dp2
               psix(k)=z
               IF (k.le.dp1) psiy(k)=z*y(j)
               z=z*xij
            END DO
            DO k=1,dp2
               bi0(k,i)=bi0(k,i)+wij*psix(k)
	    END DO
	    res=y(j)-theta(1,i)
	    DO k=2,dp1
	       res=res-theta(k,i)*psix(k)
	    END DO
            IF (aws) THEN
C    now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
               DO l=1,dp1
                 thij(l)=theta(l,j)
               END DO
               IF (dp1.gt.1) THEN
                  z=1.d0
                  zj=1.d0
                  DO l=2,dp1
                     z=-z*xij
                     DO k=1,dp1-l+1
                        thij(k)=thij(k)+cb(k+l-1,k)*z*theta(k+l-1,j)
                     END DO
                  END DO
               END IF
C     thij contains theta_{j} in the model centered at xi
C
C                 now get sij
C
	       DO l=1,dp1
                  thij(l)=theta(l,i)-thij(l)
               END DO
C
C   thats the difference between thetai and thetaij
C
               sij=0.0d0
               DO l=1,dp1
	          thijl=thij(l)
		  sij=sij+dmat(l,l)*thijl*thijl
		  IF (l.eq.dp1) CYCLE
                  DO k=l+1,dp1
                     sij=sij+2.d0*dmat(k,l)*thijl*thij(k)
                  END DO
               END DO
	       sij=sij/sigma(i)-1+dlog(sigma(i)/sigma(j))
C
C     now we have everything to compute  w_{ij}
C       
               IF (sij.gt.spmax) CYCLE
               wij=wij*dexp(-sij)
            END IF
            IF (wij.gt.eps) nwij=nwij+1
C        now compute contributions to bi(i),bi0(i),ai(i)  
            DO k=1,dp2
               bii(k)=bii(k)+wij*psix(k)
               IF (k.le.dp1) ai(k,i)=ai(k,i)+wij*psiy(k)
            END DO 
	    swijs=swijs+res*res*wij   
         END DO
C    this was the j - loop
         IF (nwij.lt.dp1.and.aws) THEN
            lambda=lambda*1.25
C    increase lambda to weaken stochastic penalization
            goto 8999
	 END IF
	 DO k=1,dp2 
	    bi(k,i)=bii(k)
         END DO	
	 IF (nwij.gt.dp1) sigma(i)=swijs/(bii(1)-dp1)
      END DO      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in bivariate local polynomial aws (gridded) 
C
C   p > 0   only  !!!! 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cpawsnbi(n1,n2,dp1,dp2,y,fix,theta,sigma,bi,bi0,ai,
     1    lam,h,kern,dmat,thij,psix,si,si0,siy,ind,wght,spmax)
C
C     n1         number of points in first dimension
C     n2         number of points in second dimension
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     y          observed values at design points
C     theta      estimates from step (k-1)   (input)
C     bi         \sum \Psi^T Wi \Psi  from step (k-1)
C     bin        \sum \Psi^T Wi \Psi  (output)
C     bi0        \sum \Psi^T Wi0 \Psi  from step (k-1)
C     ai        \sum \Psi^T Wi Y     (output)
C     lam        lambda*sigma2            (stochastic penalty parameter)
C     h          actual bandwidth         (location penalty parameter)
C     dmat       working arrays  dp1 times dp1
C     thij               vector for parameter differences
C     psix                     working memory for Psij
C     si         working array  for \sum \Psi^T Wi \Psi
C     si0        working array  for \sum \Psi^T Wi0 \Psi
C     siy        working array  for \sum \Psi^T Wi Y 
C
      implicit logical (a-z)
      integer n1,n2,dp1,dp2,i1,i2,j1,j2,k,l,info,je1,ja1,je2,ja2,
     1        j,ih1,ih2,ind(dp1,dp1),m,kern
      logical aws,fix(n1,n2)
      real*8 bi(dp2,n1,n2),bi0(dp2,n1,n2),lkern,
     1       ai(dp1,n1,n2),wght,theta(dp1,n1,n2),siy(dp1),
     2       dmat(dp1,dp1),thij(dp1),spmax,y(n1,n2),h,lam,
     3       psix(dp2),si(dp2),si0(dp2),sii,wijy,d,z11,z12,z22,ha2,
     4       lambda,wij,s0i,z1,z2,thijl,sigma(n1,n2),swijs,res,
     5       sigmai,thi1,thi2,thi3,thi4,thi5,thi6
      external lkern
C     
C     in case of dp1==1  lawsbi should be called  (p=0)
C
      ha2=h*h
      ih1=h
      aws=lam.lt.1.d20
      DO i1=1,n1
         DO i2=1,n2
	    IF (fix(i1,i2)) CYCLE
C    nothing to do, estimate in (i1,i2) is already fixed by control
            lambda=lam
	    sigmai=sigma(i1,i2)
	    thi1=theta(1,i1,i2)
	    thi2=theta(2,i1,i2)
	    thi3=theta(3,i1,i2)
	    IF (dp1.gt.3) THEN
	       thi4=theta(4,i1,i2)
	       thi5=theta(5,i1,i2)
	       thi6=theta(6,i1,i2)
	    END IF
C  first fill si and siy with 0's
C  fields are used to sum components of ai, bin and bi0
8999        DO j=1,dp2
               si(j)=0.d0
               si0(j)=0.d0
               IF (j.le.dp1) siy(j)=0.d0
            END DO
            IF (aws) THEN
               DO j=1,dp1
                  DO k=1,dp1
                     m=ind(j,k)
                     dmat(j,k)=bi(m,i1,i2)/lambda
                  END DO
               END DO
               s0i=bi0(1,i1,i2)
               sii=bi(1,i1,i2)
            END IF
	    swijs=0.d0
C
C     Prepare for loop over j's
C
            ja1=max0(1,i1-ih1)
            je1=min0(n1,i1+ih1)
            DO j1=ja1,je1
               z1=(i1-j1)
               z11=z1*z1
               ih2=dsqrt(ha2-z11)/wght
               ja2=max0(1,i2-ih2)
               je2=min0(n2,i2+ih2)
               DO j2=ja2,je2
                  z2=(i2-j2)*wght
                  z22=z2*z2
                  z12=z1*z2
C          first compute location part
                  wij=lkern(kern,(z11+z22)/ha2)
                  si0(1)=si0(1)+wij
                  si0(2)=si0(2)-z1*wij
                  si0(3)=si0(3)-z2*wij
                  si0(4)=si0(4)+z11*wij
                  si0(5)=si0(5)+z12*wij
                  si0(6)=si0(6)+z22*wij
		  res=y(j1,j2)-thi1+z1*thi2+z2*thi3
                  IF (dp1.gt.3) THEN
                     si0(7)=si0(7)-z11*z1*wij
                     si0(8)=si0(8)-z11*z2*wij
                     si0(9)=si0(9)-z1*z22*wij
                     si0(10)=si0(10)-z2*z22*wij
                     si0(11)=si0(11)+z11*z11*wij
                     si0(12)=si0(12)+z11*z12*wij
                     si0(13)=si0(13)+z11*z22*wij
                     si0(14)=si0(14)+z12*z22*wij
                     si0(15)=si0(15)+z22*z22*wij
		     res=res-z11*thi4-z11*thi5-z22*thi6
		  END IF
C          this is the location penalty, now fill si0
                  IF (aws) THEN 
C          now fill psix 
C           now translate thetaj into model centered in xi
C
                     thij(1)=theta(1,j1,j2)
                     thij(2)=theta(2,j1,j2)
                     thij(3)=theta(3,j1,j2)
                     thij(1)=thij(1)+thij(2)*z1+thij(3)*z2
                     IF (dp1.gt.3) THEN
                        thij(4)=theta(4,j1,j2)
                        thij(5)=theta(5,j1,j2)
                        thij(6)=theta(6,j1,j2)
                        thij(1)=thij(1)+thij(4)*z11+
     +                           thij(5)*z12+thij(6)*z22
                        thij(2)=thij(2)+thij(5)*z2+2.d0*thij(4)*z1
                        thij(3)=thij(3)+thij(5)*z1+2.d0*thij(6)*z2
		     END IF
C  
C           get difference of thetas
C
                     DO l=1,dp1
                        thij(l)=theta(l,i1,i2)-thij(l)
                     END DO
C
C           get stochastic penalty
C
                     d=0.d0
                     DO l=1,dp1
		        thijl=thij(l)
		        d=d+dmat(l,l)*thijl*thijl
			IF(l.eq.dp1) CYCLE
                        DO k=l+1,dp1
                           d=d+2.d0*dmat(k,l)*thijl*thij(k)
                        END DO
                     END DO
		     d=d/sigma(i1,i2)-1+dlog(sigma(i1,i2)/sigma(j1,j2))
C
C           now compute weights 
C           stochastic and extension penalty can be added because kernel is exp
C
                     IF (d.gt.spmax) CYCLE
                     wij=wij*dexp(-d)
	          END IF
C           now compute contributions to bi(i),ai(i)  
                  wijy=y(j1,j2)*wij
	          swijs=swijs+res*res*wij   
                  si(1)=si(1)+wij
                  si(2)=si(2)-z1*wij
                  si(3)=si(3)-z2*wij
                  si(4)=si(4)+z11*wij
                  si(5)=si(5)+z12*wij
                  si(6)=si(6)+z22*wij
                  siy(1)=siy(1)+wijy
                  siy(2)=siy(2)-z1*wijy
                  siy(3)=siy(3)-z2*wijy
        	  IF (dp1.le.3) CYCLE
                  si(7)=si(7)-z11*z1*wij
                  si(8)=si(8)-z11*z2*wij
                  si(9)=si(9)-z1*z22*wij
                  si(10)=si(10)-z2*z22*wij
                  si(11)=si(11)+z11*z11*wij
                  si(12)=si(12)+z11*z12*wij
                  si(13)=si(13)+z11*z22*wij
                  si(14)=si(14)+z12*z22*wij
                  si(15)=si(15)+z22*z22*wij
                  siy(4)=siy(4)+z11*wijy
                  siy(5)=siy(5)+z12*wijy
		  siy(6)=siy(6)+z22*wijy
               END DO
            END DO
C        prepare matrix to test for singularity of Bi 
C        this should be changed to SVD at some point
            DO k=1,dp1
               DO j=1,dp1
                  IF (j.gt.k) THEN
                     dmat(j,k)=0.d0
                  ELSE
                     dmat(j,k)=si(ind(j,k))
                  ENDIF
               END DO
            END DO
C
C     compute choleski decomposition
C
            call invers(dmat,dp1,info)
            IF (info.gt.0) THEN
C
C          if singular relax stochastic and extension penalty 
C
               lambda=1.5*lambda
               goto 8999
            END IF
C     
C     now fill ai, bi and bi0
C
            DO j=1,dp1
               ai(j,i1,i2)=siy(j)
            END DO
            DO j=1,dp2
               bi(j,i1,i2)=si(j)
               bi0(j,i1,i2)=si0(j)
            END DO
	 IF (si(1).gt.dp1) sigma(i1,i2)=swijs/(si(1)-dp1)
         END DO     
      END DO     
      RETURN
      END
