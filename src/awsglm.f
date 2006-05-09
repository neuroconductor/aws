CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FORTRAN CODE related to awsglm.r
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine glawsunw(n,dp1,dp2,dp3,y,fix,mfamily,theta0,theta,
     1                    bii,bi,bi2,bi0,bi02,ai,ni,lam,h,hw,kern,cb,
     2                    dmat,thij,psix,wghts,wghts0,spmax,iter,
     3                    work)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
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
      integer n,dp1,dp2,i,j,k,l,je,ja,kern,mfamily,ih,iter,
     1        dp3,info
      logical aws,fix(n)
      real*8 y(n),psix(dp2),theta(dp1,n),bi2(dp2,n),bi(dp2,n),
     1     ai(dp1),lam,spmax,dmat(dp1,dp1),lkern,thij(dp1),cb(dp1,dp1),
     2     bi0(dp2,n),lambda,h,ni(n),ha,ha2,bi02(dp2,n),wghts(n,n),
     3     wghts0(n),theta0(dp1,n),bii(dp3,n),d,work(n),hw,bcorr
      external lkern
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      aws=lam.lt.1.d20
C
C   Regularization for Poisson and Bernoulli
C
      IF(mfamily.eq.3) THEN
         bcorr=0.5d0/h
         DO i=1,n
	    if(y(i).gt.0.5d0) THEN
	        y(i)=y(i)-bcorr
	    ELSE
	        y(i)=y(i)+bcorr
	    END IF
	 END DO
      ELSE IF(mfamily.eq.2) THEN
         bcorr=0.5d0/h
         DO i=1,n
	    y(i)=y(i)+bcorr
	 END DO
      END IF
      DO i=1,n
C        loop over design points, initialization
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
	 lambda=lam
         ih=h
         ja=max0(1,i-ih)
         je=min0(n,i+ih)
         ha=h
         ha2=ha*ha
         IF (aws) THEN
C
C            prepare for computation of stochastic and extension penalty 
C
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmat
C    expand Bi in dmat
            l=1
            DO j=1,dp1
               DO k=1,j
                  d=bii(l,i)/lambda
                  dmat(k,j)=d
                  dmat(j,k)=d
                  l=l+1
	       END DO
            END DO
         END IF
C
C          Now get weights in wghts and wghts0
C
         call gluniwgt(n,i,ja,je,dp1,dp2,ha2,aws,kern,
     1              theta0,spmax,psix,thij,dmat,cb,wghts(1,i),wghts0)
            call smwghtgl(n,i,ja,je,wghts(1,i),work,ha,hw,kern)
         ni(i)=0.d0
         DO j=ja,je
            ni(i)=ni(i)+wghts(j,i)
         END DO
C
C          Now do the iterations to obtain the new estimates
C
         call gluniit(mfamily,n,dp1,dp2,i,ja,je,y,theta(1,i),psix,
     1                dmat,wghts(1,i),wghts0,thij,bi(1,i),ai,bi0(1,i),
     2                bi2(1,i),bi02(1,i),iter,info)
         if(info.gt.0) fix(i)=.TRUE.
C   singularity or failed convergence keep old estimate in this case
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine glawsuni(n,dp1,dp2,dp3,y,fix,mfamily,theta0,theta,
     1                    bii,bi,bi2,bi0,bi02,ai,ni,lam,h,hw,kern,cb,
     2                    dmat,thij,psix,wghts,wghts0,spmax,iter,
     3                    work)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
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
      integer n,dp1,dp2,i,j,k,l,je,ja,kern,mfamily,ih,iter,
     1        dp3,info
      logical aws,fix(n)
      real*8 y(n),psix(dp2),theta(dp1,n),bi2(dp2,n),bi(dp2,n),
     1     ai(dp1),lam,spmax,dmat(dp1,dp1),lkern,thij(dp1),cb(dp1,dp1),
     2     bi0(dp2,n),lambda,h,ni(n),ha,ha2,bi02(dp2,n),wghts(n),
     3     wghts0(n),theta0(dp1,n),bii(dp3,n),d,work(n),hw,bcorr
      external lkern
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      aws=lam.lt.1.d20
C
C   Regularization for Poisson and Bernoulli
C
      IF(mfamily.eq.3) THEN
         bcorr=0.5d0/h
         DO i=1,n
	    if(y(i).gt.0.5d0) THEN
	        y(i)=y(i)-bcorr
	    ELSE
	        y(i)=y(i)+bcorr
	    END IF
	 END DO
      ELSE IF(mfamily.eq.2) THEN
         bcorr=0.5d0/h
         DO i=1,n
	    y(i)=y(i)+bcorr
	 END DO
      END IF
      DO i=1,n
C        loop over design points, initialization
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
	 lambda=lam
         ih=h
         ja=max0(1,i-ih)
         je=min0(n,i+ih)
         ha=h
         ha2=ha*ha
         IF (aws) THEN
C
C            prepare for computation of stochastic and extension penalty 
C
C    get (Psi W_i Psi^T) and (Psi \tilde{W}_i Psi^T) in dmat
C    expand Bi in dmat
            l=1
            DO j=1,dp1
               DO k=1,j
                  d=bii(l,i)/lambda
                  dmat(k,j)=d
                  dmat(j,k)=d
                  l=l+1
	       END DO
            END DO
         END IF
C
C          Now get weights in wghts and wghts0
C
         call gluniwgt(n,i,ja,je,dp1,dp2,ha2,aws,kern,
     1                theta0,spmax,psix,thij,dmat,cb,wghts,wghts0)
            call smwghtgl(n,i,ja,je,wghts,work,ha,hw,kern)
         ni(i)=0.d0
         DO j=ja,je
            ni(i)=ni(i)+wghts(j)
         END DO
C
C          Now do the iterations to obtain the new estimates
C
         call gluniit(mfamily,n,dp1,dp2,i,ja,je,y,theta(1,i),psix,
     1                dmat,wghts,wghts0,thij,bi(1,i),ai,bi0(1,i),
     2                bi2(1,i),bi02(1,i),iter,info)
         if(info.gt.0) fix(i)=.TRUE.
C   singularity or failed convergence keep old estimate in this case
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    smooth univariate weights
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smwghtgl(n,i,ja,je,wghts,work,hakt,hw,kern)
      integer n,kern,i,ja,je
      real*8 wghts(n),work(n),hw,hakt
      integer j,k,ka,ke,iha,ihe,ih
      real*8 maxwght,lkern,d,d2,cc,ww
      external lkern
      d=0.d0
      DO k=ja,je
         d=d+wghts(k)
      END DO
      d=(2.d0-d/2*d0)*hw-1+d/2.d0
      d=dmax1(.1d0,dmin1(d,hw))
      cc=dmin1(d-1.d0,1.d0)
      IF(cc.gt.0.d0) THEN
         d2=d*d
         ih=d
	 iha=min0(ja-1,ih)
	 ihe=min0(n-je,ih)
	 ja=ja-iha
	 je=je+ihe
	 DO j=ja,je
	    work(j)=0.d0
	 END DO
         DO j=ja+iha,je-ihe
	    DO k=-iha,ihe
	       ww=1.d0-k*k/d2
	       if(ww.lt.1.d0) ww=cc*ww
               work(j+k)=work(j+k)+wghts(j)*ww
	    END DO
         END DO
	 maxwght=work(i)
C   scale such that wghts(i)=1
         DO j=ja,je
            wghts(j)=work(j)/maxwght
         END DO
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    one iteration for local GLM estimates
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gluniit(mfamily,n,dp1,dp2,i,ja,je,y,theta,psix,
     1           dmat,wghts,wghts0,d,bi,ai,bi0,bi2,bi02,iter,info)
      implicit logical (a-z)
      integer n,dp1,dp2,mfamily,i,ja,je,iter,info
      real*8 y(n),theta(dp1),psix(dp2),bi(dp2),ai(dp1),bi0(dp2),
     1       bi2(dp2),bi02(dp2),d(dp1),wghts(n),wghts0(n),
     2       dmat(dp1,dp1)
      integer it,j,k,l
      real*8 z,dist,xij
      DO it=1,iter
         DO l=1,dp2
            bi(l)=0.d0
            bi0(l)=0.d0
            bi2(l)=0.d0
            bi02(l)=0.d0
            IF (l.gt.dp1) CYCLE
            ai(l)=0.d0
         END DO
         DO j=ja,je
            xij=(j-i)
            z=1.d0
            DO k=1,dp2
               psix(k)=z
               z=z*xij
            END DO
C        now compute contributions to bi(i),bi0(i),ai(i)  
            call guniaibi(mfamily,dp1,dp2,psix,theta,y(j),bi,
     1           ai,bi0,bi2,bi02,wghts(j),wghts0(j))
         END DO
         call mpawsun0(dp1,dp2,ai,bi,dmat,d,info)
         IF(info.gt.0) goto 999
         dist=0.d0
         DO k=1,dp1
            DO l=1,dp1
               dist=dist+d(k)*bi(k+l-1)*d(l)
            END DO
         END DO
         DO k=1,dp1
            theta(k)=theta(k)+d(k)
         END DO
         IF(dist.lt.1.d-12) goto 999
         info=10
      END DO
999   RETURN
      END
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    compute weights for local GLM estimates
C
C   fills arrays of wghts  ! only inices ja:je are affected 
C
C   wghts contains adaptive and wghts0 nonadaptive wghts
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gluniwgt(n,i,ja,je,dp1,dp2,ha2,aws,kern,theta,
     1                    spmax,psix,thij,dmat,cb,wghts,wghts0)
      implicit logical (a-z)
      integer n,i,ja,je,dp2,dp1,kern
      real*8 ha2,wghts(n),wghts0(n),psix(dp2),theta(dp1,n),
     1       thij(dp1),dmat(dp1,dp1),cb(dp1,dp1),spmax
      logical aws
      integer j,k,l
      real*8 xij,z,sij,thijl,wij,wij0,lkern
      external lkern
      DO j=ja,je
         xij=(j-i)
         wij0=lkern(kern,xij*xij/ha2)
         wghts0(j)=wij0
         z=1.d0
         DO k=1,dp2
            psix(k)=z
            z=z*xij
         END DO
         IF (aws) THEN
C    now stochastic penalty for alpha=Inf
C    need translation of theta(l,j) to model centered in xi
            DO l=1,dp1
               thij(l)=theta(l,j)
            END DO
            IF (dp1.gt.1) THEN
               z=1.d0
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
C
C     now we have everything to compute  w_{ij}
C       
            IF (sij.le.spmax) THEN
               wij=wij0*dexp(-sij)
               wghts(j)=wij
            ELSE
               wghts(j)=0.d0
            ENDIF
         ELSE
            wghts(j)=wij0
         END IF
      END DO
      RETURN
      END     
