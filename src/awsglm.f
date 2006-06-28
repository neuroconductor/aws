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
C     dp3        number of components in bii (p+2)*(p+1)/2
C     y          observed values at design points
C     fix        index of that are alredy fixed (no convergence or memory step)
C     theta0     old estimates from last step (input)
C     theta      new estimates from last step (output)
C     bii        old variances bi %*% bi2^{-1} %*% bi (usevar = TRUE) or bi (usevar = FALSE)   (input) 
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
     1        dp3,info,skern
      logical aws,fix(n)
      real*8 y(n),psix(dp2),theta(dp1,n),bi2(dp2,n),bi(dp2,n),
     1     ai(dp1),lam,spmax,dmat(dp1,dp1),thij(dp1),cb(dp1,dp1),
     2     bi0(dp2,n),lambda,h,ni(n),ha,ha2,bi02(dp2,n),wghts(n,n),
     3     wghts0(n),theta0(dp1,n),bii(dp3,n),d,work(n),hw,bcorr
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      aws=lam.lt.1.d20
      skern=1
C
C   Regularization for Poisson and Bernoulli
C
      IF(mfamily.eq.3) THEN
         bcorr=0.25d0*dp1/dsqrt(h)
         DO i=1,n
	    if(y(i).gt.0.5d0) THEN
	        y(i)=y(i)-bcorr
	    ELSE
	        y(i)=y(i)+bcorr
	    END IF
	 END DO
      ELSE IF(mfamily.eq.2) THEN
         bcorr=1d0/dsqrt(h)
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
         call gluniwgt(n,i,ja,je,dp1,dp2,ha2,aws,kern,skern,
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
         if(info.gt.0) THEN
            fix(i)=.TRUE.
            DO k=1,dp1
               theta(k,i)=theta0(k,i)
            END DO
         END IF
C   singularity or failed convergence keep old estimate in this case
         call rchkusr()
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine glawsuni(n,dp1,dp2,dp3,y,fix,mfamily,skern,theta0,
     1                    theta,bii,bi,bi2,bi0,bi02,ai,ni,lam,h,hw,
     2                    kern,cb,dmat,thij,psix,wghts,wghts0,spmax,
     3                    iter,work)
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
     1        dp3,info,ja0,je0,counts,skern
      logical aws,fix(n)
      real*8 y(n),psix(dp2),theta(dp1,n),bi2(dp2,n),bi(dp2,n),
     1     ai(dp1),lam,spmax,dmat(dp1,dp1),thij(dp1),cb(dp1,dp1),
     2     bi0(dp2,n),lambda,h,ni(n),ha,ha2,bi02(dp2,n),wghts(n),
     3     wghts0(n),theta0(dp1,n),bii(dp3,n),d,work(n),hw,bcorr
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      aws=lam.lt.1.d20
C
C   Regularization for Poisson and Bernoulli
C
      IF(mfamily.eq.3) THEN
C         bcorr=0.5d0/dexp(dlog(h)/dp1)
         bcorr=0.25d0*dp1/dsqrt(h)
         DO i=1,n
	    if(y(i).gt.0.5d0) THEN
	        y(i)=y(i)-bcorr
	    ELSE
	        y(i)=y(i)+bcorr
	    END IF
	 END DO
      ELSE IF(mfamily.eq.2) THEN
C         bcorr=1.d0/h
         bcorr=dmin1(0.1d0,dp1/h)
         DO i=1,n
	    if(y(i).eq.0) y(i)=y(i)+bcorr
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
         ja0=ja
         je0=je
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
         call gluniwgt(n,i,ja,je,dp1,dp2,ha2,aws,kern,skern,
     1                theta0,spmax,psix,thij,dmat,cb,wghts,wghts0)
         call smwghtgl(n,i,ja,je,wghts,work,ha,hw,kern)
         call smwghtgl(n,i,ja0,je0,wghts0,work,ha,hw,kern)
         ni(i)=0.d0
         DO j=ja,je
            ni(i)=ni(i)+wghts(j)
         END DO
C
C          Now do the iterations to obtain the new estimates
C
         counts=10
9000     call gluniit(mfamily,n,dp1,dp2,i,ja,je,y,theta(1,i),psix,
     1                dmat,wghts,wghts0,thij,bi(1,i),ai,bi0(1,i),
     2                bi2(1,i),bi02(1,i),iter,info)
C         if(info.gt.0) fix(i)=.TRUE.
         if(info.gt.0) THEN 
C            call intpr("Regularization in i=",20,i,1)
            DO j=ja,je
               wghts(j)=.75d0*wghts(j)+.25d0*wghts0(j)
            END DO
            theta(1,i)=theta0(1,i)
            DO k=2,dp1
               theta(k,i)=0.d0
            END DO
            counts=counts-1
C            call intpr("failed in",9,i,1)
            if(counts.le.0) THEN
               fix(i)=.TRUE.
               DO k=1,dp1
                  theta(k,i)=theta0(k,i)
               END DO
            ELSE
               goto 9000
            END IF
         END IF
         
C   singularity or failed convergence keep old estimate in this case
         call rchkusr()
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
      integer j,k,ka,ke,iha,jan,jen,ih,jpk
      real*8 maxwght,d,d2,cc,ww
      d=0.d0
      DO k=ja,je
         d=d+wghts(k)
      END DO
C      d=(hw-d)*hw
C      d=dmax1(.1d0,dmin1(d,hw))
C      cc=dmin1(d-1.d0,1.d0/dsqrt(hakt))
C      cc=dmin1(d-1.d0,1.d0)
      d=dmax1((hw-d),0.d0)
      cc=dmin1(d/hw,1.d0)   
      IF(cc.gt.0.d0) THEN
         d2=d*d
         ih=d-1.d0
	 jan=max0(1,ja-ih)
	 jen=min0(je+ih,n)
	 DO j=jan,jen
	    work(j)=0.d0
	 END DO
         DO j=ja,je
	    DO k=-ih,ih
	       jpk=j+k
	       if(jpk.le.0.or.jpk.gt.n) CYCLE
	       ww=1.d0-k*k/d2
	       if(ww.lt.1.d0) ww=cc*ww
               work(j+k)=work(j+k)+wghts(j)*ww
	    END DO
         END DO
	 maxwght=work(i)
C   scale such that wghts(i)=1
         DO j=jan,jen
            wghts(j)=dmin1(work(j)/maxwght,1.d0)
         END DO
	 ja=jan
	 je=jen
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
      real*8 z,dist,xij,dth
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
     1           ai,bi0,bi2,bi02,wghts(j),wghts0(j),info)
            IF(info.gt.0) goto 999
         END DO
         call rchkusr()
         call mpawsun0(dp1,dp2,ai,bi,dmat,d,info)
         IF(info.gt.0) goto 999
         dist=0.d0
	 dth=0.d0
         DO k=1,dp1
            DO l=1,dp1
               dist=dist+d(k)*bi(k+l-1)*d(l)
               dth=dth+theta(k)*bi(k+l-1)*theta(l)
            END DO
         END DO
	 dist=dist/(dth+1.d-6)
         DO k=1,dp1
            theta(k)=theta(k)+d(k)
         END DO
C         IF(dist.lt.1.d-12) goto 999
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
      subroutine gluniwgt(n,i,ja,je,dp1,dp2,ha2,aws,kern,skern,theta,
     1                    spmax,psix,thij,dmat,cb,wghts,wghts0)
      implicit logical (a-z)
      integer n,i,ja,je,dp2,dp1,kern,skern
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
                     thij(k)=thij(k)+cb(k+l-1,l)*z*theta(k+l-1,j)
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
               if(skern.eq.1) then
               wij=wij0*dexp(-sij)
               else
               wij=wij0*(1.d0-sij/spmax)
               endif
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Calculate contribution of Y_j to ai and bi in univariate local polynomial aws (GLM)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine guniaibi(mfamily,dp1,dp2,psix,theta,y,bii,ai,
     1                    bi0,bi2,bi02,wij,wij0,info)
      implicit logical (a-z)
      integer mfamily,dp1,dp2,info
      real*8 psix(dp1),theta(dp1),bii(dp2),ai(dp1),bi0(dp2),wij,wij0,
     1       y,ymz1,z,bi2(dp2),bi02(dp2),wij2,wij02
      integer i
      real*8 psith,z1,z2
      info=0
      psith=0.0d0
      wij2=wij*wij
      wij02=wij0*wij0
      DO i=1,dp1
         psith=psith+psix(i)*theta(i)
      ENDDO
      IF (mfamily.eq.1) THEN
C  Gaussian case
         DO i=1,dp2
            bii(i)=bii(i)+psix(i)*wij
            bi0(i)=bi0(i)+psix(i)*wij0
            bi2(i)=bi2(i)+psix(i)*wij2
            bi02(i)=bi02(i)+psix(i)*wij02
            IF (i.le.dp1) ai(i)=ai(i)+(y-psith)*psix(i)*wij
         ENDDO 
      ELSE IF (mfamily.eq.2) THEN
C  Poisson case
         IF (dabs(psith).gt.1.d2) THEN
	    info=1 
	    return
	 END IF
         IF (psith.gt.5.d0) THEN 
             z1=psith-5.d0
             z1=148.4132*(1.d0+z1*(1.d0+.5d0*z1*(1.d0+z1/3.d0)))
         ELSE  
             z1=dexp(psith)
         ENDIF
         IF(wij.gt.0.d0) ymz1=y-z1
         DO i=1,dp2
            z=psix(i)*z1
            bi0(i)=bi0(i)+z*wij0
            bi02(i)=bi02(i)+z*wij02
            if(wij.le.0.d0) CYCLE
            bii(i)=bii(i)+z*wij
            bi2(i)=bi2(i)+z*wij2
            IF (i.le.dp1) ai(i)=ai(i)+ymz1*psix(i)*wij
         ENDDO 
      ELSE IF (mfamily.eq.3) THEN
C  Bernoulli case
         IF (dabs(psith).gt.3.6d1) THEN
	    info=1 
	    return
	 END IF
C         IF (psith.gt.1.d1) THEN 
C             z1=psith-1.d1
C             z2=22026.47d0*(1.d0+z1*(1.d0+.5d0*z1*(1.d0+z1/3.d0)))
C         ELSE  IF (psith.lt.-1.d1) THEN 
C             z1=1.d1-psith
C            z1=1.d0/22026.47d0/(1.d0+z1*(1.d0+.5d0*z1*(1.d0+z1/3.d0)))
C         ELSE 
            z2=dexp(-psith)
C	 END IF
         z1=1.d0/(1.d0+z2)
         z2=z1/(1.d0+1.d0/z2)
         DO i=1,dp2
            bii(i)=bii(i)+z2*psix(i)*wij
            bi0(i)=bi0(i)+z2*psix(i)*wij0
            IF (i.le.dp1) ai(i)=ai(i)+(y-z1)*psix(i)*wij
         ENDDO 
      ELSE IF (mfamily.eq.4) THEN
C  Exponential case
         z1=-1/psith
         z2=-z1/psith
         DO i=1,dp2
            bii(i)=bii(i)+z1*psix(i)*wij
            bi0(i)=bi0(i)+z1*psix(i)*wij0
            IF (i.le.dp1) ai(i)=ai(i)+(y-z2)*psix(i)*wij
         ENDDO 
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (univariate polynomial aws)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mpawsun0(dp1,dp2,ai,bi,dmat,d,info)
C    
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     dpm        number of components in di  (dp1+1)*dp1/2
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi
C     di         inverse of bi     
C     di0         inverse of bi0     
C     theta      new parameter estimate
C     dmat       working array
C
C      implicit logical(a-z)
      integer dp1,dp2,info
      real*8 ai(dp1),bi(dp2),dmat(dp1,dp1),d(dp1)
      integer j,k,eins
      eins=1
      DO k=1,dp1
	 DO j=k,dp1
	    dmat(k,j)=bi(j+k-1)
	 END DO
	 d(k)=ai(k)
      END DO
      call dposv("U",dp1,eins,dmat,dp1,d,dp1,info)
C     just keep the old estimate if info > 0
      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate  B_i Bsq_i^{-1} B_i
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bibi2ibi(n,dp1,dp2,dp3,bi,bi2,erg,dmat,dmat2)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     dp3        number of components in erg  (dp1+1)*dp1/2
C     bi         \sum \Psi^T Wi        
C     bi2        \sum \Psi^T Wi^2        
C     bii        factor for KL-distance
C     erg       
C     dmat       working array
C
      implicit logical(a-z)
      integer n,dp1,dp2,dp3
      real*8 bi(dp2,n),bi2(dp2,n),erg(dp3,n),dmat(dp1,dp1),
     1       dmat2(dp1,dp1)
      integer ii,i,j,k,l,m,info
      real*8 d
      DO ii=1,n
         DO k=1,dp1
	    DO j=k,dp1
	       dmat(k,j)=bi2(j+k-1,ii)
	    END DO
	 END DO
	 call dpotri("U",dp1,dmat,dp1,info)
         IF (info.ne.0) CYCLE 
         m=1
         DO i=1,dp1
            DO j=1,i
               d=0.d0
               DO k=1,dp1
                  DO l=1,dp1
                     d=d+dmat(k,l)*bi(i+k-1,ii)*bi(j+l-1,ii)
                  END DO
               END DO
               erg(m,ii)=d
               m=m+1
            END DO
         END DO
      END DO
      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (univariate polynomial aws)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine confuni(n,dp1,dp2,theta,bi,bi2,conf,dmat)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     dpm        number of components in di  (dp1+1)*dp1/2
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi
C     di         inverse of bi     
C     di0         inverse of bi0     
C     theta      new parameter estimate
C     dmat       working array
C
C      implicit logical(a-z)
      integer n,dp1,dp2
      real*8 bi(dp2,n),bi2(dp2,n),theta(n),dmat(dp1,dp1),conf(2,n)
      integer i,j,k,info
      real*8 d
      DO i=1,n
         DO k=1,dp1
	    DO j=k,dp1
	       dmat(k,j)=bi(j+k-1,i)
	    END DO
	 END DO
	 call dpotri("U",dp1,dmat,dp1,info)
         IF (info.ne.0) CYCLE 
         d=0.d0
         DO j=1,dp1
            DO k=1,dp1
               d=d+dmat(1,k)*dmat(1,j)*bi2(j+k-1,i)
            END DO
         END DO
         d=dsqrt(d)*1.96
         conf(1,i)=theta(i)-d
         conf(2,i)=theta(i)+d
      END DO
      RETURN
      END      
