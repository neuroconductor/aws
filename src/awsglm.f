CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FORTRAN CODE related to awsglm.r
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        search for sufficiently long interval to fit a GLM for Poisson
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine poisdes(y,n,p,h)
      integer n,p,h(2,n)
      real*8 y(n)
      integer i,j,k
      real*8 z
      DO i=1,n
         k=0
         h(1,i)=i
         h(2,i)=i
         IF(y(i).gt.0.d0) k=k+1
         DO j=1,n
            IF(k.gt.p) CYCLE
            z=0.d0
            IF(i.gt.j) z=y(i-j)
            IF(z.gt.0.d0) THEN
               h(1,i)=i-j
               k=k+1
            ENDIF
            z=0.d0
            IF(i+j.le.n) z=y(i+j)
            IF(z.gt.0.d0) THEN
               h(2,i)=i+j
               k=k+1
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in univariate local polynomial aws
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpawsunw(n,dp1,dp2,dp3,y,xd,fix,mfamily,theta0,theta,
     1                    bii,bi,bi2,bi0,bi02,ai,ni,lam,h,h0,kern,cb,
     2                    dmat,thij,psix,wghts,wghts0,spmax,iter,smw,
     3                    work)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     y          observed values at design points
C     xd         gridlength
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
      integer n,dp1,dp2,i,j,k,l,je,ja,kern,mfamily,h0(2,n),ih,iter,
     1        dp3,info
      logical aws,fix(n),smw
      real*8 y(n),xd,psix(dp2),theta(dp1,n),bi2(dp2,n),bi(dp2,n),
     1     ai(dp1),lam,spmax,dmat(dp1,dp1),lkern,thij(dp1),cb(dp1,dp1),
     2     bi0(dp2,n),lambda,h,ni(n),ha,ha2,bi02(dp2,n),wghts(n,n),
     3     wghts0(n),theta0(dp1,n),bii(dp3,n),d,work(n)
      external lkern
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      aws=lam.lt.1.d20
      DO i=1,n
C        loop over design points, initialization
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
	 lambda=lam
         ih=h
         ja=min0(h0(1,i),max0(1,i-ih))
         je=max0(h0(2,i),min0(n,i+ih))
         ha=max0(i-ja,je-i)
         ha=dmax1(h,ha)*xd
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
         call lpuniwgt(n,i,ja,je,dp1,dp2,xd,ha2,h0(1,i),aws,kern,
     1              theta0,spmax,psix,thij,dmat,cb,wghts(1,i),wghts0)
         IF (smw) THEN
            call smwghtun(n,i,ja,je,wghts(1,i),work,h0(1,i),kern)
C            call smwghtun(n,i,ja,je,wghts0,work,h0(1,i),kern)
         END IF
         ni(i)=0.d0
         DO j=ja,je
            ni(i)=ni(i)+wghts(j,i)
         END DO
C
C          Now do the iterations to obtain the new estimates
C
         call lpuniit(mfamily,n,dp1,dp2,i,ja,je,y,xd,theta(1,i),psix,
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
      subroutine lpawsuni(n,dp1,dp2,dp3,y,xd,fix,mfamily,theta0,theta,
     1                    bii,bi,bi2,bi0,bi02,ai,ni,lam,h,h0,kern,cb,
     2                    dmat,thij,psix,wghts,wghts0,spmax,iter,smw,
     3                    work)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (2*p+1)
C     y          observed values at design points
C     xd         gridlength
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
      integer n,dp1,dp2,i,j,k,l,je,ja,kern,mfamily,h0(2,n),ih,iter,
     1        dp3,info
      logical aws,fix(n),smw
      real*8 y(n),xd,psix(dp2),theta(dp1,n),bi2(dp2,n),bi(dp2,n),
     1     ai(dp1),lam,spmax,dmat(dp1,dp1),lkern,thij(dp1),cb(dp1,dp1),
     2     bi0(dp2,n),lambda,h,ni(n),ha,ha2,bi02(dp2,n),wghts(n),
     3     wghts0(n),theta0(dp1,n),bii(dp3,n),d,work(n)
      external lkern
C     
C     in case of dp1==1  lawsuni should be preferred  (p=0)
C
      aws=lam.lt.1.d20
      DO i=1,n
C        loop over design points, initialization
         IF (fix(i)) CYCLE
C    nothing to do, final estimate is already fixed by control 
	 lambda=lam
         ih=h
C         ja=min0(h0(1,i),max0(1,i-ih))
C         je=max0(h0(2,i),min0(n,i+ih))
         ja=max0(1,h0(1,i)-ih)
         je=min0(n,h0(2,i)+ih)
         ha=max0(i-ja,je-i)
         ha=dmax1(h,ha+h-ih)*xd
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
         call lpuniwgt(n,i,ja,je,dp1,dp2,xd,ha2,h0(1,i),aws,kern,
     1                theta0,spmax,psix,thij,dmat,cb,wghts,wghts0)
         IF (smw) THEN
            call smwghtun(n,i,ja,je,wghts,work,h0(1,i),kern)
C            call smwghtun(n,i,ja,je,wghts0,work,h0(1,i),kern)
         END IF
         ni(i)=0.d0
         DO j=ja,je
            ni(i)=ni(i)+wghts(j)
         END DO
C
C          Now do the iterations to obtain the new estimates
C
         call lpuniit(mfamily,n,dp1,dp2,i,ja,je,y,xd,theta(1,i),psix,
     1                dmat,wghts,wghts0,thij,bi(1,i),ai,bi0(1,i),
     2                bi2(1,i),bi02(1,i),iter,info)
         if(info.gt.0) fix(i)=.TRUE.
C   singularity or failed convergence keep old estimate in this case
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
         DO j=1,dp1
            DO k=1,dp1
               IF (j.gt.k) then 
                  dmat(j,k)=0.0d0
               ELSE
                  dmat(j,k)=bi(j+k-1,i)
               END IF
            END DO
	 END DO
         call invers(dmat,dp1,info)
         IF (info.ne.0) CYCLE 
         d=0.d0
         DO j=1,dp1
            DO k=1,dp1
               d=d+dmat(1,k)*dmat(1,j)*bi2(j+k-1,i)
            END DO
         END DO
         d=dsqrt(2*d)*1.96
         conf(1,i)=theta(i)-d
         conf(2,i)=theta(i)+d
      END DO
      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (univariate polynomial aws)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bibi2ibi(n,dp1,dp2,dp3,bi,bi2,erg,dmat,dmat2)
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
      implicit logical(a-z)
      integer n,dp1,dp2,dp3
      real*8 bi(dp2,n),bi2(dp2,n),erg(dp3,n),dmat(dp1,dp1),
     1       dmat2(dp1,dp1)
      integer ii,i,j,k,l,m,info
      real*8 d
      DO ii=1,n
         DO j=1,dp1
            DO k=1,dp1
               IF (j.gt.k) then 
                  dmat(j,k)=0.0d0
               ELSE
                  dmat(j,k)=bi2(j+k-1,ii)
               END IF
            END DO
	 END DO
         call invers(dmat,dp1,info)
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
C    smooth univariate weights
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smwghtun(n,i,ja,je,wghts,work,h,kern)
      integer n,kern,i,ja,je,h(2)
      real*8 wghts(n),work(n)
      integer j,k,ka,ke,ih
      real*8 maxwght,lkern,d
      external lkern
      ih=max0(i-h(1),h(2)-i)
      d=0.d0
      DO k=ja,ja+ih
         if(k.le.je) d=d+wghts(k)
      END DO
      DO j=ja,je
         work(j)=d
         ka=j-ih
         ke=j+ih
         if(ke.le.je) d=d+wghts(ke)
         if(ka.ge.ja) d=d-wghts(ka)
      END DO
      maxwght=work(i)
C      DO j=ja,je
C         if(work(j).gt.maxwght) maxwght=work(j)
C      END DO
      DO j=ja,je
         wghts(j)=dmin1(1.d0,work(j)/maxwght)
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    one iteration for local GLM estimates
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpuniit(mfamily,n,dp1,dp2,i,ja,je,y,xd,theta,psix,
     1           dmat,wghts,wghts0,d,bi,ai,bi0,bi2,bi02,iter,info)
      implicit logical (a-z)
      integer n,dp1,dp2,mfamily,i,ja,je,iter,info
      real*8 y(n),theta(dp1),psix(dp2),bi(dp2),ai(dp1),bi0(dp2),
     1       bi2(dp2),bi02(dp2),d(dp1),wghts(n),wghts0(n),xd,
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
            xij=(j-i)*xd
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lpuniwgt(n,i,ja,je,dp1,dp2,xd,ha2,h0,aws,kern,theta,
     1                    spmax,psix,thij,dmat,cb,wghts,wghts0)
      implicit logical (a-z)
      integer n,i,ja,je,dp2,dp1,kern,h0(2)
      real*8 xd,ha2,wghts(n),wghts0(n),psix(dp2),theta(dp1,n),
     1       thij(dp1),dmat(dp1,dp1),cb(dp1,dp1),spmax
      logical aws
      integer j,k,l
      real*8 xij,z,sij,thijl,wij,wij0,lkern
      external lkern
      DO j=ja,je
         xij=(j-i)*xd
         wij0=lkern(kern,xij*xij/ha2)
         wghts0(j)=wij0
         z=1.d0
         DO k=1,dp2
            psix(k)=z
            z=z*xij
         END DO
C         IF ((j.lt.h0(1).or.j.gt.h0(2)).and.aws) THEN
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Calculate contribution of Y_j to ai and bi in univariate local polynomial aws (GLM)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine guniaibi(mfamily,dp1,dp2,psix,theta,y,bii,ai,
     1                    bi0,bi2,bi02,wij,wij0)
      implicit logical (a-z)
      integer mfamily,dp1,dp2
      real*8 psix(dp1),theta(dp1),bii(dp2),ai(dp1),bi0(dp2),wij,wij0,
     1       y,ymz1,z,bi2(dp2),bi02(dp2),wij2,wij02
      integer i
      real*8 psith,z1,z2
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
         IF (psith.gt.5.d0) THEN 
             z1=psith-5.d0
C             z1=dexp(5.d0)*(1.d0+z1*(1.d0+.5d0*z1))
             z1=148.4132d0*(1.d0+z1*(1.d0+.5d0*z1*(1.d0+z1/6.d0)))
             z2=z1*dexp(-psith)
C             wij=0.d0
             IF(wij.gt.0.d0) ymz1=z2*y-z1
C             call dblepr("psith",5,psith,1)
         ELSE 
             z1=dexp(psith)
             IF(wij.gt.0.d0) ymz1=y-z1
         ENDIF
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
         z2=dexp(psith)
         z1=z2/(1+z2)
         z2=z1/(1+z2)
         DO i=1,dp2
            bii(i)=bii(i)+z1*psix(i)*wij
            bi0(i)=bi0(i)+z1*psix(i)*wij0
            IF (i.le.dp1) ai(i)=ai(i)+(y-z2)*psix(i)*wij
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
      integer j,k
      DO j=1,dp1
         DO k=1,dp1
            IF (j.gt.k) then 
               dmat(j,k)=0.0d0
            ELSE
               dmat(j,k)=bi(j+k-1)
            END IF
         END DO
      END DO
      call invers(dmat,dp1,info)
C      now dmat contains inverse of B_i 
C      now calculate theta as B_i^{-1} A_i
      DO j=1,dp1
         d(j)=0.0d0
         DO k=1,dp1
            d(j)=d(j)+dmat(j,k)*ai(k)  
         END DO
      END DO
C     just keep the old estimate if info > 0
      RETURN
      END      
