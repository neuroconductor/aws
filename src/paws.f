CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant  aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awspuni(y,fix,n,degr,hw,hakt,lambda,theta,bi,
     1        bi2,bi0,ai,model,kern,spmax,lw,w,slw,sw,ind)
C   
C   y        observed values of regression function
C   fix      logical TRUE fro points where we have nothing to do
C   n1,n2    design dimensions
C   degr     degree of polynomials 0,1 or 2
C   hw       bandwidth used to smooth weights
C   hakt     actual bandwidth in aws
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       Matrix Bi dim(n1,n2,dp2)
C   bi2      Matrix Bi dim(n1,n2,dp2) (with wij^2 instead of wij)
C   bi0      Matrix Bi0 dim(n1,n2,dp2) (with location weights only)
C   ai       \sum  Wi Y     (output) dim(n1,n2,dp1)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   lw       array of location weights dim(dlw,dlw) dlw=2*ih+1
C   w        array of weights dim(dlw,dlw)
C   sw       array of "smoothed" weights dim(dls,dls) dls=2*(ih+ihw)+1
C   
C   temporary arrays set for maximum degree 2
C
      implicit logical (a-z)
      external kldistp,lkern
      real*8 kldistp,lkern
      integer n,model,kern,degr,ind(1)
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,spmax,
     1       bi2(1),hakt,lw(1),w(1),hw,sw(1),slw(1)
      integer ih,i1,j1,k,iind,jind,jwind,dlw,clw,jw1,
     2        dp1,dp2,ihs,csw,dsw,l,dsw2
      real*8 bii(5),sij,swj(5),swj2(5),swj0(5),swjy(3),z1,wj,
     1       hakt2,thij(3),zz(5),lwj,yj,hs2,hs,z,cc
C   arrays with variable length are organized as 
C   theta(n,dp1)
C   bi(n,dp2)
C   arrays of fixed length correspond to degr=2
C   first set dimensions for arrays depending on degree
      aws=lambda.lt.1.d20
      if(degr.eq.0) THEN
         dp1=1
	 dp2=1
      ELSE IF (degr.eq.1) THEN
         dp1=2
	 dp2=3
      ELSE 
         dp1=3
	 dp2=5
      END IF
      hakt2=hakt*hakt
      ih=hakt
      dlw=2*ih+1
      clw=ih+1
      hs=hakt+hw
      hs2=hs*hs
      ihs=hs
      dsw=2*ihs+1
      csw=ihs+1
C   compute location weights first
      DO j1=1,dlw
         z1=clw-j1
         lw(j1)=lkern(kern,z1*z1/hakt2)
      END DO
      cc=0.0d0
      call smwghts1(lw,hakt,hw,slw,dlw,dsw,cc)
C  now stochastic term
      zz(1)=1.d0
      DO iind=1,n
         IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
         DO k=1,dp2
            bii(k)=bi(iind+(k-1)*n)/lambda
	 END DO
C   scaling of sij outside the loop
         DO jw1=1,dlw
            w(jw1)=0.d0
	    jind=jw1-clw+iind
	    if(jind.lt.1.or.jind.gt.n) CYCLE
            wj=lw(jw1)
            z1=jw1-clw
            zz(2)=z1
            zz(3)=z1*z1
            IF (aws) THEN
	       DO k=1,dp1
	          thij(k)=theta(jind+(k-1)*n)
               END DO
               thij(1)=thij(1)-thij(2)*z1
               IF (dp1.gt.2) THEN
                  thij(1)=thij(1)+thij(3)*zz(3)
                  thij(2)=thij(2)-2.d0*thij(3)*z1
               END IF
C  
C           get difference of thetas
C
               DO k=1,dp1
                  thij(k)=theta(iind+(k-1)*n)-thij(k)
               END DO
               sij=kldistp(dp1,thij,bii,ind)
               IF (sij.le.spmax) THEN
                  w(jwind)=wj*exp(-sij)
               ELSE
	          w(jwind)=0.d0
               END IF
	    ELSE
               w(jwind)=wj     
            END IF
         END DO
C
C      Smooth the weights
C   
         z=0.d0
	 DO jw1=1,dlw
	    if(jw1.eq.clw) CYCLE
               z=z+w(jw1)
         END DO
	 z=(2.d0-z/2.d0)*hw-1+z/2.d0
	 z=dmax1(.1d0,dmin1(z,hw))
	 cc=dmin1(z-1.d0,1.d0/hakt2)
         call smwghts1(w,hakt,z,sw,dlw,dsw,cc)
         DO k=1,dp2
            swj(k)=0.d0
            swj2(k)=0.d0
            swj0(k)=0.d0
         END DO
         DO k=1,dp1
               swjy(k)=0.d0
         END DO
	 dsw2=dsw*dsw
         DO jw1=1,dsw
	    j1=jw1-csw+i1
	    if(j1.lt.1.or.j1.gt.n) CYCLE
	    z1=jw1-csw
	    lwj=slw(jw1)
	    wj=sw(jw1)
	    if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE  
	    zz(2)=z1
	    zz(3)=z1*z1
	    IF(dp1.gt.2) THEN
	       zz(4)=z1*zz(3)
	       zz(5)=z1*zz(4)
	    END IF
	    DO k=1,dp2
               swj0(k)=swj0(k)+lwj*zz(k)
	    END DO
	    if(wj.le.0.d0) CYCLE  
	    DO k=1,dp2
               swj(k)=swj(k)+wj*zz(k)
               swj2(k)=swj2(k)+wj*wj*zz(k)
	    END DO
	    yj=y(j1)
	    DO l=1,dp1
               swjy(l)=swjy(l)+wj*zz(l)*yj
	    END DO
         END DO
         DO k=1,dp1
            ai(iind+(k-1)*n)=swjy(k)
         END DO
         DO k=1,dp2
            bi(iind+(k-1)*n)=swj(k)
            bi2(iind+(k-1)*n)=swj2(k)
            bi0(iind+(k-1)*n)=swj0(k)
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant  aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awspbi(y,fix,n1,n2,degr,hw,hakt,lambda,theta,bi,
     1        bi2,bi0,ai,model,kern,spmax,lw,w,slw,sw,ind)
C   
C   y        observed values of regression function
C   fix      logical TRUE fro points where we have nothing to do
C   n1,n2    design dimensions
C   degr     degree of polynomials 0,1 or 2
C   hw       bandwidth used to smooth weights
C   hakt     actual bandwidth in aws
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       Matrix Bi dim(n1,n2,dp2)
C   bi2      Matrix Bi dim(n1,n2,dp2) (with wij^2 instead of wij)
C   bi0      Matrix Bi0 dim(n1,n2,dp2) (with location weights only)
C   ai       \sum  Wi Y     (output) dim(n1,n2,dp1)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   lw       array of location weights dim(dlw,dlw) dlw=2*ih+1
C   w        array of weights dim(dlw,dlw)
C   sw       array of "smoothed" weights dim(dls,dls) dls=2*(ih+ihw)+1
C   
C   temporary arrays set for maximum degree 2
C
      implicit logical (a-z)
      external kldistp,lkern
      real*8 kldistp,lkern
      integer n1,n2,model,kern,degr,ind(1)
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,spmax,
     1       bi2(1),hakt,lw(1),w(1),hw,sw(1),slw(1)
      integer ih,ih1,i1,i2,j1,j2,k,n,iind,jind,jind2,jwind,jwind2,
     1        dlw,clw,jw1,jw2,dp1,dp2,ihs,csw,dsw,l,dsw2
      real*8 bii(15),sij,swj(15),swj2(15),swj0(15),swjy(6),z1,z2,wj,
     1       hakt2,thij(6),zz(15),lwj,yj,hs2,hs,z,cc
C   arrays with variable length are organized as 
C   theta(n1,n2,dp1)
C   bi(n1,n2,dp2)
C   arrays of fixed length correspond to degr=2
C   first set dimensions for arrays depending on degree
      aws=lambda.lt.1.d20
      if(degr.eq.0) THEN
         dp1=1
	 dp2=1
      ELSE IF (degr.eq.1) THEN
         dp1=3
	 dp2=6
      ELSE 
         dp1=6
	 dp2=15
      END IF
      hakt2=hakt*hakt
      ih=hakt
      dlw=2*ih+1
      clw=ih+1
      hs=hakt+hw
      hs2=hs*hs
      ihs=hs
      dsw=2*ihs+1
      csw=ihs+1
      n=n1*n2
C   compute location weights first
      DO j2=1,dlw
         z2=clw-j2
         z2=z2*z2
         ih1=dsqrt(hakt2-z2)
         jind2=(j2-1)*dlw
         DO j1=clw-ih1,clw+ih1
C  first stochastic term
            jind=j1+jind2
            z1=clw-j1
            lw(jind)=lkern(kern,(z1*z1+z2)/hakt2)
         END DO
      END DO
      cc=0.0d0
      call smwghts2(lw,hakt,hw,slw,dlw,dsw,cc)
C  now stochastic term
      zz(1)=1.d0
      DO i2=1,n2
         DO i1=1,n1
            iind=i1+(i2-1)*n1
            IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
            DO k=1,dp2
               bii(k)=bi(iind+(k-1)*n)/lambda
	    END DO
C   scaling of sij outside the loop
            DO jw2=1,dlw
               jwind2=(jw2-1)*dlw
	       DO jw1=1,dlw
		  w(jw1+jwind2)=0.d0
               END DO
	    END DO
            DO jw2=1,dlw
	       j2=jw2-clw+i2
	       if(j2.lt.1.or.j2.gt.n2) CYCLE
	       jind2=(j2-1)*n1
               jwind2=(jw2-1)*dlw
               z2=jw2-clw
C  get directional differences that only depend on i2-j2
               zz(3)=z2
	       zz(6)=z2*z2
               ih1=dsqrt(hakt2-z2*z2)
               DO jw1=clw-ih1,clw+ih1
		  j1=jw1-clw+i1
	          if(j1.lt.1.or.j1.gt.n1) CYCLE
		  jind=j1+jind2
		  jwind=jw1+jwind2
		  wj=lw(jwind)
                  z1=jw1-clw
C  get rest of directional differences 
		  zz(2)=z1
		  zz(4)=z1*z1
		  zz(5)=z1*z2
                  IF (aws) THEN
		     DO k=1,dp1
		        thij(k)=theta(jind+(k-1)*n)
		     END DO
		     thij(1)=thij(1)-thij(2)*z1-thij(3)*z2
                     IF (dp1.gt.3) THEN
                        thij(1)=thij(1)+thij(4)*zz(4)+
     1                          thij(5)*zz(5)+thij(6)*zz(6)
                        thij(2)=thij(2)-thij(5)*z2-2.d0*thij(4)*z1
                        thij(3)=thij(3)-thij(5)*z1-2.d0*thij(6)*z2
		     END IF
C  
C           get difference of thetas
C
		     DO k=1,dp1
                        thij(k)=theta(iind+(k-1)*n)-thij(k)
                     END DO
                   sij=kldistp(dp1,thij,bii,ind)
                     IF (sij.le.spmax) THEN
		        w(jwind)=wj*exp(-sij)
		     ELSE
		        w(jwind)=0.d0
		     END IF
		  ELSE
		     w(jwind)=wj		     
                  END IF
               END DO
            END DO
C
C      Smooth the weights
C   
            z=0.d0
            DO jw2=1,dlw
	       if(jw2.eq.clw) CYCLE
               jwind2=(jw2-1)*dlw
	       DO jw1=1,dlw
	          if(jw1.eq.clw) CYCLE
                     z=z+w(jw1+jwind2)
               END DO
	    END DO
	    z=(2.d0-z/2.d0)*hw-1+z/2.d0
	    z=dmax1(.1d0,dmin1(z,hw))
	    cc=dmin1(z-1.d0,1.d0/hakt2)
            call smwghts2(w,hakt,z,sw,dlw,dsw,cc)
            DO k=1,dp2
               swj(k)=0.d0
               swj2(k)=0.d0
               swj0(k)=0.d0
            END DO
            DO k=1,dp1
               swjy(k)=0.d0
            END DO
	    dsw2=dsw*dsw
            DO jw2=1,dsw
	       j2=jw2-csw+i2
	       if(j2.lt.1.or.j2.gt.n2) CYCLE
	       jind2=(j2-1)*n1
               jwind2=(jw2-1)*dsw
               z2=jw2-csw
               zz(3)=z2
	       zz(6)=z2*z2
	       IF(dp1.gt.3) THEN
		  zz(10)=z2*zz(6)
		  zz(15)=z2*zz(10)  
	       END IF
               ih1=dsqrt(hs2-z2*z2)
               DO jw1=csw-ih1,csw+ih1
		  j1=jw1-csw+i1
	          if(j1.lt.1.or.j1.gt.n1) CYCLE
		  z1=jw1-csw
	          jwind=jw1+jwind2
		  jind=j1+jind2
		  lwj=slw(jwind)
		  wj=sw(jwind)
		  if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE  
		  zz(2)=z1
		  zz(4)=z1*z1
		  zz(5)=z1*z2
		  IF(dp1.gt.3) THEN
		     zz(7)=z1*zz(4)
		     zz(8)=z1*zz(5)
		     zz(9)=z1*zz(6)
		     zz(11)=z1*zz(7)
		     zz(12)=z1*zz(8)
		     zz(13)=z1*zz(9)
		     zz(14)=z1*zz(10)
		  END IF
		  DO k=1,dp2
                     swj0(k)=swj0(k)+lwj*zz(k)
		  END DO
		  if(wj.le.0.d0) CYCLE  
		  DO k=1,dp2
                     swj(k)=swj(k)+wj*zz(k)
                     swj2(k)=swj2(k)+wj*wj*zz(k)
		  END DO
		  yj=y(jind)
		  DO l=1,dp1
                     swjy(l)=swjy(l)+wj*zz(l)*yj
		  END DO
               END DO
            END DO
            DO k=1,dp1
               ai(iind+(k-1)*n)=swjy(k)
            END DO
            DO k=1,dp2
               bi(iind+(k-1)*n)=swj(k)
               bi2(iind+(k-1)*n)=swj2(k)
               bi0(iind+(k-1)*n)=swj0(k)
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant  aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine awsptri(y,fix,n1,n2,n3,degr,hw,hakt,lambda,theta,bi,
     1        bi2,bi0,ai,model,kern,spmax,lw,w,slw,sw,ind)
C   
C   y        observed values of regression function
C   fix      logical TRUE fro points where we have nothing to do
C   n1,n2    design dimensions
C   degr     degree of polynomials 0,1 or 2
C   hw       bandwidth used to smooth weights
C   hakt     actual bandwidth in aws
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       Matrix Bi dim(n1,n2,dp2)
C   bi2      Matrix Bi dim(n1,n2,dp2) (with wij^2 instead of wij)
C   bi0      Matrix Bi0 dim(n1,n2,dp2) (with location weights only)
C   ai       \sum  Wi Y     (output) dim(n1,n2,dp1)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   lw       array of location weights dim(dlw,dlw) dlw=2*ih+1
C   w        array of weights dim(dlw,dlw)
C   sw       array of "smoothed" weights dim(dls,dls) dls=2*(ih+ihw)+1
C   
C   temporary arrays set for maximum degree 2
C
      implicit logical (a-z)
      external kldistp,lkern
      real*8 kldistp,lkern
      integer n1,n2,n3,model,kern,degr,ind(1)
      logical aws,fix(1)
      real*8 y(1),theta(1),bi(1),bi0(1),ai(1),lambda,spmax,
     1       bi2(1),hakt,lw(1),w(1),hw,sw(1),slw(1)
      integer ih,ih1,ih2,i1,i2,i3,j1,j2,j3,k,n,iind,jind,jind2,jwind,
     1        jwind2,jind3,jwind3,jw3,
     1        dlw,clw,jw1,jw2,dp1,dp2,ihs,csw,dsw,l,dsw2,dlw2
      real*8 bii(35),sij,swj(35),swj2(35),swj0(35),swjy(10),z1,z2,z3,
     1       hakt2,thij(10),zz(35),lwj,yj,hs2,hs,z,cc,wj
C   arrays with variable length are organized as 
C   theta(n1,n2,n3,dp1)
C   bi(n1,n2,n3,dp2)
C   arrays of fixed length correspond to degr=2
C   first set dimensions for arrays depending on degree
      aws=lambda.lt.1.d20
      if(degr.eq.0) THEN
         dp1=1
	 dp2=1
      ELSE IF (degr.eq.1) THEN
         dp1=4
	 dp2=10
      ELSE 
         dp1=10
	 dp2=35
      END IF
      hakt2=hakt*hakt
      ih=hakt
      dlw=2*ih+1
      dlw2=dlw*dlw
      clw=ih+1
      hs=hakt+hw
      hs2=hs*hs
      ihs=hs
      dsw=2*ihs+1
      csw=ihs+1
      dsw2=dsw*dsw
      n=n1*n2*n3
C   compute location weights first
      DO j3=1,dlw
         z3=clw-j3
	 z3=z3*z3
         ih2=dsqrt(hakt2-z3)
         jind3=(j3-1)*dlw2
         DO j2=clw-ih1,clw+ih1
            z2=clw-j2
            z2=z2*z2
            ih1=dsqrt(hakt2-z3-z2)
            jind2=(j2-1)*dlw+jind3
            DO j1=clw-ih1,clw+ih1
C  first stochastic term
               jind=j1+jind2
               z1=clw-j1
               lw(jind)=lkern(kern,(z1*z1+z2+z3)/hakt2)
            END DO
         END DO
      END DO
      cc=0.0d0
      call smwghts3(lw,hakt,hw,slw,dlw,dsw,cc)
C  now stochastic term
      zz(1)=1.d0
      DO i3=1,n3
         DO i2=1,n2
            DO i1=1,n1
               iind=i1+(i2-1)*n1+(i3-1)*n1*n2
               IF (fix(iind)) CYCLE
C    nothing to do, final estimate is already fixed by control 
               DO k=1,dp2
                  bii(k)=bi(iind+(k-1)*n)/lambda
	       END DO
C   scaling of sij outside the loop
               DO jw3=1,dlw
	          jwind3=(jw3-1)*dlw2
                  DO jw2=1,dlw
                     jwind2=(jw2-1)*dlw+jwind3
	             DO jw1=1,dlw
		        w(jw1+jwind2)=0.d0
                     END DO
                  END DO
	       END DO
               DO jw3=1,dlw
	          j3=jw3-clw+i3
	          if(j3.lt.1.or.j3.gt.n3) CYCLE
	          jind3=(j3-1)*n1*n2
                  jwind3=(jw3-1)*dlw2
                  z3=jw3-clw
                  zz(4)=z3
	          zz(10)=z3*z3
                  ih2=dsqrt(hakt2-z3*z3)
                  DO jw2=clw-ih2,clw+ih2
	             j2=jw2-clw+i2
	             if(j2.lt.1.or.j2.gt.n2) CYCLE
	             jind2=(j2-1)*n1+jind3
                     jwind2=(jw2-1)*dlw+jwind3
                     z2=jw2-clw
C  get directional differences that only depend on i2-j2
                     zz(3)=z2
	             zz(8)=z2*z2
	             zz(9)=z2*z3
                     ih1=dsqrt(hakt2-z2*z2-z3*z3)
                     DO jw1=clw-ih1,clw+ih1
		        j1=jw1-clw+i1
	                if(j1.lt.1.or.j1.gt.n1) CYCLE
		        jind=j1+jind2
		        jwind=jw1+jwind2
		        wj=lw(jwind)
                        z1=jw1-clw
C  get rest of directional differences 
		        zz(2)=z1
		        zz(5)=z1*z1
		        zz(6)=z1*z2
		        zz(7)=z1*z3
                        IF (aws) THEN
		           DO k=1,dp1
		              thij(k)=theta(jind+(k-1)*n)
		           END DO
		  thij(1)=thij(1)-thij(2)*z1-thij(3)*z2-thij(4)*z3
                           IF (dp1.gt.4) THEN
                  thij(1)=thij(1)+thij(5)*zz(5)+thij(6)*zz(6)+
     1                    thij(7)*zz(7)+thij(8)*zz(8)+thij(9)*zz(9)+
     2                    thij(10)*zz(10)
              thij(2)=thij(2)-thij(7)*z3-thij(6)*z2-2.d0*thij(5)*z1
              thij(3)=thij(3)-thij(6)*z1-thij(9)*z3-2.d0*thij(8)*z2
              thij(4)=thij(4)-thij(7)*z1-thij(9)*z2-2.d0*thij(10)*z3
		           END IF
C  
C           get difference of thetas
C
		           DO k=1,dp1
                              thij(k)=theta(iind+(k-1)*n)-thij(k)
                           END DO
                           sij=kldistp(dp1,thij,bii,ind)
                           IF (sij.le.spmax) THEN
		              w(jwind)=wj*exp(-sij)
		           ELSE
		              w(jwind)=0.d0
		           END IF
		        ELSE
		           w(jwind)=wj		     
                        END IF
                     END DO
		  END DO
               END DO
C
C      Smooth the weights
C   
               z=0.d0
               DO jw3=1,dlw
	          if(jw2.eq.clw) CYCLE
                  jwind3=(jw3-1)*dlw2
                  DO jw2=1,dlw
	             if(jw2.eq.clw) CYCLE
                     jwind2=(jw2-1)*dlw+jwind3
	             DO jw1=1,dlw
	                if(jw1.eq.clw) CYCLE
                        z=z+w(jw1+jwind2)
                     END DO
                  END DO
	       END DO
	       z=(2.d0-z/2.d0)*hw-1+z/2.d0
	       z=dmax1(.1d0,dmin1(z,hw))
	       cc=dmin1(z-1.d0,1.d0/hakt2)
               call smwghts3(w,hakt,z,sw,dlw,dsw,cc)
               DO k=1,dp2
                  swj(k)=0.d0
                  swj2(k)=0.d0
                  swj0(k)=0.d0
               END DO
               DO k=1,dp1
                  swjy(k)=0.d0
               END DO
               DO jw3=1,dsw
	          j3=jw3-csw+i3
	          if(j3.lt.1.or.j3.gt.n3) CYCLE
	          jind3=(j3-1)*n1*n2
                  jwind3=(jw3-1)*dsw2
                  z3=jw3-csw
                  zz(4)=z3
	          zz(10)=z3*z3
	          IF(dp1.gt.4) THEN
		     zz(20)=z3*zz(10)
		     zz(35)=z3*zz(20)  
	          END IF
                  ih2=dsqrt(hs2-z3*z3)
                  DO jw2=csw-ih2,csw+ih2
	             j2=jw2-csw+i2
	             if(j2.lt.1.or.j2.gt.n2) CYCLE
	             jind2=(j2-1)*n1+jind3
                     jwind2=(jw2-1)*dsw+jwind3
                     z2=jw2-csw
                     zz(3)=z2
	             zz(8)=z2*z2
		     zz(9)=z2*z3
	             IF(dp1.gt.4) THEN
		        zz(17)=z2*zz(8)
		        zz(18)=z2*zz(9)
		        zz(19)=z2*zz(10)
		        zz(31)=z2*zz(17)  
		        zz(32)=z2*zz(18)  
		        zz(33)=z2*zz(19)  
		        zz(34)=z2*zz(20)
	             END IF
                     ih1=dsqrt(hs2-z2*z2-z3*z3)
                     DO jw1=csw-ih1,csw+ih1
		        j1=jw1-csw+i1
	                if(j1.lt.1.or.j1.gt.n1) CYCLE
		        z1=jw1-csw
	                jwind=jw1+jwind2
		        jind=j1+jind2
		        lwj=slw(jwind)
		        wj=sw(jwind)
		        if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE  
		        zz(2)=z1
		        zz(5)=z1*z1
		        zz(6)=z1*z2
		        zz(7)=z1*z3
		        IF(dp1.gt.4) THEN
		           zz(11)=z1*zz(5)
		           zz(12)=z1*zz(6)
		           zz(13)=z1*zz(7)
		           zz(14)=z1*zz(8)
		           zz(15)=z1*zz(9)
		           zz(16)=z1*zz(10)
		           zz(21)=z1*zz(11)
		           zz(22)=z1*zz(12)
		           zz(23)=z1*zz(13)
		           zz(24)=z1*zz(14)
		           zz(25)=z1*zz(15)
		           zz(26)=z1*zz(16)
		           zz(27)=z1*zz(17)
		           zz(28)=z1*zz(18)
		           zz(29)=z1*zz(19)
		           zz(30)=z1*zz(20)
		        END IF
		        DO k=1,dp2
                           swj0(k)=swj0(k)+lwj*zz(k)
		        END DO
		        if(wj.le.0.d0) CYCLE  
		        DO k=1,dp2
                           swj(k)=swj(k)+wj*zz(k)
                           swj2(k)=swj2(k)+wj*wj*zz(k)
		        END DO
		        yj=y(jind)
		        DO l=1,dp1
                           swjy(l)=swjy(l)+wj*zz(l)*yj
		        END DO
                     END DO
                  END DO
                  DO k=1,dp1
                     ai(iind+(k-1)*n)=swjy(k)
                  END DO
                  DO k=1,dp2
                     bi(iind+(k-1)*n)=swj(k)
                     bi2(iind+(k-1)*n)=swj2(k)
                     bi0(iind+(k-1)*n)=swj0(k)
                  END DO
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smwghts1(w,hakt,hw,sw,dw,dsw,cc)
C
C  smooth w with epakern and bandwidth hw, result in sw
C
C     w  array of weights dim(dw,dw)   dw=2*ih+1
C     hakt aktual bandwidth in w
C     hw   bandwidth for smoothing of w
C     sw   array of smoothed weights dim(dsw,dsw)   dsw=2*(ihw+ih)+1
C     cc   dumping factor of weights
C
      implicit logical (a-z)
      integer dw,dsw,cw,csw,cdiff
      real*8 w(dw),sw(dsw),hw,hakt,cc
      integer i1,ja1,je1,j1,i10
      real*8 z,z0,z1,hw2,zmax,hakt2,hsw,hsw2,ww
      cw=(dw+1)/2
      csw=(dsw+1)/2
      cdiff=csw-cw
      hsw=hw+hakt
      hsw2=hsw*hsw
      hakt2=hakt*hakt
      hw2=hw*hw
      DO i1=1,dsw
	 sw(i1)=0.d0
      END DO
      IF(cc.le.0.d0) THEN
         DO j1=1,dw
	       sw(j1+cdiff)=w(j1)
	 END DO
      ELSE
         DO i1=1,dsw
	    z1=i1-csw
	    i10=i1-cdiff
	    ja1=max0(i1-2*cdiff,1)
	    je1=min0(i1,dw)
            z=0.d0
	    z0=0.d0
	    DO j1=ja1,je1
	       z1=(i10-j1)
	       z1=z1*z1
	       if(hw2-z1.lt.0.d0) CYCLE
	       ww=(1.d0-z1/hw2)
	       if(ww.lt.1.d0) ww=cc*ww
	        z=z+ww*w(j1)
            END DO
	    sw(i1)=z
	    zmax=dmax1(zmax,z)
         END DO
         DO i1=1,dsw
	       sw(i1)=sw(i1)/zmax
         END DO
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smwghts2(w,hakt,hw,sw,dw,dsw,cc)
C
C  smooth w with epakern and bandwidth hw, result in sw
C
C     w  array of weights dim(dw,dw)   dw=2*ih+1
C     hakt aktual bandwidth in w
C     hw   bandwidth for smoothing of w
C     sw   array of smoothed weights dim(dsw,dsw)   dsw=2*(ihw+ih)+1
C     cc   dumping factor of weights
C
      implicit logical (a-z)
      integer dw,dsw,cw,csw,cdiff
      real*8 w(dw,dw),sw(dsw,dsw),hw,hakt,cc
      integer i1,i2,id,jd,ja1,je1,ja2,je2,j1,j2,i10,i20
      real*8 z,z0,z1,z2,hw2,zmax,hakt2,hsw,hsw2,ww
      cw=(dw+1)/2
      csw=(dsw+1)/2
      cdiff=csw-cw
      hsw=hw+hakt
      hsw2=hsw*hsw
      hakt2=hakt*hakt
      hw2=hw*hw
      DO i1=1,dsw
         DO i2=1,dsw
	    sw(i1,i2)=0.d0
	 END DO
      END DO
      IF(cc.le.0.d0) THEN
         DO j1=1,dw
	    DO j2=1,dw
	       sw(j1+cdiff,j2+cdiff)=w(j1,j2)
	    END DO
	 END DO
      ELSE
         DO i1=1,dsw
	    z1=i1-csw
	    i10=i1-cdiff
	    ja1=max0(i1-2*cdiff,1)
	    je1=min0(i1,dw)
	    id=dsqrt(hsw2-z1*z1)
	    if(csw-id.lt.1) CYCLE
            DO i2=csw-id,csw+id
	       i20=i2-cdiff
               z=0.d0
	       z0=0.d0
	       DO j1=ja1,je1
	          z1=(i10-j1)
	          z1=z1*z1
	          if(hw2-z1.lt.0.d0) CYCLE
	          jd=dsqrt(hw2-z1)
	          ja2=max0(i20-jd,1)
	          je2=min0(i20+jd,dw)
	          DO j2=ja2,je2
	             z2=(i20-j2)
		     ww=(1.d0-(z1+z2*z2)/hw2)
		     if(ww.lt.1.d0) ww=cc*ww
	             z=z+ww*w(j1,j2)
                  END DO
	       END DO
	       sw(i1,i2)=z
	       zmax=dmax1(zmax,z)
            END DO
         END DO
         DO i1=1,dsw
            DO i2=1,dsw
	       sw(i1,i2)=sw(i1,i2)/zmax
            END DO
         END DO
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine smwghts3(w,hakt,hw,sw,dw,dsw,cc)
C
C  smooth w with epakern and bandwidth hw, result in sw
C
C     w  array of weights dim(dw,dw)   dw=2*ih+1
C     hakt aktual bandwidth in w
C     hw   bandwidth for smoothing of w
C     sw   array of smoothed weights dim(dsw,dsw)   dsw=2*(ihw+ih)+1
C     cc   dumping factor of weights
C
      implicit logical (a-z)
      integer dw,dsw,cw,csw,cdiff
      real*8 w(dw,dw,dw),sw(dsw,dsw,dsw),hw,hakt,cc
      integer i1,i2,i3,id,jd,ja1,je1,ja2,je2,ja3,je3,j1,j2,j3,
     1        i10,i20,i30,id2,jd2
      real*8 z,z0,z1,z2,z3,hw2,zmax,hakt2,hsw,hsw2,ww
      cw=(dw+1)/2
      csw=(dsw+1)/2
      cdiff=csw-cw
      hsw=hw+hakt
      hsw2=hsw*hsw
      hakt2=hakt*hakt
      hw2=hw*hw
      DO i1=1,dsw
         DO i2=1,dsw
            DO i3=1,dsw
	       sw(i1,i2,i3)=0.d0
	    END DO
	 END DO
      END DO
      IF(cc.le.0.d0) THEN
         DO j1=1,dw
	    DO j2=1,dw
	       DO j3=1,dw
	          sw(j1+cdiff,j2+cdiff,j3+cdiff)=w(j1,j2,j3)
	       END DO
	    END DO
	 END DO
      ELSE
         DO i1=1,dsw
	    z1=i1-csw
	    i10=i1-cdiff
	    ja1=max0(i1-2*cdiff,1)
	    je1=min0(i1,dw)
	    id=dsqrt(hsw2-z1*z1)
	    if(csw-id.lt.1) CYCLE
            DO i2=csw-id,csw+id
	       i20=i2-cdiff
	       z2=i2-csw
	       id2=dsqrt(hsw2-z1*z1-z2*z2)
               DO i3=csw-id,csw+id
	          i30=i3-cdiff
                  z=0.d0
	          z0=0.d0
	          DO j1=ja1,je1
	             z1=(i10-j1)
	             z1=z1*z1
	             if(hw2-z1.lt.0.d0) CYCLE
	             jd=dsqrt(hw2-z1)
	             ja2=max0(i20-jd,1)
	             je2=min0(i20+jd,dw)
	             DO j2=ja2,je2
	                z2=(i20-j2)
		        z2=z2*z2+z1
	                if(hw2-z2.lt.0.d0) CYCLE
	                jd2=dsqrt(hw2-z2)
	                ja3=max0(i30-jd2,1)
	                je3=min0(i30+jd2,dw)
                        DO j3=ja3,je3
	                   z3=(i30-j3)
                           ww=(1.d0-(z2+z3*z3)/hw2)
		           if(ww.lt.1.d0) ww=cc*ww
	                   z=z+ww*w(j1,j2,j3)
	                END DO
                     END DO
	          END DO
	          sw(i1,i2,i3)=z
	          zmax=dmax1(zmax,z)
	       END DO
            END DO
         END DO
         DO i1=1,dsw
            DO i2=1,dsw
               DO i3=1,dsw
	          sw(i1,i2,i3)=sw(i1,i2,i3)/zmax
               END DO
            END DO
         END DO
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function kldistp(dp1,thij,bii,ind)
C
C  search for maximum in w within bandwidth hw, result in sw
C
C     dp1  polynomial degree +1
C     thij parameter estimate in j dim(dp1*nwght) for basis in i
C     bii  XTX dim(dp2)
C     wght     weight for color channels
C     nwght    number of positive weights (<=dv)
C     ind   index matrix to access the correct elements in bii
C
      implicit logical (a-z)
      integer dp1,ind(dp1,dp1)
      real*8 thij(1),bii(1),thijl
      integer l,k
      real*8 d
      d=0.d0
      DO l=1,dp1
         thijl=thij(l)
         d=d+bii(ind(l,l))*thijl*thijl
	 IF(l.eq.dp1) CYCLE
         DO k=l+1,dp1
            d=d+2.d0*bii(ind(k,l))*thijl*thij(k)
         END DO
      END DO
      kldistp=d
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (bivariate case)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mpaws(n,dp1,dp2,ai,bi,theta,dmat,ind)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi    
C     theta      new parameter estimate
C     dmat       working arrays
C
      implicit logical (a-z)
      integer n,dp1,dp2
      real*8 ai(n,dp1),bi(n,dp2),theta(n,dp1),dmat(dp1,dp1)
      integer i,j,k,info,ind(dp1,dp1)
      real*8 d
      DO i=1,n
         DO k=1,dp1
            DO j=1,dp1
               IF (j.gt.k) then
                  dmat(j,k)=0.d0
               ELSE
                  dmat(j,k)=bi(i,ind(j,k))
               END IF
            END DO
C            dmat(k,k)=dmat(k,k)*1.001
	 END DO
         call invers(dmat,dp1,info)
         IF (info.gt.0) CYCLE  
C     just keep the old estimate
C     now dmat contains inverse of B_i 
C     now calculate theta as B_i^{-1} A_i
         DO j=1,dp1
            d=0.0d0
            DO k=1,dp1
               d=d+dmat(j,k)*ai(i,k)  
            END DO
            theta(i,j)=d
         END DO
      END DO
      RETURN
      END      
