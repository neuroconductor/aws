CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform nonadaptiv local constant smoothing
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gaws0(y,n1,n2,n3,hakt,theta,bi,bi0,
     1       thnew,kern,lw,wght)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew       \sum  Wi Y     (output)
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldistd,lkern
      real*8 kldistd,lkern
      integer n1,n2,n3,kern,skern,y(1),theta(1),thnew(1)
      logical aws
      real*8 bi(1),bi0,lambda,spmax,spmin,hakt,lw(1)
      integer ih,ih1,i1,i2,j1,j2,k,n,iind1,jind1,
     1        iind,jind,jind2,jwind2,dlw,clw,jw1,jw2
      real*8 bii,sij,swj,swj0,swjy,z1,z2,wj,hakt2,bii0,spf
      hakt2=hakt*hakt
C      spf=spmax/(spmax-spmin)
      spf=spmax/(spmax-spmin)
      ih=hakt
      dlw=2*ih+1
      clw=ih+1
      aws=lambda.lt.1d40
      n=n1*n2*n3
      bii0=bi0
C   compute location weights first
      swj0=0.d0
      DO j3=1,dlw
         DO j2=1,dlw
            z2=clw-j2
            z2=z2*z2
            ih1=dsqrt(hakt2-z2)
            jind2=(j2-1)*dlw
            DO j1=clw-ih1,clw+ih1
C  first stochastic term
            jind=j1+jind2
            z1=clw-j1
            wj=lkern(kern,(z1*z1+z2)/hakt2)
            swj0=swj0+wj
            lw(jind)=wj
            END DO
         END DO
      END DO
      bi0=swj0
      call rchkusr()
      DO i2=1,n2
         DO i1=1,n1
            iind=i1+(i2-1)*n1
            bii=bi(iind)/lambda
C   scaling of sij outside the loop
            swj=0.d0
            swjy=0.d0
            DO jw2=1,dlw
	       j2=jw2-clw+i2
	       if(j2.lt.1.or.j2.gt.n2) CYCLE
	       jind2=(j2-1)*n1
               jwind2=(jw2-1)*dlw
               z2=clw-jw2
               ih1=dsqrt(hakt2-z2*z2)
               DO jw1=clw-ih1,clw+ih1
		  j1=jw1-clw+i1
	          if(j1.lt.1.or.j1.gt.n1) CYCLE
		  jind=j1+jind2
		  wj=lw(jw1+jwind2)
                  swj=swj+wj
                  swjy=swjy+wj*y(jind)
               END DO
            END DO
            thnew(iind)=swjy/swj
            bi(iind)=swj
            call rchkusr()
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gvaws(y,n1,n2,vcoef,nvpar,meanvar,chcorr,
     1                   hakt,lambda,theta,bi,bi0,thnew,kern,skern,
     2                   spmin,spmax,wghts,lw,swjy)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew       \sum  Wi Y     (output)
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldistgc,lkern
      real*8 kldistgc,lkern
      integer n1,n2,kern,skern,nvpar,y(1),theta(1),thnew(1)
      logical aws
      real*8 bi(1),lambda,spmax,spmin,hakt,lw(1),wghts,bi0,
     2       vcoef(nvpar),meanvar
      integer ih,ih1,i1,i2,j1,j2,ja1,je1,l,k,n,info,i2n1,
     1        iind,jind,jind2,jwind2,dlw,clw,jw1,jw2,m0,thi(4)
      real*8 bii,sij,swj,swjy,z1,z2,wj,hakt2,spf,thij(4),
     1       s2i(16),si(4),swj0
C  s2i, s2ii temporay stor sigma^2_i and its inverse (nneded for KL-distance)
C  maximaum dv = 4
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      ih=hakt
      dlw=2*ih+1
      clw=ih+1
      aws=lambda.lt.1d40
      n=n1*n2
C   compute location weights first
      swj0=0.d0
      DO j2=1,dlw
         z2=clw-j2
         z2=z2*z2
         ih1=dsqrt(hakt2-z2)
         ja1=max0(1,clw-ih1)
         je1=min0(dlw,clw+ih1)
         jind2=(j2-1)*dlw
         DO j1=ja1,je1
C  first stochastic term
            jind=j1+jind2
            z1=clw-j1
            wj=lkern(kern,(z1*z1+z2)/hakt2)
            swj0=swj0+wj
            lw(jind)=wj
         END DO
      END DO
      bi0=swj0
      call rchkusr()
      DO i2=1,n2
         i2n1=(i2-1)*n1
         DO i1=1,n1
            iind=i1+i2n1
            bii=bi(iind)/lambda
C   scaling of sij outside the loop
            swj=0.d0
            swjy=0.d0
            thi=theta(iind)
            si = vcoef(1)
            if(nvpar.gt.1) THEN 
               si = si + vcoef(2) * thi
            END IF
            s2i = 1.d0/dmax1(si,0.1*meanvar)
C set small variances to  0.1 * mean variance
C  Now fill estimated Covariancematrix in pixel i
            DO jw2=1,dlw
	       j2=jw2-clw+i2
	       if(j2.lt.1.or.j2.gt.n2) CYCLE
	       jind2=(j2-1)*n1
               jwind2=(jw2-1)*dlw
               z2=clw-jw2
               ih1=dsqrt(hakt2-z2*z2)
               DO jw1=clw-ih1,clw+ih1
		  j1=jw1-clw+i1
	          if(j1.lt.1.or.j1.gt.n1) CYCLE
		  jind=j1+jind2
                  thij=thi-theta(jind)
		  wj=lw(jw1+jwind2)
                  IF (aws) THEN
                     sij=bii*thij*thij*s2ii
                     IF (sij.gt.spmax) CYCLE
		     IF (skern.eq.1) THEN
C  skern == "Triangle"
                        wj=wj*(1.d0-sij)
		     ELSE
C  skern == "Exp"
		        IF (sij.gt.spmin) wj=wj*dexp(-spf*(sij-spmin))
		     END IF
                  END IF
                  swj=swj+wj
                  swjy=swjy+wj*y(jind)
               END DO
            END DO
            thnew(iind)=swjy/swj
            bi(iind)=swj
            call rchkusr()
         END DO
      END DO
      RETURN
      END
