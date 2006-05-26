CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C   find max of x values within an ellipse with excentricity a, orientation
C   given by (th1,th2) and minimal half axis equal to h 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine rangex2D(th1,th2,a,h,xmax)
      real*8 th1,th2,a,h,xmax
      real*8 phi,nu
      IF (abs(th1).gt.1e-6) THEN
         phi=datan(th2/th1)
	 nu=datan(th2/th1/a)
	 xmax=h*dabs(a*dcos(phi)*dcos(nu)-dsin(phi)*dsin(nu))
      ELSE
         xmax=h
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C   find range of y values for given x 
C   within an ellipse with excentricity a, orientation
C   given by (th1,th2) and minimal half axis equal to h 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine rangey2D(th1,th2,a,h,x,ymin,ymax)
      real*8 th1,th2,a,h,x,ymin,ymax
      real*8 th11,th12,th22,a2,c1,p,q,z
      th11 = th1*th1
      th12 = th1*th2
      th22 = th2*th2
      a2 = a*a
      C1 = a2*th11+th22
      p = x*th12*(a2-1)/C1
      q = (x*x*(a2*th22+th11)-a2)/C1
      z = dsqrt(p*p-q) 
      ymin = p - z
      ymax = p + z
      RETURN
      END 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C   calculate Epanechnicov wght for point (x,y)
C   within an ellipse with excentricity a, orientation
C   given by (th1,th2) and minimal half axis equal to h 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      real*8 function anisowgt(th1,th2,a2,x,y,h2)
      real*8 th1,th2,a2,h2,x,y
      real*8 z1,z2
      z1 = x*th1+y*th2
      z2 = y*th1-x*th2
      anisowgt = (h2 - z1*z1/a2 - z2*z2)/h2
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant aws (gridded)
C   for anisotropic 2D situations
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawsa2D(y,fix,n1,n2,hakt,lambda,theta,bi,bi2,bi0,ai,
     1       phi,ani,model,kern,skern,spmin,spmax,lwght,wght)
C   
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   ai       \sum  Wi Y     (output)
C   model    specifies the probablilistic model for the KL-Distance
C   kern     specifies the location kernel
C   spmax    specifies the truncation point of the stochastic kernel
C   wght     scaling factor for second and third dimension (larger values shrink)
C   
      implicit logical (a-z)
      external kldist,anisowgt
      real*8 kldist,anisowgt
      integer n1,n2,model,kern,skern
      logical aws,fix(1)
      real*8 y(n1,n2),theta(n1,n2),bi(n1,n2),bi0(n1,n2),ai(n1,n2),
     1       lambda,spmax,wght,bi2(n1,n2),phi(n1,n2),ani(n1,n2),hakt,
     2       lwght(1),spmin,spf
      integer ih1,i1,i2,j1,j2,ij1,ij2,j2a,j2e
      real*8 thetai,bii,sij,swj,swj2,swj0,swjy,z1,z2,wj,hakt2,bii0,
     1        anii,phii,h1,a2,th1,th2,ymin,ymax
      hakt2=hakt*hakt
      spf=spmax/(spmax-spmin)
      ih1=hakt
      aws=lambda.lt.1d40
C
C   location wghts depend on local orientation !!!
C
      DO i1=1,n1
         DO i2=1,n2
            thetai=theta(i1,i2)
            bii=bi(i1,i2)/lambda
            bii0=bi0(i1,i2)
            swj=0.d0
	    swj2=0.d0
            swj0=0.d0
            swjy=0.d0
	    anii=ani(i1,i2)
	    a2=anii*anii
	    phii=phi(i1,i2)
	    th1=dcos(phii)
	    th2=dsin(phii)
C      get range for j1
            call rangex2D(th1,th2,anii,hakt,h1)
            ih1=h1
	    DO j1=-ih1,ih1
	       ij1=i1+j1
	       if(ij1.lt.1.or.ij1.gt.n1) CYCLE
C      get range for j2
               z1=j1
               call rangey2D(th1,th2,anii,hakt,z1,ymin,ymax)
	       j2a=ymin
	       j2e=ymax
	       DO j2=j2a,j2e
	          ij2=i2+j2
	          if(ij2.lt.1.or.ij2.gt.n2) CYCLE
                  z2=j2
C      get location wghts
                  wj=anisowgt(th1,th2,a2,z1,z2,hakt2)
                  swj0=swj0+wj
                  IF (aws) THEN
                     sij=bii*kldist(model,thetai,theta(ij1,ij2),bii0)
                     IF (sij.gt.spmax) CYCLE
	             IF (skern.eq.2) THEN
			wj=wj*(1.d0-sij)
		     ELSE
		        IF (sij.gt.spmin) wj=wj*dexp(-spf*(sij-spmin))
		     ENDIF
C   if sij <= spmin  this just keeps the location penalty
C    spmin = 0 corresponds to old choice of K_s 
C   new kernel is flat in [0,spmin] and then decays exponentially
                  END IF
                  swj=swj+wj
                  swj2=swj2+wj*wj
                  swjy=swjy+wj*y(ij1,ij2)
               END DO
            END DO
            ai(i1,i2)=swjy
            bi(i1,i2)=swj
            bi2(i1,i2)=swj2
            bi0(i1,i2)=swj0
            call rchkusr()
         END DO
      END DO
      RETURN
      END
      