CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local polynomial  aws (gridded) !D
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine hiraws1(y,fix,nfix,n,degrp1,adegr,hw,hakt,hhom,
     1        jrange,lambda,theta,bi,vi,ai,kern,spmin,lw,w,
     2        slw,sw)
C   
C    y,       observed values of regression function
C    fix,     logical TRUE fro points where we have nothing to do
C    nfix,    number of successive points with zero weights requiered
C             for fixing the estimate
C    n,       number of observations
C    degrp1,  degree of polynomial + 1
C    adegr,   degree of polynomial used in adaptation 
C    hw,      bandwidth used to smooth weights
C    hakt,    actual bandwidth in aws
C    hhom,    max radii of a circle containing only points with
C                s_ij==0 in preceeding step
C    jrange,  jrange(1,i),jrange(2,i)  defines maximal interval 
C             for positive weights. entries are forced to fulfil
C             jrange(1,i)<=jrange(1,i+1) and  jrange(2,i)<=jrange(2,i+1)
C    lambda,  lambda or lambda*sigma2 for Gaussian models
C    theta,   estimates from last step (input), theta(k,) for
C             k > adegr+1  are fixed as final estimates from last hirarchy
C    bi,      Matrix Bi dim(n,dp2)
C    vi,      individual variances of components of parameter estimate 
C    ai,      \sum  \Psi Wi Y     (output) dim(n,dp1)
C    kern,    specifies the location kernel
C    spmin,   width of plateau in statistical penalty
C    lw,      vector for location weights length 2*ih+1
C    w,       vector of weights length 2*ih+1
C    slw,     vector for smoothed location weights length 2*(ih+ihw)+1
C    sw       vector for smoothed weights length 2*(ih+ihw)+1
C
C   temporary arrays set for maximum degree 3
C
      implicit logical (a-z)
      external kldisth1,lkern
      real*8 kldisth1,lkern
      integer n,kern,degr,degrp1,adegr,nfix,jrange(2,n)
      logical fix(1),lfix
      real*8 y(1),theta(degrp1,n),bi(1),ai(1),lambda(4),spmin,
     1       vi(degrp1,n),hakt,lw(1),w(1),hw,sw(1),slw(1),hhom(n,2)
      integer ih,j1,k,iind,jind,dlw,clw,jw1,adegrp1,
     2        dp1,dp2,ihs,csw,dsw,l,idegr
      real*8 sij,swj(7),swjy(7),z1,wj,
     1       hakt2,zz(7),lwj,yj,hs2,hs,z,cc,spf,
     2       hhommax,hhommin,az1,hfixmax,hnfix,ssij,spmax,
     3       hhomimin,hhomimax
C   arrays with variable length are organized as 
C   theta(n,dp1)
C   bi(dp2,n)
C   arrays of fixed length correspond to degr=2
C   first set dimensions for arrays depending on degree
C      aws=lambda(degr).lt.1.d20
      lfix=nfix.gt.0
      hnfix=nfix
C      hnfix=max(hnfix,.2d0*hakt)
      degr=degrp1-1
      idegr=degrp1-adegr
C      call intpr("idegr",5,idegr,1)
      adegrp1=adegr+1
      spmax=degrp1
      spf=1.d0/(spmax-spmin)
      if(degr.eq.0) THEN
         dp1=1
	 dp2=1
      ELSE IF (degr.eq.1) THEN
         dp1=2
	 dp2=3
      ELSE IF (degr.eq.2) THEN
         dp1=3
	 dp2=5
      ELSE 
         dp1=4
	 dp2=7
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
         hhomimin=hhom(iind,1)
         hhomimax=hhom(iind,2)
         hhommax=hakt
         hhommin=hakt
         hfixmax=max(hhomimin,hhomimax)
C   scaling of sij outside the loop
         DO jw1=1,dlw
            w(jw1)=0.d0
	    jind=jw1-clw+iind
	    if(jind.lt.jrange(1,iind).or.jind.gt.jrange(2,iind)) CYCLE
            wj=lw(jw1)
            z1=jw1-clw
            az1=abs(z1)
            IF (z1.le.-hhomimin.or.z1.ge.hhomimax) THEN
               sij=1.d0
               DO k=0,idegr-1
                  sij=sij*kldisth1(
     1                idegr-k,theta(adegrp1+k,iind),
     2                        theta(adegrp1+k,jind),
     3                z1,vi(adegrp1+k,iind))/lambda(adegrp1+k)
               END DO
               IF (sij.le.spmax) THEN
                  hfixmax=max(hfixmax,az1)
                  IF (sij.gt.spmin) THEN
                      ssij=1.d0-spf*(sij-spmin)
		      w(jw1)=wj*ssij
                      if(z1.gt.0) THEN
                         hhommax=min(hhommax,az1)
                      ELSE
                         hhommin=min(hhommin,az1)
                      END IF
                  ELSE 
                     w(jw1)=wj
                  END IF
	       ELSE 
		  w(jw1)=0.d0
                  if(z1.gt.0) THEN
                     hhommax=min(hhommax,az1)
                  ELSE
                     hhommin=min(hhommin,az1)
                  END IF
	       END IF
	    ELSE
               w(jw1)=wj
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
	 z=max(.1d0,min(z,hw))
	 cc=min(z-1.d0,1.d0/hakt)
         call smwghts1(w,hakt,z,sw,dlw,dsw,cc)
         DO k=1,dp2
            swj(k)=0.d0
         END DO
         DO k=1,dp1
               swjy(k)=0.d0
         END DO
         DO jw1=csw,dsw
	    j1=jw1-csw+iind
	    if(j1.lt.1.or.j1.gt.n) CYCLE
	    z1=jw1-csw
	    lwj=slw(jw1)
	    wj=sw(jw1)
	    if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE  
            DO k=2,dp2
               zz(k)=z1*zz(k-1)
            END DO
	    if(wj.le.0.d0) CYCLE  
	    DO k=1,dp2
               swj(k)=swj(k)+wj*zz(k)
	    END DO
	    yj=y(j1)
	    DO l=1,dp1
               swjy(l)=swjy(l)+wj*zz(l)*yj
	    END DO
         END DO
         DO jw1=1,csw-1
	    j1=jw1-csw+iind
	    if(j1.lt.1.or.j1.gt.n) CYCLE
	    z1=jw1-csw
	    lwj=slw(jw1)
	    wj=sw(jw1)
	    if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE  
            DO k=2,dp2
               zz(k)=z1*zz(k-1)
            END DO
	    if(wj.le.0.d0) CYCLE  
	    DO k=1,dp2
               swj(k)=swj(k)+wj*zz(k)
	    END DO
	    yj=y(j1)
	    DO l=1,dp1
               swjy(l)=swjy(l)+wj*zz(l)*yj
	    END DO
         END DO
         DO k=1,dp1
            ai((iind-1)*dp1+k)=swjy(k)
         END DO
         DO k=1,dp2
            bi((iind-1)*dp2+k)=swj(k)
         END DO
         hhom(iind,1)=hhommin
         hhom(iind,2)=hhommax
         if(lfix.and.hakt-hfixmax.ge.hnfix) THEN
            jind=iind-hfixmax
            jrange(1,iind)=max0(jrange(1,iind),jind)
            jind=iind+hfixmax
            jrange(2,iind)=min0(jrange(2,iind),jind)
            fix(iind)=.TRUE.
         END IF
         call rchkusr()
      END DO
      RETURN
      END
      real*8 function  kldisth1(degr,thi,thj,xji,vi)
      implicit logical (a-z)
      integer degr
      real*8 thi(degr),thj,xji,vi(degr)
      integer k
      real*8 z,vij,zxji
      z=thi(1)-thj
      vij=vi(1)
      if(degr.gt.1) THEN
         zxji=xji
         do k=2,degr
            z=z+thi(k)*zxji
C            vij=vij+vi(k)*zxji*zxji
            zxji=zxji*xji/k
         END DO
      END IF
      kldisth1=z*z/vij
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      Generate estimates from ai and bi (univariate case)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pawsest1(n,dp1,dp2,idegr,ai,bi,theta,vi,dmat,ind)
C    
C     n          number of design points
C     dp1        number of parameters  (p+1)
C     dp2        number of components in bi  (1,6,15)
C     ai         \sum \Psi^T Wi^k Y       
C     bi         \sum \Psi^T Wi^k \Psi    
C     theta      new parameter estimate
C     dmat       working arrays
C     restricted to dp2<=20
      implicit logical (a-z)
      integer n,dp1,dp2,idegr
      real*8 ai(dp1,n),bi(dp2,n),theta(dp1,n),dmat(dp1,dp1),vi(n)
      integer i,j,k,info,ind(dp1,dp1),iii,ipiv(4),lwork,idp1
      real*8 work(4),h,cii(7),aa(7)
      lwork=4
      idp1=idegr+1
C      call intpr("n",1,n,1)
C      call intpr("dp1",3,dp1,1)
C      call intpr("dp2",3,dp2,1)
C      call intpr("idegr",5,idegr,1)
C      call intpr("ind",3,ind,dp1*dp1)
      DO i=1,n
C         call intpr("i",1,i,1)
         cii(1)=1.d0
         h=bi(1,i)
         if(dp1.gt.1) THEN
            DO k=2,dp2
               cii(k)=cii(k-1)*h
            END DO
         END IF
         DO k=1,dp1
	    DO j=k,dp1
               iii=ind(k,j)
	       dmat(k,j)=bi(iii,i)/cii(iii)
               dmat(j,k)=dmat(k,j)
	    END DO
	    aa(k)=ai(k,i)/cii(k)
	 END DO
C     now calculate theta as B_i^{-1} A_i
	 call dsytrf("U",dp1,dmat,dp1,ipiv,work,lwork,info)
         IF (info.ne.0) THEN
            call dblepr("dmat1",5,dmat,dp1*dp1)
            STOP
            CYCLE
         END IF  
         call dsytrs("U",dp1,1,dmat,dp1,ipiv,aa,dp1,info)
         IF (info.ne.0) THEN
            call dblepr("dmat1",5,dmat,dp1*dp1)
            STOP
            CYCLE
         END IF  
         DO j=1,dp1
            theta(j,i)=aa(j)/cii(j)
	 END DO
         call dsytri("U",dp1,dmat,dp1,ipiv,work,info)
         IF (info.ne.0) THEN
            call dblepr("dmat2",5,dmat,dp1*dp1)
            STOP
            CYCLE
         END IF  
         vi(i)=dmat(idp1,idp1)/cii(idp1)/cii(idp1)
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Evaluate jrange to produce consistent intervals
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine jrange(n,jr)
      implicit logical (a-z)
      integer n,jr(2,n)
      integer i
      DO i=2,n
         jr(1,i)=max(jr(1,i),jr(1,i-1))
      END DO
      DO i=n-1,1,-1
         jr(2,i)=min(jr(2,i),jr(2,i+1))
      END DO
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Evaluate jrange to produce consistent intervals
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine jrofsets(n,jr,setS,lsetS)
      implicit logical (a-z)
      integer n,jr(2,n),lsetS,setS(2,lsetS)
      integer i,j,k
      if(lsetS.le.0) RETURN
C   code information from setS in jr
      DO k=1,lsetS 
         i=setS(1,k)
         j=setS(2,k)
         if(i.le.j) THEN
            jr(2,i)=min(jr(2,i),j-1)
         ELSE
            jr(1,i)=max(jr(1,i),j+1)
         END IF
      END DO
C  now complete information in jr 
      call jrange(n,jr)
C  now reinitialize setS
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local polynomial  aws (gridded) !D
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine hiraws1b(y,fix,nfix,n,degrp1,adegr,hw,hakt,hhom,
     1        jrange,setS,lsetS,lambda,theta,bi,vi,ai,kern,spmin,
     2        lw,w,slw,sw)
C   
C    y,       observed values of regression function
C    fix,     logical TRUE fro points where we have nothing to do
C    nfix,    number of successive points with zero weights requiered
C             for fixing the estimate
C    n,       number of observations
C    degrp1,  degree of polynomial + 1
C    adegr,   degree of polynomial used in adaptation 
C    hw,      bandwidth used to smooth weights
C    hakt,    actual bandwidth in aws
C    hhom,    max radii of a circle containing only points with
C                s_ij==0 in preceeding step
C    jrange,  jrange(1,i),jrange(2,i)  defines maximal interval 
C             for positive weights. entries are forced to fulfil
C             jrange(1,i)<=jrange(1,i+1) and  jrange(2,i)<=jrange(2,i+1)
C    setS,    set of identified discontinuities (i,j) such that s_{ij}>1
C    lsetS,   number of points in setS
C    lambda,  lambda or lambda*sigma2 for Gaussian models
C    theta,   estimates from last step (input), theta(k,) for
C             k > adegr+1  are fixed as final estimates from last hirarchy
C    bi,      Matrix Bi dim(n,dp2)
C    vi,      individual variances of components of parameter estimate 
C    ai,      \sum  \Psi Wi Y     (output) dim(n,dp1)
C    kern,    specifies the location kernel
C    spmin,   width of plateau in statistical penalty
C    lw,      vector for location weights length 2*ih+1
C    w,       vector of weights length 2*ih+1
C    slw,     vector for smoothed location weights length 2*(ih+ihw)+1
C    sw       vector for smoothed weights length 2*(ih+ihw)+1
C   
C   temporary arrays set for maximum degree 3
C
      implicit logical (a-z)
      external kldisth1,lkern
      real*8 kldisth1,lkern
      integer n,kern,degr,degrp1,adegr,jrange(2,n),nfix,lsetS,
     1        setS(2,n)
      logical fix(1),lfix
      real*8 y(1),theta(degrp1,n),bi(1),ai(1),lambda(4),spmin,
     1       vi(degrp1,n),hakt,lw(1),w(1),hw,sw(1),slw(1),hhom(n,2)
      integer ih,j1,k,iind,jind,dlw,clw,jw1,adegrp1,
     2        dp1,dp2,ihs,csw,dsw,l,idegr
      real*8 sij,swj(7),swjy(7),z1,wj,
     1       hakt2,zz(7),lwj,yj,hs2,hs,z,cc,spf,
     2       hhommax,hhommin,az1,ssij,spmax,
     3       hhomimin,hhomimax,hfixmax,hnfix
      lfix=nfix.gt.0
      hnfix=nfix
      degr=degrp1-1
      idegr=degrp1-adegr
      adegrp1=adegr+1
      spmax=degrp1
      spf=1.d0/(spmax-spmin)
      if(degr.eq.0) THEN
         dp1=1
	 dp2=1
      ELSE IF (degr.eq.1) THEN
         dp1=2
	 dp2=3
      ELSE IF (degr.eq.2) THEN
         dp1=3
	 dp2=5
      ELSE 
         dp1=4
	 dp2=7
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
         hhomimin=hhom(iind,1)
         hhomimax=hhom(iind,2)
         hhommax=hakt
         hhommin=hakt
         hfixmax=max(hhomimin,hhomimax)
C   scaling of sij outside the loop
         w(clw)=1.d0
         DO jw1=clw-1,1,-1
            w(jw1)=0.d0
	    jind=jw1-clw+iind
	    if(jind.lt.jrange(1,iind)) EXIT
            wj=lw(jw1)
            z1=jw1-clw
            az1=abs(z1)
            IF (z1.le.-hhomimin) THEN
            sij=1.d0
            DO k=0,idegr-1
               sij=sij*kldisth1(
     1             idegr-k,theta(adegrp1+k,iind),
     2                        theta(adegrp1+k,jind),
     3                        z1,vi(adegrp1+k,iind))/lambda(adegrp1+k)
            END DO
               IF (sij.le.spmax) THEN
                  hfixmax=max(hfixmax,az1)
                  IF (sij.gt.spmin) THEN
                      ssij=1.d0-spf*(sij-spmin)
		      w(jw1)=wj*ssij
                      if(z1.gt.0) THEN
                         hhommax=min(hhommax,az1)
                      ELSE
                         hhommin=min(hhommin,az1)
                      END IF
                  ELSE 
                     w(jw1)=wj
                  END IF
	       ELSE 
		  w(jw1)=0.d0
                  setS(1,lsetS)=iind
                  setS(2,lsetS)=jind
                  lsetS=lsetS+1
                  EXIT
                  if(z1.gt.0) THEN
                     hhommax=min(hhommax,az1)
                  ELSE
                     hhommin=min(hhommin,az1)
                  END IF
	       END IF
	    ELSE
               w(jw1)=wj
            END IF
         END DO
         DO jw1=clw+1,dlw
            w(jw1)=0.d0
	    jind=jw1-clw+iind
	    if(jind.gt.jrange(2,iind)) EXIT
            wj=lw(jw1)
            z1=jw1-clw
            az1=abs(z1)
            IF (z1.ge.hhomimax) THEN
            sij=1.d0
            DO k=0,idegr-1
               sij=sij*kldisth1(
     1             idegr-k,theta(adegrp1+k,iind),
     2                        theta(adegrp1+k,jind),
     3                        z1,vi(adegrp1+k,iind))/lambda(adegrp1+k)
            END DO
               IF (sij.le.spmax) THEN
                  IF (sij.gt.spmin) THEN
                      ssij=1.d0-spf*(sij-spmin)
		      w(jw1)=wj*ssij
                      if(z1.gt.0) THEN
                         hhommax=min(hhommax,az1)
                      ELSE
                         hhommin=min(hhommin,az1)
                      END IF
                  ELSE 
                     w(jw1)=wj
                  END IF
	       ELSE 
		  w(jw1)=0.d0
                  setS(1,lsetS)=iind
                  setS(2,lsetS)=jind
                  lsetS=lsetS+1
                  EXIT
                  if(z1.gt.0) THEN
                     hhommax=min(hhommax,az1)
                  ELSE
                     hhommin=min(hhommin,az1)
                  END IF
	       END IF
	    ELSE
               w(jw1)=wj
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
	 z=max(.1d0,min(z,hw))
	 cc=min(z-1.d0,1.d0/hakt)
         call smwghts1(w,hakt,z,sw,dlw,dsw,cc)
         DO k=1,dp2
            swj(k)=0.d0
         END DO
         DO k=1,dp1
               swjy(k)=0.d0
         END DO
         DO jw1=csw,dsw
	    j1=jw1-csw+iind
	    if(j1.lt.1.or.j1.gt.n) CYCLE
	    z1=jw1-csw
	    lwj=slw(jw1)
	    wj=sw(jw1)
	    if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE  
            DO k=2,dp2
               zz(k)=z1*zz(k-1)
            END DO
	    if(wj.le.0.d0) CYCLE  
	    DO k=1,dp2
               swj(k)=swj(k)+wj*zz(k)
	    END DO
	    yj=y(j1)
	    DO l=1,dp1
               swjy(l)=swjy(l)+wj*zz(l)*yj
	    END DO
         END DO
         DO jw1=1,csw-1
	    j1=jw1-csw+iind
	    if(j1.lt.1.or.j1.gt.n) CYCLE
	    z1=jw1-csw
	    lwj=slw(jw1)
	    wj=sw(jw1)
	    if(lwj.le.0.d0.and.wj.le.0.d0) CYCLE  
            DO k=2,dp2
               zz(k)=z1*zz(k-1)
            END DO
	    if(wj.le.0.d0) CYCLE  
	    DO k=1,dp2
               swj(k)=swj(k)+wj*zz(k)
	    END DO
	    yj=y(j1)
	    DO l=1,dp1
               swjy(l)=swjy(l)+wj*zz(l)*yj
	    END DO
         END DO
         DO k=1,dp1
            ai((iind-1)*dp1+k)=swjy(k)
         END DO
         DO k=1,dp2
            bi((iind-1)*dp2+k)=swj(k)
         END DO
         hhom(iind,1)=hhommin
         hhom(iind,2)=hhommax
         if(lfix.and.hakt-hfixmax.ge.hnfix) THEN
            jind=iind-hfixmax
            jrange(1,iind)=max0(jrange(1,iind),jind)
            jind=iind+hfixmax
            jrange(2,iind)=min0(jrange(2,iind),jind)
            fix(iind)=.TRUE.
         END IF
         call rchkusr()
      END DO
      RETURN
      END
