      double precision function lncchi2(sigma,ni,ksi,wj,sj,L,n,
     1                                  work, x0, flb, nfb)
C
C  compute local weighted noncentral chi^2 log-likelihood * (-1)
C
C  sigma - parameter to estimate
C  ni    - sum(wj)
C  ksi   - sum(wj*Sj^2)/ni
C  wj    - local weights
C  sj    - observed values
C  L     - df/2
C  n     - number of local weights/observations
C
      implicit none
      integer n,nfb
      double precision sigma,ni,ksi,wj(n),sj(n),L,work(*),flb(nfb)
      integer j,ib
      double precision eta,z,sig2,zs,sl,lm1,za,x0,rb,eps,sl0,nfb1
      double precision bessliex
      external bessliex
      eps=1.d-16
      lm1=L-1
      nfb1 = nfb-1
C define level for use of large value approximation NIST 10.30.4
      sig2=sigma*sigma
      eta=0.d0
      sl=ksi-2.d0*L*sig2
      if(sl.lt.eps) THEN
         z = sqrt(eps)
         sl0 = eps
         sig2 = (ksi-eps)/2.d0/L
      ELSE
         z = sqrt(sl)
         sl0 = sl
      END IF
      zs=z/sig2
      DO j=1,n
         if(wj(j).gt.0.d0) THEN
            za=sj(j)*zs
            if(za.lt.x0) THEN
C calculate for small arguments za for better accuracy than interpolation
               za=log(bessliex(za,lm1,1.d0,work))
            ELSE IF (za.gt.nfb1) THEN
C  asymptotic large value approximation
               za=za-log(za*6.283185d0)/2.d0
            ELSE
C  linear interpolation between tabulated values
               ib = int(za)
               rb = za-dble(ib)
               za = (1.d0-rb)*flb(ib)+rb*flb(ib+1)
            END IF
C  taylor series approximation of log(besselI(za,lm1)) in x0
C  avoids overflow and costly evaluation for large za
         END IF
         eta=eta+wj(j)*za
      END DO
      lncchi2=ksi/sig2+log(sig2)+lm1/2.d0*log(sl0)-eta/ni
      if(sl.le.eps) lncchi2=lncchi2*(1+(eps-sl)/eps)
C penalize for domain violation for theta
      RETURN
      END
C
C    Minimization of noncentral chi2-likelihood
C    Adapted from procedure localmin in Richard Brent, Algorithms for
C    Minimization without Derivatives, Prentice-Hall, Inc. (1973)
C
      subroutine localmin(low,up,wj,sj,L,n,tol,maxit,work,xmin,fmin,
     1                    flb,nfb)
      implicit none
      integer n,maxit,nfb
      double precision low,up,wj(n),sj(n),L,tol,xmin,fmin,work(*)
      double precision goldc,a,b,d,e,eps,xm,p,q,r,eps1,eps2,u,v,w,fu,
     1       fv,fw,fx,x,x0,flb(nfb)
      double precision ni,ksi,sjj,Lm1
      integer it,j
      logical gsect
      double precision lncchi2,lgammaf,bessliex
      external lncchi2,lgammaf,bessliex
      eps=1d-8
      goldc=0.381966
C  Initialize
C  ni    - sum(wj)
C  ksi   - sum(wj*Sj^2)/ni
      Lm1=L-1
      x0 = max(3.d0*L,1.d1)
C  x0 such that interpolation of log(I(x,L-1)) has rel. error < 1e-4
      ni=0.d0
      ksi=0.d0
      DO j=1,n
C         if(wj(j).gt.0.d0) THEN
            ni=ni+wj(j)
            sjj=sj(j)
            ksi=ksi+wj(j)*sjj*sjj
C         END IF
      END DO
      ksi=ksi/ni
      a=low
      b=up
      v=a+goldc*(b-a)
      w=v
      x=v
      e=0.d0
      d=0.d0
      fx=lncchi2(x,ni,ksi,wj,sj,L,n,work,x0,flb,nfb)
      fv=fx
      fw=fx
C  Search for minimum
      DO it=1,maxit
         xm=0.5d0*(a+b)
         eps1=eps*abs(x)+tol/3.d0
         eps2=2.d0*eps1
         IF (abs(x-xm).le.(eps2-0.5d0*(b-a)) ) EXIT
         gsect = .TRUE.
         IF (abs(e) .gt. eps1) THEN
            gsect = .FALSE.
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.d0*(q-r)
            IF (q .gt. 0.d0) p=-p
            q=abs(q)
            r=e
            e=d
            IF (abs(p).ge.abs(0.5d0*q*r).or.
     1                  p.le.q*(a-x).or.p.ge.q*(b-x)) THEN
                gsect = .TRUE.
            ELSE
               d=p/q
               u=x+d
               IF ((u-a).lt.eps2.or.(b-u).lt.eps2) d=sign(eps1,xm-x)
            END IF
         END IF
         IF (gsect) THEN
            IF (x .ge. xm) THEN
               e=a-x
            ELSE
               e=b-x
            END IF
            d=goldc*e
         END IF
         IF (abs(d) .ge. eps1) THEN
            u=x+d
         ELSE
            u=x+sign(eps1,d)
         END IF
         fu=lncchi2(u,ni,ksi,wj,sj,L,n,work,x0,flb,nfb)
         IF (fu .le. fx) THEN
            IF (u .ge. x) THEN
               a=x
            ELSE
               b=x
            END IF
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
         ELSE
            IF (u.lt.x) THEN
               a=u
            ELSE
               b=u
            END IF
            IF (fu.le.fw.or.w.eq.x) THEN
               v=w
               fv=fw
               w=u
               fw=fu
            ELSE IF (fu.le.fv.or.v.eq.x.or.v.eq.w) THEN
               v=u
               fv=fu
            END IF
         END IF
      END DO
      xmin=sqrt(ni/(ni-1.d0))*x
C this seems to correct the bias
      fmin=fx
      RETURN
      END
