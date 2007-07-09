      subroutine hiraws(y,n,sigma2,degree,fw,ni,sigma2D,hakt,lambda,
     1                  kern,spmin,yhat,bi,yhatnew,binew)
      implicit logical (a-z)
      integer n,degree,kern
      real*8 y(n),sigma2,fw(degree,n),ni(degree,n),sigma2D(degree),
     1       hakt,lambda,spmin,yhat(n),bi(n),yhatnew(n),binew(n)
      integer i,j,k,ih,ia,ie
      real*8 z,xd,yhati,yhatj,wij,wijk,swij,swijy,yj
      integer lkern
      real*8 skern
      external lkern
      ih=hakt
      DO i=1,n
         ia=max0(1,i-ih)
         ie=min0(n,i+ih)
         yhati=yhat(i)
         swij=0.d0
         swijy=0.d0
         DO j=ia,ie
            xd=j-i
            wij=lkern(kern,xd*xd)
C   we have location weights
C   now test if derivatives are homogeneous
            IF(degree.gt.1) THEN
               DO k=1,degree
                  z=fw(k,j)-fw(k,i)
                  z=z*z/sigma2D(k)/lambda*ni(k,i)
                  wijk=1.d0
                  if(z.gt.spmin) wijk=(1.d0-z)/(1.d0-spmin)
                  if(wijk.lt.0.d0) wijk=0.d0
                  wij=wij*wijk
               END DO
            END IF
            if(wij.le.0.d0) CYCLE
C   Now the statistical penalty
            yhatj=yhat(j)
            yj=y(j)
            if(degree.gt.0) THEN
               z=xd
               DO k=1,degree
                  yhatj=yhati+fw(k,i)*z
                  yj=yj-fw(k,i)*z
                  z=z*xd
               END DO
            END IF
            z=yhatj-yhat(j)
            z=z*z/sigma2/lambda
            if(z.gt.1.d0) CYCLE
            if(z.gt.spmin) THEN
               wij=wij*(1.d0-z)/(1.d0-spmin)
            END IF
            swij=swij+wij
            swijy=swijy+wij*yj
         END DO
         yhatnew(i)=swijy/swij
         binew(i)=swij
      END DO
      RETURN
      END