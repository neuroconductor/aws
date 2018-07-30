      subroutine fillpat3(x,n1,n2,n3,phw,psize,pmat)
      implicit NONE
      integer n1,n2,n3,phw,psize
      double precision x(n1,n2,n3),pmat(n1,n2,n3,psize)
      integer i1,i2,i3,j1,j2,j3,k,l1,l2,l3
      DO i1=1,n1
        DO i2=1,n2
          DO i3=1,n3
            k=0
            DO j1=-phw,phw
               l1=i1+j1
               if(l1.lt.1) l1=2-l1
               if(l1.gt.n1) l1=2*n1-l1
               DO j2=-phw,phw
                  l2=i2+j2
                  if(l2.lt.1) l2=2-l2
                  if(l2.gt.n2) l2=2*n2-l2
                  DO j3=-phw,phw
                     l3=i3+j3
                     if(l3.lt.1) l3=2-l3
                     if(l3.gt.n3) l3=2*n3-l3
                     k=k+1
                     pmat(i1,i2,i3,k)=x(l1,l2,l3)
                  END DO
               END DO
            END DO
          END DO
        END DO
      END DO
      RETURN
      END
      subroutine fillpat2(x,n1,n2,phw,psize,pmat)
      implicit NONE
      integer n1,n2,phw,psize
      double precision x(n1,n2),pmat(n1,n2,psize)
      integer i1,i2,j1,j2,k,l1,l2
      DO i1=1,n1
        DO i2=1,n2
           k=0
           DO j1=-phw,phw
              l1=i1+j1
              if(l1.lt.1) l1=2-l1
              if(l1.gt.n1) l1=2*n1-l1
              DO j2=-phw,phw
                 l2=i2+j2
                 if(l2.lt.1) l2=2-l2
                 if(l2.gt.n2) l2=2*n2-l2
                 k=k+1
                 pmat(i1,i2,k)=x(l1,l2)
              END DO
           END DO
        END DO
      END DO
      RETURN
      END
      subroutine fillpat1(x,n1,phw,psize,pmat)
      implicit NONE
      integer n1,phw,psize
      double precision x(n1),pmat(n1,psize)
      integer i1,j1,k,l1
      DO i1=1,n1
        k=0
        DO j1=-phw,phw
           l1=i1+j1
           if(l1.lt.1) l1=2-l1
           if(l1.gt.n1) l1=2*n1-l1
           k=k+1
           pmat(i1,k)=x(l1)
        END DO
      END DO
      RETURN
      END
      subroutine nlmeans(x,n1,n2,n3,patch,pd,swd,tau,xhat)
      implicit none
      integer n1,n2,n3,pd,swd
      double precision x(n1,n2,n3),patch(pd,n1,n2,n3),tau,
     1                 xhat(n1,n2,n3)
      integer i1,i2,i3,j1,j2,j3
      double precision dij,swi,sywi,wi,pd2
      external enorm
      double precision enorm
      pd2=-2*pd
      pd2=pd2*tau*tau
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(x,xhat,patch,pd,n1,n2,n3,swd,tau,pd2)
C$OMP& PRIVATE(i1,i2,i3,wi,swi,sywi,dij,j1,j2,j3)
C$OMP DO SCHEDULE(GUIDED)
      DO i1=1,n1
        DO i2=1,n2
          DO i3=1,n3
            swi=0.d0
            sywi=0.d0
            DO j1=max(1,i1-swd),min(n1,i1+swd)
              DO j2=max(1,i2-swd),min(n2,i2+swd)
                DO j3=max(1,i3-swd),min(n3,i3+swd)
                  dij=enorm(patch(1,i1,i2,i3),patch(1,j1,j2,j3),pd)
                  wi=dexp(dij/pd2)
                  swi=swi+wi
                  sywi=sywi+wi*x(j1,j2,j3)
                END DO
              END DO
            END DO
            xhat(i1,i2,i3)=sywi/swi
          END DO
        END DO
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(xhat)
      RETURN
      END

      double precision function enorm(x,y,n)
      implicit none
      integer n,i
      double precision x(n),y(n),s,si
      s=0.d0
      DO i=1,n
        si=x(i)-y(i)
        s=s+si*si
      END DO
      enorm=s
      return
      end
