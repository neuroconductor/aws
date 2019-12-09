CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in patchwise local constant aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pvawsb(y,mask,nv,n1,n2,n3,hakt,lambda,theta,bi,
     1                bin,thnew,ncores,spmin,lwght,wght,swjy,
     2                np1,np2,np3,thpatch,biipatch)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew    \sum  Wi Y / bi     (output)
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit none

      integer nv,n1,n2,n3,ncores
      logical aws,mask(*)
      double precision y(nv,*),theta(nv,*),bi(*),thnew(nv,*),lambda,
     1  wght(2),hakt,lwght(*),spmin,spf,swjy(nv,ncores),bin(*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,k,thrednr
      double precision bii,biinv,sij,swj,z,z1,z2,z3,wj,hakt2,
     1        w1,w2,spmb,sijp
      integer np1,np2,np3
      double precision thpatch(nv,np1*np2*np3,*),biipatch(np1*np2*np3,*)
      integer ip1,ip2,ip3,nph1,nph2,nph3,pc,ipind,pcj,jp1,jp2,jp3,jpind
      external lkern
      double precision lkern
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      thrednr = 1
C just to prevent a compiler warning
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih1=FLOOR(hakt)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      nph1=(np1-1)/2
      nph2=(np2-1)/2
      nph3=(np3-1)/2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(2,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thnew,bi,nv,n1,n2,n3,hakt2,theta,bin,
C$OMP& ih3,lwght,wght,y,swjy,mask,thpatch,biipatch,nph1,nph2,nph3)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,bii,biinv,swj,spmb,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,z,thrednr,ip1,ip2,ip3,pc,ipind,pcj,
C$OMP& jp1,jp2,jp3,jpind,sijp)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         if(.not.mask(iind)) CYCLE
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
         pc=0
         DO ip1=i1-nph1,i1+nph1
            if(ip1.le.0.or.ip1.gt.n1) CYCLE
            DO ip2=i2-nph2,i2+nph2
               if(ip2.le.0.or.ip2.gt.n2) CYCLE
               DO ip3=i3-nph3,i3+nph3
                  if(ip3.le.0.or.ip3.gt.n3) CYCLE
                  ipind=ip1+(ip2-1)*n1+(ip3-1)*n12
                  if(.not.mask(ipind)) CYCLE
                  pc=pc+1
                  DO k=1,nv
                     thpatch(k,pc,thrednr)=theta(k,ipind)
                  END DO
                  biipatch(pc,thrednr)=bi(ipind)/lambda
               END DO
            END DO
         END DO

C   scaling of sij outside the loop
         swj=0.d0
         DO k=1,nv
            swjy(k,thrednr)=0.d0
         END DO
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  if(.not.mask(jind)) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  IF (aws) THEN
                     pcj=0
                     sij=0.d0
                     DO ip1=i1-nph1,i1+nph1
                        if(ip1.le.0.or.ip1.gt.n1) CYCLE
                        jp1=ip1+jw1
                        if(jp1.le.0.or.jp1.gt.n1) CYCLE
                        DO ip2=i2-nph2,i2+nph2
                           if(ip2.le.0.or.ip2.gt.n2) CYCLE
                           jp2=ip2+jw2
                           if(jp2.le.0.or.jp2.gt.n2) CYCLE
                           DO ip3=i3-nph3,i3+nph3
                              if(sij.gt.1.d0) CYCLE
                              if(ip3.le.0.or.ip3.gt.n3) CYCLE
                              jp3=ip3+jw3
                              if(jp3.le.0.or.jp3.gt.n3) CYCLE
                              jpind=jp1+(jp2-1)*n1+(jp3-1)*n12
                              if(.not.mask(jpind)) CYCLE
                              pcj=pcj+1
                              sijp=0.d0
                              DO k=1,nv
                                 z=thpatch(k,pcj,thrednr)-theta(k,jpind)
                                 sijp=sijp+z*z
                              END DO
                              sij=max(sij,biipatch(pcj,thrednr)*sijp)
                           END DO
                        END DO
                     END DO
                     IF (sij.ge.1.d0) CYCLE
                     IF (sij.gt.spmin) wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
                  swj=swj+wj
                  DO k=1,nv
                     swjy(k,thrednr)=swjy(k,thrednr)+wj*y(k,jind)
                  END DO
               END DO
            END DO
         END DO
         DO k=1,nv
            thnew(k,iind)=swjy(k,thrednr)/swj
         END DO
         bin(iind)=swj
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thnew,bin)
      RETURN
      END
      subroutine pvaws2b(y,mask,nv,nvd,n1,n2,n3,hakt,lambda,theta,bi,
     1                bin,thnew,invcov,ncores,spmin,lwght,wght,swjy,
     2                np1,np2,np3,thpatch,invcovp,biipatch)
C
C   y        observed values of regression function
C   n1,n2,n3    design dimensions
C   hakt     actual bandwidth
C   lambda   lambda or lambda*sigma2 for Gaussian models
C   theta    estimates from last step   (input)
C   bi       \sum  Wi   (output)
C   thnew    \sum  Wi Y / bi     (output)
C   wght     scaling factor for second and third dimension (larger values shrink)
C
      implicit none

      integer nv,n1,n2,n3,ncores,nvd
      logical aws,mask(*)
      double precision y(nv,*),theta(nv,*),bi(*),thnew(nv,*),lambda,
     1  wght(2),hakt,lwght(*),spmin,spf,swjy(nv,ncores),invcov(nvd,*),
     2  bin(*)
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,k,thrednr
      double precision bii,biinv,sij,swj,z,z1,z2,z3,wj,hakt2,
     1        w1,w2,spmb,sijp
      integer np1,np2,np3,l,m
      double precision thpatch(nv,np1*np2*np3,*),
     1        invcovp(nvd,np1*np2*np3,*),biipatch(np1*np2*np3,*)
      integer ip1,ip2,ip3,nph1,nph2,nph3,pc,ipind,pcj,jp1,jp2,jp3,jpind
      external lkern, KLdistsi
      double precision lkern, KLdistsi
!$      integer omp_get_thread_num
!$      external omp_get_thread_num
      thrednr = 1
C just to prevent a compiler warning
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
      ih1=FLOOR(hakt)
      aws=lambda.lt.1d35
C
C   first calculate location weights
C
      w1=wght(1)
      w2=wght(2)
      ih3=FLOOR(hakt/w2)
      ih2=FLOOR(hakt/w1)
      ih1=FLOOR(hakt)
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      nph1=(np1-1)/2
      nph2=(np2-1)/2
      nph3=(np3-1)/2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3*w2
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3)/w1)
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               jind2=jind3+(j2+clw2)*dlw1
            ELSE
               jind2=0
            END IF
            DO j1=-ih1,ih1
C  first stochastic term
               jind=j1+clw1+1+jind2
               z1=j1
               lwght(jind)=lkern(2,(z1*z1+z2)/hakt2)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(thnew,bi,nv,nvd,n1,n2,n3,hakt2,theta,invcov,
C$OMP& ih3,lwght,wght,y,swjy,mask,thpatch,biipatch,nph1,nph2,nph3,
C$OMP& invcovp,bin)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,aws,n12,
C$OMP& spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12,w1,w2)
C$OMP& PRIVATE(i1,i2,i3,iind,bii,biinv,swj,spmb,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,z,thrednr,ip1,ip2,ip3,pc,ipind,pcj,
C$OMP& jp1,jp2,jp3,jpind,sijp,l,m)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         if(.not.mask(iind)) CYCLE
!$         thrednr = omp_get_thread_num()+1
C returns value in 0:(ncores-1)
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
         pc=0
         DO ip1=i1-nph1,i1+nph1
            if(ip1.le.0.or.ip1.gt.n1) CYCLE
            DO ip2=i2-nph2,i2+nph2
               if(ip2.le.0.or.ip2.gt.n2) CYCLE
               DO ip3=i3-nph3,i3+nph3
                  if(ip3.le.0.or.ip3.gt.n3) CYCLE
                  ipind=ip1+(ip2-1)*n1+(ip3-1)*n12
                  if(.not.mask(ipind)) CYCLE
                  pc=pc+1
                  DO k=1,nv
                     thpatch(k,pc,thrednr)=theta(k,ipind)
                  END DO
                  DO k=1,nvd
                     invcovp(k,pc,thrednr)=invcov(k,ipind)
                  END DO
                  biipatch(pc,thrednr)=bi(ipind)/lambda
               END DO
            END DO
         END DO

C   scaling of sij outside the loop
         swj=0.d0
         DO k=1,nv
            swjy(k,thrednr)=0.d0
         END DO
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3*w2
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3)/w1)
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2*w1
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  if(.not.mask(jind)) CYCLE
                  wj=lwght(jw1+clw1+1+jwind2)
                  IF (aws) THEN
                     pcj=0
                     sij=0.d0
                     DO ip1=i1-nph1,i1+nph1
                        if(ip1.le.0.or.ip1.gt.n1) CYCLE
                        jp1=ip1+jw1
                        if(jp1.le.0.or.jp1.gt.n1) CYCLE
                        DO ip2=i2-nph2,i2+nph2
                           if(ip2.le.0.or.ip2.gt.n2) CYCLE
                           jp2=ip2+jw2
                           if(jp2.le.0.or.jp2.gt.n2) CYCLE
                           DO ip3=i3-nph3,i3+nph3
                              if(sij.gt.1.d0) CYCLE
                              if(ip3.le.0.or.ip3.gt.n3) CYCLE
                              ipind=ip1+(ip2-1)*n1+(ip3-1)*n12
                              if(.not.mask(ipind)) CYCLE
                              jp3=ip3+jw3
                              if(jp3.le.0.or.jp3.gt.n3) CYCLE
                              jpind=jp1+(jp2-1)*n1+(jp3-1)*n12
                              pcj=pcj+1
                              if(.not.mask(jpind)) CYCLE
C   need both ipind and jpind in mask, pcj only depends on i
                              sijp=KLdistsi(theta(1,jpind),
     1                                thpatch(1,pcj,thrednr),
     2                                invcovp(1,pcj,thrednr),nv)
                              sij=max(sij,biipatch(pcj,thrednr)*sijp)
                           END DO
                        END DO
                     END DO
                     IF (sij.ge.1.d0) CYCLE
                     IF (sij.gt.spmin) wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
                  swj=swj+wj
                  DO k=1,nv
                     swjy(k,thrednr)=swjy(k,thrednr)+wj*y(k,jind)
                  END DO
               END DO
            END DO
         END DO
         DO k=1,nv
            thnew(k,iind)=swjy(k,thrednr)/swj
         END DO
         bin(iind)=swj
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(thnew,bin)
      RETURN
      END
