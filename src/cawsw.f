CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Perform one iteration in local constant three-variate aws (gridded)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cawsw2io(y,fix,n1,n2,n3,hakt,hhom,lambda,theta,bi,bi2,
     1                bi0,ai,model,kern,spmin,lwght,inwght,outwght)
C   compute weights for all direct neighbors of design points
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
C
C   components of wght:
C     1: x-1  2: x+1  3: y-1  4: y+1  5: z-1  6: z+1
C
      implicit none

      external kldist,lkern
      double precision kldist,lkern
      integer n1,n2,n3,model,kern
      logical fix(*)
      double precision y(*),theta(*),bi(*),lambda,hhom(*),
     1       inwght(6,n1,n2,n3),outwght(6,n1,n2,n3),hakt,lwght(*),
     2       spmin,spf,bi0(*),ai(*),bi2(*),swj,swj2,swjy,swj0
      integer ih1,ih2,ih3,i1,i2,i3,j1,j2,j3,jw1,jw2,jw3,jwind3,jwind2,
     1        iind,jind,jind3,jind2,clw1,clw2,clw3,dlw1,dlw2,dlw3,
     2        dlw12,n12,ind
      double precision thetai,bii,sij,z1,z2,z3,wj,hakt2,hmax2,hhomi,
     1        hhommax,w1,w2,w3,hdist
      hakt2=hakt*hakt
      spf=1.d0/(1.d0-spmin)
C
C   first calculate location weights
C
      ih3=FLOOR(hakt)
      ih2=ih3
      ih1=ih3
      if(n3.eq.1) ih3=0
      if(n2.eq.1) ih2=0
      clw1=ih1
      clw2=ih2
      clw3=ih3
      dlw1=ih1+clw1+1
      dlw2=ih2+clw2+1
      dlw3=ih3+clw3+1
      dlw12=dlw1*dlw2
      n12=n1*n2
      z2=0.d0
      z3=0.d0
      hmax2=0.d0
      DO j3=-clw3,clw3
         if(n3.gt.1) THEN
            z3=j3
            z3=z3*z3
            ih2=FLOOR(sqrt(hakt2-z3))
            jind3=(j3+clw3)*dlw12
         ELSE
            jind3=0
         END IF
         DO j2=-ih2,ih2
            if(n2.gt.1) THEN
               z2=j2
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
               lwght(jind)=lkern(kern,(z1*z1+z2)/hakt2)
               if(lwght(jind).gt.0.d0) hmax2=max(hmax2,z2+z1*z1)
            END DO
         END DO
      END DO
      call rchkusr()
C$OMP PARALLEL DEFAULT(NONE)
C$OMP& SHARED(y,fix,bi,n1,n2,n3,hakt2,hmax2,theta,hhom,
C$OMP& ih3,lwght,outwght,inwght,bi0,ai,bi2)
C$OMP& FIRSTPRIVATE(ih1,ih2,lambda,n12,
C$OMP& model,spmin,spf,dlw1,clw1,dlw2,clw2,dlw3,clw3,dlw12)
C$OMP& PRIVATE(i1,i2,i3,iind,thetai,bii,swj,swj2,swjy,swj0,
C$OMP& sij,wj,j3,jw3,jind3,z3,jwind3,j2,jw2,jind2,z2,jwind2,
C$OMP& j1,jw1,jind,z1,ind,hhomi,hhommax,w1,w2,w3,hdist)
C$OMP DO SCHEDULE(GUIDED)
      DO iind=1,n1*n2*n3
         hhomi=hhom(iind)
         hhomi=hhomi*hhomi
         hhommax=hmax2
         IF (fix(iind)) CYCLE
         i1=mod(iind,n1)
         if(i1.eq.0) i1=n1
         i2=mod((iind-i1)/n1+1,n2)
         if(i2.eq.0) i2=n2
         i3=(iind-i1-(i2-1)*n1)/n12+1
         thetai=theta(iind)
         bii=bi(iind)/lambda
         swj=0.d0
         swj2=0.d0
         swj0=0.d0
         swjy=0.d0
C   scaling of sij outside the loop
         DO jw3=-clw3,clw3
            j3=jw3+i3
            if(j3.lt.1.or.j3.gt.n3) CYCLE
            jwind3=(jw3+clw3)*dlw12
            jind3=(j3-1)*n12
            z3=jw3
            z3=z3*z3
            if(n2.gt.1) ih2=FLOOR(sqrt(hakt2-z3))
            DO jw2=-ih2,ih2
               j2=jw2+i2
               if(j2.lt.1.or.j2.gt.n2) CYCLE
               jwind2=jwind3+(jw2+clw2)*dlw1
               jind2=(j2-1)*n1+jind3
               z2=jw2
               z2=z3+z2*z2
               ih1=FLOOR(sqrt(hakt2-z2))
               DO jw1=-ih1,ih1
C  first stochastic term
                  j1=jw1+i1
                  if(j1.lt.1.or.j1.gt.n1) CYCLE
                  jind=j1+jind2
                  wj=lwght(jw1+clw1+1+jwind2)
                  sij=bii*kldist(model,thetai,theta(jind))
                  IF (sij.gt.1.d0) CYCLE
C
C   use local structure of weights from previous step to modify the location kernel
C
                  if(jw1.lt.0) w1 = inwght(1,i1,i2,i3)
                  if(jw1.gt.0) w1 = inwght(2,i1,i2,i3)
                  if(jw2.lt.0) w2 = inwght(3,i1,i2,i3)
                  if(jw2.gt.0) w2 = inwght(4,i1,i2,i3)
                  if(jw3.lt.0) w3 = inwght(5,i1,i2,i3)
                  if(jw3.gt.0) w3 = inwght(6,i1,i2,i3)
                  hdist = abs(jw1)+abs(jw2)+abs(jw3)
                  wj = dmin1(wj,(w1*abs(jw1)+w2*abs(jw2)+w3*abs(jw3))/
     1                              sqrt(hdist))
                  IF (sij.gt.spmin) THEN
                     wj=wj*(1.d0-spf*(sij-spmin))
                  END IF
C
C   store local weights
C
                  if(abs(i1-j1)+abs(i2-j2)+abs(i3-j3).eq.1) THEN
                     ind = (i1-j1)+2*(i2-j2)+3*(i3-j3)
                     SELECT CASE (ind)
                        CASE(-1)
                           outwght(1,i1,i2,i3) = wj
                        CASE(1)
                           outwght(2,i1,i2,i3) = wj
                        CASE(-2)
                           outwght(3,i1,i2,i3) = wj
                        CASE(2)
                           outwght(4,i1,i2,i3) = wj
                        CASE(-3)
                           outwght(5,i1,i2,i3) = wj
                        CASE(3)
                           outwght(6,i1,i2,i3) = wj
                        CASE DEFAULT
                           ind = 0
                     END SELECT
                  END IF
                  swj=swj+wj
                  swj2=swj2+wj*wj
                  swjy=swjy+wj*y(jind)
               END DO
            END DO
         END DO
         ai(iind)=swjy
         bi(iind)=swj
         bi2(iind)=swj2
         bi0(iind)=swj0
         hhom(iind)=sqrt(hhommax)
      END DO
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C$OMP FLUSH(outwght)
      RETURN
      END
