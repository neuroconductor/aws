      subroutine initflob(y,n1,n2,h0,ahat,chat,gamma,Ni,Bmat,bvec,
     1                    Svec,Si,Pvec,Qmat)
      implicit logical (a-z)
      Integer n1,n2
      real*8 y(n1,n2),h0,ahat(n1,n2),chat(n1,n2),gamma(2,n1,n2),
     1       Ni(n1,n2),Pvec(2,n1,n2),Qmat(3,n1,n2),Bmat(3,n1,n2),
     2       bvec(2,n1,n2),Svec(2,n1,n2),Si(n1,n2)
      integer i1,i2
      real*8 z1,z2,ch,chsq,ah
      DO i1=1,n1
         DO i2=1,n2
            call aggrej20(y,n1,n2,i1,i2,h0,Bmat(1,i1,i2),
     1                    bvec(1,i1,i2),Svec(1,i1,i2),Si(i1,i2),
     2                    Ni(i1,i2))
            call estajcj0(Bmat(1,i1,i2),bvec(1,i1,i2),Svec(1,i1,i2),
     1                    Si(i1,i2),Ni(i1,i2),gamma(1,i1,i2),
     1                    ahat(i1,i2),chat(i1,i2))
            z1=gamma(1,i1,i2)
            z2=gamma(2,i1,i2)
            ch=chat(i1,i2)
            chsq=ch*ch
            Qmat(1,i1,i2)=chsq*Bmat(1,i1,i2)
            Qmat(2,i1,i2)=chsq*Bmat(2,i1,i2)
            Qmat(3,i1,i2)=chsq*Bmat(3,i1,i2)
            ah=ahat(i1,i2)
            Pvec(1,i1,i2)=ch*(Svec(1,i1,i2)-ah*bvec(1,i1,i2))
            Pvec(2,i1,i2)=ch*(Svec(2,i1,i2)-ah*bvec(2,i1,i2))
         END DO
      END DO
      RETURN 
      END
      subroutine initflow(y,n1,n2,h0,ahat,chat,gamma,Ni,Pvec,Qmat)
      implicit logical (a-z)
      Integer n1,n2
      real*8 y(n1,n2),h0,ahat(n1,n2),chat(n1,n2),gamma(2,n1,n2),
     1       Ni(n1,n2),Pvec(2,n1,n2),Qmat(3,n1,n2),Bmat(3),
     2       bvec(2),Svec(2)
      integer i1,i2
      real*8 Si,z1,z2,ch,chsq,ah
      DO i1=1,n1
         DO i2=1,n2
           call aggrej20(y,n1,n2,i1,i2,h0,Bmat,bvec,Svec,Si,Ni(i1,i2))
            call estajcj0(Bmat,bvec,Svec,Si,Ni(i1,i2),gamma(1,i1,i2),
     1                    ahat(i1,i2),chat(i1,i2))
            z1=gamma(1,i1,i2)
            z2=gamma(2,i1,i2)
            ch=chat(i1,i2)
            chsq=ch*ch
            Qmat(1,i1,i2)=chsq*Bmat(1)
            Qmat(2,i1,i2)=chsq*Bmat(2)
            Qmat(3,i1,i2)=chsq*Bmat(3)
            ah=ahat(i1,i2)
            Pvec(1,i1,i2)=ch*(Svec(1)-ah*bvec(1))
            Pvec(2,i1,i2)=ch*(Svec(2)-ah*bvec(2))
         END DO
      END DO
      RETURN 
      END
      subroutine aggrej20(y,n1,n2,j1,j2,h0,Bjmat,bjvec,Sjvec,Sj,Nj)
C  compute aggregated statistics in pixel j
      implicit logical (a-z)
      integer n1,n2,j1,j2
      real*8 y(n1,n2),h0,Bjmat(3),bjvec(2),Sjvec(2),Sj,Nj
      integer k1a,k1e,k2a,k2e,k1,k2
      real*8 x1,x2,wkj,h0sq
      h0sq=h0*h0
      k2=h0
      k2a=max(1,j2-k2)
      k2e=min(n2,j2+k2)
      Bjmat(1)=0.d0
      Bjmat(2)=0.d0
      Bjmat(3)=0.d0
      bjvec(1)=0.d0
      bjvec(2)=0.d0
      Sjvec(1)=0.d0
      Sjvec(2)=0.d0
      Sj=0.d0
      Nj=0.d0
      DO k2=k2a,k2e
         x2=k2-j2
         k1a=x2-h0+j1
         k1a=max(1,k1a)
         k1e=x2+h0+j1
         k1e=min(n1,k1e)
         DO k1=k1a,k1e
            x1=k1-j1
            wkj=max(0.d0,1.d0-(x1*x1+x2*x2)/h0sq)
            if(wkj.le.0.d0) CYCLE
            Bjmat(1)=Bjmat(1)+x1*x1*wkj
            Bjmat(2)=Bjmat(2)+x1*x2*wkj
            Bjmat(3)=Bjmat(3)+x2*x2*wkj
            bjvec(1)=bjvec(1)+x1*wkj
            bjvec(2)=bjvec(2)+x2*wkj
            Sjvec(1)=Sjvec(1)+y(k1,k2)*x1*wkj
            Sjvec(2)=Sjvec(2)+y(k1,k2)*x2*wkj
            Sj=Sj+y(k1,k2)*wkj
            Nj=Nj+wkj
         END DO
      END DO
      RETURN
      END
      subroutine estajcj0(Bmat,bvec,Svec,S,N,gamma,a,c)
C   estimate a_j and c_j
      implicit logical (a-z)
      real*8 Bmat(3),bvec(2),Svec(2),S,N,gamma(2),a,c
      real*8 mat(3,3),vec(3),z1,z2
      integer info
      mat(1,1)=N
      mat(1,2)=bvec(1)
      mat(1,3)=bvec(2)
      mat(2,2)=Bmat(1)
      mat(2,3)=Bmat(2)
      mat(3,3)=Bmat(3)
      vec(1)=S
      vec(2)=Svec(1)
      vec(3)=Svec(2)
C     now calculate theta as mat^{-1} vec
      call dposv("U",3,1,mat,3,vec,3,info)
C    if info>0 just keep the old estimate
      IF (info.gt.0) THEN
         call intpr("singular in estajcj0",20,info,1)
      ELSE
         a=vec(1)
         z1=vec(2)
         z2=vec(3)
         c=sqrt(z1*z1+z2*z2)
         if(c.gt.1d-10) THEN
            gamma(1)=z1/c
            gamma(2)=z2/c
         ELSE
            gamma(1)=1
            gamma(2)=0
         END IF
      END IF
      RETURN
      END
      
      real*8 function wtj2(dx1,dx2,a)
C   compute anisotropic local weights wtj for point t and pixel j
C   a contains 1/h^2 I + 1/eta^2 gamma_j gamma_j^T
C
      implicit logical (a-z)
      real*8 a(3)
      integer dx1,dx2
      real*8 z11,z22,z12
      z11=dx1*dx1
      z12=dx1*dx2
      z22=dx2*dx2
      wtj2=max(0.d0,1.d0-z11*a(1)-2.d0*z12*a(2)-z22*a(3))
      RETURN
      END
      subroutine aggregj2(y,n1,n2,j1,j2,gammaj,gamma,ngamma,Njold,
     1                    Qmat,Pvec,hsqinv,etasqinv,lamsig,
     1                    Bjmat,bjvec,Sjvec,Sj,Nj,mr2,sel)
C  compute aggregated statistics in pixel j
      implicit logical (a-z)
      integer n1,n2,j1,j2,sel
      real*8 y(n1,n2),gammaj(2),gamma(2,n1,n2),Njold,Qmat(3),Pvec(2),
     1       hsqinv,etasqinv,Bjmat(3),bjvec(2),Sjvec(2),Sj,Nj,lamsig,
     2       ngamma(n1,n2)
      integer k1a,k1e,k2a,k2e,k1,k2,dx1,dx2
      real*8 wtj2,vijst,nlsig,QigplPig,QuaQgam,Pvecgam,mr2
      real*8 a(3),z1,z2,x1,x2,det,wkj,dkj,vijstb
      external vijst,vijstb
      nlsig=lamsig*Njold
      QigplPig=-QuaQgam(gamma(1,j1,j2),Qmat)+
     1                 2.d0*Pvecgam(gamma(1,j1,j2),Pvec)
C      call dblepr("QigplPig",8,QigplPig,1)
C      call dblepr("nlsig",5,nlsig,1)
      z1=gammaj(1)
      z2=gammaj(2)
      a(1)=hsqinv+etasqinv*z1*z1
      a(2)=etasqinv*z1*z2
      a(3)=hsqinv+etasqinv*z2*z2
      det=(a(3)*a(1)-a(2)*a(2))
      k2=sqrt(a(1)/det)
      k2a=max(1,j2-k2)
      k2e=min(n2,j2+k2)
      Bjmat(1)=0.d0
      Bjmat(2)=0.d0
      Bjmat(3)=0.d0
      bjvec(1)=0.d0
      bjvec(2)=0.d0
      Sjvec(1)=0.d0
      Sjvec(2)=0.d0
      Sj=0.d0
      Nj=0.d0
      DO k2=k2a,k2e
         dx2=k2-j2
         x2=dx2
         z1=sqrt(a(1)-det*x2*x2)/a(1)
         z2=-x2*a(2)/a(1)
         k1a=z2-z1+j1
         k1a=max(1,k1a)
         k1e=z2+z1+j1
         k1e=min(n1,k1e)
         DO k1=k1a,k1e
            dx1=k1-j1
            x1=dx1
            dkj=ngamma(j1,j2)/ngamma(k1,k2)
            wkj=wtj2(dx1,dx2,a)
            if(x1*x1+x2*x2.ge.mr2) THEN
               if(sel.eq.1) THEN
               wkj=wkj*vijst(gamma(1,k1,k2),Qmat,
     1                      Pvec,QigplPig,nlsig,dkj)
               ELSE
               wkj=wkj*vijstb(gamma(1,k1,k2),Qmat,
     1                      Pvec,QigplPig,nlsig,dkj)
               END IF
C               wkj=wkj*sijac(Qmat,j1,j2,k1,k2,
C     1                             ahat(i1,i2),chat(i1,i2),
C     2                             gamma(1,i1,i2),ngamma(i1,i2),
C     3                             ahat(j1,j2),chat(j1,j2),
C     4                             gamma(1,j1,j2),ngamma(j1,j2),nlsig)
               if(wkj.le.0.d0) CYCLE
            END IF
            Bjmat(1)=Bjmat(1)+x1*x1*wkj
            Bjmat(2)=Bjmat(2)+x1*x2*wkj
            Bjmat(3)=Bjmat(3)+x2*x2*wkj
            bjvec(1)=bjvec(1)+x1*wkj
            bjvec(2)=bjvec(2)+x2*wkj
            Sjvec(1)=Sjvec(1)+y(k1,k2)*x1*wkj
            Sjvec(2)=Sjvec(2)+y(k1,k2)*x2*wkj
            Sj=Sj+y(k1,k2)*wkj
            Nj=Nj+wkj
         END DO
      END DO
      RETURN
      END
      subroutine aggreg2b(y,n1,n2,j1,j2,ahat,chat,gamma,ngamma,hsqinv,
     1                    etasqinv,lamsig,Bjmat,bjvec,Sjvec,Sj,Nj,mr2)
C  compute aggregated statistics in pixel j
      implicit logical (a-z)
      integer n1,n2,j1,j2,sel
      real*8 y(n1,n2),gamma(2,n1,n2),hsqinv,etasqinv,Bjmat(3),
     1       bjvec(2),Sjvec(2),Sj,Nj,lamsig,ngamma(n1,n2),
     2       ahat(n1,n2),chat(n1,n2)
      integer k1a,k1e,k2a,k2e,k1,k2,dx1,dx2
      real*8 wtj2,mr2,wijst,sijaci
      real*8 a(3),z1,z2,x1,x2,det,wkj,dkj,sjj,gamj(2)
      real*8 No,Bomat(3),bovec(2),Sovec(2),So
      external wijst,sijaci
      gamj(1)=gamma(1,j1,j2)/ngamma(j1,j2)
      gamj(2)=gamma(2,j1,j2)/ngamma(j1,j2)
      z1=gamj(1)
      z2=gamj(2)
      a(1)=hsqinv+etasqinv*z1*z1
      a(2)=etasqinv*z1*z2
      a(3)=hsqinv+etasqinv*z2*z2
      det=(a(3)*a(1)-a(2)*a(2))
      k2=sqrt(a(1)/det)
      k2a=max(1,j2-k2)
      k2e=min(n2,j2+k2)
C  first copy old content of sum statistics
      No=Nj
      So=Sj
      Bomat(1)=Bjmat(1)
      Bomat(2)=Bjmat(2)
      Bomat(3)=Bjmat(3)
      bovec(1)=bjvec(1)
      bovec(2)=bjvec(2)
      Sovec(1)=Sjvec(1)
      Sovec(2)=Sjvec(2)
      Bjmat(1)=0.d0
      Bjmat(2)=0.d0
      Bjmat(3)=0.d0
      bjvec(1)=0.d0
      bjvec(2)=0.d0
      Sjvec(1)=0.d0
      Sjvec(2)=0.d0
      Sj=0.d0
      Nj=0.d0
      sjj=sijaci(ahat(j1,j2),chat(j1,j2),gamj,So,Sovec,No,Bomat,bovec)
      DO k2=k2a,k2e
         dx2=k2-j2
         x2=dx2
         z1=sqrt(a(1)-det*x2*x2)/a(1)
         z2=-x2*a(2)/a(1)
         k1a=z2-z1+j1
         k1a=max(1,k1a)
         k1e=z2+z1+j1
         k1e=min(n1,k1e)
         DO k1=k1a,k1e
            dx1=k1-j1
            x1=dx1
            dkj=ngamma(j1,j2)/ngamma(k1,k2)
            wkj=wtj2(dx1,dx2,a)
            if(x1*x1+x2*x2.ge.mr2) THEN
               gamj(1)=gamma(1,k1,k2)/ngamma(k1,k2)
               gamj(2)=gamma(2,k1,k2)/ngamma(k1,k2)
               wkj=wkj*wijst(ahat(k1,k2),chat(k1,k2),x1,x2,
     1                       gamj,so,Sovec,No,Bomat,bovec,sjj,lamsig)
               if(wkj.le.0.d0) CYCLE
            END IF
            Bjmat(1)=Bjmat(1)+x1*x1*wkj
            Bjmat(2)=Bjmat(2)+x1*x2*wkj
            Bjmat(3)=Bjmat(3)+x2*x2*wkj
            bjvec(1)=bjvec(1)+x1*wkj
            bjvec(2)=bjvec(2)+x2*wkj
            Sjvec(1)=Sjvec(1)+y(k1,k2)*x1*wkj
            Sjvec(2)=Sjvec(2)+y(k1,k2)*x2*wkj
            Sj=Sj+y(k1,k2)*wkj
            Nj=Nj+wkj
         END DO
      END DO
      RETURN
      END
      real*8 function sijaci(a,c,gamma,s,Svec,N,Bmat,bvec)
      implicit logical (a-z)
      real*8 a,c,gamma(2),s,Svec(2),N,Bmat(3),bvec(2)
      real*8 z,cg1,cg2
      cg1=c*gamma(1)
      cg2=c*gamma(2)
      z=cg1*(a*bvec(1)-Svec(1))+cg2*(a*bvec(2)-Svec(2))-a*s
      z=2.d0*z+a*a*N+cg1*cg1*Bmat(1)+2*cg1*cg2*Bmat(2)+cg2*cg2*Bmat(3)
      sijaci=z
      RETURN
      END
      real*8 function sijacij(a,c,x1,x2,gamma,s,Svec,N,Bmat,bvec)
      implicit logical (a-z)
      real*8 x1,x2
      real*8 a,c,gamma(2),s,Svec(2),N,Bmat(3),bvec(2)
      real*8 z,cg1,cg2,aij
      cg1=c*gamma(1)
      cg2=c*gamma(2)
      aij=a-cg1*x1-cg2*x2
      z=cg1*(aij*bvec(1)-Svec(1))+cg2*(aij*bvec(2)-Svec(2))-aij*s
      z=2.d0*z+aij*aij*N
      sijacij=z+cg1*cg1*Bmat(1)+2*cg1*cg2*Bmat(2)+cg2*cg2*Bmat(3)
      RETURN
      END
      real*8 function wijst(a,c,x1,x2,gamma,s,Svec,N,Bmat,bvec,sii,
     1                      lamsig)
      implicit logical (a-z)
      real*8 x1,x2
      real*8 a,c,gamma(2),s,Svec(2),N,Bmat(3),bvec(2),sii,lamsig
      real*8 sijacij
      external sijacij
      real*8 z
      z=sijacij(a,c,x1,x2,gamma,s,Svec,N,Bmat,bvec)-sii
      if(z.gt.lamsig) THEN
         wijst=0.d0
      ELSE
         z=z/lamsig
         if(z.le.0.25) THEN
            wijst=1.d0
         ELSE
            wijst=4.d0/3.d0*(1.d0-z)
         END IF
      END IF
      RETURN
      END
      subroutine estacjb(Bjmat,bjvec,Sjvec,Sj,Nj,gamma,ngamma,aj,cj)
C   estimate a_j and c_j
      implicit logical (a-z)
      real*8 Bjmat(3),bjvec(2),Sjvec(2),Sj,Nj,gamma(2),ngamma,aj,cj
      real*8 a11,a12,a22,b1,b2,d,cg1,cg2
      cg1=gamma(1)/ngamma
      cg2=gamma(2)/ngamma
      a11=Nj
      a12=bjvec(1)*cg1+bjvec(2)*cg2
      a22=cg1*Bjmat(1)*cg1+2.d0*cg1*Bjmat(2)*cg2+cg2*Bjmat(3)*cg2
      b1=Sj
      b2=Sjvec(1)*cg1+Sjvec(2)*cg2
      d=a11*a22-a12*a12
      IF(d>1.e-10) THEN
         aj=(a22*b1-a12*b2)/d
         cj=(a11*b2-a12*b1)/d
         if(cj.lt.0.d0) THEN 
            cj=-cj
            gamma(1)=-gamma(1)
            gamma(2)=-gamma(2)
         END IF
C         cj=abs(a22*b2-a12*b1)/d
C   restict to positive c to keep information on direction
      ELSE
C         call dblepr("non-positive determinant",24,d,1)
         aj=Sj/Nj
         cj=0.d0
      END IF
      RETURN
      END
      subroutine estajcj(Bjmat,bjvec,Sjvec,Sj,Nj,gamma,aj,cj)
C   estimate a_j and c_j
      implicit logical (a-z)
      real*8 Bjmat(3),bjvec(2),Sjvec(2),Sj,Nj,gamma(2),ngamma,aj,cj
      real*8 a11,a12,a22,b1,b2,d,cg1,cg2
      cg1=gamma(1)/ngamma
      cg2=gamma(2)/ngamma
      a11=Nj
      a12=bjvec(1)*cg1+bjvec(2)*cg2
      a22=cg1*Bjmat(1)*cg1+2.d0*cg1*Bjmat(2)*cg2+cg2*Bjmat(3)*cg2
      b1=Sj
      b2=Sjvec(1)*cg1+Sjvec(2)*cg2
      d=a11*a22-a12*a12
      IF(d>1.e-10) THEN
         aj=(a22*b1-a12*b2)/d
         cj=abs(a11*b2-a12*b1)/d
C   restict to positive c to keep information on direction
      ELSE
C         call dblepr("non-positive determinant",24,d,1)
         aj=Sj/Nj
         cj=0.d0
      END IF
      RETURN
      END
      subroutine itfstep1(y,n1,n2,gamma,Niold,Qmat,Pvec,h,eta,
     1           sigma2,lambda,Bmat,bvec,Svec,Ni,chat,ahat,ngamma,mr2,
     2           sel)
C   flow estimation step 1
      implicit logical (a-z)
      integer n1,n2,sel
      real*8 y(n1,n2),gamma(2,n1,n2),Niold(n1,n2),Qmat(3,n1,n2),
     1       Pvec(2,n1,n2),h,eta,Bmat(3,n1,n2),bvec(2,n1,n2),
     2       Svec(2,n1,n2),Ni(n1,n2),chat(n1,n2),ahat(n1,n2),si,
     3       lambda,sigma2,ngamma(n1,n2),mr2
      integer i1,i2
      real*8 z1,z2,sgamma(2),hsqinv,etasqinv,lamsig
      lamsig=lambda*sigma2
      hsqinv=1.d0/h/h
      etasqinv=1.d0/eta/eta
      DO i1=1,n1
         DO i2=1,n2
            sgamma(1)=gamma(1,i1,i2)/ngamma(i1,i2)
            sgamma(2)=gamma(2,i1,i2)/ngamma(i1,i2)
C compute statistics
               call aggregj2(y,n1,n2,i1,i2,sgamma,gamma,ngamma,
     1                       Niold(i1,i2),Qmat(1,i1,i2),Pvec(1,i1,i2),
     2                       hsqinv,etasqinv,lamsig,Bmat(1,i1,i2),
     3                       bvec(1,i1,i2),Svec(1,i1,i2),si,Ni(i1,i2),
     4                       mr2,sel)
C compute estimates ahat and chat
            call estajcj(Bmat(1,i1,i2),bvec(1,i1,i2),Svec(1,i1,i2),si,
     1                   Ni(i1,i2),sgamma,ahat(i1,i2),chat(i1,i2))
         END DO
      END DO
      RETURN
      END
      subroutine itstep1b(y,n1,n2,gamma,hinv,etainv,sigma2,lambda,
     1                    Bmat,bvec,Svec,Ni,si,chat,ahat,ngamma,mr2)
C   flow estimation step 1
      implicit logical (a-z)
      integer n1,n2,sel
      real*8 y(n1,n2),gamma(2,n1,n2),ngamma(n1,n2),chat(n1,n2),
     1       ahat(n1,n2),hinv,etainv,Bmat(3,n1,n2),bvec(2,n1,n2),
     2       Svec(2,n1,n2),Ni(n1,n2),si(n1,n2),lambda,sigma2,mr2
      integer i1,i2
      real*8 z1,z2,hsqinv,etasqinv,lamsig
      lamsig=lambda*sigma2
      hsqinv=hinv*hinv
      etasqinv=etainv*etainv
      DO i1=1,n1
         DO i2=1,n2
C compute sum statistics, need statistics at (i1,i2) from last step 
C for the penalty
C need full content of ahat, chat, gamma from step (k-1)
            call aggreg2b(y,n1,n2,i1,i2,ahat,chat,gamma,ngamma,hsqinv,
     1                    etasqinv,lamsig,Bmat(1,i1,i2),bvec(1,i1,i2),
     2                    Svec(1,i1,i2),si(i1,i2),Ni(i1,i2),mr2)
C now sum statistics from step (k-1) have been replaced
C ahat, chat, gamma unchanged
         END DO
      END DO
      DO i1=1,n1
         DO i2=1,n2
C compute estimates ahat and chat
            call estacjb(Bmat(1,i1,i2),bvec(1,i1,i2),Svec(1,i1,i2),
     1             si(i1,i2),Ni(i1,i2),gamma(1,i1,i2),ngamma(i1,i2),
     2             ahat(i1,i2),chat(i1,i2))
C estimates ahat and chat replaced
         END DO
      END DO
      RETURN
      END
      subroutine itfstep2(gamma,n1,n2,h,sigma2,lambda,Pvec,Qmat,Niold,
     1                   chat,ahat,Bmat,bvec,Svec,gammanew,
     2                   ngamma,mr2)
      implicit logical (a-z)
      integer n1,n2
      real*8 gamma(2,n1,n2),h,sigma2,lambda,Pvec(2,n1,n2),
     1       Qmat(3,n1,n2),Niold(n1,n2),
     2       chat(n1,n2),ahat(n1,n2),Bmat(3,n1,n2),bvec(2,n1,n2),
     3       Svec(2,n1,n2),gammanew(2,n1,n2),ngamma(n1,n2),mr2
      integer i1,i2,j1,j2,j1a,j1e,j2a,j2e
      external Pvecgam,QuaQgam,vijst,vijstb,sijac
      real*8 Pvecgam,QuaQgam,vijst,vijstb,sijac
      real*8 vij,QigplPig,lamsig,nlsig,Qi1,Qi2,Qi3,Pi1,Pi2,z1,z2,
     1       h2,cj,cjsq,det,nPi,dji,z
      integer hj
C   flow estimation step 2
      lamsig=lambda*sigma2
      h2=h*h
      DO i1=1,n1
         DO i2=1,n2
            nlsig=lamsig*Niold(i1,i2)
            QigplPig=-QuaQgam(gamma(1,i1,i2),Qmat(1,i1,i2))+
     1               2.d0*Pvecgam(gamma(1,i1,i2),Pvec(1,i1,i2))
            Qi1=0.d0
            Qi2=0.d0
            Qi3=0.d0
            Pi1=0.d0
            Pi2=0.d0
            hj=h
            j1a=max(1,i1-hj)
            j1e=min(n1,i1+hj)
            DO j1=j1a,j1e
               z1=j1-i1
               hj=sqrt(h*h-z1*z1)
               j2a=max(1,i2-hj)
               j2e=min(n2,i2+hj)
               DO j2=j2a,j2e
                  z2=j2-i2
                  vij=max(0.d0,1.d0-(z1*z1+z2*z2)/h2)
                  if(vij.le.0.d0) CYCLE
C                  if(z1*z1+z2*z2.gt.mr2) THEN
                     dji=ngamma(i1,i2)/ngamma(j1,j2)
                     vij=vij*vijst(gamma(1,j1,j2),Qmat(1,i1,i2),
     1                             Pvec(1,i1,i2),QigplPig,nlsig,dji)
                     if(vij.le.0.d0) CYCLE
                  cj=abs(chat(j1,j2))
                  cjsq=cj*cj
                  Qi1=Qi1+cjsq*Bmat(1,j1,j2)*vij
                  Qi2=Qi2+cjsq*Bmat(2,j1,j2)*vij
                  Qi3=Qi3+cjsq*Bmat(3,j1,j2)*vij
              Pi1=Pi1+cj*(Svec(1,j1,j2)-ahat(j1,j2)*bvec(1,j1,j2))*vij
              Pi2=Pi2+cj*(Svec(2,j1,j2)-ahat(j1,j2)*bvec(2,j1,j2))*vij
               END DO
            END DO
            det=Qi1*Qi3-Qi2*Qi2
            nPi=Pi1*Pi1+Pi2*Pi2
            if(det.gt.1e-20.and.nPi.gt.1e-10) THEN
               gammanew(1,i1,i2)=(Qi3*Pi1-Qi2*Pi2)/det
               gammanew(2,i1,i2)=(Qi1*Pi2-Qi2*Pi1)/det
               Qmat(1,i1,i2)=Qi1
               Qmat(2,i1,i2)=Qi2
               Qmat(3,i1,i2)=Qi3
               Pvec(1,i1,i2)=Pi1
               Pvec(2,i1,i2)=Pi2
            ELSE
C   keep the old estimate
C               call dblepr("det",3,det,1)
C               call dblepr("nPi",3,nPi,1)
               gammanew(1,i1,i2)=gamma(1,i1,i2)
               gammanew(2,i1,i2)=gamma(2,i1,i2)
            END IF
         END DO
      END DO
C   recompute ngamma
      DO i1=1,n1
         DO i2=1,n2
               z1 = gammanew(1,i1,i2)
               z2 = gammanew(2,i1,i2)
               z = sqrt(z1*z1+z2*z2)
               ngamma(i1,i2) = z
         END DO
      END DO
      RETURN
      END
      real*8 function QuaQgam(gamma,Qmat)
      implicit logical(a-z)
      real*8 gamma(2),Qmat(3),z1,z2
C  compute gamma^T Qmat gamma
      z1=gamma(1)
      z2=gamma(2)
      QuaQgam=Qmat(1)*z1*z1+2.d0*Qmat(2)*z1*z2+Qmat(3)*z2*z2
      RETURN
      END
      real*8 function Pvecgam(gamma,Pvec)
      implicit logical(a-z)
      real*8 gamma(2),Pvec(2)
C  compute Pvec^T gamma
      Pvecgam=Pvec(1)*gamma(1)+Pvec(2)*gamma(2)
      RETURN
      END
      real*8 function vijst(gammaj,Qmat,Pvec,QigplPig,nlsig,dkj)
C   compute K_st(sij) with plateau kernel and plateau width 0.25 
      implicit logical(a-z)
      real*8 gammaj(2),Qmat(3),Pvec(2),QigplPig,nlsig,dkj
      external Pvecgam,QuaQgam
      real*8 Pvecgam,QuaQgam
      real*8 z,w1,w2
      z=QigplPig+dkj*dkj*QuaQgam(gammaj,Qmat)-2.d0*dkj*
     1                           (Pvecgam(gammaj,Pvec))
      if(z.gt.nlsig) THEN
         vijst=0.d0
      ELSE
         z=z/nlsig
         if(z.le.0.25) THEN
            vijst=1.d0
         ELSE
            vijst=4.d0/3.d0*(1.d0-z)
         END IF
      END IF
      RETURN
      END
      real*8 function vijstb(gammaj,Qmat,Pvec,QigplPig,nlsig,dkj)
C   compute K_st(sij) with plateau kernel and plateau width 0.25 
      implicit logical(a-z)
      real*8 gammaj(2),Qmat(3),Pvec(2),QigplPig,nlsig,dkj
      external Pvecgam,QuaQgam
      real*8 Pvecgam,QuaQgam
      real*8 z,w1,w2
      z=QigplPig+dkj*dkj*QuaQgam(gammaj,Qmat)-2.d0*dkj*
     1                           abs(Pvecgam(gammaj,Pvec))
      if(z.gt.nlsig) THEN
         vijstb=0.d0
      ELSE
         z=z/nlsig
         if(z.le.0.25) THEN
            vijstb=1.d0
         ELSE
            vijstb=4.d0/3.d0*(1.d0-z)
         END IF
      END IF
      RETURN
      END
      subroutine numgrad2(yx1,yx2,yy1,yy2,d1,d2,grad)
C  compute numeric gradient 
      implicit logical (a-z)
      real*8 yx1,yx2,yy1,yy2,grad(2)
      integer d1,d2
      real*8 dx,dy
      grad(1)=(yx2-yx1)/d1
      grad(2)=(yy2-yy1)/d2
      return 
      end
      subroutine ingamma2(y,n1,n2,gamma)
C  initialize gamma using numerical gradients
      implicit logical (a-z)
      integer n1,n2
      real*8 y(n1,n2),gamma(2,n1,n2)
      integer i,j
C   first corners
      call numgrad2(y(1,1),y(2,1),y(1,1),y(1,2),1,1,gamma(1,1,1))
      call numgrad2(y(n1-1,1),y(n1,1),y(n1,1),y(n1,2),1,1,
     1              gamma(1,n1,1))
      call numgrad2(y(1,n2),y(2,n2),y(1,n2-1),y(1,n2),1,1,
     1              gamma(1,1,n2))
      call numgrad2(y(n1-1,n2),y(n1,n2),y(n1,n2-1),y(n1,n2),1,1,
     1              gamma(1,n1,n2))
C   now edges
      DO i=2,n1-1
         call numgrad2(y(i-1,1),y(i+1,1),y(i,1),y(i,2),2,1,
     1                 gamma(1,i,1))
         call numgrad2(y(i-1,n2),y(i+1,n2),y(i,n2-1),y(i,n2),2,1,
     1                 gamma(1,i,n2))
      END DO
      DO j=2,n2-1
         call numgrad2(y(1,j),y(2,j),y(1,j-1),y(1,j+1),1,2,
     1                 gamma(1,1,j))
         call numgrad2(y(n1-1,j),y(n1,j),y(n1,j-1),y(n1,j+1),1,2,
     1                 gamma(1,n1,j))
      END DO
C   now the interior part
      DO i=2,n1-1
         DO j=2,n2-1
            call numgrad2(y(i-1,j),y(i+1,j),y(i,j-1),y(i,j+1),2,2,
     1                    gamma(1,i,j))
         END DO
      END DO
      RETURN
      END
