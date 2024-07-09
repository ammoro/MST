*****************************************************************
*     LS K matrix program
*     E.F. Redish, K. Stricker-Bauer (1987)
*     Formalism described in: Phys. Rev. C36, 513 (1987)
*****************************************************************

******************************************************************
*    f90 version 
*    R. Crespo, A.M. Moro (2001)
******************************************************************

      subroutine k(ecm,ifmst)
      use kl
      use parms
      use trace
      use ampnl
      use toff
      use ampold
!      use constants
      implicit real*8 (a-h,o-z),integer*4(i-n)
c      integer nang
!      integer::lmu
      complex*16 :: det
      logical :: first,ifmst
      character*10 nnpot

      namelist /system/ tcm,jmax,nnpot,pm,tm !,icase
      namelist /pw/ ll,is,jj

! AMoro (4/9/07): Remove am0 and amtau
!      namelist /out/ am0,amtau,ampp,ampn,amon,dels,mon,moff,am0,
      namelist /out/ ampp,ampn,amon,dels,mon,moff,am0,
     &               aeon,aeoff,xsec
c
c      Program to calculate partial wave K matrix by matrix inversion
c      of the Lippmann-Schwinger equation. The block data file 'dat.for'
c      the library file 'lib.for', and a 'potslj' file such as 'potrsc.for'
c      or 'potpar.for' must be linked in.
c
c        The partial wave K matrix for a given slj is stored
c        in triangular form in (k,k').  tkl contains the k values
c        with the (k,k') index given by n=(i-1)*(2*k1-i)/2+j where (i,j)
c        are the row and column numbers for the matrix.


       ndel =0
 2223 format(i1)
      delfac=57.295780/(1.+ndel*56.295780)
      amu=931.49432d0
      ch=197.33
      hbarc=ch
 
! Set some default values
       first=.true.
       tcm=0.
       jmax=-1
       nnxsec=0
       rewind(10)
       read(10,nml=system) 
       read(10,nml=out)

       if (tcm<1.e-5) then
          if (.not.ifmst) then
          write(*,*) '**WARNING **:NN program:tcm not found in namelist'
c          write(*,*) ' Using value from MSO input tcm=',ecm 
          endif
          tcm=ecm
       endif
 

c AMoro 30/09/04
      if (pm*tm < 1e-3) then
      	pm=1.0078*amu !938.93
        tm=1.0078*amu !938.93
        write(*,'(" - Nucleon masses not specified=>assuming",$)')
        write(*,'(" pm=",1f6.2," tm=",1f6.2, " MeV")')pm,tm
      else
        pm=pm*amu
        tm=tm*amu
        write(*,*)'- From NN input: pm=',pm/amu,' tm=',tm/amu,' amu'
      endif
      redm=pm*tm/(pm+tm)
      rm=redm
      ttof=ch*ch*pi/rm
      	

c     Phase-shifts
       if (dels.ne.0) then
          open(unit=20,file='dels.dat',status='unknown')
          write(20,2020) tcm*2.
 2020 format(' *********************************************'/
     & '   Lab Energy =',f6.1,' MeV'/)
       endif
c

      p=dsqrt(2.0*redm*tcm)/ch
      write(*,1030) tcm,p
1030  format(' - NNAMP: Using tcm=',1f6.2,' p=',1f6.4)

c
c          choose points and weights
c
       call qw
c
c          input parameters for s, l, j state
c
      if (dels.ne.0) write(20,1020)
 1020 format(' T  s l l" j',2x,' deltal',3x,' deltalp',3x,
     & ' epsilon',7x,' det'/)

      write(14,*)kmax

! Preread to calculate lmax
      lmu=0
91    if (jmax<=0) then
         read(10,nml=pw)
!         write(*,nml=pw)
         if (is<0)goto 93
         if (ll.gt.lmu) lmu=ll
         goto 91
      else
         lmu=jmax
      endif
93    rewind(10)   
      lmu=lmu+1
      write(*,fmt='(" - Reading NN partial waves up to l=",1i2)')lmu

      nmx=kmax*(kmax+1)/2
       write(99,'("+Allocating",1i4,"x",1i4," elements for d01")') 
     &          lmu,nmx
       allocate(d01(lmu,nmx))
       allocate(dz1(lmu,nmx))
       allocate(dm1(lmu,nmx))
      write(99,'("+Allocating",1i4,"x",1i4," elements for dmm1")') 
     &          lmu,kmax
       allocate(dmm1(lmu,kmax,kmax))
       allocate(dpp1(lmu,kmax,kmax))
       allocate(dp1(lmu,nmx))
       d01=0d0;dz1=0d0;dm1=0d0;dmm1=0d0;dpp1=0d0;dp1=0d0
      
100   if (jmax.lt.0) then
         read(10,nml=pw) 
      else 
         if (first) then
            is=0; ll=0; jj=0
            first=.false.
         else if ((jj.eq.jmax).and.(ll.eq.jj+1)) then
            is=-1; ll=0; jj=0
         else if ((ll.eq.jj).and.(is.eq.1)) then
            ll=jj+1
            if (ll.gt.jmax) is=-1
         else if ((ll.eq.jj).and.(is.eq.0)) then
            is=1; jj=jj; ll=jj
            if (jj.eq.0) ll=is
         else if (ll.eq.jj+1) then
            is=0; jj=jj+1; ll=jj
         endif
      endif

 101  format(3i2)
      it=0
      if ((-1)**(is+ll).eq.1) it=1
      ic=1
      if (jj.eq.(ll+1)) ll=ll+2
      if (jj.eq.(ll-1).and.jj.gt.0) ic=2

      k1=kmax+1
      k1c=ic*k1
      nmax=k1c*(k1c+1)/2
 
      if (allocated(tvl)) then
         deallocate (tvl)
         deallocate (tkl)
      endif
     
c *** Allocate memory for arrays
      allocate (tvl(nmax))
      allocate (tkl(nmax))

c
c        Calculate k-sub-l
c
      call klkkp(det,nnpot)
      deallocate(vl)

! IF INFORMATION ON THE K-MATRIX IS REQUIRED, THEN UNCOMMENT 
! THE FOLLOWING LINES

!!$c
!!$c        write k-sub-l`s on file for construction of t(k,q)
!!$c
!!$      write(12,1000) kmax,is,ll,jj,it,ic
!!$
!!$c *** Write k-matrix      
!!$       write(12,1010) p,redm,(qq(i),wq(i),i=1,k1)
!!$       write(12,1010)(tvl(n),n=1,nmax)
!!$       write(12,1010)(tkl(n),n=1,nmax)

c *** K-matrix -> T-matrix 
       call kt !(lmu)
       if (is.lt.0) goto 99



 1000 format(6i5)
 1010 format(8e14.6)
     
c
c          print out phase shifts as a check
c
       
      if (dels.ne.0) then
      alf=2.0*redm*p/ch
      if(ic.eq.2)go to 10
      delta=atan(-alf*tkl(1))
      delta=delta*delfac
      if (dels.ne.0) write(20,1021) it,is,ll,jj,delta,det
 1021 format(i2,2x,2i2,2x,i2,3x,f9.5,21x,2e10.2)
      endif      
      go to 100
 10   continue 
      
      if (dels.ne.0) then
      nod=k1+1
      np=k1*(3*k1+1)/2+1
      delkl=tkl(np)-tkl(1)
      sumkl=tkl(1)+tkl(np)
      ej2=atan(2.0*tkl(nod)/delkl)
      sinej2=sin(ej2)
      cosej2=cos(ej2)
      delta=atan(-alf/2.0*(sumkl-delkl/cosej2))
      deltap=atan(-alf/2.0*(sumkl+delkl/cosej2))
c      delta=delta*delfac
c      deltap=deltap*delfac
      llp=ll-2
c     print 5
c     print 4,is,ll,llp,jj,delta,deltap,sinej2
      dsum=deltap+delta
      ddif=deltap-delta
      sine2b=sinej2*sin(ddif)
      ej2bar=asin(sine2b)
      dbdif=asin(tan(ej2bar)/tan(ej2))
      dbsum=dsum
      deltb=(dbsum-dbdif)/2.0
      deltbp=(dbsum+dbdif)/2.0
      deltb=deltb*delfac
      deltbp=deltbp*delfac
      ejbar=ej2bar*delfac/2.
      write(20,1022) it,is,ll,llp,jj,deltb,deltbp,ejbar,det
      endif
      go to 100
 1022 format(i2,2x,4i2,3x,3(f9.5,1x),2e10.2)
c 5    format(' eigenphase shifts:')
c 6    format(' nuclear bar:')
 99    continue
 1033  format(' -1 0 0 0 0 0')
      end

c===================================================================
      subroutine qw
c====================================================================
c          setup gauss points and weights for integral equation.
c          weights include all factors in integral except v and t.
c          the integration region 0 to infinity is divided into 4
c          intervals.  the second is chosen to be symmetric about the
c          on shell point p and must have an even number of points.
c          the 4th region extends to infinity and is treated by
c          gauss laguerre.  the first point is the on shell point and
c          has weight zero.
c
      use toff
      use pw
      use parms
      use trace
      use kl
      implicit double precision (a-h,o-z),integer*4(i-n)
! AMoro (4/9/07): Remove b01 and qlo and use default b01=qlo=0
!      namelist /kmat/ hw,b01,b34,n1,n2,n3,n4,qlo,qhi,s,vnn
      namelist /kmat/ hw,b01,b34,n1,n2,n3,n4,qlo,qhi,s,vnn



      alf=2.0*redm*p/ch
!          set up limits of 4 intervals and the number of points in each
c      read(10,*) hw
c      read(10,*) b01,b34
c      read(10,*) n1,n2,n3,n4
c      read(10,*) s
c      read(10,*) qlo,qhi

!      Defalult values
       hw=0.5
       qlo=0.
       b01=0.
       qhi=1000.
       s=1.
       vnn=.false.

       read(10,nml=kmat)
      if (b01.ne.0.) then
      write(*,*)'** WARNING ** ignoring b01 from input and using b01=0' 
      b01=0
      endif

!     *********JAT Addition...ic undefined here***********************
      write(99,fmt='("qw:n1,n2,n3,n4,ic,kmax",6i3)')n1,n2,n3,n4,ic,kmax
      ic=0
c     nnn=ic*(n1+n2+n3+n4+1)

      if (dels.ne.0) write(20,1) p,hw
 1    format(' pcm =',f6.3,' f-1'/' half width of pv interval=',f6.3)
      b12=p-hw
      b23=p+hw
      if (dels.ne.0) write(20,3) n1,n2,n3,n4
 3    format(' number of points in the intervals =',4i4)
      if(b12.lt.b01)n1=0.
      if(b23.gt.b34)n3=0.
      if(n1.le.0)b12=b01
      if(n3.le.0)b23=b34
      if(hw.le.0)n2=0
      if(n2/2*2.ne.n2)n2=n2+1
      if (dels.ne.0) write(20,4) b12,b23,b34
 4    format(' boundary points of the regions:  0.',3f10.5)
 25   kmax=n1+n2+n3+n4
      k1=kmax+1
      
      write(99,2000) k1 
 2000 format('+Allocating',i3,' points for qq and wq...')
      allocate (qq(k1))
      allocate (wq(k1))
     
c          scale points and weights appropriately for the 4 intervals
      ni=1
c          the on shell point
      qq(1)=p
      wq(1)=0.0
      if(n1.eq.0)go to 30
c          interval 1
      do 35 i=1,n1
      ni=ni+1
      c1=(b12-b01)/2.0
      c2=(b12+b01)/2.0
      qq(ni)=pwg(n1+1,i)*c1+c2
      wq(ni)=pwg(npwg-n1,npwg-i+1)*c1
 35   continue
 30   if(n2.eq.0)go to 40
c          interval 2
      do 45 i=1,n2
      ni=ni+1
      c1=(b23-b12)/2.0
      c2=(b23+b12)/2.0
      qq(ni)=pwg(n2+1,i)*c1+c2
      wq(ni)=pwg(npwg-n2,npwg-i+1)*c1
 45   continue
 40   if(n3.eq.0)go to 50
c          interval 3
      do 55 i=1,n3
      ni=ni+1
      c1=(b34-b23)/2.0
      c2=(b34+b23)/2.0
      qq(ni)=pwg(n3+1,i)*c1+c2
      wq(ni)=pwg(npwg-n3,npwg-i+1)*c1
 55   continue
 50   if(n4.eq.0)go to 60
c          interval 4
      if (dels.ne.0) write(20,5) s
 5    format(' scale for gauss laguerre quadrature',f5.1)
      do 65 i=1,n4
      ni=ni+1
      z=pwgl(n4+1,i)
      x=z/s+b34
      qq(ni)=x
      wq(ni)=pwgl(npwgl-n4,npwgl-i+1)/s
 65   continue
 60   continue
c          multiply weights by other factors in integral
      if (dels.ne.0) write(20,7) qlo,qhi
 7    format(' low and high momentum cutoffs=',2f8.1/)
      do 70 i=2,k1
      q=qq(i)
      wq(i)=wq(i)*2.0*alf*q*q/pi/p/(p*p-q*q)
      if(q.lt.qlo)wq(i)=0.0
      if(q.gt.qhi)wq(i)=0.0
 70   continue
      go to 99

 1000 format(i5)
 1010 format(6e13.7)
 99   write(99,*)'Exiting qw subroutine..'
      return
      end
c************************************************************
      subroutine klkkp(det,pot)
        use kl
        use toff
        use wk
c
c          routine to calculate k matrix elements
c
      implicit real*8 (a-h,o-z),integer*4(i-n)
      complex*16 det
      character*10 :: pot

      write(99,*)'Entering klkkp subroutine...'

      epsck=1.e-4
c          choose points and weights
      k1=kmax+1
      k1c=ic*k1
c      mc=2*k1c
      
      write(99,2010) k1c,k1c
 2010 format('+Allocating ',i3,'x',i3,' elements for cvl,ak,bk')
      allocate(cvl(k1c,k1c))
      allocate(ak(k1c,k1c))
      allocate(bk(k1c,k1c))

c          calculate v for slj state
      if (pot.eq.'paris') then
!         print*,'- Using Paris potential...'
         call potslj
      else if (pot.eq.'bonn') then
!         print*,'- Using Bonn potential...'
         call potsljbonn
      else 
         print*,'**ERROR** Specified potential',pot,'not defined'
         stop
      endif


      do 50 i=1,k1c
      do 50 j=1,k1c
  50  cvl(i,j)=vl(i,j)
c          set up equation for matrix inversion routine
     

      do 200 i=1,k1c
      do 200 j=1,k1c
      if(j.le.k1)wwq=wq(j)
      if(j.gt.k1)wwq=wq(j-k1)
      ak(i,j)=-wwq*cvl(i,j)
      if(i.eq.j)ak(i,j)=1.0+ak(i,j)
 200  continue
c          invert matrix
      do 300 i=1,k1c
      do 300 j=1,k1c
  300 bk(i,j)=cvl(i,j)
c     call leq(ak,bk,k1c,k1c,kmx,kmx,det)
      call leq(k1c,k1c,k1c,k1c,det)
c          store k`s in triangular form
      do 100 i=1,k1c
      do 100 j=i,k1c
      n=(i-1)*(2*k1c-i)/2+j
      if (vnn) then
         tkl(n)=vl(i,j) !approximates T by NN interaction
      else
         tkl(n)=bk(i,j)
        
      endif
      tvl(n)=vl(i,j)
!      write(8,*) n,' ',tkl(n) 
c          check to be sure k is symmetric
      if(abs(bk(i,j)-bk(j,i)).gt.epsck) then
         write(*,1)i,j,bk(i,j),bk(j,i)
      endif
 100  continue
 1    format(1x,' i=',i3,2x,' j=',i3,2x,' k(i,j)=',2e12.5,2x,' k(j,i)=',
     1 2e12.5)
      
      write(99,*)'-Deallocating memory for cvl,ak,bk'
      deallocate(cvl)
      deallocate(ak)
      deallocate(bk)
      return
      end

! This is now in a module (AMoro)!
c*************************************************************
c     Gauss point data for K matrix program
c      block data
c      use pw
c      implicit double precision (a-h,o-z),integer*4(i-n)
c      parameter (npwg=21)
c      parameter (npwgl=7)
c      common/pw/pwg(npwg,npwg),pwgl(npwgl,npwgl)

c     library routines for LS programs
c**************************************************************

      subroutine ql(xs,l,qql)
c
c          subroutine to calculate the legendre function of the second
c          kind ql(x) for given x and l.
c
      implicit real*8 (a-h,o-z),integer*4(i-n)
c      real*8 xs,qql
      parameter (numx=21)
      parameter (nnmx=21)
      parameter (nlmx=2*nnmx+numx)
      integer :: i
      integer,save::init
      real*8:: fl(nlmx),dl(nlmx)
      real*8:: al2,pol(numx),dos
      dos=2d0
      al2=dlog(dos)
!AMORO 10/2/04
!      common/bfl/ fl,dl !!NOT USED????
!      data init /0/
      x=xs
      l1=l+1
      nn1=numx
      if(x.gt.5)nn1=11
      if(x.gt.10)nn1=6
      if(nn1.le.nnmx)go to 110
      print 2
 2    format(' error in ql: nnmx too small')
      return
 110  if(l1.le.numx) go to 100
      print 1
 1    format(' error in ql: l too large')
      return
 100  if(init.ne.0)go to 200
      init=1
c          fl contains the log of the factorial function, dl
c          the double factorial.
!      al2=dlog(dos)
     
      fl(1)=0d0
      dl(1)=0d0
      do 10 i=1,nlmx-1
         fl(i+1)=fl(i)+dlog(dfloat(i))
         dl(i+1)=dl(i)+dlog(dfloat(2*i+1))
 10   continue
          
 200  continue
      if(x.le.2.)goto 50
c  for x>2 use the hypergeometric series to calculate ql10 and q
      qol=0d0
      do 30 i=1,nn1
         nn=i-1
         test=dble(2*nn+l1)*dlog(x)
         if(test.gt.100.0)go to 30
         if (2*nn+l+1>nlmx) then
            write(*,*)'ql: argument of fl too large. Aborting...'
            stop
         endif
         a1=dexp(fl(2*nn+l+1)-dble(nn)*al2-fl(nn+1)-dl(l+nn+1))
         qol=qol+a1/(x**(2*nn+l+1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
!      write(94,1002)x,i,a1,al2,qol
!1002  format(1f12.6,1i3,3f12.6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 30   continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      write(95,1001)x,l,qol
!1001  format(1f12.6,1i3,1f12.6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      go to 99
c          for x<2 use the equation for ql in terms of the legendre
c          polynomials pol(x), calculated in pldp (see a&s 8.6.19).
 50   call pldp(x,l1,pol)
      qol=pol(l1)*dlog((x+1d0)/(x-1d0))/2d0
      if(l.eq.0) go to 99
      do 70 i=1,l
      qol=qol-pol(i)*pol(l1-i)/float(i)
 70   continue
 99   qql=qol
      return
      end
c***************************************************************
      subroutine pldp(x,l,pol)
c         computes legendre polynominals of order zero to order
c         (l) for a given argument (x=cos(theta))
c         also valid for x>1
      implicit real*8 (a-h,o-z),integer*4(i-n)
      real*8:: s,x,pol(*)
      pol(1)=1.0
      pol(2)=x
      if(l.le.1)goto20
      do 10 jj=2,l
      s=jj
   10 pol(jj+1)=((2.0*s-1.0)*x*pol(jj)-(s-1.0)*pol(jj-1))/s
   20 return
      end


c***************************************************
      function phasef(n)
      implicit real*8(a-h,o-z),integer*4(i-n)
      phasef = float(1-2*iabs(n-2*(n/2)))
      return
      end
c***************************************************
      function yxfct(m,n)
        use factorials
        implicit none
        integer::m,n,numax,nu
        real*8 :: yxfct,cntrl,fctor
c     computes nfact/mfact
      yxfct = 1.0
      numax = m-n
      if(numax.gt.40) go to 6
      if(numax) 2,7,1
    1 cntrl = -1.0
      fctor =float(n+1)
      go to 3
    2 cntrl = 1.0
      numax = -numax
      fctor =float(m+1)
    3 continue
      do 4 nu=1,numax
      yxfct = yxfct*fctor
      fctor = fctor+1.0
    4 continue
      if(cntrl) 5,7,7
    5 yxfct = 1.0/yxfct
      go to 7
    6 continue
!AMORO (27/11/2003)
      fctor = flog(n+1)-flog(m+1)
!     fctor = faclog(n+1)-faclog(m+1)
      yxfct = exp(fctor)
    7 continue
      return
      end
c*****************************************************************
      subroutine leq (n,m,la,lb,det)
        use wk
        use jat

      implicit complex*16(a-h,o-z),integer*4(i-n)
      real*8 t1, t2
c     leq         linear equations subroutine
c      solves the matrix equation ax=b
c       a=coefficient matrix
c       b=constant matrix, solution stored in b.
c       n=number of equations and unknowns
c       m=number of sets of equations
c       la=dimension of first subscript of a
c       lb=dimension of first subscript of b

!      real*8 vl
!       parameter (kmx=100)
!      common/wk/vl(kmx,kmx),ak(kmx,kmx),bk(kmx,kmx),cvl(kmx,kmx)
!      common/jat/ a(kmx*kmx),b(kmx*kmx)

      write(99,*)'Entering leq...'
      allocate(a(n*n))
      allocate(b(n*n))
      
      i=1
c      do 777 jt1=1,kmx
c      do 778 jt2=1,kmx
      do 777 jt1=1,n
      do 778 jt2=1,n
      a(i)=ak(jt2,jt1)
      b(i)=bk(jt2,jt1)
      i=i+1
  778 continue
  777 continue

      lna=la*n-la
      lmb=lb*m-lb

      do 70 i=2,n
      ii=i-1
      lj=-la
      do 70 j=1,ii
      lj=lj+la
      lij=i+lj
      ljj=j+lj
      t1=a(lij)*conjg(a(lij))
c      write(99,*)'t1=',t1
      if (t1.eq.0.0) go to 70
      t2=a(ljj)*conjg(a(ljj))
c      write(99,*)'t2=',t2
      if (t2.lt.t1) go to 10
      go to 40
   10 lnao=lna+1

      do 20 lko=1,lnao,la
      lk=lko-1
      ljk=j+lk
      lik=i+lk
      t=a(ljk)
      a(ljk)=a(lik)
   20 a(lik)=-t
      lmbo=lmb+1
      do 30 lko=1,lmbo,lb
      lk=lko-1
      ljk=j+lk
      lik=i+lk
      t=b(ljk)
      b(ljk)=b(lik)
   30 b(lik)=-t
   40 r=a(lij)/a(ljj)
      llj=lj+la
      do 50 lk=llj,lna,la
      lik=i+lk
      ljk=j+lk
   50 a(lik)=a(lik)-r*a(ljk)
      lmbo=lmb+1
      do 60 lko=1,lmbo,lb
      lk=lko-1
      lik=i+lk
      ljk=j+lk
   60 b(lik)=b(lik)-r*b(ljk)
   70 continue
      lmbo=lmb+1
      do 90 ljo=1,lmbo,lb
      lj=ljo-1
      lnj=n+lj
      lnn=n+lna
      b(lnj)=b(lnj)/a(lnn)
      li=lna
      do 90 ii=2,n
      li=li-la
      i=n-ii+1
      kk=i+1
      t=0.0
      lk=la*i-la
      do 80 k=kk,n
      lk=lk+la
      lik=i+lk
      lkj=k+lj
   80 t=t+a(lik)*b(lkj)
      lij=i+lj
      lii=i+li
   90 b(lij)=(b(lij)-t)/a(lii)
      det=(1.0,0.0)
      do 100 j=1,n
      ljj=j+(j-1)*la
  100 det=det*a(ljj)
      i=1
c      do 887 jt1=1,kmx
c      do 888 jt2=1,kmx
      do 887 jt1=1,n
      do 888 jt2=1,n
      ak(jt2,jt1)=a(i)
      bk(jt2,jt1)=b(i)
      i=i+1
  888 continue
  887 continue
      write(99,*)'-Deallocating memory for a,b'
      deallocate(a)
      deallocate(b)
      write(99,*)'Leaving leq...'
      return
      end


c********************************************************
      function plm1(l,cth,pol,lmx)
      implicit double precision (a-h,o-z),integer*4(i-n)
      dimension pol(lmx)
      plm1=0.
      if(abs(cth).ge.1.0) return
      if(l.eq.0)return
      sth=sqrt(1.-cth*cth)
      if(abs(sth).lt.1.e-3)go to 10
      plm1=l*(cth*pol(l+1)-pol(l))/sth
      return
 10   if(sth.lt.0.)sth=0.
      plm1=-l*(l+1.)*sth/2.
      if(cth.lt.0.)plm1=-plm1*(1.-2.*mod(l,2))
      return
      end
c*******************************************************
      function plm1p(l,cth,pol,lmx)
      implicit double precision (a-h,o-z),integer*4(i-n)
      dimension pol(lmx)
      plm1p=0.
      if(l.eq.0)return
      plm1p=-l*(l+1.)/2.
      if(cth.lt.0.)plm1p=-plm1p*(1.-2.*mod(l,2))
      if(abs(cth).ge.1.0) return
      sth=sqrt(1.-cth*cth)
      if(abs(sth).lt.1.e-3)return
      plm1p=l*(cth*pol(l+1)-pol(l))/sth/sth
      return
      end
c********************************************************
      subroutine wrt2(nx,ny,dx,dy,p,th,f,ndim,n)
c
c          writes a 2 dimensional matrix on file n
c
      implicit double precision (a-h,o-z),integer*4(i-n)
      complex*16 f
      dimension f(ndim,ndim)
      write(n,1000)nx,ny,dx,dy,p,th
      do 1 i=1,ndim
 1    write(n,*) (f(j,i),j=1,ndim)
      return
 1000 format(2i5,4e14.7)
      end
c*******************************************************
      function nnn(m,n,i,j,k1c)
c
c     this function is used to untangle
c     the k matrix for coupled states
c     produced by KLS
c
      implicit double precision (a-h,o-z),integer*4(i-n)
      k=m+n
      go to (10,20,30,40,10) k
 20   nnn=i*(2*k1c-i-1)/2+j+1
      return
 30   if(m.gt.n) go to 35
      nnn=i*(2*k1c-i-1)/2+j+k1c/2+1
      return
 35   nnn=j*(2*k1c-j-1)/2+i+k1c/2+1
      return
 40   nnn=(k1c/2+i)*(3*k1c/2-i-1)/2+k1c/2+j+1
      return
 10   print 9
 9    format(' error in nnn')
      return
      end
 
*************************************************************
      function fint2d(xtab,ytab,fxytab,xbar,ybar,nnx,nny,nord,mmx)
c
c          2 dimension aitkin interpolation routine.  note that mesh points
c          xtab, ytab must be in increasing order.
c
      implicit double precision (a-h,o-z),integer*4(i-n)
      dimension xytab(2,100),x(10),y(10),
     1 xybar(2),nbg(2),xtab(mmx),ytab(mmx),nxy(2)
      dimension fxytab(mmx,mmx),fxy(10,10)
      nnn=nord+1
c          set up arrays for loop
      do 40 j=1,nnx
 40   xytab(1,j)=xtab(j)
      do 45 j=1,nny
 45   xytab(2,j)=ytab(j)
      nxy(1)=nnx
      nxy(2)=nny
      xybar(1)=xbar
c          begin loop to determine the (order+1) x and y points over
c          which interpolation is made; 1=x, 2=y.
      xybar(2)=ybar
      do 10 i=1,2
30    num=1
      if(xybar(i).lt.xytab(i,1))go to 85
      num=nxy(i)-nnn+1
      if(xybar(i).gt.xytab(i,nxy(i)))go to 85
50     min=1
      max=nxy(i)
      num=nxy(i)/2
55    if(max-min.lt.4)goto70
      if(xytab(i,num)-xybar(i))60,82,61
60     min=num
          num=num+(max-min+1)/2
      goto55
  61   max=num
      num=(max-min+1)/2
      goto55
  70   num=max
  71  if(xytab(i,num)-xybar(i))82,82,81
 81   num=num-1
      goto71
 82   num=max0(1,num-nnn/2)
      num=min0(num,nxy(i)-nnn+1)
 85   nbg(i)=num
 10   continue
c          end loop.  set up x, y, function arrays of required
c          (order+1)**2 points
      nx=nbg(1)
      ny=nbg(2)
      do 20 i=1,nnn
      x(i)=xytab(1,nx+i-1)
      y(i)=xytab(2,ny+i-1)
      do 20 j=1,nnn
      fxy(i,j)=fxytab(nx+i-1,ny+j-1)
 20   continue
c          do interpolation
 90   do95ii=2,nnn
      do95 jj=ii,nnn
      do 95 mm=ii,nnn
      fxy(jj,mm)=((fxy(ii-1,ii-1)*(x(jj)-xybar(1))-fxy(jj,ii-1)*(x(ii-1)
     1 -xybar(1)))*(y(mm)-xybar(2))-(fxy(ii-1,mm)*(x(jj)-xybar(1))
     2 -fxy(jj,mm)*(x(ii-1)-xybar(1)))*(y(ii-1)-xybar(2)))
     3 /(x(jj)-x(ii-1))/(y(mm)-y(ii-1))
  95   continue
      fint2d=fxy(nnn,nnn)
      return
      end

 
