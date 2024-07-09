*****************************************************************
*    J. Tostevin (July 1995)
*
*    Formalism and f90 version: A.M. Moro, R. Crespo (July 2001): 
c                                  Phys. Rev. C 65, 054001 (2002)
*****************************************************************

      subroutine maboff()
*****************************************************************
*   Routine to generate a(q,k)-e(q,k)
*
*   Main changes:
*   - kmx, nmx,lnomx,nmx,np parameters no longer used 
*   - common blocks replaced by modules 
*****************************************************************
      use ampold
      use ampnew
      use nnamp
      use j
      use rc
      use jat1
      use parms
      use toff
      use ampnl
      use nnamp
      USE trace
!      use constants

      implicit none
      integer :: itj,isj,ikj,iqj,isgpj,isgj,iaj,ibj,iqp
      real*8 :: sjt,cjt1,stat,rkj,cjt2,qj,sgpj,sgj,cleb,aj,bj,u9,cjt3

      integer :: i,kth,kkmax,nqmax,jq,ik,nmax,
     $            n,m,kn,ir,kr,iprt,ith,
     $            knp,norder

       real*8 :: rm,dth,eps,angle,dmat,
     $       fkq,sthkq,cthkq,th,xq,xkk,xk,xkp,cth,
     $       x0,w0,sth,epsrc,alpha,beta,thqk,norm
     
      real*8, allocatable :: thpr(:)
      complex*16 :: aaon,bbon,ccon,ddon,eeon
      complex*16 :: aa2on,bb2on,cc2on,dd2on,ee2on
      complex*16 :: aa3on,bb3on,cc3on,dd3on,ee3on
      complex*16 :: aa4on,bb4on,cc4on,dd4on,ee4on
      complex*16 :: zi,cint2db
      complex*16 :: m1a,m2a,m3a,m4a,m5a,m6a,m1b,m2b,m3b,m4b,m5b,m6b
      real*8, dimension(2) :: fisp(2,2)
      data fisp/0.25,-0.25,0.75,0.25/

c *** Interfaces:  required to pass dummy arrays as arguments
       interface
         subroutine newgd(lmu,kmax,xk,xkp,dz1,dz)
           integer:: lmu,kmax
           real*8 :: xk, xkp
           complex*16, dimension(:,:), intent(in) :: dz1
           complex*16, dimension(:), intent(out) ::  dz
         end subroutine newgd
       end interface

       interface
         subroutine newgnd(lmu,kmax,xk,xkp,dpp1,dpp)
         integer:: lmu,kmax
         real*8 :: xk, xkp
         complex*16, dimension(:,:,:), intent(in) :: dpp1
         complex*16, dimension(:), intent(out) ::  dpp
        end subroutine newgnd
       end interface
c *** ----------------- end interfaces ------------------------------

       complex*16 :: maux1,maux2
       complex*16 ,dimension(2) :: m00,m11,m1m1,m10,m01,mss
       complex*16 , allocatable :: mst(:,:,:,:,:,:),mkq(:,:,:,:,:,:)
       complex*16 , allocatable :: mhj(:,:,:,:,:,:)
       complex*16 , allocatable :: mfoff(:,:,:,:,:,:,:)
       complex*16 , allocatable :: mfHJ(:,:,:,:,:,:,:) !Hooton-Johnson

       write(99,*)'Entering moff with lmu=',lmu
       write(*,*) '- Calculating off-shell amplitudes M(a,b)...'

!       open(unit=14,file='tout.dat',status='unknown')
!       rewind(14)
      open(unit=60,file='m0000.off',status='unknown')
      open(unit=61,file='m0111.off',status='unknown')
      open(unit=62,file='m1100.off',status='unknown')
      open(unit=63,file='m1120.off',status='unknown')
      open(unit=64,file='m1121.off',status='unknown')
      open(unit=65,file='m1122.off',status='unknown')

      open(unit=66,file='mHJ0000.off',status='unknown')
      open(unit=67,file='mHJ0111.off',status='unknown')
      open(unit=68,file='mHJ1100.off',status='unknown')
      open(unit=69,file='mHJ1120.off',status='unknown')
      open(unit=70,file='mHJ1121.off',status='unknown')
      open(unit=71,file='mHJ1122.off',status='unknown')

c      open(unit=66,file='m1121.off',status='unknown')
c      open(unit=67,file='m1122.off',status='unknown')

      zi=cmplx(0.0,1.0)
      eps=1.e-3
c
c       ifkq = 0 for q=(k"-k)/sqrt(2), kk=(k+k")/sqrt(2)
c            = 1 for q=(k"-k)        , kk=(k+k")/2*)
c     xkmax= max val. of kk, xqmax= max val. of q
c     dk   = kk step,        dq   = q step
c     theta= angle bet. q and kk
c     icase = 0 for pp, 1 for pn, 2 for average, 3 for T=0
c     lx = max L in file
c

      if(ifkq.eq.0)fkq=1.
      if(ifkq.eq.1)fkq=2.

!      write(99,*) 'theta,icase,itype,nth=',theta,icase,itype,nth
      write(99,*) 'theta,icase,nth=',theta,icase,nth
    
      call mesh(xkmax,xqmax,dk,dq,kkmax,nqmax)
      allocate(mst(0:1,0:1,-1:1,-1:1,kkmax,nqmax))
      allocate(mkq(0:1,0:1,0:2,-2:2,kkmax,nqmax))
      allocate(mhj(0:1,0:1,0:2,-2:2,kkmax,nqmax)) !Hooton-Johnson
      allocate(mfoff(0:1,0:1,0:1,0:2,-2:2,kkmax,nqmax))
      allocate(mfHJ(0:1,0:1,0:1,0:2,-2:2,kkmax,nqmax))!Hooton-Johnson

      if (.not.allocated(xkm)) then
         write(99,*) '+Allocating memory for xkm,xkmp,ctheta...'
         allocate(xkm(kkmax,nqmax))
         allocate(xkmp(kkmax,nqmax))
         allocate(ctheta(kkmax,nqmax))
      endif

      allocate(thpr(nth))     

      if (allocated(d0)) deallocate(d0,dz,dm,dmm,dpp,dp)
      allocate(d0(lmu))
      allocate(dz(lmu))
      allocate(dm(lmu))
      allocate(dmm(lmu))
      allocate(dpp(lmu))
      allocate(dp(lmu))

c
c         if theta=90deg calculate on shell quantities
c         if theta is 0deg or 180deg set c, c` to zero
c
      kth=0
      if(abs(theta-90.).lt.eps)kth=1
      if(abs(theta).lt.eps.or.abs(theta-180.).lt.eps)kth=-1
      th=theta*pi/180.
      cthkq=cos(th)
      sthkq=sin(th)
 
      dth=180./(nth-1)
      do 33 i=1,nth
 33   thpr(i)=dth*(i-1.)

      if (.not.allocated(pol)) then
         write(99,*) '+Allocating memory for pol, cs'
         allocate(pol(lmu))
         allocate(cs(nth))
      end if
      
c
c          set up k(kk,q), k`(kk,q), pl(kk,q)
c
      do 110 jq=1,nqmax
      xq=xxq(jq)
      do 110 ik=1,kkmax
      xkk=xxk(ik)
      xk=dsqrt(xq*xq/fkq/2.+xkk*xkk*fkq/2.-xkk*xq*cthkq)
      xkp=dsqrt(xq*xq/fkq/2.+xkk*xkk*fkq/2.+xkk*xq*cthkq)
      xkm(ik,jq)=xk
      xkmp(ik,jq)=xkp
      cth=1.
      if(xk.le.0..or.xkp.le.0.)go to 115
      cth=(xkk*xkk*fkq/2.-xq*xq/fkq/2.)/xk/xkp
 115  ctheta(ik,jq) = cth
c  print*,ik,jq,xkk,xq,ctheta(ik,jq)
 110  continue

c
c          read tl`s from file 14
c 
!!$       rewind(14)
!!$       read(14,*)kmax
!!$       if(kmax.lt.0)go to 160
!!$
!!$       nmax=kmax*(kmax+1)/2
!!$       write(99,*) 'Reading T-matrix from unit 14'
!!$       
!!$100    read(14,1000)is,ll,jj,it,ic
!!$1000   format (6i5)
!!$c       write(*,1000)is,ll,jj,it,ic
!!$       if(is.lt.0) go to 160
!!$c
!!$c **       spin zero case
!!$      if(is.eq.0)then
!!$         read(14,1010)(d01(jj+1,kn),kn=1,nmax)
!!$!         write(9,*) 'Reading d01'
!!$!         write(9,1010)(d01(jj+1,kn),kn=1,nmax)
!!$c **       spin one case (ic=1->uncoupled, ic=2->coupled)
!!$      else
!!$        if (ic.eq.1) then
!!$          if (jj.gt.0)then
!!$          read(14,1010)(dz1(jj+1,kn),kn=1,nmax)
!!$          else
!!$          read(14,1010)(dm1(jj+2,kn),kn=1,nmax)
!!$1010      format(8e14.7)
!!$          endif
!!$        else if (ic.eq.2) then
!!$        do 105 m=1,ic
!!$        do 105 n=1,ic
!!$          if(m.eq.1.and.n.eq.1)then
!!$          read(14,1010)(dm1(jj+2,kn),kn=1,nmax)
!!$          else if(m.eq.1.and.n.eq.2)then
!!$          read(14,1010)((dmm1(jj,kn,knp),kn=1,kmax),knp=1,kmax)
!!$          else if(m.eq.2.and.n.eq.1)then
!!$          read(14,1010)((dpp1(jj+2,kn,knp),kn=1,kmax),knp=1,kmax)
!!$!          write(98,*)'jj=',jj
!!$!               do kn=1, kmax
!!$!               do knp=1,kmax
!!$!                  write(98,'(2i3,6g12.4)')kn,knp,dpp1(jj+2,kn,knp)
!!$!               enddo
!!$!            enddo
!!$          else if(m.eq.2.and.n.eq.2)then
!!$          read(14,1010)(dp1(jj,kn),kn=1,nmax)
!!$          endif
!!$ 105    continue
!!$        endif
!!$      endif
!!$      go to 100

c
c          begin loops to interpolate t and v on kk, q mesh,
c          and add together for off-shell a and c
c
 160  continue
     
      write(99,*) 't-matrix read. begin jq,ik loop'

      do 140 jq=1,nqmax
      xq = xxq(jq)
      do 140 ik=1,kkmax
      xkk = xxk(ik)
      xk=xkm(ik,jq)
      xkp=xkmp(ik,jq)
!      thqk = acos(ctheta(ik,jq))
      thqk = dacos(ctheta(ik,jq))
      if (xk.ne.xkp) then
         norm = (xk-xkp)**2 + 4.*xk*xkp*cos(thqk/2.)**2
         norm = 1/norm
         alpha = norm*(xk*sin(thqk))**2
         beta = norm*(xk+xkp*cos(thqk))**2
! angle to convert from Redish coordinates to 
! Hooton-Johnson by rotation around y axis
         angle=dacos(sqrt(beta)) 
      endif

c **      generates all partial wave amplitudes in the new grid
      call newgd(lmu,kmax,xk,xkp,d01,d0)
      call newgd(lmu,kmax,xk,xkp,dz1,dz)
      call newgd(lmu,kmax,xk,xkp,dm1,dm)
      call newgd(lmu,kmax,xk,xkp,dp1,dp)
      call newgnd(lmu,kmax,xk,xkp,dpp1,dpp)
      call newgnd(lmu,kmax,xk,xkp,dmm1,dmm)
           
c **      calculates M00,M01,M1-1,M10,M11,MSS
      call ampcal(thqk,lmu,m11,m10,m1m1,m01,m00,mss)

c     ir=1:T=0, ir=2:T=1

      do 111 ir=1,2
c     ------------------------------------------------------------------
c     assign S=1 amplitudes to mst array, from ampcal
c     ------------------------------------------------------------------
      mst(1,ir-1,1, 1,ik,jq)=m11(ir)
      mst(1,ir-1,1, 0,ik,jq)=m10(ir)
      mst(1,ir-1,1,-1,ik,jq)=m1m1(ir)
      mst(1,ir-1,0, 1,ik,jq)=m01(ir)
      mst(1,ir-1,0, 0,ik,jq)=m00(ir)
      
c     ------------------------------------------------------------------
c     from symmetry relations
c     ------------------------------------------------------------------
      mst(1,ir-1, 0,-1,ik,jq)=-mst(1,ir-1,0, 1,ik,jq)
      mst(1,ir-1,-1, 1,ik,jq)= mst(1,ir-1,1,-1,ik,jq)
      mst(1,ir-1,-1, 0,ik,jq)=-mst(1,ir-1,1, 0,ik,jq)
      mst(1,ir-1,-1,-1,ik,jq)= mst(1,ir-1,1, 1,ik,jq)
c     ------------------------------------------------------------------
c     assign S=0 amplitudes, from ampcal
c     ------------------------------------------------------------------
      mst(0,ir-1, 0, 0,ik,jq)=mss(ir)
c     ------------------------------------------------------------------
  111 continue
!      write(9,*) 'Mkq'
c     ------------------------------------------------------------------
c     translate to Mkq form (first zero array)
c     ------------------------------------------------------------------
      do 171 itj=0,1
      do 171 isj=0,1
      do 171 ikj=0,2,1
      do 171 iqj=-2,2,1
      mkq(isj,itj,ikj,iqj,ik,jq)=(0.d0,0.d0)
  171 continue
c     ------------------------------------------------------------------
c     translate to Mkq form 
c     ------------------------------------------------------------------
      do 113 itj=0,1
      do 113 isj=0,1
      sjt=isj
      cjt1=1.d0/stat(sjt)
      cjt1=cjt1*cjt1
      do 113 ikj=0,2,1
      rkj=ikj
      cjt2=cjt1*stat(rkj)
      do 1116 iqj=-ikj,ikj,1
      qj=iqj
c     ------------------------------------------------------------------
c *** Sum over sigma' and sigma for Mkq'
c     ------------------------------------------------------------------
      do 1115 isgpj=-isj,isj,1
      sgpj=isgpj
      do 1115 isgj=-isj,isj,1
      sgj=isgj
      mkq(isj,itj,ikj,iqj,ik,jq)=mkq(isj,itj,ikj,iqj,ik,jq)+
     +      mst(isj,itj,isgpj,isgj,ik,jq)*cleb(sjt,sgpj,rkj,qj,sjt,sgj)
c     ------------------------------------------------------------------
 1115 continue
      mkq(isj,itj,ikj,iqj,ik,jq)=mkq(isj,itj,ikj,iqj,ik,jq)*cjt2
 1116 continue

c *** Hooton-Johnson
      do 186 iqj=-ikj,ikj,1
      mhj(isj,itj,ikj,iqj,ik,jq)=(0.d0,0.d0)
      do 173 iqp=-2,2,1
      mhj(isj,itj,ikj,iqj,ik,jq)=mhj(isj,itj,ikj,iqj,ik,jq)
     & + dmat(ikj,iqp,iqj,angle)*mkq(isj,itj,ikj,iqp,ik,jq)
173   continue

     
c      if ((ikj.eq.2).and.(iqj.eq.1).and.(isj.eq.1)) then
c      write(98,1024) itj,isj,ikj,iqj,mhj(isj,itj,ikj,iqj,ik,jq)
c      endif
186   continue


1024  format(4i,2e14.6)
  113 continue

      if(jq.eq.3) then
      do 150 isj=0,1
      do 150 itj=0,1
      do 150 ikj=0,2*isj
      do 150 iqj=0,ikj
!      print*,isj,itj,ikj,iqj,mkq(isj,itj,ikj,iqj,ik,jq)
  150 continue
      endif
      
c     ------------------------------------------------------------------
c     now a,b T form
c     ------------------------------------------------------------------      
      do 181 iaj=0,1
      do 181 ibj=0,1
      do 181 itj=0,1
      aj=iaj
      bj=ibj
      cjt1=stat(aj)*stat(bj)/2.d0 
      do 182 ikj=0,2
      do 183 iqj=-ikj,ikj
      mfoff(iaj,ibj,itj,ikj,iqj,ik,jq)=(0.d0,0.d0)
      
      rkj=ikj
      do 183 isj=0,1
      sjt=isj
      cjt2=cjt1*u9(0.5d0,0.5d0,sjt,0.5d0,0.5d0,sjt,aj,bj,rkj)
      cjt3=cjt2*(stat(sjt))**3
      mfoff(iaj,ibj,itj,ikj,iqj,ik,jq)=mfoff(iaj,ibj,itj,ikj,iqj,ik,jq)+
     +       cjt3*mkq(isj,itj,ikj,iqj,ik,jq)
  183 continue
    
c Hooton-Johnson a,b,T amplitudes
      do 184 iqj=-ikj,ikj
         mfHJ(iaj,ibj,itj,ikj,iqj,ik,jq)=(0.d0,0.d0)
      do 184 iqp=-ikj,ikj
         mfHJ(iaj,ibj,itj,ikj,iqj,ik,jq)= 
     &   mfHJ(iaj,ibj,itj,ikj,iqj,ik,jq)
     & + dmat(ikj,iqp,iqj,angle)*mfoff(iaj,ibj,itj,ikj,iqp,ik,jq)
  184 continue
  182 continue
  181 continue
140   continue





! Write amplitudes a,b,T to external files
      write(99,*) 'Writing off-shell amplitudes...'
      do 850 jq=1,nqmax
      xq=xxq(jq)
      do 850 ik=1,kkmax
      xkk=xxk(ik)

! mfHJ(a,b,t,ik,iq,k,q): Hooton-Jonhson amplitudes
      write(66,1023) xq,xkk,mfHJ(0,0,0,0,0,ik,jq),mfHJ(0,0,1,0,0,ik,jq)
      write(67,1023) xq,xkk,mfHJ(0,1,0,1,1,ik,jq),mfHJ(0,1,1,1,1,ik,jq)
      write(68,1023) xq,xkk,mfHJ(1,1,0,0,0,ik,jq),mfHJ(1,1,1,0,0,ik,jq)
      write(69,1023) xq,xkk,mfHJ(1,1,0,2,0,ik,jq),mfHJ(1,1,1,2,0,ik,jq)
      write(70,1023) xq,xkk,mfHJ(1,1,0,2,1,ik,jq),mfHJ(1,1,1,2,1,ik,jq)
      write(71,1023) xq,xkk,mfHJ(1,1,0,2,2,ik,jq),mfHJ(1,1,1,2,2,ik,jq)

      select case (moff)
         case (2)
! mfoff(a,b,t=0,ik,iq,k,q)
         m1a=mfoff(0,0,0,0,0,ik,jq)         
         m2a=mfoff(0,1,0,1,1,ik,jq)
         m3a=mfoff(1,1,0,0,0,ik,jq)
         m4a=mfoff(1,1,0,2,0,ik,jq)
         m5a=mfoff(1,1,0,2,1,ik,jq)
         m6a=mfoff(1,1,0,2,2,ik,jq)
! mfoff(a,b,t=1,ik,iq,k,q)
         m1b=mfoff(0,0,1,0,0,ik,jq)
         m2b=mfoff(0,1,1,1,1,ik,jq)
         m3b=mfoff(1,1,1,0,0,ik,jq)
         m4b=mfoff(1,1,1,2,0,ik,jq)
         m5b=mfoff(1,1,1,2,1,ik,jq)
         m6b=mfoff(1,1,1,2,2,ik,jq)

         case (1)
! Isoscalar:
         m1a=(mfoff(0,0,0,0,0,ik,jq)+3*mfoff(0,0,1,0,0,ik,jq))/4.
         m2a=(mfoff(0,1,0,1,1,ik,jq)+3*mfoff(0,1,1,1,1,ik,jq))/4.
         m3a=(mfoff(1,1,0,0,0,ik,jq)+3*mfoff(1,1,1,0,0,ik,jq))/4.
         m4a=(mfoff(1,1,0,2,0,ik,jq)+3*mfoff(1,1,1,2,0,ik,jq))/4.
         m5a=(mfoff(1,1,0,2,1,ik,jq)+3*mfoff(1,1,1,2,1,ik,jq))/4.
         m6a=(mfoff(1,1,0,2,2,ik,jq)+3*mfoff(1,1,1,2,2,ik,jq))/4.
! Isovector         
         m1b=(mfoff(0,0,1,0,0,ik,jq)-mfoff(0,0,0,0,0,ik,jq))/4.
         m2b=(mfoff(0,1,1,1,1,ik,jq)-mfoff(0,1,0,1,1,ik,jq))/4.
         m3b=(mfoff(1,1,1,0,0,ik,jq)-mfoff(1,1,0,0,0,ik,jq))/4.
         m4b=(mfoff(1,1,1,2,0,ik,jq)-mfoff(1,1,0,2,0,ik,jq))/4.
         m5b=(mfoff(1,1,1,2,1,ik,jq)-mfoff(1,1,0,2,1,ik,jq))/4.
         m6b=(mfoff(1,1,1,2,2,ik,jq)-mfoff(1,1,0,2,2,ik,jq))/4. 

        case (3)
! pn
         m1a=(mfoff(0,0,0,0,0,ik,jq)+mfoff(0,0,1,0,0,ik,jq))/2.
         m2a=(mfoff(0,1,0,1,1,ik,jq)+mfoff(0,1,1,1,1,ik,jq))/2.
         m3a=(mfoff(1,1,0,0,0,ik,jq)+mfoff(1,1,1,0,0,ik,jq))/2.
         m4a=(mfoff(1,1,0,2,0,ik,jq)+mfoff(1,1,1,2,0,ik,jq))/2.
         m5a=(mfoff(1,1,0,2,1,ik,jq)+mfoff(1,1,1,2,1,ik,jq))/2.
         m6a=(mfoff(1,1,0,2,2,ik,jq)+mfoff(1,1,1,2,2,ik,jq))/2.

! pp 
         m1b=mfoff(0,0,1,0,0,ik,jq)
         m2b=mfoff(0,1,1,1,1,ik,jq)
         m3b=mfoff(1,1,1,0,0,ik,jq)
         m4b=mfoff(1,1,1,2,0,ik,jq)
         m5b=mfoff(1,1,1,2,1,ik,jq)
         m6b=mfoff(1,1,1,2,2,ik,jq)
        end select

        if (moff.ne.0) then
           write(60,1023) xq,xkk,m1a,m1b
           write(61,1023) xq,xkk,m2a,m2b
           write(62,1023) xq,xkk,m3a,m3b
           write(63,1023) xq,xkk,m4a,m4b
           write(64,1023) xq,xkk,m5a,m5b
           write(65,1023) xq,xkk,m6a,m6b
        endif
           

c      write(66,1023) xq,xkk,mfoff(1,0,0,1,1,ik,jq),mfoff(1,0,1,1,1,ik,jq)
c      write(67,1023) xq,xkk,mfoff(1,1,0,1,1,ik,jq),mfoff(1,1,1,1,1,ik,jq)
c      write(68,1023) xq,xkk,mfoff(1,1,0,2,2,ik,jq),mfoff(1,1,1,2,2,ik,jq)
c      write(69,1023) xq,xkk,mfoff(0,1,0,1,0,ik,jq),mfoff(0,1,1,1,0,ik,jq)
c      write(70,1023) xq,xkk,mfoff(1,0,0,1,0,ik,jq),mfoff(1,0,1,1,0,ik,jq)
c      write(71,1023) xq,xkk,mfoff(1,1,0,1,0,ik,jq),mfoff(1,1,1,1,0,ik,jq)
 850  continue
 1023 format(2f8.4,4e14.6)
      write(99,*)'- Deallocating xkm,xkmp,ctheta'
      deallocate(xkm,xkmp,ctheta,thpr)
      deallocate(d0,dz,dm,dp,dpp,dmm)
      deallocate(mfoff,mst,mkq,mfHJ)
      close(60);close(61);close(62);close(63);close(64);close(65)
      return
      end












