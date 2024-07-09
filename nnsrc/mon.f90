!*********************************************************
!    mabon : Evaluates NN on-shell amplitudes in tensor 
!            representation
!    J. Tostevin (July 1995)
!
!    Formalism and f90 version: A.M. Moro, R. Crespo (July 2001): 
!                Phys. Rev. C 65, 054001 (2002)
!**********************************************************

!     Mixed representation tensors a,b,T,    JAT July 1995
       subroutine mabon()
        use jat1
        use jj2
        use rc2
        use ampold
        use ampnew
        use nnamp
        use parms
        use toff
        use ampnl
        use trace
        use amps
!        use constants

      implicit none
      integer,save::ncalls=0
      integer :: i,kth,
     $            kkmax,nqmax,jq,ik,nmax,
     $            n,m,kn,ir,kr,iprt,ith,
     $            knp,norder,itj,isj,ikj,iqj,isgpj,isgj,iaj,ibj
       real*8 :: rm,dth,eps,sjt,cjt1,stat,rkj,cjt2,qj,sgj,
     $       fkq,sthkq,cthkq,th,xq,xkk,xk,xkp,cth,aj,bj,u9,cjt3,cleb,
     $       x0,w0,sth,epsrc,alpha,beta,thqk,norm,zi,sgpj,raux
      complex*16 ,dimension(2) :: m00,m11,m1m1,m10,m01,mss,maux1,maux2
      complex*16 , allocatable :: mst(:,:,:,:,:),mkq(:,:,:,:,:)
      complex*16 :: m1a,m2a,m3a,m4a,m5a,m6a,m1b,m2b,m3b,m4b,m5b,m6b

! This is required to pass dummy arrays as arguments
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

c --------------------- end interfaces --------------------------------
       
       allocate(mst(0:1,0:1,-1:1,-1:1,nth))
       allocate(mkq(0:1,0:1,0:2,-2:2,nth))
       if (allocated(mfon)) deallocate(mfon)
       allocate(mfon(0:1,0:1,0:1,0:2,-2:2,nth))
       mst=(0.,0.)
       mkq=(0.,0.)
       mfon=(0.,0.)
       ncalls=ncalls+1

!! COMMENTED BY AMORO (nthin no longer used) (26/11/203)
!!       nthin=nth !compatibility with MST program
!!
!!       write(99,*)' - mon: nthin=',nthin

       if (allocated(d0)) deallocate(d0,dz,dm,dmm,dpp,dp)

       allocate(d0(lmu))
       allocate(dz(lmu))
       allocate(dm(lmu))
       allocate(dmm(lmu))
       allocate(dpp(lmu))
       allocate(dp(lmu))
       d0=0d0;dz=0d0;dm=0d0;dmm=0d0;dpp=0d0;dp=0d0
       
      write(99,*) 'Calculating coeficients in tensor representation'
      write(*,*) '- Calculating on-shell amplitudes M(a,b) lmu=',lmu

      if (mon.ne.0) then
         open(unit=22,file='m0000.on',status='unknown')
         open(unit=23,file='m0111.on',status='unknown')
         open(unit=24,file='m1100.on',status='unknown')
         open(unit=25,file='m1120.on',status='unknown')
         open(unit=26,file='m1121.on',status='unknown')
         open(unit=27,file='m1122.on',status='unknown')
      endif
c      open(unit=28,file='m1121.on',status='unknown')
c      open(unit=29,file='m1122.on',status='unknown')

      norder=5
      zi=cmplx(0.0,1.0)
      dth=180./(nth-1)
c     q=(k"-k), kk=(k+k")/2*)


 1010 format(8e14.7)

      if(allocated(xxq))deallocate(xxq)
      if(allocated(xxk))deallocate(xxk)
      
      allocate(xxq(nth))
      allocate(xxk(nth))
      allocate(xkm(nth))
      allocate(xkmp(nth))
      allocate(ctheta(nth))

      do 110 jq=1,nth
      th=(jq-1)*dth
      cth=cos(th*pi/180.d0)
      sth=sin(th*pi/180.d0)
c  xq is momentum transfer (fm -1)
      xq=p*sqrt(abs((1.-cth)*2.d0))
      xxq(jq)=xq
      xkk=p*sqrt(abs((1.+cth)/2.d0))
      xxk(jq)=xkk
      xk=sqrt(xq*xq/4.+xkk*xkk)
      xkp=xk
      xkm(jq)=xk
      xkmp(jq)=xkp
      ctheta(jq)=1.
      if(xk.le.0.d0.or.xkp.le.0.d0) go to 115
      ctheta(jq)=(xkk*xkk-xq*xq/4.d0)/xk/xkp
 115  continue
 110  continue


c     begin loops to interpolate t on q mesh,
      do 140 jq=1,nth
      xq = xxq(jq)
      xkk = xxk(jq)
      xk=xkm(jq)
      xkp=xkmp(jq)
      if(ctheta(jq).gt.0.999999999999d0) then
      thqk=0.d0
      else if(ctheta(jq).lt.-0.999999999999d0) then
      thqk=pi
      else
      thqk = acos(ctheta(jq))
      endif

c     generates all partial wave amplitudes in the new grid
      call newgd(lmu,kmax,xk,xkp,d01,d0)
      call newgd(lmu,kmax,xk,xkp,dz1,dz)
      call newgd(lmu,kmax,xk,xkp,dm1,dm)
      call newgd(lmu,kmax,xk,xkp,dp1,dp)
      call newgnd(lmu,kmax,xk,xkp,dpp1,dpp)
      call newgnd(lmu,kmax,xk,xkp,dmm1,dmm)

c     calculates M00,M01,M1-1,M10,M11,MSS
      call ampcal(thqk,lmu,m11,m10,m1m1,m01,m00,mss,ncalls)

c     calculates the NN scattering amplitudes A-E
c     ir=1:T=0, ir=2:T=1

      do 111 ir=1,2
c     ------------------------------------------------------------------
c     assign S=1 amplitudes to mst array, from ampcal
c     ------------------------------------------------------------------
      mst(1,ir-1,1, 1,jq)=m11(ir)
      mst(1,ir-1,1, 0,jq)=m10(ir)
      mst(1,ir-1,1,-1,jq)=m1m1(ir)
      mst(1,ir-1,0, 1,jq)=m01(ir)
      mst(1,ir-1,0, 0,jq)=m00(ir)

c     ------------------------------------------------------------------
c     from symmetry relations
c     ------------------------------------------------------------------
      mst(1,ir-1, 0,-1,jq)=-mst(1,ir-1,0, 1,jq)
      mst(1,ir-1,-1, 1,jq)= mst(1,ir-1,1,-1,jq)
      mst(1,ir-1,-1, 0,jq)=-mst(1,ir-1,1, 0,jq)
      mst(1,ir-1,-1,-1,jq)= mst(1,ir-1,1, 1,jq)
c     ------------------------------------------------------------------
c     assign S=0 amplitudes, from ampcal
c     ------------------------------------------------------------------
      mst(0,ir-1, 0, 0,jq)=mss(ir)
c     ------------------------------------------------------------------
  111 continue
    
c     ------------------------------------------------------------------
c     translate to Mkq form (first zero array)
c     ------------------------------------------------------------------
      do 171 itj=0,1
      do 171 isj=0,1
      do 171 ikj=0,2,1
      do 171 iqj=-2,2,1
      mkq(isj,itj,ikj,iqj,jq)=(0.d0,0.d0)
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
c     sum over sigma' and sigma for Mkq'
c     ------------------------------------------------------------------
      do 1115 isgpj=-isj,isj,1
      sgpj=isgpj
      do 1115 isgj=-isj,isj,1
      sgj=isgj
      mkq(isj,itj,ikj,iqj,jq)=mkq(isj,itj,ikj,iqj,jq)+
     +      mst(isj,itj,isgpj,isgj,jq)*cleb(sjt,sgpj,rkj,qj,sjt,sgj)
c     ------------------------------------------------------------------
 1115 continue
      mkq(isj,itj,ikj,iqj,jq)=mkq(isj,itj,ikj,iqj,jq)*cjt2

 1116 continue
  113 continue


c     ------------------------------------------------------------------
c     now a,b T form
c     ------------------------------------------------------------------      
      do 181 iaj=0,1
      do 181 ibj=0,1
      do 181 itj=0,1
      aj=iaj
      bj=ibj
      cjt1=stat(aj)*stat(bj)/2.d0 
      do 182 ikj=0,2,1
      do 182 iqj=-ikj,ikj,1
      mfon(iaj,ibj,itj,ikj,iqj,jq)=(0.d0,0.d0)
      rkj=ikj
      do 183 isj=0,1
      sjt=isj
      cjt2=cjt1*u9(0.5d0,0.5d0,sjt,0.5d0,0.5d0,sjt,aj,bj,rkj)
      cjt3=cjt2*(stat(sjt))**3
      mfon(iaj,ibj,itj,ikj,iqj,jq)=mfon(iaj,ibj,itj,ikj,iqj,jq)+
     +       cjt3*mkq(isj,itj,ikj,iqj,jq)
  183 continue
  182 continue
  181 continue



*     do 150 ir=1,2
*     aa(ir,jq)=(2.*m11(ir)+m00(ir)+mss(ir))/4. + aa(ir,jq)
*     bb(ir,jq)=(-2.*m1m1(ir)+m00(ir)-mss(ir))/4. + bb(ir,jq)
*     cc(ir,jq)=zi*(m10(ir)-m01(ir))/2./sqrt(2.) + cc(ir,jq)
*     maux1 = m11(ir)+m1m1(ir)-mss(ir)
*     maux2 = m11(ir)-m1m1(ir)-m00(ir)
*     dd(ir,jq)=(maux1-maux2/cos(thqk))/4. + dd(ir,jq)
*     ee(ir,jq)=(maux1+maux2/cos(thqk))/4. + ee(ir,jq)
*150  continue
 
 140  continue

*     ttol=0.001d0
*     call badpt(ttol,nth,1,ee)
*     call badpt(ttol,nth,2,ee)
*     call badpt(ttol,nth,1,dd) 
*     call badpt(ttol,nth,2,dd)


!      write(21,*)'# NCALLS=',ncalls
      do 850 ith=1,nth
      th=dble(ith-1)*dth
      xkk=xxk(ith)
      xq=xxq(ith)
      
      select case(mon)
        case (1) 
! Isoscalar (alpha)
         m1a=(mfon(0,0,0,0,0,ith)+3*mfon(0,0,1,0,0,ith))/4.d0
         m2a=(mfon(0,1,0,1,1,ith)+3*mfon(0,1,1,1,1,ith))/4.d0
         m3a=(mfon(1,1,0,0,0,ith)+3*mfon(1,1,1,0,0,ith))/4.d0
         m4a=(mfon(1,1,0,2,0,ith)+3*mfon(1,1,1,2,0,ith))/4.d0
         m5a=(mfon(1,1,0,2,1,ith)+3*mfon(1,1,1,2,1,ith))/4.d0
         m6a=(mfon(1,1,0,2,2,ith)+3*mfon(1,1,1,2,2,ith))/4.d0
! Isovector (beta)
         m1b=(-mfon(0,0,0,0,0,ith)+mfon(0,0,1,0,0,ith))/4.d0
         m2b=(-mfon(0,1,0,1,1,ith)+mfon(0,1,1,1,1,ith))/4.d0
         m3b=(-mfon(1,1,0,0,0,ith)+mfon(1,1,1,0,0,ith))/4.d0
         m4b=(-mfon(1,1,0,2,0,ith)+mfon(1,1,1,2,0,ith))/4.d0
         m5b=(-mfon(1,1,0,2,1,ith)+mfon(1,1,1,2,1,ith))/4.d0
         m6b=(-mfon(1,1,0,2,2,ith)+mfon(1,1,1,2,2,ith))/4.d0  

        case (2)
c       print*,'-Writing T=0,1 Tensor on-shell amplitudes'
        m1a=mfon(0,0,0,0,0,ith)
        m2a=mfon(0,1,0,1,1,ith)
        m3a=mfon(1,1,0,0,0,ith)
        m4a=mfon(1,1,0,2,0,ith)
        m5a=mfon(1,1,0,2,1,ith)
        m6a=mfon(1,1,0,2,2,ith)
        
        m1b=mfon(0,0,1,0,0,ith)
        m2b=mfon(0,1,1,1,1,ith)
        m3b=mfon(1,1,1,0,0,ith)
        m4b=mfon(1,1,1,2,0,ith)
        m5b=mfon(1,1,1,2,1,ith)
        m6b=mfon(1,1,1,2,2,ith)
        
       case (3)
c      print*,'-Writing pn Tensor on-shell amplitudes'
! pn
      m1a=(mfon(0,0,0,0,0,ith)+mfon(0,0,1,0,0,ith))/2.d0
      m2a=(mfon(0,1,0,1,1,ith)+mfon(0,1,1,1,1,ith))/2.d0
      m3a=(mfon(1,1,0,0,0,ith)+mfon(1,1,1,0,0,ith))/2.d0
      m4a=(mfon(1,1,0,2,0,ith)+mfon(1,1,1,2,0,ith))/2.d0
      m5a=(mfon(1,1,0,2,1,ith)+mfon(1,1,1,2,1,ith))/2.d0
      m6a=(mfon(1,1,0,2,2,ith)+mfon(1,1,1,2,2,ith))/2.d0

! pp (T=1) 
       m1b=mfon(0,0,1,0,0,ith)
       m2b=mfon(0,1,1,1,1,ith)
       m3b=mfon(1,1,1,0,0,ith)
       m4b=mfon(1,1,1,2,0,ith)
       m5b=mfon(1,1,1,2,1,ith)
       m6b=mfon(1,1,1,2,2,ith)      
      end select
    
      if (mon.ne.0) then
      write(22,1023) th,xq,m1a,m1b
      write(23,1023) th,xq,m2a,m2b
      write(24,1023) th,xq,m3a,m3b
      write(25,1023) th,xq,m4a,m4b
      write(26,1023) th,xq,m5a,m5b
      write(27,1023) th,xq,m6a,m6b
      end if

!!!!!!!!
!      write(21,1023) xq,xkk,mfon(0,0,0,0,0,ith),mfon(0,0,1,0,0,ith)
!!!!!!!!
      
 850  continue
      

      do 950 ith=1,nth
      th=(ith-1)*dth
      xkk=xxk(ith)
      xq=xxq(ith)
c      write(53,1023) th,xq,mfon(0,1,0,1,0,ith),mfon(0,1,1,1,0,ith)
c      write(54,1023) th,xq,mfon(1,0,0,1,0,ith),mfon(1,0,1,1,0,ith)
c      write(56,1023) th,xq,mfon(1,1,0,1,0,ith),mfon(1,1,1,1,0,ith)
 950  continue
!       write(21,*)'#------------------------------------------------'

!      deallocate(xxq,xxk)
      deallocate(xkm,xkmp,ctheta)
      deallocate(d0,dz,dm,dp,dpp,dmm)
      if (mon.ne.0)then
      close(22);close(23);close(24);close(25);close(26)
      close(28);close(29)
      endif
1023  format('  ',f5.1,f8.4,4e14.6)
      return
      end


