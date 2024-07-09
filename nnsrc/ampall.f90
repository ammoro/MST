*****************************************************************
*     E.F. Redish, K. Stricker-Bauer (1987)
*     Formalism described in: Phys. Rev. C36, 513 (1987)
*****************************************************************
c Routine to generate a(q,k)-e(q,k)

      subroutine ampall ()

*****************************************************************
*   f90 version
*   R.Crespo, A.M.Moro (2001)
*
*   Main changes:
*   - kmx, nmx,lnomx,nmx,np parameters no longer used 
*   - common blocks replaced by modules 
*****************************************************************
      use ampold
      use ampnew
      use nnamp
      use ampaux
      use j
      use rc
      use jat1
      use parms
      use toff
      use trace
      use ampnl
      use ampon
!      use constants

      implicit none
      integer,save::ncalls=0
      integer ::  i,kth,
     $            kkmax,nqmax,jq,ik,nmax,
     $            n,m,kn,ir,kr,iprt,ith,
     $            knp,norder
       real*8 :: rm,dth,eps,
     $       fkq,sthkq,cthkq,th,xq,xkk,xk,xkp,cth,
     $       x0,w0,sth,epsrc,alpha,beta,thqk,norm

      real*8, allocatable :: thpr(:)
      complex*16 :: aaon,bbon,ccon,ddon,eeon,ffon
      complex*16 :: aa2on,bb2on,cc2on,dd2on,ee2on,ff2on
      complex*16 :: aa3on,bb3on,cc3on,dd3on,ee3on,ff3on
      complex*16 :: aa4on,bb4on,cc4on,dd4on,ee4on,ff4on
      complex*16 :: aa5on,bb5on,cc5on,dd5on,ee5on,ff5on
      complex*16 :: zi,maux1,maux2,cint2db
      complex*16, dimension(2):: m00(2),m11(2),m1m1(2),
     $                           m10(2),m01(2),mss(2)
!      common/bitype/itype
      real*8, dimension(2) :: fisp(2,2)
      data fisp/0.25,-0.25,0.75,0.25/

! This is required to pass dummy dynamical arrays as arguments
       interface
         subroutine newgd(l,kmax,xk,xkp,dz1,dz)
           integer:: l,kmax
           real*8 :: xk, xkp
           complex*16, dimension(:,:), intent(in) :: dz1
           complex*16, dimension(:), intent(out) ::  dz
         end subroutine newgd
       end interface

       interface
         subroutine newgnd(l,kmax,xk,xkp,dpp1,dpp)
         integer:: l,kmax
         real*8 :: xk, xkp
         complex*16, dimension(:,:,:), intent(in) :: dpp1
         complex*16, dimension(:), intent(out) ::  dpp
        end subroutine newgnd
       end interface
     
      write(99,*)'Entering ampall with lmu=',lmu
      write(*,*) '- Calculating A-E amplitudes (Wolfestein)...'

      open(unit=14,file='tout.dat',status='unknown')

c      norder=5
      zi=cmplx(0.0,1.0)
      eps=1.e-3
      ncalls=ncalls+1
c
c       ifkq = 0 for q=(k"-k)/sqrt(2), kk=(k+k")/sqrt(2)
c            = 1 for q=(k"-k)        , kk=(k+k")/2*)
c     xkmax= max val. of kk, xqmax= max val. of q
c     dk   = kk step,        dq   = q step
c     theta= angle bet. q and kk


c      read(10,*) ifkq
      if(ifkq.eq.0)fkq=1.
      if(ifkq.eq.1)fkq=2.
!      write(99,*) 'theta,icase,itype,nth=',theta,icase,itype,nth
      write(99,*) 'theta,icase,nth=',theta,icase,nth

!      write(*,*)'ampall: p,pi=',p,pi
    
      call mesh(xkmax,xqmax,dk,dq,kkmax,nqmax)

      write(99,*) '+Allocating memory for aa,bb,...'
      write(99,*)' kkmax,nqmax=',kkmax,nqmax
      if(allocated(aa)) then
!        write(0,*)'ampall: ERROR: aa already allocated'
        deallocate(aa,bb,cc,dd,ee,ff)
!        stop
      endif
      allocate(aa(2,kkmax,nqmax))
      allocate(bb(2,kkmax,nqmax))
      allocate(cc(2,kkmax,nqmax))
      allocate(dd(2,kkmax,nqmax))
      allocate(ee(2,kkmax,nqmax))
      allocate(ff(2,kkmax,nqmax))

      write(99,*) '+Allocating memory for ampl,ampw,ampz...'
      if(allocated(ampl)) then
!        write(0,*)'ampall: ERROR: ampl already allocated'
        deallocate(ampl,ampw,ampz)
!        stop
      endif

      allocate(ampl(kkmax,nqmax))
      allocate(ampw(kkmax,nqmax))
      allocate(ampz(kkmax,nqmax))

      write(99,*) '+Allocating memory for xkm,xkmp,ctheta...'
      if(allocated(xkm)) then
!        write(0,*)'ampall: ERROR: ampl already allocated'
        deallocate(xkm,xkmp,ctheta,thpr)
!        stop
      endif
      allocate(xkm(kkmax,nqmax))
      allocate(xkmp(kkmax,nqmax))
      allocate(ctheta(kkmax,nqmax))
      allocate(thpr(nth))
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

      write(99,*) 'ampall: allocating memory for pol, cs'
!      allocate(pol(lmu))
!      allocate(cs(nth))
      
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
       rewind(14)
       read(14,*)kmax
       if(kmax.lt.0)go to 160

       nmax=kmax*(kmax+1)/2
       
       write(99,'("+Allocating",1i3," elements for d0")') lmu
       if (allocated(d0)) deallocate(d0,dz,dm,dmm,dpp,dp)
       
       allocate(d0(lmu))
       allocate(dz(lmu))
       allocate(dm(lmu))
       allocate(dmm(lmu))
       allocate(dpp(lmu))
       allocate(dp(lmu))
       

c          begin loops to interpolate t and v on kk, q mesh,
c          and add together for off-shell a and c
c
 160      continue

      do 140 jq=1,nqmax
      xq = xxq(jq)
      do 140 ik=1,kkmax
      xkk = xxk(ik)
      xk=xkm(ik,jq)
      xkp=xkmp(ik,jq)
      thqk = acos(ctheta(ik,jq))
      if (xk.ne.xkp) then
      norm = (xk-xkp)**2 + 4.*xk*xkp*cos(thqk/2.)**2
      norm = 1/norm
      alpha = norm*(xk*sin(thqk))**2
      beta = norm*(xk+xkp*cos(thqk))**2
      
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
 
c **      calculates the NN scattering amplitudes A-E
c **      ir=1:A0=(3*A(T=1)+A(T=0))/4, ir=2:Atau=(A(T=1)-A(T=0))/4
c **      kr=1:T=0, kr=2:T=1
      do 145 ir = 1,2
      aa(ir,ik,jq) = 0.d0
      bb(ir,ik,jq) = 0.d0
      cc(ir,ik,jq) = 0.d0
      dd(ir,ik,jq) = 0.d0
      ee(ir,ik,jq) = 0.d0
      ff(ir,ik,jq) = 0.d0
 145  continue

      do 150 ir=1,2
      do 155 kr=1,2
      aa(ir,ik,jq)=(2.d0*m11(kr)+m00(kr)+mss(kr))/4.d0*fisp(ir,kr)
     +              + aa(ir,ik,jq)
      bb(ir,ik,jq)=(-2.*m1m1(kr)+m00(kr)-mss(kr))/4.d0*fisp(ir,kr)
     +              + bb(ir,ik,jq)
      cc(ir,ik,jq)=zi*(m10(kr)-m01(kr))/2./dsqrt(2.d0)*fisp(ir,kr)
     +              + cc(ir,ik,jq)
      if ((abs(xk-xkp).gt.0.0001)) then
      dd(ir,ik,jq)=(2*(beta+2*alpha-1.)*m11(kr) + 2*beta*m1m1(kr)
     *     + (2*beta-1.)*m00(kr) - mss(kr)
     *     + 2*dsqrt(2*alpha*beta)*(m01(kr)+m10(kr)))/4.d0*fisp(ir,kr)
     *     +  dd(ir,ik,jq)
      ee(ir,ik,jq)=(2*(alpha+2*beta-1)*m11(kr) + 2*alpha*m1m1(kr)
     *     + (2*alpha-1.)*m00(kr) - mss(kr)
     *    + 2.*dsqrt(2.*alpha*beta)*(m01(kr)+m10(kr)))/4.*fisp(ir,kr)
     *       + ee(ir,ik,jq)
      ff(ir,ik,jq) = (dsqrt(alpha*beta)*(m11(kr)-m00(kr)-m1m1(kr))
     *      + (alpha-beta)/dsqrt(2.d0)*(m10(kr)+m01(kr)))*fisp(ir,kr)
     *      + ff(ir,ik,jq)
       else
      maux1 = m11(kr)+m1m1(kr)-mss(kr)
      maux2 = m11(kr)-m1m1(kr)-m00(kr)
      dd(ir,ik,jq)=(maux1-maux2/cos(thqk))/4.d0*fisp(ir,kr)
     +              + dd(ir,ik,jq)
      ee(ir,ik,jq)=(maux1+maux2/cos(thqk))/4.d0*fisp(ir,kr)
     +               + ee(ir,ik,jq)
      ff(ir,ik,jq) =0.
       endif
 155  continue
c     cc(ir,ik,jq) = sthkq * xkk * xq / xk /xkp * cc(ir,ik,jq)
 150  continue

c      if (abs(ctheta(ik,jq)).lt.0.01)then
c      write(22,*)ik,jq,ctheta(ik,jq),ee(1,ik,jq)
c      endif

      if (sin(thqk).ne.0)then
         do 700 kr=1,2
            maux1 = m11(kr)-m1m1(kr)-m00(kr)
            maux2 = dsqrt(2.d0)*cos(thqk)/sin(thqk)*(m10(kr)+m01(kr))  
 700  continue
      endif
      ampl(ik,jq) = ( 1.d0*(m11(2)+m1m1(2)-mss(2))
     *                - (m11(1)+m1m1(1)-mss(1)) )/4.d0
      ampw(ik,jq) = ( 1.d0*(m11(2)-m1m1(2)-m00(2))
     *               - (m11(1)-m1m1(1)-m00(1)) )/4.d0
      
c      write(22,1010) real(m11(1)),dimag(m11(1)),real(m10(1)),
c     *            dimag(m10(1)),real(m1m1(1)),dimag(m1m1(1))
c      write(22,1010) real(m01(1)),dimag(m01(1)),real(m00(1)),
c     *            dimag(m00(1)),real(mss(1)),dimag(mss(1))
 140  continue

c
c ****       calculates bad points such cos(theta)=0
       epsrc = 0.01
       call badpt(epsrc,nqmax,kkmax,1,ee)
       call badpt(epsrc,nqmax,kkmax,2,ee)

       call badpt(epsrc,nqmax,kkmax,1,dd)
       call badpt(epsrc,nqmax,kkmax,2,dd)

c **** calculates on-shell quantities
      if(kth.ne.1)go to 800

      if ((aeon.ne.0).or.(xsec.ne.0)) then      
      if (allocated(aon)) deallocate(aon,con)
      allocate(aon(nth))
      allocate(con(nth))
  
      do 850 ith=1,nth
      th=(ith-1)*dth
      cth=cos(th*pi/180.)
      sth=sin(th*pi/180.)
      if(cth.eq.0) go to 850
      xkk=p*dsqrt((1.+cth)/fkq)
      xq=p*dsqrt((1.-cth)*fkq)
      
      

  
!  ir=1: Aalpha, Balpha, Calpha,....
      call onamp(aa,xkk,xq,kkmax,nqmax,aaon,1)
      call onamp(cc,xkk,xq,kkmax,nqmax,ccon,1)
      call onamp(bb,xkk,xq,kkmax,nqmax,bbon,1)
      call onamp(dd,xkk,xq,kkmax,nqmax,ddon,1)
      call onamp(ee,xkk,xq,kkmax,nqmax,eeon,1)
      call onamp(ff,xkk,xq,kkmax,nqmax,ffon,1)

c      write(23,1013)th,real(aaon),dimag(aaon),real(ccon),
c     *               dimag(ccon),real(bbon),dimag(bbon)

c      write(22,1012)th,real(ddon),dimag(ddon),real(eeon),
c     *              dimag(eeon)


!  ir=2: Atau, Btau, Ctau,...
      call onamp(aa,xkk,xq,kkmax,nqmax,aa2on,2)
      call onamp(cc,xkk,xq,kkmax,nqmax,cc2on,2)
      call onamp(bb,xkk,xq,kkmax,nqmax,bb2on,2)
      call onamp(dd,xkk,xq,kkmax,nqmax,dd2on,2)
      call onamp(ee,xkk,xq,kkmax,nqmax,ee2on,2)
      call onamp(ff,xkk,xq,kkmax,nqmax,ff2on,2)

c      write(33,1013)th,real(aa2on),dimag(aa2on),real(cc2on),
c     *               dimag(cc2on),real(bb2on),dimag(bb2on)
c      write(32,1012)th,real(dd2on),dimag(dd2on),real(ee2on),
c     *              dimag(ee2on)

!  aa3on = Apn = A0-Atau = [A(T=1)+A(T=0)]/2 
      aa3on = aaon - aa2on
      cc3on = ccon - cc2on
      bb3on = bbon - bb2on
      dd3on = ddon - dd2on
      ee3on = eeon - ee2on
      ff3on = ffon - ff2on

      
!  aa4on = App = Ann = A0+Atau = A(T=1) 
      aa4on = aaon + aa2on
      cc4on = ccon + cc2on
      bb4on = bbon + bb2on
      dd4on = ddon + dd2on
      ee4on = eeon + ee2on
      ff4on = ffon + ff2on

!  aa5on = A0 - Atau = A(T=0) 
      aa5on = aaon - 3*aa2on
      cc5on = ccon - 3*cc2on
      bb5on = bbon - 3*bb2on
      dd5on = ddon - 3*dd2on
      ee5on = eeon - 3*ee2on
      ff5on = ffon - 3*ff2on

      select case (aeon)
      case (1)
!  alpha
         write(220,777) th,
     *               real(aaon), dimag(aaon),
     *               real(bbon), dimag(bbon),
     *               real(ccon), dimag(ccon),
     *               real(ddon), dimag(ddon),
     *               real(eeon), dimag(eeon),
     *               real(ffon), dimag(ffon)
       aon(ith)=aaon
       con(ith)=ccon
!  tau
        write(230,777) th,
     *               real(aa2on),dimag(aa2on),
     *               real(bb2on), dimag(bb2on),
     *               real(cc2on), dimag(cc2on),
     *               real(dd2on), dimag(dd2on),
     *               real(ee2on), dimag(ee2on),
     *               real(ff2on), dimag(ff2on)
      aon(ith)=aa2on
      con(ith)=cc2on

      case(3)
!  pn 
        write(320,777) th,real(aa3on),dimag(aa3on),
     *               real(bb3on), dimag(bb3on),
     *               real(cc3on), dimag(cc3on),
     *               real(dd3on), dimag(dd3on),
     *               real(ee3on), dimag(ee3on),
     *               real(ff3on), dimag(ff3on) 
      aon(ith)=aa3on
      con(ith)=cc3on

      case(2)
!  T=1 (=pp=nn)
        write(330,777) th,real(aa4on),dimag(aa4on),
     *               real(bb4on), dimag(bb4on),
     *               real(cc4on), dimag(cc4on),
     *               real(dd4on), dimag(dd4on),
     *               real(ee4on), dimag(ee4on),
     *               real(ff4on), dimag(ff4on) 
      aon(ith)=aa4on
      con(ith)=cc4on
!  T=0     
        write(340,777) th,
     *               real(aa5on), dimag(aa5on),
     *               real(bb5on), dimag(bb5on),
     *               real(cc5on), dimag(cc5on),
     *               real(dd5on), dimag(dd5on),
     *               real(ee5on), dimag(ee5on),
     *               real(ff5on), dimag(ff5on)  
      aon(ith)=aa5on
      con(ith)=cc5on
      end select

 850  continue
      endif
 800  continue

 777  format(14e11.3)

      write(99,*)'-Deallocating xkm,xkmp,ctheta,thpr'
      deallocate(xkm,xkmp,ctheta,thpr)
!      deallocate(d01,dz1,dm1,dmm1,dpp1,dp1,d0,dz,dm,dmm,dpp,dp)
      write(99,*) 'Writing amplitudes'
      call wrt2p(kkmax,nqmax,dk,dq,p,theta,111,113)
c     call wrtall(kkmax,nqmax,dk,dq,p,theta,np,112)
!      if(kth.ne.1)stop
      if(kth.ne.1) return
c
      go to 99
c2000 print 2001,kmax
 2000  continue
 2001 format('mesh pts. on which t is stored=',i3,' too big')
c     print 2002,kmx
 2002 format('must be less than',i3)
c
 1000 format(6i5)
 1010 format(8e14.7)
 1011 format(' Re A',8x,' Re A int',4x,' Im A',8x,' Im A int',4x,
     & ' Re C',8x,' Re C int',4x,' Im C',8x,' Im C int'/)
 1012 format(8e12.4)
 1013 format(7e11.3)
 1050 format('On-shell  A and C amplitudes from',2x,a5,' potl.')
 1051 format('Off-shell  A and C amplitudes from',2x,a5,' potl.')
 1052 format('E. F. Redish and K.  Stricker-Bauer')
 1053 format('0-180 degrees in',f5.1,' degree intervals')
 1060 format('proton-proton case, pcm=',f10.4,' fm-1')
 1061 format('proton-neutron case, pcm=',f10.4,' fm-1')
 1062 format('average pp+pn case, pcm=',f10.4,' fm-1')
 1063 format('Theta        Re A       Im A       Re C       Im C')
 4000 format(' too many partial waves')
c4010 print 4000
 4010 continue
 99   end



c********************************************************
      subroutine wrt2p(nx,ny,dx,dy,p,th,n,n2)
        use nnamp
        use ampaux
        use block
        use trace
        use ampnl
        use amps
        implicit none
c       implicit double precision (a-h,o-z),integer*4(i-n)
       integer ::nx,ny,n,n2,i,j
       real*8 :: dx,dy,p,th
       complex*16, target :: a,b,c,d,e,f

       if (allocated(app)) deallocate(app,cpp,apn,cpn)
       write(99,'("Allocating",2i4,"elmts for app,cpp...")')nx,ny
       allocate(app(nx,ny))
       allocate(cpp(nx,ny))
       allocate(apn(nx,ny))
       allocate(cpn(nx,ny))
       nkmx1=nx
       nqmx1=ny
       qq=>Bq
       q=>Sq

!      write(n,1000)nx,ny,dx,dy,p,th
!      write(n,1001)nx,ny

      do 1 i=1,ny
         do 2 j=1,nx,1
            ampl(j,i) = aa(1,j,i) 
            ampw(j,i) = cc(1,j,i)

! Punch on unit 22 Alpha, Balpha,... components 
          if (aeoff.eq.1) then
               write(22,777) Sq(i),Bq(j),
     $    real(aa(1,j,i)),dimag(aa(1,j,i)),
     $    real(bb(1,j,i)),dimag(bb(1,j,i)),
     $    real(cc(1,j,i)),dimag(cc(1,j,i)),
     $    real(dd(1,j,i)),dimag(dd(1,j,i)),
     $    real(ee(1,j,i)),dimag(ee(1,j,i)),
     $    real(ff(1,j,i)),dimag(ff(1,j,i))

!       write(50,777) Sq(i),Bq(j),
!     $    real(ee(1,j,i)),dimag(ee(1,j,i))


! Punch on unit 23 Abeta, Bbeta,... components
          write(23,777) Sq(i),Bq(j),
     $    real(aa(2,j,i)),dimag(aa(2,j,i)),
     $    real(bb(2,j,i)),dimag(bb(2,j,i)),
     $    real(cc(2,j,i)),dimag(cc(2,j,i)),
     $    real(dd(2,j,i)),dimag(dd(2,j,i)),
     $    real(ee(2,j,i)),dimag(ee(2,j,i)),
     $    real(ff(2,j,i)),dimag(ff(2,j,i))
       endif

! Punch on unit 32 Apn,Bpn,Cpn,Dpn,Epn       
      
          a=aa(1,j,i) - aa(2,j,i)
          b=bb(1,j,i) - bb(2,j,i)
          c=cc(1,j,i) - cc(2,j,i)
          d=dd(1,j,i) - dd(2,j,i)
          e=ee(1,j,i) - ee(2,j,i)
          f=ff(1,j,i) - ff(2,j,i)

          apn(j,i)=a
          cpn(j,i)=c

     
        if (aeoff.eq.3) then
          write(32,777) Sq(i),Bq(j),
     $    real(a),dimag(a),
     $    real(b),dimag(b),
     $    real(c),dimag(c),
     $    real(d),dimag(d),
     $    real(e),dimag(e),
     $    real(f),dimag(f)
       endif

! Punch on unit 33 App,Bpp,Cpp,Dpp,Epp 
      
          a=aa(1,j,i) + aa(2,j,i)
          b=bb(1,j,i) + bb(2,j,i)
          c=cc(1,j,i) + cc(2,j,i)
          d=dd(1,j,i) + dd(2,j,i)
          e=ee(1,j,i) + ee(2,j,i)
          f=ff(1,j,i) + ff(2,j,i)
          app(j,i)=a
          cpp(j,i)=c
       
        if (aeoff.eq.2) then
          write(33,777) Sq(i),Bq(j),
     $    real(a),dimag(a),
     $    real(b),dimag(b),
     $    real(c),dimag(c),
     $    real(d),dimag(d),
     $    real(e),dimag(e),
     $    real(f),dimag(f)
        endif

! Punch on unit 34 T=0 components 
          a=aa(1,j,i) - 3*aa(2,j,i)
          b=bb(1,j,i) - 3*bb(2,j,i)
          c=cc(1,j,i) - 3*cc(2,j,i)
          d=dd(1,j,i) - 3*dd(2,j,i)
          e=ee(1,j,i) - 3*ee(2,j,i)
          f=ff(1,j,i) - 3*ff(2,j,i)

        if (aeoff.eq.2) then
          write(34,777) Sq(i),Bq(j),
     $    real(a),dimag(a),
     $    real(b),dimag(b),
     $    real(c),dimag(c),
     $    real(d),dimag(d),
     $    real(e),dimag(e),
     $    real(f),dimag(f)
       endif

 2     continue
 1    continue

c ***  redish.data & ampnn.data
!!$      write(n)nx,ny
!!$      write(n)Sq,Bq,ampl,ampw
!!$      if (itype.eq.2)then
!!$      do 10 i=1,ny
!!$       do 20 j=nx,1,-1
!!$       ampl(j,i) = aa(1,j,i) - aa(2,j,i)
!!$       ampw(j,i) = cc(1,j,i) - cc(2,j,i)
!!$c       write(128,777) Sq(i),Bq(j),real(ampl(j,i)),dimag(ampl(j,i)),
!!$c     $ real(ampw(j,i)),dimag(ampw(j,i))
!!$ 20    continue
!!$ 10   continue
!!$       write(n2)nx,ny
!!$       write(n2)Sq,Bq,ampl,ampw
!!$      end if
!!$      close(n); close(n2)
!!$      deallocate(aa,bb,cc,dd,ee,ff)
!!$      deallocate(ampl,ampw,ampz)
      
      return
 1000 format(2i5,4e14.7)
 1001 format(2i5)
 777  format(14e11.3)
 778  format(2e11.3)
      end




c *** Calculattes x-section for NN scattering
c *** using previously calculated on-shell amplitudes
       subroutine nnxsec
         use ampon
         use ampnl
         implicit none
         integer:: ith
         real*8::xsect,dth,th
!         if (.not.allocated(aon)) then
!            write(*,*) 'Cannot calculate xsec.Set aeon>0'
!            stop
!         endif
         open(20,file='nnxsec.out',status='unknown')
         dth=180./(nth-1)
         do ith=1,nth
            th=(ith-1)*dth
            xsect=abs(aon(ith))**2 + abs(con(ith))**2
            xsect=xsect*10 ! convert fm^2 -> mb
            write(20,*) th,xsect
         enddo
         close(20)
         return
       end subroutine nnxsec




