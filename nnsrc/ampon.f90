*****************************************************************
*     E.F. Redish, K. Stricker-Bauer (1987)
*     Formalism described in: Phys. Rev. C36, 513 (1987)
*****************************************************************
!      program ampon
       subroutine amponsh
        use jat1
        use j
        use rc
        use ampold
        use ampnew
        use nnamp
        use parms
        use constants

      implicit real*8(a-h,o-z)
      integer,save:: ncalls=0
c     program to generate a(q)-e(q) on shell
!      parameter (nth=91)
!      parameter (lmx=20)
!      parameter (lnomx=50)
!      parameter (kmx=50)
!      parameter (nmx=2*kmx*(2*kmx+1)/2)
!      parameter (np=91)
!      complex*16 ttl,tl,zi,aa,bb,cc,dd,ee,ff,t,tlon,cint2d
!      complex*16 d01,dz1,dm1,dp1,dpp1,dmm1
!      complex*16 d0,dz,dm,dp,dpp,dmm
!      complex*16 m00,m11,m1m1,m10,m01,mss,maux1,maux2
!      common/jat1/ vl(nmx),tl(nmx),xx(kmx),wx(kmx),ttl(kmx,kmx)
!      common/j/ xkm(np),xkmp(np),xxk(np),xxq(np),pol(lmx),cs(nth)
       complex*16:: m00(2),m11(2),m1m1(2),m10(2),m01(2),mss(2)
!      dimension lscan(lnomx)
!      common/rc/ctheta(np)
!      common/ampold/d01(lmx,nmx),dz1(lmx,nmx),dp1(lmx,nmx),dm1(lmx,nmx),
!     1              dpp1(lmx,kmx,kmx),dmm1(lmx,kmx,kmx)
!      common/ampnew/d0(lmx),dz(lmx),dp(lmx),dm(lmx),dpp(lmx),dmm(lmx)
!      common/nnamp/aa(2,np),bb(2,np),cc(2,np),dd(2,np),ee(2,np),ff(2,np)
!      common/parms/pi,ttof
 
      open(unit=14,file='tout.dat',status='unknown')
      open(unit=22,file='ampi.a',status='unknown')
      open(unit=23,file='ampi.b',status='unknown')
      open(unit=24,file='ampi.c',status='unknown')
      open(unit=25,file='ampi.d',status='unknown')
      open(unit=26,file='ampi.e',status='unknown')
 
      norder=5
      ch=197.33
      
	if (pm*tm.lt.1e-3) then
        write(*,*)'** WARNING ** ampon:pm=',pm,'tm=',tm
        stop
      endif

      rm=pm*tm/(pm+tm)
     
!      rm=938.9/2.
      ttof=ch*ch*pi/rm
      write(*,*)'ampon: pm,tm,ttof=',pm,tm,ttof
      zi=cmplx(0.0,1.0)
! nth=2 ????????????????
      dth=2.d0
c     q=(k"-k), kk=(k+k")/2*)

      read(14,*)kmax
      read(14,1010)p,redm,x0,w0,(xx(i),wx(i),i=1,kmax)
 1010 format(8e14.7)

      do 110 jq=1,nth
      th=(jq-1)*dth
      cth=cos(th*pi/180.d0)
      sth=sin(th*pi/180.d0)
c     xq is momentum transfer (fm -1)
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

c     read tl's from file 14'
      lno = 0
 100  read(14,1000)is,ll,jj,it,ic
 1000 format(6i5)  
      if(is.lt.0)go to 160
      lno = lno+1
!      lscan(lno) = ll
      if(kmax.gt.kmx) stop 'kmax'
      if(kmax.lt.0)go to 160
      nmax=kmax*(kmax+1)/2

c     spin zero case
      if(is.eq.0)then
      read(14,1010)(d01(jj+1,kn),kn=1,nmax)

c     spin one case (ic=1->uncoupled, ic=2->coupled)
      else
      if (ic.eq.1) then
      if (jj.gt.0)then
      read(14,1010)(dz1(jj+1,kn),kn=1,nmax)
      else
      read(14,1010)(dm1(jj+2,kn),kn=1,nmax)
      endif
      else if (ic.eq.2) then
      do 105 m=1,ic
      do 105 n=1,ic
      if(m.eq.1.and.n.eq.1)then
      read(14,1010)(dm1(jj+2,kn),kn=1,nmax)
      else if(m.eq.1.and.n.eq.2)then
      read(14,1010)((dmm1(jj,kn,knp),kn=1,kmax),knp=1,kmax)
      else if(m.eq.2.and.n.eq.1)then
      read(14,1010)((dpp1(jj+2,kn,knp),kn=1,kmax),knp=1,kmax)
      else if(m.eq.2.and.n.eq.2)then
      read(14,1010)(dp1(jj,kn),kn=1,nmax)
      endif
 105  continue
      endif
      endif
      go to 100

 160  continue
c     calculates maximum number of partial waves read in
      iflag =1
!      do 130 i=1,lnomx
!      aux = lscan(iflag)-lscan(i)
!  130 if(aux.lt.0)iflag=i
!      lmu = lscan(iflag)+1

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
      ncalls=ncalls+1
      call ampcal(thqk,lmu,m11,m10,m1m1,m01,m00,mss,ncalls)

c     calculates the NN scattering amplitudes A-E
c     ir=1:T=0, ir=2:T=1
      do 145 ir = 1,2
      aa(ir,jq) = 0.d0
      bb(ir,jq) = 0.d0
      cc(ir,jq) = 0.d0
      dd(ir,jq) = 0.d0
      ee(ir,jq) = 0.d0
 145  continue

      do 150 ir=1,2
      aa(ir,jq)=(2.*m11(ir)+m00(ir)+mss(ir))/4. + aa(ir,jq)
      bb(ir,jq)=(-2.*m1m1(ir)+m00(ir)-mss(ir))/4. + bb(ir,jq)
      cc(ir,jq)=zi*(m10(ir)-m01(ir))/2./sqrt(2.) + cc(ir,jq)
      maux1 = m11(ir)+m1m1(ir)-mss(ir)
      maux2 = m11(ir)-m1m1(ir)-m00(ir)
      dd(ir,jq)=(maux1-maux2/cos(thqk))/4. + dd(ir,jq)
      ee(ir,jq)=(maux1+maux2/cos(thqk))/4. + ee(ir,jq)
 150  continue
 140  continue

      ttol=0.001d0
      call badpt(ttol,nth,1,ee)
      call badpt(ttol,nth,2,ee)
      call badpt(ttol,nth,1,dd) 
      call badpt(ttol,nth,2,dd)

      do 850 ith=1,nth
      th=(ith-1)*dth
      xkk=xxk(ith)
      xq=xxq(ith)
      write(22,1023) th,xq,aa(1,ith),aa(2,ith)
      write(23,1023) th,xq,bb(1,ith),bb(2,ith)
      write(24,1023) th,xq,cc(1,ith),cc(2,ith)
      write(25,1023) th,xq,dd(1,ith),dd(2,ith)
      write(26,1023) th,xq,ee(1,ith),ee(2,ith)
 850  continue

 1023 format('  ',f5.1,f8.4,4e14.6)
      end
