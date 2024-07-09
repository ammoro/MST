c ***************************************************************
c    program LPOTPS
c    R. Crespo, R. C. Johnson, and J. A. Tostevin
c    Phys. Rev. C 41, 2257 (1990)

c ** f90 version by R. Crespo and A.M. Moro  (spring 2002)
c ***************************************************************

c **********************************************************
      subroutine vopt2nl(v2l,k,kp,ik,ikp,n1,lmas2)
c **********************************************************
c *   this subroutine calculates the second order coupling *
c *   potentials for oxygen                                *
c **********************************************************
      implicit real*8(a-h,k,o-z)
      real*8 mn,mp
      real*8 mnuc,mr
      integer option6, kcall
      dimension kgama(20),den(50)
      dimension ri0(0:2,0:2,0:30,0:30)
      dimension v2l(0:99,2,6)
      dimension kk(50),wt(48)
      dimension vaux1(2),vaux2(2)
c     dimension fa(10),zexp(0:110),onexp(0:110)
      complex*16 v1c,v2c,z,vaux
      complex*16 vcon(20),vc
      complex*16 vaux1,vaux2
      complex*16 z1,z2,z3,z4,zf,z14,z23
      complex*16 zexp
      complex*16 v2,vrc1(50),vrc2(50)
      complex*16 y1,y2,dy1,dy2
      complex*16  tmatt
      common/params/ hbarc, pi, mp, mn, nz, na
      common/block6/option6
      common/bgrid/kk,wt
      common/bwf/grho(0:4,201)
      common/bprop/deng(50,10)
      common/vcng/vc(2,50,20)
      common/brit1/rit1(0:30,0:30,0:2,0:2,50)
      common/brit2/rit2(0:30,0:30,0:2,0:2,50)
      common/qnumb/lamx,lbmx,l1mx,lfmx
      common/bkcall/kcall
      common/clebma/faclog(500)
c     common/bv2/v2l(0:99,2,6)
c     common/btaux/t1(100000),t2(100000)
      common/btaux/t1(100000)
      common/btmatt/tmatt(0:7,10,50,50)
      common/ffch/ich
      common/bsflip/isflip
      common/blcount/lcount(4,0:99)
      common/bfast/fa(10),zexp(0:110),onexp(0:110)
c
      data g1,h1,k1,g2,h2,k2/6*0.d0/

      rewind 15
      z = cmplx(0.,1.)
      mr = (mn*na)*mp/(mn*na+mp)
      mnuc = mn*na
      h3 = hbarc**3
      k0=kk(n1)
      if (isflip.eq.0)then
      jamx = 1
      else if (isflip.eq.1)then
      jamx = 10
      endif
c *** generates the radial integrals
      call rint2i(k,kp,ik,ikp,n1)

c *** initialize partial wave 2nd order optical potential
      do 5 ld=0,lmas2
      v2l(ld,1,1)=0.
      v2l(ld,2,1)=0.
 5    continue
c
      do 6 l1=0,lmas2+lfmx+2
      zexp(l1) = z**l1
      onexp(l1) = (-1)**l1
 6    continue
c
      do 7 jamp=1,10
      if(jamp.le.5)then
      fa(jamp)=1.
      if(jamp.eq.3)fa(jamp)=fa(jamp)*2.
      else
      fa(jamp)=3.
      if(jamp.eq.8)fa(jamp)=fa(jamp)*2.
      endif
 7    continue
      ishell=0
c *** calculates on-shell shift
      do 8 la=0,lamx
      do 8 lb=0,lbmx
      ishell = ishell + 1
      kgama(ishell)=sqrt(k0*k0+eshell(la)-eshell(lb))
 8    continue
c
c *** loop over partial wave
c     lmas2n = l1mx+lfmx
c      lmas2n = lmas2
      lmas2n =l1mx
      do 20 ld=0,lmas2n
      if(kcall.eq.1)then
c *** calculates and stores all cleb and racah
      im = 0
      ishell = 0
      do 10 la=0,lamx
      do 10 lb=0,lbmx
      ishell = ishell + 1
      lcount(ishell,ld) = im
      do 10 lc=abs(la-lb),la+lb,2
         f1=clebz(real(la),real(lb),real(lc))
     #      * hat(la)*hat(lb)/hat(ld)
      do 10 jc=abs(ld-lc),ld+lc
      do 15 lf=0,lfmx
      do 15 l1=abs(ld-lf),ld+lf,2
      do 15 l2=abs(l1-lc),l1+lc,2
         f2=clebz(real(lc),real(l2),real(l1))
     #      *  hat(l1)*hat(l2)**2
         im = im+1
         g1 = rac(real(lc),real(l2),real(ld),real(lf),real(l1),
     #         real(jc))
         h1 = 0
         k1 = 0
         if(g1.ne.0.)then
         h1 = clebz(real(lf),real(l1),real(ld))
         if(h1.ne.0.)then
         k1 = clebz(real(l2),real(lf),real(jc))
         end if
         end if
         t1(im) = ( hat(lf)**2
     #             *f1*f2*g1*h1*k1*(-1)**(l1+l2))

c        g2 = rac(real(lc),real(l2),real(jc),real(lf),real(l1),
c    #         real(ld))
c        h2 = 0
c        k2 = 0
c        if(g2.ne.0.)then
c        h2 = clebz(real(lf),real(l2),real(ld))
c        if(h2.ne.0.)then
c        k2 = clebz(real(l1),real(lf),real(jc))
c        endif
c        end if
c        t2(im) = ( hat(lf)**2
c    #              *f1*f2*g2*h2*k2)
 15      continue
 10      continue
      write (15) im
      write (15) (t1(i),i=1,im)
      else
      read (15) im
      read (15) (t1(i),i=1,im)
      endif

      vaux = 0.d0
      vaux1(1) = 0.d0
      vaux2(1) = 0.d0
c *** loop in the intermediate momentum grid
      do 30 ip=1,n1
      p = kk(ip)
      icoup=0
      ishell=0
c *** loop in all combinations of quantum numbers
      im=0
      do 40 la=0,lamx
      do 40 lb=0,lbmx
      ishell=ishell+1
      vc(1,ip,ishell)=0.
      do 40 jamp =1,jamx
      im = lcount(ishell,ld)
      do 40 lc=abs(la-lb),la+lb,2
      do 40 jc=abs(ld-lc),ld+lc
      icoup=icoup+1
      v1c = 0.
      v2c = 0.
        do 50 lf=0,lfmx
        z1 = tmatt(lf,jamp,ik,ip)
        z2 = tmatt(lf,jamp,ip,ikp)
        do 50 l1=abs(ld-lf),ld+lf,2
        z3 = zexp(l1)
        z4 = z3*onexp(l1)
        z14 = z1*z4
        z23 = z2*z3
        do 50 l2=abs(l1-lc),l1+lc,2
        im = im+1
        if (t1(im).ne.0.)then
        zf = t1(im) * zexp(l2)
        v1c = v1c +
     #                  zf*rit1(l2,l1,lb,la,ip)
     *                  *z14
c       v2c(ip,icoup) = v2c(ip,icoup) +
        v2c = v2c +
     #                   zf*rit2(l1,l2,lb,la,ip)
     *                   *onexp(l2)*z23
        end if
 50     continue
         vc(1,ip,ishell)=vc(1,ip,ishell)+v1c*v2c
     *                    * fa(jamp)
c        kgama(icoup) = sqrt(k0**2+eshell(la)-eshell(lb))
         vaux1(1) = vaux1(1) + (-1.)*deng(ip,ishell)
     *              *v1c*v2c
 40     continue
 30    continue

        vaux1(1)=0.
        vaux2(1) = 0.
       i=0
       do 80 la=0,lamx
       do 80 lb=0,lbmx
       i = i+1
       if (option6.eq.2.and.abs(la-lb).ne.0)then
c      write(6,*)'option6=',option6
       call vcinter(n1,i)
       endif
        do 90 ip=1,n1
        vaux1(1)=vaux1(1)+(-1.)*deng(ip,i)*vc(1,ip,i)
  90    continue
       vcon(i) = vc(1,n1,i)
       vaux2(1)=vaux2(1) - z*2.*mr* kgama(i)*vcon(i)
 80    continue
c      write(6,*)ld,vaux1(1),vaux2(1)
       if (ld.eq.0)then
c      write(6,*)ld,k,kp,vaux1(1),vaux2(1)
       end if
       vaux = vaux1(1)+vaux2(1)
       v2l(ld,1,1) = - (na-1)/(na*1.)*real(vaux)*4./h3/h3
       v2l(ld,2,1) = -(na-1.)/(na*1.)*dimag(vaux)*4./h3/h3
       v2l(ld,1,2)=v2l(ld,1,1)
       v2l(ld,2,2)=v2l(ld,2,1)
 20    continue
c
      return
      end


c *********************************************************
      subroutine denomg(kk,wt,n1)
c **********************************************************
c *   this subroutine calculates the coupling propagators  *
c *   for oxygen                                           *
c **********************************************************
      implicit real*8(a-h,k,o-z)
      real*8 mn,mp,mr
      integer option6
      dimension kk(50),wt(48),kks(50),wts(50)
      common/params/ hbarc, pi, mp, mn, nz, na
      common/block6/option6
      common/qnumb/lamx,lbmx,l1mx,lfmx
      common/bprop/deng(50,10)

      mr = (mn*na)*mp/(mn*na+mp)
      k02 = kk(n1)*kk(n1)

      icode=022
      b=21000.
      ngp=n1-1
      icoup=0
c *** loop in all combinations of quantum numbers
      do 20 la=0,lamx
      do 20 lb=0,lbmx
      icoup=icoup+1

c***  default values
      k0s2 = k02
      do 25 ip=1,ngp
      kks(ip)=kk(ip)
      wts(ip)=wt(ip)
 25   continue
c***
        if(option6.eq.2.and.abs(la-lb).ne.0)then
        k0s2=k02+eshell(la)-eshell(lb)
        k0s=sqrt(k0s2)
        call gauss2(ngp,icode,k0s,b,kks,wts)
        endif
c *** loop in the k' grid
       sum = 0.d0
       do 30 ikp=1,n1-1
       kp= kks(ikp)
c ***  calculates the coupling propagators
       k2 = kp*kp
       deng(ikp,icoup)=+(4*mr/pi)*k2*wts(ikp)
     #                  /(k2-k0s2)
       sum = sum + wts(ikp)/(k2-k0s2)
 30    continue
c ***  end of loop in k'
       deng(n1,icoup)=-(4*mr/pi)*k0s2*sum
 20   continue
      return
      end


c **************************************************************
      function eshell(la)
c***************************************************************
c *** this function calculates the energy of the single particle
c *** levels for oxigen
      use blocks
      use params
      implicit real*8 (a-h,o-z)
!      real*8 mn,mp,mr
!      integer option6
!     common/block6/option6
!      common/params/ hbarc, pi, mp, mn, nz, na

      mr = (mn*na)*mp/(mn*na+mp)
      eshell = 0.0
      if (option6.eq.2.and.la.eq.1)then
      eshell = 20.0*2.*mr
c     eshell = 0.
      end if
      if(la.gt.1)stop
      return
      end






c ****************************************************************
       subroutine matwf
c ******************************************************************
      use blocks
      use qnumb
      use radgr
      use bgauss3
      use bwf
      use bwfnum
      use ffch
      implicit real*8(a-h,o-z)
!      integer option5
      complex*16 wf(220),wxxi
      dimension rr(220)
!      common/block5/option5
!      common/qnumb/lamx,lbmx,l1mx,lfmx
!      common/radgr/radxis(201),rad(201)
!      common/bgauss3/rmaxr,quin,mquadi,mquado,irmx
!      common/bwf/grho(0:4,201)
!      common/bwfnum/wfs(220),wfp(220),densr(220),drx,nramax
!      common/ffch/ich
c
      if (option5.ne.4)then
      do 20 la=0,lamx
      do 20 ik=1,irmx
      grho(la,ik) = rho(la,radxis(ik))
      if(ich.eq.1) grho(la,ik) = rhochk(la,radxis(ik))
 20   continue
      else
      do 30 i=1,nramax
      rr(i) = (i-1)*drx
  30  continue
      do 40 il=0,1
      do 45 i=1,nramax
      if (il.eq.0) wf(i)=wfs(i)
      if (il.eq.1) wf(i)=wfp(i)
  45  continue
      do 40 ik=1,irmx
      grho(il,ik)=wxxi(radxis(ik),rr,wf,220,nramax)
  40  continue
      endif
      return
      end

c ****************************************************************
       subroutine mbess(kk,n1)
c ****************************************************************
       use params
       use sizes
       use radgr
       use bgauss3
       use bsdess
       use qnumb

       implicit real*8 (a-h,o-z)
       integer :: kout=16
       real*8 k,kp,kk(n1)
!       real*8 mp,mn
       complex*16,allocatable,dimension(:):: xj,xjp,xh1,xh1p
       complex*16 ::x
       integer code,xk
!       dimension kk(50)
!       common/params/ hbarc, pi, mp, mn, nz, na
!       common/sizes/  achp, acmp, wsp, achn, acmn, wsn
!       common/radgr/  radxis(201), radwt(201)
!       common/bgauss3/rmaxr,quin,mquadi,mquado,irmx
!       common/bsdess/sbess(49,201,0:30)
!       common/qnumb/lamx,lbmx,l1mx,lfmx

       lmx=l1mx+lfmx+2
       h2 = hbarc*hbarc
       lmx2=lmx+2

       write(*,*) 'Allocating',lmx2,'elements for xj, xjp...'
       allocate(xj(0:lmx2))
       allocate(xjp(0:lmx2))
       allocate(xh1(0:lmx2))
       allocate(xh1p(0:lmx2))
       allocate(sbess(n1,201,0:lmx2))
!       allocate(radxis(201))
!       allocate(radwt(201))
       
       
       nr = irmx
c      a = 0.001
c      b = 8.001
c      code = 11
c      call gauss2( nr, code, a, b, radxis, radwt )

!       write(*,*)'calling gauss3 from mbess',
!     &   rmaxr,quin,mquadi,mquado
!       write(*,*) 'sizes=',size(radxis),size(radwt)
       call gauss3(rmaxr,quin,mquadi,mquado,radxis,radwt)
!       call gauss3(rmaxr,quin,mquadi,mquado)
       write(*,*)'exiting gauss3 from mbess'
!       write(*,*)'radxis=',radxis
        write(*,*) 'n1,nrm,hbarc=',n1,nr,hbarc
!      write(*,*)'kk=',kk
       do 80 ik=1,n1
       write(kout,*)'mbess, ik=',ik
       do 80 xk=1,nr
       if (kk(ik).gt.1e5)then
       do 82 lb=0,lmx2
       sbess(ik,xk,lb)=0.
 82    continue
       go to 80
       endif
       x = cmplx(kk(ik)*radxis(xk)/hbarc)
       
       call sbesjh(x,lmx2,xj,xjp,xh1,xh1p,ifail)
        do 85 lb=0,lmx2
        sbess(ik,xk,lb) = real(xj(lb))
 85     continue
 80    continue
       return
       end




c *******************************************************************
      subroutine mtmatt(kk,n1,lmas2,ir)
c ********************************************************************
c this subroutine calculates the partial wave decomposition of the
c nucleon-nucleon transition matrix
c *********************************************************************
      implicit real*8(a-h,k,o-z)
      real*8 mn,mp
      complex*16 z,tmatt
      complex*16 amp0,amp1
      integer code,xi
      integer option1, option2
      dimension xis(48),wts(48),kk(50)
      dimension pl(0:99), plcth(0:99,48)
      dimension tapb1(4),te1(4),tapb2(4),te2(4)
      common/params/hbarc,pi,mp,mn,nz,na
      common/ampl/q(200),qq(200),amp0(200,200),amp1(200,200)
      common/btmatt/tmatt(0:7,10,50,50)
      common/block1/option1
      common/block2/option2

      z = (0.,1.)
      lmx = lmas2
      k0 = kk(n1)
      twopi2=2.*pi*pi
      h3 = hbarc**3
c *** forward approximation
      if (lmas2.eq.0) then
      write(6,*)'forward approximation in vopt2nl'
      call redish(k0,k0,k0,1.d0,mn,tapb1,te1)
      do 10 ik=1,n1
      do 10 ip=1,n1
      tmatt(0,ir,ip,ik)=twopi2*(tapb1(1)+z*tapb1(3))*h3
      tmatt(0,ir+5,ip,ik)=twopi2*(te1(1)+z*te1(3))*h3
 10   continue
      return
      endif
c ***
      a = -1
      b = 1
      ncos = 32
      write(*,*)'rc3: ncos=',ncos
      code = 11
      call gauss2(ncos,code,a,b,xis,wts)

c *** initialize nn transition matrix
      do 30 ik=1,n1
      do 30 ip = 1,n1
      do 30 l=0,lmx
      tmatt(l,ir,ip,ik)=0.
      tmatt(l,ir+5,ip,ik)=0.
 30   continue


      do 60 xi=1,ncos
      cthnuc = xis(xi)
      call legpol(cthnuc,pl,lmx+1)
      do 80 ik = 1,n1
      do 80 ip = 1,ik
      k = kk(ik)
      p = kk(ip)
      call redish(k,p,k0,cthnuc,mn,tapb1,te1)
      do 80 l=0,lmx
      tmatt(l,ir,ip,ik)=tmatt(l,ir,ip,ik)
     #    +twopi2*1/2.*pl(l)*wts(xi)*(tapb1(1)+z*tapb1(3))*h3
      tmatt(l,ir+5,ip,ik)=tmatt(l,ir+5,ip,ik)
     #    +twopi2*1/2.*pl(l)*wts(xi)*(te1(1)+z*te1(3))*h3
 80   continue
 60   continue
      do 90 ik=1,n1
      do 90 ip=1,ik
      do 90 l=0,lmx
      tmatt(l,ir,ik,ip)=tmatt(l,ir,ip,ik)
      tmatt(l,ir+5,ik,ip)=tmatt(l,ir+5,ip,ik)
 90   continue
      return
      end
c*************************************************************
      subroutine newgrid(k0,n1)
c*************************************************************
      implicit real*8(a-h,o-z)
      real*8 k0,k0s,kk
      common/bgrid/kk(50),wt(48)
      common/xsng/xs(50,20)
      common/qnumb/lamx,lbmx,l1mx,lfmx

      kode=022
      b=21000.
      ngp=n1-1
      ishell=0
      do 40 la=0,lamx
      do 40 lb=0,lbmx
      ishell=ishell+1
      k0s = sqrt(k0*k0 + eshell(la)-eshell(lb))
        do 60 ip=1,n1
        xs(ip,ishell)=(kk(ip)-k0s)/(k0s+kk(ip))
 60     continue
 40   continue
      return
      end


c***********************************************************
       function rho(l,r)
c**************************************************************
c ***  calculates radial ho wave functions normalised to one
c ***  for oxigen
c**************************************************************
       implicit real*8 (a-h,o-z)
       real*8 mn,mp
       common/sizes/achp,acmp,wsp,achn,acmn,wsn
       common/params/hbarc,pi,mp,mn,nz,na

       exp1 = exp(-(r/sqrt(2.)/achp)**2)
       if (l.eq.0) then
       sn10 = 2.*(1./(pi*achp**6))**(1/4.)
c      sn10 = 1.
       rho= sn10*exp1
       else if (l.eq.1)then
       sn11 = 2.*(4./(9.*pi*achp**10))**(1/4.)
       rho= sn11*r*exp1
       else
       stop
       end if
       return
       end


c ****************************************************************
       subroutine rint2i(k,kp,ik,ikp,n1)
c ****************************************************************
       implicit real*8 (a-h,o-z)
       real*8 k,kp,kk,kkmx
       real*8 mp,mn
       complex*16 xj(0:30),xjp(0:30),xh1(0:30),xh1p(0:30),x
       integer code,xk
       dimension kk(50)
       dimension wt(48)
       common/params/ hbarc, pi, mp, mn, nz, na
       common/sizes/  achp, acmp, wsp, achn, acmn, wsn
       common/radgr/  radxis(201), radwt(201)
       common/bgauss3/rmaxr,quin,mquadi,mquado,irmx
       common/bsdess/sbess(49,201,0:30)
       common/bwf/grho(0:4,201)
       common/qnumb/lamx,lbmx,l1mx,lfmx
       common/brit1/rit1(0:30,0:30,0:2,0:2,50)
       common/brit2/rit2(0:30,0:30,0:2,0:2,50)
       common/bgrid/kk,wt

       ibl = lamx+lbmx
       lmx = l1mx+lfmx+ibl
       lmx2 = lmx
       h2 = hbarc*hbarc
       nr=irmx
c      kkmx=4000.
       kkmx=2300.




       if(k.gt.kkmx.or.kp.gt.kkmx) go to 60
       do 15 ip=1,n1
       if(kk(ip).gt.kkmx)go to 15
c      do 20 la=0,lamx
c      do 20 lb=0,lbmx
       do 20 i=0,lmx
       j1=abs(i-ibl)
       if (i.lt.2)j1=0
c      j1 = 0
       j2 = i+ibl
       do 20 j=j1,j2
c      do 20 ip=1,n1
       do 25 la=0,lamx
       do 25 lb=0,la
       r1 = 0.
       r2 = 0.
       do 30 xk=1,nr
       r1 = r1
     *   +grho(la,xk)*grho(lb,xk)
     *    *sbess(ik,xk,i)*sbess(ip,xk,j)*radwt(xk)*radxis(xk)*radxis(xk)

       r2=r2
     *   +grho(la,xk)*grho(lb,xk)
     *   *sbess(ip,xk,i)*sbess(ikp,xk,j)*radwt(xk)*radxis(xk)*radxis(xk)
 30    continue
       rit1(j,i,lb,la,ip) = r1
       rit2(j,i,lb,la,ip) = r2
 25    continue
c      if (lamx.gt.0.and.lbmx.gt.0) then
       do 40 lb=1,lbmx
       do 40 la=0,lb-1
       rit1(j,i,lb,la,ip) = rit1(j,i,la,lb,ip)
       rit2(j,i,lb,la,ip) = rit2(j,i,la,lb,ip)
  40   continue
c       end if
  20   continue
  15   continue

 60    return
       end


c**********************************************************
      subroutine vcinter(n1,ishell)
c**********************************************************
      implicit real*8(a-h,o-z)
      complex*16 vc,wf,wxxi
      dimension rr(50),wf(50)

      common/qnumb/lamx,lbmx,l1mx,lfmx
      common/vcng/vc(2,50,20)
      common/xsng/xs(50,20)

      ndim = 50
        do 60 ip=1,n1
        rr(ip) = xs(ip,ishell)
        wf(ip) = vc(1,ip,ishell)
c       write(6,*)rr(ip),vc(1,ip,ishell)
 60     continue
        do 80 ip=1,n1-1
        r = xs(ip,1)
        vc(1,ip,ishell)=wxxi(r,rr,wf,ndim,n1-1)
c       write(6,*)'vc=',vc(1,ip,ishell)
 80     continue
        vc(1,n1,ishell)=wxxi(xs(n1,ishell),rr,wf,ndim,n1-1)
c       write(6,*)'vc(1,n1,ishell)=',vc(1,n1,ishell)
      return
      end


 
       subroutine rint2c(k,kp,ik,ikp,n1)
c ****************************************************************
       implicit real*8 (a-h,o-z)
       real*8 k,kp,kk,kkmx
       real*8 mp,mn
       complex*16 xj(0:30),xjp(0:30),xh1(0:30),xh1p(0:30),x
       integer code,xk
       dimension kk(50)
       dimension wt(48)
       dimension fj(0:2,0:2,201)
       common/params/ hbarc, pi, mp, mn, nz, na
       common/sizes/  achp, acmp, wsp, achn, acmn, wsn
       common/radgr/  radxis(201), radwt(201)
       common/bgauss3/rmaxr,quin,mquadi,mquado,irmx
       common/bsdess/sbess(49,201,0:30)
       common/bwf/grho(0:4,201)
       common/qnumb/lamx,lbmx,l1mx,lfmx
       common/brit1/rit1(0:30,0:30,0:2,0:2,50)
       common/brit2/rit2(0:30,0:30,0:2,0:2,50)
       common/bgrid/kk,wt

       ibl = 2
       h2 = hbarc*hbarc
       lmx = l1mx+ibl+lfmx
       lmx2 = lmx
       nr=irmx
c      kkmx = 2300.
       kkmx = 1500.
c      del = 5.0*hbarc
       del = 2.0*hbarc

       if(na.eq.16)then
       do 5 xk=1,nr
       fj(0,0,xk) = (grho(0,xk)**2 - grho(1,xk)**2)  /sqrt(8.)
       fj(0,1,xk) = grho(0,xk)**2/sqrt(6.)
       fj(0,2,xk) = grho(1,xk)**2/sqrt(2.)
       fj(1,0,xk) = grho(0,xk)*grho(1,xk)
       fj(2,0,xk) = grho(1,xk)*grho(1,xk)
c      g02 = grho(0,xk)**2
c      g12 = grho(1,xk)**2
c      fj(0,0,xk) = (g02 - g12) /sqrt(8.)
c      fj(0,1,xk) = g02/sqrt(6.)
c      fj(0,2,xk) = g12/sqrt(2.)
c      fj(1,0,xk) = grho(0,xk)*grho(1,xk)
c      fj(2,0,xk) = g12
 5     continue
       endif

        do 10 ip=1,n1
        do 10 jc=0,2
        do 10 jsf=0,2
        do 10 i=0,lmx
        do 10 j=abs(i-jc),i+jc,2
c       do 10 j=0,lmx
        rit1(j,i,jc,jsf,ip)=0.
        rit2(j,i,jc,jsf,ip)=0.
 10     continue


       if(k.gt.kkmx.or.kp.gt.kkmx) go to 60
       do 15 ip=1,n1
       dk = abs(kk(ip)-k)
       dkp = abs(kk(ip)-kp)
       if(dk.gt.del.and.dkp.gt.del)go to 15
       if(kk(ip).gt.kkmx)go to 15
       do 20 jc=0,2
       jsfmx = 0
       if(jc.eq.0)jsfmx=2
       do 20 jsf=0,jsfmx
       do 20 i=0,lmx
       do 20 j=abs(i-jc),i+jc,2
c      do 20 j=0,lmx
       r1 = 0.
       r2 = 0.
       do 30 xk=1,nr
       raux = fj(jc,jsf,xk)*radwt(xk)*radxis(xk)*radxis(xk)
       r1 = r1 + raux
     *    *sbess(ik,xk,i)*sbess(ip,xk,j)

       r2 = r2 + raux
     *   *sbess(ip,xk,i)*sbess(ikp,xk,j)
 30    continue
       if (dk.gt.del) r1=0.
       if (dkp.gt.del) r2=0.
       rit1(j,i,jc,jsf,ip) = r1
       rit2(j,i,jc,jsf,ip) = r2
 20    continue
 15    continue

 60    return
       end

c **********************************************************
      subroutine vopt2nlc(v2l,k,kp,ik,ikp,n1,lmas2)
c **********************************************************
c *   this subroutine calculates the second order coupling *
c *   potentials for oxygen using closure                  *
c **********************************************************
      implicit real*8(a-h,k,o-z)
      real*8 mn,mp
      real*8 mnuc,mr
      integer option6
      integer kcall
      dimension v2l(0:99,2,6)
      dimension kk(50),wt(48)
      dimension vaux1(2),vaux2(2)
c     dimension f2(10),zexp(0:110),onexp(0:110)
      complex*16 v1c,v2c,z,vaux
      complex*16 z1,z2,z3,z4,zf,z14,z23
      complex*16 zexp
      complex*16 vcon(20),vc
      complex*16 vaux1,vaux2
      complex*16  tmatt
      common/params/ hbarc, pi, mp, mn, nz, na
      common/bdeltav/ideltav
      common/block6/option6
      common/bgrid/kk,wt
      common/bprop/deng(50,10)
      common/brit1/rit1(0:30,0:30,0:2,0:2,50)
      common/brit2/rit2(0:30,0:30,0:2,0:2,50)
      common/qnumb/lamx,lbmx,l1mx,lfmx
      common/bkcall/kcall
      common/clebma/faclog(500)
      common/btaux/t1(100000)
      common/btmatt/tmatt(0:7,10,50,50)
      common/ffch/ich
      common/bsflip/isflip
      common/bjcount/jcount(0:3,0:99)
      common/bfast/f2(10),zexp(0:110),onexp(0:110)
c
      data g1,h1,k1,g2,h2,k2/6*0.d0/

      rewind 15
      z = cmplx(0.,1.)
      mr = (mn*na)*mp/(mn*na+mp)
      mnuc = mn*na
      h3 = hbarc**3
c     del = 5.0*hbarc
      del = 2.0*hbarc
c     kkmx = 2300.
      kkmx = 1500.
      k0=kk(n1)
      ideltav=0
      ipst=1
      if(ideltav.eq.1)ipst=n1
      tfact = -(na-1.)
      tfact1 = tfact * 3./2.
      tfact2 = tfact * 24./(na*1.)
      vaux1(2) = 0.
      vaux2(2) = 0.
c ***   mulligan factor
c      if(ich.eq.1)tfact1 = tfact1*(0.8)**2

c *** generates the radial integrals
      call rint2c(k,kp,ik,ikp,n1)

c *** initialize partial wave 2nd order optical potential
      do 5 ld=0,lmas2
      v2l(ld,1,1)=0.
      v2l(ld,2,1)=0.
      v2l(ld,1,2)=0.
      v2l(ld,2,2)=0.
 5    continue
c
      do 6 l1=0,lmas2+lfmx+2
      zexp(l1) = z**l1
      onexp(l1) = (-1)**l1
 6    continue
c
      do 7 jamp=2,10
      if(jamp.le.5)then
      f2(jamp)=1.
      if(jamp.eq.3)f2(jamp)=f2(jamp)*2.
      else
      f2(jamp)=3.
      if(jamp.eq.8)f2(jamp)=f2(jamp)*2.
      endif
 7    continue
c *** loop over partial wave
c **** at the moment only l1mx partial waves are calculated
      lmas2n = l1mx
c     lmas2n = lmas2
      do 20 ld=0,lmas2n
      if (kcall.eq.1)then
c *** calculates and stores all cleb and racah
      im = 0
      do 10 jc=0,2
      jcount(jc,ld) = im
      do 10 lc=abs(jc-ld),(jc+ld)
      do 15 lf =0,lfmx
      c1 = 2*lf+1.
      do 15 l1=abs(ld-lf),ld+lf,2
      do 15 l2=abs(l1-jc),l1+jc,2
      f1=clebz(real(jc),real(l2),real(l1))
     #        *  hat(l1)*hat(l2)**2/hat(ld)
         im = im+1
         g1 = rac(real(jc),real(l2),real(ld),real(lf),real(l1),
     #         real(lc))

         h1 = 0
         k1 = 0
c          if(g1.ne.0.)then
           if(abs(g1).gt.1e-6)then
           h1 = clebz(real(lf),real(l1),real(ld))
c          if(h1.ne.0.)then
           if(abs(h1).gt.1e-6)then
           k1 = clebz(real(lf),real(l2),real(lc))
           end if
           end if
         t1(im) = c1*f1*g1*h1*k1
 17      continue
 15      continue
 10      continue

       write(15)im
       write(15)(t1(i),i=1,im)
       else
       read(15)im
       read(15)(t1(i),i=1,im)
       endif
      vaux1(1) = 0.d0
      vaux2(1) = 0.d0
      vaux1(2) = 0.d0
      vaux2(2) = 0.d0
c *** loop in the intermediate momentum grid
      do 30 ip=ipst,n1
      p = kk(ip)
      if (p.gt.kkmx)go to 30
      dk=abs(k-p)
      dkp=abs(kp-p)
      if(dk.gt.del.or.dkp.gt.del)go to 30
c *** loop in all combinations of quantum numbers
      im=0
      do 40 jc=0,2
      do 40 lc=abs(jc-ld),jc+ld
      v1c=0.d0
      v2c=0.d0
        do 50 lf=0,lfmx
        z1 = tmatt(lf,1,ik,ip)
        z2 = tmatt(lf,1,ip,ikp)
        do 50 l1=abs(ld-lf),ld+lf,2
c       z3 = z**l1
        z3 = zexp(l1)
        z4 = z3*onexp(l1)
        z14 = z1*z4
        z23 = z2*z3
        do 50 l2=abs(l1-jc),l1+jc,2
        im = im+1
c       if (t1(im).ne.0.)then
        if (abs(t1(im)).gt.1e-6)then
c       zf = t1(im) * z**l2
        zf = t1(im) * zexp(l2)
        v1c = v1c +
     #                  zf*rit1(l2,l1,jc,0,ip)
     *                  *z14
        v2c = v2c +
     #                   zf*rit2(l1,l2,jc,0,ip)
     *                   *onexp(l2)*z23
        end if
 50     continue
       vc=v1c*v2c
       vaux1(1) = vaux1(1) + (-1.)*deng(ip,1) * vc
       if (ip.eq.n1)then
       vaux2(1)=vaux2(1) - z*2.*mr* k0*vc
       endif
 40    continue
c
c ***    introduction of the spin-flip term
      if (isflip.eq.1)then
      do 60 jsf =0,2
        if (jsf.eq.0)then
        jci = 1
        jcf = 2
        else if (jsf.gt.0)then
        jci = 0
        jcf = 0
        endif
      do 60 jamp=2,10
      im = 0
      do 60 jc=jci,jcf
      im = jcount(jc,ld)
      do 60 lc=abs(jc-ld),jc+ld
      v1c=0.d0
      v2c=0.d0
        do 70 lf=0,lfmx
        z1 = tmatt(lf,jamp,ik,ip)
        z2 = tmatt(lf,jamp,ip,ikp)
        do 70 l1=abs(ld-lf),ld+lf,2
        z3 = zexp(l1)
        z4 = z3*onexp(l1)
        z14 = z1*z4
        z23 = z2*z3
        do 70 l2=abs(l1-jc),l1+jc,2
        im = im+1
c       if (t1(im).ne.0.)then
        if (abs(t1(im)).gt.1.e-6)then
        zf = t1(im) * zexp(l2)
        v1c = v1c +
     #                  zf*rit1(l2,l1,jc,jsf,ip)
     *                  *z14
        v2c = v2c +
     #                   zf*rit2(l1,l2,jc,jsf,ip)
     *                   *onexp(l2)*z23
        end if
 70     continue
       vc=v1c*v2c
       vaux1(2) = vaux1(2) + (-1.)*deng(ip,1) * vc*f2(jamp)
       if (ip.eq.n1)then
       vaux2(2)=vaux2(2) - z*2.*mr* k0*vc*f2(jamp)
       endif
 60    continue
      endif
c ***
c
 30    continue

c      if (ld.lt.20)then
c      write(6,*)ld,k,kp,vaux1(1),vaux2(1)
c      end if
       if(ideltav.lt.1)then
       vaux = (vaux1(1)+vaux2(1))*tfact1
     *      + (vaux1(2)+vaux2(2))*tfact2
        else
        vaux=vaux2(1)*tfact1 + vaux2(2)*tfact2
        endif

       v2l(ld,1,1)= real(vaux)/h3/h3
       v2l(ld,2,1)= dimag(vaux)/h3/h3
       v2l(ld,1,2)=v2l(ld,1,1)
       v2l(ld,2,2)=v2l(ld,2,1)
 20    continue
c
      return
      end





c **********************************************************
      subroutine vnl2(v2l,k,kp,ik,ikp,n1,lmas2)
c **********************************************************
c *   this subroutine calculates the second order coupling *
c *   potentials for oxygen without using closure          *
c *   second component
c **********************************************************
      implicit real*8(a-h,k,o-z)
      real*8 mn,mp
      real*8 mnuc,mr
      integer option6
      dimension kgama(20),den(50)
      dimension ri0(0:2,0:2,0:30,0:30)
      dimension v2l(0:99,2,6)
      dimension kk(50),wt(48)
      complex*16 v1c,v2c,z,vaux
      complex*16 vcon(20),vc
      complex*16 vaux1,vaux2
      complex*16 v2,vrc1(50),vrc2(50)
      complex*16  tmatt
      common/params/ hbarc, pi, mp, mn, nz, na
      common/block6/option6
      common/bgrid/kk,wt
      common/bwf/grho(0:4,201)
      common/bprop/deng(50,10)
      common/rcv2/v1c(50,20),v2c(50,20)
      common/vcng/vc(2,50,20)
      common/brit1/rit1(0:30,0:30,0:2,0:2,50)
      common/brit2/rit2(0:30,0:30,0:2,0:2,50)
      common/qnumb/lamx,lbmx,l1mx,lfmx
      common/clebma/faclog(500)
c     common/bv2/v2l(0:99,2,6)
c     common/btaux/t1(100000),t2(100000)
      common/btaux/t1(100000)
      common/btmatt/tmatt(0:7,10,50,50)
c
      data g1,h1,k1,g2,h2,k2/6*0.d0/

      z = cmplx(0.,1.)
      mr = (mn*na)*mp/(mn*na+mp)
      mnuc = mn*na
      h3 = hbarc**3
      k0=kk(n1)

c *** initialize partial wave 2nd order optical potential
      do 5 ld=0,lmas2
      v2l(ld,1,1)=0.
      v2l(ld,2,1)=0.
 5    continue
c
      ishell=0
c *** calculates on-shell shift momentum
      do 7 la=0,lamx
      do 7 lb=0,lbmx
      ishell = ishell + 1
      kgama(ishell)=sqrt(k0*k0+eshell(la)-eshell(lb))
 7    continue
c
c *** loop over partial wave
c     lmas2n = l1mx+lfmx
c     lmas2n = lmas2
      lmas2n = l1mx
      do 200 ld=0,lmas2n

c *** calculates and stores all cleb and racah
      im = 0
      do 10 la=0,lamx
      do 10 lf=0,lfmx
      do 10 l1=abs(ld-lf),ld+lf,2
c     do 10 l1=0,l1mx
         t1aux = 0.
         do 15 lc=abs(la-l1),la+l1
         t1aux = t1aux+clebz(real(la),real(l1),real(lc))**2
 15      continue
c        do 20 lf=0,lfmx
         im = im+1
         t1(im) = t1aux * (2*la+1.)*(2*l1+1)*(2*lf+1.)*
     *              clebz(real(l1),real(lf),real(ld))**2
     *              /(2*ld+1.)
 20      continue
 10      continue

      vaux = 0.d0
      vaux1 = 0.d0
      vaux2 = 0.d0
c *** loop in the intermediate momentum grid
      do 30 ip=1,n1
      p = kk(ip)
c *** loop in all combinations of quantum numbers
      im=0
      do 40 la=0,lamx
      v1c(ip,la)=0.d0
      v2c(ip,la)=0.d0
        do 50 lf=0,lfmx
        do 50 l1=abs(lf-ld),lf+ld,2
c       do 50 l1=0,l1mx
c       do 50 lf=0,lfmx
        im = im+1
        if (t1(im).ne.0.)then
        v1c(ip,la) = v1c(ip,la) +
     #                  t1(im)*rit1(l1,l1,la,la,ip)
     *                  *tmatt(lf,1,ik,ip)

        v2c(ip,la) = v2c(ip,la) +
     #                  t1(im)*rit2(l1,l1,la,la,ip)
     *                  *tmatt(lf,1,ip,ikp)
        end if
 50      continue
 40     continue
 30    continue

       i=0
       do 80 la=0,lamx
       do 80 lb=0,lbmx
       i = i+1
         do 85 ip=1,n1
         vc(1,ip,i)=v1c(ip,la)*v2c(ip,lb)
  85     continue
c        if (option6.eq.2.and.abs(la-lb).ne.0)then
c        call vcinter(n1,i)
c        endif
       do 90 ip=1,n1
       vaux1=vaux1+(-1.)*deng(ip,1)*vc(1,ip,i)
  90   continue
c      vaux2=vaux2 - z*2.*mr* kgama(i)*vc(1,n1,i)
       vaux2=vaux2 - z*2.*mr* k0*vc(1,n1,i)
 80    continue
       if (ld.eq.0)then
c      write(6,*)ld,k,kp,vaux1,vaux2
       end if
       vaux = vaux1+vaux2
       v2l(ld,1,1)=(na-1.)/(na*1.)/(na*1.)*real(vaux)*16./h3/h3
       v2l(ld,2,1)=(na-1.)/(na*1.)/(na*1.)*dimag(vaux)*16./h3/h3
       v2l(ld,1,2)=v2l(ld,1,1)
       v2l(ld,2,2)=v2l(ld,2,1)
 200   continue
c
      return
      end

c***********************************************************
       function rhochk(l,r)
c**************************************************************
c ***  calculates radial ho wave functions normalised to one
c ***  for oxigen
c**************************************************************
       implicit real*8 (a-h,o-z)
       real*8 mn,mp
       common/sizes/achp,acmp,wsp,achn,acmn,wsn
       common/params/hbarc,pi,mp,mn,nz,na

        achk = 1.3*achp
       exp1 = exp(-(r/sqrt(2.)/achk)**2)
       if (l.eq.0) then
       sn10 = 2.*(1./(pi*achk**6))**(1/4.)
c      sn10 = 1.
       rhochk= sn10*exp1
       else if (l.eq.1)then
       sn11 = 2.*(4./(9.*pi*achk**10))**(1/4.)
       rhochk= sn11*r*exp1
       else
       stop
       end if
       return
       end

