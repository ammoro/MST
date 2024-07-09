c ***************************************************************
c    program LPOTPS
c    R. Crespo, R. C. Johnson, and J. A. Tostevin
c    Phys. Rev. C 41, 2257 (1990)

c ** f90 version by R. Crespo and A.M. Moro  (spring 2002)
c ***************************************************************


c ** Calculates elastic scattering wavefunction as partial wave 
c ** expansion
c ***********************************************************
      subroutine wavfn(l,nsp,n1,xgama,lmax)
c ***********************************************************
      use params
      use bscoul
      use switch
      use ranges
      use bgrid
      use bcwfn
      use tcomn
      use rcomn
      use foptp
      use radwf
      use btaux
      use bsdess

      implicit real*8(a-h,k,o-z)
      real*8 ::jj
      complex*16 ::znorm, zomega, zi, zaux ! tcwf
      complex*16 ::rmat,fmat,ton
      complex*16,allocatable,dimension(:)::xj,xjp,xh1,xh1p

!      dimension nifty(20)
!      dimension scoul(30)
!      dimension tbr(100,12),tbi(100,12),tr(100,12),ti(100,12)

      real*8,allocatable,dimension(:,:)::reul,aimul

      allocate(reul(npoints+1,n1))
      allocate(aimul(npoints+1,n1))
      
      
!      dimension reul(50,30,2),aimul(50,30,2)
     


!      common/params/hbarc,pi,mp,mn,nz,na,nes,nwaves
!      common/bscoul/scoul
!      common/switch/nifty
!      common/ranges/rcoul,rcut
!      common/bgrid/kk(50),wt(50)
!      common/bcwfn/sa(30),fc(30),fpc(30),gc(30),gpc(30)
!      common/tcomn/tr,ti,tbr,tbi
!      common/rcomn/rr,ri,rrb,rib
!      common/foptp/fr(98,98),fi(98,98)
!      common/radwf/igwf,lstore,rmax,stepr
!      common/btaux/tcwf(25000,2)
!      common/bsdess/sbess(49,161,0:30)

      if (allocated(tcwf)) deallocate(tcwf)
      allocate(tcwf(0:lmax,2))
      allocate(sbess(n1,npoints+1,0:lmax))
!      allocate(fc(0:ldum),gc(0:ldum))
!      allocate(xj(0:ldum),xjp(0:ldum),xh1(0:ldum),xh1p(0:ldum))

      zi = (0.,1.)
      npoints = (rmaxwf-stepr)/stepr+1
      if (l.eq.1)then
      write(12,700)lstore,npoints,stepr
      endif
      eta = xgama
      if(nifty(10).lt.3)eta=0.
      kkmx=10.*hbarc
c **  calculates on-shell t matrix extracting the coulomb phase
      sigl = scoul(l)
      ton = tr(l,nsp) + zi*ti(l,nsp)
      if(nifty(10).ge.3)then
      ton = ton*(cos(2*sigl)-zi*sin(2*sigl))
      endif

! Better in main????
c *** calculates all spherical bessel and coulomb functions
      if (l.eq.1)then
      do 30 npt=1,n1
      do 30 i=1,npoints+1
      do 30 ldum=0,lmax
      sbess(npt,i,ldum)=0.
 30   continue
      do 40 npt=1,n1
      do 40 i=1,npoints+1
      r = stepr*i
      if (i.gt.npoints)r=rcut
      rhl = kk(npt)*r/hbarc
      zaux = cmplx(rhl)
      call sbesjh(zaux,lmax,xj,xjp,xh1,xh1p,ifail)
      do 50 ldum=0,lmax
      sbess(npt,i,ldum) = real(xj(ldum))
      if(ldum.eq.0.and.i.ge.npoints)then
       endif
 50   continue
 40   continue
      rhl=rcut*kk(n1)/hbarc
      call coulfn(eta,rhl,stepr,lmax)
      do 60 ldum =1,lmax
      tcwf(ldum,1)=fc(ldum)
      tcwf(ldum,2)=gc(ldum)+zi*fc(ldum)
 60   continue
      endif
c
       rmat = rr + zi * ri
      if(nifty(10).ge.3)then
      znorm = 0.
c *** calculates denominator
      do 70 npt=1,n1
      fj = sbess(npt,npoints+1,l-1)
      fmat = fr(n1,npt) + zi*fi(n1,npt)
      zomega = na*fmat/(1-zi*rmat)/(na-1.)
      if (npt.eq.n1)zomega = zomega-1/(na-1.)
      znorm = znorm+fj*zomega
 70   continue

c *** calculates numerator
      rhl=kk(n1)*rcut/hbarc
      znorm = ( tcwf(l,1) + ton*tcwf(l,2))/rhl
     *              /znorm
c **  multiplies by the coulomb phase
      znorm = znorm*(cos(sigl)+zi*sin(sigl))
      else
      znorm = (1.,0.)
      endif
c     write(12,*)'znorm=',znorm
c
c *** calculates wave function for l up to lstore
c *** do loop over grid points
      do 130 i=1,npoints
      r = stepr*i
      if(r.le.rcut)then
      zaux=0.
      do 150 npt=1,n1
      fmat = fr(n1,npt) + zi*fi(n1,npt)
      fj = sbess(npt,i,l-1)
      zomega = na*fmat/(1-zi*rmat)/(na-1.)
      if (npt.eq.n1)zomega = zomega-1/(na-1.)
      zomega=zomega*fj
      zaux=zaux+zomega
 150  continue
      zaux = zaux*znorm
      else
      rhl=kk(n1)*r/hbarc
      call coulfn(eta,rhl,stepr,l)
      zaux = fc(l) + ton*(gc(l)+zi*fc(l))
      zaux = zaux/rhl
      endif
!      reul(i,l,nsp)=r*kk(n1)/hbarc*real(zaux)
!      aimul(i,l,nsp)=r*kk(n1)/hbarc*dimag(zaux)
      reul(i,nsp)=r*kk(n1)/hbarc*real(zaux)
      aimul(i,nsp)=r*kk(n1)/hbarc*dimag(zaux)
 130  continue

      ll=l-1
      if(nsp.eq.1)jj=ll+1/2.
      if(nsp.eq.2)jj=ll-1/2.
      write(12,730)ll,jj
c     do 160 i=1,npoints
c     r=stepr*i
c     write(12,777)r, reul(i,l,nsp),aimul(i,l,nsp)
c160  continue
!      write(12,750)(reul(i,l,nsp),aimul(i,l,nsp),i=1,npoints)
      write(12,750)(reul(i,nsp),aimul(i,nsp),i=1,npoints)
c ***
c     write(12,*)'asymptotic wave functions'
c     do 200 i=1,npoints
c     r = rcut+i*stepr
c     rhl=kk(n1)*r/hbarc
c     call coulfn(eta,rhl,stepr,l)
c     zaux = fc(l) + (tr(l,nsp)+zi*ti(l,nsp))*(gc(l)+zi*fc(l))
c     zaux = zaux/rhl
c     zaux=zaux*r*kk(n1)/hbarc
c     write(12,777)r,real(zaux),dimag(zaux)
c200  continue
 700  format(2i4,f8.4)
 730  format(i4,f8.4)
 750  format(6e12.4)
 777  format(f7.2,2e12.4)

      return
      end

      subroutine coulfn (eta,rhoz,drho,lmax)
      implicit real*8 ( a-h , o-z )
      common/bcwfn/ sa(30),fa(30),fpp(30),ga(30),gpp(30)
      dimension gd(5)
      data cg0, cg1, cg2, cg3, cg4, cg5/
     +1.223404016d0,4.959570165d-2,8.888888889d-3,
     +2.455199181d-3,9.108958061d-4,2.534684115d-4/
      data cgp,cp1,cp2,cp3,cp4,cp5 /
     +-7.078817734d-1,1.728260369d-1,3.174603174d-4,
     +3.581214850d-3,3.117824680d-4,9.073966427d-4/
      data gd(1),gd(2),gd(3),gd(4),gd(5) / 12.,-360.,1260.,-1680.,1188./
      nr=1
      epsi=0.d0
      eps=0.1*epsi
      eps=max (eps, 1.0d-11)
      eta2=eta+eta
      etas=eta*eta
      lp=max0 (lmax,12)+1
      t=lp
      u=t*t+etas
      v=t/u
      w=eta/u
      x=v*v-w*w
      y=2.*v*w
      u=sqrt (u)
      sig0=eta*(dlog(u)-1.)+(t-0.5)* datan(eta/t)
      do 20 i=1,5
      sig0=sig0-w/gd(i)
      t=v*x-w*y
      w=v*y+w*x
      v=t
   20 continue
   30 if (lp .le. lmax+1) sa(lp)=sig0
      lp=lp-1
      if (lp .le. 0) go to 100
      t=lp
      sig0=sig0- datan (eta/t)
      go to 30
  100 emax=(1.0e-5/eps)**0.16666667
      if (eta .lt. emax) go to 200
      r=eta2
      t=6.0
      t=eta**(1.0/t)
      w=eta*t
      u=t-t*(cg2+cg4/etas)/etas
      v=(cg1+(cg3+cg5/etas)/etas)/w
      g=cg0*(u+v)
      t=1./t
      w=eta*t
      u=t+t*(cp2+cp4/etas)/etas
      v=(cp1+(cp3+cp5/etas)/etas)/w
      gp=cgp*(u-v)
      go to 300
  200 r=max0 (nr,1)-1
      r=rhoz+r*drho
      t=12.+1.4*eta
      if (t .gt. r) r=t
      fk=1.
      f=1.
      gk=0.
      g=0.
      fsk=0.
      fp=0.
      gsk=1.-eta/r
      gp=gsk
      epss=eps*eps
      n=r+r
      do 210 kp=1,n
      t=kp+kp
      u=t*r
      ak=(t-1.)*eta/u
      v=kp*(kp-1)
      bk=(etas-v)/u
      t=ak*fk-bk*gk
      gk=ak*gk+bk*fk
      fk=t
      t=ak*fsk-bk*gsk-fk/r
      gsk=ak*gsk+bk*fsk-gk/r
      fsk=t
      f=f+fk
      g=g+gk
      fp=fp+fsk
      gp=gp+gsk
      test=fk*fk+gk*gk+fsk*fsk + gsk*gsk
      if (test .lt. epss) go to 220
  210 continue
  220 t=r-eta*dlog(r+r)+sig0
      u=cos (t)
      v=sin (t)
      g=f*u-g*v
      gp=fp*u-gp*v
  300 rs=r
      rho=rhoz
      f=g
      fp=gp
      is=0
      ir=1
      t=r-rho
      if (t) 600,700,310
  310 if (nr .le. 1) go to 320
      is=t/drho
      is=min0 (is+1,nr)
  320 t=is
      rho=rhoz+t*drho
      gg=rho
      is=is+1
      ir=is
  330 rho=rho-drho
      ir=ir-1
      if (ir .gt. 0) go to 600
      ir=max0 (is,1)
      r=rs
      rho=gg
      f=g
      fp=gp
      go to 350
  340 rho=rho+drho
      ir=ir+1
  350 if (ir .gt. nr) return
  600 h=0.5
      w=r-eta2
      if (r-1.0) 601,602,602
  601 h=0.5*r
  602 if (w) 603,605,605
  603 t=sqrt (-r/(w+w))
      if (t-h) 604,605,605
  604 h=t
  605 last=0
      t=rho-r
      if (t) 606,700,607
  606 h=-h
  607 u=t-h
      if (u*h) 608,608,609
  608 h=t
      last=1
  609 u=0.0
      t=1.0
      b1=0.0
      b2=f
      b3=h*fp
      f=f+b3
      v=0.0
  610 it=0
  620 v=-h*(h*b1+w*b2+u*v)/(r*t)
      fp=fp+v
      u=t
      t=t+1.0
      b1=b2
      b2=b3
      b3=h*v/t
      f=f+b3
      test=b3
      testp=v
      if (w) 630,640,640
  630 test=b3/f
      testp=v/fp
  640 if (abs(test)+abs(testp)-eps) 650,610,610
  650 if (it) 660,660,670
  660 it=1
      go to 620
  670 r=r+h
      if (last) 600,600,700
  700 k=lmax+1
      x=f
      y=fp
      do 710 j=1,k
      ga(j)=x
      al=j
      t=j*j
      u=t/rho+eta
      v=sqrt (t+etas)
      w=(u*x-al*y)/v
      y=(v*x-u*w)/al
      x=w
  710 continue
      lp=rho
      lp=max0 (lp+10, lmax+20)
      b3=0.
      b2=1.0e-20
      w=1.0/rho
      al=lp+1
      v=eta/al
      u=0.
      do 840 j=1,lp
      k=lp+1-j
      al=k
      t=k+k+1
      b1=t*(v/al+w)*b2-u*b3
      v=eta/al
      u=sqrt (1.0+v*v)
      b1=b1/u
      b3=b2
      b2=b1
      if (k-lmax-1) 810,810,820
  810 fa(k)=b1
      go to 840
  820 test=b1
      if (abs(test)-1.) 840,840,830
  830 b2=b2*1.0d-20
      b3=b3*1.0d-20
  840 continue
      t=(w+eta)*b2-u*b3
      u=1./(t*f-b1*fp)
      k=lmax+1
      do 850 j=1,k
      fa(j)=u*fa(j)
  850 continue
      do 400 l=1,k
      fl=l
      flsq=fl*fl
      fac1=eta/fl+fl/rhoz
      fac2=dsqrt(eta*eta+flsq)/fl
      fpp(l)=fac1*fa(l)-fac2*fa(l+1)
      gpp(l)=fac1*ga(l)-fac2*ga(l+1)
  400 continue
      if (ir-is) 330,340,340
      end
