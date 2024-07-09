c ** These procedures have disappeared in the f90 version
c ** or should disappear in forthcoming versions


c *** energies, momenta, and masses in mev units throughout program

c *** input switches - nifty(1) thru nifty(20) follow
c
c     nifty(1)     = -2, (p,n) cex after p,n - not used
c                  =  5, protons
c                  =  9, neutron

c     nifty(2)*5mev = binding energy (-5. * nifty(2) used) for 3-body
c                     energy

c     nifty(3)     = 0, default case--takes no action.
c                    1, single scatt (when n(6)used else)
c                    2, set t=r (nonunitarize) but with ms
c                    3, double scattering for 1/2*1/2

c     nifty(4)     = 0, not used by the proton code

c     nifty(5)     = 0, use softh(k0)
c                    1, sin(theat)
c                    2, so(on shell),no angle tf
c                    5, e3b energy
c                    6, e3b energy with folded tmatrice
c                    7, e3b with aay magic vector prescription for theta
c                    9, soth(ko), aay

c     nifty(6)     = 0, nospin flip and no calc of ds and ts (0 x 0)
c                    1, calculate r matrix up to and include ds
c                    2, calculate r matrix up to and include ts
c                    3, include spin flip terms (0 x 1/2)
c                    4, no scatt (u=0) as check of wf
c                    5, r matrix for only ss
c                    6, ustrong = 0
c                    7, imaginary part of u = 0
c                    8, spin 1/2 x 1/2 with 6 spin states - not used

c     nifty(7)     = 0, set h, i, j, and k nn-partial waves to zero
c                  = 1, retain all nn-partial waves.

c     nifty(8)     = 0, don't read in phases, use stored t - not used
c                  = 1, read in phases - not used'

c     nifty(9)     = 0, not used by the proton code.

c     nifty(10)    = 0, no coulomb amplitude
c                    1, coulomb amplitude and phase change included
c                    2, coulomb amplitude, but not phase change
c                    4, exact coulomb via uniform charge sphere
c                    3, exact coulomb with realistic charge distribution

c *** nifty(11) may be used as temporary switch
c     nifty(11)    = 0, schroedinger theory
c                  = 1, relativistic (dirac equation) theory

c          (12)    = 0, not used for proton code
c          (13)    = 0, not used for proton code
c          (14)    = 0, not used for proton code

c     nifty(15)    = 0, divide out proton ff from he ff
c                    1, dont divide out proton ff

c     nifty(16)    = 0, modified gaussian rho and rho**2 & for p* z=n
c                       nucleus will also use t(spin)*rho(matter)
c                    1, wood saxon rho,..but gaussian rho**2
c                    2, wood saxon rho and rho**2
c                    3, special form for he or deuteron
c                    4, he3..fnsp=fnmat,fpsp=0
c                    5, he3..fnsp=fnmat=fpmt, fpsp=0
c                    6, use fmatter*tspin for 0*.5(non-gaussian nucleus)
c                    7, uspin=(tsp*fsp+tspn*fmt) for 0*.5

c     nifty(17)    = 0, kmt potential
c                    1, watson potential, no a/(a-1)

c     nifty(18)    = 0, dont turn off any 1/2*1/2 amps
c                    1, no e terms
c                    2, no a+b terms
c                    3, no e and no a+b terms

c     nifty(19)    = 0, nifty(16) controls ff's '
c                    1, hadjimichael impulse approximation form factors
c                       (no mec)
c                    2, full hadjimichael (i.e. with exchange currents)

c     nifty(20)    = 0, reduced output, no plot file to tape10
c                  = 1, full output, no plot file to tape10
c                  = 2, reduced output, write plot file to tape10
c                  = 3, full output, write plot file to tape10

c     nsex =0-normal/ =1-1st time thru/ =2-2nd time thru/ =3 diff cex

c ** plot up lab cross sections if nang was read in as a negative number
c ** u(nk,1) = uplus for spin flip,  u(nk, 2)= uminus
c ** set = for nonflip (same code for tr,ti,trc,tic,tbr,tbi)
c
cc** new options in lpotps
c** option1 = 1     love-franey (not yet)
c**         = 2     redish
c**         = 3     martin
c** option2 = 1     on-shell
c**         = 2     off-shell
c** option3 = 1     vopt1
c**         = 2     vopt1 + vopt2(local)
c**         = 3     vopt1 + vopt2 (non local)
c**         = 4     vopt1 + vopt2 (full nonlocal)
c** option4 = 1     factorization approximation
c**         = 2     full-folding
c** option5 = 1     ho density for factorized approximation
c**         = 2     ws2pf density for factorized approximation
c**         = 3     ho with different s and p parameters
c**         = 4     numerical density
c** option6 = 1     closure in second order term
c**         = 2     no closure in 2nd order term
c** n.b. option 6 is to be used if option4.ge.3
c** n.b. option 4 can be used only with redish tnn
c       ----------------------------------------------


c @(#)logplt.f 2.2 8/25/86 12:28:25

!!$      subroutine logplt (x,y,ye,n,ne,chart,yminx)
!!$************************************************************************
!!$
!!$c *** semi-log graphing on printer,from pirk,modified by rhl
!!$c *** x(i),y(i) are cartesian coordinates to be plotted
!!$c *** on standard size 5cycle graph paper
!!$c *** ne gt 0 also causes ye to be plotted on same graph
!!$c *** ne lt o for non log polar plot
!!$
!!$      implicit real*8 (a-h, o-z)
!!$
!!$      dimension x(1), y(1), chart(61,61), xf(11), ye(1)
!!$      dimension plot(10), nplot(30)
!!$
!!$c *** ymin must be lowest power of 10, e.g. -2.,or 3.
!!$c *** but not smallest value of y
!!$c ***
!!$
!!$       character*1 :: blank=" ",point=".",aminus="-",ex="0"
!!$       real*8 :: 
!!$       integer :: nit=0
!!$!,point=' ',aminus=' ',ex='0',plus='+'
!!$
!!$!      data plus, pt, blank, ex      /1h+, 1h., 1h , 1h0/
!!$!      data aminus, point            /1h-, 1h+/
!!$!      data nit                      /0/
!!$
!!$!      if (nit .ne. 0) go to 20
!!$
!!$      do 10 i=1,29
!!$         j = i - 7
!!$         if (i .le. 7) plot(i) = chart(i,1)
!!$         if (i .gt. 7) nplot(j) = chart(i,1)
!!$ 10   continue
!!$
!!$ 20   nit = 1
!!$      xmin = 0
!!$      xmax = 180.
!!$      ymin = yminx
!!$      ymax = ymin + 5.
!!$
!!$c *** polarization
!!$
!!$      if (ne .lt. 0) ymin = -1.
!!$      if (ne .lt. 0) ymax = 1.5
!!$      rs = 60./(ymax - ymin)
!!$      mm = 6
!!$
!!$c *** construct graph
!!$
!!$      do 30 i=1,61
!!$         do 30 j=1,61
!!$!            chart(j,i) = blank
!!$ 30   continue
!!$
!!$      do 40 i=1,61,10
!!$         do 40 j=1,61,1
!!$!            chart(j,i) = pt
!!$ 40   continue
!!$
!!$      do 50 k=1,mm
!!$         i = (k - 1) * rs + 1.1
!!$         if (i .gt. 61) i = 61
!!$         do 50 j=1,61,2
!!$            chart(i,j) = pt
!!$ 50   continue
!!$
!!$c *** place data on graph
!!$
!!$      xlen = (xmax - xmin)/60.
!!$      ylen = (ymax - ymin)/60.
!!$
!!$      do 70 i=1,n
!!$         rhl1 = x(i)
!!$         if (rhl1 .gt. xmax) rhl1 = xmax
!!$         if (rhl1 .lt. xmin) rhl1 = xmin
!!$         jx = int((rhl1 - xmin)/xlen)
!!$         if (jx .gt. 60) jx = 60
!!$         point = plus
!!$         if (y(i) .le. 0.) point = aminus
!!$         rhl = y(i)
!!$         if (ne .lt. 0) go to 60
!!$         y(i) = abs(y(i))
!!$         if (y(i) .eq. 0.) y(i) = 10.**ymin
!!$         rhl = log10( y(i) )
!!$ 60      continue
!!$         if (rhl .lt. ymin) rhl = ymin
!!$         if (rhl .gt. ymax) rhl = ymax
!!$         iy = int((rhl - ymin)/ylen)
!!$         if (iy .gt. 60) iy = 60
!!$         chart(61-iy,1+jx) = point
!!$         if (ne .le. 0) go to 70
!!$         if (ye(i) .le. 0.0) ye(i) = 10.**ymin
!!$         rhl = log10( ye(i) )
!!$         if (rhl .lt. ymin) rhl = ymin
!!$         if (rhl .gt. ymax) rhl = ymax
!!$         iy = int((rhl - ymin)/ylen)
!!$         if (iy .gt. 60) iy = 60
!!$         point = chart(61-iy,1+jx)
!!$         if ((point .eq. plus) .or. (point .eq. aminus)) go to 70
!!$         chart(61-iy,1+jx) = ex
!!$ 70   continue
!!$
!!$c *** print graph
!!$
!!$      if (ne .ge. 0) then
!!$         write (6,115)
!!$         write (6,110) plot(1), (nplot(ii),ii=1,22)
!!$      endif
!!$      if (ne .lt. 0) then
!!$         write(kout,115)
!!$         write (6,120) plot(1), (nplot(ii),ii=1,22)
!!$      endif
!!$      xlen = xlen * 100./10.
!!$      ylen = ylen *  60./10.
!!$      ymin = ymax
!!$      if (ne .ge. 0) ymin = 10.**ymax
!!$
!!$      do 90 i=1,55,6
!!$         write (6,130) ymin, (chart(i,j),j=1,61,1)
!!$         if (i .eq. 1) write (6,150) (plot(ii),ii=2,7)
!!$         ymax = ymax - ylen
!!$         if (ne .ge. 0) ymin = 10.**ymax
!!$         ip1 = i + 1
!!$         ip5 = i + 5
!!$
!!$         do 80 i1=ip1,ip5
!!$            write (6,140) (chart(i1,j),j=1,61,1)
!!$ 80      continue
!!$ 90   continue
!!$
!!$      write (6,130) ymin, (chart(61,j),j=1,61)
!!$
!!$      do 100 i=1,7
!!$         xf(i) = xmin
!!$         xmin = xmin + xlen
!!$ 100  continue
!!$
!!$      return
!!$
!!$ 110  format (1h+,7hsig(th),f6.1,2i4,4x,i1,i2,8i1,10i3)
!!$ 115  format (1h1)
!!$ 120  format (1h+,6h polar,f6.1,2i4,4x,i1,i2,8i1,10i3)
!!$
!!$c *** 1h_ inserted in next line for triumf printers (1x removed after
!!$c *** e8.2)
!!$c *** same insertion (deletion) necessary for the ridge.
!!$
!!$ 130  format (1h ,e8.2,101a1)
!!$ 140  format (9x,101a1)
!!$ 150  format (1h+,20x,6f7.3)
!!$      end





!!$
!!$      subroutine totsig( k0, scoul, sigl, xgam, lmax, nspin, nsex )
!!$************************************************************************
!!$      use params
!!$      use switch
!!$      use inputs
!!$      use tcomn
!!$
!!$      implicit real*8 ( a-h , o-z )
!!$
!!$!      real*8   k0 , mp , mn
!!$      
!!$
!!$c      dimension tbr(100,12), tbi(100,12), tr(100,12),ti(100,12)
!!$c      dimension scoul(30)
!!$c      dimension nifty(20)
!!$c      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
!!$c      common /switch/   nifty
!!$c      common /inputs/   tlab , b , ymin1 , ymin2 ,
!!$c     $                  kode , lxmax , nang , ngp , nr
!!$c      common /tcomn/    tr , ti , tbr , tbi
!!$
!!$************************************************************************
!!$
!!$      if ( nifty(6) .eq. 8 ) go to 542
!!$      sigel = 0.
!!$      sigtot = 0.
!!$      brnel = 0.
!!$      brntot = 0.
!!$
!!$c *** total cross section summations
!!$c *** remove coulomb phase from  trhs
!!$c *** n.b. the t*s still contain the short range part of the coulomb
!!$c *** integral so they are not pure nuclear
!!$
!!$!      if (nifty(6) .eq. 8) go to 555
!!$       if (nifty(6) .ne. 8) then
!!$
!!$ 542  do 560 ldum = 1,lmax
!!$         if ((nifty(10) .eq. 0) .or. (nifty(10) .eq. 2)) go to 538
!!$         sigl = scoul(ldum)
!!$         rhl5 = cos(2 * sigl)
!!$         rhl6 = sin(2 * sigl)
!!$
!!$         do 530 nspin = 1,2
!!$            trc = rhl5 * tr(ldum,nspin) + rhl6 * ti(ldum,nspin)
!!$            tic = rhl5 * ti(ldum,nspin) - rhl6 * tr(ldum,nspin)
!!$            tr(ldum,nspin) = trc
!!$            ti(ldum,nspin) = tic
!!$ 530     continue
!!$            endif
!!$!538      continue
!!$ 538     sigtot=sigtot + ldum * ti(ldum,1) + (ldum - 1.) * ti(ldum,2)
!!$         brntot=brntot+ldum * tbi(ldum,1) +
!!$     $                     (ldum - 1.) * tbi(ldum,2)
!!$         sigel=sigel+ldum * (tr(ldum,1)**2 + ti(ldum,1)**2) +
!!$     $                   (ldum - 1.) * (tr(ldum,2)**2 + ti(ldum,2)**2)
!!$         brnel=brnel+ldum * (tbr(ldum,1)**2 + tbi(ldum,1)**2) +
!!$     $                  (ldum - 1.) * (tbr(ldum,2)**2 + tbi(ldum,2)**2)
!!$! 541     if (nsex .ne. 2) go to 560
!!$       endif
!!$! 555      if (nsex .eq. 2) then
!!$         if (nsex .eq. 2) then
!!$c *** 2nd pass on charge exch,redefine t
!!$
!!$         l = ldum
!!$         np2m = 2
!!$         if (nifty(6) .eq. 8) np2m = 6
!!$         
!!$         
!!$         do 550 n = 1,np2m
!!$            np2 = n + np2m
!!$            tbr(l,n) = tbr(l,np2)
!!$            tbi(l,n) = tbi(l,np2)
!!$            tr(l,n) = tr(l,np2)
!!$            ti(l,n) = ti(l,np2)
!!$ 550     continue
!!$            endif
!!$         endif
!!$! 560  continue
!!$
!!$
!!$      if (nifty(6) .eq. 8) then
!!$         return
!!$      endif
!!$      rhl = 10. * 4. * pi * (hbarc/k0)**2
!!$      sigtot = rhl * sigtot
!!$      sigel = sigel * rhl
!!$      brnel = brnel * rhl
!!$      brntot = brntot * rhl
!!$      rhl = sigtot - sigel
!!$      write (6,1130) sigtot, sigel, rhl, tlab, xgam, brntot, brnel
!!$
!!$      return
!!$
!!$c *** formats ***
!!$
!!$ 1130 format (1h0,22htotal cross-section  =,f10.3,9x,8helast , ,
!!$     $20hinel cross-section =,2f10.3//19h cross-sections  in,22h milliba
!!$     $rns  ---------,25h energy(kinetic),gamcoul=,f10.1,e11.3/20h born a
!!$     $pprox total =,f10.3,28h elastic              sigma=,f10.3)
!!$
!!$      end


      subroutine ffacthe (q2, ff, nff)
************************************************************************
      implicit real*8 (a-h, o-z)
      real*8 mp, mn

      dimension ff(4)
      dimension nifty(20)

      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
      common /switch/   nifty
      common /sizes/    achp, acmp, wsp, achn, acmn, wsn

c***   achp for core and achn for valence neutrons
c >>>  first executable statement
       h2 = hbarc*hbarc
c ***  core
       rhl = q2*achp*achp/4./h2
       ff(1) = 0.
       if (rhl.lt.150.0)ff(1)=exp(-rhl)
c ***  valence neutrons
       rhl = q2*achn*achn/4./h2
       rhl1 = q2*achn*achn/6./h2
       ff(2) = 0.
       if(rhl.lt.150.0)ff(2)=(1.-rhl1)*exp(-rhl)
       return
       end
c
      subroutine ffactbe11(q2,formf)
c*********************************************************************
c **  calculate the transformed density rho(q) for berilium
      implicit real*8(a-h,o-z)
      dimension formf(4)
      dimension ffbe10(4)

      call ffactbe10(q2,ffbe10)
      call halodnbe(q2,ffhalo,ffccm)
c *** core
      formf(1) = ffbe10(1)*ffccm
      formf(2) = ffbe10(2)*ffccm
      formf(3) = ffbe10(3)*ffccm
c *** halo
      formf(4) = ffhalo
      return
      end
      
      subroutine ffactli11(q2,formf)
c*********************************************************************
c **  calculate the transformed density rho(q) for lithium
      implicit real*8(a-h,o-z)
      dimension formf(4)
      dimension ffli9(4)

      call ffactli9(q2,ffli9)
      call halodn11(q2,ffhalo,ffccm)
c *** core
      formf(1) = ffli9(1)*ffccm
      formf(2) = ffli9(2)*ffccm
      formf(3) = ffli9(3)*ffccm
c *** halo
      formf(4) = ffhalo
      return
      end



      subroutine ffactbe7(q2,formf)
c*********************************************************************
c **  calculate the  density rho(q) for berilium 7
      use params
      use sizes
      use bhalodn
      use brms

      implicit real*8(a-h,k,m,o-z)
      integer :: kout=16
      dimension formf(4),ff(4)
      
!      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
!      common /sizes/    achp, acmp, wsp, achn, acmn, wsn
!      common/bhalodn/itydn9,itydn11,itycmdn
!      common/brms/rms9,rms11

      h2 = hbarc*hbarc
      q = q2/hbarc/hbarc

       do 20 i=1,4
       formf(i) = 0.d0
 20    continue

       call ffactb8(q2,ff)
       formf(3)=ff(3)
       formf(2)=ff(2)
       formf(4)=0. 
      if (q.lt.001)write(kout,*)'rhoq',formf(2)+formf(3)


c       a9 = rms9/sqrt(3.)
c       rhl = q2*a9*a9/hbarc/hbarc/2.
c       if (rhl.lt.150.0) then
c       formf(1) = 3.*  exp(-rhl) 
c       formf(2) = 4.*  formf(1)
c       formf(3) = 0.
c       formf(4) = 0. 
c       endif


      return
      end

      subroutine martin(kp,k,k0,cthnuc,xmn,tapb,te)
         use params
         use switch
         use blocks
c ****************************************************
c
       implicit real*8 (a-h, i, k, m, o-z)
!       integer option2

!       dimension   nifty(20)
       dimension a(4),b(4),c(4), tapb(4), te(4)

!       common /params/   hbarc, pi, mp, xmn, nz, na, nes, nwaves
!       common /switch/   nifty
!       common/block2/option2

c **   qqq is the square of the transfer momentum
       if (option2.eq.1)then
       qqq = 2.*k0*k0*(1.-cthnuc)/hbarc/hbarc
       else
       qqq = (kp**2+k**2-2.*kp*k*cthnuc)/hbarc/hbarc
       end if
       if(qqq.lt.0)qqq=-qqq
c **   calculations a la martin
       sthnuc0 = 1. - cthnuc*cthnuc
       sthnuc0 = 2.*sqrt(sthnuc0)
c**    parametrizations for elab=156mev
       a(1)=0.47570
       a(2)=0.40160
       a(3)=0.10801
       a(4)=0.41842
       b(1)=-0.35493
       b(2)=-0.33218
       b(3)=0.87967
       b(4)=0.21862
       c(1)=0.37200
       c(2)=0.32240
       c(3)=1.19499
       c(4)=0.37691
       tapb(1) = a(1)*(1.+b(1)*qqq)*exp(-c(1)*qqq)
       tapb(2) = tapb(1)
       tapb(3) = a(2)*(1.+b(2)*qqq)*exp(-c(2)*qqq)
       tapb(4) = tapb(3)
       te(1) = a(3)*(1.+b(3)*qqq)*exp(-c(3)*qqq)
       te(2) = te(1)
       te(3) = a(4)*(1.+b(4)*qqq)*exp(-c(4)*qqq)
       te(4) = te(3)
c **   change from martin to landau conventions
c **   change to calculations a la martin
       do 300 j=1,4
       tapb(j)=-tapb(j)/2./pi/pi/hbarc/mp
       te(j)=-te(j)/2./pi/pi/hbarc/mp*sthnuc0
       te(j) = te(j)*na*kp*k/(na+1.)/k0/k0
c      te(j) = te(j)*0.
 300    continue
c **   reduce amplitudes by a common factor as a check to
c **   the born approximation
c      do 350 j=1,4
c      tapb(j) = tapb(j)*0.01
c      te(j) = te(j)*0.01
c350    continue
       return
       end



      subroutine xsrecoil()
c==============================================================
      use params
      implicit real*8  (a-h, k, m, o-z)
      integer :: kout=16
      complex*16 zz,aux
      real*8, allocatable :: theta(:)
      dimension   formf(4) !theta(181)
!      common /params/ hbarc, pi, mp, mn, nz, na, nes, nwaves
      ffccm = 0.
      zz = cmplx(0.d0,1.d0)
      h2 = hbarc*hbarc
      rewind 77
      read(77,*)noangs, kcore
      allocate(theta(noangs))

      do 50 i=1,50
      q = (i-1)*0.1
      q2 = q*q*hbarc*hbarc
      call halodnbe(q2,ffhalo,ffccm)
      write(kout,*)'rhocm',q, ffhalo,ffccm
 50   continue

      do 100 i=1,noangs
      read(77,*) theta(i),tncre,tncimag
      write(kout,*)'tna', tncre,tncimag
      x = cos(theta(i)*pi/180.)
      q2 = (2*kcore*kcore*(1-x))
      q2 = q2*h2
      if (na.eq.10) then 
      call halodnbe(q2,ffhalo,ffccm)
      endif
      aux = (tncre + zz*tncimag)*ffccm 
      xssh = aux*conjg(aux)
      xssh = xssh*10.
      write(78,*) theta(i), xssh
 100  continue
      return
      end







c===================================================================
c *** calculates gaussian form factor for harmonic oscillator nuclei
c===================================================================
      subroutine ffhmo( q2 , nff , ff1 , ff2 )
      use params
      use switch
      use sizes
      use inputs
     

      implicit real*8 (a-h, o-z)

!      integer :: nz,nn
!               real*8 mp, mn
!      dimension nifty(20)
!      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
!      common /switch/   nifty
!      common /sizes/    achp, acmp, wsp, achn, acmn, wsn

c >>> first executable statement <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      h2 = hbarc * hbarc
      
      achp=ahop
      acmp=ahop

      achn=ahon
      acmn=ahon 

      nz=nzclus
      nn=nnclus
 
!      write(*,*) 'nff=',nff; stop
!      write(*,*)'achp,achn',achp,achn     
!      write(*,*)'nz,nn=',nz,nn ; stop

c *** protons

      rhl = q2 * achp * achp/4./h2
      rhl1 = acmp * acmp * q2/(6.0 * nz)/h2
      ff1 = 0.
      if (rhl .lt. 150.0) ff1 = (1.-(nz-2.)*rhl1)*exp(-rhl)
      write(91,*) 'q2,ff1=',q2,ff1

c *** divide out the proton form factor

! commented by AMoro
!      if (nifty(15) .eq. 0) ff1 = ff1 * (1.0 + q2/18.2/h2)**2
      if (nff .eq. 1) return

c *** neutron matter form factor

      rhl = q2 * achn * achn/4./h2
!      rhl1 = acmn * acmn * q2/(6.0 * (na - nz))/h2
      rhl1 = acmn * acmn * q2/(6.0 * nn)/h2
      ff2 = 0.
!      if (rhl .lt. 150.0) ff2 = (1.-(na-nz-2.)*rhl1)*exp(-rhl)
       if (rhl .lt. 150.0) ff2 = (1.-(nn-2.)*rhl1)*exp(-rhl)

c *** divide out the neutron form factor

!      if (nifty(15) .eq. 0) ff2 = ff2 * (1.0 + q2/18.2/h2)**2
      return
      end



c
c
      subroutine voptbe( vll, k, kp, neles, k0 )
************************************************************************


c *** calculates the optical potential for nucleons scattering from a
c     Lithium  nucleus and does a full t*rho projection

      use switch
      use inputs
      use wr
      use halo
      use params

      implicit real*8  (a-h, k, m, o-z)
      integer  code, xi, maxls
      real*8   imucl, imusl, imvc, imvs
      real*8   imu1,imu2
      real*8   imucl2
!      complex*16  ucentl( 0:99), uspinl( 0:99 ), zi
      complex*16, allocatable::  ucentl(:), uspinl(:)
      complex*16  ucl, usl, uc, us, zi

!      dimension   vll( 0:99, 2, 6 )
      real*8 :: vll(0:neles,2,6)
      dimension   xis( 32 ), wt( 32 )
!      dimension   pl( 0:99 ), plp( 0:99 )
      real*8, allocatable ::  pl(:), plp(:)
      dimension   tapb(4), tamb(4), tcpd(4), tcmd(4), te(4)
      dimension   tapbm(4), tem(4)
      dimension   tapb2(4),te2(4)
!      dimension   nifty(20)
      dimension   formf(4)

!!$      common   /params/  hbarc, pi, mp, mn, nz, na, nes, nwaves
!!$      common   /switch/ nifty
!!$      common  /inputs/tlab
!!$      common  /wr/iwrite,iread,istrong
!!$      common/halo/ioptls,ioptcnt

      allocate(ucentl(0:neles))
      allocate(uspinl(0:neles))
      allocate(pl(0:neles))
      allocate(plp(0:neles))

      data  maxls / 100 /
c       nlicav = Be10 alpha core
c       nlicn =  Be10 4 neutron valence
c       nlicp =  Be10 2 proton valence
c       nlivls = Be11 halo (spin orbit)
c       nlivc =  Be11 halo (central)
  
        nlicav = 1
        nlicn = 1
        nlicp = 1
        nliv = 1
        nlivls = 1
        nlivc = 1
        if (ioptls.eq.2) nlivls = 0.
        if (ioptcnt.eq.1) nlivc = 0.
        ncos = 32
        nn = na - nz
        amass = na * mn
        aovera = (na -1.)/(na*1.)
c ***   set ustrong=0
        if (istrong.eq.1)aovera=0.
         twopi2 = 2. * pi * pi
         if ( nifty(17) .eq. 1 ) aovera = 1.
         a = -1.
         b = 1.
         code = 11
         call gauss2( ncos, code, a, b, xis, wt )

      zi = ( 0. , 1.)

      do 10 l = 0, maxls - 1
         ucentl(l) = ( 0., 0. )
         uspinl(l) = ( 0., 0. )
  10  continue

      do 20 xi = 1, ncos
         x = xis( xi )
         cthnuc = x
         theta = acos( cthnuc ) * 180./pi
!         call amp2(kp, k, k0, cthnuc, mn, tapb, te)
         call amp2(kp, k, k0, cthnuc, tapb, te)
         q2 = kp**2 + k**2 - 2.* kp * k * cthnuc
         if ( q2 .lt. 0 ) q2 = -q2
         if (na.eq.10) call ffactbe10( q2, formf)
         if (na.eq.11) call ffactbe11( q2, formf)

         call legpol( x, pl, neles )
         call plprme( x, plp, neles )

         do 30 l = 0, neles
            if ( l .eq. 0 ) then
               uspinl(l) = ( 0., 0. )
            else
               reusl=nlicav*te(3)*formf(1) + nlicn * te(4) * formf(2)
     $               + nlicp * (2*te(3)-te(4)) * formf(3)
     $               +  nlivls * te(4) * formf(4)
               reusl = 1./( (2.*l) * (l + 1.) ) * reusl
               reusl = reusl * sqrt(1.- x**2) * plp(l) * wt(xi)
               imusl=nlicav*te(1)*formf(1) + nlicn * te(2) * formf(2)
     $                +  nlicp *(2*te(1) - te(2)) * formf(3)
     $                + nlivls * te(2) * formf(4)       
               imusl = -1./( (2.*l) * (l + 1.) ) * imusl
               imusl = imusl * sqrt(1.- x**2) * plp(l) * wt(xi)
               uspinl(l) = uspinl(l) + cmplx( reusl, imusl )
            endif
            reucl=nlicav * tapb(1) * formf(1) + nlicn*tapb(2)*formf(2)
     $                 + nlicp*(2*tapb(1) - tapb(2))*formf(3)
     $                 + nlivc*tapb(2)*formf(4)
            reucl = 1./2. * reucl
            reucl = reucl * pl(l) * wt(xi)
            imucl=nlicav * tapb(3) * formf(1) + nlicn*tapb(4)*formf(2)
     $                  + nlicp*(2*tapb(3) - tapb(4))*formf(3)   
     $                  + nlivc*tapb(4)*formf(4)     
            imucl = 1./2. * imucl
            imucl = imucl * pl(l) * wt(xi)
            ucentl(l) = ucentl(l) + cmplx( reucl, imucl )
  30     continue
         if (( kp .eq. k) .and. ( kp .eq. k0 ) ) then
            if (xi .eq. 1 ) then
c              *** set up headers for write of resum ***
               write(6,1000)
               write(6,1001)
               write(6,*)
            endif
            ucl = 0.
            usl = 0.
            do 60 l = 0, neles
               ucl = ucl + ucentl(l) * pl(l) * (2*l+1.)
               usl = usl + uspinl(l) * plp(l) * (2*l+1.)
  60        continue
            revc = nhe4 * tapb(1) * formf(1) + nhev * tapb(2) * formf(2)
            imvc = nhe4 * tapb(3) * formf(1) + nhev * tapb(4) * formf(2)
            revs = nhe4 * te(3) * formf(1) + nhev * te(4) * formf(2)
            imvs = -(nhe4 * te(1) * formf(1) + nhev * te(2) * formf(2))
            uc = cmplx( revc, imvc )
            us = cmplx( revs, imvs )
            write(6,1002) theta, q2, abs(ucl), abs(uc),
     $                               abs(usl), abs(us)
         endif
  20  continue


      do 40 l = 0, neles
         if ( l .eq. 0 ) then
            vll(l,1,1) = twopi2 * ucentl(l)
            vll(l,2,1) = twopi2 * dimag( ucentl(l) )
            vll(l,1,2) = vll(l,1,1)
            vll(l,2,2) = vll(l,2,1)
         else
            vll(l,1,1) = twopi2 * ( ucentl(l) + l * uspinl(l) )
            vll(l,2,1) = twopi2 * dimag( ucentl(l) + l * uspinl(l) )
            vll(l,1,2) = twopi2 * ( ucentl(l) - (l+1.) * uspinl(l) )
            vll(l,2,2) = twopi2 * dimag( ucentl(l) - (l+1.) * uspinl(l))
         endif
         do 50 ir = 1,2
            do 50 nspin = 1, 2
               vll(l, ir, nspin) = aovera * vll( l, ir, nspin )
  50     continue
  40  continue
c ***************************************
      return

c     *** formats ***

 1000 format('1')
 1001 format(' ',' cos  ',9x,'q2',13x,'resum for vcent',11x,'vcent',
     $           11x,'resum for vspin',11x,'vspin')
 1002 format(' ',f5.1,5x,f12.3,9x,e13.6,3(8x,e13.6))
 1003 format(' ',31x,e13.6,3(8x,e13.6))
      end



      subroutine voptboro( vll, k, kp, neles, k0 )
************************************************************************


c *** calculates the optical potential for nucleons scattering from a
c     b8  nucleus and does a full t*rho projection
      use params
      use switch
      use inputs
      use wr
      use halo

      implicit real*8  (a-h, k, m, o-z)
      integer  code, xi, maxls
      real*8   imucl, imusl, imvc, imvs
      real*8   imu1,imu2
      real*8   imucl2
!      complex*16  ucentl( 0:99), uspinl( 0:99 ), zi
      complex*16, allocatable::  ucentl(:), uspinl(:)
      complex*16  ucl, usl, uc, us

!      dimension   vll( 0:99, 2, 6 )
      real*8::   vll( 0:neles, 2, 6 )
      dimension   xis( 32 ), wt( 32 )
!      dimension   pl( 0:99 ), plp( 0:99 )
      real*8, allocatable :: pl(:),plp(:)
      dimension   tapb(4), tamb(4), tcpd(4), tcmd(4), te(4)
      dimension   tapbm(4), tem(4)
      dimension   tapb2(4),te2(4)
!      dimension   nifty(20)
      dimension   formf(4)

!      common   /params/  hbarc, pi, mp, mn, nz, na, nes, nwaves
!      common   /switch/ nifty
!      common  /inputs/tlab
!      common  /wr/iwrite,iread,istrong
!      common/halo/ioptls,ioptcnt

      allocate(ucentl(0:neles))
      allocate(uspinl(0:neles))
      allocate(pl(0:neles))
      allocate(plp(0:neles))

      data  maxls / 100 /

        ncos = 32
        nn = na - nz
        amass = na * mn
        aovera = (na -1.)/(na*1.)
        nval = 1
        if (nz.eq.4.and.na.eq.6) nval = 0.
c         if (k.eq.kp)write(6,*)'in voptboro'     
   
c ***   set ustrong=0
        if (istrong.eq.1)aovera=0.
         twopi2 = 2. * pi * pi
         if ( nifty(17) .eq. 1 ) aovera = 1.
         a = -1.
         b = 1.
         code = 11
         call gauss2( ncos, code, a, b, xis, wt )

      zi = ( 0. , 1.)

      do 10 l = 0, maxls - 1
         ucentl(l) = ( 0., 0. )
         uspinl(l) = ( 0., 0. )
  10  continue

      do 20 xi = 1, ncos
         x = xis( xi )
         cthnuc = x
         theta = acos( cthnuc ) * 180./pi
!        call amp2(kp, k, k0, cthnuc, mn, tapb, te)
         call amp2(kp, k, k0, cthnuc, tapb, te)
         q2 = kp**2 + k**2 - 2.* kp * k * cthnuc
         if ( q2 .lt. 0 ) q2 = -q2         
         call ffactb8( q2, formf)
         call legpol( x, pl, neles )
         call plprme( x, plp, neles )

         do 30 l = 0, neles
            if ( l .eq. 0 ) then
               uspinl(l) = ( 0., 0. )
            else
               reusl= te(3)*formf(1) +  te(4) * formf(2)
     $                +  (2*te(3)-te(4)) * formf(3)
     $                +  nval*  (2*te(3)-te(4)) * formf(4)
               reusl = 1./( (2.*l) * (l + 1.) ) * reusl
               reusl = reusl * sqrt(1.- x**2) * plp(l) * wt(xi)
               imusl=  te(1)*formf(1) +  te(2) * formf(2)
     $                +    (2*te(1) - te(2)) * formf(3)
     $                +   nval*  (2*te(1) - te(2))    * formf(4)       
               imusl = -1./( (2.*l) * (l + 1.) ) * imusl
               imusl = imusl * sqrt(1.- x**2) * plp(l) * wt(xi)
               uspinl(l) = uspinl(l) + cmplx( reusl, imusl )
            endif
            reucl=  tapb(1) * formf(1) +  tapb(2)*formf(2)
     $                 +  (2*tapb(1) - tapb(2))*formf(3)
     $                 +  nval* (2*tapb(1) - tapb(2)) *formf(4)
            reucl = 1./2. * reucl
            reucl = reucl * pl(l) * wt(xi)
            imucl=  tapb(3) * formf(1) + tapb(4)*formf(2)
     $                  + (2*tapb(3) - tapb(4))*formf(3)   
     $                  + nval*  (2*tapb(3) - tapb(4)) *formf(4)     
            imucl = 1./2. * imucl
            imucl = imucl * pl(l) * wt(xi)
            ucentl(l) = ucentl(l) + cmplx( reucl, imucl )
  30     continue
         if (( kp .eq. k) .and. ( kp .eq. k0 ) ) then
            if (xi .eq. 1 ) then
c              *** set up headers for write of resum ***
               write(6,1000)
               write(6,1001)
               write(6,*)
            endif

         endif
  20  continue


      do 40 l = 0, neles
         if ( l .eq. 0 ) then
            vll(l,1,1) = twopi2 * ucentl(l)
            vll(l,2,1) = twopi2 * dimag( ucentl(l) )
            vll(l,1,2) = vll(l,1,1)
            vll(l,2,2) = vll(l,2,1)
         else
            vll(l,1,1) = twopi2 * ( ucentl(l) + l * uspinl(l) )
            vll(l,2,1) = twopi2 * dimag( ucentl(l) + l * uspinl(l) )
            vll(l,1,2) = twopi2 * ( ucentl(l) - (l+1.) * uspinl(l) )
            vll(l,2,2) = twopi2 * dimag( ucentl(l) - (l+1.) * uspinl(l))
         endif
         do 50 ir = 1,2
            do 50 nspin = 1, 2
               vll(l, ir, nspin) = aovera * vll( l, ir, nspin )
  50     continue
  40  continue
c ***************************************
      return

c     *** formats ***

 1000 format('1')
 1001 format(' ',' cos  ',9x,'q2',13x,'resum for vcent',11x,'vcent',
     $           11x,'resum for vspin',11x,'vspin')
 1002 format(' ',f5.1,5x,f12.3,9x,e13.6,3(8x,e13.6))
 1003 format(' ',31x,e13.6,3(8x,e13.6))
      end


