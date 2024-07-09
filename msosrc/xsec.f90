      subroutine totsig( k0, sigl, xgam, lmax, nspin, nsex )
c ***************************************************************
c *** calculation of r-matrix and t-matrix by matrix inversion with
c     potential in momentum space for spin 0 and 1/2, at any energy.
c     
c     This subroutine solves the LS equation in momentum space

c ** Original program: lpotp nucleon elasting scattering from spin 0:
c     M. J. Paez, M. E. Sagen and R. H. Landau, 
c     Comput. Phys. Commun. 52, 141 (1988)  

c ** LPOTPS is an extensively revised and extended version of the 
c    program LPOTP     
c    R. Crespo, R. C. Johnson, and J. A. Tostevin
c    Phys. Rev. C 41, 2257 (1990)

c ** f90 version by R. Crespo and A.M. Moro  (spring 2002)
c ***************************************************************
      use params
      use inputs
      use tcomn
      use bscoul
      implicit real*8 ( a-h , o-z )
      real*8   k0 
!       if ( nifty(6) .ne. 8 ) then
            sigel = 0.
            sigtot = 0.
!       end if
c *** total cross section summations
c *** remove coulomb phase from  trhs
c *** n.b. the t*s still contain the short range part of the coulomb
c *** integral so they are not pure nuclear

!      if (nifty(6) .eq. 8) go to 541

 542  do 560 ldum = 0,lmax
         if ((coulomb .eq. 0) .or. (coulomb .eq. 2)) go to 540
         sigl = scoul(ldum)
         rhl5 = cos(2 * sigl)
         rhl6 = sin(2 * sigl)

         do 530 nspin = 1,2
            trc = rhl5 * tr(ldum,nspin) + rhl6 * ti(ldum,nspin)
            tic = rhl5 * ti(ldum,nspin) - rhl6 * tr(ldum,nspin)
            tr(ldum,nspin) = trc
            ti(ldum,nspin) = tic
 530     continue

 540     continue
         sigtot =sigtot+(ldum + 1) * ti(ldum,1) + ldum * ti(ldum,2)
         sigel = sigel + (ldum + 1)*(tr(ldum,1)**2 + ti(ldum,1)**2)+
     $                   ldum * (tr(ldum,2)**2 + ti(ldum,2)**2)
 541     if (nsex .ne. 2) go to 560

c *** 2nd pass on charge exch,redefine t

         l = ldum
         np2m = 2
!         if (nifty(6) .eq. 8) np2m = 6

         do 550 n = 1,np2m
            np2 = n + np2m
            tr(l,n) = tr(l,np2)
            ti(l,n) = ti(l,np2)
 550     continue
 560  continue


!      if (nifty(6) .eq. 8)    return
      rhl = 10. * 4. * pi * (hbarc/k0)**2
      sigtot = rhl * sigtot
      sigel = sigel * rhl
      rhl = sigtot - sigel
      write (kout,1130) sigtot, sigel, rhl, tlab, xgam

      return

c *** formats ***

 1130 format (1h0,22htotal cross-section  =,f10.3,9x,8helast , ,
     $20hinel cross-section =,2f10.3//19h cross-sections  in,22h milliba
     $rns  ---------,25h energy(kinetic),gamcoul=,f10.1,e11.3)

      end

c--------------------------------------------------------------------
      subroutine xsects(alph, xgam, sig0, ldmax,nspin, nsex)
*********************************************************************
      use kinemt
      use params
      use switch
      use sizes
      use ranges
      use inputs
      use rcomn
      use spins
      use tcomn
      use bscoul
      use tatheta !AMoro
 
      implicit real*8 ( a-h, m , o-z )
      integer :: ldmax,nspin,nsex
      real*8::ff(4)
      real*8:: k02, kone, ktwo,x
      real*8:: alph,xgam,sig0
      real*8,allocatable:: pl(:), plp(:)
      real*8, allocatable,dimension(:):: dsig,sigbrn,
     $        dflip,dnof,polar,spinrot
      complex*16:: zi , zcul , zifase
      complex*16:: za , zb , zc , zd , ze
      complex*16:: zan , zbn , zcn , zdn , zen
      
      dimension sefz(15)


c-----calculate angular distribution -----------------------------------

      zi = ( 0. , 1. )
      k02 = k0 * k0
      noangs = (180./nang) + 1
      cosold = 2.
    
       write (77,*)noangs,k0/hbarc
       write(99,*)' + Allocating',noangs,'angles for theta,dnof,...'
       allocate(polar(noangs))
       allocate(spinrot(noangs))
       allocate(sigbrn(noangs))
       allocate(dflip(noangs))
       allocate(dsig(noangs))
       allocate(dnof(noangs))
       allocate(theta(noangs)) 
       allocate(tnct(noangs)) 

       allocate(pl(0:ldmax))
       allocate(plp(0:ldmax))
      
       do 520 iang = 1,noangs
         theta(iang) = (iang - 1) * nang
         x = cos(theta(iang) * pi/180.)

c *** rhl chanmge of paez
         call legpol (x, pl, ldmax)

         gbr = 0.e00
         gbi = 0.e00
         gr =  0.e00
         gi =  0.e00
!         if (nifty(6) .eq. 3) then
            call plprme (x, plp, ldmax)
!         endif
         sumbr = 0.
         sumbi = 0.
         sumr = 0.
         sumi = 0.
         
         do 450 ldum = 0, ldmax
            sag = ldum + 1
            sumr = sumr   + (sag * tr(ldum,1)  +
     $             ldum * tr(ldum,2))  * pl(ldum)
            sumi = sumi   + (sag * ti(ldum,1)  +
     $             ldum * ti(ldum,2))  * pl(ldum)
!               if (nifty(6) .eq. 3) then
            gr  = gr  + (tr(ldum,1)  - tr(ldum,2))  * plp(ldum)
            gi  = gi  + (ti(ldum,1)  - ti(ldum,2))  * plp(ldum)
!               endif
!            write(15,*) ldum,tr(ldum,1),tr(ldum,2)
 450        continue
            rhl = (sumr * hbarc/k0)
            rhl5 = (sumi * hbarc/k0)
            rhl6 = 0.
            rhl7 = (hbarc/k0)**2
            dflip(iang) = rhl7 * (1. - x * x) * (gr**2 + gi**2)
            dnof(iang) = rhl7 * (sumr**2 + sumi**2)
            dsig(iang) = dflip(iang) + dnof(iang)
            dsigma = dsig(iang)

!           write(kout,*)'dflip=',dflip(iang)
!           write(kout,*)'dnof=',dnof(iang)
c           *** polarization

            if (dsig(iang) .ne. 0) then
               polar(iang) = 2. * rhl7 * sqrt(1. - x * x) *
     $                        (sumi*gr-sumr*gi)/dsig(iang)
               spinrot(iang) = -2.*rhl7*sqrt(1.-x*x) *
     $                        (sumr*gr+sumi*gi)/dsig(iang)
               r = (x * (dnof(iang) - dflip(iang)) - 2. * (1. - x * x) *
     $                  rhl7 * (sumr * gr + sumi * gi))/dsig(iang)
            endif


 453     coscm = x
         rhl7 = (hbarc)**2
         tivt = +2. * (k0**2) * (1. - coscm)
         qtf = 0.
         q2f = tivt/rhl7
         if (q2f .ge. 0.) qtf = sqrt(q2f)
c        if (coulomb .eq. 0) go to 510

c *** include coulomb ampl + phase changes

         ff(1) = 1
         if (coulomb .ge. 3) go to 500
c ***     reset raddii in case equiv wiped them out
         if (na .eq. 2) go to 480
!         if ((na .le. 16) .and. (na .gt. 4)) go to 490
!         if (nifty(16) .eq. 0) go to 490
         goto 490

         go to 500
 480     qr = qtf * rcoul
         if (qr .eq. 0.) go to 500
         ff(1) = (3./qr/qr) * (sin(qr)/qr - cos(qr))
         go to 500
 490     ff(1) = (1. - alph * (qtf * acmp/1.)**2) *
     $            exp(-.25 * (qtf * achp/1.)**2)
c *** important: hbarc changed to 1. error in lpotp
c         write(kout,*)'ff(1)=',ff(1)
c *** important: ff(1)=1 for coulomb=1. error in lpotp
 500      if (coulomb.eq.1.or.coulomb.ge.3)ff(1)=1.

c *** the value of theta(iang) from e-10 to e-06 necessary for vax

         if (theta(iang) .eq. 0) theta(iang) = 1.e-06
         phasec = sin(theta(iang) * pi/360.)
         fcoul = (-xgam * ff(1)/2./k0/phasec/phasec) * hbarc
         dsigcl = fcoul * fcoul * 10.
         phasec = -2. * xgam * log(phasec)
         rhl6 = phasec + 2 * sig0
         phasec = phasec + 2. * sig0
         zifase = zi * phasec
         zcul = fcoul * (cos(phasec) + zi * sin(phasec))

c *** regular cases of either 0 x 0 or 0 x 1/2

 501     rhl = (sumr * hbarc/k0) + fcoul * cos(phasec)
         rhl5 = (sumi * hbarc/k0) + fcoul * sin(phasec)
         dsigma = rhl * rhl + rhl5 * rhl5 + dflip(iang)

c ***    now coulomb in polarization
         if (dsigma>1.e-15)then
         polar(iang) = 2.*(hbarc/k0)*sqrt(1.-x*x)*
     *               (rhl5*gr-rhl*gi)/dsigma
         spinrot(iang) = -2.*(hbarc/k0)*sqrt(1.-x*x)*
     *                (rhl*gr+rhl5*gi)/dsigma
         endif
         dsi0ma = dsigma
         rhl = (sumr * hbarc/k0)
         rhl5 = (sumi * hbarc/k0)
         sumbr = sumbr * hbarc/k0 + fcoul * cos(rhl6)
         sumbi = sumbi * hbarc/k0 + fcoul * sin(rhl6)
         rhl6 = fcoul/sqrt(dsig(iang))
         dsigb = sumbr * sumbr + sumbi * sumbi + dflipb
         dnof(iang) = rhl * rhl + rhl5 * rhl5
         dnofb = sumbr**2 + sumbi**2

c *** calc of t.f. to lab from cm, see goldberg+watson
c ***

c510     if (coulomb.eq. 0) then
c           call expera( 1.d0, can1, ldmax, za, zb,
c    $                   zc, zd, ze, tr, ti, sefz)
c        endif

         ektwo = mnuc + plab * plab * mnuc * (1. - coscm)/s
         ktwo = sqrt(ektwo**2 - mnuc2)
         ekone = el + mnuc - ektwo
         kone = sqrt(abs(ekone**2 - mp2))
         coslab = sqrt(abs(kone**2 - k02 * (1. - coscm**2)))/kone
    
         if (coslab .gt. cosold) coslab = -coslab
         cosold = abs(coslab)
         dot = (mnuc2 + mp2 + 2. * ekone * ektwo - s)/2.
         cmtol = kone * kone * s/
     $                  (plab * mnuc *(ektwo * kone - dot * ekone/kone))

c *** convert cross sections from fm**2 to mb

!         if (nifty(6) .eq. 3) dsi0ma = dsigma
         dsi0ma = dsigma !AMoro
         dsi0ma = dsi0ma * 10.
 513     dnof(iang) = dnof(iang) * 10.
         dsigma = dsigma * 10.
         dsig(iang) = dsig(iang) * 10.
         dflip(iang) = dflip(iang) * 10.
         dsigl = dsigma * cmtol

c *** print out of cross sections vs theta

c *** scattering p_A amplitudes
         rc1 = (hbarc/k0)
         rc2 = (hbarc/k0*(1-x*x))
         tnct(iang)=cmplx(sumr*rc1,sumi*rc1)
!         write(77,4050)theta(iang),sumr*rc1,sumi*rc1,
!     $             gr*rc2,gi*rc2
         write(25,1130)theta(iang),dsigma,polar(iang),spinrot(iang)
!         write(25,1130)theta(iang),dsigma,abs(fcoul)**2
         
 1986    dsig(iang) = dsigma

 520  continue

c *** calculate total cross-sections ***
      call totsig( k0, sigl, xgam, ldmax, nspin, nsex )
     
      deallocate(pl,plp)
      deallocate(spinrot,sigbrn,dflip,dsig,dnof,polar)

      return

c *** formats ***

 1093 format(1h ,f7.0,2f8.3,3e11.3,f8.3,2e11.3,3f9.3,e11.3)
 1100 format(1h ,f7.0,f8.3,2e11.3,2f8.1,f10.3,2e11.3,2f9.3,e12.4
     1,e10.2)
 1120 format (1h ,18x,10e11.3)
 1130 format (f7.0,3e13.3)
 1921 format(1h ,' sigot= ',e13.5,' sigit= ',e13.5,
     $           ' sig2t= ',e13.5,' e= ',f8.2,' xgam= ',e13.5)
 4000 format(' ', f7.2, 3e11.3)
 4010 format(1h ,f7.2,2e11.3)
 4050 format(f7.2,4e11.3)

      end


