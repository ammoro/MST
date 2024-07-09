
       subroutine amp2(kp,k,k0,cthnuc,tapb,te)
c ***************************************************************
c *** calculation of r-matrix and t-matrix by matrix inversion with
c     potential in momentum space for spin 0 and 1/2, at any energy.

c ** Original program: lpotp nucleon elasting scattering from spin 0:
c     M. J. Paez, M. E. Sagen and R. H. Landau, 
c     Comput. Phys. Commun. 52, 141 (1988)  

c ** LPOTPS is an extensively revised and extended version of the 
c    program LPOTP     
c    R. Crespo, R. C. Johnson, and J. A. Tostevin
c    Phys. Rev. C 41, 2257 (1990)

c ** f90 version by R. Crespo and A.M. Moro  (spring 2002)
c 
c    Main changes
c    - common blocks replaced by modules
c    - loops on partial waves start on l=0 (not l=1) 
c    - vll changed from real*8 to complex*16 (1 index dropped)
c    - namelist input
c ***************************************************************
       use params
       use switch
       use blocks
       use halo
!       use ampl
       use amps
       use amphe
       implicit real*8 (a-h, i, k, m, o-z)
C       integer option2
!       integer ioptls,ioptcnt
!       complex*16 an,ap,cn,cp
       complex*16 taux
       complex*16 cint2d
       dimension tapb(4),te(4)

!      dimension   nifty(20)

!      common /params/   hbarc, pi, mp, xmn, nz, na, nes, nwaves
!      common /switch/   nifty
!      common/block2/option2
!      common/halo/ioptls,ioptcnt
!      common/ampl/q(200),qq(200),ap(200,200),cp(200,200),nkmx1,nqmx1
!      common/amphe/an(200,200),cn(200,200)

c ***  qp and qqp are transfer and total n-n momentum
       if (option2.eq.1)then
       qp = k*k + kp*kp -2.*k*kp*cthnuc
       if(qp.lt.0)qp=-qp
       qqp = ((na+1)*k0/(na*1.))**2 - qp
       if(qqp.lt.0)qqp=-qqp
       qp = sqrt(qp)/hbarc
c ***  no recoil effects
       qqp = sqrt(qqp)/2./hbarc
       else if (option2.eq.3) then
       qp = 0.
       qqp = (na+1)*k0/(2.*na)/hbarc
       else
       qp = k**2 +kp**2 -2.*k*kp*cthnuc
       if(qp.lt.0)qp=-qp
       qqp = (k*k + kp*kp)/2. - qp/4.
       if(qqp.lt.0.)qqp=-qqp
       qp = sqrt(qp)/hbarc
c ***  no recoil effects
       qqp = sqrt(qqp)/2./hbarc
       end if
       taux=cint2d(qq,q,ap,qqp,qp,nkmx1,nqmx1,5,200)
       tapb(1)=real(taux)
       tapb(3)=dimag(taux)
       taux=cint2d(qq,q,an,qqp,qp,nkmx1,nqmx1,5,200)
       tapb(2)=real(taux)
       tapb(4)=dimag(taux)
       taux = cint2d(qq,q,cp,qqp,qp,nkmx1,nqmx1,5,200)
       te(1) = real(taux)
       te(3) = dimag(taux)
       taux = cint2d(qq,q,cn,qqp,qp,nkmx1,nqmx1,5,200)
       te(2)=real(taux)
       te(4)=dimag(taux)
c ** change from scattering amplitude to transtion amplitude
c ** and to landau conventions
       do 300 j=1,4
       tapb(j) = -tapb(j)/2./pi/pi/mp/hbarc
       te(j) = -te(j)/2./pi/pi/mp/hbarc
  300  continue
       if (ioptls.eq.1)then
c      write(*,*)'ioptls=',ioptls
       do 305 j=1,4
       te(j) = 0.
 305   continue
       endif
       return
       end
c


       function rofc(kp,k,k0,cthnuc,qqt)
c*******************************************************
       use params
       use sizes
       implicit real*8 (a-h,o-z)
       real*8 k,kp,k0
       real*8 iss, ipp, j0, j1, i2ss,idd
!       real*8 mp,mn
       complex*16 z, zj0, zj1
!       common/params/ hbarc, pi, mp, mn, nz, na,nes,nwaves
!       common/sizes/  achp, acmp, wsp, achn, acmn, wsn
c ***  statement functions now follow
       zj0(z)=sin(z)/z
       zj1(z)=sin(z)/z/z - cos(z)/z

       z = cmplx(0.,1.)
       h2 = hbarc*hbarc
c***  dealing with computer limitations
       sw1 = k*k/h2 + kp*kp/h2 + 8.*qqt
       sw2 = 270.*k0*k0/h2
       if(sw1.gt.sw2)then
       rofc = 0.d0
c ***
       else
       sn10 = 4.*achp**3/sqrt(pi)
       sn11 = 8.*achp**5/3./sqrt(pi)
       sn12 = 16.*achp**7/15./sqrt(pi)
       sn20 = 8.*achp**3/3./sqrt(pi)
       a = 9./4.
       b = achp**4
       c = 3.*achp*achp
       oscp = achp*achp/2.
       u = 2.*oscp*2.*qqt*sqrt((k*k +kp*kp +2.*k*kp*cthnuc)/h2)
       u1 = sqrt((k*k + kp*kp + 2.*k*kp*cthnuc)/h2)
       u2 = (k*k + k*kp*cthnuc)/h2
       u3 = (kp*kp + k*kp*cthnuc)/h2
       qqt2 = 2.*qqt*hbarc
       j0 = real(zj0(z*u))
       j1 = real(z*zj1(z*u))
       exp1 = -oscp*(k*k + kp*kp + 2.*qqt2*qqt2)/h2
       exp1 = exp(exp1)
       iss = 4.*pi*sn10*exp1*j0
       ipp = 2.*oscp*(k*kp*cthnuc+qqt2*qqt2)/h2*j0 + u*j1
       ipp = sn11*2.*pi*exp1/oscp*ipp

       if (na.eq.16)then
       rofc  = (iss + 3.*ipp)/4./pi
       rofc  = 8.*2.*rofc
       else if (na.eq.40)then
        coef4 = k*kp*cthnuc + qqt2*qqt2
        coef4 = coef4/h2
       i2ss = a *j0
     *        - c*qqt2*u1/hbarc*j1
     *        - c*(k*k+kp*kp+2.*qqt2*qqt2)/2./h2*j0
       i2ss = i2ss * sn20*4.*pi*exp1
       idd = (3.*qqt2*u1/hbarc*(2.*coef4+1./oscp))*j1
     *      +(3.*coef4*coef4
     *       +3.*qqt2*qqt2*u1*u1/h2)*j0
       idd = idd * sn12*2.*pi*exp1
       rofc  = (iss + 3.*ipp + 5.*idd + i2ss)/4./pi
       rofc  = 2.*8.*rofc
       else
       stop
       end if
       end if
       return
       end

c


       function rofls(kp,k,k0,cthnuc,qqt)
c*******************************************************
       use params
       use sizes 
       implicit real*8 (a-h,o-z)
       real*8 k,kp,k0
       real*8 iss, ipp, j0, j1, i2ss,idd
!       real*8 mp,mn
       complex*16 z, zj0, zj1
!       common/params/ hbarc, pi, mp, mn, nz, na,nes,nwaves
!       common/sizes/  achp, acmp, wsp, achn, acmn, wsn
c ***  statement functions now follow
       zj0(z)=sin(z)/z
       zj1(z)=sin(z)/z/z - cos(z)/z

       z = cmplx(0.,1.)
       h2 = hbarc*hbarc
c ***  dealing with computer limitations
       sw1 = k*k/h2 + kp*kp/h2 +8.*qqt
       sw2 = 270.*k0*k0/h2
       if(sw1.gt.sw2)then
       rofls = 0.d0
c ***
       else
       sn10 = 4.*achp**3/sqrt(pi)
       sn11 = 8.*achp**5/3./sqrt(pi)
       sn12 = 16.*achp**7/15./sqrt(pi)
       sn20 = 8.*achp**3/3./sqrt(pi)
       a = 9./4.
       c = 3.*achp*achp
       oscp = achp*achp/2.
       u = 2.*oscp*2.*qqt*sqrt((k*k +kp*kp +2.*k*kp*cthnuc)/h2)
       u1= sqrt((k*k +kp*kp +2.*k*kp*cthnuc)/h2)
c      qqt2 is 2.*q^ = x in r.c. notes
       qqt2 = 2.*qqt*hbarc
       coef =(k*kp*cthnuc+qqt2*qqt2)/h2
       j0 = real(zj0(z*u))
       j1 = real(z*zj1(z*u))
       exp1 = -oscp*(k*k + kp*kp + 2.*qqt2*qqt2)/h2
       exp1 = exp(exp1)
       iss = -4.*pi*sn10*exp1*j1
       ipp = (k*kp*cthnuc/h2 + qqt2*qqt2/h2 + 1./oscp)*j1
     *       + qqt2/hbarc*u1*j0
       ipp = -sn11*4.*pi*exp1*ipp
       tfact = k*kp*dsqrt(1.-cthnuc*cthnuc)*2./u1/h2
       tfact = tfact/sqrt((k*k +kp*kp -2.*k*kp*cthnuc)/h2)
c      tfact =1.
       if (na.eq.16)then
       rofls= (iss + 3.*ipp)/4./pi
       rofls= 8.*2.*rofls*tfact
       else if (na.eq.40)then
       r0 = -4.*pi*j1
       r1 = 4.*pi/oscp*(j1 + u*j0/2.)
c ***  sum = 5.*idd + i2ss
       sum = sn20*(a-c/2.*(k*k+kp*kp+2.*qqt2*qqt2)/h2)*r0
     * + sn20*c*r1 + sn12*15.*coef*coef*r0/2.
     * -15.*sn12*coef*r1
     * -15.*sn12*pi/oscp/oscp*( (3.+u*u/2.)*j1 + u*j0)
       sum = sum*exp1
       rofls= (iss + 3.*ipp + sum)/4./pi
       rofls= 2.*8.*rofls*tfact
       else
       stop
       end if
       end if
       return
       end

c
c
