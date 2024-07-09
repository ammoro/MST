      subroutine vopt0a( vll, k, kp, neles, k0,i1,i2)
c----------------------------------------------------------------------
c *** Calculates optical potential in t*rho approximation
c *** for each cluster
c     R. Crespo and A.M. Moro (2002)
c----------------------------------------------------------------------
         

c----------------------------------------------------------------------
c *** Calculates optical potential in t*rho approximation
c *** for each cluster
c----------------------------------------------------------------------
      
c *** calculates the optical potential for nucleons scattering from a
c     spin 0 nucleus and does a full t*rho projection
***********************************************************************
      use params
      use switch
      use inputs
!      use blocks
!      use wr
      use rhos

      implicit real*8  (a-h, k, m, o-z)
!      logical::  first
      integer::  code, xi, maxls,test
      real*8:: imucl, imusl, imvc, imvs
      real*8:: imu1,imu2
      real*8:: imucl2
      complex*16:: ucl, usl, uc, us,zi, vc,vs
      complex*16:: vll(0:neles,2)
      real*8, allocatable :: pl(:),plp(:)
      complex*16, allocatable:: ucentl(:),uspinl(:)
      parameter (ncos=48)
      dimension   xis( ncos ), wt( ncos )
!      dimension   tapb(4), tamb(4), tcpd(4), tcmd(4), te(4)
!      dimension   tapbm(4), tem(4)
!      dimension   tapb2(4),te2(4)
!      real*8 :: formf(4)
      real*8 :: rhop,rhon
      real*8:: gffr=1,gffi=1
      complex*16:: tapb(2), tamb(2), tcpd(2), tcmd(2), te(2)
      complex*16:: tapbm(2), tem(2),tapb2(2),te2(2)
      
      data  maxls / 100 /

      allocate(pl(0:neles))
      allocate(plp(0:neles))

      zi = (0.,1.)
      
      twopi2 = 2. * pi * pi

!      amass = na * mn
c ** wrong!!! this is not the kmt factor. look inconsistency in
c ** the units !!!!!!
c     aovera = ( amass - 1. )/amass
      aovera = (na -1.)/na

!      if ( nifty(17) .eq. 1 ) aovera = 1.
      if (.not.kmt) aovera=1.
      a = -1.
      b = 1.
      code = 11
      
      write(99,*)'gauss2: ncos=',ncos
      call gauss2( ncos, code, a, b, xis, wt )

       allocate(ucentl(0:neles))
       allocate(uspinl(0:neles))
       ucentl(0:neles)=(0.,0.)
       uspinl(0:neles)=(0.,0.)

!!$       write(90,*)'ncluster=',ncluster
!!$       write(90,*) ffp(1,i1,i2,1), ffn(1,i1,i2,1)
!!$       write(90,*) ffp(2,i1,i2,1), ffn(2,i1,i2,1)
       
              
       do 20 xi = 1, ncos
          x = xis( xi )
          cthnuc = x
          theta = acos( cthnuc ) * 180./pi
          call redish(kp, k, k0, cthnuc, mn, tapb, te)
          q2 = kp**2 + k**2 - 2.* kp * k * cthnuc
          call legpol( x, pl, neles )
          call plprme( x, plp, neles )
!          write(*,*)usnr(1:2),usni(1:2)
          do ncl=1,ncluster
             rhop=ffp(ncl,i1,i2,xi)
             rhon=ffn(ncl,i1,i2,xi)
!             write(90,*) ncl,rhop,rhon
             do 30 l = 0, neles
c *** gaussian formmactor
             if (br>0) gffr=exp(-br*abs(q2)/hbarc/hbarc)
             if (bi>0) gffi=exp(-bi*abs(q2)/hbarc/hbarc)
             ucl=tapb(1)*rhop + tapb(2)*rhon
             ucl = (1./2.) * pl(l) * wt(xi) * ucl
             ucl=ucnr(ncl)*real(ucl) + zi*ucni(ncl)*dimag(ucl)
             if (ucnr(ncl)<1.0) then
             write(*,*)'ncl,ucnr(ncl)=',ncl,ucnr(ncl)
             write(*,*)'ncl,ucni(ncl)=',ncl,ucni(ncl)
             endif
             ucentl(l) = ucentl(l) + ucl
!             if(ncl.eq.1.and.l<1.and.x<1e-3)then !on-shell
! if(xi.eq.1.and.l<1.and.(i1.eq.10).and.(i2.eq.10))then
!                write(90,*)'Cluster:',ncl
!             write(90,'(6g12.3)')theta,tapb(1),tapb(2),rhop,rhon
!                if (k<200) then
!             write(90,'(1f10.2,2x,1g12.5,2x,1g12.5)')k, rhop,rhon
!             endif
!             endif

             if (l.gt.0) then
                usl=-zi*(te(1)*rhop + te(2)*rhon)
                usl=1./((2.*l) * (l + 1.))*usl
                usl= sqrt(1.- x**2) * plp(l) * wt(xi) *usl
                usl=usnr(ncl)*real(usl)+zi*usni(ncl)*dimag(usl)
                uspinl(l) = uspinl(l) + usl
             endif

30             continue !loop in partial waves 
            end do ! end loop clusters

c *** Check: can comment----------------------------------
           if (( kp .eq. k) .and. ( kp .eq. k0 ) ) then
            if (xi .eq. 1 ) then
c              *** set up headers for write of resum ***
               write(kout,1000)
               write(kout,1001)
               write(kout,*)
            endif
            ucl = 0.
            usl = 0.
            do 60 l = 0, neles
               ucl = ucl + ucentl(l) * pl(l) * (2*l+1.)
               usl = usl + uspinl(l) * plp(l) * (2*l+1.)
  60        continue
           write(96,*) k,real(ucl),dimag(ucl)
!!$            revc = tapb(1) * rhop + tapb(2) * rhon
!!$            imvc = tapb(3) * rhop + tapb(4) * rhon
!!$            revs = te(3) * rhop + te(4) * rhon
!!$            imvs = -(te(1) * rhop + te(2) * rhon)
!!$            uc = cmplx( revc, imvc )
!!$            us = cmplx( revs, imvs )
            uc = tapb(1) * rhop + tapb(2) * rhon
            us = -zi*(te(1) * rhop + te(2) * rhon)

            write(kout,1002) theta, q2, abs(ucl), abs(uc),
     $                               abs(usl), abs(us)
         endif
c *** end test ---------------------------------------------
  20  continue

         vll(0,1)=twopi2 * ucentl(0)
         vll(0,2)=vll(0,1)

      do 40 l = 1, neles
            vll(l,1) = twopi2 * ( ucentl(l) + l * uspinl(l) )
            vll(l,2) = twopi2 * ( ucentl(l) - (l+1.) * uspinl(l) )
   40  continue

c *** introduce overall factor
      vll=aovera*vll 

      deallocate(pl)
      deallocate(plp)
      return

 1000 format('1')
 1001 format(' ',' cos  ',9x,'q2',13x,'resum for vcent',11x,'vcent',
     $           11x,'resum for vspin',11x,'vspin')
 1002 format(' ',f5.1,5x,f12.3,9x,e13.6,3(8x,e13.6))
 1003 format(' ',31x,e13.6,3(8x,e13.6))
      end





c ***
c *** ffext: calculates the transformed density rho(q)
c *** 
      subroutine ffext(q2,formf)
      use params
      use bwfnum

      implicit real*8(a-h,k,m,o-z)
      dimension formf(4),bs(220)
!      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
!      common /bwfnum  /wfs(220),wfp(220),densr(220),drx,nramax

      q = sqrt(q2)/hbarc
      do 502 i=1,nramax
      r=(i-1)*drx
      qr=q*r
      if(qr.lt.1.d-10) then
      bs(i)=densr(i)*r*r
      else
      bs(i)=densr(i)*r*r*sin(qr)/qr
      endif
  502 continue
      call sim(bs,dens,1,nramax,drx)
      do 503 i=1,4
      formf(i)=dens
  503 continue
      return
      end



c  @(#)ffact.f 2.3



      subroutine denstgr()
c******************************************************************
      use bdenstg
      implicit real*8(a-h,o-z)
      integer ::kout=16
!      common/bdenstg/nrmxtg,drxtg,rtg(5002),denstg(5002,3)  
       mromx=5002
       read(17,*)rmin,rmax,drxtg
       nrmxtg = (rmax-rmin)/drxtg + 1
! nramax is used but never set??????
c        if (nramax.gt.mromx)then !commented by AMoro
        if (nrmxtg.gt.mromx)then
       write(kout,*)'nramax.gt.mromx'
       stop
       endif
       do 10 i=1,nrmxtg
       read(17,*)rtg(i), denstg(i,1),denstg(i,2),denstg(i,3)
 10    continue
      return
      end

      subroutine ffactb8(q2,formf)
c*********************************************************************
c **  calculate the  density rho(q) for boro 8
      use params
      use bdenstg
      implicit real*8(a-h,k,m,o-z)
      integer :: kout=16
      real*8 bs(5002),formf(4),rhotg(3)
!      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
!      common/bdenstg/nrmxtg,drxtg,rtg(5002),denstg(5002,3)

      allocate(rtg(nrmxtg))
      allocate(denstg(nrmxtg,3))
        
c *******************************************************************
c      Target density normalized to number of nucleons
c >>>  formf(1)-neutrons formf(2)-protons formf(3)-valence proton
      mromx=10000
      
      h2 = hbarc*hbarc
      q = q2/hbarc/hbarc
      do 20 i=1,4
      formf(i) = 0.d0
 20   continue
      do 53 j=1,3
      do 52 i=1,nrmxtg
      qr=q*rtg(i)
      if(qr.lt.1.d-10) then
      bs(i)=denstg(i,j)*rtg(i)*rtg(i)
      else
      bs(i)=denstg(i,j)*rtg(i)*rtg(i)*sin(qr)/qr
      endif

  52  continue
      call sim(bs,res,1,nrmxtg,drxtg)
      rhotg(j) = res*4*pi
      if (rhotg(j).lt.0.00001)rhotg(j)=0.d0
  53  continue
      formf(2)=  3/7. * rhotg(2)
      formf(3)=  4/7. * rhotg(2)
      formf(4) = rhotg(3)
      
      if (q.lt.001)write(kout,*)'rhoq',formf(4),formf(2)+formf(3)

      return
      end


      subroutine ffactli9(q2,formf)
c*********************************************************************
c **  calculate the  density rho(q) for lithium 9
      use params
      use sizes
      use bhalodn
      use brms
      implicit real*8(a-h,k,m,o-z)
      integer :: kout=16
      dimension formf(4)
!      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
!      common /sizes/    achp, acmp, wsp, achn, acmn, wsn
!      common/bhalodn/itydn9,itydn11,itycmdn
!      common/brms/rms9,rms11

c >>>  first executable statement
       h2 = hbarc*hbarc

       do 20 i=1,4
       formf(i) = 0.d0
 20    continue

       if (itydn9.eq.1)then
c ***  core
       rhl = q2*achp*achp/4./h2
       if (rhl.lt.150.0)formf(1)=4*exp(-rhl)
c ***  valence neutrons/proton
       rhl = q2*achn*achn/4./h2
       rhl1 = q2*achn*achn/6./h2
       if(rhl.lt.150.0)formf(2)=4*(1.-rhl1)*exp(-rhl)
       if(rhl.lt.150.0)formf(3)=(1.-rhl1)*exp(-rhl)
       formf(4) = 0.
       else if (itydn9.eq.0)then
       a9 = rms9/sqrt(3.)
       rhl = q2*a9*a9/hbarc/hbarc/2.
       if (rhl.lt.150.0) then
       aux =  9.* exp(-rhl) 
       formf(1) = 6/9. * aux
       formf(2) = 3/9.* aux
       formf(3) = 0.
       formf(4) = 0. 
       end if
       else
       write(kout,*)'wrong itydens'
       STOP 
       end if

      return
      end
      
      subroutine ffactbe10(q2,formf)
c*********************************************************************
c **  calculate the  density rho(q) for lithium 9
      use params
      use sizes
      use bhalodn
      use brms

      implicit real*8(a-h,k,m,o-z)
      dimension formf(4)
!      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
!      common /sizes/    achp, acmp, wsp, achn, acmn, wsn
!      common/bhalodn/itydn9,itydn11,itycmdn
!      common/brms/rms9,rms11
       integer :: kout=16

c >>>  first executable statement
       h2 = hbarc*hbarc
       rmsc = rms9
       rmscha = rms11
       itydnc = itydn9
       itydnha = itydnha
      

       do 20 i=1,4
       formf(i) = 0.d0
 20    continue

       if (itydnc.eq.1)then
c ***   core
        rhl = q2*achp*achp/4./h2
        if (rhl.lt.150.0)formf(1)=4*exp(-rhl)
c ***   valence neutrons/proton
        rhl = q2*achn*achn/4./h2
        rhl1 = q2*achn*achn/6./h2
        if(rhl.lt.150.0)formf(2)=4*(1.-rhl1)*exp(-rhl)
        if(rhl.lt.150.0)formf(3)=2*(1.-rhl1)*exp(-rhl)
        formf(4) = 0.
       else if (itydnc.eq.0)then
        ac = rmsc/sqrt(3.)
        rhl = q2*ac*ac/hbarc/hbarc/2.
        if (rhl.lt.150.0) then
        aux =  10.* exp(-rhl) 
        formf(1) = 8/10. * aux
        formf(2) = 2/10.* aux
        formf(3) = 0.
        formf(4) = 0. 
        end if
       else if (itydnc.eq.2)then
c ***   core
        rhl = q2*achp*achp/4./h2
        if (rhl.lt.150.0)formf(1)=8*exp(-rhl)
c ***   valence neutrons/proton
        rhl = q2*achn*achn/4./h2
        rhl1 = q2*achn*achn/6./h2
        if(rhl.lt.150.0)formf(2)=2*(1.-rhl1)*exp(-rhl)
        formf(3)=0.
        formf(4) = 0.

       else
       write(kout,*)'wrong itydens'
       STOP 
       end if

      return
      end



      subroutine halodn11(q2,ffhalo,ffccm)
c*********************************************************************
      use params
      use bhalodn
      use brms
      implicit real*8(a-h,m,o-z)
      dimension rnmodel(5)
      dimension ap(5), rgp(5), acmp(5), rgcmp(5)
      dimension ah(5), rgh(5), acmh(5), rgcmh(5)
      dimension aso(5), rgso(5), acmso(5), rgcmso(5)
      dimension ai1(5), rgi1(5), acmi1(5), rgcmi1(5)
      dimension ai2(5), rgi2(5), acmi2(5), rgcmi2(5)
!      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
!      common/bhalodn/itydn9,itydn11,itycmdn
!      common/brms/rms9,rms11
      integer :: kout=16

      data rnmodel/5,5,5,5,5/
      data ap/ 9.4551717017698D-02, 0.71304972331458, 1.5060923656690,
     *  -0.30713593782402, -7.2203628101574D-03/
      data rgp/ 0.25000000000000,0.43058676498168, 0.74161984870957,
     *      1.2773267660082, 2.2000000000000/   
      data acmp/ 5.0806336202514D-02, -3.5647883298848D-02,
     *  0.81396573552984, -4.0514984059165D-02, 0.21176273957505/
      data rgcmp/ 1.0000000000000, 1.6817928305074, 2.8284271247462,
     *  4.7568284600109, 8.0000000000000/

      data ai1/ 0.33880679813338, 1.0519101493118,0.64676748481913,
     *   -0.14262661286842, 7.8684879514324D-02 /
      data rgi1/ 0.20000000000000, 0.35565588200778, 0.63245553203368,
     *   1.1246826503807, 2.0000000000000 /
      data acmi1/ 5.7424594236150D-02, 0.48549845843248,
     *  0.41812422896829, 2.8375117376032D-02,9.8898531562224D-03/
      data rgcmi1/1.0000000000000, 1.6817928305074, 2.8284271247462,
     *   4.7568284600109, 8.0000000000000 /

      data ah/ 0.63300296428133, 0.83738369338481,0.52876727485225,
     *   -1.8523641037583D-02, 6.3813961359542D-03 /
      data rgh/ 0.30000000000000, 0.53348382301168, 0.94868329805051,
     *   1.6870239755710, 3.0000000000000/
      data acmh/ 1.9650383889859D-02, 0.18699167232550,0.47233392336042,
     *   0.29997032680215, 2.0628361721867D-02 /
      data rgcmh/1.0000000000000, 1.6817928305074,2.8284271247462,
     *  4.7568284600109, 8.0000000000000/

      data aso/ 0.20017502444560, 0.76135258933939,1.2973806563704,
     *  -0.24003570520186, -1.7023657017945D-02/
      data rgso/ 0.25000000000000, 0.44456985250973,0.79056941504209,
     *   1.4058533129759, 2.5000000000000/
      data acmso/ -1.0807747301582D-03, 0.11115380152594,
     *   0.47926895415049, 0.45672706661538,-4.6525854195310D-02/
      data rgcmso/ 1.0000000000000, 1.6817928305074,
     *   2.8284271247462, 4.7568284600109, 8.0000000000000/

      data ai2/ 0.17535568131644, 0.76816156692909,1.1704145313703,
     *  -0.11387290957603, -9.4603264730216D-04/
      data rgi2/ 0.22000000000000,0.39122147020856,
     * 0.69570108523704, 1.2371509154188, 2.2000000000000/
      data acmi2/-6.3909865402492D-03, 0.25347948777479,
     * 0.53776763872864, 0.17953177621300, 3.5382144677764D-02/
      data rgcmi2/ 1.0000000000000, 1.6817928305074,
     * 2.8284271247462, 4.7568284600109, 8.0000000000000/

c      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
c      common/bhalodn/itydn9,itydn11,itycmdn
c      common/brms/rms9,rms11

c >>>  first executable statement
       h2 = hbarc*hbarc

c ***  itydn11=1 ---- pairing(p)
c ***  itydn11=2 ---- hansen(h)
c ***  itydn11=3 ---- spin-orbit(so)
c ***  itydn11=4 ---- intruderI(i1)
c ***  itydn11=5 ---- intruderII(i2)
c ***  itydn11=6 ---- Zhukov (COSMA)

       if (itydn11.lt.6)then

       imx = rnmodel(itydn11)
       aux = 0.
       auxcm = 0.
       do 20 i = 1,imx
       if (itydn11.eq.1)then
       rhl = (-q2/rgp(i)/rgp(i)/h2)
       rhl1 = 0.
       if(-rhl.lt.150.0)  rhl1 = exp(rhl)
       aux = aux + ap(i)*rhl1
       rhlcm = (-q2/rgcmp(i)/rgcmp(i)/h2)
       rhl1cm = 0.
       if(-rhlcm.lt.150.0)  rhl1cm = exp(rhlcm)
       auxcm = auxcm + acmp(i)*rhl1cm

       elseif (itydn11.eq.2)then
       aux = aux + ah(i)*exp(-q2/rgh(i)/rgh(i)/h2)
       auxcm = auxcm + acmh(i)*exp(-q2/rgcmh(i)/rgcmh(i)/h2)

       elseif (itydn11.eq.3)then
       aux = aux + aso(i)*exp(-q2/rgso(i)/rgso(i)/h2)
       auxcm = auxcm + acmso(i)*exp(-q2/rgcmso(i)/rgcmso(i)/h2)

       elseif (itydn11.eq.4)then
       aux = aux + ai1(i)*exp(-q2/rgi1(i)/rgi1(i)/h2)
       auxcm = auxcm + acmi1(i)*exp(-q2/rgcmi1(i)/rgcmi1(i)/h2)

       elseif  (itydn11.eq.5) then
       aux = aux + ai2(i)*exp(-q2/rgi2(i)/rgi2(i)/h2)
       auxcm = auxcm + acmi2(i)*exp(-q2/rgcmi2(i)/rgcmi2(i)/h2) 
       endif      
 20    continue

       elseif (itydn11.eq.6)then
       b = sqrt( 11.*rms11*rms11/5. - 9.*rms9*rms9/5.) 
       rhl = b*b*q2/h2/4.
       rhl1 = b*b*q2/h2/6.
       aux = 2*(1 -  rhl1)*exp(-rhl) 

       else
       write(kout,*)'wrong type of density for Li11'
       stop
       end if


       ffhalo = aux
       ffccm = auxcm
       if (itycmdn.eq.1)ffccm = 1.
       if (itydn11.eq.6)ffccm = 1.
       return
       end


      subroutine halodnbe(q2,ffhalo,ffccm)
c*********************************************************************
      use params
      use bhalodn
      use brms
      implicit real*8(a-h,m,o-z)
      integer :: kout=16
      dimension rnmodel(2)
      dimension ap(8), rgp(8), acmp(6), rgcmp(6)
      dimension ah(8), rgh(8), acmh(6), rgcmh(6)
!      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
!      common/bhalodn/itydn9,itydn11,itycmdn
!      common/brms/rms9,rms11

      data rnmodel/8,8/
      data ap/ 0.61431992635579, -0.32030496033444, 0.91865405085961,
     *  -9.1035470067696D-02, -0.47385749654132, 0.35560120478844,
     *  -5.2088345323200D-03, -3.8703090858109D-03/
      data rgp/ 0.35000000000000, 0.48632342303060, 0.67574420510914,
     *  0.93894352834790, 1.3046578021102, 1.8128161377309,
     *  2.5188998555040, 3.5000000000000/   
      data acmp/ -1.0113672423068D-02, 4.1626830874519D-02,
     *  3.6491646602639D-02, 0.32702409329388, 0.38645627838139,
     *  0.21817099557942 /
      data rgcmp/ 1.0000000000000, 1.5157165665104, 2.2973967099941,
     *  3.4822022531845, 5.2780316430916, 8.0000000000000/

      data ah/ 0.75797653920326, -0.40816684404758, 1.0836294120775,
     * -0.45885625241802, -0.24148301479338, 0.25008491295767,
     *  1.7851636859113D-02, -8.0535687034861D-03 /
      data rgh/ 0.35000000000000, 0.48632342303060, 0.67574420510914,
     *  0.93894352834790, 1.3046578021102, 1.8128161377309,
     *  2.5188998555040, 3.5000000000000 /
      data acmh/ -9.0352977482308D-03, 3.9604538152411D-02,
     *  6.6663125825300D-02, 0.37341214510238, 0.48433921814049,
     *  4.4691305768633D-02 /
      data rgcmh/ 1.0000000000000, 1.5157165665104, 2.2973967099941,
     *  3.4822022531845, 5.2780316430916, 8.0000000000000/



c >>>  first executable statement
       h2 = hbarc*hbarc
       rmsc = rms9
       rmsha = rms11
       itydnha = itydn11
       itydnc = itydn9

c ***  itydnha=1 ---- Filomena(FN)
c ***  itydnha=2 ---- Ian(IT)

       if (itydnha.lt.3)then

       imx = rnmodel(itydnha)
       aux = 0.
       auxcm = 0.
       
       do 20 i = 1,imx
       if (itydnha.eq.1)then
       rhl = (-q2/rgp(i)/rgp(i)/h2)
       rhl1 = 0.
       if(-rhl.lt.150.0)  rhl1 = exp(rhl)
       aux = aux + ap(i)*rhl1
       elseif (itydnha.eq.2)then
       aux = aux + ah(i)*exp(-q2/rgh(i)/rgh(i)/h2)
       endif
 20    continue
 
       do 25 i = 1,6
       if (itydnha.eq.1)then
       rhlcm = (-q2/rgcmp(i)/rgcmp(i)/h2)
       rhl1cm = 0.
       if(-rhlcm.lt.150.0)  rhl1cm = exp(rhlcm)
       auxcm = auxcm + acmp(i)*rhl1cm
       elseif (itydnha.eq.2)then
       auxcm = auxcm + acmh(i)*exp(-q2/rgcmh(i)/rgcmh(i)/h2)
       endif
 25    continue

       elseif (itydnha.eq.3)then
c ***  not avaiable yet       
       STOP
       b = sqrt( 11.*rmsha*rmsha/5. - 9.*rmsc*rmsc/5.) 
       rhl = b*b*q2/h2/4.
       rhl1 = b*b*q2/h2/6.
       aux = 2*(1 -  rhl1)*exp(-rhl) 

       else
       write(kout,*)'wrong type of density for Be11'
       stop
       end if


       ffhalo = aux
       ffccm = auxcm
       if (itycmdn.eq.1)ffccm = 1.
       if (itydn11.eq.6)ffccm = 1.
       return
       end




      subroutine ffhesw (q2, ff, nff)
************************************************************************
      use params
      use switch
      use sizes
      implicit real*8 (a-h, o-z)
!      real*8 mp, mn

      dimension ff(4)
!      dimension nifty(20)

!      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
!      common /switch/   nifty
!      common /sizes/    achp, acmp, wsp, achn, acmn, wsn

c ***  statement functions now follow
       rj1(r)=sin(r)/r/r - cos(r)/r

c >>>  first executable statement
       h2 = hbarc*hbarc
       rad = 3.26
       argu = sqrt(q2/h2)*rad
       if (argu.lt.0.00000001)then
       aux = 1.
       else
       aux = 3 * rj1(argu)/argu
       end if
       ff(1) = aux
       ff(2) = ff(1)

       return
       end
c
 

      subroutine ffact1(q2,ff,nff)
************************************************************************
      use params
      use sizes

      implicit real*8(a-h,o-z)
!      real*8 mp,mn
      dimension ff(4)

!      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
!      common/sizes/achp,acmp,wsp,achn,acmn,wsn

      h2 = hbarc*hbarc
      ff(1)=0.
      if(na.eq.16)then
        bs = q2/h2/4./0.310
        bp = q2/h2/4./0.336
        if(bs.lt.150.and.bp.lt.150)then
        ff(1) = 4.*exp(-bs) + 12.*(1-2./3.*bp)*exp(-bp)
        end if
        else
        stop
      end if
      ff(2) = ff(1)
      ff(1)=ff(1)/2.
      ff(2)=ff(2)/2.
      return
      end



      subroutine ffacthe (q2, ff, nff)
      end subroutine ffacthe

       subroutine ffactli11(q2,formf)
       end subroutine ffactli11

       subroutine ffhoca(q2,nff,ff1,ff2)
c***********************************************************************
c *** does not divide by n/p form factor
c *** raquel crespo addition
       use params
       use switch
       use sizes

       implicit real*8 (a-h,o-z)
!       real*8 mp, mn
!       dimension nifty(20)
!       common/params/ hbarc, pi, mp, mn, nz, na, nes, nwaves
!       common/switch/nifty
!       common/sizes/achp, acmp, wsp, achn, acmn, wsn
c >>>  first executable statement
       h2 = hbarc*hbarc
c ***  protons
       rhl = q2*achp*achp/4./h2
       rhl1 = q2*acmp*acmp/4./h2
       ff1 = 0.
       if (rhl.lt.150.0)ff1=(1.-rhl1+rhl1*rhl1/5.)*exp(-rhl)
c ***  neutrons
       rhl = q2*achn*achn/4./h2
       rhl1 = q2*acmn*acmn/4./h2
       ff2 = 0.
       if(rhl.lt.150.0)ff2=(1.-rhl1+rhl1*rhl1/5.)*exp(-rhl)
       return
       end



       function f2(q2)
c******************************************************
c*** raquel crespo addition ... *************
c*** this function calculates the square of the density
c*** in momentum space to be used in the calculation of
c*** the second order local optical potential.
   

       use sizes
       use params
       implicit real*8(a-h,o-z)
c       real*8 mp,mn
c       common/sizes/ achp, acmp, wsp, achn, acmn, wsn
c       common/params/ hbarc, pi, mp, mn, nz, na, nes, nwaves
       h2 = hbarc*hbarc
       alpha = (nz-2.)/3.
       ach2=achp*achp
       ach3=ach2*achp
       acm2=acmp*acmp
c
      ro0=2.*acm2/sqrt(pi*pi*pi)/ach3/(2.*acm2+3.*ach2*alpha)
      gama = sqrt(2.)/achp
      gama2 = gama*gama
      a2 = 2.*alpha/acm2
      a3 = alpha*alpha/acm2/acm2
      f2 = ( gama2 + a2*3./2. + a3*15./4./gama2
     #   - q2/h2/gama2*(a2/4. + a3*15./12./gama2)
     #   + q2*q2/h2/h2*a3/16./gama2**3 )
     #   *sqrt(pi*pi*pi)*ro0*ro0/gama2/gama2/gama
      f2 = f2*exp(-q2/h2/4./gama2)
      return
      end

c *******************************************************************
       subroutine ftvoff(ld,nsp,n1)
c *******************************************************************
c      this subroutine calculates fourier transform of the
c      partial wave of the optical potential
c
      use params
      use radgr
      use bgauss3
      use bsdess
      use nlspfl
      use bgrid

      implicit real*8(a-h,k,o-z)
c      real*8 mp,mn
      complex*16 voffr(201,201),zi,sum,sum1,vaux
c      dimension ur(49,49,6), ui(49,49,6)
c      dimension kk(50), wt(48)
c      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
c      common/radgr/radxis(201),radwt(201)
c      common/bgauss3/rmaxr,quin,mquadi,mquado,irmx
c      common/bsdess/sbess(49,201,0:30)
c      common/nlspfl/ur,ui
c      common/bgrid/kk,wt
c
       del = 2
      zi = (0.,1.)
      do 10 i1=1,irmx
      do 10 i2=1,i1
      sum = 0.
        do 20 ik=1,n1-6
        sum1 = 0.
        do 25 jk=1,n1-6
!        vaux = ur(ik,jk,nsp)+zi*ui(ik,jk,nsp)
        vaux = u(ik,jk,nsp)
        sum1 = sum1 + vaux*wt(jk)
     *        *(kk(jk)/hbarc)**2*sbess(jk,i2,ld)
  25    continue
       sum = sum + sum1*wt(ik)*(kk(ik)/hbarc)**2*sbess(ik,i1,ld)
  20  continue
      voffr(i1,i2)=sum
  10  continue
      do 30 i1=1,irmx
      do 30 i2=1,i1
       voffr(i2,i1)=voffr(i1,i2)
  30   continue
       do 40 i1=1,mquadi
       do 40 i2=1,mquadi
c      if(abs(radxis(i1)-radxis(i2)).gt.del)go to 40
       write(14,777)radxis(i1),radxis(i2),real(voffr(i1,i2)),
     *             dimag(voffr(i1,i2))
  40   continue
 777   format(4e12.4)
c
       return
       end





      subroutine vopt0af( vll, k, kp, neles, k0 )
************************************************************************

c *** calculates the full folding optical potential for nucleons
c *** scattering from a double closed shell nucleus
c *** (off-shell calculations for the moment)
c
      use params
      use switch
      use blocks
!      use ampl
      use amps

      implicit real*8  (a-h, k, m, o-z)

      integer  code, xi, maxls
      integer  codep, xip
!      integer  option3

      real*8   imucl, imusl, imvc, imvs
      real*8   imucl2


      complex*16,allocatable::  ucentl(:), uspinl(:)
      complex*16  ucl, usl, uc, us, zi 
      complex*16:: an, cn
!      complex*16 ap, cp
      complex*16 cint2d
      

!      dimension   vll( 0:99, 2, 6 )
      real*8:: vll(0:neles,2,2)
      dimension   xis( 48 ), wt( 48 )
      dimension   xisp(48), wtp(48)
!      dimension   pl( 0:99 ), plp( 0:99 )
      real*8, allocatable :: pl(:),plp(:)
      dimension   tapb(4), tamb(4), tcpd(4), tcmd(4), te(4)
      dimension   tapbf(4), tef(4)
      dimension   tapb2(4),te2(4)
!      dimension   nifty(20)
      dimension   formf(4)

      allocate(pl(0:neles))
      allocate(plp(0:neles))
      allocate(ucentl(0:neles))
      allocate(uspinl(0:neles))

!!$      common   /params/  hbarc, pi, mp, mn, nz, na, nes, nwaves
!!$      common   /switch/ nifty
!!$      common  /block3/option3
!!$      common/ampl/q(200),qq(200),ap(200,200),cp(200,200),nkmx1,nqmx1

      data  maxls / 100 /

         nff = 1
         ncos = 48
c ***    ncos changed to 48 to increase accuracy
c        ncos =48
         ncosp = 48
         nn = na - nz
         aovera = (na -1.)/(na*1.)
         twopi2 = 2. * pi * pi
         if ( nifty(17) .eq. 1 ) aovera = 1.
         a = -1.
         b = 1.
         aa = 0.d0
c        bb = max(k,kp)/hbarc
c        bb = 3.0*k0/hbarc
         bb = 3.0*k0/hbarc
         code = 11
         call gauss2( ncos, code, a, b, xis, wt )
         call gauss2( ncosp, code, aa, bb, xisp, wtp )



       tfact = (na)/(na+1.)*mn*hbarc**2*((2*pi)**3)*1.38/k0
       tfact = tfact*na*na/4.
       if (option3.eq.2)then
       t2 = 0.d0
       tt2 = (na+1.)*k0/(2.*na)
       tt2 = tt2/hbarc
       tapb2(1) = real(cint2d(qq,q,ap,tt2,t2,nkmx1,nqmx1,5,200))
       tapb2(2) = tapb2(1)
       tapb2(3) = dimag(cint2d(qq,q,ap,tt2,t2,nkmx1,nqmx1,5,200))
       tapb2(4) = tapb2(3)
       do 5 i=1,4
       tapb2(i) = -tapb2(i)/2./pi/pi/mp/hbarc
  5   continue
       endif
c
      zi = ( 0. , 1.)

      do 10 l = 0, maxls - 1
         ucentl(l) = ( 0., 0. )
         uspinl(l) = ( 0., 0. )
  10  continue

      do 20 xi = 1, ncos
         x = xis( xi )
         cthnuc = x
         q2 = kp**2 + k**2 - 2.* kp * k * cthnuc
         if ( q2 .lt. 0 ) q2 = -q2
         call ffact(q2, formf, nff)
         call legpol( x, pl, neles )
         call plprme( x, plp, neles )
c
c ***  qt, qqt, qqp are in units fm-1
c *** (transfer and total nucleon-nucleon momentum)
c *** qqp used as a check to this subroutine
       qqp = k*k + kp*kp + 2.*k*kp*cthnuc
       if(qqp.lt.0)qqp=-qqp
       qqp = sqrt(qqp)/2./hbarc
c *** no recoil effects
c      qqp = qqp*(na+1.)/(na*1.)/2.
       qqp = qqp/2.
       qt = sqrt(q2)/hbarc
c
       do 24 i=1,4
       tapbf(i)=0.d0
       tef(i)=0.d0
  24   continue
        do 25 xip = 1,ncosp
        qqt = xisp(xip)
        qqp = qqt

       tapb(1)=real(cint2d(qq,q,ap,qqp,qt,nkmx1,nqmx1,5,200))
c      tapb(2)=real(cint2d(qq,q,an,qqp,qt,nkmx2,nqmx2,5,200))
       tapb(2)=tapb(1)
       tapb(3)=dimag(cint2d(qq,q,ap,qqp,qt,nkmx1,nqmx1,5,200))
c      tapb(4)=dimag(cint2d(qq,q,an,qqp,qt,nkmx2,nqmx2,5,200))
       tapb(4)=tapb(3)
       te(1)=real(cint2d(qq,q,cp,qqp,qt,nkmx1,nqmx1,5,200))
c      te(2)=real(cint2d(qq,q,cn,qqp,qt,nkmx2,nqmx2,5,200))
       te(2)=te(1)
       te(3)=dimag(cint2d(qq,q,cp,qqp,qt,nkmx1,nqmx1,5,200))
c      te(4)=dimag(cint2d(qq,q,cn,qqp,qt,nkmx2,nqmx2,5,200))
       te(4)=te(3)
c **   change from scattering amplitude to transition amplitude
cc **   and to landau conventions
       do 300 i=1,4
       tapb(i) = -tapb(i)/2./pi/pi/mp/hbarc
       te(i) = -te(i)/2./pi/pi/mp/hbarc
 300   continue
c
       do 27 i=1,4
       tapbf(i) = tapbf(i) + tapb(i)*rofc(kp,k,k0,cthnuc,qqt)
     *                       *qqt*qqt*wtp(xip)
       tef(i) = tef(i) + te(i)*rofls(kp,k,k0,cthnuc,qqt)
     *                   *qqt*qqt*wtp(xip)
  27    continue
  25    continue

         do 30 l = 0, neles
            if ( l .eq. 0 ) then
               uspinl(l) = ( 0., 0. )
            else
       reusl = (tef(3) + tef(4))
c      reusl = nz*te(3)*formf(1) + nn*te(4)*formf(1)
               reusl = 1./( (2.*l) * (l + 1.) ) * reusl
               reusl = reusl * sqrt(1.- x**2) * plp(l) * wt(xi)
       imusl = (tef(1) + tef(2))
c      imusl = nz*te(1)*formf(1) + nn*te(2)*formf(1)
               imusl = -1./( (2.*l) * (l + 1.) ) * imusl
               imusl = imusl * sqrt(1.- x**2) * plp(l) * wt(xi)
               uspinl(l) = uspinl(l) + cmplx( reusl, imusl )
            endif
       reucl = (tapbf(1) + tapbf(2))
c      reucl = nz*tapb(1)*formf(1) + nn*tapb(2)*formf(1)
            reucl = 1./2. * reucl
            reucl = reucl * pl(l) * wt(xi)
       imucl = (tapbf(3) + tapbf(4))
c      imucl = nz*tapb(3)*formf(1) + nn*tapb(4)*formf(1)
            imucl = 1./2. * imucl
            imucl = imucl * pl(l) * wt(xi)
          if (option3.eq.2)then
         reucl2 = -2.*(tapb2(1)*tapb2(3))*tfact*f2(q2)
         reucl2 = 1./2. * reucl2
         reucl2 = reucl2 * pl(l) * wt(xi)
         imucl2 = (tapb2(1)**2-tapb2(3)**2)*tfact*f2(q2)
         imucl2 = 1./2. * imucl2
         imucl2 = imucl2 * pl(l) * wt(xi)
         reucl=reucl+reucl2
         imucl=imucl+imucl2
         end if
            ucentl(l) = ucentl(l) + cmplx( reucl, imucl )
  30     continue
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
      deallocate(pl,plp)
      return
      end




      subroutine vopthe8( vll, k, kp, neles, k0 )
************************************************************************


c *** calculates the optical potential for nucleons scattering from a
c     Helium 8  nucleus and does a full t*rho projection

      use params
      use inputs
      use wr
      use halo
      use switch

      implicit real*8  (a-h, k, m, o-z)
      integer  code, xi, maxls
      real*8   imucl, imusl, imvc, imvs
      real*8   imu1,imu2
      real*8   imucl2
      real*8 ::vll(0:neles,2,6)

!      complex*16  ucentl( 0:99), uspinl( 0:99 ), zi
      complex*16,allocatable:: ucentl(:), uspinl(:)
      complex*16  ucl, usl, uc, us ,zi
      dimension   xis( 48 ), wt( 48 )


      real*8::   pl(0:neles), plp(0:neles)
      dimension   tapb(4), tamb(4), tcpd(4), tcmd(4), te(4)
      dimension   tapbm(4), tem(4)
      dimension   tapb2(4),te2(4)
      dimension   formf(4)

!      dimension   nifty(20)
!      dimension   vll( 0:99, 2, 6 )
!      dimension   pl( 0:99 ), plp( 0:99 )

!!$      common   /params/  hbarc, pi, mp, mn, nz, na, nes, nwaves
!!$      common   /switch/ nifty
!!$      common  /inputs/tlab
!!$      common  /wr/iwrite,iread,istrong
!!$      common /halo/ioptls,ioptcnt

     
      allocate(ucentl(0:neles))
      allocate(uspinl(0:neles))
!      allocate(pl(0:neles))
!      allocate(plp(0:neles))

      data  maxls / 100 /
        nhe4 = 4
        nhev = 4
        ncos = 48
        nn = na - nz
        amass = na * mn
        aovera = (na -1.)/(na*1.)
c ***   set ustrong=0
!        if (istrong.eq.1)aovera=0.
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
         call ffacthe( q2, formf, 2 )
         call legpol( x, pl, neles )
         call plprme( x, plp, neles )

         do 30 l = 0, neles
            if ( l .eq. 0 ) then
               uspinl(l) = ( 0., 0. )
            else
               reusl=nhe4 * te(3) * formf(1) + nhev * te(4) * formf(2)
               reusl = 1./( (2.*l) * (l + 1.) ) * reusl
               reusl = reusl * sqrt(1.- x**2) * plp(l) * wt(xi)
               imusl=nhe4 * te(1) * formf(1) + nhev * te(2) * formf(2)
               imusl = -1./( (2.*l) * (l + 1.) ) * imusl
               imusl = imusl * sqrt(1.- x**2) * plp(l) * wt(xi)
               uspinl(l) = uspinl(l) + cmplx( reusl, imusl )
            endif
            reucl=nhe4 * tapb(1) * formf(1) + nhev*tapb(2)*formf(2)
            reucl = 1./2. * reucl
            reucl = reucl * pl(l) * wt(xi)
            imucl=nhe4 * tapb(3) * formf(1) + nhev*tapb(4)*formf(2)
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


      subroutine voptli( vll, k, kp, neles, k0 )
************************************************************************


c *** calculates the optical potential for nucleons scattering from a
c     Lithium  nucleus and does a full t*rho projection

      use params
      use switch
      use inputs
      use halo
      use wr

      implicit real*8  (a-h, k, m, o-z)
      integer  code, xi, maxls
      real*8   imucl, imusl, imvc, imvs
      real*8   imu1,imu2
      real*8   imucl2
!      complex*16  ucentl( 0:99), uspinl( 0:99 )
      complex*16  ucl, usl, uc, us, zi

!      dimension   vll( 0:99, 2, 6 )
      real*8 ::  vll( 0:neles, 2, 6 )
      complex*16,allocatable:: ucentl(:), uspinl(:)
      dimension   xis( 48 ), wt( 48 )
!      dimension   pl( 0:99 ), plp( 0:99 )
      real*8, allocatable:: pl(:), plp(:)
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

        nlicav = 1
        nlicn = 1
        nlicp = 1
        nliv = 1
        nlivls = 1
        nlivc = 1
        if (ioptls.eq.2) nlivls = 0.
        if (ioptcnt.eq.1) nlivc = 0.
        ncos = 48
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
         if (na.eq.9) call ffactli9( q2, formf)
         if (na.eq.11) call ffactli11( q2, formf)

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
! nhev is used but never set ?????
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












c
      subroutine vopt2l( v2l, k, kp, neles, k0 )
************************************************************************
      use params
      use switch
      use inputs
      use blocks 
      use wr
      use bsflip
      use btmatt

      implicit real*8  (a-h, k, m, o-z)

      integer  code, xi, maxls
!      integer option1,option5

      real*8   imucl, imusl, imvc, imvs
      real*8 imu1,imu2
      real*8   imucl2

      logical  first

      complex*16, allocatable ::  ucentl(:), uspinl(:)
      complex*16  ucl, usl, uc, us, zi
!      complex*16  tamp2,tmatt
      complex*16:: tamp2

!      real*8, allocatable::   v2l(:,:,:)
      real*8 :: v2l(0:neles,2,2)
      dimension   xis( 48 ), wt( 48 )
!      dimension   pl( 0:99 ), plp( 0:99 )
      real*8,allocatable::   pl(:), plp(:)
      dimension   tapb(4), tamb(4), tcpd(4), tcmd(4), te(4)
      dimension   tapbm(4), tem(4)
      dimension   tapb2(4),te2(4)
!      dimension   nifty(20)
      dimension   formf(4)
      dimension   famp(10)

      allocate(pl(0:neles))
      allocate(plp(0:neles))
      allocate(ucentl(0:neles))
      allocate(uspinl(0:neles))


!!$      common   /params/  hbarc, pi, mp, mn, nz, na, nes, nwaves
!!$      common   /switch/ nifty
!!$      common /inputs/   tlab , b , ymin1 , ymin2 ,
!!$     $                  kode , lxmax , nang , ngp , nr
!!$      common  /block3/option3
!!$      common  /block1/option1
!!$      common  /wr/iwrite,iread,istrong
!!$      common  /bsflip/isflip
!!$      common  /btmatt/tmatt(0:7,10,50,50)

      data  first / .true. /
      data  maxls / 100 /

      if ( first ) then
         first = .false.
         ncos = 48
         nn = na - nz
         amass = na * mn
         aovera = (na -1.)/(na*1.)
c ***  set ustrong=0
         if (istrong.eq.1)aovera=0.
         twopi2 = 2. * pi * pi
         if ( nifty(17) .eq. 1 ) aovera = 1.
         a = -1.
         b = 1.
         code = 11
         call gauss2( ncos, code, a, b, xis, wt )

      endif

       tfact = (na)/(na+1.)*mn*hbarc**2*((2*pi)**3)*1.38/k0
       tfact = tfact*na*na/4.
        call redish(k0,k0,k0,1.d0,mn,tapb2,te2)

      zi = ( 0. , 1.)
      n1 = ngp+1

      do 10 l = 0, maxls - 1
         ucentl(l) = ( 0., 0. )
         uspinl(l) = ( 0., 0. )
  10  continue
c
      do 15 ir = 1,10
        if(ir.le.5)then
        famp(ir) = 1.
        if(ir.eq.3)famp(ir)=famp(ir)*2.
        else
        famp(ir) = 3.
        if(ir.eq.8)famp(ir)=famp(ir)*2.
        endif
 15   continue
      twopi2 = 2.*pi*pi
      h3 = hbarc**3
      c2 = 1./twopi2/twopi2/h3/h3
      if(isflip.eq.0)then
        tamp2 = tmatt(0,1,n1,n1)*tmatt(0,1,n1,n1)*c2
      else
        tamp2 = 0.
        do 17 ir=1,10
        tamp2 = tamp2 + tmatt(0,ir,n1,n1)*tmatt(0,ir,n1,n1)*c2
     *                 * famp(ir)
 17     continue
      endif
      ret2 = real(zi*tamp2)
      aimgt2 = dimag(zi*tamp2)


      do 20 xi = 1, ncos
         x = xis( xi )
         cthnuc = x
         theta = acos( cthnuc ) * 180./pi
         q2 = kp**2 + k**2 - 2.* kp * k * cthnuc
         if ( q2 .lt. 0 ) q2 = -q2
         call legpol( x, pl, neles )
         call plprme( x, plp, neles )

c ***
           if(k.eq.kp.and.k.eq.k0)then
            call ffact(q2,formf,1)
          write(6,*)'cthnuc= ',cthnuc,'  f2=  ',f2(q2)
         write(6,*)'cthnuc= ',cthnuc,'  rho=  ',formf(1)
           endif
         do 30 l = 0, neles
c        reucl2 = -2.*(tapb2(1)*tapb2(3))*tfact*f2(q2)
         reucl2 = ret2 * tfact *f2(q2)
         reucl2 = 1./2. * reucl2
         reucl2 = reucl2 * pl(l) * wt(xi)
c        imucl2 = (tapb2(1)**2-tapb2(3)**2)*tfact*f2(q2)
         imucl2 = aimgt2*tfact*f2(q2)
         imucl2 = 1./2. * imucl2
         imucl2 = imucl2 * pl(l) * wt(xi)
            ucentl(l) = ucentl(l) + cmplx( reucl2, imucl2 )
  30     continue
  20  continue


      do 40 l = 0, neles
            v2l(l,1,1) = twopi2 * ucentl(l)
            v2l(l,2,1) = twopi2 * dimag( ucentl(l) )
            v2l(l,1,2) = v2l(l,1,1)
            v2l(l,2,2) = v2l(l,2,1)
         do 50 ir = 1,2
            do 50 nspin = 1, 2
               v2l(l, ir, nspin) = aovera * v2l( l, ir, nspin )
  50     continue
  40  continue
c ***************************************
      deallocate(pl,plp)
      return

c     *** formats ***

 1000 format('1')
 1001 format(' ',' cos  ',9x,'q2',13x,'resum for vcent',11x,'vcent',
     $           11x,'resum for vspin',11x,'vspin')
 1002 format(' ',f5.1,5x,f12.3,9x,e13.6,3(8x,e13.6))
 1003 format(' ',31x,e13.6,3(8x,e13.6))
      end



       subroutine vcoul(k,kp,vl,neles,type)
c***************************************************************
c
c ***  this subroutine determines partial wave decomposition of
c ***  coulomb potential, with a realistic charge form factor
c ***  raquel crespo addition

c **** type=0 point (screened)
c          =3 realistic (HO)
c          =4 uniform sphere 

       use params
       use ranges
       use sizes
       use inputs

         
       implicit real*8 (a-h,k,m,o-z)
       integer :: code, xi, type
       dimension xis(48), wt(48)

        real*8:: vl(0:neles)  
! 
!        real*8:: vl(:) !AMoro 21/05/04
!        real*8:: vl(*) ! with this line coulomb=3 crashes! 
       real*8, allocatable :: pl(:)
       real*8:: formf
       allocate(pl(0:neles))

       h2 = hbarc*hbarc
       twopi2 = 2.*pi*pi
       fact = nz/2./137.
       ra = 1.3*na**(1./3.)

!       gama = (nz-2.)/(6.*nz)

       nff = 1
       ncos = 48
       a = -1.
       b = 1.
       code = 011
       vl(0:neles)=0.

       call gauss2(ncos, code, a, b, xis, wt)
       do 20 xi=1,ncos
          cthnuc = xis(xi)
          q2 = kp**2 + k**2 - 2.*kp*k*cthnuc
          if (q2.lt.0) q2 =-q2
          q = sqrt(q2)
          call legpol(cthnuc,pl,neles)

          select case(type)
          case (0) ! point coulomb
             formf = 1.
          case (3) ! charge form factor
             rhl1=1.-gama*q2/h2-alfa*(q2/h2)**2
!             rhl = q2*achp*achp/4./h2
             rhl=beta*q2/h2
             formf = 0.
             if(rhl.lt.150.) formf=rhl1*exp(-rhl)
          case (4) ! uniform spherical charge density
             if (q.lt.1)then
                formf = 3.
             else
                x = q*ra/hbarc
                formf = 3.*(sin(x)/x - cos(x))/x/x
             end if
          end select


       if (q.lt.1.)then
          select case(type)
          case(0)
             vofq=fact*rcut**2/2.
          case(3)
!             vofq=fact*(-achp*achp/4.-
!     &            gama*achp*achp + rcut*rcut/2.)
              vofq=fact*(-beta-gama+rcut**2/2.)
          case(4)
             vofq=fact*(-ra**2/10. + rcut**2/2.)
          end select
       else
          vofq = fact * (formf - cos(q*rcut/hbarc))/q2
       end if

       do i=0,neles
          vl(i) = vl(i) + vofq*pl(i)*wt(xi)
       enddo
20     continue
       deallocate(pl)
       return
       end



       subroutine vcoul2(k,kp,vl,neles)
c***************************************************************
c
c ***  this subroutine determines partial wave decomposition of
c ***  coulomb potential, with a realistic charge form factor
c ***  raquel crespo addition
c
       use params
       use ranges
       use switch
       use wr
       use bhalodn
       use brms
       use sizes
       use inputs

         
       implicit real*8 (a-h,k,m,o-z)
       integer code, xi
       dimension xis(48), wt(48)
       real*8:: vl(0:neles)
       real*8, allocatable :: pl(:)
       dimension formf(4)
       dimension ff9(4),ffc(4)

       allocate(pl(0:neles))

!!$       common/params/ hbarc, pi, mp, mn, nz, na, nes, nwaves
!!$       common/sizes / achp, acmp, wsp, achn, acmn, wsn
!!$       common/ranges/ rcoul, rcut
!!$       common/switch/ nifty
!!$       common/wr/iwrite,iread,istrong
!!$       common/bhalodn/itydn9,itydn11,itycmdn
!!$       common/brms/rms9,rms11

       h2 = hbarc*hbarc
       twopi2 = 2.*pi*pi
       fact = nz/2./137.
       ra = 1.3*na**(1./3.)
       rmsc = rms9
       rmsha = rms11
       itydnc = itydn9
       itydnha = itydn11

       if (nz.eq.2)then
       gama = 0.
       else if (nz.gt.2)then
       if (na.le.39)then
       gama = (nz-2.)/(6.*nz)
       else if (na.eq.40)then
       gama = 1/4.
       else
       stop
       end if
       end if

       nff = 1
       ncos = 48
       a = -1.
       b = 1.
       code = 011
       do 10 i=0,neles
       vl(i) = 0.d0
 10    continue
       call gauss2(ncos, code, a, b, xis, wt)

       do 20 xi=1,ncos
       cthnuc = xis(xi)
       q2 = kp**2 + k**2 - 2.*kp*k*cthnuc
       if (q2.lt.0) q2 =-q2
       q = sqrt(q2)
       call legpol(cthnuc,pl,neles)
c ***  charge form factor
       if(coulomb.eq.3)then
        rhl = q2*achp*achp/4./h2

        if (nz.eq.2)then
           if (na.eq.4) then
           formf(1) = 0.
           rhl1 = (acmp*acmp*q2/h2)**6
           if(rhl.lt.150.0) formf(1)=(1.-rhl1)*exp(-rhl)
           if(iwrite.eq.1) formf(1)=1.
           else if (na.eq.8)then
           formf(1) = 0.
           if(rhl.lt.150.0) formf(1)=exp(-rhl)
!           if(iwrite.eq.1) formf(1)=1. 
           end if  
        else if (nz.eq.3) then
          call ffactli9(q2,ff9)
          formf(1) = (ff9(1)+ff9(2)+ff9(3))/9.
          if(iwrite.eq.1) formf(1)=1. 
        else if (nz.eq.4) then
          call ffactbe10(q2,ffc)
          formf(1) = (ffc(1)+ffc(2)+ffc(3))/10.
!          if(iwrite.eq.1) formf(1)=1.          
                 
        else if (nz.gt.4)then
          if (na.le.39)then
          rhl1 = achp*achp*q2/(6.0*nz)/h2
          formf(1) = 0.
          if(rhl.lt.150.0) formf(1)=(1.-(nz-2.)*rhl1)*exp(-rhl)
          if(iwrite.eq.1) formf(1)=1.
          else if(na.eq.40)then
          formf(1) = 0.
          rhl1=(-achp*achp*q2/h2/4.+(achp*achp*q2/h2)**2/80.)
          if(rhl.lt.150.)formf(1)=(1.+rhl1)*exp(-rhl)
!          if (iwrite.eq.1) formf(1)=1.
          else
          stop
          end if
        end if

       else if(coulomb.eq.4)then
        if (q.lt.1)then
        formf(1) = 3.
        else
        x = q*ra/hbarc
        formf(1) = 3.*(sin(x)/x - cos(x))/x/x
!        if (iwrite.eq.1) formf(1)=1.
        end if
       else
       stop
       end if
c ***  coulomb potential
       if (q.lt.1.)then
        if(coulomb.eq.3)then

        if (nz.eq.3)then
          if (itydn9.eq.0)then
          a9 = rms9/sqrt(3.)
          vofq = (-a9*a9 + rcut*rcut)/2.
          else if  (itydn9.eq.1)then
          vofq = -(achp**2/9. + 5*achn**2/(9*4.) + 5*achn**2/(9*6.))
     *            + rcut*rcut/2.
          end if
          else if (nz.eq.4)then
          if (itydnc.eq.0)then
          ac = rmsc/sqrt(3.)
          vofq = (-ac*ac + rcut*rcut)/2.
          else if  (itydnc.eq.1)then
          vofq = -(achp**2/10. + 6*achn**2/(10*4.) + achn**2/(10.))
     *            + rcut*rcut/2.
          end if         
          
        else
        vofq=fact*(-achp*achp/4.-gama*achp*achp+rcut*rcut/2.)
        end if

        else if(coulomb.eq.4)then
        vofq = fact*(-ra**2/10. + rcut**2/2.)
        else
        stop
        end if
!        if (iwrite.eq.1)vofq=fact*rcut**2/2.
       else
       vofq = fact * (formf(1) - cos(q*rcut/hbarc))/q2
       end if
        do 40 i=0,neles
        vl(i) = vl(i) + vofq*pl(i)*wt(xi)
 40     continue
 20    continue
        deallocate(pl)
       return
       end






      subroutine wsax(q2,formf)
        use params
c*****************************************************************
c***  calculates fourier transform of 2pf for oxygen
c
      implicit real*8 (a-h,k,m,o-z)
      dimension formf(4)
!      common /params/hbarc, pi, mp, mn, nz, na, nes, nwaves

      q = sqrt(q2)/hbarc
      if (q.eq.0)then
      do 5 i=1,4
      formf(i) = 1.
  5   continue
      return
      endif
c     a = 0.55
c     r = 3.25
      a=0.389
      r=1.1*(na*1.)**(1./3.)
      exp3 = exp(-r/a)
      v0 = 4.*pi*(r**3)/3.
      v0 = v0*(  1. + pi*pi*(a/r)*(a/r)
     *           + 6.*(a/r)**3*exp3*(1-exp3/8.) )
      exp1 = exp(-pi*a*q)
      exp2 = exp(-2.*pi*a*q)
      v1= 8.*pi*pi*a*a*exp1/q/(1-exp2)/(1-exp2)
      v1= v1*(pi*(1.+exp2)*sin(q*r) - r*(1.-exp2)*cos(q*r)/a)
      v2= 8.*pi*a**3*exp3/(1.+(q*a)**2)**2
      v2= v2*(  1.-2.*exp3
     *     *( (1.+(q/a)**2)/(4.+(q/a)**2))**2  )
      do 10 i=1,4
      formf(i) = (v1+v2)/v0
  10  continue
      return
      end

c-------------------------------------------------------------------
c *** Calculates all form factors rho(cluster,q)
c-------------------------------------------------------------------
      subroutine formfactors(n1,n2) 
      use rhos
      use bgrid
      use inputs
      use params

      implicit real*8  (a-h, k, m, o-z)
      integer :: code,xi,type,shape
      integer :: nzclus,nnclus
      real*8  ::  ucrnorm,ucinorm,usrnorm,usinorm
      parameter (ncos=48)
!      real*8 ::xis( ncos ),ws(ncos) !, wt( 32 )
      dimension   formf(4)
            
      namelist /cluster/type,shape,nzclus,nnclus,ahop,ahon,
     &        ucrnorm,ucinorm,usrnorm,usinorm

      write(99,*) '+ Allocating',ncluster,n1,n1,ncos,'for ffp' 
 
      allocate(ff(ncluster,n1,n1,ncos,2))
      if (ncluster>10) then
         write(*,*)'** formfactors ** Please, increase dimension'
         write(*,*)'of usnr and usni in modules.f90'
         stop
      endif
      ff=0d0
      ucnr(:)=1.d0
      ucni(:)=1.d0
      usnr(:)=1.d0
      usni(:)=1.d0
      do ncl=1,ncluster
         ucrnorm=0d0;ucinorm=0d0
         usrnorm=0d0;usinorm=0d0
         read(kin,nml=cluster)
!         write(kout,nml=cluster)
         if (ucrnorm>1e-6) ucnr(ncl)=ucrnorm
         if (ucinorm>1e-6) ucni(ncl)=ucinorm
         if (usrnorm>1e-6) usnr(ncl)=usrnorm
         if (usinorm>1e-6) usni(ncl)=usinorm
!         write(*,*)'norms:',ucrnorm,ucinorm,usrnorm,usinorm
         
         select case(type)
            case (0)
               write(*,'("- Cluster",1i2,"->numerical density")')ncl
               if (shape.eq.5) then 
                  call ffactextr(ncl,n1,n2)
               else if (shape.eq.6) then
                  call ffactextq(ncl,n1,n2)
               else
                  call ffactbs(ncl,shape,n1,n2)
               endif
            case (1)
               select case(shape)
                  case(0) 
                     write(*,'(" -Cluster",1i2,"->HO density")') ncl
                     call ffactho(ncl,n1,n2)
                  case(1)
                     write(*,'(" -Cluster",1i2,"->Fermi density")')ncl
!                     call ffactfermi(ncl,n1,n2,nzclus,nnclus)
                      call ffactfermi2(ncl,n1,n2,nzclus,nnclus)
                  case(2)
                     write(*,'(" -Cluster",1i2,"->G3 density")')ncl
                     call ffactg3(ncl,n1,n2,nzclus,nnclus)
                  case default
                     write(*,*) 'Shape',shape,
     &            'not implemented for type',type,'potential'  
              end select
           case default
              write(*,*) 'Type',type,'undefined!'
              stop
           end select
      enddo
      return
      end

c--------------------------------------------------------------------
c *** Calculates numerical density by calculating the single-particle
c *** orbitals using an analytical binding potential
c--------------------------------------------------------------------
      subroutine ffactbs(ncl,ib,n1,n2)
        use params
        use rhos
        implicit real*8 (a-h,k,o-z)
      
        integer :: norba
        integer :: kin=13
        integer :: ia,ib,ic,is2,j2a,lmoma,nodd,nramax,jhw,kcheck
        real*8  :: cmass,vmass,wzz,wal,wr0,wls,vdepth,bengy,dmat
        real*8  :: wrz,drx,pnloc,wr0ls,wals
        real*8 ::chis(10)
        real*8,allocatable ::ffr(:,:)
        real*8,allocatable ::bs(:),u(:,:)
        real*8,allocatable :: densr(:,:),denstot(:),dd(:)
        parameter (ic=1,jhw=1,pnloc=0.d0)

        namelist /bsparm/ cmass,vmass,zc,zv,
     &                  nramax,drx,dmat,nshell

        namelist /bshell/ ia,bengy,vdepth,
     &                  wr0,wal,wls,
     &                  j2a,lmoma,is2,nodd,norba


*----------------------------------------------------------------------
*     read single particle bound state formfactors
*----------------------------------------------------------------------
       open(19,file='wfn.out',status='unknown')
     
       do 887 iso=1,2 
          bengy = 0.d0
          vdepth = 0.d0
          cmass=0.
          zv=0.
          read(kin,nml=bsparm)
          allocate(densr(nshell,nramax))
          allocate(u(nshell,nramax))
          allocate(ffr(nramax,4))
          allocate(bs(nramax))
          allocate(dd(nramax))
          allocate(denstot(nramax))

          u=0d0
          write(99,nml=bsparm)

          if (cmass<1.e-5) cycle
          denstot=0.
c *** loop from outer to inner shells
          write(*,*) 'Reading',nshell,'nshells...'
          do 888 ibs =1,nshell
            read(kin,nml=bshell)
*----------------------------------------------------------------------
*          ' --------------------------------------'
*          ' Bound states in a two body potential  '
*          ' ------------------------------------  '
*     search well depth for ia=0, search energy for ia=1
*     separation energy: bengy
*     potential depth: vdepth
*     core and bound particle masses: cmass,vmass
*     core and bound particle charges: zc,zv
*     potential shape ib=
*                       0 = Woods-Saxon
*                       1 = Gauss
*                       2 = Yukawa
*                       3 = Hulthen
*                       4 = Cosh
*     potential radius and diffuseness: wr0,wal
*     spin-orbit potential strength (~6 MeV): wls
*     2*j value, lmom, and nodes (from 0): j2a,lmoma,nodd
*     2*valence partice spin: is2
*     no. integ steps and step: nramax,drx
*     dmat (fine tuning of interior-exterior matching. Usually input 0. 
*           If no state found try increase of dmat to of order 10.
*     wzz is zc*zv
*--------------------------------------------------------------------------
      wzz=zc*zv
      wr0ls=wr0
      wals=wal
      wrz=wr0
*------------------------------------------------------------------------
      call bound(ia,ib,ic,is2,j2a,lmoma,nodd,nramax,jhw,cmass,
     1 vmass,wzz,wal,wr0,wls,ffr,vdepth,bengy,dmat,kcheck,wrz,
     2 chis,drx,pnloc,wr0ls,wals)
*------------------------------------------------------------------------
!      write(6,21) vdepth,bengy,wal,wr0,wzz,wls,ia,ib,ic,
!     1 lmoma,nodd,j2a,is2,dmat
   21 format('  vdepth= ',f10.5,'  bengy= ',f10.5,'  wa= ',f10.5,
     + /,'  wr0= ',f10.5,'  wzz= ',f10.5,'  wls= ',f10.5,/,
     +   '  ia= ',i5,'  ib= ',i5,'  ic= ',i5,'  lmom= ',i5,
     +   '  nod= ',i5,/,'  2*j=',i5,' 2*s= ',i5,' dmat= ',f10.5)
*------------------------------------------------------------------------
*     wave function at the origin
*------------------------------------------------------------------------
      u(ibs,1)=chis(1)
      if(ibs.lt.2) u(ibs,1)=0.d0
*------------------------------------------------------------------------
      do 44 i=1,nramax
      u(ibs,i+1)= ffr(i,1)
   44 continue
*------------------------------------------------------------------------
*     print, check normalizations, and calculate rms of relative motion
*------------------------------------------------------------------------
      do 86 i=1,nramax
      r=(i-1)*drx
      bs(i)=(r*u(ibs,i))**2
      densr(ibs,i)=bs(i)*norba
      dd(i) = densr(ibs,i)
   86 continue
      call sim2(bs,rnorm,1,nramax,drx,nramax)
      call sim2(dd,rnd,1,nramax,drx,nramax)
      write(6,fmt='(" Norm=",2f6.4)') rnorm,rnd 
      do 90 i=1,nramax
      r=(i-1)*drx
      denstot(i)=denstot(i)+norba*u(ibs,i)**2
      bs(i)=bs(i)*r*r
   90 continue
      call sim2(bs,rnorm,1,nramax,drx,nramax)
      write(6,'(" Wave function relative rms",1f8.3)') sqrt(rnorm) 
      

 888  continue
 
*---------------------------------------------------------------------
*     print bound wave functions and density
*---------------------------------------------------------------------
      write(19,*)'wave functions'
      write(14,*) '#', nramax,drx
      do 507 i=1,nramax
      rad = (i-1)*drx
      write(14,*) rad,denstot(i)
      write(19,508) rad,u(1,i),u(2,i)
 507  continue
 508  format(e12.5,3x,e12.5,3x,e12.5)

c *** transforms nuclear density to momentum space
      write(8,*) '#Numerical potential'
     
      if (zv<1.e-5) then 
         np=2 !neutrons
          write(8,*) '#neutrons'
      else 
         np=1 !protons
          write(8,*) '#protons'
      endif
      call rtoq(denstot,nramax,drx,ncl,n1,n2,np)
      deallocate(densr,u,ffr,bs,denstot,dd)
887   continue
!!$      if (associated(ffp).and.(.not.associated(ffn))) then
!!$         write(*,*) '- Neutron density not read. 
!!$     &              Assuming proton density'
!!$         ffn=>ffp
!!$      else if (associated(ffn).and.(.not.associated(ffp))) then
!!$          write(*,*) '- Proton density not read. 
!!$     &              Assuming neutron density'
!!$         ffp=>ffn
!!$      endif
      return
      end


      subroutine bound(ia,ib,ic,is2,j2,lmom,nod,nramax,jhw,cmass,vmass,
     1wzz,wa,wr0,wls,ffr,vdepth,bengy,dmat,kcheck,wrz,chis,drx,pnloc,
     1wr0ls,wals)
      implicit real*8(a-h,o-z)
      integer:: ia,ib,ic,is2,j2,lmom,nod,nramax,jhw,kcheck
      real*8 :: chis(*),drd(4)
      real*8 :: cmass,vmass ,wzz,wa,wr0,wls,ffr(nramax,4),
     1         vdepth,bengy,dmat,wrz,drx,pnloc,wr0ls,wals
      real*8,allocatable ::wfcr(:),wfsr(:),wfc(:)
      allocate(wfcr(nramax))
      allocate(wfsr(nramax))
      allocate(wfc(nramax))

      korec=0
      do 6000 i=1,nramax
      ffr(i,jhw)=0.0
 6000 continue
      niter=0
      incr=0
      eps7=1.d-6
      drd(jhw)=drx
      lmom1=lmom+1
      test=1.0e+20
      lmom2=lmom1+1
      drz=wrz
      fnod=nod
      flmom=lmom
      flmom1=lmom1
      radi=wr0*(cmass)**0.333333333d0
      radls=wr0ls*(cmass)**0.333333333d0
      yyyy = cosh(radi/wa)
      radz=radi
      wr=drd(jhw)
      nr1=nramax+1
      do 8 i=1,nramax
      yyy=(wr-radi)/wa
      zzz=(wr-radls)/wals
      if(yyy.gt.90.d0) yyy=90.d0
      if(zzz.gt.90.d0) zzz=90.d0
      ex=exp(yyy)
      exls=exp(zzz)
      ib1=ib+1
      select case (ib)
         case (0)
            wfcr(i)=1.0/(1.0+ex)
         case (1)
            wfcr(i)=exp(-(wr/wa)*(wr/wa))
         case (2)
            wfcr(i)=exp(-wa*wr)/wr
         case (3)
            wfcr(i)=exp(-wa*wr)/(exp(-wr0*wr)-exp(-wa*wr))
         case (4)
            wfcr(i) = (1.0+yyyy)/(cosh(wr/wa)+yyyy)
         case default
            write(*,*) 'Shape',ib,'not implemented!'
            stop
      end select
      
      wfsr(i)=exls/(1.0+exls)/(1.0+exls)
      if(wr-radi) 3,3,4
    3 wfc(i)=0.7199262*wzz*(3.0-wr*wr/(radi*radi))/radi
      go to 5
    4 wfc(i)=1.4398523*wzz/wr
    5 wr=wr+drd(jhw)
    8 continue
      if(ia-1) 10,20,20
   10 vdepth=bengy+(3.1415926*(fnod+0.5*flmom1))**2/(0.048228*vmass*
     1 (radi+drz)**2)
      go to 30
   20 bengy=vdepth-(3.1415926*(fnod+0.5*flmom1))**2/(0.048228*vmass
     1*radz*radz)
      if(bengy-eps7) 25,25,30
   25 radz=radz+drz
      incr=incr+1
      if(incr-20) 20,20,27
   27 kcheck=11
      go to 7400
   30 is2mn=abs(j2-2*lmom)
      is2mx=j2+2*lmom
      if(is2.lt.is2mn.or.is2.gt.is2mx) go to 40
      fjs=0.5*real(j2)
      fis=0.5*real(is2)
      flns=fjs*(fjs+1.)-flmom*(flmom+1.)-fis*(fis+1.)
      go to 70
   40 lsq=j2-2*lmom
      if(lsq) 42,44,46
   42 flns=-lmom1
      go to 70
   44 flns=0.
      go to 70
   46 flns=lmom
   70 match=radi/drd(jhw)+dmat
   80 fmu=vmass*cmass/(vmass+cmass)
      wk=0.2195376d0*sqrt(fmu*bengy)
      wrhon=wk*radi
      wrhoc=wrhon
      wrhoz=wk*radz
      wrhocs=wrhoc*wrhoc
      weta=0.7199262*wzz*wk/bengy
      wetac=weta/wrhoc
      wdrho=drd(jhw)*wk
      wvs=2.*wls*wk/(bengy*wals)
      drhosq=wdrho*wdrho
      dr56=0.8333333333*drhosq
      dr12=0.1*dr56
      fl1=lmom*lmom1
  100 wrho=wdrho
      wvc=vdepth/bengy
      zer=1.0
      do 180 j=1,lmom2
      a1=-wvs*flns*wfsr(j)/(flmom1+flmom1)
      b1=1.0-wvc*wfcr(j)+3.0*wetac
      b2=wvs*flns*wfsr(j)
      a2=(b1-b2*a1)/(4.0*flmom1+2.0)
      a3=(b1*a1-b2*a2)/(6.0*flmom1+6.0)
      wrhosq=wrho*wrho
      b3=weta/(wrhoc*wrhocs)
      a4=(b1*a2-b2*a3-b3)/(8.0*flmom1+12.0)
      a5=(b1*a3-b2*a4-b3*a1)/(10.0*flmom1+20.0)
      a6=(b1*a4-b2*a5-b3*a2)/(12.0*flmom1+30.0)
      ffr(j,jhw)=(wrho**lmom1)*(1.0+a1*wrho+a2*wrho*wrho+a3*wrhosq*wrho
     1+a4*wrhosq*wrhosq+a5*wrho*wrhosq*wrhosq+a6*wrhosq**3)
  180 wrho=wrho+wdrho
      mat1=match+1
      x1=wdrho*flmom1
      x2=x1+wdrho
      x3=x2+wdrho
      do 200 i=lmom1,mat1
      fac1=1.0-dr12*(fl1/(x1*x1)+1.0-wvc*wfcr(i)-wvs*flns*wfsr(i)/x1
     1+wfc(i)/bengy)
      fac2=2.0+dr56*(fl1/(x2*x2)+1.0-wvc*wfcr(i+1)-wvs*flns*wfsr(i+1)/x2
     1+wfc(i+1)/bengy)
      fac3=1.0-dr12*(fl1/(x3*x3)+1.0-wvc*wfcr(i+2)-wvs*flns*wfsr(i+2)/x3
     1+wfc(i+2)/bengy)
      ffr(i+2,jhw)=(ffr(i+1,jhw)*fac2-ffr(i,jhw)*fac1)/fac3
      if(dabs(ffr(i+2,jhw))-test) 195,185,185
  185 ip2=i+2
      do 190 ip=1,ip2
      ffr(ip,jhw)=ffr(ip,jhw)/test
  190 continue
      zer=zer/test
  195 x1=x2
      x2=x3
      x3=x3+wdrho
  200 continue
      nodes=0
      do 202 i=lmom1,match
      if(ffr(i,jhw)*ffr(i+1,jhw)) 210,215,202
  210 nodes=nodes+2
      go to 202
  215 nodes=nodes+1
  202 continue
      nn=nodes/2
      fnn=nn
      if(nod-nn) 225,240,225
  225 korec=korec+1
      if(korec-10) 228,228,226
  226 kcheck=10
      go to 7400
  228 vcor=(wrhoz*wrhoz+9.86959*(fnod+0.5*flmom1)**2)/(wrhoz*wrhoz
     1+9.86959*(fnn+0.5*flmom1)**2)
      vcor=sqrt(vcor)
      if(ia-1) 230,235,235
  230 vdepth=vcor*vdepth
      go to 100
  235 bengy=bengy/vcor
      go to 80
  240 dffr1=((ffr(match+3,jhw)-ffr(match-3,jhw))/60.0+3.0*(ffr(match-2,
     1jhw)-ffr(match+2,jhw))/20.0+3.0*(ffr(match+1,jhw)-ffr(match-1,jhw)
     2)/4.0)/wdrho
      korec=0
      rhoa=wk*drd(jhw)*real(nramax)
      wrho=rhoa
      jrho=nramax
  325 ex=wrho+weta*log(wrho+wrho)
      ffr(jrho,jhw+2)=exp(-ex)
      if(jrho-nramax) 340,340,350
  340 jrho=jrho+1
      wrho=wrho+wdrho
      go to 325
  350 x1=rhoa-wdrho
      x2=rhoa
      x3=x2+wdrho
      imax=nramax-match+3
      do 360 i=1,imax
      k=nramax-i
      fac1=1.0-dr12*(fl1/(x1*x1)+1.0-wvc*wfcr(k)-wvs*flns*wfsr(k)/x1
     1+wfc(k)/bengy)
      fac2=2.0+dr56*(fl1/(x2*x2)+1.0-wvc*wfcr(k+1)-wvs*flns*wfsr(k+1)/x2
     1+wfc(k+1)/bengy)
      fac3=1.0-dr12*(fl1/(x3*x3)+1.0-wvc*wfcr(k+2)-wvs*flns*wfsr(k+2)/x3
     1+wfc(k+2)/bengy)
      ffr(k,jhw+2)=(ffr(k+1,jhw+2)*fac2-ffr(k+2,jhw+2)*fac3)/fac1
      if(dabs(ffr(k,jhw+2))-test) 358,352,352
  352 nrten=nramax-k+2
      do 356 iten=1,nrten
      kten=iten+k-1
      ffr(kten,jhw+2)=ffr(kten,jhw+2)/test
  356 continue
  358 x3=x2
      x2=x1
      x1=x1-wdrho
  360 continue
      dffr2=((ffr(match+3,jhw+2)-ffr(match-3,jhw+2))/60.0+3.0*(ffr(match
     1-2,jhw+2)-ffr(match+2,jhw+2))/20.0+3.0*(ffr(match+1,jhw+2)
     2-ffr(match-1,jhw+2))/4.0)/wdrho
      ratio=ffr(match,jhw)/ffr(match,jhw+2)
      tlogd1=dffr1/ffr(match,jhw)
      tlogd2=dffr2/ffr(match,jhw+2)
      difnce=dabs(tlogd1-tlogd2)
      if(difnce-eps7) 510,510,400
  400 niter=niter+1
      if(niter-100) 410,410,405
  405 kcheck=12
      go to 7400
  410 fnum=ffr(match,jhw+2)*dffr2*ratio*ratio-ffr(match,jhw)*dffr1
      sum=0.0
      do 480 i=1,nramax,2
      if(i-1) 420,420,430
  420 sum1=0.0
      go to 445
  430 if(i-match) 440,440,450
  440 sum1=ffr(i-1,jhw)*ffr(i-1,jhw)
  445 sum2=ffr(i,jhw)*ffr(i,jhw)
      sum3=ffr(i+1,jhw)*ffr(i+1,jhw)
      go to 455
  450 sum1=ffr(i-1,jhw+2)*ffr(i-1,jhw+2)*ratio*ratio
      sum2=ffr(i,jhw+2)*ffr(i,jhw+2)*ratio*ratio
      sum3=ffr(i+1,jhw+2)*ffr(i+1,jhw+2)*ratio*ratio
  455 if(ia-1) 460,470,470
  460 if(i-1) 462,462,465
  462 sum1=0.0
      go to 467
  465 sum1=-sum1*wfcr(i-1)*wvc
  467 sum2=-sum2*wfcr(i)*wvc
      sum3=-sum3*wfcr(i+1)*wvc
  470 sum=sum+sum1+4.0*sum2+sum3
  480 continue
      denom=sum*wdrho/3.0
      incr=0
      ram1=fnum/denom
  482 ramda=1.0+ram1
      if(ramda-eps7) 485,485,488
  485 ram1=0.5*ram1
      incr=incr+1
      if(incr-10) 482,482,486
  486 kcheck=13
      go to 7400
  488 if(ia-1) 489,500,500
  489 if(ramda-2.0) 490,495,495
  490 vdepth=ramda*vdepth
      go to 100
  495 drz=drz-0.08*radi
      go to 10
  500 bengy=bengy*ramda
      go to 80
  510 nr1=nramax+1
      do 520 i=match,nr1
      ffr(i,jhw)=ratio*ffr(i,jhw+2)
  520 continue
      sum=0.0
      do 570 i=1,nramax,2
      if(i-1) 540,540,550
  540 sum1=0.0
      go to 560
  550 sum1=ffr(i-1,jhw)*ffr(i-1,jhw)
  560 sum2=ffr(i,jhw)*ffr(i,jhw)
      sum3=ffr(i+1,jhw)*ffr(i+1,jhw)
      sum=sum+sum1+4.0*sum2+sum3
  570 continue
      sum=sum*drd(jhw)/3.0
      znorm=1.0/sqrt(sum)
      wlss=2.*wls*flns/wals
      r=drd(jhw)
      dr=r
c     *****( non-local correction )*****
      ipnl=0
      if(pnloc.lt.0.0) ipnl=dabs(pnloc)/dr+1
      fact=0.048196758*fmu*pnloc**2/8.0
      sum=0.0
      do 599 i=1,nr1
      ffr(i,jhw)=znorm*ffr(i,jhw)
      if(fact.lt.1.e-10) go to 599
      ffr(i,3)=vdepth*wfcr(i)+wlss*wfsr(i)/r-wfc(i)
      ffr(i,jhw)=ffr(i,jhw)*exp(-fact*ffr(i,3))
      if(i.lt.ipnl) ffr(i,jhw)=0.0
      sum=sum+(ffr(i,jhw))**2
  599 r=r+dr
      chis(jhw)=znorm*zer*wk**lmom1
      znorm=1.0
      if(fact.lt.1.e-10) go to 8611
      chis(jhw)=chis(jhw)*exp(-fact*ffr(1,3))
      if(ipnl.gt.0) chis(jhw)=0.0
      znorm=1.0/sqrt(sum*dr)
 8611 r=dr
      do 600 i=1,nr1
      go to (575,580,585),ic
  575 ffr(i,jhw)=znorm*ffr(i,jhw)/r
      ffr(i,3)=vdepth*wfcr(i)+wlss*wfsr(i)/r-wfc(i)      
      go to 600
  580 ffr(i,3)=vdepth*wfcr(i)+wlss*wfsr(i)/r
      go to 590
  585 ffr(i,3)=vdepth*wfcr(i)+wlss*wfsr(i)/r-wfc(i)
  590 ffr(i,jhw)=znorm*ffr(i,jhw)*ffr(i,3)
  600 r=r+dr
      chis(jhw)=znorm*chis(jhw)
      if(ic.ne.1) chis(jhw)=chis(jhw)*ffr(1,3)
      go to 8000
 7400 write (*,7500) kcheck,ramda
      bengy=10.d0
 7500 format(19h0 subroutine bound ,8h kcheck=,i3,7h ramda=,f10.6)

      deallocate(wfcr,wfsr,wfc)
 8000 return
      end



c =====================================================================
c *** Calculates density rho(q) within different models for cluster ncl:
c *** 0.- HO
c *** 1.- WS
c *** 2.- External
c =====================================================================
      subroutine ffact (q2, ff, nff)
c *** see r.h. landau''s program lpott, 1981
c *** calculates the form factor
c *** via special direct formula

c *** can now handle rhop .ne. rhon and spin,
c *** nff = f.f.*s needed, ff(i)= rho p,n,psp,nsp for i=1,2,3,4

      use params
      use switch
      use sizes
      use inputs
      implicit real*8 (a-h, o-z)
      
      integer :: type
      real*8 :: ff(4)

!      common /params/   hbarc, pi, mp, mn, nz, na, nes, nwaves
!      common /sizes/    achp, acmp, wsp, achn, acmn, wsn

      data ndata  /0/

c >>> first executable statement <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


c *** calls to different form factor models
         TYPE=0 !to avoid errors
         select case(type)
            case (0) ! HO
!               write(*,*) 'calling ffhmo'
               call ffhmo (q2, nff, ff(1), ff(2),ncl)
            case (1) !WS
               call wsax(q2, ff, ncl)
            case (3) !Read external
               call ffext(q2,ff, ncl)
               if (q2.gt.800000) ff(1)=0.
         end select

 130  format (29h xxxxbuggy in ffact nff gt2,=,i5)
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



c *** -----------------------------------------------------
c *** Calculates for 3 parameter Fermi density distribution 
c *** normalized to number of particles
c *** -----------------------------------------------------
      subroutine ffactfermi(ncl,n1,n2,nzclus,nnclus)
        use params
        implicit real*8 (a-h,k,o-z)
       integer :: ncl,n1,n2,i,m,kout=16
       integer :: nramax,kin=13, nzclus,nnclus,norm
       real*8  :: w3pf,z3pf,c3pf
       real*8  :: drx,rnorm,densq,yyy,r
       real*8,allocatable,target ::bs(:,:),densr(:), bsr2(:,:)
       real*8, pointer :: bp(:),bpr2(:)

       namelist /bs3pf/ w3p,z3p,c3p,w3n,z3n,c3n,nramax,drx
       
       read(kin,nml=bs3pF)
       write(99,nml=bs3PF)
       allocate(densr(nramax))
       allocate(bs(nramax,2))
       allocate(bsr2(nramax,2))
       

       do 86 i=1,nramax
          r=(i-1.)*drx
          yyy = 1. + exp( (r-c3p)/z3p )
          bs(i,1)= (1. + w3p*(r/c3p)**2)/yyy !rho(r) protons 
          bsr2(i,1)=bs(i,1)*r*r ! rho(r)*r^2
          yyy = 1. + exp( (r-c3n)/z3n )
          bs(i,2)= (1. + w3n*(r/c3n)**2)/yyy !rho(r) neutrons
          bsr2(i,2)=bs(i,2)*r*r ! rho(r)*r^2
 86    continue

       do iso=1,2
          bp=>bs(:,iso)
          bpr2=>bsr2(:,iso)
          
          
          if (iso.eq.1) then
             norm=nzclus
             write(*,80)"protons: ",nzclus
             write(8,*) '#Fermi: protons' 
          else
             norm=nnclus
             write(*,80)"neutrons: ",nnclus
             write(8,*) '&'
              write(14,*) '&'
             write(8,*) '#Fermi: neutrons'
          endif
          write(14,*) '#',nramax,drx
80        format(a8,i3)
          call sim2(bpr2,rnorm,1,nramax,drx,nramax)
c          write(6,*) 'normalization =',rnorm
          bp=bp/rnorm
          densr=dble(norm)*bp
      
       do 95 i=1,nramax
          r=(i-1)*drx
          bp(i)=bp(i)*r**4     !rho*r^4
c          densr(i)= densr(i)/rnorm
          write(14,*) r,densr(i)
!          if (densr(i)<1e-10) densr(i)=0.       
95     continue

          call sim2(bp,rms,1,nramax,drx,nramax)
          write(6,'(" rms=",f6.3)') sqrt(rms)
          

c *** transforms density to momentum space
c        write(8,*) '#3 parameter Fermi density'
       call rtoq(densr,nramax,drx,ncl,n1,n2,iso)
       enddo ! end loop p/n
       deallocate(densr,bs,bsr2)
       return
       end 



c *** -----------------------------------------------------
c *** Calculates for 3 parameter Fermi density distribution 
c *** normalized to number of particles
c *** -----------------------------------------------------
      subroutine ffactfermi2(ncl,n1,n2,nzclus,nnclus)
        use params
        implicit real*8 (a-h,k,o-z)
       integer :: ncl,n1,n2,i,m,kout=16
       integer :: nramax,kin=13, nzclus,nnclus,norm
       real*8  :: w3pf,z3pf,c3pf
       real*8  :: drx,rnorm,densq,yyy,r
       real*8,allocatable,target ::bs(:,:),densr(:), bsr2(:,:)
       real*8, pointer :: bp(:),bpr2(:)

       namelist /bs3pf/ w3p,z3p,c3p,w3n,z3n,c3n,nramax,drx
       
       read(kin,nml=bs3pF)
       write(*,nml=bs3PF)
       allocate(densr(nramax))
       allocate(bs(nramax,2))
       allocate(bsr2(nramax,2))
       

       
       do 86 i=1,nramax
          r=(i-1.)*drx
          yyy = 1. + exp( (r-c3p)/z3p )
          bs(i,1)= (1. + w3p*(r/c3p)**2)/yyy !rho(r) protons 
          bsr2(i,1)=bs(i,1)*r*r ! rho(r)*r^2
          yyy = 1. + exp( (r-c3n)/z3n )
          bs(i,2)= (1. + w3n*(r/c3n)**2)/yyy !rho(r) neutrons
          bsr2(i,2)=bs(i,2)*r*r ! rho(r)*r^2
 86    continue

       do iso=1,2
          bp=>bs(:,iso)
          bpr2=>bsr2(:,iso)
                    
          if (iso.eq.1) then
             norm=nzclus
             write(*,80)"protons",nzclus
             write(8,*) '#Fermi: protons' 
          else
             norm=nnclus
             write(*,80)"neutrons",nnclus
             write(8,*) '&'
              write(14,*) '&'
             write(8,*) '#Fermi: neutrons'
          endif
          write(14,*) '#',nramax,drx
80        format(a8,i3)
          call sim2(bpr2,rnorm,1,nramax,drx,nramax)
c          write(6,*) 'normalization =',rnorm
          bp=bp/rnorm
          densr=dble(norm)*bp
      
       do 95 i=1,nramax
          r=(i-1)*drx
          bp(i)=bp(i)*r**4     !rho*r^4
c          densr(i)= densr(i)/rnorm
          write(14,*) r,densr(i)
          
!          if (densr(i)<1e-10) densr(i)=0.       
95     continue

          call sim2(bp,rms,1,nramax,drx,nramax)
          write(6,'(" rms=",f6.3)') sqrt(rms)
          

c *** transforms density to momentum space
c        write(8,*) '#3 parameter Fermi density'
       call rtoq(densr,nramax,drx,ncl,n1,n2,iso)
       enddo ! end loop p/n
       deallocate(densr,bs,bsr2)
       return
       end 



c *** -----------------------------------------------------
c *** Calculates for G3 density distribution 
c *** normalized to number of particles
c *** -----------------------------------------------------
      subroutine ffactg3(ncl,n1,n2,nzclus,nnclus)
        use params
        implicit real*8 (a-h,k,o-z)
       integer :: ncl,n1,n2,i,m,kout=16
       integer :: nramax,kin=13, nzclus,nnclus,norm
       real*8  :: w3pf,z3pf,c3pf
       real*8  :: drx,rnorm,densq,yyy,r
       real*8,allocatable,target ::bs(:,:),densr(:), bsr2(:,:)
       real*8, pointer :: bp(:),bpr2(:)

       namelist /bs3pf/ w3p,z3p,c3p,w3n,z3n,c3n,nramax,drx
       
       read(kin,nml=bs3pF)
       write(kout,nml=bs3PF)
       allocate(densr(nramax))
       allocate(bs(nramax,2))
       allocate(bsr2(nramax,2))
       

       do 86 i=1,nramax
          r=(i-1.)*drx
          x=(r-c3p)/z3p
          yyy = 1. + exp(x**2)
          bs(i,1)= (1. + w3p*(r/c3p)**2)/yyy !rho(r) protons 
          bsr2(i,1)=bs(i,1)*r*r ! rho(r)*r^2
          x= (r-c3n)/z3n 
          yyy = 1. + exp(x**2)
          bs(i,2)= (1. + w3n*(r/c3n)**2)/yyy !rho(r) neutrons
          bsr2(i,2)=bs(i,2)*r*r ! rho(r)*r^2
 86    continue

       do iso=1,2
          bp=>bs(:,iso)
          bpr2=>bsr2(:,iso)
          
         
          if (iso.eq.1) then
             norm=nzclus
             write(*,'(" protons=",i3)') nzclus
             write(8,*) '#G3: protons' 
          else
             norm=nnclus
             write(*,'(" neutrons=",i3)') nnclus
             write(8,*) '&'
             write(14,*) '&'
             write(8,*) '#G3: neutrons'
          endif
           write(14,*) '#',nramax,drx
          call sim2(bpr2,rnorm,1,nramax,drx,nramax)
c          write(6,*) 'normalization =',rnorm
          bp=bp/rnorm
          densr=dble(norm)*bp
      
       do 95 i=1,nramax
          r=(i-1)*drx
          bp(i)=bp(i)*r**4     !rho*r^4
c          densr(i)= densr(i)/rnorm
          write(14,*) r,densr(i)
!          if (densr(i)<1e-10) densr(i)=0.       
95     continue

          call sim2(bp,rms,1,nramax,drx,nramax)
          write(6,'(" rms=",f6.3)') sqrt(rms)
          

c *** transforms density to momentum space
c        write(8,*) '#G3 density'
       call rtoq(densr,nramax,drx,ncl,n1,n2,iso)
       enddo ! end loop p/n
       deallocate(densr,bs,bsr2)
       return
       end 




c *** --------------------------------------------------------
c *** Calculates density in momentum space for HO distribution
c *** (normalized to number of nucleons)
c *** --------------------------------------------------------
      subroutine ffactho(ncl,n1,n2)
        use rhos
        use params
        use bgrid
!        implicit real*8 (a-h,k,o-z)
        implicit none
        integer ncos
        parameter (ncos=48)
        integer :: xi
        integer::ncl,n1,n2,nshell,ibs,isos,ir,i2,m,iq,i1
        integer :: kin=13
        integer :: code,norba,zorba
        real*8  ::densq,rho,norm,drx,densr,theta,cthnuc
        real*8 :: rmin,rmax,r,aap,aan !zorba,norba
        real*8,allocatable :: aa1(:,:)
        integer,allocatable :: norba1(:,:)
        real*8 :: xis(ncos), ws(ncos)
        real*8:: q2,k,kp,yyy,qmin,qmax,dq,aa,x,a,b,q
        character*2 shname(4)
       
        
        namelist /hoparm/ nshell
        namelist /hoshell/aap,norba,aan,zorba
        shname(1)='1s'
        shname(2)='1p'
        shname(3)='1d'
        shname(4)='2s'
        
     
        read(kin,nml=hoparm)
        write(99,nml=hoparm)
        allocate(aa1(nshell,2))
        allocate(norba1(nshell,2))

!     Loop from ** inner to outer ** shells
        write(*,'(3a10)')' Shell ',' Protons ',' Neutrons '
        do ibs=1,nshell
           norba=0; zorba=0.
           read(kin,nml=hoshell)
           write(*,'(4x,1a2,2i8)')shname(ibs),zorba,norba
           write(99,nml=hoshell)
           aa1(ibs,1)=aap
           aa1(ibs,2)=aan
           norba1(ibs,1)=zorba
           norba1(ibs,2)=norba
        enddo


c *** Density in configuration space
c *** (not used, just for check)
      rmin=0.00
      rmax=5.
      drx=0.1
      m=(rmax-rmin)/drx+1
      
      do isos=1,2 ! loop in p/n
         if (isos.eq.2) write(14,*)'&' 
         write(14,*) '#',m,drx
      do ir=1,m
         r = rmin + (ir-1)*drx
         densr = 0.d0
         do ibs=1,nshell
            a=aa1(ibs,isos)
            yyy = exp(-r*r/a/a)
            select case(ibs)
            case(1)  !    1s shell
               norm=4./dsqrt(pi*a**6)
               rho=1.
            case(2) !     1p shell
               norm=4.*dsqrt(4./9./pi/a**10)
               rho=r**2
            case(3) !     1d shell
               norm=16./dsqrt(9./25./pi/a**14)
               rho=r**4
            case(4) !     2s shell
               norm=dsqrt(2**6/9./pi/a**6)
               rho= (1.5-r**2/a/a)**2
            case default ! ibs>4
               write(*,*)'HO only up to 2s shell implemented'
               STOP
            end select
            rho=rho*norba1(ibs,isos)*norm * yyy
            densr=densr+rho
!            write(*,*) isos,shname(ibs),rho
         enddo !end loop shells
         
         write(14,*) r,densr
      enddo !end loop p/n
      enddo

!     Grid and loop in q for printing
        qmin=0.05
        qmax=5.
        dq=0.1
        m=(qmax-qmin)/dq+1
        write(8,*) '#HO densities for protons/neutrons'
       
        do 950 isos=1,2 !p/n
           if (isos.eq.1) then 
              write(8,*) '#protons'
           else
              write(8,*) '&'
              write(8,*)'#neutrons'
           end if
              
        do 900 iq=1,m
           q = qmin + (iq-1)*dq
           densq = 0.d0
           rho=0.
           do 951 ibs=1,nshell
              aa=aa1(ibs,isos)
              yyy = exp(-q*q*aa*aa/4.)
              select case(ibs)
              case(4) !     2s shell
                 rho=yyy/6.*(6.-2.*(q*aa)**2+(q*aa)**4/4.)
              case(3) !     1d shell
                 rho = yyy/5.*(5 -10/6.*(q*aa)**2+
     &                 (q*aa)**4/12.) 
              case(2) !     1p shell
                 rho = yyy*2./3.*(3/2. - 1/4.*(q*aa)**2)
              case(1)  !    1s shell
                 rho = yyy 
              case default ! ibs>4
                 write(*,*)'HO only up to 2s shell implemented'
                 STOP
              end select
              rho=rho*norba1(ibs,isos)
              densq=densq+rho
 951  continue
            write(8,*)q,densq    
 900  continue
 950   continue ! p/n
c      write(8,*) '#---------------------------------------------'
      write(8,*) '&'

      
!     Grid and loop in q for MSO and storage
      a = -1.
      b = 1.
      code = 11
      call gauss2( ncos, code, a, b, xis, ws )

!     Grid and loop in q 
      do 240 i1 = 1,n1
         do 240 i2 = 1,i1
            k = kk(i2)
            kp = kk(i1)
            do  xi = 1, ncos
               x = xis( xi )
               cthnuc = x
               theta = acos( cthnuc ) * 180./pi
               q2 = kp**2 + k**2 - 2.* kp * k * cthnuc
               q2=q2/hbarc/hbarc !convert to natural units
               if ( q2 .lt. 0 ) q2 = -q2
               q=dsqrt(q2)              
               do isos=1,2
                  densq = 0.d0
               do 952 ibs=1,nshell
                  aa=aa1(ibs,isos) 
                  yyy = exp(-q*q*aa*aa/4.)
                  select case(ibs)
                  case(4) !     2s shell
                     rho = yyy/6.*(6-2*(q*aa)**2+
     &                       (q*aa)**4/4.)
                  case(3) !     1d shell
                     rho = yyy/5.*(5 -10/6.*(q*aa)**2+
     &                       (q*aa)**4/12.) 
                  case(2) !     1p shell
                     rho = yyy*2./3.*(3/2. - 1/4.*(q*aa)**2)
                  case(1)  !    1s shell
                     rho = yyy 
                  case default ! ibs>4
                     write(*,*)'HO only up to 2s shell implemented'
                     STOP
                  end select
                  rho=norba1(ibs,isos)*rho
                  densq=densq+rho
952              continue ! end loop ibs
!                 if ((i1.eq.i2).and.(isos.eq.1).and.q<20) 
!     &            write(98,*) q,densq
!                 if ((i1.eq.i2).and.(isos.eq.2).and.q<20) 
!     &             write(97,*) q,densq
                  ff(ncl,i1,i2,xi,isos)=densq
                 enddo ! end loop p/n
                enddo ! end loop xi               
240             continue ! end loop k,kp
                ffp=>ff(:,:,:,:,1) !protons
                ffn=>ff(:,:,:,:,2) !neutrons
                deallocate(aa1,norba1)
      return
      end

c-----------------------------------------------------------------
c *** Reads from external file radial matter density and transforms
c *** to momentum space
c-----------------------------------------------------------------
      subroutine ffactextr(ncl,n1,n2)
        implicit real*8 (a-h,k,o-z)
        integer :: xi,np
        integer ::norba,kin=13
        integer :: code,ncos,nramax
        real*8, allocatable :: densr(:)

c *** protons
        read(4,*) nramax,drx
        allocate(densr(nramax))
        do i=1,nramax
           read(4,*)  xr,densr(i)
        end do

c *** transforms nuclear density to momentum space
        write(8,*) '#External density'
        write(8,*) '#protons'
        np=1
        call rtoq(densr,nramax,drx,ncl,n1,n2,np)

c *** neutrons
        read(4,*) nramax,drx
        do i=1,nramax
           read(4,*)  xr,densr(i)
        end do

        write(8,*) '#neutrons'
        np=2
        call rtoq(densr,nramax,drx,ncl,n1,n2,np)
        deallocate(densr)
        return
        end

c--------------------------------------------------------------
c *** fourier transforms the radial density to momentum space
c--------------------------------------------------------------
      subroutine rtoq(densr,nramax,drx,ncl,n1,n2,np)
        use rhos
        use params
        use bgrid 
!        implicit real*8 (a-h,k,o-z)
        implicit none
        integer:: xi,m,i1,i2,ncl,n1,n2,np,i,iq
        integer:: kin=13
        integer:: code,ncos,nramax
        real*8 :: densq,rho,x,rnorm,drx,r,qr,a,b
        real*8 :: k,kp,q2,q,cthnuc,theta
        real*8 :: qmin,qmax,dq
        real*8, allocatable :: dd(:)
        real*8, intent(in) :: densr(nramax) 
        parameter (ncos=48)
        real*8 :: xis(ncos), ws(ncos) 
        allocate(dd(nramax))

        

c *** check normalization of density
        do i=1,nramax
           r=(i-1)*drx
           dd(i)=densr(i)*r*r
        end do
        call sim2(dd,rnorm,1,nramax,drx,nramax)
        write(*,'(5x," density normalized to: ",1f6.2)') rnorm

c ***  Grid and loop in q for printing
        qmin=0.0001
        qmax=5.
        dq=0.1
        m=(qmax-qmin)/dq+1
        do 900 iq=1,m
          q = qmin + (iq-1)*dq
          do 602 i=1,nramax
             r=(i-1)*drx
             qr=q*r
!!!!             write(*,*)i,densr(i)
             if(qr>1.d-10) dd(i)=r*r*densr(i)*sin(qr)/qr
  602     continue
          call sim2(dd,densq,1,nramax,drx,nramax) 
          write(8,*)q,densq
  900  continue


c ***  Grid and loop in q for MSO and storage
       a = -1.
       b = 1.
       code = 11 
       call gauss2( ncos, code, a, b, xis, ws )
       do 240 i1 = 1,n1
         do 240 i2 = 1,i1
            k = kk(i2)
            kp = kk(i1)
            do  xi = 1, ncos
               x = xis( xi )
               cthnuc = x
               theta = acos( cthnuc ) * 180./pi
               q2 = kp**2 + k**2 - 2.* kp * k * cthnuc
               q2=q2/hbarc/hbarc !convert to natural units
               if ( q2 <0 ) q2 = -q2
               q=dsqrt(q2)              
               densq = 0.d0
               do i=1,nramax
                  r=(i-1)*drx
                  qr=q*r
                  if(qr > 1.d-10) dd(i)=r*r*densr(i)*sin(qr)/qr
               enddo
! This should be improved!!!!!!!!????????????
               if (qr<100.) call sim2(dd,densq,1,nramax,drx,nramax) 
               ff(ncl,i1,i2,xi,np)=densq

cc  TEST: print on-shell values at quadrature points
!                if ((i1.eq.i2).and.(np.eq.1).and.q<20) 
!     &            write(98,*) q,densq
!                 if ((i1.eq.i2).and.(np.eq.2).and.q<20) 
!     &             write(97,*) q,densq


            enddo ! end loop xi
240         continue ! end loop k,kp
            if (np.eq.1) then
               ffp=>ff(:,:,:,:,1)
            else
               ffn=>ff(:,:,:,:,2)
            endif
          deallocate(dd) 
        return
        end
        
        
c-----------------------------------------------------------------
c *** Reads from external file radial matter density 
c *** in MOMENTUM space
c-----------------------------------------------------------------
      subroutine ffactextq(ncl,n1,n2)
        use rhos
        use params
        use bgrid 
        implicit real*8 (a-h,k,o-z)
        integer :: xi,np
        integer ::norba,kin=13
        integer :: code,ncos,nqp,nqn
        real*8, allocatable :: qp(:),qn(:)
        complex*16,allocatable::rhoqp(:),rhoqn(:)
        complex*16 wxxi,zaux
        real*8::raux,cero=0.

        parameter (ncos=48)
        dimension xis(ncos),wts(ncos)

c *** protons
        read(8,*) nqp
        write(*,'("- Protons: Reading",1i4," points")') nqp
        allocate(rhoqp(nqp))
        allocate(qp(nqp))
        do i=1,nqp
           read(8,*,err=800) qp(i),raux
           rhoqp(i)=cmplx(raux,cero)
!           write(*,*)  qp(i),rhoqp(i)
        end do

        
c *** neutrons
        read(8,*) nqn
        write(*,'("- Neutrons: Reading",1i4," points")') nqn
        allocate(rhoqn(nqn))
        allocate(qn(nqn))
        do i=1,nqn
           read(8,*,err=800) qn(i),raux
           rhoqn(i)=cmplx(raux,cero)
        end do
        
!        write(*,*) cmplx(rhoqn,0)

c ***  Grid and loop in q for MSO and storage
       a = -1.
       b = 1.
       code = 11 
       call gauss2( ncos, code, a, b, xis, wts )
       do 240 i1 = 1,n1
          do 240 i2 = 1,i1
             k = kk(i2)
             kp = kk(i1)
             do  xi = 1, ncos
                x = xis( xi )
                cthnuc = x
                theta = acos( cthnuc ) * 180./pi
                q2 = kp**2 + k**2 - 2.* kp * k * cthnuc
                q2=q2/hbarc/hbarc !convert to natural units
                if ( q2 <0 ) q2 = -q2
                q=dsqrt(q2)
                zaux=cmplx(0,0)
                if (q<qp(nqp)) zaux=wxxi(q,qp,rhoqp,nqp,nqp)
                ff(ncl,i1,i2,xi,1)=real(zaux)
            
                zaux=cmplx(0,0)
                if (q<qn(nqn)) zaux=wxxi(q,qn,rhoqn,nqn,nqn)
                ff(ncl,i1,i2,xi,2)=real(zaux)
!                if ((k.eq.kp).and.abs(zaux)>1) then
!                   write(*,*)k,kp,q, ff(ncl,i1,i2,xi,2)
!                endif
               
             enddo ! end loop xi
240     continue ! end loop k,kp
               
        ffp=>ff(:,:,:,:,1)
        ffn=>ff(:,:,:,:,2)

        deallocate(rhoqp,rhoqn,qp,qn)
        return
800     write(*,*)'Error reading external densities.Aborting';stop
        end

        
