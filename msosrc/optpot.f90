      subroutine optp ( n1, lmaxin)
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
c 
c    Main changes
c    - common blocks replaced by modules
c    - loops on partial waves start on l=0 (not l=1) 
c    - vll changed from real*8 to complex*16 (1 index dropped)
c    - namelist input
c ***************************************************************

c *** need to reset dimensions and equivalences as max # grid increase
c     should use n = ngp + 1 for scattering.
c     ff(ndimf/6, 6) also needed ,and equivalenced to f(1)
c     for equivalence need length + 1 to overlay
c     also see initialize statements which now change with dimensions

c *** for spin 0 * 0.5 calculate the partial wave decomposition
c     of rho (the form factor).
c     rho(nff,i) is the i-1 th partial wave amplitude of the form factor
c     for a value of k and kp.

c *** nff determines which f.f. is used
c        nff = 1 proton matter     = 2 neutron matter
c            = 3 proton spin       = 4 neutron spin
c *** revisions ***
c
c        ms (12/4/87)
c           changed nuc0 to a logical type variable called spin0
c
c ***************************************************************
      use bandt
      use params
      use switch
      use nlspfl
      use foptp
      use sizes
      use ranges
      use blocks
      use bgrid
      use bprop
      use clebma
      use qnumb
      use bkcall
      use inputs
      use rhos

      implicit real*8 (a-h, k, o-z)

      logical :: spin0
      integer :: xi,code
      real*8  :: k, kp, k1
      real*8, allocatable :: vcoull(:)
      real*8, allocatable :: vcpoint(:)
      complex*16 :: v2(2)
      complex*16, allocatable :: vll(:,:)
      complex*16, allocatable :: v2l(:,:)
      complex*16, allocatable :: v2l2(:,:)
      real*8, pointer :: kpts(:)
      dimension vr(6), vi(6)
      real*8::vc
      
      vc=0d0
     
!      data  nota        /1/
!      data  nint        /0/
      data  maxls/ 100 /
      data (v2(j),j=1,2)/2*0.d0/

c >>> first executable statement <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      kpts=>kk

      nz=zt
      na=at
    
!      if (lpia .ne. 1) go to 250

c *** 1st call to optp, calculate all u's
c *** initialize (can't use data for this due to commoms)

      comp2=0.
!      nint = 1
      k = -200.
      spin0 = .false.

c *** special 0 * 0.5 spin case, use old t * rholl (the partial wave
c *** projection of the form factor).

!      if (nifty(6) .eq. 3 .or. nifty(6) .eq. 0) then
         spin0 = .true.
!      endif
      if (kpts(n1) .gt. 0.) then
         k1 = kpts(n1)
      else
         k1 = kpts(1)
      endif
      write(kout,410) k1
      n2 = 2 * n1
      nn = na - nz

      rhl3 = (achp/hbarc)**2/2.
      rhl4 = (achn/hbarc)**2/2.

c *** added terms for c-13 spin...may want to change size here

      rhl5 = (wsn/hbarc)**2/2.
      as = 1./6.

c *** do loop over all grid points
!      if (nifty(6) .eq. 6) nifty(17) = 1
      aovera = (na - 1.)/(na*1.)
!      if (nifty(17) .eq. 1) aovera = 1.
      if (na .le. 2) aovera = 1.

!      if (nifty(6) .eq. 4)  aovera = 0.
      xhl = aovera * 2. * pi * pi
      aparam = 1.
!      if (nifty(1) .eq. 9)  aparam = 0.
       if (zp .eq. 1)  aparam = 0. !neutrons
      if (rcoul .lt. 1.e-10 ) rcoul = 1.e-6
      aparam = (1.5 * nz * aparam/rcoul) * hbarc * 7.29735e-3
      xhl1 = rcoul/hbarc
      xhl2 = rcut/hbarc
      nif = coulomb

      kcall = 0
      write(99,*) '+ Allocating ul  with',lmaxin,n1,n1
      allocate(ul(0:lmaxin,n1,n1,2))
      allocate(ucoul(0:lmaxin,n1,n1,2))
      ul=0d0;ucoul=0d0;

c *** formfactors (densities in momentum space)
      call formfactors(n1,n2)
      
      do 240 i1 = 1,n1
         do 240 i2 = 1,i1
         kcall = kcall + 1
         k = kpts(i2)
         kp = kpts(i1)
         lmax = lmaxin + 1
         
         if(coulomb>2)then
             allocate(vcpoint(0:lmax))
             allocate(vcoull(0:lmax))
             call vcoul(k,kp,vcpoint,lmax,0)
             call vcoul(k,kp,vcoull,lmax,coulomb)
         end if

c *** new version....outer do loop over kp, k; inner over all l''s

         xmn = mn

c ***
c *** vll(neles, nspin-mp) 
c ***
         lmas2 = lmaxin + 2
         allocate(vll(0:lmas2,2))
         allocate(v2l(0:lmas2,2))
         allocate(v2l2(0:lmas2,2))

            
c **** first order calculations
c ***  vopt0a should contain all possible form factors
          call vopt0a( vll, k, kp, lmas2, k1,i1,i2)

c *** second order calculations
         select case(second)
         case(2)
            write(99,*) 'Entering vopt2l for second order calc.'
            call vopt2l(v2l,k,kp,lmas2,k1)
         case(3,4)
            if (closure)then
               comp2 = 0.
         call vopt2nlc(v2l,k,kp,i2,i1,n1,lmas2)
c ***                 check
c                     call vopt2nl(v2l,k,kp,i2,i1,n1,lmas2)
c                     call vnl2(v2l2,k,kp,i2,i1,n1,lmas2)
c                     comp2=1.
c ***                 end of the check
                   else
                      comp2 = 1.
                      call vopt2nl(v2l,k,kp,i2,i1,n1,lmas2)
                      call vnl2(v2l2,k,kp,i2,i1,n1,lmas2)
                   endif
                end select

c *** do loop over l ***

            do 230 ldum = 0,lmaxin
               j1 = ldum
               urc = 0.
               vc = 0.
               uic = 0.
               urs = 0.
               uis = 0.

               if (coulomb.ge.3)then
!                  write(*,*) 'test1'
                  vc = vcoull(ldum)
c *** ucoul(ldum,k,kp,spin)
!! Can we drop the argument on spin????
                  ucoul(ldum,i1,i2,1)=vcpoint(ldum)
                  ucoul(ldum,i1,i2,2)=vcpoint(ldum)
!                   write(*,*) 'test2'
               end if

c *** set up potential needed for this l in integral equation
c *** n.b. the f is + imv as stored.

               npot1 = 2 * ldum

               if(second.gt.1)then
                  v2(1) = v2l(ldum,1) + comp2*v2l2(ldum,1)
                  v2(2) = v2l(ldum,2) + comp2*v2l2(ldum,2)
               else
                  v2=0d0
               endif

c *** this is for singlet


c *** url(ldum,k,kp,spin)
               ul(ldum,i1,i2,1)=vll(ldum, 1)+vc+v2(1)
               ul(ldum,i1,i2,2)=vll(ldum, 2)+vc+v2(2)
               

               if (   kp  .eq.  k1
     $                    .and.
     $                kp  .eq.  k
     $                    .and.
     $              ldum  .eq.  1) write(kout,410) k1

c *** print check for kp=k1 and ldum=0,1,2
               if (kp .ne. k1 .or. ldum .gt.2) go to 330
               write(kout,400) k, real(ul(ldum,i1,i2,1)), 
     &                   vc,dimag(ul(ldum,i1,i2,1))


c *** set up potential with u(kp,k) = u(k,kp)..symmetric u

  330       if (i1 .eq. i2) go to 230
            ul(ldum,i2,i1,1) = ul(ldum,i1,i2,1)
            ul(ldum,i2,i1,2) = ul(ldum,i1,i2,2)

            if (coulomb.ge.3) then
!               write(*,*)'test3'
               ucoul(ldum,i2,i1,1)= ucoul(ldum,i1,i2,1)
               ucoul(ldum,i2,i1,2)= ucoul(ldum,i1,i2,2)
!                write(*,*)'test4'
            endif
              

c *** write out potential matrices on printer (half off-shell)

c *** do loop over ldum ends
 230         continue

c *** vll no longer used; can release memory
            if (allocated(vll)) deallocate(vll)
            if (allocated(v2l)) deallocate(v2l)
            if (allocated(v2l2)) deallocate(v2l2)
            if (allocated(vcpoint)) deallocate(vcpoint)
            
            if (allocated(vcoull)) deallocate(vcoull)
!             write(*,*)'Deallocate vll,v2l...'
            if ((i1 .eq. i2) .and. (i1 .eq. n1)) npot1m = 2 * lmaxin

c *** do loop over kp and k ends

 240  continue
 340  continue


      deallocate(ff)
      return


 350  format(' c13 neutron spin ff for k = kp = ',e15.5/10(10x,6e15.6/))
 400  format (3x, 3e19.4, 22x, e19.4)
 401  format(1x, i2, 19x, e19.4, 41x, e19.4)
 410  format('0','optical potenial (optp) for k = k0 =',
     $            f13.3,' mev  and l = 0',
     $           /, 16x, 'kp', 17x, 'rev', 16x, 'vc', 38x, '-imv', /)
 999   format(6e12.4)
      return
      end



c-----------------------------------------------------------------------
c     F matrix
c-----------------------------------------------------------------------
      subroutine fmatrx( den , n1 , nmax , nspin , nptmax , ldum )
      use nlspfl
      use foptp
!      use switch
      use spins

      implicit real*8 ( a-h, o-z )
      integer :: kout=16, fwrite,ldum,nptmax,nmax,n1
      real*8 :: den(:) !den(nmax)

c *** assign pointer u to array ul for this ldum
      u=>ul(ldum,:,:,:)
      uc=>ucoul(ldum,:,:,:)
88    format(1i3,6g12.4)

      do 240 i2 = 1,nmax
         do 240 i3 = 1,nmax

c *** regular case  j=l or j=l +/- 1/2
      if (nspine.gt.2)  then
         write(*,*) 'nspine>2 not implemented. Aborting!!!'
         stop
      endif

c *** nuclear + coulomb
         fr(i2,i3) = real(u(i2,i3,nspine)) * den(i3)
         fi(i2,i3) = dimag(u(i2,i3,nspine)) * den(i3)

c *** coulomb
         fcr(i2,i3) = uc(i2,i3,nspine) * den(i3)
         fci(i2,i3) = 0.
c          write(kout,*) 'fr: ',i2,i3,fr(i2,i3),fi(i2,i3)
c           write(kout,*) 'frc: ',i2,i3,fcr(i2,i3)

c *** diagonal term
         if (i2 .eq. i3) then 
            fr(i2,i3) = fr(i2,i3) + 1.
            fcr(i2,i3) = fcr(i2,i3) + 1.
         endif
            
 240  continue

c      if (ldum.eq.1)then
c      do 235 ii=1,nmax
c        do 235 jj=1,nmax
c        write(kout,*)'(fr,fi)=',fr(ii,jj),fi(ii,jj)
c235    continue
c         endif

c *** write out at least part of f matrix

      fwrite=0
      if (fwrite.eq.0) goto 208

      nn = nmax
      if (ldum .ne. 1) go to 208
c     if( nifty(20) .eq. 1 .or. nifty(20) .eq. 3 ) then
        write(kout,928)
c     endif
      i = 1
      if (nn .gt. 10) i = 6
      if (nn .gt. 32) i = nn/6

c     if( nifty(20) .eq. 1 .or. nifty(20) .eq. 3 ) then
        do 217 ii = 1,nn,i
217        write(kout,227) ii, (fr(ii,ij), ij=1,nn)
c     endif

      ii = nn
c     if( nifty(20) .eq. 1 .or. nifty(20) .eq. 3 ) then
        write(kout,227) ii, (fr(ii,ij), ij=1,nn)
c     endif

c     if( nifty(20) .eq. 1 .or. nifty(20) .eq. 3 ) then
c        write(kout,928)

        do 218 ii = 1,nn,i
218        write(kout,227) ii, (fi(ii,ij), ij=1,nn)
c     endif

      ii = nn
c     if( nifty(20) .eq. 1 .or. nifty(20) .eq. 3 ) then
c        write(kout,227) ii, (fi(ii,ij), ij=1,nn)
c     endif
 208  continue
      nfdim = nptmax * 2

c *** cmatin does matrix inversion for complex matrix
      call cmatin( fr , fi , nmax , nfdim )
      call cmatin( fcr , fci , nmax , nfdim )
      return

c 250  nms = 1 + nifty(6)
  250  nms = 1 + 3
      write(kout,820) nms

      return

c *** formats ***

 227  format(i3, 8e15.6 / 200(3x, 8e15.6 / ))
 820  format (' r matrix calculated with ms series up to order ',i3)
 928  format('0',' re/im f matrix before inversion ')

      end

c-------------------------------------------------------------------
c     R-matrix
c-------------------------------------------------------------------
      subroutine rmatrx ( psfac , n1, n2 , nmax, nspin )
      use nlspfl
      use foptp
      use rcomn
      use spins

      implicit real*8 (a-h, o-z)
      integer ::kout=16
      rr = 0.   ; rcr=0.
      ri = 0.   ; rci=0.
      rrb = 0.  ; rcrb=0.
      rib = 0.  ; rcib=0.

c *** calculation  of  r-matrix   (on  shell)
c ***
c *** do matrix mult   r= finv * u      f always square r not necly so
c ***

      do 270 j3 = 1,nmax

c *** nspin =1 or 2    i.e. spin 0x0, 0 x 1/2, .5 x .5 (singlet+j=l)
c *** uncoupled channels
c ***

         if (j3 .gt. n1) write(kout,242) j3
         rr = rr + fr(n1,j3) * real(u(j3,n1,nspine))
     $           - fi(n1,j3) * dimag(u(j3,n1,nspine))
         ri = ri + fr(n1,j3) * dimag(u(j3,n1,nspine))
     $           + fi(n1,j3) * real(u(j3,n1,nspine))

c *** coulomb 
         rcr = rcr + fcr(n1,j3) * uc(j3,n1,nspine)
         rci = rci + fci(n1,j3) * uc(j3,n1,nspine)


 270  continue

c       write(kout,970) rr, ri
c      write(kout,972) rcr, rci
      
      rr = -psfac * rr
      ri = -psfac * ri
      rrb = -psfac * rrb
      rib = -psfac * rib
c      write(15,*) ldum-1,nspine,rr,ri

c *** point coulomb (screened)
      rcr = -psfac * rcr
      rci = -psfac * rci
      rcrb = -psfac * rcrb
      rcib = -psfac * rcib
      
      write(kout,2423) rr, ri, rrb, rib
      write(kout,2423) rcr, rci, rcrb, rcib
      return

c *** format statements ***

 242  format(' ********** j3 too big in rcalc at 242 ************',i8/)
 970  format(1h ,37h  unnormalized r matrix(re,im)       ,13x,4e14.7)
 972  format(1h ,37h  unnormalized rc matrix(re,im)       ,13x,4e14.7)
 2423 format(1h ,39h normalized re and im ra,rb(ko,ko,ko)=  ,
     $ 11x,4e14.7)

      end

c====================================================================
c *** calculation  of    t-matrix  (on - shell),
c *** t normalized to exp(i del)*sin(del)
c===================================================================
      subroutine tmatrx( sigl, aovera, ld, nspin, njump )
      use kinemt
      use params
!      use switch
      use sizes
      use ranges
      use inputs
      use rcomn
      use spins
      use tcoul
      use tcomn
      use tcoulb
!      use wr
      use bscoul
      use etapbl  
      implicit real*8 (a-h, o-z)

      real*4 :: jspin
      complex*16:: zi , ztan
      complex*16:: zr34 , zr56 , zrmix , ztjp, ztjm , ztj
      complex*16:: tcpl, tcplc
      complex*16:: tfact,smatrx
      n = nspin
      zi = ( 0. , 1. )

c *** nspin=1,2  usual, old case
c *** using born for high l, this  unitarizes,therefore differ from ss

c *** coulomb + nuclear
 280  cnst1 = ri * ri
      cnst2 = rr * rr
      cnst3 = 1. + cnst1 - ri - ri
      cnst4 = cnst3 + 4. * ri

      aovera = at/(at - 1.)

      if (.not.kmt) then 
         avera=1.
         write(*,*) 'Watson->aovera=1'
      endif
      if (at .eq. 2) aovera = 1.

      tr(ld,n) = aovera*rr/(cnst4 + cnst2) 
      ti(ld,n) = aovera*(ri + cnst1 + cnst2)/(cnst4 + cnst2) 

      trc = tr(ld,n)
      tic = ti(ld,n)      
      tcpl=cmplx(trc,tic)
      

c *** coulomb      
      if (coulomb.ge.3)then
         cnst1 = rci * rci
         cnst2 = rcr * rcr
         cnst3 = 1. + cnst1 - rci - rci
         cnst4 = cnst3 + 4. * rci

         if ((cnst4 + cnst2).lt.1.e-10) write(*,*)' tmatrx: error!'
         trcoul(ld,n) = aovera*rcr/(cnst4 + cnst2)
         ticoul(ld,n) = aovera*(rci + cnst1 + cnst2)/(cnst4 + cnst2)

         tcplc = cmplx(trcoul(ld,n),ticoul(ld,n))
         tcpl = (tcpl - tcplc)/(2*zi*tcplc + 1.)
         trc = real(tcpl)
         tic = dimag(tcpl)
      endif
      
     


c *** do not unitarize for single scattering or if nifty(3) = 2
  284  write(kout,1030) tr(ld,n), ti(ld,n)
       write(kout,1031) trcoul(ld,n), ticoul(ld,n)
       write(kout,1032) trc,tic
c  284   continue

c *** calculation  of phase  shift+  eta

      eta = sqrt(1. + 4. *(tr(ld,n)**2 + ti(ld,n)**2 - ti(ld,n)))
      delta = 0.5 * atan2(2. * tr(ld,n),
     $                     (1. - 2. * ti(ld,n))) * 180./pi

c     if (eta .gt. 1.005) write(6,1020)
c     if (eta.gt.1.005) write(6,1010) eta, delta
c      deltap=atan2(dimag(2.*zi*tcpl+1.),real(2.*zi*tcpl+1.))
c      deltap = deltap/2.
c      write(*,*)'phase',ld,n,deltap
      
      etap(ld,n) = conjg(2.*zi*tcpl+1.)*(2.*zi*tcpl+1.)
c      write(*,*)ld,n,tcpl
      etap(ld,n)=sqrt(etap(ld,n))
c      write(kout,*) 'etap(ld,n)',etap(ld,n)
      smatrx = 2.*zi*tcpl+1.
      
      if (nspin.eq.1) then
         jspin=float(ld)+0.5
      else
         jspin=abs(float(ld)-0.5)
      endif
      write(11,777)real(smatrx),dimag(smatrx),ld,jspin

 283  continue
!      write(kout,*)ld,n,trc,tic
!      if (coulomb .eq. 0) go to 300

c *** include coulomb phase factor in nuclear amplitude
      if (coulomb > 0) then
         sigl = scoul(ld)
         if (coulomb .eq. 2) sigl = 0
         rhl5 = cos(2 * sigl)
         rhl6 = sin(2 * sigl)
         tr(ld,n) = rhl5 * trc - rhl6 * tic
         ti(ld,n) = rhl6 * trc + rhl5 * tic
         write(kout,1050) tr(ld,n), ti(ld,n)
      endif
      return

c *** formats

 777  format (2f12.6,1i5,2x,1f4.1)
 1010 format (1h+,85x,4heta=,f10.4,8h  delta=,f10.4)
 1020 format (41h0********** eta greater than 1 **********,//)
 1030 format (1h ,15h  t(k,k,k)=    ,35x,4e14.7)
 1031 format (1h ,15h  tcoul(k,k,k)=    ,35x,4e14.7)
 1032 format (1h ,15h  t_sub(k,k,k)=    ,35x,4e14.7)
 1050 format (20x,40hcoulomb phase modified nuclear amplitude/16h   t(k,
     $k.k) =   ,4e14.7)
 1055  format(2i4,2e14.6)
 3335 format(1h ,19h tr ti (nspin=3-6)=,8e13.5)

      end




c *** c.f. Raquel''s thesis, pag. 40 

      subroutine denom( k0, mnuc, enuck0, epk0, psfac, kg, wt, den )
        use params
        use inputs
        implicit real*8 ( a-h , m , o-z )

        integer code
        real*8 :: k , k2  , ko2 , k0
!        real*8 :: den(*),wt(*) !kg(ngp)
 
! AMORO 6/2/04     
!       dimension   den(1), kg(1) , wt(1)
        real*8::den(:),kg(:)
        real*8::wt(*)

!      common /params/   hbarc , pi , mp , mn , nz , na , nes , nwaves
!      common /inputs/   tlab , b , y1 , y2 , code , lmax , nang ,
!     $                  ngp , nr

************************************************************************

      mnuc2 = mnuc * mnuc
      mp2 = mp * mp
      mr = mnuc * mp / ( mnuc + mp )
      ko2 = k0 * k0

      n1 = ngp + 1
      piinv = 1. / pi
      rhl1 = mr + mr
      rhl2 = ( rhl1 + rhl1 ) * piinv
      rhl3 = piinv + piinv
      rhl4 = rhl3 + rhl3

      sum = 0.

      if ( nr .le. 0 ) then
         psfac = k0 * rhl1
      else
         psfac = 2. * k0 * epk0 * enuck0/( epk0 + enuck0 )
      endif

      write(99,*)'denom: ngp=',ngp
      do 10 i = 1 , ngp
         k = kg(i)
         k2 = k * k
         if ( nr .le. 0 ) then
            den(i) = k2 * wt(i) * rhl2/( k2 - ko2 )
         else
            enuck = sqrt( k2 + mnuc2 )
            epk = sqrt( k2 + mp2 )
            den(i) = rhl3 * k2 * wt(i)/( enuck + epk - enuck0 - epk0 )
            if ( nr .gt. 4 ) then
               den(i) = rhl3 * k2 * wt(i)/( (k2 - ko2)/(psfac/k0) )
            endif
         endif
         sum = sum + wt(i)/( k2 - ko2 )
  10  continue

      if ( nr .le. 0 ) then
         den(n1) = -sum * rhl2 * ko2
      else
         den(n1) = -sum * rhl4 * ko2 * enuck0 * epk0/( enuck0 + epk0 )
      endif

      return
      end


