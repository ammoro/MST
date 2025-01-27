      subroutine lptps(elab,ifmst)
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

      use nlspfl
      use foptp
      use kinemt
      use params
      use switch
      use sizes
      use ranges
      use inputs
      use rcomn
      use spins
      use tcoul
      use tcomn
      use tcoulb
      use wr
      use bgrid
      use radgr
      use bprop
      use clebma
      use bscoul
      use qnumb
      use bkcall
      use blocks
      use radwf
      use etapbl
      use btmatc
      use bnoyes
      
      parameter(lfact=100)
      parameter(npts=200)
      implicit real*8 (a-h, o-z)

      character*80 sccsid
      integer :: klog=99
      real*8::  ko , ko2 
      real*8,allocatable::  theta(:)
      logical::ifmst
      real*8:: gp(npts),ws(npts)

c ***  -------------------------------------------------------------
c ***  Interaces: They are required to pass dynamic arrays 
c ***  as dummy arguments in subroutines
       interface 
          subroutine legpol(x, plofx, n)
            integer:: n
            real*8 :: x
            real*8 ,intent(out):: plofx(:)
          end subroutine legpol
       end interface

       interface
          subroutine plprme(x, plp, n)
          integer:: n
          real*8 :: x
          real*8, dimension(:),intent(out) ::plp
        end subroutine plprme
       end interface


! AMORO 5/4/2004 Commented (otherwise does not work in ifc)
!!$       interface 
!!$          subroutine gauss2(n, kode, a, b, gp, wt)
!!$            integer::n,kode
!!$            real*8:: a,b
!!$!            real*8, dimension(:),intent(out)::gp,wt
!!$            real*8::gp(200),wt(200)
!!$!             real*8, dimension(:),intent(out)::xx(200),ww(200)
!!$          end subroutine gauss2
!!$       end interface

       interface
          subroutine gauss3(rmaxr,quin,mquadi,mquado,xri,wri)
            integer:: mquadi,mquado
            real*8::rmaxr,quin
            real*8::xri(201),wri(201)
          end subroutine gauss3
       end interface

        interface
           subroutine vopt0a( vll, k, kp, neles, k0,i1,i2 )
           integer :: neles,i1,i2
           real*8 :: k,kp,k0
           complex*16 :: vll(:,:)
         end subroutine vopt0a
      end interface

      interface 
         subroutine vopt2l(v2l,k,kp,lmas2,k1)
           integer :: lmas2
           real*8:: k,kp,k1
           complex*16  ::v2l(:,:)
         end subroutine vopt2l
      end interface

      interface 
         subroutine vopt2nl(v2l,k,kp,i2,i1,n1,lmas2)
           integer :: i1,i2,n1,lmas2
           real*8:: k,kp
           complex*16  ::v2l(:,:)
         end subroutine vopt2nl
      end interface

      interface
         subroutine vopt2nlc(v2l,k,kp,i2,i1,n1,lmas2)
           integer :: i1,i2,n1,lmas2
           real*8 :: k,kp
           complex*16  ::v2l(:,:)
         end subroutine vopt2nlc
      end interface

      interface
         subroutine xsects(alph,xgam,sig0,ldmax,nspin,nsex)
         integer::ldmax,noangs,nspin,nsex
         real*8 :: ff(4),alph,xgam,sig0
        end subroutine xsects
      end interface

      interface
         subroutine cmatin(a, b, n)
           integer :: n
           real*8 :: a(:,:),b(:,:)
         end subroutine cmatin
      end interface
           
      interface
       subroutine vcoul(k,kp,vl,neles,type)
         integer::neles,type
         real*8:: k,kp,vl(:)
       end subroutine vcoul
      end interface

      interface
          subroutine formfactors(n1,n2)
            real*8:: k,kp
           integer :: n1,n2
         end subroutine formfactors
      end interface

      interface
         subroutine ffactbs(ncl,n1,n2)
            integer :: ncl,n1,n2
          end subroutine ffactbs
       end interface

      interface
         subroutine bound(ia,ib,ic,is2,j2a,lmoma,
     1   nodd,nramax,jhw,cmass,
     2   vmass,wzz,wal,wr0,wls,ffr,vdepth,bengy,dmat,kcheck,wrz,
     3   chis,drx,pnloc,wr0ls,wals)
         implicit real*8 (a-h,o-z)
         integer :: ia,ib,ic,is2,j2a,lmoma,nodd,nramax,jhw,kcheck
         real*8  :: ffr(:,:),chis(:)
         real*8  ::cmass,vmass,wzz,wal,wr0,wls,vdepth,bengy,dmat,
     1   wrz,drx,pnloc,wr0ls,wals
         end subroutine bound
      end interface

      interface
         subroutine sim2(fa,res,m,n,dr,nramax)  
           real*8  ::res,dr,fa(:)
           integer ::m,n,nramax
         end subroutine sim2
      end interface

      interface
         subroutine rtoq(densr,nramax,drx,ncl,n1,n2)
           integer:: ncl,n1,n2,nramax
           real*8 :: densr(:)
         end subroutine rtoq
      end interface

      interface
         subroutine sigcl(eta,scou,lmax)
           integer::lmax
           real*8 ::eta
           real*8:: scou(:)
         end subroutine sigcl
      end interface

      interface
          subroutine denom(k0,mnuc,enuck0,epk0,psfac,kg,wt,den)
            real*8 :: k0,mnuc,enuck0,epk0,psfac
            real*8 :: den(:),kg(:),wt(200)
          end subroutine denom
       end interface

       interface
          subroutine fmatrx(den,n1,nmax,nspin,nptmax,ldum )
            integer :: n1,nmax,nspin,nptmax,ldum
            real*8 :: den(:)
          end subroutine fmatrx
       end interface

c **** end interfaces ----------------------------------------------

c *** dimension of u depends on ngp (only need change in main)
c *** double size for spin 0 case,so can up n grid pts

      
       real*8,allocatable:: den(:) !,gp(:)
       real*8,allocatable:: v2l(:,:,:), v2l2(:,:,:)
     
      hbarc=197.3289e+0
      pi=acos(-1.d0)
!      pi=3.141593e+0
      kin=13
      kout=16

      open(unit=kin,file='mso.in',status='unknown')
      open(unit=kout,file='mso.out',status='unknown')    
      open(unit=11,file='smat.out',status='unknown')
      open(unit=25,file='xsec.out')
     

c DEPRECATED UNITS
c      open(unit=9,file='fort9.data',form='unformatted')
c      open (unit=77, file='tna.dat',status='unknown')
c      open (unit=78, file='xssha.dat',status='unknown')
c      open(unit=10,file='tape10',status='unknown')
c      open(unit=7,file='fort7.data',status='unknown')
c      open(unit=30,file='redish.data',form='unformatted')
c      open(unit=32,file='ampnp.data',form='unformatted')
c      open(unit=35,file='aaee.data',form='unformatted',status='unknown')
c      open(unit=40,file='twrcoul.data',form='unformatted')      
c      open(unit=45,file='tcoul.data',form='unformatted')
c      open(unit=50,file='wsf.data',status='unknown')
c      open (unit=52 , file='densli.dat',status='unknown')


      icex = 1


      call getins(elab)
      write(*,*)'elab=',elab
!      call denstgr()

      allocate(kk(ngp+1))
!      allocate(wt(ngp+1))
      allocate(den(ngp+1))
!      allocate(gp(ngp+1))
      gp=0d0;den=0d0;ws=0d0;kk=0d0 !initialize
!      write(*,2030)ngp+1
!2030  format('lptps: allocating',i3,' elements for kk,wt,den,gp')
      
      nang = iabs(nang)

      amu=931.49432
!      amass = na
      amass=masst
      zcm = acmp
      zch = achp
      zws = wsp
      nsex = 0
      nifty1 = nifty(1)

 60   continue

c *** cex(proton) mod  (m.s.)

!       if (nifty(1) .eq. 5) mp = 938.2796
!       if (nifty(1) .eq. 9) mp = 939.5731

c *** mnuc here refers to nucleus mass, mn is nucleon

!       mn = (nz * 938.2796 + (na - nz) * 939.5731)/na
       mn = (zt * 938.2796 + (at - zt) * 939.5731)/at
!      mnuc = mn * amass
       mnuc=masst*amu
       mp=massp*amu
!      write(kout,780) amass, mp, mn
       amass=masst
      write(kout,780) amass, mp, mn
      mnuc2 = mnuc * mnuc
      mp2 = mp * mp

 80   continue

c *** now k0 can have nonrelativistic and relativistic value

      el = tlab + mp
      ecm = sqrt( mnuc2 + mp2 + 2. * mnuc * el )
      if (nr.le.0)then
!         k0 = (amass/(amass+1.))**2 * 2.*mp*tlab
          k0 = (massp*amass/(amass+massp)) * 2.*amu*tlab
         k0 = sqrt(k0)
      else
         k0 = sqrt( (el * el - mp2) * mnuc2/ecm/ecm )
      end if
      plab = sqrt( el * el - mp2 )
      write(kout,960) tlab, el, ecm, k0, plab
      write(kout,980)
      write(kout,940)
      write(kout,950) tlab, k0
      n1 = ngp + 1
      n2 = n1 + n1
      a = k0

      write(klog,*) '+ Allocating fr,fi with',n1,n1
      allocate(fr(n1,n1)) 
      allocate(fi(n1,n1))
      allocate(fcr(n1,n1)) 
      allocate(fci(n1,n1))

c *** special check if b=0, use different distrb of points

      if (b .eq. 0.) a = 200.

      
      call gauss2 (ngp, kode, a, b, gp, ws)
      do 90 i1 = 1 , ngp
         kk(i1) = gp(i1)
         write(kout,*)'i1 , gp  = ',i1,kk(i1),ws(i1)
  90  continue
c      write(7,*)(ws(i1),i1=1,n1)

      kk(n1) = k0
      write(kout,*)'i1 , gp  = ',n1,kk(n1),ws(n1)

c *********
c       call wavfn(1,1,n1,xgama,10)
c        stop
c***************
c ** read numerical NN scattering amplitudes
!      call readnn(ifmst)
       write(99,*)'MSO: calculating factorials for l<=',lfact
       call logfac(lfact)
!       write(*,*)'MSO: calling calcnn'
       call calcnn(ifmst)
!      call mbess(kk,n1)

c **  2nd order calculations
      if(option3.ge.3)then
      call mbess(kk,n1)
      call matwf
      call newgrid(k0,n1)
      call denomg(kk,ws,n1)
c AMORO (27/11/2003)

!!      call factor
c ************** check
c      call rint2c(k0,k0,n1,n1,n1)
c      stop
c      write(kout,*)'call voptnl'
c     call vopt2nl(v2l,kk(n1),kk(n1),n1,n1,n1,2)
c     write(kout,*)'nonlocal1  ',v2l(0,1,1),v2l(0,2,1)
c     write(kout,*)'nonlocal1  ',v2l(1,1,1),v2l(1,2,1)
c     call vnl2(v2l2,kk(n1),kk(n1),n1,n1,n1,2)
c     write(kout,*)'nonlocal2  ',v2l2(0,1,1),v2l2(0,2,1)
c     write(kout,*)'nonlocal2  ',v2l2(1,1,1),v2l2(1,2,1)
c     sum1=v2l(0,1,1)+v2l2(0,1,1)
c     sum2=v2l(0,2,1)+v2l2(0,2,1)
c     sum3=v2l(1,1,1)+v2l2(1,1,1)
c     sum4=v2l(1,2,1)+v2l2(1,2,1)
c     write(kout,*)'sum=  ',sum1,sum2
c     write(kout,*)'sum= ',sum3,sum4
c     write(kout,*)'call voptnlc'
c     kcall = 1
c      call vopt2nlc(v2l,kk(n1),kk(n1),n1,n1,n1,25)
c     stop
c     write(kout,*)'non local  ',v2l(0,1,1),v2l(0,2,1)
c     write(kout,*)'nonlocal   ',v2l(1,1,1),v2l(1,2,1)
c     write(kout,*)'nonlocal   ',v2l(10,1,1),v2l(10,2,1)
c
c     write(kout,*)'call voptnlc'
c     kcall = 2
c      call vopt2nlc(v2l,kk(n1),kk(n1),n1,n1,n1,25)
c     write(kout,*)'non local  ',v2l(0,1,1),v2l(0,2,1)
c     write(kout,*)'nonlocal   ',v2l(1,1,1),v2l(1,2,1)
c     write(kout,*)'nonlocal   ',v2l(10,1,1),v2l(10,2,1)
c     call vopt2l(v2l,kk(n1),kk(n1),10,kk(n1))
c     write(kout,*)'local',v2l(0,1,1),v2l(0,2,1)
c     do 771 i=0,lfmx
c     write(kout,*)'tmatc= ',tmatc(n1,n1,i)
c771  continue
c     stop
c ************** end of check
      end if
c **
      ko = k0
      ko2 = ko * ko
      emnos = sqrt(ko2 + mnuc2)
      epios = sqrt(ko2 + mp2)
      if (nr.eq.0) then
         emnos = mnuc
         epios = mp
      endif
      s = (emnos + epios)**2

c *** now set up the denominator for either relativistic or
c *** non-relativistic case
      call denom( k0, mnuc, emnos, epios, psfac, kk, ws, den)

c ****** check
c       write(kout,*)'check to propagator'
c      do 105 i=1,n1
c     write(kout,*)deng(i,1),den(i)
c105  continue
c     stop
c *******************************

c ***
c *** set lmax if needed
c ***

      lclass = ko * 1.33 * (na**(.333333))/hbarc
      ldmax = lxmax
      if (lxmax .eq. 0) ldmax = 2.5 * lclass + 2.5
      xgam = 0.
      icex = 2
      write(kout,730) ldmax, lxmax
      write(kout,4020)


c *** stop if doing dirac.  this is temporary.

      if (nifty(11) .eq. 1) stop

c ***

c *** include coulomb, set uop ***

      if (coulomb .ne. 0)  then
!         zp = 1
!         if ( nifty(1) .eq. 9 ) then
!            zp = 0.
!         endif

c  ***   now xgam can have relativistic and nonrelativistic value
         if(nr.eq.0)then
             if (k0.lt.1e-6) then
                write(*,*)'** ERROR ** In lptps k0=',k0
                stop
             endif
!            xgam = zp*zt*mn*na/(na+1.)/137.036/k0
             if (zp*zt.lt.1e-6) then
                xgam=0
                write(*,*)'** WARNING ** In lptps xgam=0'
             else
                xgam = zp*zt*mn*at/(at+1.)/137.036/k0
             endif
         else
            write(*,*)'*** In lptps using RELATIVISTIC kinematics'
            xgam = zp * zt/137.036/(ko * (emnos + epios)/emnos/epios)
         end if
         write(kout,*)' MSO: Sommerfeld parameter=',xgam
         write(kout,*)' MSO: Reduced mass',mn*at/(at+1.)
         write(kout,*)' MSO: kcore=',k0
         alph = (zt - 2)/3.
         alph = alph/2./(2. + 3. * alph)
!         if (na .lt. 4) alph = 0.
          if (at .lt. 4) alph = 0.

c        *** set up sigl ***

         rhoc = ko * rcut/hbarc
         lcmax = ldmax+1!Amoro
         if (lcmax .eq. 1) lcmax = 2
         eta = xgam

         
! Commented 21/05/2004
         if (lcmax<0) then
            write(*,*)'*** ERROR *** In subroutine COUL lmax=',lcmax
            stop
         endif
         allocate(scoul(0:lcmax))
!�         call sigcl (eta, scoul, lcmax)
! AMoro (24/09/06): lmax => lcmax, cph => scoul
!         call coul(eta,cph,lmax)
          call coul(eta,scoul,lcmax)
!
         sig0 = scoul(0)
         if (nsex .eq. 2)  go to 430
      endif

c ***
c *** do loop over l for this e and  state
c ***

      nspinm = 2
      if (nifty(6) .eq. 8) nspinm = 6

      write(klog,*)'+ Allocating for fnoyes:',ldmax,nspinm
      allocate(fnoyes(0:ldmax,nspinm))
      allocate(tbr(0:ldmax,nspinm))
      allocate(tbi(0:ldmax,nspinm))
      allocate(tr(0:ldmax,nspinm))
      allocate(ti(0:ldmax,nspinm))
      allocate(trcoul(0:ldmax,nspinm))
      allocate(ticoul(0:ldmax,nspinm))
      allocate(etap(0:ldmax,nspinm))
      trcoul=0d0;ticoul=0d0;
!      tr=0d0;ti=0d0;
!      tbr=0d0;tbi=0d0;
!      etap=0d0
!      fnoyes=0d0
  
     

! This has to dissappear at some point (AMoro)
      nptmax = ldmax
      nspdim =  nspinm
      call optp (n1, ldmax)

      do 420 ldum = 0,ldmax
         ld = ldum

c ***
c ----------   do loop over 2 - 6 spin sates(if necessary)--------
c ***
         do 410 nspin = 1,nspinm
            n = nspin
            if (n .eq. 3 .or. n .eq. 5) go to 410
            nspine = nspin

c *** set spin indices needed for 1/2 x 1/2

!!$            if (nifty(6) .eq. 8) then
!!$               nspina = nspin - 1
!!$               nspinb = nspin
!!$               nspinc = nspina + 2 * (5 - nspin)
!!$               nspind = nspinc + 1
!!$            endif

c *** special storage r0(11)-nspin=5 unreachable with our l sum scheme
c *** store in ro(00)-nspin=2 for ldum=1(l=0)

            nmax = n1
            if (nspin .gt. 2) nmax = n2

c *** set up the f-matrix ***
            call fmatrx( den , n1 , nmax , nspin , nptmax, ldum )

 260        continue

c *** calculate the on-shell r matrix from inverse of f * u ***
            call rmatrx ( psfac , n1 , n2 , nmax , nspin )

c *** calculation  of    t-matrix  (on - shell),
c *** t normalized to exp(i del)*sin(del)

 2424       continue
            call tmatrx(sigl, aovera, ldum, nspin, njump )
            call noyes(kk,psfac,n1,ld,nspin)
c            if(ldum.le.2.and.n.eq.1)then
c           call ftvoff(ldum,n,n1)
c            end if
            if(igwf.eq.1.and.ld.lt.lstore)then
               call wavfn(ld,nspin,n1,xgam,lxmax)
            endif
 410     continue
 420    continue !end loop in ldum



  430  continue


! 430  if (nsex .ne. 3) go to 435
!      write(kout,740) nsex
!
!      do 434 ll = 0,ldmax
! 434     write(kout,1031) ll, (tr(ll,np2), ti(ll,np2), np2=1,6)

! 435  continue
!      write (kout,1110)

c-----calculate angular disrtbn    -----------------------------------

       write(kout,*)'rcut=',rcut
       call xsects(alph, xgam, sig0, ldmax, nspin, nsex)
c *** Release memory for (fr,fi,tr,ti,ul) 
        write(klog,*) '- Deallocating fr,fi,tr,ti'
        if(.not.allocated(fi)) then
           write(klog,*) 'fr not allocated!!!'
        end if
        deallocate(fr,fi)
        deallocate(tr,ti)
        deallocate(ul)


 620  continue


       

 650  write(kout,841)
      return

cc ****  prints eta in unit 10
      do 655 ld=0,ldmax
         write(kout,*)ld,real(etap(ld,1)),real(etap(ld,2))
      write(10,*)ld,real(etap(ld,1)),real(etap(ld,2))
 655  continue
      stop
c     _______________     formats
c ***
 227  format(i3, 8e15.6 / 200(3x, 8e15.6 / ))
 700  format (15h p, n, then cex)
 730  format (1h0,'lmax(calculated) = ',i4,5x,'lborn = ',i4,5x,
     $        'lmax(read) = ',i4)
 740  format (1h1,5i4)
 780  format (' ********** a, mp, m(nucleon) =', 3f10.4,
     $        ' **********', /, /, 7x,
     $        ' tlab', 9x, 'eplab ', 8x, 'ecm(sqrts)', 4x,
     $        'kcm', 11x, 'plab')
 820  format (47h r matrix calculated with ms series up to order,i3)
 841  format(12h eof in main ,i3)
 928  format(1h0,33h re/im f matrix before inversion )
 930  format (1h ,26h  normalized r(k,k..k)=   ,24x,2e14.7)
 940  format (1h0,10x,8h  energy,10x,15h    momentum   ,/)
 950  format (1h ,7x,f14.4,10x,f14.4)
 960  format (1h ,5f14.4)
 980  format (1h0,31h      relativistic  calculation)
 1010 format (1h+,85x,4heta=,f10.4,8h  delta=,f10.4)
 1020 format (41h0********** eta greater than 1 **********,//)
 1030 format (1h ,15h  t(k,k,k)=    ,35x,4e14.7)
 1031 format (1h ,16h t(k,k,k),for l=,i3,/,10(12e10.3,/))
 1040 format (38h now with exact coulomb after matching)
 1050 format (20x,40hcoulomb phase modified nuclear amplitude/16h   t(k,
     1k.k) =   ,4e14.7)
 1070 format (1h0,30h----- angular   momentum  =   ,i4,5h-----)
 1080 format (1h ,8x,23h  born  approximation  )
 1090 format (19h0  th-c.m.  cos(cm),3x,10hdsig/dw-cm,2x,11hds/dw(b.a.),
     13x,1ht,4x,6hth-lab,2x,8hcos(lab),2x,5hds/dw,4x,11hds/dw(b.a.),2x,
     210href(th)-cm,1x,3himf,4x,6hsignuc,4x,7hsigcoul//)
 1091 format(17h0th-c.m.  cos(cm),3x,4hq2fm,7x,10hdsig/dw-cm  ,
     1         1x,10hds/dw(nuc) ,1x,11hds/dw(b.a.) ,2x,5hpolar ,4x,
     2 7hnoflip ,4x,4hflip ,5x,4hpnoo ,4x,5hdlomo ,4x,
     35hclmoo ,2x,7hsigcoul  /1h ,18x,
     420hza ,zb ,zc ,zd ,ze  /1h ,18x,
     520hzab,zbb,zcb,zdb,zeb  //)
 1093 format(1h ,f7.0,2f8.3,3e11.3,f8.3,2e11.3,3f9.3,e11.3)
 1100 format(1h ,f7.0,f8.3,2e11.3,2f8.1,f10.3,2e11.3,2f9.3,e12.4
     1,e10.2)
 1110 format (1h+,21x,6hnoflip,7x,8hnoflip-b,5x,5hpolar,9x,4hflip,9x,
     18hflip-b,r,//)
 1120 format (1h ,18x,10e11.3)
 1921 format(1h ,7h sigot=,e13.5,7h sigit=,e13.5,
     27h sig2t=,e13.5,2he=,f8.2,5hxgam=,e13.5)
 1979  format(1h ,17h ur ui(k0,k0,k0)=   ,2e13.5)
 2211  format(1h ,' ..........   nspin=   ',i4,' ...... ')
 3335 format(1h ,19h tr ti (nspin=3-6)=,8e13.5)
 4000 format(' ', f7.2, 3e11.3)
 4010 format(1h ,f7.2,2e11.3)
 4020 format(1h1)
      end


c ***
c *** Read in input parameters
c *** 
      subroutine getins (elab)
************************************************************************
      use params
      use switch
      use sizes
      use ranges
      use inputs
      use blocks
      use halo
      use qnumb
      use bgauss3
      use bwfnum
      use radwf
      use wr
      use ffch
      use bsflip
      use bion2
      use bhalodn
      use brms
  

      implicit real*8  ( a-h , o-z )
    
c *** namelists for input/output
      namelist /mso/ kmt, offshell,tlab,igwf,lstore,nang,
     &              nr,ioptls,ioptcnt,lxmax,ngp,
     &               kode,b,ymin1,ymin2,
     &               ucnr, ucni, usnr, usni, br,bi
!      namelist /coul/ coulomb,rcut, achp,  achn
      namelist /coul/ coulomb,rcut, alfa,beta,gama
      namelist /proj/ massp,zp
      namelist /targ/ ncluster,masst, at,zt
      namelist /highorder/ second,isflip,lamx,lbmx,l1mx,lfmx,
     &                     rmaxr, quin, mquadi, mquado,closure

c-----------------------------------------------------------------
c  &mso/ namelist:
c     tlab=kinetic energy
c     kmt=T (default) kmt potential
c        =F watson potential, no a/(a-1)
c     
c     coulomb: 
c     rcut=
c
c     igwf=T calculates wavefunction up to lstore
c     lstore=max. partial wave for wf
c     nr=0
c     ucnr, ucni=central renormalization constants
c     usnr, usni=s.o.            "            " 
   
c  &coul / namelist:
c     coul=0  No coulomb
c     coul=1  Coulomb amplitude and phase change included
c     coul=2  Coulomb amplitude but no phase change included
c     coul=3  Realistic charge distribution (HO)
c     coul=4  Uniform sphere of radius 1.3*A**(1/3)
c     alfa, beta, gama= HO parameters


c  &proj/
c     massp=mass number of projectile (a.m.u.)
c     zp=charge


c  &targ/ namelist:
c     ncluster=number of clusters
    
c  &cluster/ namelist:
c     type=0 HO model for density
c         =1 WS
c         =2 Read external density
c     nz,na= number of protons/neutrons
c     achp= size parameter (charge) for protons
c     achn= size parameter (charge) for neutrons
c     acmp= size parameter (mass)   for protons
c     acmn= size parameter (mass)   for neutrons


c  &highorder/ namelist (by now 2nd order):
c     second: former option3
c     closure=true/false
c     isflip=
c     ich=
c     ion2=
c     lamx=
c     lbmx=
c     l1mx=
c     lfmx=
c     rmaxr=maximum radius
c     quin= maximum radius for inner region
c     mquadi=number of inner quad points (multiple of 6)
c     mquado= number of outer quad points (multiple of 6)


c------------------------------------------------------------------
       real*8 imap,imcp,iman,imcn
       complex*16 ap,cp,an,cn
       complex*16 aack,ccck,cint2d
       dimension formf(4)

c *** open statements are standard fortran 77 and may not be necessary
c *** depending upon the operating system of the machine in which the
c *** code is run.

c     open(3,file='tape3',status='unknown',form='unformatted')
c     open(5,file='tape5',status='unknown')
c     open(6,file='tape6',status='unknown')
c     open(9,status='scratch',form='unformatted')

c * call to ctchsig necessary for ridge to catch and ignore underflows *

c      call ctchsig()

c *** mnuc here refers to nucleus mass, mn (in common) to nucleon
c *** read  momenta  in  centre _of_mass
c *** nr <= 0  is  nonrelativistic  nr >= 0 relativistic  case
c ***    >  4  for approximate klein-gordon

c *** set some default values before reading
      tlab=0.
      coulomb=0
      offshell=.true.
      kmt=.true.
      igwf=0
      lstore=0
      nr=0
      ucnr=1.; ucni=1.
      usnr=1.; usni=1.
      bi=0.; br=0.
c *** -------------------------------------

      read(kin,nml=mso)
      if (offshell) then
         write(*,*) '- Using off-shell NN amplitudes'
      else
         write(*,*) '- Using on-shell NN amplitudes'
      endif
      if (tlab<1.e-5) then
         tlab=elab
      else
         elab=tlab
      endif
      write(99,*)'getins:elab=',tlab
      read(kin,nml=coul)
      read(kin,nml=highorder)
      read(kin,nml=proj)
      read(kin,nml=targ)


       irmx = mquadi + mquado
       write(kout,*)'ioptls=',ioptls
       write(kout,*)'itydn9=', itydn9,'itydn11=',itydn11
       write(kout,*)'itycmdn=',itycmdn
       write(kout,*)'qnumbrs=',lamx,lbmx,l1mx,lfmx
       write(kout,*)'radial grid  ',rmaxr,quin,mquadi,mquado
       write(kout,741) nr, lxmax
       write(kout,770) tlab

      
      if (nang .eq. 0) then
         nang = 3
      endif
      if ((ymin1 + ymin2) .eq. 0.) then
         ymin1 = -2.
         ymin2 = -4.
      endif
      write(kout,750) ngp, kode, b, nang, ymin1, ymin2
      if (ngp .le. 0) then
         write(kout,*) 'number of gauss points cannot be less than one'
         stop
      endif
      write(kout,790) ngp
      write(kout,880) achp, acmp, wsp, achn, acmn, wsn
      write(kout,830) nz, na, (nifty(n),n=1,20)
      write(kout,910) nes, nwaves

        if(na.ne.16)then
        write(kout,*)'read wave functions for s and p only'
        endif
        write(kout,*)'reading from 15'
      return

c *** format statements ***

 741  format (1h0,5i4)
 750  format (1h ,i4,i3,e10.2,i4,2f5.2)
 770  format (bz,f10.4)
 790  format (1h ,18hno of grid points=,i4)
 830  format (1h ,2i5,1x,i2,i2,8i1,10i3)
 880  format (bz,8f10.4)
 910  format (' nes, nwaves = ',2i5)

      end





c***********************************************************************
      subroutine noyes(kk,psfac,n1,ld,nspin)
c***********************************************************************
c *** calculates the kowalski-noyes f ratio for the n-a transition
c *** matrix
      use params
      use nlspfl
      use foptp
      use tcomn
      use bnoyes

      implicit real*8 (a-h,o-z)
      real*8:: kk(n1), kkr,kmt
      complex*16 z
      complex*16 taux, rhalf

      z = (0.,1.)
      kmt = na/(na-1.)
c      write(*,*) 'kk=',kk(n1)
c      write(*,*) hbarc
      lg = kk(n1)/hbarc*1.3*(na*1.)**0.333333
      lg = lg+1
c     kmt = 1
c      write(*,*) 'fnoyes',size(fnoyes)
c      write(*,*)' tr,ti',size(tr),size(ti)
      taux = tr(ld,nspin) + z*ti(ld,nspin)

      do 270 j=1,n1
      rhalf = 0.
         do 280 j3=1,n1
         rhalf = rhalf + (fr(j,j3)+z*fi(j,j3))*u(j3,n1,nspin)
!     *            (ur(j3,n1,nspin)+z*ui(j3,n1,nspin))
  280    continue
      rhalf = -psfac*rhalf
      fnoyes(j,nspin) = rhalf*(kmt + z*taux)
c     fnoyes(j,nspin) = fnoyes(j,nspin)/taux
  270 continue
c     if(ld.eq.1.or.ld.eq.lg)then
!      if (nspin.eq.2) then
c      write(20,*)ld
!      do 290 j=1,n1
!      kkr = kk(j)/hbarc
!      write(20,777) kkr, real(fnoyes(j,1)),real(fnoyes(j,2)),
!     *              dimag(fnoyes(j,1)), dimag(fnoyes(j,2))
c      write(*,777) kkr, real(fnoyes(j,1)),real(fnoyes(j,2)),
c     *              dimag(fnoyes(j,1)), dimag(fnoyes(j,2))
!  290 continue
!      endif
c     endif
  777 format (5e12.3)


      return
      end

      

      subroutine calcnn(ifmst)
      use params
      use switch
      use sizes
      use ranges
      use inputs
      use bgrid
      use blocks
      use qnumb
      use wr
      use ffch
      use bion2
      use amps
      use amphe
  
      implicit real*8  ( a-h , o-z )

      real*8 imap,imcp,iman,imcn
      complex*16 aack,ccck,cint2d
      logical::ifmst,ifmso
      real*8::sqmax
      real*8::bqmax
      integer::ifkq=1

      ifmst=.true.
      sqmax=0d0; bqmax=0d0

       tcm=tlab/2.
       if (tcm<1.e-5) then
          write(*,*) '**ERROR*** Energy not specified!. Aborting'
          stop
       endif
       write(kout,*)'-Calling ampnn from MSO with tcm=',tcm
       call ampnn(tcm,ifmst,sqmax,bqmax,ifkq)
!      check
!       write(*,'(4g12.6)')cpp(:),cpn(:)
!       stop
       return
      end subroutine calcnn


c *** Readnn: old subroutine to read NN amplitudes from redish.data and 
c *** ampnn.data 
      subroutine readnn (ifmst)
************************************************************************
      use params
      use switch
      use sizes
      use ranges
      use inputs
      use bgrid
      use blocks
      use qnumb
      use wr
      use ffch
      use bion2
!      use ampl
      use amps
      use amphe
  
      implicit real*8  ( a-h , o-z )

      real*8 imap,imcp,iman,imcn
      complex*16 aack,ccck,cint2d
      logical::ifmst
      real*8::sqmax=0.
      real*8::bqmax=0.
      integer::ifkq=1

      if(option3.ge.2)then
      ion1=option2
      if (option3.eq.4.and.ion2.eq.1)option2=1
      if (option3.eq.4.and.ion2.eq.2)option2=2
        n1 = ngp+1
        read(35)nkmx1,nqmx1
        allocate(q(nqmx1))
        allocate(qq(nkmx1))
        read(35)q,qq
        do 100 ir=1,5
        write(kout,*)'ir=',ir
        read(35)ap,cp
!        write(kout,*)'calling mtmatt'
        call mtmatt(kk,n1,lfmx,ir)
 100    continue
        option2=ion1
      endif
c
      read(30)nkmx1,nqmx1
      allocate(q(nqmx1))
      allocate(qq(nkmx1))
      write(kout,*)'in 30 nkmx1,nqmx1=',nkmx1,nqmx1

      allocate(app(nkmx1,nqmx1))
      allocate(cpp(nkmx1,nqmx1))
      read(30)q,qq,app,cpp
     
      
!      allocate(ap(nkmx1,nqmx1))
!      allocate(cp(nkmx1,nqmx1))
!      read(30)q,qq,ap,cp
111   format (1f14.6)
      
      write(kout,*)'- Reading nn amplitudes'
      read(32)nkmx1,nqmx1
      write(kout,*)'- In 32 nkmx1 x nqmx1=',nkmx1,nqmx1
! This will have to be changed to right bounds (ngpxngp ???)
     
c      allocate(an(nkmx1,nqmx1))
c      allocate(cn(nkmx1,nqmx1)) 
c      read(32)q,qq,an,cn

      allocate(apn(nkmx1,nqmx1))
      allocate(cpn(nkmx1,nqmx1)) 
      
      read(32)q,qq,apn,cpn
c      endif
c
c ***** this is a check to the interpolation of the scattering
c *****  amplitudes
c      nth = 46
c      dth=180./(nth-1.)
c      pch = 1.55283
c      do 100 ith=1,nth
c       th =(ith-1.)*dth
c       aack = 0.
c       ccck = 0.
c       cth = cos(th*pi/180.)
c       xkk = pch*sqrt((1+cth)/2.)
c       xq = pch*sqrt((1.-cth)*2.)
c       aack = cint2d(qq,q,ap,xkk,xq,nkmx1,nqmx1,5,200)
c       ccck = cint2d(qq,q,cp,xkk,xq,nkmx1,nqmx1,5,200)
c       write(kout,107) real(aack),dimag(aack),real(ccck),dimag(ccck)
c100    continue
c107     format (4e12.4)
c *************************************************
c
      return

c *** format statements ***

 741  format (1h0,5i4)
 750  format (1h ,i4,i3,e10.2,i4,2f5.2)
 770  format (bz,f10.4)
 790  format (1h ,18hno of grid points=,i4)
 830  format (1h ,2i5,1x,i2,i2,8i1,10i3)
 880  format (bz,8f10.4)
 910  format (' nes, nwaves = ',2i5)

      end


       subroutine redish(kp,k,k0,cthnuc,xmn,tapb,te)
c*** raquel crespo addition ... ************************
c
       use params
       use switch
       use blocks
!       use ampl
       use amps
       use halo
       use inputs
       use blocks
 
       implicit real*8 (a-h, i, k, m, o-z)

       real*8:: kp,k,k0
       complex*16 :: taux,cn
       complex*16 :: cint2d
       complex*16 :: tapb(2),te(2)

c ***  qp and qqp are transfer and total n-n momentum
!       if (option2.eq.1)then
       if (.not.offshell) then         
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

!       if (.not.associated(qq)) then
!          qq=>Bq
!          q=>Sq
!          write(*,*) "qq=",qq(1),"q=",q(1)
!       endif

       if (.not.allocated(app)) then
          write(*,*)'**WARNING** In redish, app not allocated!'
          stop
!       else
!          write(*,*) "app",app
       endif
       
       if (zp>0.) then
c          tapb(1) = cint2d(qq,q,app,qqp,qp,nkmx1,nqmx1,5,200) !App
c          tapb(2) = cint2d(qq,q,apn,qqp,qp,nkmx1,nqmx1,5,200) !Apn
           tapb(1) = cint2d(qq,q,app,qqp,qp,nkmx1,nqmx1,5,nqmx1) !App
          tapb(2) = cint2d(qq,q,apn,qqp,qp,nkmx1,nqmx1,5,nqmx1) !Apn

c          te(1) = cint2d(qq,q,cpp,qqp,qp,nkmx1,nqmx1,5,200) !Cpp
c          te(2) = cint2d(qq,q,cpn,qqp,qp,nkmx1,nqmx1,5,200) !Cpn
          te(1) = cint2d(qq,q,cpp,qqp,qp,nkmx1,nqmx1,5,nqmx1) !Cpp
          te(2) = cint2d(qq,q,cpn,qqp,qp,nkmx1,nqmx1,5,nqmx1) !Cpn
       else
!          tapb(1) = cint2d(qq,q,apn,qqp,qp,nkmx1,nqmx1,5,200) !App
!          tapb(2) = cint2d(qq,q,app,qqp,qp,nkmx1,nqmx1,5,200) !Ann=App          
          tapb(1) = cint2d(qq,q,apn,qqp,qp,nkmx1,nqmx1,5,nqmx1) !App
          tapb(2) = cint2d(qq,q,app,qqp,qp,nkmx1,nqmx1,5,nqmx1) !Ann=App
!          te(1) = cint2d(qq,q,cpn,qqp,qp,nkmx1,nqmx1,5,200) !Cpn
!          te(2) = cint2d(qq,q,cpp,qqp,qp,nkmx1,nqmx1,5,200) !Cnn        
          te(1) = cint2d(qq,q,cpn,qqp,qp,nkmx1,nqmx1,5,nqmx1) !Cpn
          te(2) = cint2d(qq,q,cpp,qqp,qp,nkmx1,nqmx1,5,nqmx1) !Cnn
       endif
!       tapb(2)=tapb(1)

       

c ** change from scattering amplitude to transtion amplitude
c ** and to landau conventions
       tapb=-tapb/2./pi/pi/mp/hbarc
       te= -te/2./pi/pi/mp/hbarc
       return
       end




