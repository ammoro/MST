      subroutine ampnn(tcm,ifmst,sqmax,bqmax,fkq)
      use toff
      use parms
      use trace
      use ampnl
      use ampold
!      use constants
      
!      implicit real*8 (a-h,o-z),integer*4(i-n)
      implicit none
      integer lmax,fkq
      integer,save::ncall=1
      complex*16 :: det
      logical :: first,uu,ifmst
      character*10 nnpot
      real*8::sqmax,bqmax,tcm
      integer itype
      
!      namelist /amp/ ifkq,xkmax,xqmax,dk,dq,theta,nth,itype,icase
      namelist /amp/ ifkq,xkmax,xqmax,dk,dq,theta,nth,
     X         itype,icase ! deprecated; just for backward compatibility (AMM Jul'24)
      
c---------------------  Constants -------------------
      pi=acos(-1d0)
      ch=197.33
      pm=0d0 !projectile mass; default=938.93
      tm=0d0 !target mass; default=938.93
!      redm=pm*tm/(tm+pm)
c AMoro 30/09/04
!!      rm=938.9/2.
!      rm=redm
!      ttof=ch*ch*pi/rm
      

      write(99,*)ncall,'call to NNAMP'
      write(99,*)'ampnn: tcm=',tcm, ' ifmst=',ifmst
      ncall=ncall+1

c--------------------   I/O files  -------------------
c     input parameters
      inquire(file='nnamp.in',exist=uu)
      if(uu) then
         write(*,*) ' ***  USING AS INPUT FILE nnamp.in,',
     X          ' NOT stdin  ***'
         open(10,file='nnamp.in',status='old')
      endif      

c *** T-matrix
      open(14,file='tout.dat',status='unknown')

c     write np off-shell amplitudes on file 28
!      open(unit=113,file='ampnp.data',form='unformatted',
!     1               status='unknown')



c----------------------------------------------------------------

c Build factorials
!! AMORO (27/11/2003). Not necessary. it is done in main.f90 or in mstmain.f90
!!      call fctrls
 

c Calculate k-matrix & T-matrix
       write(99,*) 'Entering k with ifmst=',ifmst
       write(99,*)'nnamin: bqmax=',bqmax
       call k(tcm,ifmst)
       write(99,*)'bqmax=',bqmax

! Set some default values, in case they are not entered in
! the input file
       nth=46
       dk=0.1
       dq=0.1
       ifkq=1
       
       
       read(10,nml=amp)
       if (bqmax>1e-3) then
        write(99,*)'NN:Ignoring xkmax and using bqmax=',bqmax
        xkmax=bqmax
       endif
       if (sqmax>1e-10) then
        write(99,*)'NN:Ignoring xqmax and using=',sqmax
          xqmax=sqmax
       endif
       if (fkq>0) ifkq=fkq

       write(99,*)'NN: ifkq=',ifkq

!  Calculate scattering amplitudes: A,B...E
       close(20)
       if ((aeoff.ne.0).or.(aeon.ne.0).or.ifmst.or.(xsec.ne.0)) 
     & call ampall() ! (lmu)

!  Calculate M(a,b;alpha,beta) (tensor representation)
cc on-shell
      if ((mon.ne.0).or.ifmst) then
         call mabon() ! (lmu)
      endif
      write(99,*)  '<-Leaving mabon...'

cc off-shell
      if (moff.ne.0) then
         write(99,*) '-> Entering maboff...'
         call maboff() ! (lmu)
         write(99,*) 'Leaving maboff...'
      end if
      if (xsec.ne.0)  then
         write(99,*) '-> Entering xsec...'
         call nnxsec
          write(99,*) '<- Leaving xsec...'
       endif
      deallocate(qq,wq)
      deallocate(d01,dz1,dm1,dmm1,dpp1,dp1)
!      close(10); close(14)
      return
400   write(*,*)'Error opening NN input file.Aborting';stop
      end

************************************************************
*  END OF MAIN PROGRAM
************************************************************






