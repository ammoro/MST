      subroutine readsmat(ns,tmat,kk,z,mu,jj)
c*****************************************************************************
c      Reads S-matrix for nucleon-cluster scattering. 
c      to be used within the Multiple Scattering of the T-matrix formalism, 
c      described in:
c 
c      R.Crespo, A.M.Moro, I.J.Thompson, Nucl.Phys. A771, 26 (2006)
c******************************************************************************
       use scattering
       use constants
       implicit real*8(a-h,o-z)
       character*80 line
        integer npw,lcmax,l,istat,nspin,i,ns
        complex*16::tl, fc,tl1,tl2
        real*8:: sre,sim,sin2,jspin,kk,mu,z,jj,dsdt,t
        real*8, allocatable :: pl(:),delc(:)!,sc(:)
        complex*16,allocatable:: score(:,:)      
        complex*16:: tmat(nangles)
        
       npw=0
       lcmax=0
       zz=cmplx(0d0,1d0)
       write(99,*)'Readsmat: jj=',jj
       write(99,*)'**WARNING**: This version of MST does not'
       write(99,*)'consider spin-orbit part in N-core T-matrix'
c **** Preread file to determine number of partial waves
       rewind(ns)
310    read(ns,*,end=320,err=320) sre,sim,l,jspin
!       write(*,*)'sre,sim,l,jspin=',sre,sim,l,jspin
       npw=npw+1
       if (l>lcmax) lcmax=l
       goto 310

320    if(npw.eq.0) then
          write(*,*)'**ERROR**: Error reading S-matrix from unit',ns
          stop
       endif
       write(99,fmt='(" - readsmat: read",1i3," partial waves")')lcmax
     
       allocate(score(0:lcmax,2))
       score=cmplx(0.,0.)
c **** Read trully now       
       rewind(ns)
       do i=1,npw
          read(ns,*)sre,sim,l,jspin 
          if (abs(l+jj-jspin)<1.e-3) then !j=l+jj
             score(l,1)=cmplx(sre,sim)
             if (abs(jj)<1e-3) score(l,2)=score(l,1)
             if (l.eq.0) score(l,2)=score(l,1)
          else if (abs(abs(l-jj)-jspin)<1.e-5) then !j=l-jj
             score(l,2)=cmplx(sre,sim)
          else
             write(*,*)'**ERROR** reading unit',ns,'. This version' 
             write(*,*)' only allows 0x1/2 spins and 0x0 spins'
             write(*,*)'l,j,jspin=',l,jj,jspin
             write(*,*)' Aborting...'
             stop
          end if         
       enddo
       

       write(99,*)' External S-matrix:S(l+1/2) S(l-1/2)'
       do l=0,lcmax
       write(99,330)score(l,1),score(l,2),l
330    format(4f8.4,2x,1i3,2x,f6.3)
       enddo
       write(99,*)'-------------------------------------'
      
       allocate(pl(0:lcmax),stat=istat)
       if (istat>0) then
          write(*,*) 'readsmat: allocating memory for PL failed!';
          stop
       endif
340    format(2x,"# zp=",1f6.3," z=",1f6.2," mu=",
     & f8.2," kk=",f6.2," spin=",f4.1)
 
     
       if (coulmst.eq.0) then
          write(99,*)'*** CLUSTER T-MATRIX WITHOUT COULOMB ***'
          xgam=0.d0
       else
          xgam = zp*z*mu/137.036/kk/hbarc
       endif
       write(99,*)'readsmat: Sommerfeld param=',xgam
       allocate(delc(0:lcmax+1),stat=istat)
!       allocate(sc(0:lcmax),stat=istat)
   
       if (istat>0) then
          write(*,*) 'readsmat: allocating memory for DELC failed!';
          stop
       endif
       write(99,*)'Entering coul() with lmax=',lcmax+1
       call coul (xgam,delc,lcmax+1) !coulomb phase-shifts
       write(99,*)'Exiting coul'
       write(99,*)'delc(0:1)',delc(0:1)
 
       write(99,*)'readsmat:hbarc,kk,delc(0)=',hbarc,kk,delc(0)
       write(64,340) zp,z,mu/amu,kk,jj
       write(64,350)"Theta","dsigma/dw","t","dsigma/dt"
       write(64,350)"(deg)","(mb/sr)",
     &       "(GeV/c)^2)","(mb/(GeV/c)^2)"
350    format(3x,"#",a5,a10,a5,a10)
       do ith=1,nangles
          th = thmin + dth*(ith-1)
          thrad = th*pi/180.
          cthna = cos(thrad)
          call legpol( cthna, pl, lcmax )

c      Point Coulomb
          if (thrad<1e-6) thrad=1.e-6
          sin2=sin(thrad/2.)**2
          au=-xgam*log(sin2)+ pi + 2.*delc(0)
          fc=xgam*exp(zz*au)/(2.*kk*sin2)
                 
          tmat(ith)=fc !sum Coulomb amplitude 
!          tmat(ith)=0.  !do not sum Coulomb amplitude

         
          do l=0,lcmax
             tl1=(score(l,1)-1d0)/2./zz
             tl2=(score(l,2)-1d0)/2./zz
             tl=((dble(l+1)*tl1 + dble(l)*tl2))/kk
             tl=tl*exp(2*zz*delc(l)) ! coulomb phase
             tmat(ith)=tmat(ith) + tl*pl(l)
          enddo !partial waves
          dsdw=10*abs(tmat(ith))**2
          t=4*kk**2*sin2 ! momentum transfer**2
          t=t*hbarc**2/1e6 !convert to (GeV/c)^2
          dsdt=10*abs(tmat(ith))**2*(pi/kk/kk)*1e6/hbarc**2
          write(64,777) th,dsdw,t,dsdt
       enddo !angles
777    format(f8.3,1e12.3,1f10.4,1e12.3)
       write(64,*)'&' ! new set
       deallocate(pl)
       deallocate(delc)
       deallocate(score)
      end subroutine readsmat
