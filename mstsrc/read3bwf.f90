	subroutine readwf(onlyelastic)          
c*****************************************************************************
c      Reads 3-body wavefuncion in {r,R} coordinates or in the hyperradius.
c      to be used within the Multiple Scattering of the T-matrix formalism, 
c      described in:
c 
c      R. Crespo and R. C. Johnson, Phys. Rev. C 60, 034007 (1999)
c      R.Crespo, I.J.Thompson, A.A.Korsheninnikov, Phys.Rev. C66, 021002 (2002)
c      R.Crespo, A.M.Moro, I.J.Thompson, Nucl.Phys. A771, 26 (2006)
c      R.Crespo, I.J.Thompson, A.M.Moro, Phys.Rev. C 74, 044616 (2006)
c      A.M.Moro,M.Rodriguez-Gallardo,R.Crespo,I.J.Thompson,Phys.Rev.C 75,017603 (2007)
c******************************************************************************



C  coupling order |[(S1,S2)S, (LX,LY)LL, Jnn], JC; JTOT>
C
c 	where S1=S2=1/2 and JC = core spin (zero for 6Li and 6He).
c       (IWF.EQ.1) OR (IN=IWF-1.EQ.0)   elastic
c       (IWF.GT.1) OR (IN=IWF-1.GT.1)   inelastic  
c
	use wfns
        use scattering
        use trdens
        
      IMPLICIT REAL*8(A-H,O-Z)
      logical:: onlyelastic
      real*8:: jlast=-1
      real*8:: iparlast=-2,elast=0.0
      integer:: ijpi=-1
      integer:: jimax=0,jfmax=0

      real*8:: x2,y2,x2avg,y2avg,rhorms
      real*8:: cx,cy
      
      CX=SQRT((m2*m3)/(m2+m3))
      CY=SQRT((m4*(m2+m3))/(m4+m2+m3))

C
c --------------------------------------------------------------
c  Dummy read to set dimensions:
c -------------------------------------------------------------
        write(99,*)'Prereading wf to set dimensions'
        write(*,*)
        write(*,*)'------------- J Pi components --------------'
        rewind(nfl)
        IWF=0
        nma=0
        mllmx=0
11        READ(NFL,2,end=998) NRXY,RSTEP,II,en,TN 
2       format(i4,f10.5,i4,f10.5,f10.6)
c        write(*,*)'IWF=',IWF,'ens(in)=',ens0  
	if(NRXY.gt.MRXY) stop 'MRXY'
!	WN = 0.0
!	VN = 0.0
! AMoro (25/11/2003) 
        if (onlyelastic.and.(iwf.eq.1)) then 
           write(*,*)!'========================================='
           write(*,'(A)')" Onlyel=T =>Skipping non-elastic components"
           write(*,*)!'========================================='
           goto 998
        endif
        
	IWF=IWF+1

        IF(IWF.EQ.1)ITYTR=elastic
        IF(IWF.GT.1)ITYTR=inelastic
      DO 201 IVERT=0,0
      READ(NFL,*) 
      READ(NFL,*) 
     
      ima=0
      DO 111 IA=1,ma
      READ(NFL,915) LX,LY,LL,S,JNN,IC,JC,JTOT 
!      write(*,915) LX,LY,LL,S,JNN,IC,JC,JTOT    
      if (max(lx,ly)>mllmx) mllmx=2*max(lx,ly)
      
      
      if(LX.lt.0) go to 121
      ima=ima+1
      IPAR = (-1)**(LX+LY)
 !     if ((ia.eq.1).and.iwf<3) write(*,916)JTOT,PSIGN(IPAR+2),EN
      if ((jlast.ne.jtot).or.(ipar.ne.iparlast).or.
     & (elast<0.and.abs(elast-en)>0.001)) then
         ijpi=ijpi+1
         write(*,916)ijpi,JTOT,PSIGN(IPAR+2),EN
      endif
      
      if (ijpi>0) then
         jfmax=max(jnn,jfmax) 
      else
         jimax=max(jimax,jnn)
      endif

!       write(*,916)JTOT,PSIGN(IPAR+2),EN
916    format('  Component: J,PI,=',i3,f5.1,a1,' at',f8.3,' MeV')
      jlast=jtot
      iparlast=ipar
      elast=en
      if(IC>1.or.JC>0.) then
	write(6,*) 'This program not written for core excitation!'
	stop
	endif
C
!40      RV(1) = 0.0
        DO 51 IX=1,NRXY
!	RV(IX) = (IX-1)*RSTEP
	read(NFL,*) IIX
51	READ(NFL,*) (WF0,IY=1,NRXY)
	XN = 0.0
        DO 54 IX=1,NRXY
        DO 54 IY=1,NRXY
54      continue
 	read(NFL,*) PNORM
81    FORMAT(' Input ch.',I3,' ( lx,ly,L,S=',4i3,') norm =',f9.5)
111   CONTINUE
121   if (ima>nma) nma=ima
      read(NFL,*) TNORM
201   CONTINUE
      indx = IWF-1

c AMoro addition
      if ((nwfmax>0).and.(iwf.eq.nwfmax)) goto 998      
      
      go to 11
998    NWF=IWF
      

       write(*,*)'  Read',nwf,'states'
       write(*,*)'  Max number of components =',nma
       write(*,*)'  MLLMX =',mllmx
      write(*,*)'-----------------------------------------------------'

c *** Allocate according to determined bounds     
      dmax=jimax+jfmax
      mll=dmax+isnx
      write(99,*)'isnx,ndel,dmax,mll',isnx,ndel,dmax,mll

      allocate(wf(nrxy,nrxy,nma,2))
      wf=0d0

      allocate(tnq(8,nma,nwf))
!      allocate(tqn(8,nma,nwf))
      allocate(na(0:nwf))
      allocate(ens(0:nwf))
      allocate(rv(mrxy))
      allocate(lvalmax(0:nwf))
     
      write(99,*)'+ Allocating ffcore: mdelta,mll,nwf=',mdelta,mll,nwf
      allocate(ffcore(mdelta,0:mll,0:nwf))
      allocate(ffcth(nangles,0:mll,0:nwf))
      allocate(ffval(mdelta,0:mll,0:isnx,0:dmax,0:nwf))
      allocate(ffvth(nangles,0:mll,0:isnx,0:dmax,0:nwf))
      ffcore=0d0
      ffval=0d0
      ffcth=0d0
      ffvth=0d0

c ----------------------------------------------------------
c NOW READ WF FOR THE SECOND TIME      
c------------------------------------------------------------
       ijpi=-1
       jlast=-1
       iparlast=-2
       elast=0.0
       REWIND NFL    
       do IWF=1,nwf
        
!1      READ(NFL,2,end=999) NWF,NRXY,RSTEP,II,ens(IWF),TN
1      READ(NFL,2) NRXY,RSTEP,II,en,TN 
          ens(iwf-1)=en
                 write(*,*)'IWF=',IWF-1,'ens(in)=',ens(IWF-1)  
        IF(IWF.EQ.1)ITYTR=elastic
        IF(IWF.GT.1)ITYTR=inelastic
        WN=0.0
        XRM=0.0
      DO 200 IVERT=0,0
      READ(NFL,*) 
      READ(NFL,*) 
      DO 110 IA=1,ma !channels (lx,ly,ll,s,jnn)
      READ(NFL,915) LX,LY,LL,S,JNN,IC,JC,JTOT
	if(LX.lt.0) go to 120
915	format(3x,6i3,2f4.1)
	NA(IWF) = IA
	TNQ(1,IA,IWF) = LX
	TNQ(2,IA,IWF) = LY
	TNQ(3,IA,IWF) = LL
        TNQ(4,IA,IWF) = S
	TNQ(5,IA,IWF) = JNN
!	TNQ(6,IA,IWF) = IC
!	TNQ(7,IA,IWF) = nint(2.*JC)
!	TNQ(8,IA,IWF) = nint(2.*JTOT)

!        TQN(1:5,IA,IWF)=TNQ(1:5,IA,ITYTR)
        IPAR = (-1)**(LX+LY)
      if ((jlast.ne.jtot).or.(ipar.ne.iparlast).or.
     & (elast<0.and.abs(elast-en)>0.001)) then
         ijpi=ijpi+1
         write(*,916)ijpi,JTOT,PSIGN(IPAR+2),EN
      endif
      jlast=jtot
      iparlast=ipar
      elast=en
      if(IA.eq.1) write(KO,25) JTOT,PSIGN(IPAR+2),en
25	format(//' Wave function for J,pi =',F5.1,a1,f5.1/)

      IF(IPC.GE.5) WRITE(KO,*) 'Have:',LX,LY,LL,S,JNN,IC,JC,JTOT
C
      IF(IPC.GE.4)WRITE(KO,*) IA,LX,LY,LL,S,JNN,IC,JC,JTOT
30    FORMAT(' Input ch.',I3,' ( nos.',6i3,2f4.1,') ')
     
      if(IC>1.or.JC>0.) then
	write(6,*) 'This program not written for core excitation!'
	stop
	endif
C
40      RV(1) = 0.0
        DO 50 IX=1,NRXY
	RV(IX) = (IX-1)*RSTEP
	read(NFL,*) IIX
50	READ(NFL,*) (WF(IX,IY,IA,ITYTR),IY=1,NRXY)
	XN = 0.0
        
        DO 53 IX=1,NRXY
           x2=(float(ix)*rstep)**2
        DO 53 IY=1,NRXY
           y2=(float(iy)*rstep)**2
           RHO2 = (ix*rstep*CX)**2 + (iy*rstep*CY)**2

	if(IVERT.eq.0) then
           T=WF(IX,IY,IA,ITYTR)**2 * RSTEP**2
           XN = XN +  T
!! AMoro addition (26/01/05)
           XRM = XRM + T * RHO2
           rho2=rho2+ T*RHO2
        endif
53	if(IVERT.gt.0) XN = XN + WF(IX,IY,IA,ITYTR) * 
     X		RV(IX)*RV(IY)*RSTEP**2 *FPI

 	read(NFL,*) PNORM
      IF(IPC.GE.4) WRITE(KO,*) 'Input partial integral:',real(XN)
      IF(IPC.eq.3)WRITE(KO,80) IA,LX,LY,LL,S,XN
80    FORMAT(' Input ch.',I3,' ( lx,ly,L,S=',4i3,') norm =',f9.5)
      WN = WN + XN
110   CONTINUE
120 	read(NFL,*) TNORM
        XRM = sqrt(XRM/WN)
      IF(IPC.GE.3) WRITE(KO,*) 'Input summed integrals:',real(TNORM)
      IF(IPC.GE.3) WRITE(KO,*) 'Input rms hyperradius:',real(XRM)
!      write(*,'(/," Sqrt[<rho^2>]=",1f8.4,/)') sqrt(x2avg+y2avg)
200   CONTINUE
      indx = IWF-1
      
      if (ijpi.eq.0) then
         call dens(indx)
      else if (quais(ijpi)>0) then
         write(99,*)'Calling dens for ijpi=',ijpi
         call dens(indx)
      else
         write(99,*)'For ijpi=',ijpi,'skipping dens'
      endif
      enddo

      
c************************************************
c      CHECK FOR READING THE ENERGIES PROPERLY
c************************************************
      write(99,*) "---------------------------------------"
      write(99,*) 'CHECK FOR READING THE ENERGIES PROPERLY'
      do 500  in=0,NWF-1
      write(99,*)'in',in, 'ens(in)', ens(in)
 500  continue
      write(99,*) "-------------------------------"
   

      RETURN
      END



c *** ----------------------------------------------------------
c *** Read 3B wave function in hyperspherical coordinates
c *** Uses format from Efaddy code
c *** ----------------------------------------------------------
      subroutine readrowf(onlyelastic)
      use trdens
      use rhowf
      use wfns
      use scattering
      USE realloc_mod
      use time
      implicit none  
     
      logical:: onlyelastic
      real*8:: co,pr,en
      integer :: iread,i,j,m,kout,ncomp,iwf,ima,nchan,pty
      character*40 filewf
      character*80 line
      character*20 head1
      character*7 word
      integer:: ir,k,lx,ly,ll,ichan,ip
      real*8::dxy
      real*8:: aux
      real*8:: y(0:20)
      real*8:: ener,wfnorm
      integer njpi,nma,jnn
     
      real*8:: jlast=-1
      real*8::iparlast=-2,elast=0.0
      integer:: ijpi=-1,ipar
      integer:: jimax=0,jfmax=0

!!!
!!!      integer lxmax, lymax
      integer nchmax
!!!
     
      namelist /wfrho/ filewf,dxy,smallchan,lxmax,lymax
      
      kout=6
      nma=0
      mllmx=0
      smallchan=0d0
      
!!! PROVISIONAL     
      ca=1.d0/(2.d0*sqrt(0.95637d0/20.7502d0)) !6He
!!!

      
     
      write(*,'(/,A,/)')'Dummy read to setup dimensions:'
      do iread=1,2
         cumtime1=0
         cumtime2=0
         rewind(18)
c Initialize some variables 
         iwf=0
         jlast=-1; iparlast=-2; elast=0d0
         ijpi=-1
         jimax=0; jfmax=0

      if (iread.eq.2) then
         write(99,*)'+Allocating ffcore,ffcth,ffval,tnq...'
         write(99,*)' with mll,nwf,nma=',mll,nwf,nma
         allocate(ffcore(mdelta,0:mll,0:nwf))
         allocate(ffcth(nangles,0:mll,0:nwf))
         allocate(ffval(mdelta,0:mll,0:isnx,0:dmax,0:nwf))
         allocate(ffvth(nangles,0:mll,0:isnx,0:dmax,0:nwf))
!!! CHANGE 28/11/05
!!         allocate(tnq(8,nma,2))
         allocate(tnq(8,nma,nwf))
         allocate(tqn(8,nma,nwf))
         allocate(na(0:nwf))     
         allocate(ens(0:nwf))
         allocate(lvalmax(0:nwf))
         allocate(chnorm(0:nwf,nma))
         nrxy=nint(nrho*drho/rstep)+1
         write(99,*)'+Allocating wf with nrxy=',nrxy
         allocate(wf(nrxy,nrxy,nma,2))
         ffcore=0d0
         ffval=0d0
         ffcth=0d0
         ffvth=0d0
         na=0d0
         ens=0d0
      endif   

50       dxy=0d0
         smallchan=0.0
         lxmax=0
         lymax=0
         read(18,nml=wfrho) !main input: mst.in
         write(kout,*)' '
         write(kout,*)'- Reading 3B wf from file: ',filewf
!         write(*,nml=wfrho)
!         write(*,*)'iread=',iread
         if (dxy.lt.1e-6) goto 998
         rstep=dxy
         open(11,file=filewf,status='old',err=90)


      rewind(11)
      nchmax=0
10    read(11,385,end=995) drho,nrho,nchan,JTOT,PTY,en
385   format('#',f8.4,2i4,f4.1,i2,f10.5)
!      write(*,*) drho,nrho,nchan,JTOT,PTY,en
      nchinc=0
      
      if (nchan.gt.nma) nma=nchan
      if(allocated(rr)) deallocate(rr)
      allocate(rr(nrho))
      read(11,'(a)') line
!      write(*,*)line
      read(11,'(a)') line
!      write(*,*)line
      if (onlyelastic.and.(iwf.eq.1)) then
         write(*,*)!'========================================='
         write(*,*)'Onlyel=T =>Skipping non-elastic components'
         write(*,*)!'========================================='
         goto 995
      endif
      
      IWF=IWF+1      
      IF(IWF.EQ.1)ITYTR=elastic
      IF(IWF.GT.1)ITYTR=inelastic

      write(99,*)'+Allocating memory for ',nchan,' channels'
      allocate(qn(nchan,10))
      allocate(wf0(nchan,nrho))
!      allocate(wfro(njpi,nchan,nrho))      
      allocate(chknorm(nchan)) !channel normalization
      chknorm=0d0
      wf0=0.d0
      rr(1)=0.d0     

      do i=1,nchan,10 !K  L Sx lx ly jp Iy Ix Iz JT
         ncomp=10
         if (nchan-i.lt.10) ncomp=nchan-i+1
!         write(*,*)'Reading wfrom i=',i,'to i=',i+ncomp-1
         do j=i,i+ncomp-1
            read(11,200) word,(qn(j,m),m=1,7)
!           read(11,200) word,(qn(j,m),m=1,10)
!            write(*,'(11i3)') j,(qn(j,m),m=1,7)
 200        format(1a7,10i3)
            k=qn(j,1)
            lx=qn(j,4)
            ly=qn(j,5)
            ll=qn(j,2)
            S=qn(j,3)
            jnn=qn(j,6)
!           jt=qn(j,10)
            if((lxmax.eq.0).and.(lymax.eq.0)) then
               nchinc=nchinc+1
            else if (lx.le.lxmax.and.ly.le.lymax) then
               nchinc=nchinc+1
            endif
!            if ((lxmax.gt.0).and.(lymax.gt.0).and.
!     &      (lx.le.lxmax.and.ly.le.lymax)) nchinc=nchinc+1
            write(99,121) j,k,lx,ly,ll,s 
121        format("Chan:",1i3,": (k,lx,ly,L,S)=",5i3)

            if (max(lx,ly)>mllmx) mllmx=2*max(lx,ly)
            IPAR = (-1)**(LX+LY)

            if ((jlast.ne.jtot).or.(ipar.ne.iparlast).or.
     &      (elast<0.and.abs(elast-en)>0.001)) then
             ijpi=ijpi+1
!             write(*,*)'en,elast,jlast=',en,elast,jlast
!             write(*,916)ijpi,JTOT,PSIGN(IPAR+2),EN
916          format('  Component: J,PI,=',i3,f5.1,a1,' at',f8.3,' MeV')
            endif
 
            if (ijpi>0) then
               jfmax=max(jnn,jfmax)
            else
               jimax=max(jimax,jnn)
            endif
            dmax=jimax+jfmax
            mll=dmax+isnx
                 
            jlast=jtot
            iparlast=ipar
            elast=en
         enddo
         
         do j=1,nrho
            read(11,395) rr(j),(wf0(m,j),m=i,i+ncomp-1)
395         format(f7.3,1p,10e12.4)
            chknorm(i:i+ncomp-1)=chknorm(i:i+ncomp-1)+
     &                      wf0(i:i+ncomp-1,j)**2*drho
         enddo 
         read(11,*) line ! This line should contain "&" character
      enddo !ichan
      wfnorm=0d0
      do ichan=1,nchan
         write(99,*)'  Ch.',ichan,'  Norm=',chknorm(ichan)
         wfnorm=wfnorm+chknorm(ichan)
      enddo
      write(99,*)'   State:',iwf,'  Total norm=',wfnorm
      write(*,'(2x,"o State ",i3," at E=",1f6.3," MeV  J/pi="
     &,F5.1,a1," :", i3," chans.  Norm=",1f8.5 )')
     & iwf,en,jtot,psign(ipar+2),nchan,wfnorm

      if (nchinc.gt.nchmax) nchmax=nchinc
!!      write(*,*)'nma=',nma,' nchmax=',nchmax 
!!! TEST OJO
      nma=nchmax
!!! 
ccc     call wfrho2p(nchan)  !NOT USED NOW

ccc      goto 111

      if (iread.eq.2) then
!         write(*,*)'ijpi=',ijpi
         if (ijpi.eq.0) then
            call rho2xy(drho,dxy,nchan,iwf)
            call groupk(nchan,iwf)
            call dens(iwf-1)
         else if (quais(ijpi)>0) then
            call rho2xy(drho,dxy,nchan,iwf)
            call groupk(nchan,iwf)
            call dens(iwf-1)
         else
            write(*,*)'For ijpi=',ijpi,'skipping dens'
         endif
      endif
ccc111   continue

      if (allocated(qn)) then
         deallocate(qn)
      else
         write(*,*)'qn not allocated!'
         stop
      endif


      if (allocated(wf0)) then
         write(99,*)'deallocating wf0'
         deallocate(wf0)
      else
         write(*,*)'wf0 not allocated!'
         stop
      endif
      
       if (allocated(chknorm)) then
          write(99,*)'-Deallocate chknorm'
         deallocate(chknorm)
      else
         write(*,*)'chknorm not allocated!'
         stop
      endif
      goto 10 !iwf
995   close(11)
      goto 50 !try to read a new WF file

998   nwf=iwf
      
      write(*,'(/)')
      write(*,*)' - Read',nwf,'states'
      write(*,*)' - J/pi combinations:',ijpi
      write(*,*)' - Max. number of channels:',nma
      write(*,*)' - Max. number of channels actually included:',nchmax
      write(*,'(/)') 
      enddo !iread=1,2
      call write_time_stat()
      
      return
90    write(*,*) "Could not open file ",filewf
      stop
      end subroutine readrowf
      



c *** ------------------------------------------------------------
c *** Transforms the 3B wf from hyperradial coord. to Jacobi coord.
c *** {rho,alpha}-> {r,R} 
c *** ------------------------------------------------------------
      subroutine rho2xy(hro,dxy,nchan,iwf)
        use rhowf
        use wfns
        use scattering
        use time
        implicit none
        integer(4)::ic2,crate,cmax
        integer(4)::ic4,ticks,secs
        integer::i,j,n,njpi,nma,iwf,ichan,ix,iy,k,lx,ly,nchan,ll,jnn
        real*8::tw,x,y,ro,alphaxy,fival,wfalpha
        real*8::ax,ay,hx,hy,hro,ang,chh,dxy,cnorm,wfnorm,T
        real*8::rx,ry,xrm
        real*8,pointer::w(:)         
        tw=0.d0
        ang=0.d0

        ax=sqrt(m2*m3/(m2+m3))
        ay=sqrt(m4*(m2+m3)/(m2+m3+m4))
      
        rstep=dxy ! compatibility with the rest of the code
        nrxy=nint(nrho*hro/dxy)+1
        call system_clock(count=ic2,count_rate=crate,count_max=cmax)
        write(*,'(/)')
        print*,'   rho2xy: {rho,alpha} => {r,R}'
        write(99,*)'hro,nchan,nrxy,rstep,iwf',hro,nchan,nrxy,rstep,iwf
         
        if (allocated(wfk)) deallocate(wfk)
        write(99,*)'rho2xy:allocating wfk with nrxy,nchan=',nrxy,nchan
        allocate(wfk(nrxy,nrxy,nchan,2))
       
        IF(IWF.EQ.1)ITYTR=elastic
        IF(IWF.GT.1)ITYTR=inelastic
        wfnorm=0d0
        xrm=0d0
        do ichan=1,nchan
           if (smallchan>0.and.chknorm(ichan)<smallchan) cycle
!            if (smallchan>0.and.chknorm(ichan)<smallchan) then
!            write(*,*)'(skipping remaining channels)'
!            cycle
!                write(99,*)'rho2xy: skipping channel',ichan
!               cycle
!            endif
!            write(*,'(a1,$)')'.'
           write(*,'(i3,",",$)')ichan
            cnorm=0d0
            k=qn(ichan,1)
            lx=qn(ichan,4)
            ly=qn(ichan,5)
            ll=qn(ichan,2)
            S=qn(ichan,3)
            jnn=qn(ichan,6)
            write(99,*)ichan,k,lx,ly,ll,jnn
            w=>wf0(ichan,:)
            do ix=1,nrxy
              do iy=1,nrxy
                 rx=dble(ix)*rstep
                 ry=dble(iy)*rstep
                 x=ax*rx
                 y=ay*ry
                 ro=sqrt(x**2  + y**2)
                 tw=fival(ro,rr,w,nrho,ang)     
                 if (y<1e-6) y=1e-6
                 alphaxy=atan(x/y)
                 wfk(ix,iy,ichan,itytr)=tw*wfalpha(k,lx,ly,alphaxy)*
     &            sqrt(ax*ay/ro)
                 T=WFK(IX,IY,ichan,ITYTR)**2 * RSTEP**2
                 cnorm=cnorm+T
                 WFNORM=WFNORM+T                    
                 XRM = XRM + T * ro*ro
                 end do !iy
              end do !ix
!              write(99,*)'Channel:',ichan,' Norm=',cnorm
              write(99,122) j,k,lx,ly,ll,s,cnorm
122           format("Chan:",1i3,": (k,lx,ly,L,S)=",5i3,' Norm=',1f6.4)
           end do ! ichan
           XRM = sqrt(XRM/WFNORM)
          
           call system_clock(count=ic4)
           ticks = ic4-ic2
!           ticks = mod(ticks+cmax4, cmax4)! reset neg. numbers
           secs = float(ticks)/float(crate)

           write(*,*)'=> ',secs,' secs'
           cumtime1=cumtime1+secs
!           write(*,'(/)')
!           write(*,*) 'State iwf=',iwf,' Norm=',wfnorm
!           write(*,*)  'Input rms hyperradius:',real(XRM)
           write(*,125) wfnorm,xrm
125        format(5x,"Norm=",1f8.5," RMS hyperradius=",1f8.3," fm")
           write(*,'(/)')
           return
      end subroutine rho2xy


c ***
c *** Group K-channels with the same lx,ly,..
c ***
      subroutine groupk(nchan,iwf)
        use realloc_mod
        use rhowf
        use wfns
        implicit none
        logical comp5
        integer lx1,ly1,ll1,st1,jnn1,k
        integer lx2,ly2,ll2,st2,jnn2
        integer ich,nchan,ima,ia,iwf
        real*8,allocatable:: wfaux(:,:,:)

        IF(IWF.EQ.1)ITYTR=elastic
        IF(IWF.GT.1)ITYTR=inelastic
        
        if (allocated(wfaux)) deallocate(wfaux)
       
        
        write(99,*)'+Allocating wfaux,chnorm with',nrxy,nchan
!        write(*,*)'+Allocating wfaux,chnorm with',nrxy,nchan
        allocate(wfaux(nrxy,nrxy,nchan))
!        allocate(chnorm(nchan))

        do 100 ich=1,nchan
           k=qn(ich,1)
           ll1=qn(ich,2)
           ST1=qn(ich,3)
           lx1=qn(ich,4)
           ly1=qn(ich,5)
           jnn1=qn(ich,6)
           if (ich.eq.1) then
!              write(*,*)'groupk(2):nchan,iwf',nchan,iwf
!              write(*,*)k,ll1,st1,lx1,ly1,jnn1
              wfaux(:,:,1)= wfk(:,:,1,itytr)
              tnq(1,1,iwf)=lx1
              tnq(2,1,iwf)=ly1
              tnq(3,1,iwf)=ll1
              tnq(4,1,iwf)=st1
              tnq(5,1,iwf)=jnn1
!              write(*,*)tnq(1:5,1,iwf)
              ima=1
              write(99,'("new chan:",6i3)')ima,lx1,ly1,ll1,st1,jnn1
              chnorm(iwf,1)=chknorm(1)
              goto 100
           endif

           write(99,*)'----------------------'
           write(99,*)'kchan:',ich,lx1,ly1,ll1,st1,jnn1
!           write(*,*)k,lx,ly,ll,jnn
150           format(5x,5i4)
           do ia=1,ima
              lx2 = TNQ(1,ia,iwf) 
              ly2 = TNQ(2,ia,iwf)
              ll2 = TNQ(3,ia,iwf)
              ST2  = TNQ(4,ia,iwf)
              jnn2 =TNQ(5,ia,iwf)
             
!              if((lx1.eq.lx2).and.(ly1.eq.ly2).and.(ll1.eq.ll2).and.
!     &      (st1.eq.st2).and.(jnn1.eq.jnn2)) then
              if(comp5(lx1,ly1,ll1,st1,jnn1,lx2,ly2,ll2,st2,jnn2))then
                 wfaux(:,:,ia)= wfaux(:,:,ia)+ wfk(:,:,ich,itytr)
                 write(99,*)'coincides with ia=',ia
!!!!
              chnorm(iwf,ia)=chnorm(iwf,ia)+chknorm(ich)
!!!!
                 goto 100 !go to next ichan
              endif
           enddo
           ima=ima+1
           wfaux(:,:,ia)= wfk(:,:,ich,itytr)
           tnq(1,ima,iwf)=lx1
           tnq(2,ima,iwf)=ly1
           tnq(3,ima,iwf)=ll1
           tnq(4,ima,iwf)=st1
           tnq(5,ima,iwf)=jnn1
           chnorm(iwf,ima)=chknorm(ich)
           write(99,'("new chan:",6i3)')ima,lx1,ly1,ll1,st1,jnn1
100        continue

        write(99,'(i3," K-channels grouped into",i3," beta channels")')
     &    nchan,ima
        na(iwf)=ima
!        tqn(1:5,1:ima,iwf)=tqn(1:5,1:ima,itytr)
        do ia=1,ima
           lx1 = TNQ(1,ia,iwf) 
           ly1 = TNQ(2,ia,iwf)
           ll1 = TNQ(3,ia,iwf)
           ST1 = TNQ(4,ia,iwf)
           jnn1 = TNQ(5,ia,iwf)
           write(99,120) ia,lx1,ly1,ll1,st1,chnorm(iwf,ia) 
120        format(8x,"Ch.",1i3,2x,"(lx,ly,L,S)=",4i3,"  Norm=",1f8.6)
        enddo
        
        wf(:,:,1:ima,itytr)= wfaux(:,:,1:ima)
!        deallocate(chnorm)
        return
      end subroutine groupk





c ***
c *** comp5=T if i1=j1, i2=j2, etc
c ***
      function comp5(i1,i2,i3,i4,i5,j1,j2,j3,j4,j5)
        implicit none
        integer i1,i2,i3,i4,i5,j1,j2,j3,j4,j5
        logical comp5
        if((i1.eq.j1).and.(i2.eq.j2).and.(i3.eq.j3).and.
     &      (i4.eq.j4).and.(i5.eq.j5)) then
          comp5=.true.
!          write(*,*)'comp5=T'
!          write(*,150)i1,i2,i3,i4,i5
!          write(*,150)j1,j2,j3,j4,j5
150           format(5x,5i4)
        else
!           write(*,*)'comp5=F'
!           write(*,*)i1,i2,i3,i4,i5
!           write(*,*)j1,j2,j3,j4,j5
           comp5=.false.
        endif
      end function comp5
        




c -------------------------------------------------------------------
c *** Calculate Fourier transform of wf0 functions
c     NOT NEEDED
c *** Fourier trasnform of wf rho -> hypermomentum
c --------------------------------------------------------------------
      subroutine wfrho2p(nchan)
        use rhowf
        implicit none 
        integer mb,nrmax,ip,ir,ichan,k,nchan
        real*8,allocatable:: wfp(:,:)
        parameter (mb=50,nrmax=5000)
        real*8::  p,dp,pmax,aux,pr,bessj
        real*8::  normr,normq,normt
        pmax=10
        dp=pmax/nrho
        write(*,*)'dp=',dp 

! K  L Sx lx ly jp Iy Ix Iz JT
     
!         lx=qn(ichan,4)
!         ly=qn(ichan,5)  
      
      allocate(wfp(nchan,nrho))
      normt=0d0
      do ichan=1, nchan
         normr=0
         do ir=1,nrho
            normr=normr+wf0(ichan,ir)**2*drho
         enddo
         normt=normt+normr
         normq=0d0
         do ip=1,nrho !linear hypermomentum
            p=ip*dp
            k=qn(ichan,1)
            aux=0d0
            do ir=1,nrho
               pr=rr(ir)*p   
               if (pr<1e-8) pr=1e-8
                aux=aux+bessj(k+2,pr)*wf0(ichan,ir)*drho *sqrt(pr)
            enddo
            wfp(ichan,ip)=aux
            normq=normq+wfp(ichan,ip)**2*dp
            if (ichan.lt.10) write(90,'(3e16.6)') p,aux
         enddo
         write(90,*)'&'
         write(*,777) k,normr,normq,normr/normq
777      format(' - K:',i3,' Normr=',1f12.6,' Normq=',2f12.6 )
      enddo
      write(*,*)'- Total norm=',normt

      end subroutine wfrho2p
