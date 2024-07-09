       program mstamp
c*****************************************************************************
c      Calculates nucleon-nucleus scattering for a 3-body target (6He,11Li,etc)
c      using the Multiple Scattering of the T-matrix formalism, described in:
c 
c      R. Crespo and R. C. Johnson, Phys. Rev. C 60, 034007 (1999)
c      R.Crespo, I.J.Thompson, A.A.Korsheninnikov, Phys.Rev. C66, 021002 (2002)
c      R.Crespo, A.M.Moro, I.J.Thompson, Nucl.Phys. A771, 26 (2006)
c      R.Crespo, I.J.Thompson, A.M.Moro, Phys.Rev. C 74, 044616 (2006)
c******************************************************************************
	use wfns
	use trdens
        use scattering
        use bt
        use constants
        
      	implicit real*8(a-h,o-z)
       

        logical:: dry,onlyel,ifrho
        integer ::jnn,icb,party,ttype,typev,typec,kinem
        real*8:: massp,mtclus,spin
!       real*8:: masst
        character*40::word
        
        

        interface
          function wxxi(qrr,qcore,tnct,noangs,noangs2)
            integer::noangs,noangs2
            real*8::qrr,qcore(:)
            complex*16::tnct(:)
            complex*16::wxxi
          end function wxxi
        end interface


        namelist /mst/qmax,tlab,thmin,thmax,dth,dry,
     &                onlyel,ifrho,kinem,rel,coulmst
        namelist /proj/ massp,zp,jp
        namelist /targ/ masst,zt,ncl,inelcb,nustates,quais,irho,nwfmax,
     &                  lxymax
        namelist /jpi/ jnn,party,icb
        namelist /kapas/ k0000,k1100,k0111,k1120,k1121,k1122
        namelist /tclus/ mtclus,spin,ztclus,ttype,sigtot,beta,alpha
        namelist /quad1/ qmaxr,quin,mquadi,mquado
        namelist /quad2/ qmaxrd,quind,mquadid,mquadod 
        namelist /quad3/ rmaxr,rin,mrquadi,mrquado
 

        write(*,*)'MSTamp program. Version 1.1'
      
       
      	call logfac(lfact)
        call setconstants()

        open(unit=18,file='mst.in',status='old')
        open(unit=99,file='mst.log',status='unknown')
!        open(unit=10,file='output.data',status='unknown')
!        open(unit=15,file='gammin.data',status='old') 

!        open(20,file='trcmall.den',status='unknown')
!        open(21,file='trvalall.den',status='unknown')
!        open(23,file='trval.den', status='unknown')
!        open(24,file='trcm.den', status='unknown')
!        open(25,file='shakeoff.den', status='unknown')

!        open(unit=35,file='tnn.data', status='unknown')
!        open(unit=36,file='tnnon.data', status='unknown')
!        open(unit=40,file='dens.data', status='unknown')

!        open(unit=52,file='M0000',status='unknown')
!        open(unit=53,file='M0111',status='unknown')
!!       open(unit=54,file='M1011',status='unknown')
!        open(unit=55,file='M1100',status='unknown')
!!       open(unit=56,file='M1111',status='unknown')
!        open(unit=57,file='M1120',status='unknown')
!        open(unit=58,file='M1121',status='unknown')
!        open(unit=59,file='M1122',status='unknown')

!        open(unit=65,file='doublexs.data',status='unknown')
!        open(unit=70,file='xsection11.data',status='unknown')
!        open(unit=73,file='xs11soff.data',status='unknown')
!        open(unit=75,file='xsection9.data',status='unknown')
!        open(unit=71,file='xsecthex.data',status='unknown')
!        open(unit=96,file='xsecex.data', status='unknown')

!        open(unit=77, file='tncBA.data', status='unknown')
!        open(unit=80, file='tna.data', status='old')
!        open(unit=85, file='tnaint.data', status='unknown')

!        open(unit=92, file='scratch.data', status='unknown')
!        open(unit=95, file='warning.txt', status='unknown')
        


c=======================================================================
c       reads input parameters
c=======================================================================
        itnn=0;itnnav=0 !PROVISIONAL
        quais=1
        nwfmax=0 !max. number of components read (1=only elastic)
        dry=.false.
        coulmst=1
        onlyel=.false.
        ifrho=.false.
        kinem=0 !"Chew" kinematics
        read(18,nml=mst)
        select case(kinem)
        case(0) ! MST kinematics
           word="MST  "
        case(1)!Chew kinematics
           word="Chew "
        case(2)!Rihan
           word="Rihan"
        case(3)!Kujawski & Lambert
           word="Kujawski-Lambert"
        case(4,5) ! Adiabatic  (FSA)
           word="FSA"
        case default
           write(*,*)'kinem=',kinem,'not used!';stop
        end select
        write(*,'(2x,/,"** Using ",a10," kinematics **",/)')word
        if (rel) then
           write(*,*) 'with RELATIVISTIC correction'
        else
           write(*,*) 'WITHOUT relativistic correction'
        endif
        
        read(18,nml=quad1)
        read(18,nml=quad2)
        read(18,nml=quad3)
        read(18,nml=proj)
        m1=massp
        z1=zp
        read(18,nml=targ)
        if (ncl>3) then
           write(*,*) 'Only up to three clusters allowed';stop
        endif
        if(nustates>50) then
            write(*,*) 'Too many states. Increase dimension of quais'
            stop
         endif

        read(18,nml=kapas)
        kapa(0,0,0,0)=k0000
        kapa(1,1,0,0)=k1100
        kapa(0,1,1,1)=k0111
        kapa(1,1,2,0)=k1120
        kapa(1,1,2,1)=k1121
        kapa(1,1,2,2)=k1122
 
       

c *** memory allocation here!!
        write(*,'(2x,"Reading",1i2," clusters")') ncl
        ttype=0
        read(18,nml=tclus) ! valence 1 
!        write(*,nml=tclus) 
        m2=mtclus
        z2=ztclus
        j2=spin
        typev=ttype
        ttype=0
        read(18,nml=tclus) ! valence 2
!         write(*,nml=tclus)
        m3=mtclus
        z3=ztclus
        j3=spin
        typev=ttype
        ttype=0 !cluster T-matrix calculated on-the-fly
        read(18,nml=tclus) ! core
!         write(*,nml=tclus)
        m4=mtclus
        z4=ztclus
        j4=spin
        typec=ttype
!        close(18)
        

 
!!$        do icl=1,ncl
!!$           read(15,nml=tclus)  ! read t-cluster 
!!$           if ((ztclus.eq.0).or.(ztclus.eq.1)) then
!!$              
!!$        call tntensor()
!!$           else
!!$              call tncore()
!!$           endif
!!$        enddo
        
 



c OLD INPUT-----------------------------------------
!        read(15,*)qmax
!        read(15,*)tlab
!        read(15,*)m1,m2,m3,m4
!        read(15,*)s2,s3
!        read(15,*)thmin, thmax,dth
!        read(15,*)qmaxr,quin,mquadi,mquado
!        read(15,*)qmaxrd,quind,mquadid,mquadod 
c        read(15,*)irho
c        read(15,*)itnn,itnnav
c        read(15,*)rhomx,steprho
c        read(15,*)inelcb
!        read(15,*)rmaxr,rin,mrquadi,mrquado
c        read(15,*) nustates,itstates 
c        read(15,*)((jnnpin(i),partyin(i),istatcb(i)),i=1,nustates)
c        read(15,*)k0000,k1100,k0111,k1120,k1121,k1122 
c------------------------------------
        zz = cmplx(0.d0,1.d0)
        mn = 939.5731
        hbarc = 197.3289
        s1 = 0.5
        iqmx = mquadi + mquado
        iqmxd = mquadid + mquadod
        mq = iqmx
        mdelta = iqmxd
        dq=0.1
        
        qmxnn = qmaxr
        nq=qmxnn/dq
        nangles = (thmax-thmin)/dth + 1
        call gauss3(qmaxr,quin,mquadi,mquado,radxis,radwt)
        call gauss3(qmaxrd,quind,mquadid,mquadod,radxisd,radwtd)
        call gauss3(rmaxr,rin,mrquadi,mrquado,radxisr,radwtr)
c        if (iqmx.gt.mrho.and.iqmxd.gt.mrho)then
c        write(10,*)'iqmx gt mrho'
c        stop
c        endif

        m234 = m2+m3+m4
        m1234 = m234+m1
        m23=m2+m3
        m34=m3+m4
        m14=m1+m4
        m12=m1+m2
        amass = m234
! AMoro (5/10/04)
!        muna = amass/(amass+massp)*mn
!       mucore=amu*massp*m234/(m234+massp)
!        muv=amu*massp*m2/(m2+massp)  
!        mucore = mn*m4/(m4+1.) !Changed 19/03/2003
        muna=m1*masst/(m1+masst)*amu !MeV
        mucore=m1*m4/(m1+m4)*amu
        muv=m1*m2/(m1+m2)*amu
       
       

        if (rel) then
c  Mandelstam          
!           s0=(m1+m234)**2 + 2*m234*tlab/amu !in amu^2        
!           k0=amu*sqrt(trian(s0,m1**2,m234**2))/2d0/sqrt(s0)/hbarc
!           w0=(dsqrt(s0)-m1-m234)*amu

           s0=(m1+masst)**2 + 2*masst*tlab/amu !in amu^2        
           k0=amu*sqrt(trian(s0,m1**2,masst**2))/2d0/sqrt(s0)/hbarc
           w0=(dsqrt(s0)-m1-masst)*amu
           tcm=w0
        else
!           tcm=tlab*m234/m1234
           tcm=tlab*masst/(m1+masst)
           k0 = dsqrt(2*muna*tcm)/hbarc
        endif
        write(99,*)'- In MST: k0=',k0
c energy parameters       
        mucore = amu*massp*m4/(m4+massp)

        select case(kinem) 
           case(4) !adiabatic (FSA)
              e14=tlab
              e12=tlab
              w14=tcm !In MeV
              w12=w14
              kv=k0
              kcore=k0
           case default
              e12=tlab !???*(1-m1*m34/m12/m234) !lab
              e14=tlab !???*(1-m1*m23/m14/m234) !lab
              if (rel) then
                 s12=(m1+m2)**2 + 2*m2*tlab/amu !in amu^2
                 s14=(m1+m4)**2 + 2*m4*tlab/amu !in amu^2
c Note: nu,nu14,nu12 reduces to mu, mu14, mu12 in the non-rel limit
                 nu=(s0+m1**2-masst**2)*(s0-m1**2+masst**2)/
     &        4/s0/sqrt(s0)
                 nu12=(s12+m1**2-m2**2)*(s12-m1**2+m2**2)/
     &        4/s12/sqrt(s12)
                nu14=(s14+m1**2-m4**2)*(s14-m1**2+m4**2)/
     &        4/s14/sqrt(s14)
                  
                 kv=amu*sqrt(trian(s12,m1**2,m2**2))/
     &                2/sqrt(s12)/hbarc
                 w12=(sqrt(s12)-m1-m2)*amu
                 
                 kcore=amu*sqrt(trian(s14,m1**2,m4**2))/
     &               2/sqrt(s14)/hbarc
                 w14=(sqrt(s14)-m1-m4)*amu

!                 write(*,*)'KL: s12, w12=',s12,w12
!                 write(*,*)'KL: s14, w14=',s14,w14
               else
                  w14=tcm*(1-m1*m23/m14/m234)  !cm kin energy
                  w12=tcm*(1-m1*m34/m12/m234)  

                  kv=sqrt(2*muv*w12)/hbarc         
                  kcore = dsqrt(2*mucore*w14)/hbarc

                  nu=muna
                  nu14=mucore
                  nu12=muv
               endif
        end select

        kv0=kv
        kc0=kcore
        



        write(*,'(/,"The following clusters are read:")')
        write(*,140) 
        write(*,*)'----------------------------------------'
        write(*,142)'Projectile:',m1,jp,z1,k0,tcm
        write(*,142)'Valence  1:',m2,j2,z2,kv,w12
        write(*,142)'Valence  2:',m3,j3,z3,kv,w12
        write(*,142)'Core      :',m4,j4,z4,kcore,w14
        write(*,*)'----------------------------------------'
        
140     format(13x,'Mass  ','J     Z     k    ','Ecm (w) ')
142     format(a11,1f6.3,1f4.1,1f6.1,1f6.3,1f8.2)

        write(99,*)'+Allocating',mdelta,'pts for tvalence,tcore'
        allocate(tvalence(mdelta))
        allocate(tcore(mdelta))

        
        tvalence(1:mdelta)=(0.d0,0.d0)
        tcore(1:mdelta) = (0.d0,0.d0)
      	stepdel = qmax/ndel

c====================================================================
c check
c        call statescb(1,0,1)
c        call statescb(0,0,1)
c        call statescb(1,1,1)
c        stop
c======================================================================
c       reads and calculates the density functions in appropriate grid
c======================================================================
        write(99,*) '-> Entering readwf'
        if (dry) then
           write(*,*) '**WARNING**: DRY RUN. NO WF READ'
        else
           if (ifrho) then
              write(*,*)'readrowf...'
              call readrowf(onlyel)
           else

              call readwf(onlyel)
           endif
        endif

c===================================================================
c       calculates the scattering
c==================================================================     

c **  reads and interpolates NA ( proton/neutron - core) scattering matrix
      write(99,*) '-> Entering tncore'
      call tncore(e14,typec,kinem)
      

c ** calculates scattering by valence particles
       write(99,*) '-> Entering tntensor w12=',w12
      call tntensor(w12,typev,kinem)
     
       write(99,*) '<- Exiting  tntensor'
c **  calculates cross section
       write(99,*) '-> Entering drotat'
       call drotat()
      
       write(99,*) '-> Entering xsecinel'
c        call xsecinel(kinem)
         call xsecinelnew(kinem)
        stop
       
       write(99,*) '-> Entering test'
       if (dry) then
          write(*,*)'DRY=T => Stopping now';stop
       endif
c       call test()  
c       write(99,*) '-> Entering sigex'
c       call sigex()
c       write(99,*) '<- Exiting sigex'       
       end program mstamp




c-------------------------------------------------------------
	subroutine cmdens(in)
c-------------------------------------------------------------
	use parameters
        use wfns
	use trdens
        use scattering
        use rhowf
        implicit none
!      	implicit real*8(a-h,o-z)
        integer,save::ierr
        integer :: ia,ib,ix,iy,in,nrmax,iwf
        integer:: ll,lx,ly,idel,lyp,lxp,lllp,lll,jnnp,jnn
        real*8 :: aux
        real*8 :: cleb6,racah2,hat
        real*8 :: alg,alg1,alg2,alg3,alg4,alg5
        real*8 :: delta,deltaren,sumcm,t,trr
	real*8 ::bes(nrxy),xn(ndel)
        complex*16:: taux(mrxy,mrxy), cint2d
!	real*8 hat(i) = sqrt(real(2*I+1))

        lvalmax(in)=-1
        ierr=0
        iwf=in+1
   
	do 100 ll=0,mll
	rhocm(:,ll) = 0.0
        nrmax = mrquadi+mrquado
        write(99,*)'Entering cmdens:mllmx,lxmax,lymax=',
     &  mllmx,lxmax,lymax
        write(99,*)'na(elastic)=',na(elastic)
	
	do 60 ia=1,na(elastic)
        if(smallchan>0.and.chnorm(iwf,ia)<smallchan) cycle    
	lx = TNQ(1,ia,elastic) 
	ly = TNQ(2,ia,elastic)
	lll = TNQ(3,ia,elastic)
	S = TNQ(4,ia,elastic)
	jnn = TNQ(5,ia,elastic)
!	ic = TNQ(6,ia,elastic)
!	jc2 = TNQ(7,ia,elastic)
!	jtot2 = TNQ(8,ia,elastic)

!29/11/05
        if ((lxmax.gt.0).and.(lymax.gt.0).and.
     &      (lx.gt.lxmax.or.ly.gt.lymax)) cycle
!          write(*,*)'dens:skip channels with lx,ly=',lx,ly
!          cycle
!        endif
!29/11/05




	do 60 ib=1,na(iwf)
!           write(*,*)'cmdens:smallchan=',smallchan
!           write(*,*)'cmdens:chnorm(iwf,ib)=',chnorm(iwf,ib)
           
        if(smallchan>0.and.chnorm(iwf,ib)<smallchan) then
           write(99,*)'cmdens: skipping iwf,ib=',iwf,ib
           write(99,*)'smallchan,chnorm=',smallchan,chnorm(iwf,ib)
           goto 60
        endif
	lxp = TNQ(1,ib,iwf) 
!       RC modification here
               if (lxp.ne.lx) go to 60
	lyp = TNQ(2,ib,iwf)
	lllp = TNQ(3,ib,iwf)
	Sp = TNQ(4,ib,iwf)
!       RC modification here
               if (Sp.ne.S) go to 60
	jnnp = TNQ(5,ib,iwf)
!	icp = TNQ(6,ib,iwf)
!	jc2p = TNQ(7,ib,iwf)
!	jtot2p = TNQ(8,ib,iwf)
! 	if(ic/=icp) go to 60

        if ((lxmax.gt.0).and.(lymax.gt.0).and.
     &      (lxp.gt.lxmax.or.lyp.gt.lymax)) then
         write(*,*)'dens:skip channels with lxp,lyp=',lxp,lyp
          cycle
        endif

        if(s/=sp) go to 60
        if(lx/=lxp) go to 60
        if (mllmx.eq.0) stop 'cmdens: mllmx=0!'
	
	alg1 = hat(jnnp) * hat(ll) * hat(ly)*hat(lll)*hat(lllp)
	alg2 = (-1)**(lll+jnn+2*jnnp+s + 2*lx+ly+lyp) 
	alg3 = racah2(lll+z,jnn+z,lllp+z,jnnp+z,s+z,ll+z)
        if (abs(ll)>1) write(99,*)'alg3',alg3,'ll',abs(ll)
	alg4 = racah2(ll+z,ly+z,lllp+z,lx+z,lyp+z,lll+z)
	alg5 = cleb6(ll+z,z,ly+z,z,lyp+z,z)
		if(abs(alg5)<1e-10) goto 60
	alg = alg1*alg2*alg3*alg4*alg5 

!!$        if (abs(nustates-itstates).gt.0.0001)then
!!$        call statescb(Jnnp,lxp,lyp)
!!$        if(jnncb.eq.0)alg=0.
!!$        endif

!! AMoro (25/11/2003) Send this to log file instead of stdout
	write(99,10) in,ll,ia,ib,alg1,alg2,alg3,alg4,alg5,alg
10	format(' in,ll=',2i2,' #',2i3,':',5f8.4,f10.5)

	if(abs(alg)<1e-10) goto 60
	do 50 idel=1,ndel
 	 delta = (idel-1)*stepdel
         deltaren = delta*m23/m234
	 if(ll>0.and.idel==1) go to 50
         if (mllmx.eq.0) stop'cmdens(2): mllmx=0!'
	 call bessr(ll,deltaren,bes,nrxy,rstep)! Bessel or spherical bessel?`
	XN(idel) = 0.0
        sumcm = 0.0
        aux=0d0
        DO 20 IX=1,NRXY
        DO 20 IY=1,NRXY
!	RR = RV(IY)
	T = WF(IX,IY,IA,elastic)*WF(IX,IY,IB,ITYTR) 
	trr = T * bes(IY) * rstep**2 
        taux(ix,iy) = trr
        aux=aux + WF(IX,IY,IB,ITYTR)**2*rstep**2
!       RC modification here
!	XN(idel) = XN(idel) + trr
20	sumcm  = sumcm + trr

!       Calculates radial integral more accurately
!       sumcm =  0.0
!       do 25 i=1,nrmax
!       do 25 j=1,nrmax
!       qr = radxisr(i)
!       qrr = radxisr(j)
!       tauxin = dreal(cint2d(RV,RV,taux,qr,qrr,NRXY,NRXY,5,mrxy))
!       sumcm = sumcm + tauxin * radwtr(i)*radwtr(j)
!       if (icount.lt.1) write(*,*)tauxin, radwtr(i),radwtr(j)
! 25    continue

        sumcm = sumcm * alg    
        if (abs(aux)<1e-8) then
        write(99,*)'aux=',aux,' for in,lx,ly,lxp,lyp=',in,lx,ly,lxp,lyp
           ierr=ierr+1
           if (ierr.eq.10) then
              write(*,'(6f12.6)') WF(1:5,1:5,IB,ITYTR)
              write(*,*)'IB,ITYTR=',IB,ITYTR
              
           endif
        endif
        rhocm(idel,ll) = rhocm(idel,ll) + sumcm
50	continue
60	continue

!       RC modification here	
!	rhocm(:,ll,in) = XN(:) * alg
   	lvalmax(in)=ll        
100	continue
   	if(lvalmax(in)>=0) then
        parityf = (-1)**(lxp+lyp)
	write(20,101) in, jnnp, parityf
101	format('# CM transition density to state no.',i3, i3, i3 )
	do 105 idel=1,ndel
 	 delta = (idel-1)*stepdel
105	write(20,106) delta,(rhocm(idel,ll),ll=0,lvalmax(in) )
106	format(f8.3,3f12.8)
	write(20,*) '&' 

c         if (in.eq.3)then
c         do 110 idel=1,ndel
c         delta = (idel-1)*stepdel  
c         write (9,*) 'inel',delta,rhocm(idel,0), rhocm(idel,1)
c 110     continue
c         endif 
	endif
	return
	end

	subroutine valdens(in)
c***************************************************************************
	use wfns
	use trdens
        use scattering
        use rhowf
      	implicit real*8(a-h,o-z)
	real*8 bes1(nrxy),bes2(nrxy),xn(ndel,ndel)
	integer b,c,d,phi1,phi2
	complex*16 alg1
        complex*16 alg,sumval
        integer,allocatable::cdmax(:,:)
        integer::bmax
	hat(i) = sqrt(real(2*I+1))
	hats(i) = real(2*I+1)

 	s2 = 0.5  
 	s3 = 0.5  
!        write(*,*)'valdens: s2,s3=',s2,s3
!        if (.not.allocated(cdmax))
        allocate(cdmax(2,0:mll))
        
        iwf=in+1
	bmax=-1
        rhoval(:,:,:,:)=0d0
!        write(99,*)'valdens: mll=',mll

! Corrected 24/11/05
!        do 100 b=0,mll
	do 100 b=0,isnx
           cdmax(1:2,b)=-1
           do 100 c=0,mll
              do 100 d=0,dmax
                 if (mllmx.eq.0) stop'valdens: mllmx=0!'	
                 do 60 ia=1,na(elastic)
                    if(smallchan>0.and.chnorm(iwf,ia)<smallchan) cycle
                    lx = TNQ(1,ia,elastic) 
                    ly = TNQ(2,ia,elastic)
                    lll = TNQ(3,ia,elastic)
                    S = TNQ(4,ia,elastic)
                    jnn = TNQ(5,ia,elastic)
                    if ((lxmax.gt.0).and.(lymax.gt.0).and.
     &      (lx.gt.lxmax.or.ly.gt.lymax)) cycle
!	ic = TNQ(6,ia,elastic)
!	jc2 = TNQ(7,ia,elastic)
!	jtot2 = TNQ(8,ia,elastic)
                    iwf=in+1
                    do 60 ib=1,na(iwf)
                       if(smallchan>0.and.chnorm(iwf,ib)<smallchan)cycle
                       lxp = TNQ(1,ib,iwf) 
                       lyp = TNQ(2,ib,iwf)
                       lllp = TNQ(3,ib,iwf)
                       Sp = TNQ(4,ib,iwf)
                       jnnp = TNQ(5,ib,iwf)
                        if ((lxmax.gt.0).and.(lymax.gt.0).and.
     &      (lxp.gt.lxmax.or.lyp.gt.lymax)) cycle
!	icp = TNQ(6,ib,iwf)
!	jc2p = TNQ(7,ib,iwf)
!	jtot2p = TNQ(8,ib,iwf)
! 	 if(ic/=icp) go to 60
!	do 60 l=0,mllt
!	do 60 lp=0,mllt
                       do 60 l=0,mllmx
                          do 60 lp=0,mllmx
	
 	alg1 =  (0.,1.)**(-l-lp)
	alg2 = hats(l) * hats(lp) * hats(d)
     X          * hat(b)*hat(s)*hat(sp)
     X		* hat(lll)*hat(lllp)*hat(lx)
     x          * hat(ly) * hat(jnnp) / hat(c)
	 phi1 = 2*s - sp + nint(s2 - s3) + b
	 phi2 = lll+lllp-jnnp+c-l+lp+d+jnn+2*(d-b+Sp-lyp-lx)
	alg3 = (-1)**(phi1 + phi2)
	alg4 = cleb6(l+z,z,lx+z,z,lxp+z,z) 
	alg5 = cleb6(l+z,z,lp+z,z,c+z,z) 
	alg6 = cleb6(lp+z,z,ly+z,z,lyp+z,z)
		if(abs(alg4*alg5*alg6)<1e-10) goto 60

	alg7 = racah2(s+z,s2,sp+z,s2,s3,b+z)
	alg8 = wig9j(b+z,  d+z,   c+z,
     X               sp+z, jnnp+z,lllp+z,
     X               s+z,  jnn+z, lll+z)
	alg9 = wig9j(l+z,  lp+z,  c+z,
     X               lx+z, ly+z,  lll+z,
     X               lxp+z,lyp+z, lllp+z)
	alg =alg1*alg2*alg3*alg4*alg5*alg6*alg7*alg8*alg9 

!        if (abs(nustates-itstates).gt.0.0001)then
!        call statescb(Jnnp,lxp,lyp)
!        if(jnncb.eq.0)alg=0.
!        endif

	if(abs(alg)<1e-10) goto 60
	XN(:,:) = 0.0
	do 50 idel1=1,ndel
           if(l>0.and.idel1==1) go to 50
           delta1 = (idel1-1)*stepdel
           delta1ren = delta1*m3/m23
!!!!!           bes1(:)=0.001
           call bessr(l,delta1ren,bes1,nrxy,rstep) 

           if(lp>0.and.idel1==1) go to 48
           delta2 = (idel1-1)*stepdel
           delta2ren = delta1*m4/m234
!!!           bes2(:)=0.001
        call bessr(lp,delta2ren,bes2,nrxy,rstep)  
!       XN(idel1,idel2) = 0.0
        sumval = (0.0,0.0)
        DO 20 IX=1,NRXY
        DO 20 IY=1,NRXY
!	 RR = RV(IY)
	T = WF(IX,IY,IA,elastic)*WF(IX,IY,IB,ITYTR) 
	trr = T * bes1(IX)*bes2(IY)*rstep**2
c        if (trr>1e10) then 
c           write(99,*)'T,bes1,bes2,trr=',t,bes1(IX),bes2(IY),trr
c        endif
!20	XN(idel1,idel2) = XN(idel1,idel2) + trr
 20     sumval  =  sumval + trr
        sumval  = sumval * alg
        rhoval(idel1,b,c,d)=rhoval(idel1,b,c,d) + sumval
48      continue
50	continue
60	continue
	
!	rhoval(:,b,c,d,in) = XN(:,:) * alg
	 bmax = max(b,bmax)
   	 cdmax(1,b)=c
   	 cdmax(2,b)=d
100	continue
	 do 110 b=0,bmax
   	 if(cdmax(1,b)>=0.and.cdmax(2,b)>=0) then
	 write(19,101) in,jnnp,b
101	format('# Valence transition density to state no',
     x      i3,' jnnp:',i3,' b:',i3)
        write(19,*)'real components'
	do 105 idel1=1,ndel
 	 delta1 = (idel1-1)*stepdel
105	write(19,106) delta1,((dreal(rhoval(idel1,b,c,d)),
     X    d=0,max(2,cdmax(2,b)) ),
     X    c=0,max(2,cdmax(1,b)) )
106	format(f8.3,9f10.6)
	write(19,*) '&'

        write(19,*)'imag components'
	do 115 idel1=1,ndel
 	 delta1 = (idel1-1)*stepdel
115	write(19,116) delta1,((dimag(rhoval(idel1,b,c,d)),
     X    d=0,max(2,cdmax(2,b)) ),
     X    c=0,max(2,cdmax(1,b)) )
116	format(f8.3,9f10.6)
c	write(22,*) '&'

	endif
110	continue
        deallocate(cdmax)
	return
	end


      subroutine dens(in)
c***********************************************************************
c     Calculates halo density and core cm density distribution
c     in the momentum space grid for the scattering and in theta grid
c------------------------------------------------------------------------
      use parameters
      use wfns
      use scattering
      use trdens
      use constants
      use time
!      implicit real*8 (a-h,o-z)
      implicit none
      integer(4)ic1,crate,cmax,ic2,ticks,secs
      integer ll,i,iq,ith,ic,id,ib,in
      real*8::ren,q,th,thrad,qq,a1,a2,cthna
      real*8 ::  t, p0
      integer ierr
      complex*16 wxxi

! CORE:    from rhocm in standard grid evaluates:
!               *ffcore* in q space and *ffcth* it theta space
! VALENCE: from rhoval in standard grid evaluates:
!               *ffval* in q space and *ffvth* it theta space

      call system_clock(count=ic1,count_rate=crate,count_max=cmax)
      write(99,*)'Calling dens for iwf=',in+1
      write(*,120)in+1
120   format(/,3x,"- Calling dens for  IWF=",i3,$)
      write(99,*)'dens:allocating rhoval,rhodens with'
      write(99,*)'ndel,mll,dmax=',ndel,mll,dmax
      allocate(rhoval(ndel,0:isnx,0:mll,0:dmax))
      allocate(rhocm(ndel,0:mll))
      rhocm=0d0
      rhoval=0d0
      
      write(99,*)'calling cmdens:in,mllmx=',in,mllmx
      call cmdens(in) 

!      if (in.eq.0) rhocm0=rhocm(1,0)
       
      if (irho.eq.0)then
        allocate(rhoaux(ndel))
        allocate(qdelta(ndel))
        do 20 ll=0,mll
        do 5 i=1,ndel
        rhoaux(i) = rhocm(i,ll)
        qdelta(i) = (i-1)*stepdel
        if(in.eq.0.and.ll.eq.0)then
!!        ren = 1/rhocm(1,0,0)
           ren = 1/rhocm(1,0)
           rhoaux(i) = ren * rhoaux(i)
        endif
 5      continue   
        
	do 10 iq=1,mdelta
	q = radxisd(iq)
        ffcore(iq,ll,in) = wxxi(q,qdelta,rhoaux,ndel,ndel)
 10	continue

        
        do 15  ith = 1,nangles
        th = thmin + dth*dble(ith-1)
        thrad = th*pi/180.
        cthna = cos(thrad)
        qq = sqrt(2*k0*k0*(1-cthna))
        ffcth(ith,ll,in)  = wxxi(qq,qdelta,rhoaux,ndel,ndel)
 15     continue

 20     continue
        write(99,*)'calling valdens:in,mllmx=',in,mllmx
 	call valdens(in)
        write(99,*)'leaving valdens'
        write(99,*)'+Allocate rhoauxp with ndel,mdelta=',ndel,mdelta
        allocate(rhoauxp(max(ndel,mdelta)),STAT=ierr)
        if (ierr/=0) stop 'Could not allocate rhoauxp'
        do 30 ic=0,mll
        do 30 ib=0,isnx
        do 30 id=0,dmax        
        do 35 i=1,ndel
        rhoaux(i) = dreal(rhoval(i,ib,ic,id))
        rhoauxp(i) = dimag(rhoval(i,ib,ic,id))
        qdelta(i) = (i-1)*stepdel
 35     continue 
    
        
	do 40 iq=1,mdelta
	q = radxisd(iq)
        a1 = wxxi(q,qdelta,rhoaux,ndel,ndel)
        a2 = wxxi(q,qdelta,rhoauxp,ndel,ndel)
        ffval(iq,ic,ib,id,in) = a1 + zz * a2
!        write(99,*)'ffval',iq,ic,ib,id,in, ffval(iq,ic,ib,id,in)
 40	continue


        do 45  ith = 1,nangles
        th = thmin + dth*dble(ith-1)
        thrad = th*pi/180.
        cthna = cos(thrad)
        qq = sqrt(2*k0*k0*(1-cthna))
        a1 = wxxi(qq,qdelta,rhoaux,ndel,ndel)
        a2 = wxxi(qq,qdelta,rhoauxp,ndel,ndel)
        ffvth(ith,ic,ib,id,in)= a1 + zz * a2
!        write(99,*)'ffvth:',ith,ic,ib,id,in,ffvth(ith,ic,ib,id,in)
 45     continue

 30     continue

      write(24,*)'Valence and core densities in momentum space in fm-1'
      do 60 iq=1,mdelta
	q = radxisd(iq)
        write(24,155)q, dreal(ffval(iq,0,0,0,0))*sqrt(2.),ffcore(iq,0,0)
 60   continue

      write(24,*)'Valence and core densities in momentum space in Gev/c'
       do 65  ith = 1,nangles
       th = thmin + dth*dble(ith-1)
       p0=hbarc*k0/1e3                 !linear momentum in GeV/c
       t=-(2*p0*sin(th*pi/180./2))**2 
       write(24,155)-t,dreal(ffvth(ith,0,0,0,0))*sqrt(2.),ffcth(ith,0,0)
  65   continue  

      else 
        write(99,*)'irho gt 0'
        stop
      end if


106	format(f8.3,9f8.3)
155     format(f8.3,2f10.6)
166     format(f8.3,e10.6)
      write(99,*)'-Deallocate rhoaux,rhoauxp,qdelta,rhoval,rhocm'
      deallocate(rhoaux,rhoauxp)
      deallocate(qdelta)
      deallocate(rhoval,rhocm)
      call system_clock(count=ic2)
      ticks = ic2-ic1
      secs = float(ticks)/float(crate)
      cumtime2=cumtime2+secs
      write(*,'("=>",i4," secs",/)') secs
      return
      end



       subroutine ptheta(in)
c************************************************************************
       use wfns
       use scattering
       use legendre
       use factorials
       implicit real*8(a-h,o-z)
       real*8, allocatable::aux(:,:)

       allocate(pleg(nangles))
       jnnp = TNQ(5,na(in+1),in+1)
       mnnt = 2*jnnp + 1
       NJ = 2*jnnp
       MJ = 2*NJ + 1

       allocate(aux(nj,mj))
       do 50 nc = 0, 2*jnnp
       nct = 2*nc + 1
       do 55 l = 1, nct
       ngama = - nc + (l-1)
       sum1 = 0.d0 
       do 100 i=1,mnnt
       mnn = - jnnp + (i-1)
       do 150 j=1,mnnt
       mnnp = - jnnp + (j-1)
       alg1 = cleb6(jnnp+z, -mnn+z, jnnp+z, mnnp+z, nc+z, ngama+z)
       phase = - mnn - 2*mnnp + (2*jnnp)
       alg1 = alg1 * (-1)**phase 
       sum1 = sum1 + alg1
 150   continue
 100   continue
       aux(nc,ngama) = sum1
       faux = exp((flog(nc-ngama+1) -flog(nc +ngama+1))*0.5) 
 55    continue
 50    continue

       do 200 ith = 1,nangles
       th = thmin + dth*(ith-1)
       thrad = th*pi/180.
       cthna = cos(thrad)
       call PLM(cthna,NJ,MJ,NAMT)
       sum2 = 0.d0
       do 250 nc = 0, 2*jnnp
       nct = 2*nc + 1
       do 255 l = 1, nct
       ngama = - nc + (l-1)
       faux = exp((flog(nc-ngama+1) -flog(nc +ngama+1))*0.5)  
       faux = faux  * (-1)**ngama 
       alg2 = cleb6(jnnp+z, z, jnnp+z, z, nc+z, z)
       if (ngama.lt.0) then
       ngamapl = - ngama
       plaux = PL(nc+1,ngamapl+1)*(-1)**ngamapl
       plaux = plaux*exp((flog(nc-ngamapl+1)-flog(nc +ngamapl+1)))
       else
       plaux = PL(nc+1,ngama+1)
       endif
       sum2 = sum2 + alg2*aux(nc,ngama)*faux*plaux
 255   continue
 250   continue
       pleg(ith) = sum2
 200   continue

       do 400 ith = 1,nangles
       th = thmin + dth*(ith-1)
       thrad = th*pi/180.
       cthna = cos(thrad)
        call PLM(cthna,NJ,MJ,NAMT)
       do 450 nc = 0, 2*jnnp
       nct = 2*nc + 1
       do 455 l = 1, nct
       ngama = - nc + (l-1)
       if (ngama.lt.0) then
       ngamapl = - ngama
       plaux = PL(nc+1,ngamapl+1)*(-1)**ngamapl
       else
       plaux = PL(nc+1,ngama+1)
       endif
 455   continue
 450   continue
 400   continue
       deallocate(aux)
       return
       end



       function tctens(ib,ibb,i,j,ith)
         use scattering
         complex*16:: tctens
         if (ib<=nint(2*s3)) then
            tctens=tnath(ith)
         else
            write(*,*)'tctens: error in tensor arguments!. Aborting'
            write(*,*)'ib,ibb,i,j,ith=',ib,ibb,i,j,ith
            stop
         end if
       end function tctens

        


       subroutine tncore(elab,type,kin)
****************************************************************************
       use wfns
       use scattering
       use tatheta
       use constants
       implicit none
       integer::type,kin,nang,i,iqrr,istat,ith,ns
       real*8:: xscore,qrr,th,thrad,cthna,qq,x
!      real*8:: nu,nu14
       real*8,allocatable::qcore(:)
       real*8 :: e14,factna,factnc,mu,elab
       complex*16 ::wxxi
       complex*16 ::tncchk
       logical::ifmst=.true.
       zz = cmplx(0.d0,1.d0)
      
c **** calculates tNA (**SCATTERING AMPLITUDES**) from MST - U on-shell
c ***  ( calls MSOamp program )
c *** or reading S-matrix externally

       select case(type)
c--------------------------------------------------------------
         case(0) !T-matrix calculated by lptps subroutine
c--------------------------------------------------------------      
       nang=nangles
       write(*,1001)elab
1001   format(/," o Calling MSO with e14=",1f6.2," MeV")
       call lptps(elab,ifmst)
       nang=noangs
       if (allocated(qcore))deallocate(qcore)
       allocate(qcore(noangs))
       do i=1,noangs      
          x = cos(theta(i)*pi/180.)
          qcore(i) = kcore*sqrt(2 - 2*x)
       enddo
c **** interpolates to obtain the tmatrix on the **q** grid for MST - T
       allocate(ssctna(mdelta))
       do 200 iqrr = 1, mdelta
       qrr = radxisd(iqrr)
       ssctna(iqrr) = wxxi(qrr,qcore,tnct,noangs,noangs)
 200   continue
!------------------------------------------------------       

c--------------------------------------------------------------
       case (1) !analytical expression for scattering amplitude
! WARNING: This option has NOT been TESTED already!!!!
c--------------------------------------------------------------   
       nang=nangles
       if (allocated(tnct)) deallocate(tnct)
       allocate(tnct(nangles),stat=istat)
       if (allocated(qcore))deallocate(qcore)
       allocate(qcore(nangles))
       do i=1,nangles
          th = thmin + dth*dble(i-1)
          x = cos(th*pi/180d0)
          qcore(i) = kcore*sqrt(2 - 2*x)
          tnct(i)=zz*(kcore/4/pi)*sigtot*(1d0-zz*alpha)*
     &           exp(-qcore(i)**2*beta/2)
          tnct(i)=tnct(i)/sqrt(10d0) !sigtot is in mb
       enddo  

c------------------------------------------------------------
       case default !T-matrix from S-matrix read externally
c------------------------------------------------------------
       ns=type
       nang=nangles
       if (allocated(tnct)) deallocate(tnct)
       allocate(tnct(nangles),stat=istat)
       if (allocated(qcore))deallocate(qcore)

       allocate(qcore(nangles))
       do i=1,nangles
          th = thmin + dth*dble(i-1)
          if (kin.eq.2) call rihan(th)
          x = cos(th*pi/180d0)
          qcore(i) = kcore*sqrt(2 - 2*x)
       enddo
       if (istat>0) stop 'Could not allocate memory for tnath!'
       write(99,*)'- Reading external S-mat jp=',jp
c *** tnct will be calculated in qcore() grid by readsmat
       if (kin.eq.4) then
          mu=muna
       else
          mu=mucore
       endif
       write(*,'("- Calling readsmat with:",$)')
       write(*,*)'  kcore=',kcore,' mu (amu)=',mu/amu
       call readsmat(ns,tnct,kcore,z4,mu,jp)
       end select 

c AMoro 11/10/04: kin=5 uses the adiabatic projectile-core scattering 
c amplitude, as in kin=4, but this is calculated from the physical
c scattering amplitude making the high energy assumption
c f(q)~k g(q) and hence f_ad (q) = f_phys (q)*(k0/kcore)
       if (kin.eq.5) tnct=tnct*k0/kcore

c ** KL
       if (kin.eq.3) then
           tnct=tnct*(nu14/nu)*k0/kcore
           write(*,*)'Renormalizing Tcore by ',(nu14/nu)*k0/kcore
        endif

c **** interpolates to obtain the tmatrix on the **th** grid for MST - T
       if (allocated(tnath)) deallocate(tnath)
       allocate(tnath(nangles),stat=istat)
       if (istat>0) then
          write(*,*) 'tncore: allocating TNATH failed!.Aborting';
          stop
       endif
       do 300  ith = 1,nangles
       th = thmin + dth*dble(ith-1)
       if (kin.eq.2) call rihan(th)
       thrad = th*pi/180.
       cthna = cos(thrad)
      
       select case (kin)
       case(0,1,4,5) ! MST, Chew, Rihan, FSA kinematics
c      (the momentum transfer is the same  for 3B and 2B
          qq = sqrt(2*k0*k0*(1-cthna)) 

       case (2) !Rihan
          qq = sqrt(2*k0*k0*(1-cthna)) 
c  to correct for the angle-dependence of the energy parameter
c  within  Rihan approximation. Here we make use of  the energy 
c  relation f(q)~k g(q) at
          tnct(ith)=tnct(ith)*kcore/kc0

       case(3) !Kujawski & Lambert (KL)
          if (rel) then
             qq=(nu14/nu)*sqrt(2*k0*k0*(1-cthna))             
          else
             qq=(mucore/muna)*sqrt(2*k0*k0*(1-cthna))
          endif
          if (ith.eq.1) then
             write(*,*)'mucore/mu=',mucore/mu,'nucore/nu',nu14/nu
!             write(*,*)'kcore/k0=',kcore/k0,'nu14/nu',nu14/nu
!             write(*,*)'kcore=',kcore,'(nu14/nu)k0',(nu14/nu)*kcore
          endif

       case default
          write(*,*)'tncore: kin=',kin,'not used!';stop
       end select
       write(88,'(3f12.6)') qq,tnct(ith)
       tnath(ith)=wxxi(qq,qcore,tnct,nang,nang)
!       write(*,*)ith,qq,tnath(ith),tnct(ith)
 300   continue

c--------------------------------------------------------------
c ***  calculates cross section for the core as a check      
       do 400 ith = 1,nangles
       th = thmin + dth*dble(ith-1)
       if (kin.eq.2) call rihan(th)
       thrad = th*pi/180.
       cthna = cos(thrad)
!! Check kcore here for each kinematics!!!!!!???
       qq = sqrt(2*kcore*kcore*(1-cthna))
       if (type.eq.0) then
          tncchk = wxxi(qq,qcore,tnct,noangs,noangs)
       else
          tncchk = tnath(ith)
       endif
       xscore = tncchk * conjg(tncchk)*10
       write(75,'(f8.3,2e12.3)')th,xscore
 400   continue
101    format(7f12.6)
       deallocate(tnct)
       return
       end


c ***  Calculates s12,s14 for Rihan aproximation at angle theta
       subroutine rihan(theta)
         use scattering
         use constants
         implicit none
         real*8:: a12,c12,d12,e1,e4,e2,a14,c14,d14
         real*8:: theta,thrad,p0,p2,p4,costh
         real*8:: trian

         thrad = theta*pi/180.
         costh = cos(thrad)

         p0=hbarc*k0/amu ! amu
         e1=sqrt(m1**2 + p0**2) !amu

c ** Valence
         a12=m2/masst
         c12=(m3+m4)/masst

         p2=sqrt((0.5*c12+a12)**2 + 0.25*c12**2
     &      -c12*(0.5*c12+a12)*costh)*p0
         e2=sqrt(m2**2 + p2**2) !amu
         s12=(e1+e2)**2-(c12*p0)**2*(1d0+costh)/2d0 !amu**2
         nu12=e1*e2/(e1+e2) !amu

         kv=amu*sqrt(trian(s12,m1**2,m2**2))/
     &      2/sqrt(s12)/hbarc
         w12=(sqrt(s12)-m1-m2)*amu

c ** Core:
         a14=m4/masst
         c14=(m2+m3)/masst

         p4=sqrt((0.5*c14+a14)**2 + 0.25*c14**2
     &      -c14*(0.5*c14+a14)*costh)*p0

         e4=sqrt(m4**2 + p4**2) !amu
         s14=(e1+e4)**2-(c14*p0)**2*(1d0+costh)/2d0 !amu**2
         nu14=e1*e4/(e1+e4) !amu         
                 
         kcore=amu*sqrt(trian(s14,m1**2,m4**2))/
     &         2/sqrt(s14)/hbarc
         w14=(sqrt(s14)-m1-m4)*amu
       
         
         if (theta<0.4.or.theta>34.5) then
!         write(*,*)'Rihan: theta, a14,c14,d14=',theta,a14,c14,d14
            write(*,*)'theta:',theta
         write(*,*)'s12,w12:',s12,w12
         write(*,*)'s14,w14:',s14,w14
         write(*,*)'k1,kv,kcore=',k0,kv,kcore
         endif
         
        end subroutine rihan


          

!===============================================================
c ***  Calculates TNN tensor amplitudes by calling NNAMP program
!===============================================================
       subroutine tntensor(tcm,type,kin)
         use amps
         use parameters
         use wfns
         use scattering
         use nnamps
         use jj2
         use constants
         implicit real*8(a-h,o-z) 
         logical:: ifmst=.true.
         integer::ifkq=1,type,i,kin
         real*8::tcm,rr,tnnaux,k12
         real*8:: ms1, ms1p,sqmax,bqmax
         real*8,allocatable::th(:)
         real*8,pointer::xq(:)
         complex*16,pointer:: tnaux(:)
         complex*16,allocatable,target::mf(:,:,:,:,:)
         complex*16 ::aux, taux, auxp,wxxi
         dimension sred(0:1)
         complex*16,allocatable::tabkq(:,:,:,:,:),tnmat(:)
        
!         hat(i) = sqrt(real(2*I+1))
         namelist /amp/ ifkq,xkmax,xqmax,dk,dq,theta,nth,itype,icase
         write(99,*)'Entering tntensor with tcm=',tcm
         zz = cmplx(0.d0,1.d0)
         rr=0d0
         sred(0) = 1.
         sred(1) = sqrt(3.)
         itkqopt = 1.

cc AMORO: 16/12/03
!          s1 = 0.5
          s1=jp
        

c Type of S-matrix
c type=0: calculate NN T-matrix calling ampnn
c type>0: Reads S-matrix from unit type
         select case(type)
c        ------------------------------------------------
         case(0) ! calculate NN T-matrix calling ampnn   
c        -----------------------------------------------   
c  We need to know 'nth' beforehand for memory allocation
         close(10) 
         open(10,file='nnamp.in',status='old')         
         read(10,nml=amp) !nth read in this namelist
!         nthin=nth
         close(10)
        
         if (allocated(mf))deallocate(mf)
         allocate(mf(0:isnx,0:isnx,0:ikqx,-ikqx:ikqx,nth))  
         mf=0 !initialize
         write(99,*)'+ Allocating th with',nth,' angles'
         allocate(th(nth))

         if (tcm<1.e-5) then
            write(*,*) '**ERROR*** Energy not specified!. Aborting'
            stop
         endif
         

    
         sqmax=sqrt(2.*k0*k0*(1d0-cos(nth*pi/180)))
         bqmax=k0
         write(99,*)'tntensor: bqmax=k0=',bqmax
         
c AMoro (Oct/04)
c The NN code does not allow for relativistic kinematics
c so we send as argument a modified c.m. energy, such that the
c associated non-rel momentum coincides with the rel momentum
         if (kin.eq.4) then !FSA
            tnnaux=(hbarc*kv)**2/2/muna
         else
            tnnaux=(hbarc*kv)**2/2/muv
         endif
         write(*,fmt='(/," o Calling ampnn with tcm=",f8.4)')tcm
         write(99,*)' kv,tnnaux=',kv,tnnaux
         call ampnn(tnnaux,ifmst,sqmax,bqmax,ifkq)
            
         mfon(0,0,:,0,0,:)=mfon(0,0,:,0,0,:)*k0000
         mfon(1,1,:,0,0,:)=mfon(1,1,:,0,0,:)*k1100
         mfon(1,1,:,0,0,:)=mfon(1,1,:,0,0,:)*k1100
         mfon(0,1,:,1,1,:)=mfon(0,1,:,1,1,:)*k0111
         mfon(1,1,:,2,0,:)=mfon(1,1,:,2,0,:)*k1120
         mfon(1,1,:,2,1,:)=mfon(1,1,:,2,1,:)*k1121
         mfon(1,1,:,2,2,:)=mfon(1,1,:,2,2,:)*k1122
         
         mfon(1,0,0,1,1,:) = mfon(0,1,0,1,1,:)
         mfon(1,0,1,1,1,:) = mfon(0,1,1,1,1,:)


c        ----------------------------------------------------            
         case default !read external  nucleon-nucleus S-matrix
c        ----------------------------------------------------
            if (allocated(mf))deallocate(mf)
            nth=nangles
! Use nangles from mst.in
            allocate(mf(0:isnx,0:isnx,0:ikqx,-ikqx:ikqx,nangles))  
            mf=0 !initialize
            write(99,*)'+ Allocating th,xq with',nangles,' angles'
            allocate(th(nangles))

            if (allocated(tnmat)) deallocate(tnmat)
            allocate(tnmat(nangles),stat=istat)
            if (istat>0) then
               write(*,*)'tncore:could not allocate memory for tnmat!'
               write(*,*)'Aborting...'
               stop
            end if
            write(99,*)'Calling readsmat from tntensor:type=',type
            
            call readsmat(type,tnmat,kv,z2,muv,jp)

            write(99,*)'+ Allocating xxq with',nangles,'angles'
            if (allocated(xxq)) then
               deallocate(xxq)
            else
               allocate(xxq(nangles))
            endif

           
            do iqrr = 1, nangles
               thet = thmin + dth*dble(iqrr-1)
               thrad = thet*pi/180.
               cthna = cos(thrad)
               qrr = sqrt(2*kv*kv*(1-cthna))
c     xq is momentum transfer (fm -1)
                xxq(iqrr)=qrr
             enddo
             
             write(61,*)'tnmat, as read from readsmat'
             do i=1,nth
!               write(*,*)i,tnmat(i)
!               if (abs(tnmat(i))>1e3) 
!                write(61,*)i,10*abs(tnmat(i))**2
                write(61,*)thmin + dth*dble(i-1),10*abs(tnmat(i))**2
             enddo


            if (allocated(mfon)) then 
               deallocate(mfon)
            else
               write(99,*)'+ Allocating mfon with',nth,'angles'
               allocate(mfon(0:1,0:1,0:1,0:2,-2:2,nth),stat=istat)
               if(istat>0) then 
                  write(*,*) 'Could not allocate memory for mfon'
                  stop
               end if
            endif
            mfon=0
c *** make A(T=0)=A(T=1)
            mfon(0,0,0,0,0,1:nth)=tnmat  !T=0
            mfon(0,0,1,0,0,1:nth)=tnmat  !T=1
            deallocate(tnmat)
            
c        ------------------------------------------------------
         end select
c        ------------------------------------------------------         
     
   
c ***** --------------------------------------------------------             
         if (zp.eq.z2) then
            mf=mfon(:,:,1,:,:,:) ! pp / nn
            write(99,*)'- Tntensor for pp/nn'
         else 
            mf=(mfon(:,:,0,:,:,:)+mfon(:,:,1,:,:,:))/2. !pn
            write(99,*)'- Tntensor for pn'
         endif
         
         if (allocated(mfon)) then 
            write(99,*)'- Deallocating mfon'
            deallocate(mfon)
         else
            write(*,*)' ERROR (tntensor): mfon not allocated!'
            stop
         endif


         write(99,*)'+ Allocating tabkq',isnx,ikqx,ikqx,nangles  
         
         allocate(tabkq(0:isnx,0:isnx,0:ikqx,-ikqx:ikqx,nangles),
     &           stat=istat)
         if(istat>0) then 
            write(*,*) 'Could not allocate memory for tabkq';stop
         end if
         tabkq=0d0
            
         xq=>xxq 
c--------------------------------------------------
c **** interpolates to obtain the tmatrix on the required 
c **** theta grid for MST 
       do 191 iaj=0,2*s1 !isnx
       do 191 ibj=0,2*s2 !isnx
       do 191 ikj=abs(iaj-ibj),iaj+ibj 
       do 191 iqj=0,ikj
       if(iaj.eq.1.and.ibj.eq.1.and.ikj.eq.1)go to 191
       ikq = ikj + iqj
       if(ikq.eq.1)go to 191
       if (kapa(iaj,ibj,ikj,iqj).eq.0) cycle
       tnaux=>mf(iaj,ibj,ikj,iqj,1:nangles)
 
c AMoro: 12/10/04:
c With kin=5, the projectile-valence scattering amplitude
c is calculated from the physical amplitude by introducing appropriate
c mass/momentum factors (valid at high energies)
       if (kin.eq.5) tnaux=(k0/kv)*tnaux

c ** KL
       if (kin.eq.3) then
           tnaux=tnaux*(nu12/nu)*k0/kv
           write(*,*)'Renormalizing Tv by ',(nu12/nu)*k0/kv
        endif
     
       write(61,*)'tnmat:mf' 
       write(61,*)'iaj,ibj,ikj,iqj',iaj,ibj,ikj,iqj
       write(61,*)'nangles',nangles
       do 111 nnn=1,nth
          write(61,*)xq(nnn), 10*abs(tnaux(nnn))**2
111    continue    
          write(61,*)'interpolated tnmat: tabkq' 
          do 200 iqrr = 1, nangles
             thet = thmin + dth*dble(iqrr-1)
             thrad = thet*pi/180.
             if (kin.eq.2) call rihan(th)
             cthna = cos(thrad)

             select case (kin)
             case(0,1) ! MST,Chew,Rihan kinematics
                qrr = sqrt(2*k0*k0*(1-cthna))

             case(2)
                 qrr = sqrt(2*k0*k0*(1-cthna))
c  to correct for the angle-dependence of the energy parameter
c  within  Rihan approximation. Here we make use of  the energy
c  relation f(q)~k g(q) at
                tnaux(ith)=tnaux(ith)*kv/kv0

             case(3)!Kujawski & Lambert
c ** commented 22/12/04
!                qrr=sqrt(2*kv*kv*(1-cthna))
                if (rel) then
                   qqr=(nu12/nu)*sqrt(2*k0*k0*(1-cthna))
                else
                   qqr=(muv/muna)*sqrt(2*k0*k0*(1-cthna))
                endif

             case(4,5)
                qrr=sqrt(2*k0*k0*(1-cthna))

             case default
                write(*,*)'tntensor: kin=',kin,'not used!';stop
             end select


             tabkq(iaj,ibj,ikj,iqj,iqrr)=wxxi(qrr,xq,tnaux,nth,nth)
              write(61,*) qrr,10*abs(tabkq(iaj,ibj,ikj,iqj,iqrr))**2
200       continue
          nullify(tnaux)
191       continue


192    continue

       write(99,*)'+ Allocating memory for tnnu,tauvv,tauvc'
       write(99,*)'mllmx,isnx,nangles=',mllmx,isnx,nangles
       allocate(tnnu(0:isnx,-isnx:isnx,0:1,0:1,nangles))
       allocate(tauvv(nangles,0:mllmx,-mllmx:mllmx,0:isnx,0:isnx))  
       allocate(tauvc(nangles,0:isnx,-isnx:isnx))

       write(61,*)'tnnu-----------------------'

!! RCRESPO 6/01/04
!! change to allow spin zero of projectile
!! do 250 ibj = 0,isnx
       do 250 ibj = 0,2*s2
       do 250 ibjp = -ibj,ibj
!! AMORO: Change 17/12/03: 
!! Generalization from spin=0.5 to spin generic spin s1 
!!  do 280 i = 0,1    
!!  ms1 = - 0.5 + i
!!  do 280 j = 0,1
!!  ms1p = - 0.5 + j
       do 280 i = 0,2*s1  
        ms1 = - s1 + i     
        do 280 j = 0,2*s1  
        ms1p = - s1 + j 
         do 300 iqrr = 1, nangles
          aux = (0.d0,0.d0)
!! RCRESPO 6/01/04
!! change to allow spin zero of projectile
!! do 350 iaj=0,isnx
          do 350 iaj=0,2*s1
          do 350 iajp = -iaj,iaj
          do 350 ikj=abs(iaj-ibj),iaj+ibj 
          do 350 iqj=-ikj,ikj
          taux = tabkq(iaj,ibj,ikj,iqj,iqrr)

          if (iqj.lt.0) taux = 
     *    tabkq(iaj,ibj,ikj,-iqj,iqrr)*(-1)**(ikj-iqj)  
!! sred(iaj)*sred(ibj) spin reduced matrix elements for projectile (sred(iaj))
!! and struck nucleon (sred(ibj))
          aux = aux + (-1)**iqj  
     &    *sqrt(2*s2+1)/sqrt(2*ibj+1.) * sred(iaj)*sred(ibj)
     *    *cleb6(s1+rr,ms1+rr,iaj+rr,iajp+rr,s1+rr,ms1p+rr)
     *    *cleb6(iaj+rr,iajp+rr,ibj+rr,ibjp+rr,ikj+rr,-iqj+rr)*taux 

 350   continue       
       tnnu(ibj,ibjp,i,j,iqrr) = aux 
!!  in the limit of central interaction, the sum over the spin of the
!!  projectile of tnnu is normalized
!!  to (2*s1+1)*(2*s2+1)*T00
       write(61,332)iqrr,ibj,ibjp,i,j,10*abs(aux)**2!/(2*s2+1)
332    format(5i3,2g12.6)
!       if (abs(aux)>0) write(*,*)ibj,ibjp,i,j,iqrr,aux       
 300   continue
 280   continue
 250   continue

c***   Calculates the tensors to be used in the valence-valence and
c***   valence-core contributions
       tauvv=0d0
        write(61,*)'tauvv----------------------'
       do 390 iqrr = 1, nangles 
       do 400 if = 0,mllmx
       do 400 ifp = -if,if
!! RCRESPO 6/01/04
!! do 400 ib1j = 0,isnx
!! do 400 ib2j = 0,isnx
       do 400 ib1j = 0,2*s2
       do 400 ib2j = 0,2*s2
       do 400 ib1jp = -ib1j,ib1j
       do 400 ib2jp = -ib2j,ib2j
!! AMORO
!!$       do 400 i = 0,1
!!$       ms1 = - 0.5 + i
!!$       do 400 j = 0,1
!!$       ms1p = - 0.5 + j
!! The xs term valence-valence term is averaged over the spin of the projectile
!! that is it is devided by (2*s1+1) CHECK!
        do 400 i = 0,2*s1
       ms1 = - s1 + i
       do 400 j = 0,2*s1
       ms1p = - s1 + j
       auxp = tnnu(ib2j,ib2jp,i,j,iqrr)
       tauvv(iqrr,if,ifp,ib1j,ib2j) = tauvv(iqrr,if,ifp,ib1j,ib2j)
     &     + (-1)**(ifp+ib1jp) * tnnu(ib1j,ib1jp,i,j,iqrr) 
     &     * conjg(auxp)/(2*s1+1)
     &     * cleb6(ib1j+rr,ib1jp+rr,ib2j+rr,-ib2jp+rr,if+rr,ifp+rr)

       if (abs(tauvv(iqrr,if,ifp,ib1j,ib2j))>0) 
     & write(61,*)xq(iqrr),if,ifp,ib1j,ib2j,
     & 10*abs(tauvv(iqrr,if,ifp,ib1j,ib2j))
 400   continue
      
!! RCRESPO 6/01/04
!! do 450 ib1j = 0,isnx
       do 450 ib1j = 0,2*s2
       do 450 ib1jp = -ib1j,ib1j
       tauvc(iqrr,ib1j,ib1jp)  = (0.d0, 0.d0)
!!$       do 450 i = 0,1
!!$       ms1 = - 0.5 + i
!!$       do 450 j = 0,1
!!$       ms1p = - 0.5 + j
       do 450 i = 0,2*s1
       ms1 = - s1 + i
       do 450 j = 0,2*s1
       ms1p = - s1 + j
       tauvc(iqrr,ib1j,ib1jp) = tauvc(iqrr,ib1j,ib1jp) 
     *              +  tnnu(ib1j,ib1jp,i,j,iqrr)/(2*s1+1)
 450   continue
 390   continue

!!$       do 500 iqrr = 1, nangles 
!!$       qrr = radxisd(iqrr) 
!!$       write(93,777)qrr,tauvv(iqrr,0,0,1,1) 
!!$       write(94,777)qrr,tauvv(iqrr,1,0,1,1),tauvv(iqrr,1,-1,1,1)
!!$       write(95,777)qrr,tauvc(iqrr,0,0), tauvc(iqrr,1,0)
!!$ 500   continue
 777   format(f7.3,4e12.4)

 1023  format('  ',f5.1,f8.4,4e14.6)
         return
       end subroutine tntensor
         




       subroutine statescb(jnnp,lxp,lyp)
c************************************************************************
c      This subroutine verifies if a particular Jnn+ state will
c      contribute to the scattering. Presently up to eight states

       use scattering
       implicit real*8(a-h,o-z)
       integer parityf, parin

       parityf = (-1)**(lxp+lyp)
       jnncb = 1
       eps = 0.1
c      go to 60
       write(*,*)'states',jnnp, parityf
       do 50 j=1,nustates
       parin = (-1)**(1+partyin(j))      
       aux1 = abs(Jnnp-Jnnpin(j))
       aux2 = abs(parityf - parin)
       write(*,*)'states (j)',Jnnpin(j),parin
       if(aux1.lt.eps.and.aux2.lt.eps)then
       jnncb=istatcb(j)
       write(*,*) 'jnncb=',jnncb
       go to 60
       endif
 50    continue
c 60    write(*,*)'jnncb=',jnncb
 60    continue
       return
       end
       


      subroutine xsecinel(kin)
c****************************************************************************
c     This subroutine evaluates the single scattering inelastic cross section
c*****************************************************************************
       use parameters
       use wfns
       use scattering
       use trdens
       use nnamps
       use constants
!       implicit real*8(a-h,o-z)
       implicit none
       integer:: i,ic,if,im,ib,mif,mib,jnnp,ith,id,ibp,icp,in
!       parameter(isnx=1)
       integer :: llif,jnn,kin
       real*8 :: th,racah2,cleb6,factnn,factnc,factna,hat2
       real*8 :: ecm,ratio,zt,sigruth,aux1,rr
       real*8:: t, dsdt,sinv,cte,p0
       complex*16 wxxi
       complex*16 aux,auxp,aux2
       real*8:: nvv,nvc,ncc,a12,a14,d12,d14,c12,c14
       complex*16::doublex1,doublex2
       complex*16,allocatable::fstcore(:,:,:),fstval(:,:,:,:,:)
       complex*16,allocatable::fstcth(:,:,:),fstvth(:,:,:,:,:)!not used??
       complex*16,allocatable::fstxth(:,:,:,:,:)
       complex*16,allocatable::testval(:),testcore(:)
       complex*16,allocatable::testvc1(:),testvc2(:)
       real*8, allocatable::cgeomc(:,:)
       real*8, allocatable::cgeomv(:,:,:,:,:,:),cgeomx(:,:,:)
       real*8,allocatable::xscore(:),xsval(:),xsx(:),xstot(:),xsel(:)
       real*8, allocatable:: doublec(:,:),doublev(:,:)
       real*8, allocatable:: doublex(:,:) !,doublet(:,:)
       hat2(i)=(2*i+1)
       rr=0d0
       pi=acos(-1.)
       
       allocate(fstcore(nangles,0:mll,0:nwf))
       allocate(fstval(nangles,0:mll,0:mll,0:mll,0:nwf))
       allocate(fstxth(nangles,0:mll,0:isnx,0:mll,0:nwf))
       allocate(xscore(nangles))
       allocate(xsval(nangles))
       allocate(xsx(nangles))
       allocate(xstot(nangles))
       allocate(xsel(nangles))
       allocate(testval(nangles),testcore(nangles))
       allocate(testvc1(nangles),testvc2(nangles))
       allocate(doublec(nangles,0:nwf))
       allocate(doublev(nangles,0:nwf))
       allocate(doublex(nangles,0:nwf))
       allocate(doublet(nangles,0:nwf))
      
       
       allocate (cgeomc(0:mllmx,0:mll))
       allocate (cgeomv(0:mllmx,0:mll,0:isnx,0:dmax,0:mll,0:isnx))
       allocate (cgeomx(0:mll,0:isnx,0:dmax))

       factnn = -hbarc*hbarc*(m1+m2)/4/pi/pi/mn/(m1*m2)
       factnc = -hbarc*hbarc*(m4+m1)/4/pi/pi/mn/(m1*m4)
       factna = -hbarc*hbarc*(m234+m1)/4/pi/pi/mn/(m1*m234)
       jnn = TNQ(5,na(1),1)
      
       write(99,*)'xsecinel: factnn,factnc,factna',
     &           factnn,factnc,factna    


c ***  calculates the geometric coeficients: cgeomc(if,ic),
c ***  cgeomv(if,ic,ib,id,ic',ib'), cgeomx(ic,ib,id)           
       do 20 ic=0,mll
       do 20 if=0,mllmx
       cgeomc(if,ic) = 0.
           do 25 im = -ic,ic
	   cgeomc(if,ic) = cgeomc(if,ic) + 
     *         cleb6(ic+rr,rr,ic+rr,rr,if+rr,rr) 
     *         *(-1)**im * cleb6(ic+rr,-im+rr,ic+rr,im+rr,if+rr,rr)
 25        continue
       do 30 ib=0,isnx
       do 30 id=0,dmax
       cgeomx(ic,ib,id)=(-1)**id *cleb6(ic+rr,rr,id+rr,rr,ib+rr,rr)
     *                  *hat2(ic)/sqrt(hat2(ib)*hat2(id))
       do 30 ibp=0,isnx
       do 30 icp=0,mll
       cgeomv(if,ic,ib,id,icp,ibp) = (-1)**(ibp+ic-if) 
     *         * (-1)**( 2*(-ic+icp+id) ) 
     *         * hat2(ic)/hat2(id) * sqrt(hat2(ic)*hat2(icp))
     *         *cleb6(ic+rr,rr,icp+rr,rr,if+rr,rr)
     *         *racah2(ib+rr,id+rr,if+rr,icp+rr,ic+rr,ibp+rr)        
 30    continue
 20    continue        

      

c ***  calls the wigner rotation function
       call drotat() 
     
       write(61,*)'doublevv -----------------------------------'
c ***  Calculates the cross section contributions in th space
       doublev=0d0
       do 40 ith = 1,nangles
       th = thmin + dth*dble(ith-1)
      
c ***  Calculates mass factors according to kinematics
       if (rel) then
          nu=(s0+m1**2-masst**2)*(s0-m1**2+masst**2)/
     &        4/s0/sqrt(s0)
          nu12=(s12+m1**2-m2**2)*(s12-m1**2+m2**2)/
     &        4/s12/sqrt(s12)
          nu14=(s14+m1**2-m4**2)*(s14-m1**2+m4**2)/
     &        4/s14/sqrt(s14)     
        else
          nu=(m1*masst)/(m1+masst)
          nu12=m1*m2/(m1+m2)
          nu14=m1*m4/(m1+m4)
       endif

       
 
!!$        if (ith.eq.1) then
!!$            write(*,*)' TEST:'
!!$            write(*,*)'m1,m4,masst=',m1,m4,masst
!!$            write(*,*)'mu14=',m1*m4/(m1+m4),' nu14=',nu14
!!$            write(*,*)' ** nu/nu14',nu/nu14,'k0/kcore=',k0/kcore
!!$            write(*,*)'(nu/nu14)^2(kcore/k0)^2=',(nu*kcore/nu14/k0)**2
!!$         endif


c *** Relativistic
       select case(kin)
       case(0,1,3)! MST, Chew, KL
          nvv=(nu/nu12)**2
          ncc=(nu/nu14)**2
          nvc=sqrt(nvv*ncc)
          
!          if (ith.eq.1) then
!          write(*,*)'** nvv=',sqrt(nvv),muna/muv
!          write(*,*)'** ncc=',sqrt(ncc),muna/mucore
!          endif

       case(2) ! Rihan
          call rihan(th)
          a12=nu12/m2
          c12=(m3+m4)/masst
          d12=-a12*c12*(1d0-0.5*a12*c12) ! {\cal C} in paper
          
          a14=nu14/m4
          c14=(m2+m3)/masst
          d14=-a14*c14*(1d0-0.5*a14*c14)

c    Normalization constants for:
c    valence-valence 
           nvv=1d0/(1.+ d12 + d12*cos(th*pi/180.))
c    core-core
           ncc=1d0/(1.+ d14 + d14*cos(th*pi/180.))
c    valence-core
          nvc=sqrt(nvv*ncc)

       case(4,5) !FSA
          nvv=1d0
          ncc=1d0
          nvc=1d0
          
       case default 
          write(*,*)'Kinematics',kin,'not used!';stop
       end select
        
c ***  Calculates the cross section contribution: *core-core*,
       xscore(ith) = 0.d0 
       xsval(ith) = 0.d0 
       testval(ith) = (0.d0,0d0)
       xsx(ith) = 0.d0
       testvc1(ith) = (0.d0,0.d0) 
       testvc2(ith) = (0.d0,0.d0)
c ***  loop over excited states 
!! CHANGE 3/10/2003: amoro
       do 50 in=0,NWF-1 
       doublec(ith,in) = 0.d0
       doublev(ith,in) = 0.d0
       doublex(ith,in) = 0.d0
       jnnp = TNQ(5,na(in+1),in+1)         
c ***  loop over quantum numbers
       xscin(ith,in) = 0.d0

              
       do 60 ic=0,mll

!AMORO 20/2/04
       fstcore(ith,ic,in) = 
     *          tnath(ith)*ffcth(ith,ic,in)*sqrt(ncc)! factnc/factna

       do 60 if=0,mllmx
       llif = 0
       xscore(ith) = xscore(ith) +  dtheta(ith,if,llif)*cgeomc(if,ic)*
     *    fstcore(ith,ic,in)*conjg(fstcore(ith,ic,in))*10.     

       testcore(ith) = testcore(ith) + dtheta(ith,if,llif)*cgeomc(if,ic)
     *    *fstcore(ith,ic,in)*conjg(fstcore(ith,ic,in))*10.
       
       xscin(ith,in)=xscin(ith,in) + dtheta(ith,if,llif)*cgeomc(if,ic)*
     *    fstcore(ith,ic,in)*conjg(fstcore(ith,ic,in))*10.

       doublec(ith,in)=doublec(ith,in)+dtheta(ith,if,llif)
     *    * cgeomc(if,ic) *
     *    fstcore(ith,ic,in)*conjg(fstcore(ith,ic,in))*10.

!          if (abs(doublec(ith,in))>1.e-10)
!     &   write(60,'(1i4, 4g14.6,2i4)')  ith,doublec(ith,in)

 60    continue
c 50    continue

        
!       write(99,*)'mll,dmax,mllmx',mll,dmax,mllmx
       
c ***  Calculates the cross section contribution: *valence-valence*,
c       xsval(ith) = 0.d0 
c       testval(ith) = (0.d0,0d0)
c ***  loop over excited states 
c      do 150 in=1,NWF-1
       jnnp = TNQ(5,na(in+1),in+1)
c ***  loop over quantum numbers
       do 160 ic=0,mll
       do 160 ib=0,2*s2 !isnx
       do 160 id=0,dmax
       do 160 if=0,mllmx
       do 160 icp=0,mll
       do 160 ibp=0,2*s2 !isnx 
       do 160 mif =-if,if
        aux = conjg(ffvth(ith,icp,ibp,id,in))
     *        *cgeomv(if,ic,ib,id,icp,ibp)
     
        testval(ith) =  testval(ith) + aux*ffvth(ith,ic,ib,id,in)
     *      *dtheta(ith,if,mif) * tauvv(ith,if,mif,ib,ibp)
     *       *10.*2*nvv !(factnn/factna)**2

        xsval(ith) = real( testval(ith) )
        aux2=aux*ffvth(ith,ic,ib,id,in)
     *      *dtheta(ith,if,mif) * tauvv(ith,if,mif,ib,ibp)
     *       *10.*2*nvv !(factnn/factna)**2

        doublev(ith,in) = doublev(ith,in)+aux*ffvth(ith,ic,ib,id,in)
     *      *dtheta(ith,if,mif) * tauvv(ith,if,mif,ib,ibp)
     *       *10.*2*nvv !(factnn/factna)**2

!         if (abs(tauvv(ith,if,mif,ib,ibp))>10) 
!     &    write(61,*)ith,if,mif,ib,ibp,tauvv(ith,if,mif,ib,ibp)

!        if (abs(dtheta(ith,if,mif))>10) 
!     &   write(61,*)ith,if,mif,dtheta(ith,if,mif)

!      if (abs(doublev(ith,in))>1.e-10)
!     &   write(60,'(1i4, 4g14.6,2i4)')  ith,doublev(ith,in)
 160   continue
c 150   continue

c ***  Calculates the cross section contribution: valence-core,
c       xsx(ith) = 0.d0
c       testvc1(ith) = (0.d0,0.d0) 
c       testvc2(ith) = (0.d0,0.d0)
       doublex1 = (0.d0,0.d0)
       doublex2 = (0.d0,0.d0)
c ***  loop over excited states 
c       do 250 in=1,NWF-1
       jnnp = TNQ(5,na(in+1),in+1)
c ***  loop over quantum numbers
       do 260 ib=0,isnx
        aux = (0.d0,0.d0)      
        do 265 mib = -ib,ib   

        aux = aux + tauvc(ith,ib,mib)*
     *        dtheta(ith,ib,mib) 

 265   continue  
       

       do 260 ic=0,mll
       do 260 id=0,dmax
       fstxth(ith,ic,ib,id,in) = 
     &  conjg(tnath(ith)*ffcth(ith,id,in))*
     &  ffvth(ith,ic,ib,id,in)*cgeomx(ic,ib,id)
!!     &  *factnc/factna


        testvc1(ith) =  testvc1(ith) +
     &     fstxth(ith,ic,ib,id,in)*aux*10.*2
     &    *nvc   !factnn/factna

        testvc2(ith) =  testvc2(ith) + 
     *     conjg(fstxth(ith,ic,ib,id,in)*aux)*10.*2
     &   * nvc   !factnn/factna
!        xsx(ith) = real(  testvc1(ith) +  testvc2(ith)   )

       doublex1 = doublex1 +  
     &     fstxth(ith,ic,ib,id,in)*aux*10.*2
     &    *nvc !factnn/factna

       doublex2 = doublex2 +  
     *     conjg(fstxth(ith,ic,ib,id,in)*aux)*10.*2
     *    *nvc !factnn/factna
 260   continue
        doublex(ith,in) = real (doublex1 + doublex2)
!        doublex(ith,in) = doublex1*conjg(doublex1)
       
c 250   continue
        doublet(ith,in)= doublec(ith,in)+doublev(ith,in)
     *   +       doublex(ith,in)
 50     continue

!        xsx(ith) =  testvc1(ith) + conjg(testvc1(ith)) 
!       xstot(ith) = 0*xscore(ith)+ xsval(ith)+xsx(ith)
!        write(60,777)th,xstot(ith),xscore(ith)  
 40    continue

        
       write(99,*)'Writing elastic & total inelastic xsec...'
       do 300 ith = 1,nangles
       th = thmin + dth*dble(ith-1)
       xsel(ith)=0.0
       xstot(ith) = 0.d0
       xscore(ith) = 0.d0
c Elastic
        xsel(ith) = xsel(ith)   
     &            + doublev(ith,0) ! valence-valence
     &            + doublec(ith,0) ! core-core
     &            + doublex(ith,0) ! core-valence

c Inelastic
       do 350 in=1,NWF-1
       xstot(ith) = xstot(ith) + 
     *              doublec(ith,in)+doublev(ith,in)+doublex(ith,in)
       xscore(ith) = xscore(ith) +  doublec(ith,in)
 350   continue
       if (nwf>1) write(70,777)th,xstot(ith),xscore(ith) 
         
!       write(71,777)th,xsel(ith)
c *** If projectile and target are both charged, print ratio to Rutherford
       zt=z2+z3+z4
       if (zt*zp>0.) then
          ecm=tlab*m1/(m1+masst)
          aux1=sigruth(zp,zt,ecm,th)
       else
          aux1=1
       endif
         
c cte=> deprecated; no longer used
       cte=(s0*amu**2-(amu*masst+amu*m1)**2)*
     &     (s0*amu**2-(amu*masst-amu*m1)**2)/2d0/s0/amu/amu
       cte=cte/1e6 !convert to (GeV/c)^2
!       write(0,*) th,k0,sqrt(cte/2)*1e3/hbarc
       ratio=xsel(ith)/aux1

       p0=hbarc*k0/1e3 !linear momentum in GeV/c
       t=-(2*p0*sin(th*pi/180/2))**2 
       dsdt=xsel(ith)*(pi/p0/p0)
       write(72,777)th,xsel(ith),ratio,-t,dsdt
       ratio=doublev(ith,0)/aux1
       write(73,777)th,doublev(ith,0),ratio
       ratio=doublec(ith,0)/aux1
c Note that for the core we also use k0 ad not kcore
c because doublec is the core contribution to the scattering
c of the 3-BODY system (not 2-BODY system)
        dsdt=doublec(ith,0)*(pi/p0/p0)  
        write(74,777)th,doublec(ith,0),ratio,-t,dsdt

 300   continue    


 555   format(i5, 10e12.3)
 777   format(f8.3,4e12.3)
 999   format(f8.5,5e12.3)
       return
       end





      subroutine xsecinelnew(kin)
c****************************************************************************
c     This subroutine evaluates the single scattering inelastic cross section
c*****************************************************************************
       use parameters
       use wfns
       use scattering
       use trdens
       use nnamps
       use constants

       implicit none
       integer::  i,j,ic,if,im,ib,mif,mib,jnnp,ith,id,ibp,icp,in
       integer::  mnn,mnnp,iwf
       integer::  ibb, idd, icc
       integer::  llif,jnn,kin
       integer::  itest
       integer::  stout,stoutx,stin,stinx,is1nx
       real*8 ::  th,racah2,cleb6,factnn,factnc,factna,hat
       real*8 ::  clebaux,phase
       real*8 ::  ecm,ratio,zt,sigruth,aux1,rr
       real*8 ::  t, dsdt,sinv,cte,p0
       real*8 ::  nvv,nvc,ncc
       real*8 ::  a12,a14,d12,d14,c12,c14
       complex*16 tctens
       complex*16 tnaux1,tnaux2
       real*8 ::  dthrad, thrad, sinth, res, res2, xs2pth
       integer::  icount

       complex*16,allocatable::tmatt(:,:,:,:,:,:)
       complex*16,allocatable::tmatc(:,:,:,:,:,:)
       complex*16,allocatable::tmatv(:,:,:,:,:,:)
c      real*8,    allocatable::vgeom(:,:,:,:,:,:,:,:,:)
c      real*8,    allocatable::cgeom(:,:,:,:,:,:,:,:,:)
       real*8,    allocatable::xsel(:),xsinelt(:),xsinelc(:),xsinelv(:)
       real*8,    allocatable::xstotal(:,:),xscore(:,:),xsval(:,:)
       real*8,    allocatable:: xsecex(:), bs(:)
      
!-------------------------------------------------------------------------
!     Geom factors do not depend on in, and may be calcculated in the Tmat
!     and not need to be stored in an array
!--------------------------------------------------------------------------

!       hat2(i)=(2*i+1)
       hat(i)=sqrt(2.*i+1) 
       rr=0d0
       pi=acos(-1.)

c ***  Lim of quantum nbr defined in parameters ! DEFINE THIS IN PARAMETERS
c ***  isnx = 2*S2                              ! max of b !ASSUMED  isnx=1 
c ***  mllmx = 2*maximum of bound partial waves ! 
c ***  dmax = J0 + J                            ! max of d
c ***  mll=dmax + isnx                          ! max of c
c ***  These limits follow from the Racah coefficient
c ***  !  b  d    c
c ***  !  S' Jnn' L'
c **   !  J  Jnn  L
c ***  !  tnnu(b,beta,ms1,ms1',ith)) -> tnnu(0:isnx,-isnx:isnx,0:1,0:1,nangles)
       stout = 2                                ! Jnnp 
       stoutx = 2 * stout + 1                   ! 2 * Jnnp + 1
       stin = 0                                 ! Jnn
       stinx =  2 * stin + 1                    ! 2 * Jnn  + 1
       is1nx = 2*s1
       allocate (tmatc(
     &  0:is1nx,0:is1nx,-stout:stout,-stin:stin,nangles,0:nwf))
       allocate(tmatv(
     &  0:is1nx,0:is1nx,-stout:stout,-stin:stin,nangles,0:nwf))
       allocate(tmatt(
     &  0:is1nx,0:is1nx,-stout:stout,-stin:stin,nangles,0:nwf))
       allocate(xsel(nangles),xsinelt(nangles))           
       allocate(xstotal(nangles,0:nwf)) 
       allocate(xscore(nangles,0:nwf),xsval(nangles,0:nwf)) 
       allocate(xsinelc(nangles),xsinelv(nangles))
c       allocate(cgeom(0:isnx,-isnx:isnx,-stout:stout,-stin:stin,
c       0:dmax,-dmax:dmax,0:mllmx,-mllmx:mllmx,0:nwf))
c       allocate(vgeom(0:isnx,-isnx:isnx,-stout:stout,-stin:stin,
c        0:dmax,-dmax:dmax,0:mllmx,-mllmx:mllmx,0:nwf))
     
       factnn = -hbarc*hbarc*(m1+m2)/4/pi/pi/mn/(m1*m2)
       factnc = -hbarc*hbarc*(m4+m1)/4/pi/pi/mn/(m1*m4)
       factna = -hbarc*hbarc*(m234+m1)/4/pi/pi/mn/(m1*m234)
       jnn = TNQ(5,na(1),1)
       write(99,*)'xsecinel: factnn,factnc,factna',
     &           factnn,factnc,factna
       write(99,*)'max q nbrs: isnx,dmax, icmx',isnx,dmax,mll         
       
c ***  calls the wigner rotation function
       call drotat() 
     
       write(61,*)'Calculating Tmat  ---------------------------------'

c ***  Calculates the transition amplitudes and cross sections
       do 40 ith = 1,nangles
       th = thmin + dth*dble(ith-1)
      
c ***  Calculates mass factors according to kinematics
       if (rel) then
          nu=(s0+m1**2-masst**2)*(s0-m1**2+masst**2)/
     &        4/s0/sqrt(s0)
          nu12=(s12+m1**2-m2**2)*(s12-m1**2+m2**2)/
     &        4/s12/sqrt(s12)
          nu14=(s14+m1**2-m4**2)*(s14-m1**2+m4**2)/
     &        4/s14/sqrt(s14)     
        else
          nu=(m1*masst)/(m1+masst)
          nu12=m1*m2/(m1+m2)
          nu14=m1*m4/(m1+m4)
       endif

c *** Relativistic
c *** Calculates nvv ! ({\cal N}_12)**1/2 in paper
c *** Calculates ncc ! ({\cal N}_14)**1/2 in paper
       select case(kin)
       case(0,1,3)! MST, Chew, KL
          nvv=(nu/nu12)  
          ncc=(nu/nu14)  
          
!          if (ith.eq.1) then
!          write(*,*)'** nvv=',nvv,muna/muv   
!          write(*,*)'** ncc=',ncc,muna/mucore  
!          endif

       case(2) ! Rihan
          call rihan(th)
          a12=nu12/m2
          c12=(m3+m4)/masst
          d12=-a12*c12*(1d0-0.5*a12*c12) ! {\cal C} in paper
          
          a14=nu14/m4
          c14=(m2+m3)/masst
          d14=-a14*c14*(1d0-0.5*a14*c14)

c    Normalization constants for:
c    valence-valence 
           nvv=1d0/(1.+ d12 + d12*cos(th*pi/180.))
           nvv=sqrt(nvv) ! ANTO CHECK HERE
c    core-core
           ncc=1d0/(1.+ d14 + d14*cos(th*pi/180.))
           ncc=sqrt(ncc) ! ANTO CHECK HERE

       case(4,5) !FSA
          nvv=1d0
          ncc=1d0
          
       case default 
          write(*,*)'Kinematics',kin,'not used!';stop
       end select

c ***  loop over excited states 
       do 150 in=0,NWF-1
       iwf=in+1
       xstotal(ith,in)=0.d0
       xscore(ith,in)=0.d0
       xsval(ith,in)=0.d0
       jnnp = TNQ(5,na(in+1),in+1)

c ***  loop over spin of the projectile and spin of the target
       do 155 i = 0,2*s1        ! projectile
       do 155 j = 0,2*s1        ! projectile
       do 155 mnn=-jnn,+jnn     ! target
       do 155 mnnp=-jnnp,+jnnp  ! target

       tmatt(i,j,mnnp,mnn,ith,in)=(0.d0,0.d0)
       tmatc(i,j,mnnp,mnn,ith,in)=(0.d0,0.d0)
       tmatv(i,j,mnnp,mnn,ith,in)=(0.d0,0.d0)
c ***  loop over quantum numbers
       itest = 0
       do 160 ib=0,isnx 
       do 160 ibb=-ib,+ib
       do 160 id=0,dmax      
       do 160 idd=-id,+id
       do 160 ic=0,mll       
       do 160 icc=-ic,+ic
       itest = itest + 1
       phase = (-1)**(mnn+mnnp-idd)
       clebaux = cleb6(jnnp+rr,mnnp+rr,id+rr,idd+rr,jnn+rr,mnn+rr) *
     &           cleb6(ib+rr,ibb+rr,id+rr,idd+rr,ic+rr,icc+rr)
c       cgeom(ib,ibb,id,idd,ic,icc,mnnp,mnn,in) = hat(ic)*clebaux
c       vgeom(ib,ibb,id,idd,ic,icc,mnnp,mnn,in) = hat(ic)*phase*clebaux
        tmatc(i,j,mnnp,mnn,ith,in) = tmatc(i,j,mnnp,mnn,ith,in)  + 
     &           ncc * dtheta(ith,ic,icc)
     &         * hat(ic)*clebaux
     &         * tctens(ib,ibb,i,j,ith) * ffcth(ith,ic,in)
     &         /sqrt(2*s1+1)  ! Attention modify tctens tctens=0 if (i.neq.j)
        tmatv(i,j,mnnp,mnn,ith,in) = tmatv(i,j,mnnp,mnn,ith,in)  + 
     &         2 * nvv * dtheta(ith,ic,icc)   ! 2 val neutrons
     &        * hat(ic)*phase*clebaux
     &        * tnnu(ib,ibb,i,j,ith) * ffvth(ith,ic,ib,id,in) 
             

        tmatt(i,j,mnnp,mnn,ith,in) = 
     &     tmatc(i,j,mnnp,mnn,ith,in)+ tmatv(i,j,mnnp,mnn,ith,in)

   
 160   continue
 155   continue

!      ELASTIC AND INELASTIC CROSS SECTIONS

       write(61,*)'calculating Cross Sections -----------------------'
c ***  Sums the contributions over the spin of the projectile and target
c ***  Average over spin and multiply by 10 to convert fm to mb 
!      --->  ATTENTION CHECK AVERAGE OVER SPIN TARGET
       xscore(ith,in)  = 0.d0
       xsval(ith,in)   = 0.d0
       xstotal(ith,in) = 0.d0
       do 165 i = 0,2*s1        !projectile
       do 165 j = 0,2*s1        !projectile
       do 165 mnn=-jnn,+jnn     ! target
       do 165 mnnp=-jnnp,+jnnp  ! target
       xscore(ith,in) = xscore(ith,in) +
     &         tmatc(i,j,mnnp,mnn,ith,in)
     &         *conjg(tmatc(i,j,mnnp,mnn,ith,in))
     &         *10/(2.*s1+1)/hat(jnn)**2
!     &         *10/(2.*s1+1)/hat(jnnp)**2/hat(jnn)**2
       xsval(ith,in)  = xsval(ith,in) +
     &         tmatv(i,j,mnnp,mnn,ith,in)
     &         *conjg(tmatv(i,j,mnnp,mnn,ith,in))
     &         *10/(2.*s1+1)/hat(jnn)**2
!     &         *10/(2.*s1+1)/hat(jnnp)**2/hat(jnn)**2
       xstotal(ith,in)= xstotal(ith,in) +
     &         tmatt(i,j,mnnp,mnn,ith,in)
     &        *conjg(tmatt(i,j,mnnp,mnn,ith,in))
     &         *10/(2.*s1+1)/hat(jnn)**2
!     &         *10/(2.*s1+1)/hat(jnnp)**2/hat(jnn)**2

       if(ith.eq.11)then
       write(999,*)i,j,mnn,mnnp
       write(999,*)'ith=',ith,'in=',in
       write(999,*)'xscore(ith,in)',xscore(ith,in)
       write(999,*)'xsval(ith,in)',xsval(ith,in)
       endif
 165   continue

 150   continue !loop en in
 40    continue

       write(99,*)'Writing elastic &  total inelastic xsec...'
       do 300 ith = 1,nangles
       th = thmin + dth*dble(ith-1)
       xsel(ith)    = 0.d0
       xsinelt(ith) = 0.d0
       xsinelc(ith) = 0.d0
       xsinelv(ith) = 0.d0
c ***  Elastic
        xsel(ith) =   xstotal(ith,0)    
c ***  Inelastic
       do 350 in=1,NWF-1
       xsinelt(ith) = xsinelt(ith) +  xstotal(ith,in) 
       xsinelc(ith) = xsinelc(ith) +  xscore(ith,in) 
       xsinelv(ith) = xsinelv(ith) +  xsval(ith,in) 
 350   continue
c ***  Write inel cross section in theta unit 71
       if (nwf>1) write(70,777)th,xsinelt(ith),xsinelc(ith),xsinelv(ith)       
!       write(71,777)th,xsel(ith)
c ***  If projectile and target are both charged, print ratio to Rutherford
       zt=z2+z3+z4   
       if (zt*zp>0.) then
!          write(99,*)'Writing ratios to Rutherford'
          ecm=tlab*m1/(m1+masst)
          aux1=sigruth(zp,zt,ecm,th)
       else
          aux1=1
       endif
         
       ratio=xsel(ith)/aux1
       p0=hbarc*k0/1e3                 !linear momentum in GeV/c
       t=-(2*p0*sin(th*pi/180/2))**2 
       dsdt=xsel(ith)*(pi/p0/p0)
       write(72,777)th,xsel(ith),ratio,-t,dsdt
       ratio= xsval(ith,0)/aux1
       dsdt= xsval(ith,0)*(pi/p0/p0) 
       write(73,777)th, xsval(ith,0),ratio,-t,dsdt
       ratio= xscore(ith,0)/aux1
c ***  Note that for the core we also use k0 ad not kcore
c ***  because xscore is the core contribution to the scattering
c ***  of the 3-BODY system (not 2-BODY system)
        dsdt= xscore(ith,0)*(pi/p0/p0)  
        write(74,777)th, xscore(ith,0),ratio,-t,dsdt
        dsdt= tnath(ith)*conjg( tnath(ith))*ncc*ncc*10*(pi/p0/p0)
        write(744,777)th,-t,dsdt
 300   continue    

c ***  Write double cross section in unit 65, 66 
        write(65,*)thmin,thmax,dth
        write(65,*)NWF-1
!       write(98,fmt='("#Double: 4i4)')thmin,thmax,dth,NWF-1
        do 450 ith = 1,nangles
        do 450 in=NWF-1,1,-1
        write(65,*)ens(in),xstotal(ith,in)
        th = thmin + dth*dble(ith-1)
        write(66,*)ens(in),th,xstotal(ith,in)
 450    continue


      
       write(99,*)' - In sigex:'
       write(99,fmt='("+ Allocating",i4," angles for bs")')nangles
       allocate(bs(nangles))
       allocate(xsecex(0:nwf))

 
        dthrad = dth*pi/180.
        icount = 0.
        do 550 in=NWF-1,1,-1
        icount = icount + 1
        sum = 0.d0
        do 560 ith = 1,nangles
        th = thmin + dth*dble(ith-1)
        thrad = th*pi/180.
        sinth = sin(thrad)
        bs(ith) =  xstotal(ith,in)*sinth
        sum = sum + bs(ith)
 560    continue
        call sim(bs,res,1,nangles,dthrad,mthmx)
        xsecex(in) = res * 2 * pi
        call sim(bs,res2,46,90,dthrad,mthmx)
        xs2pth = res2 * 2 * pi
        write(96,*) ens(in), xsecex(in),xs2pth
 550     continue
        deallocate(xsecex)



 555   format(i5, 10e12.3)
 777   format(f8.3,4e12.3)
 999   format(f8.5,5e12.3)
       return
       end




       subroutine xsection(inelcb)
c************************************************************************
c      old subroutine adds final states incoherently
       use wfns
       use scattering
       use bfst
       use legendre
      implicit real*8(a-h,o-z)
       complex*16 tmat1,snd
       complex*16 cx1q(150), cx2q(150)
       complex*16 wxxi
       nangles = (thmax-thmin)/dth + 1
!       call  ptheta(inelcb,pleg) 
        call  ptheta(inelcb) 
       do 100 iqrr = 1, mdelta
       qrr = radxisd(iqrr)
       qdelta(iqrr) = qrr
       tcore(iqrr) = (0.d0,0.d0)
       tvalence(iqrr) = (0.d0,0.d0)
       snd = tcore(iqrr) + tvalence(iqrr)
       if (ival.eq.1)  snd = tcore(iqrr)
       xsect1 = fst(1,iqrr) * conjg(fst(1,iqrr)) * 10.
       cx1q(iqrr) = xsect1
       write(9,*)'old cx1q', cx1q(iqrr),fst(1,iqrr)
       xsect2 = (fst(1,iqrr)+snd) * conjg(fst(1,iqrr)+snd) * 10.
       cx2q(iqrr) = xsect2
 100   continue

       do 150 ith = 1,nangles
       th = thmin + dth*dble(ith-1)
       thrad = th*pi/180.
       cthna = cos(thrad)
       qq = sqrt(2*k0*k0*(1-cthna))
       xsect1 = wxxi(qq,qdelta,cx1q,200,mdelta)
       xsect2 = wxxi(qq,qdelta,cx2q,200,mdelta)
       xsect1 = pleg(ith)*xsect1
       write(9,*)th,xsect1
 150   continue
       return
       end



       subroutine sigex()
******************************************************************************
*     This subroutine calculates the integrated differential cross section for
*      all the excited states from the core contribution.
*******************************************************************************
       use parameters
       use wfns
       use scattering
       use trdens
       use constants
       implicit real*8(a-h,o-z)
       real*8,allocatable:: bs(:) !(mthmx)
       real*8,allocatable:: xsecex(:)

       
       write(99,*)' - In sigex:'
       write(99,fmt='("+ Allocating",i4," angles for bs")')nangles
       allocate(bs(nangles))
       allocate(xsecex(0:nwf))

 
        dthrad = dth*pi/180.
        icount = 0.
        do 20 in=NWF-1,1,-1
        icount = icount + 1
        sum = 0.d0
        do 40 ith = 1,nangles
        th = thmin + dth*dble(ith-1)
        thrad = th*pi/180.
        sinth = sin(thrad)
c       bs(ith) =  xscin(ith,in)*sinth
        bs(ith) =  doublet(ith,in)*sinth
        sum = sum + bs(ith)
 40     continue
        call sim(bs,res,1,nangles,dthrad,mthmx)
        xsecex(in) = res * 2 * pi
        call sim(bs,res2,46,90,dthrad,mthmx)
        xs2pth = res2 * 2 * pi
        write(96,*) ens(in), xsecex(in),xs2pth
c        write(96,*) ens(icount), xsecex(in),xs2pth
 20     continue

c **    Write double cross section 
        write(65,*)thmin,thmax,dth
        write(65,*)NWF-1
!        write(98,fmt='("#Double: 4i4)')thmin,thmax,dth,NWF-1
        do 60 ith = 1,nangles
        do 60 in=NWF-1,1,-1
c       write(65,*)xscin(ith,in)
        write(65,*)ens(in),doublet(ith,in)
        th = thmin + dth*dble(ith-1)
        write(66,*)ens(in),th,doublet(ith,in)
 60     continue

!       write(98,*)'#differential cross section'
        do 75 ith = 1,nangles
        th = thmin + dth*dble(ith-1)
        auxth = 0.d0
        do 70 in=NWF-1,1,-1
        auxth = auxth + doublet(ith,in)
 70     continue
!        write(98,*)th,auxth
 75     continue
        deallocate(xsecex)
        return

1000    write(*,*)'Could not open input file mst.in!';stop
        end





      subroutine shakeoff()
*********************************************************************
c     Calculates the inelstic cross section using the shakeoff
c     mechanism

       use wfns
       use scattering
       implicit real*8(a-h,o-z)
       complex*16 tnct(181)
       complex*16 wxxi
       complex*16 tncchk
       dimension theta(181),qcore(181)

       zz = cmplx(0.d0,1.d0)

       rewind 80
c **** reads tNA (**SCATTERING AMPLITUDES**) from MST - U on-shell
       read(80,*) noangs,kcore
       do 100 i=1,noangs      
       read(80,*)theta(i),tncre,tncimag
       tnct(i) = tncre+zz*tncimag
       x = cos(theta(i)*pi/180.)
       qcore(i) = kcore*sqrt(2 - 2*x)
 100   continue
       

c **** interpolates to obtain the tmatrix on the **th** grid for MST - T
       do 300  ith = 1,nangles
       th = thmin + dth*dble(ith-1)
       thrad = th*pi/180.
       cthna = cos(thrad)
       qq = sqrt(2*k0*k0*(1-cthna))
       tnath(ith) = wxxi(qq,qcore,tnct,181,noangs)
 300   continue


c ***  calculates cross section 

       do 400 ith = 1,nangles
       th = thmin + dth*dble(ith-1)
       thrad = th*pi/180.
       cthna = cos(thrad)
       qq = sqrt(2*kcore*kcore*(1-cthna))
       tncchk = wxxi(qq,qcore,tnct,181,noangs)
       aux = 1 - ffcth(ith,0,0)**2
       xsoff = tncchk * conjg(tncchk)*10*aux
       write(73,*)th,xsoff, xsoff*0.5
 400   continue
 

       return
       end




      subroutine test()
c***********************************************************************

      use parameters
      use wfns
      use scattering
      use trdens
      implicit real*8 (a-h,o-z)
      complex*16 wxxi


!***    write interpolated densities in q and th grid
        write(40,*)'core interpolated in q grid'
 	do 140 iq=1,mdelta
	q = radxisd(iq)
  	write(40,155) q,ffcore(iq,1,3), ffcore(iq,0,0)
 140    continue
        write(40,*)'valence  interpolated in q grid'
 	do 145 iq=1,mdelta
	q = radxisd(iq)       
  	write(40,106) q,(    (dreal(ffval(iq,ic,1,id,2)),
     X    id=0,2 ),  ic=0,2 )
 145  continue
      
        do 150  ith = 1,nangles
        th = thmin + dth*dble(ith-1)
	write(23,106) th,(    (dreal(ffvth(ith,ic,1,id,2)),
     X    id=0,2 ),  ic=0,2 )
        write (24,155)th,ffcth(ith,1,3), ffcth(ith,0,0)
        aux1 = ffcth(ith,1,3)**2
        aux2 = 1 - ffcth(ith,0,0)**2
        write (25,155)th, aux1,aux2
150     continue 

106	format(f8.3,9f8.3)
155     format(f8.3,2f10.6)


      return
      end



