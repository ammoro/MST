c----------------------------------------------------------------------
c *** Calculates optical potential in t*rho approximation
c *** for each cluster
c     R. Crespo and A.M. Moro (2002)
c----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)      
      call ffactbs()    
      end

**********************************************************************
      subroutine ffactbs()
      implicit real*8 (a-h,o-z)
*--------------------------------------------------------------------- 
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
             subroutine sim(fa,res,m,n,dr,nramax)  
             real*8  ::res,dr,fa(:)
             integer ::m,n
             end subroutine sim 
      end interface
*----------------------------------------------------------------------- 
      integer :: ia,ib,ic,is2,j2a,lmoma,nodd,nramax,jhw,kcheck
      real*8  :: cmass,vmass,wzz,wal,wr0,wls,vdepth,bengy,dmat,
     1           wrz,drx,pnloc,wr0ls,wals
      integer ::norba
!     real*8 ::chis(10),ffr(2000,4),bs(2000),dd(2000)
!     real*8 ::densr(10,2000),u(20,2000)
      real*8 ::chis(10)
      real*8,allocatable ::ffr(:,:)
      real*8,allocatable ::bs(:),dd(:),densr(:,:),u(:,:)
      parameter (ic=1,jhw=1,pnloc=0.d0)
*----------------------------------------------------------------------
*     read single particle bound state formfactors
*----------------------------------------------------------------------
      open(unit=5,file='tape5.data',status='unknown')
!      open(unit=6,file='output.data',status='unknown')
      open(19,file='wfn.out',status='unknown')
      bengy = 0.d0
      vdepth = 0.d0
      namelist /bsparm/ cmass,vmass,zc,zv,
     &                  nramax,drx,dmat,nshell

      namelist /bshell/ ia,bengy,vdepth,
     &                  ib,wr0,wal,wls,
     &                  j2a,lmoma,is2,nodd,norba

      read(5,nml=bsparm)
      write(6,nml=bsparm)
      allocate(densr(nshell,nramax))
      allocate(u(nshell,nramax))
      allocate(ffr(nramax,4))
      allocate(bs(nramax))
      allocate(dd(nramax))

!     loop from outer to inner shells
      do 888 ibs =1,nshell
      read(5,nml=bshell)
*----------------------------------------------------------------------
*          ' --------------------------------------'
*          ' Bound states in a two body potential  '
*          ' ------------------------------------  '
*     search well depth for ia=0, search energy for ia=1
*     separation energy: bengy
*     potential depth: vdepth
*     core and bound particle masses: cmass,vmass
*     core and bound particle charges: zc,zv
*     potential type ib=
*                       0 = Woods-Saxon
*                       1 = Gauss
*                       2 = Yukawa
*                       3 = Hulthen
*                       4 = Cosh
*     potential radius and diffuseness: wr0,wal
*     spin-orbit potential strength (~6 MeV): wls
*     2*j value, lmom, and nodes (from 0): j2a,lmoma,nodd
*     2*valence particle spin: is2
*     no. integ steps (max 2000) and step: nramax,drx
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
      write(6,21) vdepth,bengy,wal,wr0,wzz,wls,ia,ib,ic,
     1 lmoma,nodd,j2a,is2,dmat
   21 format('  vdepth= ',f10.5,'  bengy= ',f10.5,'  wa= ',f10.5,
     + //,'  wr0= ',f10.5,'  wzz= ',f10.5,'  wls= ',f10.5,//,
     +   '  ia= ',i5,'  ib= ',i5,'  ic= ',i5,'  lmom= ',i5,
     +   '  nod= ',i5,//,'  2*j=',i5,' 2*s= ',i5,' dmat= ',f10.5)
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
!     densr(ibs,i) = r*r*u(ibs,i)**2 * norba
      do 86 i=1,nramax
      r=(i-1)*drx
      bs(i)=r*r*u(ibs,i)**2
      densr(ibs,i)=bs(i)*norba
      dd(i) = densr(ibs,i)
   86 continue
      call sim(bs,rnorm,1,nramax,drx,nramax)
      call sim(dd,rnd,1,nramax,drx,nramax)
      write(6,*),'normalization =',rnorm,rnd 
      do 90 i=1,nramax
      r=(i-1)*drx
      bs(i)=bs(i)*r*r
   90 continue
      call sim(bs,rnorm,1,nramax,drx,nramax)
      write(6,'(" Wave function relative rms",1f8.3)')sqrt(rnorm) 

 888  continue
*---------------------------------------------------------------------
*     print bound wave functions 
*---------------------------------------------------------------------
      write(19,*)'wave functions'
      do 507 i=1,nramax
      rad = (i-1)*drx
      write(19,508)rad,u(1,i),u(2,i)
 507  continue
 508  format(e12.5,3x,e12.5,3x,e12.5)
*----------------------------------------------------------------------
*     calculate the transformed density
*----------------------------------------------------------------------
!     Grid and loop in q for printing
      qmin=0.05
      qmax=2.
      dq=0.05
      m=(qmax-qmin)/dq+1
      write(6,*)'rhoq for protons/neutrons'
      do 900 iq=1,m
      q = qmin + (iq-1)*dq
      formf=0.d0
      do 605 ibs=1,nshell
      do 602 i=1,nramax
      r=(i-1)*drx
      qr=q*r
      dd(i)=densr(ibs,i)
      if(qr.gt.1.d-10) dd(i)=dd(i)*sin(qr)/qr
  602 continue
      call sim(dd,densq,1,nramax,drx,nramax)      
      formf = formf + densq
!     end loop in shell
 605  continue
!     end loop in q 
      write(6,*)q,formf
 900  continue

!     Grid and loop in q for MSO and storage

      deallocate(densr,u,ffr,bs,dd)

      return
      end
*************************************************************************
      subroutine sim(fa,res,m,n,h,nramax)
*------------------------------------------------------------------------
*     subroutine does the integral of fa stored
*     in the array of the same name using simpsons rule. the step length
*     is h and the integral is between the elements m and n of the arrays
*     only. resulting integral is placed in res.
*------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
!     dimension fa(2000),dq(2000)
      real*8 ::res,fa(nramax),dq(nramax),h
      integer m,n,nramax
      do 90 i=m,n
      dq(i)=fa(i)
   90 continue
      rq1=dq(m+1)
      rq2=dq(m+2)
      i=m+3
   98 continue
      if(i.ge.n) go to 99
      rq1=rq1+dq(i)
      rq2=rq2+dq(i+1)
      i=i+2
      go to 98
   99 continue
      res=0.33333333333d0*h*(dq(m)+4.d0*rq1+2.d0*rq2-dq(n))
      return
      end
***********************************************************************
      subroutine bound(ia,ib,ic,is2,j2,lmom,nod,nramax,jhw,cmass,vmass,
     1wzz,wa,wr0,wls,ffr,vdepth,bengy,dmat,kcheck,wrz,chis,drx,pnloc,
     1wr0ls,wals)
      implicit real*8(a-h,o-z)
!     dimension ffr(2000,4),chis(10),drd(4)
!     dimension wfcr(2000),wfsr(2000),wfc(2000)
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
!      do 8 i=1,800
      do 8 i=1,nramax
      yyy=(wr-radi)/wa
      zzz=(wr-radls)/wals
      if(yyy.gt.90.d0) yyy=90.d0
      if(zzz.gt.90.d0) zzz=90.d0
      ex=exp(yyy)
      exls=exp(zzz)
      ib1=ib+1
      go to (6100,6200,6300,6350,6355),ib1
 6100 wfcr(i)=1.0/(1.0+ex)
      go to 6400
 6200 wfcr(i)=exp(-(wr/wa)*(wr/wa))
      go to 6400
 6300 wfcr(i)=exp(-wa*wr)/wr
      go to 6400
 6350 wfcr(i)=exp(-wa*wr)/(exp(-wr0*wr)-exp(-wa*wr))
      go to 6400
 6355 wfcr(i) = (1.0+yyyy)/(cosh(wr/wa)+yyyy)
 6400 wfsr(i)=exls/(1.0+exls)/(1.0+exls)
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

