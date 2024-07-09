      Module parameters
      	parameter (mqx=200,mromx=5000,mpotmx=500,mqbmx=1000)
      end module parameters

      module scattering
	use parameters
        implicit real*8(a-h,o-z)
        real*8 pi,hbarc,mn,mp,mav,k0,k0p
        real*8 mpj, mtg
        integer  nztg, natg, nzpj, napj
        integer iqmx
        real*8 tlab,mu,mup
        real*8 rmspj, rmstg
        integer idenpj,identg 
        real*8 Fact,alpha,gama
        real*8 qmin,qmax,step
        real*8 rrmin,rrmax,drr
        integer nq,nrr
        complex*16 z
        real*8 qbin(mqbmx)
        real*8 rotg(mqbmx), ropj(mqbmx)
        real*8 pot(4,mpotmx)
      end module scattering

      module nnamp
        use parameters
        implicit real*8(a-h,o-z)
        integer  nkmx1,nqmx1
        real*8 q(mqx),qq(mqx)
        complex*16 ap(mqx,mqx),cp(mqx,mqx)
        complex*16 an(mqx,mqx),cn(mqx,mqx)
        complex*16 tnn(mqbmx,3),tnnls(mqbmx,3)
      end module nnamp

!     program dfold.f

      use nnamp
      use scattering
      implicit real*8(a-h,o-z)
!      open(unit=9,file='outdfold.dat',status='unknown')
!      open(unit=20,file='pot.dat',status='unknown')
!      open(unit=25,file='potP.dat',status='unknown')
!      open(unit=10,file='getin.dat',status='unknown')
!      open(unit=50,file='denstg.dat',status='unknown')
!      open(unit=60,file='denspj.dat',status='unknown')
       open(unit=30,file='redish.data',form='unformatted')


       mp = 938.2796
       mn = 939.5731
       mav = (mp + mn)/2.
       pi = 4d0*atan(1d0)
       z = cmplx(0.d0,1.d0)
       hbarc = 197.389

       call getins()
       mpj = napj * mav
       mtg = natg * mav
       mu = mpj*mtg/(mpj+mtg)
       mup = mp*mtg/(mp+mtg)

       k0 = 2 * mu*mu * tlab /mpj
       k0 = sqrt(k0)/hbarc
       k0p = k0 *mup/mu
       write(9,*)'k0',k0,'k0p',k0p
       Fact = 2*natg/(natg+1.)/k0p/k0p

       alpha = (nztg-2.)/3.
c **   Fourier transform coefficient + plane wave coefficinet
       gama = 4.*pi/(2*pi)**3
c **   Scattering to transition factor
       gama = - gama * hbarc*hbarc/(2.*pi*pi)/mav
c **   Densities are calculated normalized to unity
       gama = gama * natg *napj 

       iqmx = (qmax)/step
       if (iqmx.gt.mqbmx)then
       write(9,*)'iqmx.gt.mqbmx - check input grids'
       stop
       endif
       do 20 iq=1,iqmx
       qbin(iq) = step + (iq-1)*step
 20    continue

c ***  reads NN transition amplitudes
       call readnn()
       call gtnn()
c ***  reads densities for projectile and target in configuration space
c      and evaluates the Fourier transform in momentum space
       write(*,*)'denspt....'
       call denspt()
c ***  calculates double folding in configuration space
       write(*,*)'dbl...'
       call dbl()
 
       close(unit=50)
       close(unit=60)    
       end


      function cint2d(xtab,ytab,fxytab,xbar,ybar,nnx,nny,nord,
     1 mmx)
c***************************************************************
c
c          2 dimensional complex aitkin interpolation routine.  Note that
c          mesh points must be in increasing order.
c
      implicit double precision(a-h,o-z),integer*4(i-n)
      complex*16 cint2d
      dimension xytab(2,200),x(10),y(10),
     1 xybar(2),nbg(2),xtab(mmx),ytab(mmx),nxy(2)
      complex*16 fxytab(mmx,mmx),fxy(10,10)
      nnn=nord+1
c          set up arrays for loop
      do 40 j=1,nnx
 40   xytab(1,j)=xtab(j)
      do 45 j=1,nny
 45   xytab(2,j)=ytab(j)
      nxy(1)=nnx
      nxy(2)=nny
      xybar(1)=xbar
c          begin loop to determine the (order+1) x and y points over whi
c          interpolation is made; 1=x, 2=y.
      xybar(2)=ybar
      do 10 i=1,2
30    num=1
      if(xybar(i).lt.xytab(i,1))go to 85
      num=nxy(i)-nnn+1
      if(xybar(i).gt.xytab(i,nxy(i)))go to 85
50     min=1
      max=nxy(i)
      num=nxy(i)/2
55    if(max-min.lt.4)goto70
      if(xytab(i,num)-xybar(i))60,82,61
60     min=num
          num=num+(max-min+1)/2
      goto55
  61   max=num
      num=(max-min+1)/2
      goto55
  70   num=max
  71  if(xytab(i,num)-xybar(i))82,82,81
 81   num=num-1
      goto71
 82   num=max0(1,num-nnn/2)
      num=min0(num,nxy(i)-nnn+1)
 85   nbg(i)=num
 10   continue
c          end loop.  set up x, y, function arrays of required
c          (order+1)**2 points.
      nx=nbg(1)
      ny=nbg(2)
      do 20 i=1,nnn
      x(i)=xytab(1,nx+i-1)
      y(i)=xytab(2,ny+i-1)
      do 20 j=1,nnn
      fxy(i,j)=fxytab(nx+i-1,ny+j-1)
 20   continue
c          do interpolation
 90   do95ii=2,nnn
      do95 jj=ii,nnn
      do 95 mm=ii,nnn
      fxy(jj,mm)=((fxy(ii-1,ii-1)*(x(jj)-xybar(1))-fxy(jj,ii-1)*(x(ii-1)
     1 -xybar(1)))*(y(mm)-xybar(2))-(fxy(ii-1,mm)*(x(jj)-xybar(1))
     2 -fxy(jj,mm)*(x(ii-1)-xybar(1)))*(y(ii-1)-xybar(2)))
     3 /(x(jj)-x(ii-1))/(y(mm)-y(ii-1))
  95   continue
      cint2d=fxy(nnn,nnn)
      return
      end
c

      subroutine dbl()
c***********************************************************************
       use parameters
       use scattering  
       use nnamp
       implicit real*8 ( a-h,o-z )
       dimension taux(4)
       dimension bs(mqbmx)

       nrr = (rrmax)/drr
       do 50 ir=1,nrr
       r = drr + (ir-1)*drr 
       do 70 nifty=1,4
       do 80 iq=1,iqmx
       qbb = qbin(iq)
       taux(1) = real(tnn(iq,3))
       taux(2) = imag(tnn(iq,3))
       taux(3) = real(tnnls(iq,3))
       taux(4) = imag(tnnls(iq,3))

       if (nifty.eq.1.or.nifty.eq.2)then
       bs(iq) =  qbb*sin(qbb*r)*taux(nifty)*rotg(iq)*ropj(iq)/r
       else if (nifty.eq.3) then
       bs(iq) = -qbb*(qbb*cos(qbb*r)-sin(qbb*r)/r)*taux(4)*Fact*
     *    rotg(iq)*ropj(iq)/r/r
       else
       bs(iq) = qbb*(qbb*cos(qbb*r)-sin(qbb*r)/r)*taux(3)*Fact*
     *    rotg(iq)*ropj(iq)/r/r
       end if
 80    continue
       call sim(bs,res,1,iqmx,step)
       pot(nifty,ir) = res * gama 
       if (napj.gt.1) then
       pot(nifty,ir) = pot(nifty,ir)*sqrt( (2.*pi)**3)
       endif
 70    continue
 50    continue

       write(20,*)nrr,rrmin,drr
       do 100 ir=1,nrr  
       r = rrmin + (ir-1)*drr 
       write(20,*)pot(1,ir)+ z*pot(2,ir)
       write(25,*) r, pot(1,ir),pot(2,ir)
 100   continue

 777   format(f8.5,4e12.3) 

       return
       end
 
      subroutine getins ()                                              
c*********************************************************************
      use parameters
      use scattering 
      use nnamp
      implicit real*8  ( a-h , o-z ) 
      read(10,*)nztg,natg, nzpj,napj
      read(10,*)tlab
      read(10,*) qmin,qmax,step
      read(10,*)rrmin,rrmax,drr
      read(10,*)identg, idenpj, rmstg, rmspj

      write(9,*)'input parameters'
      write(9,*)nztg,natg, nzpj,napj
      write(9,*)tlab
      write(9,*) qmin,qmax,step
      write(9,*)rrmin,rrmax,drr
      write(9,*)identg, idenpj, rmstg, rmspj
 
      return
      end 
   

       subroutine gtnn()
c**************************************************************************
       use parameters
       use nnamp
       use scattering
       implicit real*8 (a-h, o-z)
       dimension tapb(4),te(4)
       complex*16 taux
       complex*16 cint2d
   
c*** qp and qqp are transfer and total n-n momentum
c *** qqp evaluated on the energy shell

      do 100 iq=1,iqmx 
      qp = qbin(iq) 
      if(qp.lt.0)qp=-qp
      qqp = k0**2  - qp/4.
      if(qqp.lt.0.)qqp=-qqp
      qp = sqrt(qp)
      qqp = sqrt(qqp)/2.
  
      taux=cint2d(qq,q,ap,qqp,qp,nkmx1,nqmx1,5,mqx)
      tapb(1)=real(taux)
      tapb(3)=dimag(taux)
      taux=cint2d(qq,q,an,qqp,qp,nkmx1,nqmx1,5,mqx)
      tapb(2)=real(taux)
      tapb(4)=dimag(taux)
      taux = cint2d(qq,q,cp,qqp,qp,nkmx1,nqmx1,5,mqx)
      te(1) = real(taux)
      te(3) = dimag(taux)
      taux = cint2d(qq,q,cn,qqp,qp,nkmx1,nqmx1,5,mqx)
      te(2)=real(taux)
      te(4)=dimag(taux)

c ***  average
       tnn(iq,3) = tapb(1) + z*tapb(3)
       tnnls(iq,3) = te(1) + z*te(3)
c ***  p-n
       tnn(iq,2) = tapb(2) + z*tapb(4)
       tnnls(iq,2) = te(2) + z*te(4)
c ***  p-p
       tnn(iq,1) = tnn(iq,3)*2 - tnn(iq,2)
       tnnls(iq,1) = tnnls(iq,3)*2 - tnnls(iq,2)
 100   continue

      return
      end


       subroutine denspt()
c-----------------------------------------------------------------------
c      Evaluated the density at the required grid normalized to unity
c-----------------------------------------------------------------------
       use parameters
       use scattering
       implicit real*8(a-h,o-z)

       dimension bs(mromx),densr(mromx)

c *******************************************************************
c      Target normalized to unity
c      identg=0 gaussian,  identg=1 read file 
c*******************************************************************
       if (identg.eq.0)then
       do 40 iq=1,iqmx
       btg = rmstg*sqrt(2/3.)
       btg2 = btg*btg
       rotg(iq) = exp(-btg2*qbin(iq)**2/4.)
c      rotg(iq) = rotg(iq)/(natg*1.)
 40    continue
       else if (identg.eq.1)then       
c***   reads target density in configuration space
       read(50,*)rmin,rmax,drx
       nramax = (rmax-rmin)/drx + 1
       if (nramax.gt.mromx)then
       write(9,*)'nramax.gt.mromx'
       stop
       endif
       do 10 i=1,nramax
       read(50,*)r, densr(i)   
 10    continue
c***  calculates the transformed density rho(q)
      do 53 iq=1,iqmx 
      do 52 i=1,nramax
      r=rmin+(i-1)*drx
      qr=qbin(iq)*r
      if(qr.lt.1.d-10) then
      bs(i)=densr(i)*r*r
      else
      bs(i)=densr(i)*r*r*sin(qr)/qr
      endif
  52  continue
      call sim(bs,res,1,nramax,drx)
      rotg(iq) = res*4*pi
      rotg(iq) = rotg(iq)/(natg*1.)
      if (rotg(iq).lt.0.00001)rotg(iq)=0.d0
      if (iq.eq.1)write(9,*)'rho(q) of the target'
      write(9,*)qbin(iq),rotg(iq)
  53  continue
      endif 

      do 100 iq=1,iqmx
      ropj(iq) = 1
 100  continue

      if (napj.gt.1)then
c********************************************************************
c     Projectile normalized to number of nucleons
c      idenpj=0 gaussian,  idenpj=1 read file 
c*******************************************************************
       if (idenpj.eq.0)then
       do 140 iq=1,iqmx
       btg = rmspj*sqrt(2/3.)
       btg2 = btg*btg
       ropj(iq) = napj*exp(-btg2*qbin(iq)**2/4.)
c       ropj(iq) = ropj(iq)/(napj*1.)
       write(9,*)-qbin(iq),ropj(iq)
 140   continue
       else if (idenpj.eq.1)then    
c***   reads projectile density in configuration space
       read(60,*)rminp,rmaxp,drxp
       nramaxp = (rmaxp-rminp)/drxp + 1
       if (nramaxp.gt.mromx)then
       write(9,*)'nramaxp.gt.mromx'
       stop
       endif
       do 150 i=1,nramaxp
       read(60,*)r, densr(i)   
 150   continue
c***  calculates the transformed density rho(q)
      do 203 iq=1,iqmx 
      do 202 i=1,nramaxp
      r=rmin+(i-1)*drx
      qr=-qbin(iq)*r
      if(qr.lt.1.d-10) then
      bs(i)=densr(i)*r*r
      else
      bs(i)=densr(i)*r*r*sin(qr)/qr
      endif
  202 continue
      call sim(bs,res,1,nramax,drx)
      ropj(iq) = res*4*pi
      ropj(iq) = ropj(iq)/(napj*1.)
      if (rotg(iq).lt.0.00001)rotg(iq)=0.d0
      write(9,*)-qbin(iq),ropj(iq)
  203 continue
      endif     
      endif     

       return
       end



      subroutine sim(fa,res,m,n,h)
***************************************************************************
*     subroutine does the integral of fa stored
*     in the arrays of the same name using simpsons rule. the step length
*     is h and the integral is between the elements m and n of the arrays
*     only. resulting integral is placed in res.
*------------------------------------------------------------------------
      use parameters
      implicit real*8(a-h,o-z)
      dimension fa(mromx),dq(mromx)
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


      subroutine readnn ()
************************************************************************
      use parameters
      use nnamp
      implicit real*8  ( a-h , o-z )

      read(30)nkmx1,nqmx1
      if (nkmx1.gt.mqx.or.nqmx1.gt.mqx)then
      write(9,*)'nkmx1.gt.mqx.or.gt.nqmx1.mqx'
      stop
      endif     
      read(30)q,qq,ap,cp

      return
      end









