************************************************************************
*     Subroutines of general use, shared by several source files
************************************************************************

************************************************************************
*     REAL 4-point lagrange interpolation routine.
*     interpolates thr FUNCTION value fival at point r from an
*     array of points stored in fdis(ndm). this array is assumed
*     to be defined such that the first element fdis(1) CONTAINS
*     the FUNCTION value at r=xv(1) and xv(2 .. ndm) are monotonically
*     increasing.
************************************************************************
      FUNCTION fival(r,xv,fdis,ndm,alpha)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 fdis(ndm),y1,y2,y3,y4
      DIMENSION xv(ndm)
      IF(r.GT.xv(ndm)) go to 9
      DO 5 k=1,ndm-2
 5    IF(r.LT.xv(k)) go to 6
      k=ndm-2
 6    nst=MAX(k-1,1)
      x1=xv(nst)
      x2=xv(nst+1)
      x3=xv(nst+2)
      x4=xv(nst+3)
      y1=fdis(nst+0)
      y2=fdis(nst+1)
      y3=fdis(nst+2)
      y4=fdis(nst+3)
      pii1=(x1-x2)*(x1-x3)*(x1-x4)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)
      xd1=r-x1
      xd2=r-x2
      xd3=r-x3
      xd4=r-x4
      pi1=xd2*xd3*xd4
      pi2=xd1*xd3*xd4
      pi3=xd1*xd2*xd4
      pi4=xd1*xd2*xd3
      fival=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4
      RETURN
 9    fival=fdis(ndm) * EXP(alpha*(xv(ndm)-r))
      RETURN
      END





*******************************************************************
* Subroutines shared by different subprograms
* f90 version A.Moro, R. Crespo (2002)
*******************************************************************


!****************************************************************
      function cint2d(xtab,ytab,fxytab,xbar,ybar,nnx,nny,
     1 nord,mmx)
c
c          2 dimensional complex aitkin interpolation routine.  note that
c          mesh points must be in increasing order.
c
      implicit none
      integer :: nnx,nny,nord,mmx,nxy(2),nbg(2),nnn,j,i,num,min,max,
     1   nx,ny,ii,jj,mm
      real*8 :: xbar,ybar,xybar(2)
      real*8 :: xytab(2,100),x(10),y(10)
      real*8 :: xtab(*),ytab(*)
      complex*16, dimension(nnx,nny) :: fxytab
      complex*16 :: fxy(10,10),cint2d
      nnn=nord+1
   
c      write(*,*) 'Entering cint2d',xtab(1),xtab(nnx),ytab(1),ytab(nny)
c      write(*,5000) xtab,ytab
5000  format(10e12.3)

c       write(3,*) ' fxytab=',fxytab
c      print*,'cint2db:', xbar,ybar,nnx,nny,nord,mmx
c        write(4,*) ' xtab=',xtab(1),xtab(kmax)

c          set up arrays for loop
      do 40 j=1,nnx
c          write(3,*) 'j=',j,'xtab(j)',xtab(j)
 40   xytab(1,j)=xtab(j) 
     
      do 45 j=1,nny
c          write(3,*) 'j=',j,'ytab(j)',ytab(j)
 45   xytab(2,j)=ytab(j)

c          write(3,*) '--------------------'
      nxy(1)=nnx
      nxy(2)=nny
        
      xybar(1)=xbar
c          begin loop to determine the (order+1) x and y points over which
c          interpolation is made; 1=x, 2=y.
      xybar(2)=ybar

c      write(3,*) 'xtab=',xtab,'ytab=',ytab
      do 10 i=1,2
30    num=1
      if(xybar(i).lt.xytab(i,1))go to 85
      num=nxy(i)-nnn+1
c xytab(i,nxy(i))
      if(xybar(i).gt.xytab(i,nxy(i)))go to 85
50     min=1
      max=nxy(i)
      num=nxy(i)/2

c       write(3,*) 'xbar=',xbar,'ybar=',ybar
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
c      write(3,*) 'i, j, fxy=',i,j,fxy(i,j)
 20   continue
c          do interpolation
 90   do95 ii=2,nnn
      do95 jj=ii,nnn
      do 95 mm=ii,nnn
      fxy(jj,mm)=((fxy(ii-1,ii-1)*(x(jj)-xybar(1))-fxy(jj,ii-1)*(x(ii-1)
     1 -xybar(1)))*(y(mm)-xybar(2))-(fxy(ii-1,mm)*(x(jj)-xybar(1))
     2 -fxy(jj,mm)*(x(ii-1)-xybar(1)))*(y(ii-1)-xybar(2)))
     3 /(x(jj)-x(ii-1))/(y(mm)-y(ii-1))
c      write(*,*) 'ii,jj,mm,fxy=',ii,jj,mm,fxy(jj,mm)
  95   continue
      cint2d=fxy(nnn,nnn)
c      write(3,*) 'cint2d=',cint2db
c      print*,'Saliendo de cint2db'
      return
      end




c  @(#)pldblp.f  2.2

      subroutine pldblp(xx,y,n)
************************************************************************

c     calculates the second derivative of the legendre pol,all l
c      taking derivative of eqs 3.28 3.29 jackson elimhna
c      tingdpl/dx

      implicit real*8 (a-h, o-z)

!      dimension y(100)

      real*8::y(1:n-1)

c >>> first executable statement <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      y(1)=0.
      y(2)=0.
      y(3)=3.
      if(n.le.3) return
        nn=n-1
      do 2 i=3,nn
   2    y(i+1)=((2*i-1)*xx*y(i)-(i+1)*y(i-1))/(i-2)
      return
c**  this program valid on ftn4 and ftn5 **
      end



c   @(#)plprme.f  2.2

      subroutine plprme (x, plp, n)
c***********************************************************************

c *** plprme calculates the derivative of the legrendre polynomials.
c *** it now uses the formulae of jackson directly by dimensioning the
c *** array plp from 0 to 49.

c *** input arguments
c ***    x   = cos(theta)
c ***    plp = array containing the pl primes evaluated at x
c ***    n   = number of pl primes the calling program wants
c ***
      implicit real*8 (a-h, o-z)
      real*8 :: plp(0:n)
!      dimension plp(0:99)


c >>> first executable statement <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      plp(0) = 0.0
      plp(1) = 1.0
      if (n .gt. 2) then
!         nn = n - 2
         do 10 l = 1,n-1
            plp(l+1) = ( (2*l+1)*x*plp(l) - (l+1)*plp(l-1) )/l
  10     continue
      endif
      return
      end






c   @(#)legpol.f 2.2

      subroutine legpol(x, plofx, n)
************************************************************************
     
c *** uses recurrence formula given in jackson`s e&m (pg. 90)
c *** directly by dimensioning plofx from l = zero to l = 49
c *** and thus differs significantly from the old code
c ***
c *** x = cos(theta)
c *** plofx`s are the legendre polynomials
c *** n = # of pl`s wanted


      implicit real*8 (a-h, o-z)
      integer :: l,n,nn
      real*8 :: plofx(0:n)
!       dimension   plofx(0:99)
!      if (allocated(plofx)) deallocate(plofx)
!      allocate(plofx(0:n))


c>>>> first executable statement <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!      if (n.gt.50) then
!        write(*,*)'legpol: degree must be < 50, but you have',n 
!        stop
!      endif

      plofx(0) = 1d0
      plofx(1) = x
      if (n .gt. 2) then
         do 100 l = 1,n-1
            plofx(l+1) = dble(2 * l + 1) * x *  plofx(l) - 
     &                    dble(l) * plofx(l-1)
            plofx(l+1) = plofx(l+1)/dble(l+1)
 100     continue
      endif
!      write(*,*) 'Dentro legpol: x,plofx',x,plofx(1),plofx(n-1)
      return
      end




      subroutine cmatin(a, b, n)
************************************************************************
c
c
c *** finds the inverse of complex matrix h=(rea,aima)=(a,b)
c
c *** code provided by norman bardsley to kwon+tabakin s prog bopit
c *** n=no of elements in a assigned(and in the sub=dimension)
c *** matrix a,b(real+imag parts) destroyed by call
c *** the other arguments are dummy just to save space

      implicit real*8 (a-h, o-z)
      integer, allocatable :: indexa(:),indexb(:),ipivot(:)
      real*8 :: a(n,n),b(n,n)

      equivalence (irow,jrow), (icolum,jcolum)

      allocate(indexa(n),indexb(n),ipivot(n))

c ***
c *** detr=1.0
c *** rhl add to define tempi for real potential

      tempi = 0.

c *** deti=0.0

      do 10 j = 1,n
  10     ipivot(j) = 0

      do 70 i = 1,n
         xmax = 0.0
         do 35 j = 1,n
            if (ipivot(j) - 1) 15, 35, 15

  15        do 30 k=1,n
               if (ipivot(k) - 1) 20, 30, 90
  20           if (abs(a(j,k)) .lt. 1.0e-20) go to 30
               xjk = abs(a(j,k)) + abs(b(j,k))
               if (xmax - xjk)  25,30,30
  25           irow = j
               icolum = k
               xmax = xjk
  30        continue

  35  continue

      ipivot(icolum) = ipivot(icolum) + 1
      if (irow .eq. icolum) go to 50

c *** detr=-deti
c *** deti=-deti

      do 45 l = 1,n
         swap = a(irow,l)
         swapi = b(irow,l)
         a(irow,l) = a(icolum,l)
         b(irow,l) = b(icolum,l)
         b(icolum,l) = swapi
  45     a(icolum,l) = swap

  50  indexa(i) = irow
      indexb(i) = icolum
      pivotr = a(icolum,icolum)
      pivoti = b(icolum,icolum)

c *** temp=detr*pivotr-deti*pivoti
c *** deti=detr*pivoti+deti*pivotr
c *** detr=temp

      a(icolum,icolum) = 1.0
      b(icolum,icolum) = 0.0
      if (pivoti .eq. 0.0) tempr = 1.0/pivotr
      if (pivoti .ne. 0.0) then
         tempr =  pivotr/(pivotr * pivotr + pivoti * pivoti)
         tempi = -pivoti/(pivotr * pivotr + pivoti * pivoti)
      endif

      do 55 l = 1,n
         temp = a(icolum,l) * tempr - b(icolum,l) * tempi
         b(icolum,l) = a(icolum,l) * tempi + b(icolum,l) * tempr
  55     a(icolum,l) = temp

      do 70 l1 = 1,n
         if (l1 - icolum) 60, 70, 60
  60     tempa = a(l1,icolum)
         tempb = b(l1,icolum)
         a(l1,icolum) = 0.0
         b(l1,icolum) = 0.0

         do 65 l = 1,n
            b(l1,l) = b(l1,l) - a(icolum,l) * tempb
     $                - b(icolum,l) * tempa
            a(l1,l) = a(l1,l) - a(icolum,l) * tempa
     $                + b(icolum,l) * tempb
  65     continue
  70  continue

      do 85 i = 1,n
         l = n + 1 - i
         if (indexa(l) - indexb(l)) 75, 85, 75
  75     jrow = indexa(l)
         jcolum = indexb(l)

         do 80 k = 1,n
            swap = a(k,jrow)
            swapi = b(k,jrow)
            a(k,jrow) = a(k,jcolum)
            b(k,jrow) = b(k,jcolum)
            a(k,jcolum) = swap
            b(k,jcolum) = swapi
  80     continue
  85  continue

            
  90  deallocate(indexa,indexb,ipivot)
      return
      end


c
c

c @(#)lagrng.f 2.1 8/24/86 18:11:44

      subroutine lagrng(x,arg,y,val,ndim,nfs,nptx,maxarg,maxfs)
************************************************************************

c     lagrange interpolation,unequally spaced points
c     npts=2,3,4,5,6, nfs functions(y*s) simultaneously interpltd
c     x= value of argument, arg is the tabulated x*s
c     y= a vector of interpolad functions, from tabulted val
c     ndim= dimension of table
c     nfs= = of functions simult interpolated
c     maxarg and maxfs are maximum values of subscripts,ie dimensions
      implicit real*8 (a-h, o-z)
      dimension arg(maxarg),val(maxfs,maxarg),y(maxfs)
c-----find x0 the closest point to x


c >>> first executable statement <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      npts=nptx
      ni=1
      nf=ndim
   10 if((x .le. arg(ni)) .or.  (x .ge. arg(nf)))go to 30
      if((nf-ni+1) .eq. 2) go to 40
      nmid=(nf+ni)/2
      if( x .gt. arg(nmid)) go to 20
      nf=nmid
      go to 10
   20 ni=nmid
      go to 10
c------ x is one of the tabltd values
   30 if( x .le. arg(ni)) go to 35
      nn=nf
   31 nused=0
      do 32 n=1,nfs
   32 y(n)=val(n,nn)
      return
   35 nn=ni
      go to 31
c------- 2 pts left, chpose smaller one
   40 n0=ni
   50 nn=npts-1
      go to (200,300,400,500,600),nn
  600 continue
      if(((n0+3) .gt. ndim) .or. ((n0-2) .lt. 1)) go to 500
      nused=6
      go to 1000
  500 continue
      if((n0+2) .gt. ndim) go to 300
      if((n0-2) .lt. 1) go to 400
      nused=5
      go to 1000
  400 continue
      if(((n0+2).gt. ndim) .or. ((n0-1) .lt. 1)) go to 300
      nused=4
      go to 1000
  300 if((n0+1) .lt. ndim) go to 305
c    n0=ndim,special case
  302 nn=ndim
      go to 31
  305 nused=3
      if((n0-1) .lt. 1) nused=2
      go to 1000
  200 if((n0+1) .gt. ndim) go to 302
      nused=2
 1000 continue
c     at least 2 pts left
      y0=x-arg(n0)
      y1=x-arg(n0+1)
      y01=y1-y0
      c0=y1/y01
      c1=-y0/y01
      if(nused .eq. 2) go to 2000
c     at least 3 pts
      ym1=x-arg(n0-1)
      y0m1=ym1-y0
      ym11=y1-ym1
      cm1=-y0*y1/y0m1/ym11
      c0=c0*ym1/y0m1
      c1=-c1*ym1/ym11
      if(nused .eq. 3) go to 3000
c------at least 4 pts
      y2=x-arg(n0+2)
      ym12=y2-ym1
      y02=y2-y0
      y12=y2-y1
      cm1=cm1*y2/ym12
      c0=c0*y2/y02
      c1=c1*y2/y12
      c2=-ym1*y0*y1/ym12/y02/y12
      if(nused .eq. 4) go to 4000
c     at least 5 pts
      ym2=x-arg(n0-2)
      ym2m1=ym1-ym2
      ym20=y0-ym2
      ym21=y1-ym2
      ym22=y2-ym2
      cm2=ym1*y0*y1*y2/ym2m1/ym20/ym21/ym22
      cm1=-cm1*ym2/ym2m1
      c0=-c0*ym2/ym20
      c1=-c1*ym2/ym21
      c2=-c2*ym2/ym22
      if( nused .eq. 5) go to 5000
c     at least 6 pts
      y3=x-arg(n0+3)
      ym23=y3-ym2
      ym13=y3-ym1
      y03=y3-y0
      y13=y3-y1
      y23=y3-y2
      cm2=cm2*y3/ym23
      cm1=cm1*y3/ym13
      c0=c0*y3/y03
      c1=c1*y3/y13
      c2=c2*y3/y23
      c3=ym2*ym1*y0*y1*y2/ym23/ym13/y03/y13/y23
      go to 6000
 2000 continue
      do 2100 n=1,nfs
 2100  y(n)=c0*val(n,n0)+c1*val(n,n0+1)
      go to 7000
 3000 continue
      do 3100 n=1,nfs
 3100 y(n)= cm1*val(n,n0-1)+c0*val(n,n0)+c1*val(n,n0+1)
      go to 7000
 4000 continue
      do 4100 n=1,nfs
 4100 y(n)=cm1*val(n,n0-1)+c0*val(n,n0)+c1*val(n,n0+1)+c2*val(n,n0+2)
      go to 7000
 5000 continue
      do 5100 n=1,nfs
 5100 y(n)=cm2*val(n,n0-2)+cm1*val(n,n0-1)+c0*val(n,n0)+c1*val(n,n0+1)
     1+c2*val(n,n0+2)
      go to 7000
 6000 continue
      do 6100 n=1,nfs
 6100 y(n)=cm2*val(n,n0-2)+cm1*val(n,n0-1)+c0*val(n,n0)+c1*val(n,n0+1)
     1 +c2*val(n,n0+2)+c3*val(n,n0+3)
 7000 return
      end







       subroutine gauss3(rmaxr,quin,mquadi,mquado,xri,wri)
c--------------------------------------------------------------------
c     calculates quadrature points in steps of equal size
c     from 0 to quin, and quadrature points in steps of linerly
c     size from quin to rmaxr
c     j.a. tostevin
c-------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      real*8:: xg(6),wg(6),xri(201),wri(201)
      integer ::kout=16
*---------------------------------------------------------------------
*     presently works with 6 point quadrature. for more change
*     nww, input the quadrature weights and pivots
*---------------------------------------------------------------------
      data nww,xg(4),xg(5),xg(6),wg(4),wg(5),wg(6)/3,
     1   .2386191861d0,.6612093865d0,.9324695142d0,
     2   .4679139346d0,.3607615730d0,.1713244924d0/
c
      ramp(x)=x*(sign(0.5d0,x)+0.5d0)
      if(mquadi+mquado.gt.201)then
      write(kout,*)'two many points in gauss3'
      stop
      endif
      nw=2*nww
      do 1 n=1,nww
      nn=nw-n+1
      xg(n)=-xg(nn)
    1 wg(n)=wg(nn)
*---------------------------------------------------------------------
c      print*
c      print*,'input rmaxr  : the maximum radius'
c      print*,'      quin   : the maximum radius of the inner region'
c      print*,'      mquadi : the number of inner quad points'
c      print*,                (must be a multiple of 6)
c      print*,'      mquado :  the number of outer quad points'
c      print*,                (must be a multiple of 6)
*---------------------------------------------------------------------
c     read *,rmaxr,quin,mquadi,mquado
      quon=rmaxr-quin
      mquad =mquadi+mquado
*---------------------------------------------------------------------
      write(kout,820)quin,mquadi,quon,mquado
  820 format(' quadrature :: inner space of ',f6.2,' fm.  divided into'
     # ,i4,' steps',
     #/'            :: outer space of ',f6.2,' fm.  divided into'
     # ,i4,' steps of linearly increasing size.'/)
*-------------------------------------------------------------------
      ic3=mquadi/nw
      ci=quin/ic3
      r1=0.5*ci
      r2=r1+0.d0
      i=0
      do 17 j=1,ic3
      do 16 nn=1,nw
      i=i+1
      wri(i)=r1*wg(nn)
   16 xri(i)=r1*xg(nn)+r2
   17 r2=r2+ci
      r1=quin+quon
      t=(rmaxr-r1)/quon**2
      ic3=mquado /nw
      ci=quon/ic3
      r1=0.5*ci
      r2=r1+quin
      do 162 j=1,ic3
      do 161 nn=1,nw
      i=i+1
      x=r1*xg(nn)+r2
      wri(i)=r1*wg(nn)*(1.+2*t*ramp(x-quin))
  161 xri(i)=x+t*(x-quin)*ramp(x-quin)
  162 r2=r2+ci
*---------------------------------------------------------------------
      write(kout,*)'inner quadrature points'
      write(kout,163)(xri(i),i=1,mquadi)
  163 format(' gaussian quadrature points at',8f8.4)
      write(kout,*)'outer quadrature points'
      write(kout,163) (xri(i),i=mquadi+1,mquad)
*---------------------------------------------------------------------
      return
      end



c @(#)gauss2.f 2.1    8/24/86 17:45:56

      subroutine gauss2(npts, kode, a, b, xs, wts)
************************************************************************

c *** see r.h. landau's program lpott, 1981
c'
c     subroutine to return a scalled set of gaussian-legendre points
c     *npts* is the number of points and must be .gt. 1
c     *kode* is a three digit  (decimal) number -- htu
c      *u* = 0 means *npts* must be in the set *numbrs*, else syserr is
c              called.  *npts* may be either a constant or variable.
c          = 1 means if *npts* isn*t in the set *numbrs*,use the largest
c              number .lt. *npts* and change *npts* accordingly.
c              *npts* must be a variable.
c          = 2 means if *npts* isn*t in the set *numbrs*, use the
c              smallest member .gt. *npts* and change *npts* accordingly
c              however if *npts*.gt.numbrs(istop) then use numbrs(istop)
c              *npts* must be a variable.
c          = 5 means use equal spaced points.       for odd *npts* they
c              will be simpson points and weights.  for even *npts* they
c              will have equal weights.
c              note that in this case *npts* may be any value .ge. 1.
c              note that they are still scalled according to *t*.
c      *t*        determines how the internal          points (on
c                 (-1, +1)) are scalled.  see the computed go to
c                 after stmnt no 109 for all the details.
c      *h* = 0 means scale as per t above.
c          = 1 means scale as per *t* above but also remove square
c              root at lower bound.  this removes integrable
c              1/sqrt(x-lower) discontinuities and gives better
c              values for  sqrt(x-lower) f(x) type integrands.
c     *a* and *b* are used in the scalling
c     *xs* will be the set of points,  must be real dimensioned .gt. n
c     *ws* will be the set of weights  must be real dimensioned .gt. n
c

      implicit real*8 (a-h, o-z)
      integer numbrs(13)
      real*8:: xs(npts),wts(npts)
      real*8 ::hh(2)
      integer ::t, u, base
      real*8 ::x1(59), w1(59), x2(59), w2(59)
c AMORO 4/2/04
!      real*8 ::gx(118), wt(118)
!      equivalence (x1(1),gx(1)),(w1(1),wt(1)),(x2(1),gx(60))
!     1 ,(w2(1),wt(60))

      real*8:: gx(118),wt(118)

      data istop /13/
      data numbrs / 2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 32, 40, 48 /

      data x1/.5773502691896261d0,.861136311594053d0,.339981043584856d0,
     *.9324695142031520d+0, .6612093864662651d+0, .2386191860831970d+0,
     *.9602898564975360d+0, .7966664774136270d+0, .5255324099163290d+0,
     *.1834346424956500d+0, .9739065285171720d+0, .8650633666889850d+0,
     *.6794095682990241d+0, .4333953941292470d+0, .1488743389816310d+0,
     *.9815606342467190d+0, .9041172563704750d+0, .7699026741943050d+0,
     *.5873179542866171d+0, .3678314989981800d+0, .1252334085114690d+0,
     *.9862838086968120d+0, .9284348836635741d+0, .8272013150697650d+0,
     *.6872929048116851d+0, .5152486363581540d+0, .3191123689278900d+0,
     *.1080549487073440d+0, .9894009349916501d+0, .9445750230732330d+0,
     *.8656312023878320d+0, .7554044083550030d+0, .6178762444026441d+0,
     *.4580167776572270d+0, .2816035507792590d+0, .0950125098376370d+0,
     *.9931285991850949d+0, .9639719272779138d+0, .9122344282513259d+0,
     *.8391169718222189d+0, .7463319064601507d+0, .6360536807265151d+0,
     *.5108670019508270d+0, .3737060887154195d+0, .2277858511416450d+0,
     *.0765265211334973d+0, .9951872199970213d+0, .9747285559713094d+0,
     *.9382745520027326d+0, .8864155270044010d+0, .8200019859739029d+0,
     *.7401241915785542d+0, .6480936519369756d+0, .5454214713888396d+0,
     *.4337935076260451d+0, .3150426796961633d+0, .1911188674736163d+0,
     *.0640568928626056d+0, .9972638618494816d+0/

      data x2/ .9856115115452683d0, .9647622555875064d0,
     *.9349060759377397d+0, .8963211557660522d+0, .8493676137325699d+0,
     *.7944837959679424d+0, .7321821187402896d+0, .6630442669302153d+0,
     *.5877157572407623d+0, .5068999089322293d+0, .4213512761306353d+0,
     *.3318686022821276d+0, .2392873622521370d+0, .1444719615827964d+0,
     *.0483076656877383d+0, .9982377097105592d+0, .9907262386994570d+0,
     *.9772599499837742d+0, .9579168192137917d+0, .9328128082786765d+0,
     *.9020988069688742d+0, .8659595032122595d+0, .8246122308333115d+0,
     *.7783056514265194d+0, .7273182551899270d+0, .6719566846141796d+0,
     *.6125538896679803d+0, .5494671250951283d+0, .4830758016861787d+0,
     *.4137792043716050d+0, .3419940908257584d+0, .2681521850072536d+0,
     *.1926975807013710d+0, .1160840706752552d+0, .0387724175060508d+0,
     *.9987710072524261d+0, .9935301722663507d+0, .9841245837228269d+0,
     *.9705915925462472d+0, .9529877031604309d+0, .9313866907065542d+0,
     *.9058791367155696d+0, .8765720202742478d+0, .8435882616243934d+0,
     *.8070662040294426d+0, .7671590325157403d+0, .7240341309238146d+0,
     *.6778723796326640d+0, .6288673967765137d+0, .5772247260839728d+0,
     *.5231609747222331d+0, .4669029047509584d+0, .4086864819907167d+0,
     *.3487558862921607d+0, .2873624873554555d+0, .2247637903946890d+0,
     *.1612223560688917d+0, .0970046992094627d+0, .0323801709628694d+0/

      data w1/1.0000d+0, 0.3478548451374540d+0,  0.6521451548625461d+0,
     *.1713244923791700d+0, .3607615730481390d+0, .4679139345726910d+0,
     *.1012285362903760d+0, .2223810344533740d+0, .3137066458778870d+0,
     *.3626837833783620d+0, .0666713443086880d+0, .1494513491505810d+0,
     *.2190863625159820d+0, .2692667193099960d+0, .2955242247147530d+0,
     *.0471753363865120d+0, .1069393259953180d+0, .1600783285433460d+0,
     *.2031674267230660d+0, .2334925365383550d+0, .2491470458134030d+0,
     *.0351194603317520d+0, .0801580871597600d+0, .1215185706879030d+0,
     *.1572031671581940d+0, .1855383974779380d+0, .2051984637212960d+0,
     *.2152638534631580d+0, .0271524594117540d+0, .0622535239386480d+0,
     *.0951585116824930d+0, .1246289712555340d+0, .1495959888165770d+0,
     *.1691565193950030d+0, .1826034150449240d+0, .1894506104550690d+0,
     *.0176140071391521d+0, .0406014298003869d+0, .0626720483341091d+0,
     *.0832767415767047d+0, .1019301198172404d+0, .1181945319615184d+0,
     *.1316886384491766d+0, .1420961093183820d+0, .1491729864726037d+0,
     *.1527533871307258d+0, .0123412297999872d+0, .0285313886289337d+0,
     *.0442774388174198d+0, .0592985849154368d+0, .0733464814110803d+0,
     *.0861901615319533d+0, .0976186521041139d+0, .1074442701159656d+0,
     *.1155056680537256d+0, .1216704729278033d+0, .1258374563468282d+0,
     *.1279381953467521d+0, .0070186100094701d+0/

      data w2/              .0162743947309057d+0, .0253920653092620d+0,
     *.0342738629130214d+0, .0428358980222267d+0, .0509980592623762d+0,
     *.0586840934785355d+0, .0658222227763618d+0, .0723457941088485d+0,
     *.0781938957870703d+0, .0833119242269467d+0, .0876520930044038d+0,
     *.0911738786957639d+0, .0938443990808046d+0, .0956387200792748d+0,
     *.0965400885147278d+0, .0045212770985332d+0, .0104982845311528d+0,
     *.0164210583819079d+0, .0222458491941669d+0, .0279370069800234d+0,
     *.0334601952825478d+0, .0387821679744720d+0, .0438709081856733d+0,
     *.0486958076350722d+0, .0532278469839368d+0, .0574397690993916d+0,
     *.0613062424929289d+0, .0648040134566010d+0, .0679120458152339d+0,
     *.0706116473912868d+0, .0728865823958040d+0, .0747231690579683d+0,
     *.0761103619006262d+0, .0770398181642480d+0, .0775059479784248d+0,
     *.0031533460523058d+0, .0073275539012763d+0, .0114772345792345d+0,
     *.0155793157229438d+0, .0196161604573555d+0, .0235707608393244d+0,
     *.0274265097083569d+0, .0311672278327981d+0, .0347772225647704d+0,
     *.0382413510658307d+0, .0415450829434647d+0, .0446745608566943d+0,
     *.0476166584924905d+0, .0503590355538545d+0, .0528901894851937d+0,
     *.0551995036999842d+0, .0572772921004032d+0, .0591148396983956d+0,
     *.0607044391658939d+0, .0620394231598927d+0, .0631141922862540d+0,
     *.0639242385846482d+0, .0644661644359501d+0, .0647376968126839d+0/

c
c     this *1 - epsilon* is for simpson points so that when we scale
c     to infinity we won*t divide by zero.
c

      data one / 1.00000000001d+0 /
c

      gx(1:59)=x1(1:59)
      wt(1:59)=w1(1:59)
      gx(60:118)=x2(1:59)
      wt(60:118)=w2(1:59)

!! TEST
      wts(:)=0d0
!! TEST

      t = kode/10
      u = kode-10*t
      h = t/10
      t = t-10*h
      if ((h .lt. 0) .or. (h .gt. 1)) go to 90
      if (npts .lt. 1) go to 95
      if (u .eq. 5) go to 400
      if ((u .lt. 0) .or. (u .gt. 2)) go to 90
      base = 0
      do 30 i = 1, istop
      if (numbrs(i) .eq. npts) go to 100
      if (numbrs(i) .gt. npts) go to 50
   30 base = base + (numbrs(i)+1)/2
c     npts .gr. last
      if (u .eq. 0) go to 95
      npts = numbrs(istop)
      base = base - (npts+1)/2
      go to 100
c
 50   if ((npts .le. 1) .or. (u    .eq. 0)) go to 95
      go to (55, 60), u
c       call gotoer
 55   npts = numbrs(i-1)
      base = base - (npts+1)/2
      go to 100
 60   npts = numbrs(i)
      go to 100
c
 90   print 91, kode
 91   format (32h0***gauss*** invalid ""kode"" --   , i20)
      stop
 95   print 96, npts, kode
 96   format(28h0gauss,invalid npts,ok kode=  ,2i20)

      stop
c
 100  do 300 n = 1, npts
      if (u .eq. 5) go to 108
      if (n .gt. (npts+1)/2) go to 105
      x = -gx(n+base)
      w = wt(n+base)
      go to 109
  105 nn = npts-n+1
      x = gx(nn+base)
      w = wt(nn+base)
      go to 109
  108 x = xs(n)
      w = wts(n)
c     are we removeing sqrae root at lower.
  109 if (h .eq. 1) w = (x+1) * w
      if (h .eq. 1) x = .5 * (x*(x+2)-1)
      go to(110, 120, 130, 140, 150,160      ), t
c      call gotoer
  160 go to 90
c
c     1, scale to (a, b) uniformly
  110 xs(n) = (b+a)/2d0 + (b-a)/2d0*x
      wts(n) = (b-a)/2d0 * w
!      write(*,*)'gauss2: n,xs(n),wts(n)',n,xs(n),wts(n)
      go to 300
c
c     2, scale to (0, infinity)
  120 xs(n) = a*(1+x)/(1-x)
      wts(n) = 2*a*w/(1-x)**2
!      write(*,*)'gauss2: n,xs(n),wts(n)',n,xs(n),wts(n)
      go to 300
c
c     3, scale to (-infinity, infinity)
  130 xs(n) = a*x/(1-x**2)
      wts(n) = a*(1+x**2)/(1-x**2)**2
!      write(*,*)'gauss2: n,xs(n),wts(n)',n,xs(n),wts(n)
      go to 300
c
c     4 scale to (b, infinity) so that for b = 0 we have case 2
  140 xs(n) = (a+2*b + a*x)/(1-x)
      wts(n) = 2*(b+a) * w/(1-x)**2
!       write(*,*)'gauss2: n,xs(n),wts(n)',n,xs(n),wts(n)
      go to 300
c
c     5, scale to (0,b) so that for b = infinity we have case 2
c     (almost)
  150 xs(n) = a*b * (1+x)/(b+a - (b-a)*x)
      wts(n) = 2*a*b**2 * w/(b+a - (b-a)*x)**2
!      write(*,*)'gauss2: n,xs(n),wts(n)',n,xs(n),wts(n)
      go to 300
  300 continue
       do 320 ii=1,npts
  320  continue
      return
c
c     uniform spaced points for any n.
c
  400 if ((npts/2)*2 .eq. npts) go to 450
c
c     odd number = sumpsons
      if (npts .gt. 1) go to 410
      xs(1) = 0
      wts(1) = 2
      go to 100
  410 h = 2. / (npts-1)
      hb3 = h/3
      xs(1) = -1
      xs(npts) = 1
      if (t .eq. 3) xs(1) = -one
      if ((t .ge. 2) .and. (t .le. 4)) xs(npts) = one
      wts(1) = hb3
      wts(npts) = hb3
      hh(1) = 4d0*hb3
      hh(2) = 2d0*hb3
      nb2 = npts/2
      ih = 1
      do 430 i = 1, nb2
      xs(i+1) = -1 + i*h
      xs(npts-i) = 1-i*h
      wts(i+1) = hh(ih)
      wts(npts-i) = hh(ih)
!       write(*,*)'gauss2: i,xs(i),wts(i)',i,xs(i),wts(i)
  430 ih = 3-ih
      go to 100
c
c     even number = equal weight
  450 h = 2./npts
      do 480 i = 1, npts
      xs(i) = -1d0-h/2d0 + i*h
480      wts(i) = h
!  480 write(*,*)'gauss2: i,xs(i),wts(i)',i,xs(i),wts(i)
!  480 wts(i) = h
      go to 100
      end



*************************************************************************
      subroutine sim2(fa,res,m,n,h,nramax)
*------------------------------------------------------------------------
*     subroutine does the integral of fa stored
*     in the array of the same name using simpsons rule. the step length
*     is h and the integral is between the elements m and n of the arrays
*     only. resulting integral is placed in res.
*------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      integer :: m,n,nramax
      real*8 ::res,fa(nramax),dq(nramax),h
      
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



c----------------------------------------------------------------------
c calculates Coulomb phase-shifts (by C. Bertulani)
c----------------------------------------------------------------------
      subroutine coul(eta,cph,lmax)
      implicit real*8(a-h,o-z)
!      include 'crossec.dim'
      real*8:: cph(0:lmax)
     
c
      sto=16.d0+eta*eta
c      write(*,*)'coul: eta,lmax=',eta,lmax
     
      if (lmax>10000) stop 'utils.f90=> In coul() lmax too large!'
 
      cph(0)=-eta+(eta/2.d0)*dlog(sto)+3.5d0*datan(eta/4.d0)-(
     1 datan(eta)+datan(eta/2.d0)+datan(eta/3.d0))-(eta/(12.d0*
     2 sto))*(1.d0+(1.d0/30.d0)*(eta**2-48.d0)/sto**2+(1.d0/105.d0)
     3 *(eta**4-160.d0*eta**2+1280.d0)/sto**4)
c
      do 1 ii=1,lmax
         fi=dble(ii)
         cph(ii)=cph(ii-1)+datan(eta/fi)
!         write(*,*)ii,cph(ii)
1     continue
      return
      end



      subroutine sigcl (eta,scou,lmax)
************************************************************************

c *** see r.h. landau''s program lpott, 1981
c
c     nb lmax must be variable

      implicit real*8 (a-h, o-z)
      parameter(ldim=30)
      dimension scou(ldim), ft(5)

c >>> first executable statement <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      abeta = abs(eta)
      if (lmax.gt.ldim) lmax = ldim
      if (abeta-1.e-06) 10,10,30
 10   do 20 i=1,lmax
         scou(i) = 0.e0
 20   continue
      go to 100
 30   c = 1.e0/(eta*eta+2601.e0)
      ch = sqrt(c)
      theta = atan(eta/51.e0)
      do 40 m=1,5
         ft(m) = sin(theta*(2*m-1))
 40   continue
      sig2 = 50.5 * theta - eta * ( log( ch ) + 1. )
      sig2 = sig2-ch*(ft(1)/12.-c*(ft(2)/360.-c*(ft(3)/1260.-c*(ft(4)/
     11680.-c*ft(5)/1188.))))
      lmax1 = 50-lmax
      if (lmax1) 70,70,50
 50   do 60 i=1,lmax1
         l = 51-i
         sig2 = sig2-atan(eta/l)
 60   continue
 70   lmax1 = lmax1+2
      scou(lmax) = sig2-atan(eta/lmax)
      if (lmax1-50) 80,80,100
 80   do 90 i=lmax1,50
         l = 51-i
         scou(l) = scou(l+1)-atan(eta/l)
 90   continue
 100  continue
      return
c**  this program valid on ftn4 and ftn5 **
      end


      subroutine sigcl3 (eta,scou,lmax)
************************************************************************

c *** see r.h. landau''s program lpott, 1981
c
      implicit real*8 (a-h, o-z)
      real*8 :: scou(0:lmax+1), ft(5)
      integer:: kout=16

c >>> first executable statement <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      abeta = abs(eta)
      scou=0.
!      if (lmax.gt.ldim) lmax = ldim
      if (abeta<1.e-06) return
 30   c = 1.e0/(eta*eta+2601.e0)
      ch = sqrt(c)
      theta = atan(eta/51.e0)

      do 40 m=1,5
         ft(m) = sin(theta*(2*m-1))
 40   continue

      sig2 = 50.5 * theta - eta * ( log( ch ) + 1. )
      sig2 = sig2-ch*(ft(1)/12.-c*(ft(2)/360.-c*(ft(3)/1260.-c*(ft(4)/
     11680.-c*ft(5)/1188.))))
      lmax1 = 50-lmax
      if (lmax1<=0) then
         do i=1,lmax1
            l = 51-i
            sig2 = sig2-atan(eta/l)
         enddo
      endif

 70   lmax1 = lmax1+2
      scou(lmax) = sig2-atan(eta/lmax)
!      if (lmax1-50) 80,80,100
      if (lmax1-50>0) return

 80   do 90 i=lmax1,50
         l = 51-i
         scou(l) = scou(l+1)-atan(eta/l)         
 90   continue
      do l=0,lmax-1
         scou(l)=scou(l+1)
      end do
      return
      end



******************************************************************
!     Coulomb S-matrix
******************************************************************
      subroutine scoul(som,sc,lmax)
        IMPLICIT NONE
        REAL*8 :: som,di,j,dj
        COMPLEX*16 :: sc(0:lmax),aux,cgamma,z,z1,z2
        INTEGER :: i,l,lmax
        
        do l=0,lmax
        dj=DBLE(l)
        z1=cmplx(1.d0+dj,som)
        z2=cmplx(1.d0+dj,-som)
        sc(l)=cgamma(z1)/cgamma(z2)
        enddo
        RETURN
      END 


*******************************************************************
C     ALGORITHM 404 COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN COMM. ACM, VOL. 14, NO. 01,
C     P. 048.
      FUNCTION CGAMMA(Z)
      implicit real*8 (a-h,o-z)
      real*8:: pir=3.14159265358
      COMPLEX*16 Z,ZM,T,TT,SUM,TERM,DEN,CGAMMA,PI,A
      DIMENSION C(12)
      LOGICAL REFLEK
C SET IOUT FOR PROPER OUTPUT CHANNEL OF COMPUTER SYSTEM FOR
C ERROR MESSAGES
      IOUT = 3
      PI = (3.14159265358979311,0.d0)
      zero=0.d0
      X = REAL(Z)
      Y = AIMAG(Z)
C TOL = LIMIT OF PRECISION OF COMPUTER SYSTEM IN SINGLE PRECISI
      TOL = 1.0E-7
      REFLEK = .TRUE.
C DETERMINE WHETHER Z IS TOO CLOSE TO A POLE
C CHECK WHETHER TOO CLOSE TO ORIGIN
      IF(X.GE.TOL) GO TO 20
C FIND THE NEAREST POLE AND COMPUTE DISTANCE TO IT
      XDIST = X-INT(X-.5)
      ZM = CMPLX(XDIST,Y)
      IF(CDABS(ZM).GE.TOL) GO TO 10
C IF Z IS TOO CLOSE TO A POLE, PRINT ERROR MESSAGE AND RETURN
C WITH CGAMMA = (1.E7,0.0E0)
      WRITE(IOUT,900) Z
      CGAMMA = (1.E7,0.E0)
      RETURN
C FOR REAL(Z) NEGATIVE EMPLOY THE REFLECTION FORMULA
C GAMMA(Z) = PI/(SIN(PI*Z)*GAMMA(1-Z))
C AND COMPUTE GAMMA(1-Z).  NOTE REFLEK IS A TAG TO INDICATE THA
C THIS RELATION MUST BE USED LATER.
10    IF(X.GE.0.0) GO TO 20
      REFLEK = .FALSE.
      Z = (1.0,0.0)-Z
      X = 1.0-X
      Y = -Y
C IF Z IS NOT TOO CLOSE TO A POLE, MAKE REAL(Z)>10 AND ARG(Z)<P
20    M = 0
40    IF(X.GE.10.) GO TO 50
      X = X + 1.0
      M = M + 1
      GO TO 40
50    IF(ABS(Y).LT.X) GO TO 60
      X = X + 1.0
      M = M + 1
      GO TO 50
60    T = CMPLX(X,Y)
      TT = T*T
      DEN = T
C COEFFICIENTS IN STIRLING*S APPROXIMATION FOR LN(GAMMA(T))
      C(1) = 1.d0/12.d0
      C(2) = -1.d0/360.d0
      C(3) = 1.d0/1260.d0
      C(4) = -1.d0/1680.d0
      C(5) = 1.d0/1188.d0
      C(6) = -691.d0/360360.d0
      C(7) = 1.d0/156.d0
      C(8) = -3617.d0/122400.d0
      C(9) = 43867.d0/244188.d0
      C(10) = -174611.d0/125400.d0
      C(11) = 77683.d0/5796.d0
!      SUM = (T-(.5d0,0.0d0))*CDLOG(T)-T+
!     &      CMPLX(.5d0*ALOG(2.*pir),zero)
      SUM = (T-(.5d0,0.0d0))*CDLOG(T)-T+
     &      CMPLX(.5d0*LOG(2.*pir),zero)
      J = 1
70    TERM = C(J)/DEN
C TEST REAL AND IMAGINARY PARTS OF LN(GAMMA(Z)) SEPARATELY FOR
C CONVERGENCE.  IF Z IS REAL SKIP IMAGINARY PART OF CHECK.
      IF(ABS(REAL(TERM)/REAL(SUM)).GE.TOL) GO TO 80
      IF(Y.EQ.0.0) GO TO 100
      IF(ABS(AIMAG(TERM)/AIMAG(SUM)).LT.TOL) GO TO 100
80    SUM = SUM + TERM
      J = J + 1
      DEN = DEN*TT
C TEST FOR NONCONVERGENCE
      IF(J.EQ.12) GO TO 90
      GO TO 70
C STIRLING*S SERIES DID NOT CONVERGE.  PRINT ERROR MESSAGE AND
C PROCEDE.
90    WRITE(IOUT,910) Z
C RECURSION RELATION USED TO OBTAIN LN(GAMMA(Z))
C LN(GAMMA(Z)) = LN(GAMMA(Z+M)/(Z*(Z+1)*...*(Z+M-1)))
C = LN(GAMMA(Z+M)-LN(Z)-LN(Z+1)-...-LN(Z+M
100   IF(M.EQ.0) GO TO 120
      DO 110 I = 1,M
      A = CMPLX(I*1.d0-1.d0,0.0d0)
110   SUM = SUM-CDLOG(Z+A)
C CHECK TO SEE IF REFLECTION FORMULA SHOULD BE USED
120   IF(REFLEK) GO TO 130
      SUM = CDLOG(PI/CDSIN(PI*Z))-SUM
      Z = (1.0d0,0.0d0) -Z
130   CGAMMA = CDEXP(SUM)
      RETURN
900   FORMAT(1X,2E14.7,10X,49HARGUMENT OF GAMMA FUNCTION IS TOO CLOSE TO
     1 A POLE)
910   FORMAT(44H ERROR - STIRLING*S SERIES HAS NOT CONVERGED/14X,4HZ = ,
     12E14.7)
      END





************************************************************************
      subroutine sim(fa,res,m,n,h)
*-----------------------------------------------------------------------
*     subroutine does the integral of fa stored
*     in the arrays of the same name using simpsons rule. the step lengt
*     is h and the integral is between the elements m and n of the array
*     only. resulting integral is placed in res.
*-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension fa(5002),dq(5002)
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



!******************************************************************
!*    Complex Bessel Function of the 1st Kind of integer order    *
!* -------------------------------------------------------------- *
!* SAMPLE RUN:                                                    *
!*                                                                *
!* Complex Bessel Function of the 1st Kind of integer order       *
!*                                                                * 
!* Input complex argument (real imaginary): 1 2                   *
!* Input integer order: 1                                         *
!*                                                                *
!* Function value:  (1.291847,1.010488)                           *
!*                                                                *
!*                                                                *
!*                          F90 Release 1.0 By J-P Moreau, Paris. *
!******************************************************************
!!$Program Test_CBESSJ
!!$
!!$integer nu
!!$complex z,z1
!!$real*8  x,y
!!$
!!$  print *,' '
!!$  print *,' Complex Bessel Function of the 1st Kind of integer order'
!!$  print *,' '
!!$  write(*,10,advance='no'); read *, x, y
!!$  z = CMPLX(x,y)
!!$  write(*,20,advance='no'); read *, nu
!!$  print *,' '
!!$
!!$  call CBESSJ(z,nu,z1)
!!$
!!$  print *,' Function value: ', z1
!!$  print *,' '
!!$
!!$10 format('  Input complex argument (real imaginary): ')
!!$20 format('  Input integer order: ')
!!$
!!$END

      real*8 Function Fact(K)
        Integer i
        Real*8  f
        F=1.d0
        do i=2, k 
           f=f*dfloat(i)
        end do
        Fact=f
        return
      End Function Fact

!*******************************************
!*           FUNCTION  GAMMA(X)            *
!* --------------------------------------- *
!* Returns the value of Gamma(x) in double *
!* precision as EXP(LN(GAMMA(X))) for X>0. *
!*******************************************
      real*8 Function Gamma(xx)
        parameter(ONE=1.d0,FPF=5.5d0,HALF=0.5d0)
        real*8 xx
        real*8 cof(6)
        real*8 stp,x,tmp,ser
        integer j
        cof(1)=76.18009173d0
        cof(2)=-86.50532033d0
        cof(3)=24.01409822d0
        cof(4)=-1.231739516d0
        cof(5)=0.120858003d-2
        cof(6)=-0.536382d-5
        stp=2.50662827465d0

        x=xx-ONE
        tmp=x+FPF
        tmp=(x+HALF)*LOG(tmp)-tmp
        ser=ONE
        do j=1, 6
           x=x+ONE
           ser=ser+cof(j)/x
        end do
        Gamma = EXP(tmp+LOG(stp*ser))
        return
      End Function Gamma






c  FROM NUMERICAL RECIPES
      FUNCTION bessj(n,x) 
      INTEGER n,IACC 
      REAL*8 bessj,x,BIGNO,BIGNI 
      PARAMETER (IACC=40,BIGNO=1.e10,BIGNI=1.e-10) 
C     USES bessj0,bessj1 
c     Returns the Bessel function Jn(x) for any real x and n.gt.2. 
      INTEGER j,jsum,m 
      REAL*8 ax,bj,bjm,bjp,sum,tox,bessj0,bessj1 
      if(n.lt.2)pause 'bad argument n in bessj' 
      ax=abs(x) 
      if(ax.eq.0.)then 
         bessj=0. 
      else if(ax.gt.float(n))then !Upwards recurrence from J0 and J1. 
         tox=2./ax 
         bjm=bessj0(ax) 
         bj=bessj1(ax) 
      do 11 j=1,n-1 
         bjp=j*tox*bj-bjm 
         bjm=bj 
         bj=bjp 
11    enddo
      bessj=bj 
      else 
!     Downwards recurrence from an even m here computed. 
!     Make IACC larger to increase accuracy .
         tox=2./ax 
         m=2*((n+int(sqrt(float(IACC*n))))/2) 
         bessj=0. 
         jsum=0 
!     jsum will alternate between 0 and 1; when it is 1, we 
!     accumulate in sum the even terms in (5.5.16). 
         sum=0.
         bjp=0. 
         bj=1. 
      do 12 j=m,1,-1 !The downward recurrence. 
         bjm=j*tox*bj-bjp 
         bjp=bj 
         bj=bjm 
         if(abs(bj).gt.BIGNO)then !Renormalize to prevent overflows. 
            bj=bj*BIGNI
            bjp=bjp*BIGNI 
            bessj=bessj*BIGNI 
            sum=sum*BIGNI 
         endif 
         if(jsum.ne.0)sum=sum+bj !Accumulate the sum. 
         jsum=1-jsum            !Change 0 to 1 or vice versa. 
         if(j.eq.n)bessj=bjp    !Save the unnormalized answer. 
12      enddo 
      sum=2.*sum-bj 
!Compute (5.5.16) and use it to normalize the answer. 
      bessj=bessj/sum 
      endif 
      if(x.lt.0..and.mod(n,2).eq.1)bessj=-bessj 
      return 
      END 


c *** -------------------------------------------------------
c *** BESSJ_1(x) (from Numerical Recipes)
c *** -------------------------------------------------------
      FUNCTION bessj0(x) 
      REAL*8 bessj0,x 
c     Returns the Bessel function J0(x) for any real x. 
      REAL*8 ax,xx,z 
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4, 
     *  r5,r6,s1,s2,s3,s4,s5,s6,y 
! We will accumulate polynomials in double precision. 
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6, 
     * s1,s2,s3,s4,s5,s6 
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4, 
     * -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1, 
     * .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/ 
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     * 651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/, 
     * s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0, 
     * 9494680.718d0,59272.64853d0,267.8532712d0,1.d0/ 
      if(abs(x).lt.8.)then      !Direct rational function fit 
         y=x**2 
         bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) 
     *  /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))))) 
      else !Fitting function (6.5.9). 
         ax=abs(x) 
         z=8./ax 
         y=z**2 
         xx=ax-.785398164 
         bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     *  *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5))))) 
      endif 
      return 
      END 


c *** -------------------------------------------------------
c *** BESSJ_1(x) (from Numerical Recipes)
c *** -------------------------------------------------------
      FUNCTION bessj1(x) 
      REAL*8 bessj1,x ,uno
c      Returns the Bessel function J1(x) for any real x. 
      REAL*8 ax,xx,z 
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,
     &  q5,r1,r2,r3,r4, 
     &  r5,r6,s1,s2,s3,s4,s5,s6,y 
c We will accumulate polynomials in double precision. 
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6, 
     * s1,s2,s3,s4,s5,s6 
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     *  242396853.1d0, 
     * -2972611.439d0,15704.48260d0,-30.16036606d0/, 
     * s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0, 
     * 18583304.74d0,99447.43394d0,376.9991397d0,1.d0/ 
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     *  .2457520174d-5, 
     * -.24033702d-6/
      DATA q1,q2,q3,q4,q5/.04687499995d0,
     * -.2002690873d-3, 
     *  .8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then !Direct rational approximation. 
      y=x**2 
      bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) 
     *   /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))))) 
      else !Fitting function (6.5.9). 
         ax=abs(x) 
         z=8./ax 
         y=z**2 
         xx=ax-2.356194491
         bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y 
     *  *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5))))) 
     *  *sign(uno,x) 
      endif 
      return 
      END 





      SUBROUTINE sphbes(n,x,sj,sy,sjp,syp)
      INTEGER n
      REAL*8 sj,sjp,sy,syp,x
C     USES bessjy
c      Returns spherical Bessel functions jn(x), yn(x), 
c      and their derivatives j0n(x), y0n(x) for integer n.
      REAL*8 factor,order,rj,rjp,ry,ryp,RTPIO2
      PARAMETER (RTPIO2=1.2533141)
      if(n.lt.0.or.x.le.0.)pause 'bad arguments in sphbes'
      order=n+0.5
c     all bessjy(x,order,rj,ry,rjp,ryp)
      factor=RTPIO2/sqrt(x)
      sj=factor*rj
      sy=factor*ry
      sjp=factor*rjp-sj/(2.*x)
      syp=factor*ryp-sy/(2.*x)
      return
      END


c **************************************************************
      function hat(lh)
c *********************************************************
      implicit real*8(a-h,o-z)
      hat = sqrt(2.*lh+1.)
      return
      end



      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp) 
      INTEGER MAXIT 
      REAL rj,rjp,ry,ryp,x,xnu,XMIN 
      DOUBLE PRECISION EPS,FPMIN,PI 
      PARAMETER (EPS=1.e-10,FPMIN=1.e-30,MAXIT=10000,XMIN=2., 
     *     PI=3.141592653589793d0) 
C     USES beschb 
c      Returns the Bessel functions rj = J , ry = Y and their derivatives rjp = J0
c      , ryp = Y0, 
c      for positive x and for xnu =   0. The relative accuracy is within one or two significant 
c      digits of EPS, except near a zero of one of the functions, where EPS controls its absolute 
c      accuracy. FPMIN is a number close to the machine''s smallest floating-point number. All 
c      internal arithmetic is in double precision. To convert the entire routine to double precision, 
c      change the REAL declaration above and decrease EPS to 10−16. Also convert the subroutine 
c      beschb. 
      INTEGER i,isign,l,nl 
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli, 
     *     dr,e,f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h, 
     *     p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1, 
     *     rymu,rymup,rytemp,sum,sum1,temp,w,x2,xi,xi2,xmu,xmu2 
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessjy' 
      if(x.lt.XMIN)then 
!nl is the number of downward recurrences of the J's and 
!      upward recurrences of Y 's. xmu lies between −1=2 and 
!      1/2 for x < XMIN, while it is chosen so that x is greater 
!      than the turning point for x ge XMIN. 
      nl=int(xnu+.5d0) 
      else
         nl=max(0,int(xnu-x+1.5d0)) 
      endif 
      xmu=xnu-nl 
      xmu2=xmu*xmu 
      xi=1.d0/x 
      xi2=2.d0*xi 
      w=xi2/PI The Wronskian
      isign=1 
!Evaluate CF1 by modified Lentz''s method.isign keeps 
!      track of sign changes in the denominator
      h=xnu*xi 
      if(h.lt.FPMIN)h=FPMIN 
      b=xi2*xnu 
      d=0.d0 
      c=h 
      do 11 i=1,MAXIT 
      b=b+xi2 
      d=b-d 
      if(abs(d).lt.FPMIN)d=FPMIN 
      c=b-1.d0/c 
      if(abs(c).lt.FPMIN)c=FPMIN 
      d=1.d0/d 
      del=c*d 
      h=del*h 
      if(d.lt.0.d0)isign=-isign 
      if(abs(del-1.d0).lt.EPS)goto 1 
11      enddo 
      pause 'x too large in bessjy; try asymptotic expansion' 
 1    continue 
      rjl=isign*FPMIN 
!Initialize J and J0
!      for downward recurrence. 
      rjpl=h*rjl 
      rjl1=rjl ! Store values for later rescaling. 
      rjp1=rjpl 
      fact=xnu*xi 
      do 12 l=nl,1,-1 
      rjtemp=fact*rjl+rjpl 
      fact=fact-xi 
      rjpl=fact*rjtemp-rjl 
      rjl=rjtemp 
12      enddo  
      if(rjl.eq.0.d0)rjl=EPS 
      f=rjpl/rjl !Now have unnormalized J and J0 
      if(x.lt.XMIN) then !Use series. 
      x2=.5d0*x 
      pimu=PI*xmu 
      if(abs(pimu).lt.EPS)then 
         fact=1.d0 
      else
         fact=pimu/sin(pimu) 
      endif 
      d=-log(x2) 
      e=xmu*d 
      if(abs(e).lt.EPS)then 
         fact2=1.d0 
      else
         fact2=sinh(e)/e 
      endif 
!Chebyshev evaluation of Gamma_1 and Gamma_2. 
      call beschb(xmu,gam1,gam2,gampl,gammi) 
      ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d) !f0. 
      e=exp(e) 
      p=e/(gampl*PI) !p0. 
      q=1.d0/(e*PI*gammi) !q0. 
      pimu2=0.5d0*pimu 
      if(abs(pimu2).lt.EPS)then 
         fact3=1.d0 
      else
         fact3=sin(pimu2)/pimu2 
      endif 
      r=PI*pimu2*fact3*fact3 
      c=1.d0 
      d=-x2*x2 
      sum=ff+r*q 
      sum1=p
      do 13 i=1,MAXIT 
      ff=(i*ff+p+q)/(i*i-xmu2) 
      c=c*d/i 
      p=p/(i-xmu) 
      q=q/(i+xmu) 
      del=c*(ff+r*q) 
      sum=sum+del 
      del1=c*p-i*del 
      sum1=sum1+del1 
      if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2 
 13     enddo  
      pause 'bessy series failed to converge' 
 2    continue 
      rymu=-sum 
      ry1=-sum1*xi2 
      rymup=xmu*xi*rymu-ry1 
      rjmu=w/(rymup-f*rymu) !Equation (6.7.13). 
      else !Evaluate CF2 by modified Lentz''s method 
         a=.25d0-xmu2 !(x5.2). 
         p=-.5d0*xi 
         q=1.d0 
         br=2.d0*x 
         bi=2.d0 
         fact=a*xi/(p*p+q*q) 
      cr=br+q*fact 
      ci=bi+p*fact 
      den=br*br+bi*bi 
      dr=br/den 
      di=-bi/den 
      dlr=cr*dr-ci*di 
      dli=cr*di+ci*dr 
         temp=p*dlr-q*dli 
         q=p*dli+q*dlr 
         p=temp 
      do 14 i=2,MAXIT 
         a=a+2*(i-1) 
         bi=bi+2.d0 
      dr=a*dr+br 
      di=a*di+bi 
         if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN 
         fact=a/(cr*cr+ci*ci) 
      cr=br+cr*fact 
      ci=bi-ci*fact 
         if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN 
      den=dr*dr+di*di 
      dr=dr/den 
      di=-di/den 
      dlr=cr*dr-ci*di 
      dli=cr*di+ci*dr 
         temp=p*dlr-q*dli 
         q=p*dli+q*dlr 
         p=temp 
         if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3 
 14     enddo  
      pause 'cf2 failed in bessjy' 
 3    continue 
      gam=(p-f)/q !Equations (6.7.6) { (6.7.10). 
      rjmu=sqrt(w/((p-f)*gam+q)) 
      rjmu=sign(rjmu,rjl) 
      rymu=rjmu*gam 
      rymup=rymu*(p+q/gam) 
      ry1=xmu*xi*rymu-rymup 
      endif 
      fact=rjmu/rjl
      rj=rjl1*fact !Scale original J and J0 
      rjp=rjp1*fact 
      do 15 i=1,nl !Upward recurrence of Y 
      rytemp=(xmu+i)*xi2*ry1-rymu 
      rymu=ry1 
      ry1=rytemp 
15    enddo
      ry=rymu 
      ryp=xnu*xi*rymu-ry1 
      return 
      END 
      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi) 
      INTEGER NUSE1,NUSE2 
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x 
      PARAMETER (NUSE1=5,NUSE2=5) 
C     USES chebev 
!      Evaluate Gamma_1 and Gamma_2 by Chebyshev expansion for |x|.le.1/2. 
!      Also returns 1/Gamma(1 + x) and 
!      1/Gamma(1 − x). If converting to double precision, set NUSE1=7, NUSE2=8. 
      REAL xx,c1(7),c2(8),chebev 
      SAVE c1,c2 
      DATA c1/-1.142022680371168d0,6.5165112670737d-3, 
     *     3.087090173086d-4,-3.4706269649d-6,6.9437664d-9, 
     *     3.67795d-11,-1.356d-13/ 
      DATA c2/1.843740587300905d0,-7.68528408447867d-2, 
     *     1.2719271366546d-3,-4.9717367042d-6,-3.31261198d-8, 
     *     2.423096d-10,-1.702d-13,-1.49d-15/ 
      xx=8.d0*x*x-1.d0 
!Multiply x by 2 to make range be −1 to 1, and then 
!      apply transformation for evaluating even Chebyshev series. 
      gam1=chebev(-1.,1.,c1,NUSE1,xx) 
      gam2=chebev(-1.,1.,c2,NUSE2,xx) 
      gampl=gam2-x*gam1 
      gammi=gam2+x*gam1 
      return 
      END 



      FUNCTION chebev(a,b,c,m,x)
        INTEGER m
        REAL*8 chebev,a,b,x,c(m)
c Chebyshev evaluation: All arguments are input. c(1:m) is an array of Chebyshev 
c coeficients , the first m elements 
c of c output from chebft (which must have been called with
c the same a and b). The Chebyshev polynomial 
c Pmk=1 ckTk 1(y)  c1=2 is evaluated at a
c point y = [x-(b+a)=2]=[(b- a)=2], and the result is returned as the function value.
        INTEGER j
        REAL*8 d,dd,sv,y,y2
        if ((x-a)*(x-b).gt.0.) pause 'x not in range in chebev'
        d=0.
        dd=0.
        y=(2.*x-a-b)/(b-a) !Change of variable.
        y2=2.*y
        do 11 j=m,2,-1 !Clenshaw  recurrence.
           sv=d
           d=y2*d-dd+c(j)
           dd=sv
11      enddo
        chebev=y*d-dd+0.5*c(1) !Last step is diferent.
        return
       END 














