******************************************************************
*
*  Bonn potential: includes obnn, obnna, obnnt
*
******************************************************************
      subroutine potsljbonn
        use toff
        use crdwrt
        use cpot
        use cstate
c
c        program calls the simplest, energy-independent potential
c        Bonn OBEP-q  (O83)
c
       implicit real*8 (a-h,o-z)
c      parameter (kmx=100)
c      dimension vl(kmx,kmx),v0(12,7,2),xmu(12,2),
       real*8, allocatable :: q(:)
c        arguments and values of potential
c

c        this has been the end of the arguments and values of the
c        potential subroutine
c
c        dimensions and specifications for this interface
c        routine with Joe Redish's k-matrix code '
c

      logical switch
c
c
      data switch/.false./
      data pih/1.570796326794897d0/
      data uf/197.3286d0/
c
c
      allocate(q(kmax))


      kread=5
      kwrite=6
      open(unit=kread,file='do83t.dat',status='unknown')
c        the following parameter is actually not needed
      kpunch=7
c
      ufd=1.d0/uf
      uf3=uf*uf*uf
c
c
c        the following statement means, that the potential
c        will be given in the helicity-formalism or not:
c        potential will be called in lsj states
      heform=.false.
 
      sing=.true.
      trip=.true.
      coup=.true.
c        the logical parameter endep is set by the program obnnnn
c        according to the fact, that the potential is either energy-
c        dependent or not (in case of hm1 or hm3: endep=.false.)
c
c
      n1=kmax+1
      n=n1-1
      n2=n2*2

            allocate(vl(2*n1,2*n1))
c
      j=jj
c
c
c        get momenta in mev
      do 100 i=1,n1
  100 q(i)=qq(i)*uf
c
c        factor for correct dimension of V (fm2)
       fac=uf*uf*pih
c
c
      do 395 ix=1,n1
      xmev=q(ix)
c
c
      do 395 iy=ix,n1
      ymev=q(iy)
c
      call obnnt
c
c
      if (is.eq.0) then
c
      vl(ix,iy)=v(1)*fac
c
      else if (is.eq.1.and.ic.eq.1) then
c
      vl(ix,iy)=v(2)*fac
c        special case of 3p0
      if (j.eq.0) vl(ix,iy)=v(3)*fac
c
      else if (is.eq.1.and.ic.eq.2) then
c
c        storage of channels:
c        v(3)=++ / v(4)=--
c        v(5)=+- / v(6)=-+
c
      vl(ix,iy)=v(3)*fac
      vl(ix+n1,iy+n1)=v(4)*fac
      vl(ix,iy+n1)= v(5)*fac
      vl(ix+n1,iy)= v(6)*fac
c
      endif
 
  395 continue
c
c        fill up remaining half of matrix
c
      if (ic.eq.1) then
c
      do 500 ix=1,n1
      do 500 iy=ix+1,n1
  500 vl(iy,ix)=vl(ix,iy)
c
      else if (ic.eq.2) then
c
      do 550 ix=1,n1
      do 550 iy=ix+1,n1
      vl(iy,ix)=vl(ix,iy)
      vl(iy+n1,ix+n1)=vl(ix+n1,iy+n1)
      vl(iy,ix+n1)=vl(ix+n1,iy)
  550 vl(iy+n1,ix)=vl(ix,iy+n1)
c
      endif
c
      return
      end
 
*****************************************************************
      subroutine obnnt
      use crdwrt
      use cpot
      use cstate
      use cpoted
c
c
c        this is the code for the bonn potential -
c        model type: momentum space one-boson-exchange potential
c
c        [changed for Sparc/Sun Feb. 1991, ansi-standard]
c
c
c
      implicit double precision (a-h,o-z)
!      common /crdwrt/ kread,kwrite,kpunch,kda(9)
!      common /cpot/ v(6),xmev,ymev
!      common /cstate/ j,heform,sing,trip,coup,endep,label
!      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,
!     1                znrl,zrel,smev,noced
!      logical heform,sing,trip,coup,endep
!      logical noced
!      character*4 label
      logical nnoced
      dimension vv(6)
      logical index,indexa
      data index/.false./,indexa/.false./
c
c
      nnoced=noced
c
c
c        t = 1 potential
c
c
      call obnn
c
c
      do 105 iv=1,6
  105 vv(iv)=v(iv)
      if (index) go to 200
      index=.true.
      write (kwrite,11000)
11000 format (//' these were the parameters for the t=1 potential.')
c
c
  200 continue
c
c
c        t = 0 potential
c
c
      noced=nnoced
c
      call obnna
c
      if (indexa) go to 300
      indexa=.true.
      write (kwrite,12000)
12000 format (//' these were the parameters for the t=0 potential.')
c
c
  300 if (mod(j,2).eq.1) go to 350
c
c
c        j is even
c
      v(1)=vv(1)
      do 305 iv=3,6
  305 v(iv)=vv(iv)
      go to 1000
c
c
c        j is odd
c
  350 v(2)=vv(2)
c
c
 1000 return
      end


**********************************************************************
      subroutine obnn
        use crdwrt
        use cpot
        use cstate
        use cob
c
c
c        one-boson-exchange nn-nn interaction;
c        version which uses numerical integration
c
c        [changed for Sparc/Sun Feb. 1991, ansi-standard]
c
c
      implicit double precision (a-h,o-z)
c
c
!      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c        arguments and values of this subroutine
c
!      common /cpot/   v(6),xmev,ymev
!      common /cstate/ j,heform,sing,trip,coup,endep,label
c
!      character*4 label
c
c        this has been the end of the common-blocks containing
c        the arguments and values of this subroutine in the case of
c        no energy-dependence of the potential;
c        in case of energy-dependence look for the common-block /cped/
c        in obai and obaa.
c
c        specifications for these two common blocks
c
!      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
!      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
!     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
!     2                ez1,ez2,ct(96),wt(96),
!     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
!     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
!     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
!     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
!      logical indc,indxy,indpar
c
c
c        further specifications
c
      logical index,index2
      logical indmg(12)
      character*4 mesong(12)
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1                   '1-  ','1-t ','1-tt','1-st','1-ss',
     2                    '1+  ','2+  '/
c
      data index/.false./,index2/.false./
      data indmg/12*.false./
      data pi/3.141592653589793d0/
      data  d3/0.3333333333333333d0/
      data td3/0.6666666666666667d0/
      data fd3/1.3333333333333333d0/
c
c
c        statement functionc needed for cray fortran
c
c      sqrt(x)=sqrt(x)
c
c
c
c
      inter=1
c
c
c
c
c        call subroutine obpar once and only once
c
c
      if (index) go to 50
      index=.true.
c      if (index2) go to 40
c      index2=.true.
c      indpar(inter)=.false.
c   40 if (indpar(inter)) go to 45
c
c
      call obpar
c
c
c   45 wdd(inter)=0.d0
      wdd(inter)=0.d0
      iftgo=ift(inter)+1
      dwn=1.d0/wnn(inter)
      iman=imaa(inter)
      imen=imea(inter)
c
c
c        prepare constant over-all factor
c
      fac=1.d0/(2.d0*pi)*dwn*dwn
c     --------------------------
c
c
c
c
c
c
c
c        prepare expressions depending on x and y
c        ----------------------------------------
c        ----------------------------------------
c
c
c
c
   50 xa=xmev*dwn
      ya=ymev*dwn
      indxy=.false.
      x=xa
      xx=x*x
      y=ya
      yy=y*y
      xy2=x*y*2.d0
      xxpyy=xx+yy
      ex=sqrt(1.d0+xx)
      ey=sqrt(1.d0+yy)
      eem12=(ex*ey-1.d0)*2.d0
c
c
c
c
   55 xy=xy2*0.5d0
c     dxy=1.d0/xy
      ee=ex*ey
      ree=sqrt(ee)
      eem1=ee-1.d0
      eme=ex-ey
      emeh=eme*0.5d0
      emehq=emeh*emeh
      eep1=ee+1.d0
       epe=ex+ey
      xxyy=xx*yy
c
c
c
c
c        prepare over-all factor
c
c
      go to (70,71,72),iftgo
c
c        no additional factor
c
   70 fff=fac
      go to 90
c
c        minimal relativity
c
   71 fff=fac/ree
      go to 90
c
c        factor m/e*m/e
c
   72 fff=fac/ee
c
c
c
c
c
c
   90 do 93 iv=1,6
   93 v(iv)=0.d0
      do 95 il=iman,imen
      do 95 iv=1,32
   95 vj(iv,il)=0.d0
c
c
c
c
c        contributions of mesons
c        -----------------------
c        -----------------------
c
c
c
c
      do 1995 img=1,mge
      mg=mggo(img,inter)
      if (mg.eq.0) go to 2000
      me=mgg(mg,inter)
      go to (100,200,300,400,500,600,600,600,900,1000,1100,1200),mg
c
c
c
c        0-  , pseudo-scalar mesons
c        --------------------------
c
c
c
c
  100 mc=1
c
      ff=1.d0
      f(1)=eem1
      f(2)=-xy
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-eme
      f(8)=-f(7)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c        0-t , pseudo-vector mesons
c        --------------------------
c
c
c
c
  200 mc=1
c
      ff=1.d0
      f(1)=eem1+emehq*(ee+3.d0)
      f(2)=-xy+emehq*xy
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-eme-eme*(emehq+eem1)
      f(8)=-f(7)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c        0-st, pseudo-scalar mesons in static limit
c        ------------------------------------------
c
c
c
c
  300 mc=1
c
      ff=1.d0
      f(1)=xxpyy*0.5d0
      f(2)=-xy
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-(xx-yy)*0.5d0
      f(8)=-f(7)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c
c        0+  , scalar mesons
c        -------------------
c
c
c
c
  400 mc=1
c
      ff=1.d0
      f(1)=-eep1
      f(2)=xy
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=epe
      f(8)=f(7)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c
c        0+st, scalar mesons in static limit
c        -----------------------------------
c
c
c
c
  500 mc=1
c
c
c
c        central term  '1'  only
c
c
c
      ff=1.d0
      f(1)=-2.d0
      f(2)=0.d0
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-f(1)
      f(8)=f(7)
      f(9)=f(2)
      f(10)=f(2)
      f(11)=f(2)
c
c
c
c        central term  ' k**2 + p**2 '
c
c
c
      f(2)=f(2)+xy
      f(10)=f(10)+xy
      f(11)=f(11)-xy
c
c
c
c        spin orbit
c
c
c
      f(4)=f(4)+xy
      f(5)=f(5)+xy
      f(10)=f(10)-xy
      f(11)=f(11)+xy
c
      call obstr(2,1,me)
c
c
c
c        case of additional sigma-l
c
c
c
c     mc=-1
c     xxyy8=xxyy/8.d0
c     f(1)=-xxyy8
c     f(2)=0.d0
c     f(3)=f(1)
c     f(4)=f(2)
c     f(5)=f(2)
c     f(6)=xxyy8
c     f(7)=f(1)
c     f(8)=f(7)
c     f(9)=f(6)*2.d0
c
c     call obstr(4,1,me)
      go to 1995
c
c
c
c
c        1-  , vector mesons
c        -------------------
c
c
c
c
c        vector coupling
c
c
c
c
  600 mc=1
c
      ff=2.d0
      f(1)=eem1+ee
      f(2)=0.d0
      f(3)=ee
      f(4)=xy
      f(5)=xy2
      f(6)=1.d0
      f(7)=-ey
      f(8)=-ex
c
      call obstr(1,1,me)
      if (mg.eq.7) go to 720
      if (mg.eq.8) go to 810
c
c
c
c
c        tensor coupling
c
c
c
c
      mc=2
c
c
      ff=0.5d0
        e1=-xxpyy+xxyy
      f(1)=eem1*(xxpyy+6.d0)+e1
      f(2)=-(xxpyy+4.d0)*xy
      f(3)=eem1*(xxpyy+2.d0)+e1
        e2=xxpyy+ee
      f(4)=-(1.d0+e2)*xy
      f(5)= (3.d0-e2)*xy
      f(6)=eem1*(xxpyy-2.d0)-xxpyy
        e3=2.d0*eme
        e4=yy*(ey-2.d0*ex)+xx*(ex-2.d0*ey)
      f(7)= e3-e4
      f(8)=-e3-e4
c        factors for additional terms
      f(9)=-xxyy
      f(10)=(1.d0+ee)*xy
      f(11)=-epe*xy
c
      call obstr(2,1,me)
c
c
c
c
c        vector-tensor coupling
c
c
c
c
      mc=3
c
      ff=2.d0
      f(1)=3.d0*eem1-xxpyy
      f(2)=-xy
      f(3)=eem1-xxpyy
      f(4)=-f(2)
      f(5)=3.d0*xy
      f(6)=-(eem1+xxpyy)
        e1=yy*ex+xx*ey
      f(7)= eme+e1
      f(8)=-eme+e1
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c
c        1-t , vector mesons a la gtg
c        ----------------------------
c
c
c
c
c        tensor coupling
c
c
c
c
  720 mc=2
c
      ff=0.25d0
      f(1)=(3.d0*ee+1.d0)*xxpyy
      f(2)=-(6.d0*ee+2.d0-xxpyy)*xy
      f(3)=eem1*xxpyy+4.d0*xxyy
      f(4)=-(4.d0*ee+xxpyy)*xy
      f(5)=(4.d0-3.d0*xxpyy)*xy
      f(6)=6.d0*xxyy-(ee+3.d0)*xxpyy
      f(7)=(ex+3.d0*ey)*xx+eme*yy
      f(8)=(ey+3.d0*ex)*yy-eme*xx
c        factors for additional terms
      f(9)=-2.d0*xxyy
      f(10)=eep1*xy2
      f(11)=-epe*xy2
c
      call obstr(2,1,me)
c
c
c
c
c        vector-tensor coupling
c
c
c
c
      mc=3
c
      ff=1.d0
      f(1)=xxpyy
      f(2)=-xy2
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=6.d0*xy
      f(6)=3.d0*f(3)
      f(7)=(ex*yy+3.d0*ey*xx)
      f(8)=(ey*xx+3.d0*ex*yy)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c
c        1-tt, vector mesons with consequent ignoration of retardation
c        -------------------------------------------------------------
c
c
c
c
c        vector coupling (additional term)
c
c
c
c
  810 mc=0
c
      f(1)=eep1
      f(2)=xy
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-epe
      f(8)=f(7)
        e1=eme*eme
c
      do 815 mx=1,me
      ff=e1/c(4,ima(mx,mg,inter))
  815 call obstr(1,mx,mx)
c
c
c
c
c
c
c        tensor coupling and vector-tensor coupling
c
c
c
c
      go to 720
c
c
c
c
c        1-st, vector mesons in static limit
c        -----------------------------------
c
c
c
c        vector coupling
c
c
c
c
  900 mc=1
c
c
c
c        central term  '1'  only
c
c
c
      ff=1.d0
      f(1)=2.d0
      f(2)=0.d0
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-f(1)
      f(8)=f(7)
      f(9)=f(2)
      f(10)=f(2)
      f(11)=f(2)
      xxpyyh=xxpyy*0.5d0
c
c
c**** goto 990
c
c
c        central term  ' k**2 + p**2 '
c
c
c
      f(1)=f(1)+xxpyyh
      f(2)=f(2)+xy2
      f(3)=f(3)+xxpyyh
      f(6)=f(6)+xxpyyh
      f(7)=f(7)-xxpyyh
      f(8)=f(8)-xxpyyh
      f(10)=f(10)+xy2
      f(11)=f(11)-xy2
c
c
c
c        spin - spin
c
c
c
      f(1)=f(1)+xxpyyh*3.d0
      f(2)=f(2)-xy*3.d0
      f(3)=f(3)-xxpyyh
      f(6)=f(6)-xxpyyh
      f(7)=f(7)+xxpyyh
      f(8)=f(8)+xxpyyh
      f(10)=f(10)+xy
      f(11)=f(11)-xy
c
c
c
c        tensor  with  spin-spin-contribution
c
c
c
      f(1)=f(1)-xxpyyh
      f(2)=f(2)+xy
      f(3)=f(3)+xxpyyh
      f(4)=f(4)-xy
      f(5)=f(5)+xy
      f(6)=f(6)-xxpyyh
      f(7)=f(7)+(xx-yy)*0.5d0
      f(8)=f(8)-(xx-yy)*0.5d0
c
c
c
c        spin - orbit
c
c
c
      xy3=xy*3.d0
      f(4)=f(4)+xy3
      f(5)=f(5)+xy3
      f(10)=f(10)-xy3
      f(11)=f(11)+xy3
c
  990 continue
c
      call obstr(2,1,me)
c
c**** goto 1995
c
c
c        case of additional sigma-l
c
c
c
c     mc=-1
c     xxyy8=xxyy/8.d0
c     f(1)=xxyy8
c     f(2)=0.d0
c     f(3)=f(1)
c     f(4)=f(2)
c     f(5)=f(2)
c     f(6)=-xxyy8
c     f(7)=f(1)
c     f(8)=f(7)
c     f(9)=f(6)*2.d0
c
c     call obstr(4,1,me)
c
c
c
c
c        tensor coupling
c
c      ( 1-ss , rho-meson in static limit, tensor coupling only )
c
c
 1000 mc=2
c
      if (mg.eq.10) mc=1
c
c
c        spin - spin
c
c
c
      f(1)=xxpyyh*3.d0
      f(2)=-xy*3.d0
      f(3)=-xxpyyh
      f(4)=0.d0
      f(5)=0.d0
      f(6)=-xxpyyh
      f(7)=xxpyyh
      f(8)=xxpyyh
      f(9)=0.d0
      f(10)=xy
      f(11)=-xy
c
c
c
c        tensor  with  spin-spin-contribution
c
c
c
      f(1)=f(1)-xxpyyh
      f(2)=f(2)+xy
      f(3)=f(3)+xxpyyh
      f(4)=f(4)-xy
      f(5)=f(5)+xy
      f(6)=f(6)-xxpyyh
      f(7)=f(7)+(xx-yy)*0.5d0
      f(8)=f(8)-(xx-yy)*0.5d0
c
      call obstr(2,1,me)
      if (mg.eq.10) go to 1995
c
c
c
c        case of additional sigma-l
c
c
c
c     f(1)=xxyy
c     f(2)=0.d0
c     f(3)=f(1)
c     f(4)=f(2)
c     f(5)=f(2)
c     f(6)=-xxyy
c     f(7)=f(1)
c     f(8)=f(7)
c     f(9)=f(6)*2.d0
c
c     call obstr(4,1,me)
c
c
c
c
c        vector-tensor coupling
c
c
c
c
      mc=3
c
c
c
c        central term  ' k**2 + p**2 '
c
c
c
      f(1)=-xxpyy
      f(2)=xy2
      f(3)=f(1)
      f(4)=0.d0
      f(5)=0.d0
      f(6)=f(1)
      f(7)=-f(1)
      f(8)=f(7)
      f(9)=0.d0
      f(10)=f(2)
      f(11)=-f(2)
c
c
c
c        spin - spin
c
c
c
      f(1)=f(1)+xxpyy*3.d0
      f(2)=f(2)-xy2*3.d0
      f(3)=f(3)-xxpyy
      f(6)=f(6)-xxpyy
      f(7)=f(7)+xxpyy
      f(8)=f(8)+xxpyy
      f(10)=f(10)+xy2
      f(11)=f(11)-xy2
c
c
c
c        tensor  with  spin-spin-contribution
c
c
c
      f(1)=f(1)-xxpyy
      f(2)=f(2)+xy2
      f(3)=f(3)+xxpyy
      f(4)=f(4)-xy2
      f(5)=f(5)+xy2
      f(6)=f(6)-xxpyy
      f(7)=f(7)+(xx-yy)
      f(8)=f(8)-(xx-yy)
c
c
c
c        spin - orbit
c
c
c
      xy4=xy*4.d0
      f(4)=f(4)+xy4
      f(5)=f(5)+xy4
      f(10)=f(10)-xy4
      f(11)=f(11)+xy4
c
      call obstr(2,1,me)
c
c
c
c        case of additional sigma-l
c
c
c
c     f(1)=xxyy
c     f(2)=0.d0
c     f(3)=f(1)
c     f(4)=f(2)
c     f(5)=f(2)
c     f(6)=-xxyy
c     f(7)=f(1)
c     f(8)=f(7)
c     f(9)=f(6)*2.d0
c
c     call obstr(4,1,me)
      go to 1995
c
c
c
c
c
c        1+  , a1-meson
c        --------------
c
c
c
c
 1100 mc=1
c
      ff=2.d0
      f(1)=eep1+ee
      f(2)=0.d0
      f(3)=-ee
      f(4)=-xy
      f(5)=xy2
      f(6)=-1.d0
      f(7)=ey
      f(8)=ex
c
      call obstr (1,1,me)
      go to 1995
c
c
c
c
c
c        2+  , f-meson
c        -------------
c
c
c
c
c        g1**2 coupling
c
c
c
c
 1200 mc=1
c
      ff=-8.d0
      ee4=ee*4.d0
      xxyy4=xxyy*4.d0
      eep143=eep1*fd3
        e1=2.d0*(xxpyy+td3)
      f(1)= eep143 +(3.d0*ee+2.d0)*xxpyy+xxyy4
      f(2)=(ee4+td3+xxpyy)*xy
      f(3)=3.d0*xxyy+eep1*e1
      f(4)=(2.d0*xxpyy+3.d0*eep1-d3)*xy
      f(5)=(ee4+11.d0*d3+3.d0*xxpyy)*xy
      f(6)= eep143 +(ee+2.d0)*xxpyy+xxyy4
        e2=-epe*e1
      f(7)=e2+xx*ex
      f(8)=e2+yy*ey
c        factors for additional terms
      f(9)=xxyy
      f(10)=ee*xy
      f(11)=xy
      f(12)=-ey*xy
      f(13)=-ex*xy
c
      call obstr(3,1,me)
      go to 1995
c
c
c
c
c        this has been the end of the contributions of mesons
c        ----------------------------------------------------
c
c
c
c
c        errors and warnings
c        -------------------
c
c
c
c
 9000 if (indmg(mg)) go to 1995
      write (kwrite,19000) mesong(mg)
19000 format(1h0////'0warning in obnn: meson-group  ',a4,'  does not exi
     1st in this program.'/'0contribution ignored. execution continued.'
     2////)
      indmg(mg)=.true.
c
c
c
c
 1995 continue
c
c
c
c
c        add up contributions of mesons
c        ------------------------------
c
c
c
c
 2000 do 2005 il=iman,imen
      do 2005 iv=1,6
 2005 v(iv)=v(iv)+vj(iv,il)
c
c
c
c
      return
      end
      subroutine obpar
        use crdwrt
        use cstate
        use cob
c
c        obpar reads writes and stores the parameter for all ob-subrout.
c
c
      implicit double precision (a-h,o-z)
c
c
!      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
!      common /cstate/ j,heform,sing,trip,coup,endep,label
!      logical heform,sing,trip,coup,endep
c
!      character*4 label
c
c        common block for all ob-subroutines
c
!      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
!     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
!     2                ez1,ez2,ct(96),wt(96),
!     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
!     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
!     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
!     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
!      logical indc,indxy,indpar
c
c
c        further specifications
c
      dimension cc(5)
      character*4 name(3),nname(15)
      integer imga(3)
      character*4 cut,cutg,end
      logical index
      character*4 mesong(12)
      logical zerocp,indcut
c
      data cut/'cut '/,cutg/'cutg'/,end/'end '/
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1                   '1-  ','1-t ','1-tt','1-st','1-ss',
     2                    '1+  ','2+  '/
      data index/.false./
      data zerocp/.true./,indcut/.false./
      data uf/197.3286d0/
c
c
c
c
c        statement functions needed for cray fortran
c
c
c      sqrt(x)=sqrt(x)
c      exp(x)=exp(x)
c      datan(x)=atan(x)
c      log(x)=log(x)
c
c
c
c
10000 format (2a4,a2,15a4)
10001 format (1h1)
10002 format (1h0/' jp  name      g**2      f/g       mass    iso-spin
     1iprop/spe'/9x,'cut typ     c u t - o f f   p a r a m e t e r s')
10003 format (2a4,a2,5f10.4)
10004 format (1h0,2a4,a2,2f10.4,f9.2,1x,2(f7.1,3x))
10005 format (1h ,2a4,a2,f3.1,f11.1,f9.4,f14.4,f10.4)
10006 format (2a4,a2,3i3)
10007 format (1h0,2a4,a2,3i3)
10008 format (1h ,61(1h-))
10011 format (1h1  //' obnn:  one-boson-exchange nn-nn interaction (nume
     1r. integra.)')
10012 format (1h1  //' obnd:  one-boson-exchange nn-nd interaction (nume
     1r. integra.)')
10013 format (1h1  //' obdd:  one-boson-exchange nn-dd interaction (nume
     1r. integra.)')
10015 format ('0input-parameter-set:'/1h ,20(1h-))
10016 format (1h0,2a4,a2,15a4)
c
c
10020 format (2a4,a2,f10.7,4f10.4)
10021 format (1h0,2a4,a2,f10.7,f10.4,f9.2,1x,2(f7.1,3x))
c
c
      if (index) go to 50
      index=.true.
c
      x=-1.d0
      y=-1.d0
c
c
c
c
c        maxima of certain indices related to the dimension as follows:
c        dimension c(mme,imee),ic(mice,imee),indc(mindce,imee),
c                  mgg(mge,3),mggo(mge,3),mesong(mge),vj(32,imee),
c                  ima(mee,mge,3)
c
      mge=12
      mee=5
      mme=10
      mice=10
      mindce=2
      imb=1
      ime=0
      imee=15
      imec=0
c        mme always ge mice, mindce
c
c        set all meson-parameters and indices to zero or .false.
c
      do 1 int=1,3
      imga(int)=0
      indpar(int)=.false.
      do 1 mgx=1,mge
      mgg(mgx,int)=0
    1 mggo(mgx,int)=0
c
c
      do 2 il=1,imee
      do 2 mm=1,mme
      if (mm.le.mindce) indc(mm,il)=.false.
      if (mm.le.mice) ic(mm,il)=0
    2 c(mm,il)=0.d0
      endep=.false.
c
c
c        indc(2,il) is reserved for information concerning the eikonal
c        form-factor
c
c
c
c
c
c
c        reading and writing of first 4,5 cards
c        --------------------------------------
c        --------------------------------------
c
c
c
c        write headline and read and write name of parameter set
c
   50 go to (51,52,53),inter
   51 write (kwrite,10011)
      go to 55
   52 write (kwrite,10012)
      go to 55
   53 write (kwrite,10013)
   55 write (kwrite,10008)
      write (kwrite,10015)
      read  (kread, 10000) name,nname
      write (kwrite,10016) name,nname
      if (inter.eq.1) label=name(1)
      indpar(inter)=.true.
c
c        read and write index-parameter concerning the factor of the
c        potential
c
      read  (kread, 10006) name,ift(inter)
      write (kwrite,10007) name,ift(inter)
      iftyp=ift(inter)
c**** if (iftyp.lt.0.or.iftyp.gt.4.or.inter.eq.1.and.iftyp.gt.2)goto9003
c
c        read and write parameters for numerical integration
c
      read  (kread, 10006) name,mint(inter),maxt(inter)
      write (kwrite,10007) name,mint(inter),maxt(inter)
c
c        read and write mass of nucleon
c
      read  (kread, 10003) name,wn
      write (kwrite,10004) name,wn
      wnq=wn*wn
      dwn=1.d0/wn
      dwnq=dwn*dwn
      wnn(inter)=wn
      if (inter.lt.2) go to 60
c
c        read and write mass of n-star
c
      read  (kread, 10003) name,wdd(inter)
      write (kwrite,10004) name,wdd(inter)
c
c        write headline for meson parameters
c
   60 write (kwrite,10002)
      write (kwrite,10008)
c
c
c
c
c        read, write and store meson parameters
c        --------------------------------------
c        --------------------------------------
c
c
c
   61 read  (kread, 10020) name,cc
c
c        check if data-card just read contains cut-off parameters
c
      if (name(1).eq.cut.or.name(1).eq.cutg) go to 70
c
c        check if end of mesons
c
      if (name(1).eq.end) go to 2000
c
c
c
c
c        write meson-parameters, which are no cut-off parameters
c        -------------------------------------------------------
c
c
c
c
      indcut=.false.
c
      write (kwrite,10021) name,cc
c
c        check if coupling constants are zero
c
      if (cc(1).ne.0.d0) go to 62
      zerocp=.true.
      go to 61
c
   62 zerocp=.false.
c
c        find out number of meson-group mg
c
      do 63 mg=1,mge
      if (name(1).eq.mesong(mg)) go to 64
   63 continue
      go to 9000
c
c
c
c
c        store meson parameters, which are no cut-off parameters
c        -------------------------------------------------------
c
c
c
c
   64 ime=ime+1
      if (ime.gt.imee) go to 9011
      mgg(mg,inter)=mgg(mg,inter)+1
      m=mgg(mg,inter)
      if (m.gt.mee) go to 9001
      ima(m,mg,inter)=ime
      if (m.ne.1) go to 65
      imga(inter)=imga(inter)+1
      mggo(imga(inter),inter)=mg
   65 continue
c
c        store coupling constant g**2
      c(1,ime)=cc(1)
c        store coupling constant f*g
      c(3,ime)=cc(1)*cc(2)
      if (inter.eq.2.and.mg.eq.10) c(1,ime)=c(1,ime)+c(3,ime)
c        store coupling constant f**2
      c(2,ime)=cc(2)*c(3,ime)
      if (inter.eq.1.and.mg.eq.10)
     1  c(1,ime)=c(1,ime)+c(3,ime)*2.d0+c(2,ime)
c        store meson mass sqare in units of nucleon mass square
      c(4,ime)=cc(3)*cc(3)*dwnq
c
c        test iso-spin
      icc=cc(4)
      if (icc.ne.0.and.icc.ne.1) go to 9004
c         store isospin as logical constant
      if (icc.eq.1) indc(1,ime)=.true.
c        store and test iprsp for meson/delta/nucleon
      icc=cc(5)
      iccc=mod(icc,100)
c        iprsp for nucleons
      ic(1,ime)=mod(iccc,10)
c        ispd for deltas
      ic(2,ime)=iabs(iccc/10)
c        ispm for mesons
      ic(3,ime)=iabs(icc/100)
      if (iabs(ic(1,ime)).gt.8) go to 9005
      if (iabs(ic(1,ime)).ge.2.and.iabs(ic(1,ime)).le.7) endep=.true.
c
c        index values for further storing
      mi=4
      mm=5
      go to 61
c
c
c
c
c        write cut-off parameters
c        ------------------------
c
c
c
c
   70 write (kwrite,10005) name,cc
c
c
c        check if individuel cut or general cut
c
      if (name(1).eq.cut) go to 73
c        case of general cut-off
      if (indcut) go to 90
      if (imec.ge.ime) go to 61
      imac=imec+1
      imec=ime
      if (imac.lt.imb) imac=imb
      go to 90
c        case of individuel cut-off
   73 imac=ime
      imec=ime
      if (zerocp) go to 61
c
c        save present values of indices
c
   90 indcut=.true.
      if (cc(1).eq.0.d0) go to 61
      mix=mi
      mmx=mm
c
c        start loop of mesons, which present cut-off refers to
c
      do 1095 im=imac,imec
      mi=mix
      mm=mmx
c
c
c
c
c        store cut-off parameters
c        ------------------------
c
c
c
c
c        store typ of cut-off
      ic(mi,im)=cc(1)
      ityp=ic(mi,im)
      if (ityp.lt.1.or.ityp.gt.9) go to 9002
c        store and test typ of propagator of cut-off
      ic(mi+1,im)=cc(2)
      if (ic(mi+1,im).lt.0.or.ic(mi+1,im).gt.8) go to 9006
      if (ic(mi+1,im).ge.2.and.ic(mi+1,im).le.7) endep=.true.
      go to (100,100,300,400,400,400,700,800,900),ityp
c
c
c        cut-off of dipole type
c        **********************
c
c
c        store and test exponent of cut-off
  100 ic(mi+2,im)=cc(3)
      if (ic(mi+2,im).lt.0) go to 9009
      if (ic(mi+2,im).gt.0) go to 101
c        exponent is zero, omit cut-off
      ic(mi,im)=0
      ic(mi+1,im)=0
      go to 1000
c        store cut-off mass for denominator
  101 c(mm+1,im)=cc(4)*cc(4)*dwnq
c        store numerator of cut-off
      c(mm,im)=c(mm+1,im)
      if (ityp.eq.2)     c(mm,im)=c(mm,im)-c(4,im)
      mi=mi+3
      mm=mm+2
      go to 1000
c
c
c        cut-off of regge type /schierholz/
c        **********************************
c
c
  300 c(mm,im)=2.d0/sqrt(4.d0-cc(3)*cc(3)*dwnq)
      c(mm+1,im)=cc(4)-1.d0
      c(mm+2,im)=cc(5)*wnq*1.d-6
      mi=mi+2
      mm=mm+3
      go to 1000
c
c
c        eikonal form factor
c        *******************
c
c
c        store gamma as -4.*gamma
  400 c(mm,im)=-cc(3)*4.d0
      ieik=ityp-3
c
c        compute and store normalization factor of t-form
      d=-c(4,im)
      d1=sqrt(-d)
      d2=sqrt(4.d0+d)
c     c(mm+1,im)=exp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*datan(d1/d2))
      c(mm+1,im)=exp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*atan2(d1,d2))
      print*,atan(d1/d2),atan2(d1,d2)
c
      go to (490,412,412),ieik
c
c        compute and store normalization factor of u-form
c        and t- * u-form
  412 ieik2=cc(4)
      if (ieik2.eq.0) ieik2=1
      ic(mi+2,im)=ieik2
      go to (421,422,423,421,422,423),ieik2
c
  421 ess=4.d0*(1.d0+cc(5)*dwn)**2
      go to 450
c
  422 endep=.true.
c
  423 go to 490
c
  450 d=-(4.d0-ess)
      if (ieik2.le.3) d=d+c(4,im)
      if (d.eq.0.d0) go to 460
      d2=sqrt(4.d0+d)
      if (d.lt.0.d0) go to 470
      d1=sqrt(d)
      c(mm+2,im)=exp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*log(0.5d0*(d1+d2)))
      go to 480
c
  460 c(mm+2,im)=exp(2.d0*cc(3))
      go to 480
c
  470 d1=sqrt(-d)
c     c(mm+2,im)=exp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*datan(d1/d2))
      c(mm+2,im)=exp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*atan2(d1,d2))
      print*,atan(d1/d2),atan2(d1,d2)
c
  480 if (ieik.eq.3) c(mm+2,im)=c(mm+2,im)*c(mm+1,im)
c
  490 mi=mi+3
      mm=mm+3
      go to 1000
c
c
c        exponential form factor
c        ***********************
c
c
c        check exponent
  700 if (cc(3).lt.0.d0) go to 9009
      if (cc(3).gt.0.d0) go to 701
c        exponent is zero, omit cutoff
      ic (mi,im)=0
      ic (mi+1,im)=0
      go to 1000
c        compute constant factor for argument of exponential function
  701 c(mm,im)=cc(3)*wnq/(cc(4)*cc(4))
      mi=mi+2
      mm=mm+1
      go to 1000
c
c
c        cloudy bag form factor
c        ***********************
c
c
c        check exponent
  800 if (cc(3).lt.0.d0) go to 9009
      if (cc(3).gt.0.d0) go to 801
c        exponent is zero, omit cutoff
      ic (mi,im)=0
      ic (mi+1,im)=0
      go to 1000
c        store exponent
  801 ic(mi+2,im)=cc(3)
c        store cutoff radius
      c(mm,im)=cc(4)*wn/uf
      mi=mi+3
      mm=mm+1
      go to 1000
c
c
c        propagator of mass-distributed meson
c        ************************************
c
c
c        spin of meson
  900 icc=cc(3)
      ic(mi+2,im)=icc
c        full width
      c(mm,im)=cc(4)*dwn
c        2*pion-mass squared
      c(mm+1,im)=cc(5)*cc(5)*dwnq
c        recalculate width
      d=c(4,im)-c(mm+1,im)
      c(mm,im)=c(mm,im)*sqrt(c(4,im))/sqrt(d)
      if (icc.lt.1) go to 910
      do 905 i=1,icc
  905 c(mm,im)=c(mm,im)/d
  910 mi=mi+3
      mm=mm+2
      go to 1000
c
c
c
c
c        end cut-offs
c        ************
c
c        test dimensions
 1000 if (mi.gt.mice.or.mm-1.gt.mme) go to 9010
c
c
 1095 continue
      go to 61
c
c
c
c
c        last card
c        ---------
c        ---------
c
c
c
c
c        write end mesons
 2000 imaa(inter)=imb
      imea(inter)=ime
      imb=ime+1
      write (kwrite,10004) name
      write (kwrite,10008)
      write (kwrite,10008)
c
c
c
c
      return
c
c
c
c        errors
c        ------
c        ------
c
c
c
c
 9000 write (kwrite,19000) name(1)
19000 format (1h0/////'0error in obpar:  meson-group   ',a4,'   does not
     1 exist in this program.'/'0execution terminated.'////)
      go to 9999
c
c
 9001 write (kwrite,19001)
19001 format (1h0/////'0error in obpar: too many mesons within a meson-g
     1roup with respect to'/'0the given dimensions. execution terminatee
     2.'////)
      go to 9999
c
c
 9002 write (kwrite,19002) cc(1)
19002 format (1h0/////'0error in obpar: cut-off typ',f10.4,'  does not e
     1xist in this program.'/'0execution terminated.'////)
      go to 9999
c
c
 9003 write (kwrite,19003) iftyp
19003 format (1h0/////'0error in obpar: factor typ has the non-permissib
     1le value',i4,' .'/'0execution terminated.'////)
      go to 9999
c
c
 9004 write (kwrite,19004) cc(4)
19004 format (1h0/////'0error in obpar: isospin has the non-permissible
     1value',f10.4,'  .'/'0execution terminated.'////)
      go to 9999
c
c
 9005 write (kwrite,19005) cc(5)
19005 format (1h0/////'0error in obpar: iprop/spe has the non-permissibl
     1e value',f10.4,'  .'/'0execution terminated.'////)
      go to 9999
c
c
 9006 write (kwrite,19006) cc(2)
19006 format (1h0/////'0error in obpar: the index for the propagator of
     1the cut-off has the'/'0non-permissible value',f10.4,'  . execution
     2 terminated.'////)
      go to 9999
c
c
 9009 write (kwrite,19009)
19009 format (1h0/////'0error in obpar: the exponent of the cut-off is l
     1ess than zero.'/'0execution terminated.'////)
      go to 9999
c
c
 9010 write (kwrite,19010)
19010 format (1h0/////'0error in obpar: too many cut-off parameters with
     1 respect to the given'/'0dimensions. execution terminated.'////)
      go to 9999
c
c
 9011 write (kwrite,19011)
19011 format (1h0/////'0error in obpar:  too many mesons with respect to
     1 the dimensions given'/'0to this program. execution terminated.'
     2////)
      go to 9999
c
c
 9999 stop
      end
      subroutine obstr (icase,max,mex)
        use crdwrt
        use cstate
        use cob
c
c        obstr computes the structure of one-boson-exchanges
c
c
      implicit double precision (a-h,o-z)
c
c
c        common blocks
c
!      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
!      common /cstate/ j,heform,sing,trip,coup,endep,label
!      logical heform,sing,trip,coup,endep
c
!      character*4 label
c
c        common block for all ob-subroutines
c
!      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
!     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
!     2                ez1,ez2,ct(96),wt(96),
!     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
!     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
!     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
!     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
!      logical indc,indxy,indpar
      logical index
c
c     further specifications
c
      dimension vv(32)
      dimension tt(2,3)
      logical indiso
      data jj/-1/
      data index/.false./
c
c
c
c
      if (index) go to 50
      index=.true.
c
c
      tt(1,1)=1.d0
      tt(2,1)=-3.d0
c
      do 1 ii=2,3
      do 1 i=1,2
    1 tt(i,ii)=1.d0
c
c
c
c
c
   50 do 1095 m=max,mex
      im=ima(m,mg,inter)
c
c
      if (mc.ne.1) go to 60
c
c
c
c
c        call integrals
c        --------------
c
c
c
c
      call obai
c
c
c
c
   60 if (mc.lt.1) mc=1
c
      if (c(mc,im).eq.0.d0) go to 1095
c
c
c
c
      go to (100,200,300),inter
c
c
c
c
c        nn-nn helicity amplitudes /combinations/
c        ----------------------------------------
c
c
c
c
c        ground structure (a factor of 2 is included in v5 and v6)
c
c
  100 ive=6
c
      vv(1)=f(1)*ai(1,m)+f(2)*ai(2,m)
      vv(2)=f(3)*ai(1,m)+f(4)*ai(3,m)
      vv(3)=f(5)*ai(1,m)+f(6)*ai(2,m)
      vv(4)=f(4)*ai(1,m)+f(3)*ai(3,m)
      vv(5)=f(7)*ai(4,m)
      vv(6)=f(8)*ai(4,m)
c
c
      write(0,*)'icase=',icase
      go to (1000,120,130,140),icase
c
c
c        additional terms  in case of tensor coupling
c
c
  120 vv(1)=vv(1)+f(9)*ai(5,m)
      vv(2)=vv(2)+f(10)*ai(2,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(10)*ai(5,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(10)*ai(6,m)
         e1=f(11)*ai(7,m)
      vv(5)=vv(5)+e1
      vv(6)=vv(6)+e1
      go to 1000
c
c
c        additional terms in case of 2+ mesons
c
c
  130 vv(2)=vv(2)+f(10)*ai(2,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(11)*ai(5,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(10)*ai(6,m)
      vv(5)=vv(5)+f(12)*ai(7,m)
      vv(6)=vv(6)+f(13)*ai(7,m)
      go to 1000
c
c
c        additional terms in case of sigma-l in static limit
c
c
  140 vv(1)=vv(1)+f(6)*ai(5,m)
      vv(2)=vv(2)+f(1)*ai(5,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(1)*ai(11,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(1)*ai(12,m)
      vv(5)=vv(5)+f(6)*ai(13,m)
      vv(6)=vv(6)+f(6)*ai(13,m)
      go to 1000
c
c
c
c
c        nn-nd helicity amplitudes
c        -------------------------
c
c
c
c
  200 ive=16
c
      ai3m1=ai(3,m)-ai(1,m)
      ai6m2=ai(6,m)-ai(2,m)
      ai3p1=ai(3,m)+ai(1,m)
      ai6p2=ai(6,m)+ai(2,m)
c
c
      vv( 1)= f( 1)* ai( 4,m) + f( 2)* ai( 7,m)
      vv( 2)= f( 3)* ai( 1,m) + f( 4)* ai( 2,m) + f( 5)* ai( 5,m)
      vv( 3)= f( 6)* ai( 4,m) + f( 7)* ai( 7,m)
      vv( 4)= f( 8)* ai( 8,m)
      vv( 5)= f( 9)* ai( 8,m)
      vv( 6)=-f(10)* ai( 4,m) + f(11)* ai( 7,m)
      vv( 7)= f(12)* ai( 1,m) - f(13)* ai( 2,m) + f(14)* ai( 5,m)
      vv( 8)=-f(15)* ai( 4,m) + f(16)* ai( 7,m)
c
      vv( 9)= f(17)* ai3m1    + f(18)* ai6m2
      vv(10)= f(19)* ai( 4,m) + f(20)* ai( 7,m)
      vv(11)= f(21)* ai3p1    + f(22)* ai6p2
      vv(12)= f(23)* ai( 9,m)
      vv(13)=-f(24)* ai(10,m)
      vv(14)=-f(25)* ai3m1    + f(26)* ai6m2
      vv(15)=-f(27)* ai( 4,m) + f(28)* ai( 7,m)
      vv(16)=-f(29)* ai3p1    + f(30)* ai6p2
      go to 1000
c
c
c
c        nn-dd helicity amplitudes
c        -------------------------
c
c
c
c
  300 ive=32
c
      ai31p=ai( 3,m)+ai( 1,m)
      ai31m=ai( 3,m)-ai( 1,m)
      ai62p=ai( 6,m)+ai( 2,m)
      ai62m=ai( 6,m)-ai( 2,m)
      aic5p=ai(12,m)+ai( 5,m)
      aic5m=ai(12,m)-ai( 5,m)
c
c
      vv( 1)= f( 1)*(ai(1,m)+ai(2,m)) + f( 2)*(ai(2,m)+ai(5,m)) +
     1        f( 3)*(ai(5,m)+ai(11,m))
c
      vv( 2)= f( 4)* ai( 4,m) + f( 5)* ai( 7,m) + f( 6)* ai(13,m)
      vv( 3)= f( 7)* ai( 8,m) + f( 8)* ai(14,m)
      vv( 4)= f( 9)* ai(17,m)
      vv( 5)=vv( 2)
      vv( 6)= f(10)* ai( 1,m) + f(11)* ai( 2,m) + f(12)* ai( 5,m) +
     1        f(13)* ai(11,m)
      vv( 7)= f(14)* ai( 4,m) + f(15)* ai( 7,m) + f(16)* ai(13,m)
      vv( 8)=-f(17)* ai( 8,m) + f(18)* ai(14,m)
      vv( 9)=vv( 3)
      vv(10)=vv( 7)
      vv(11)=-f(19)* ai( 1,m) + f(20)* ai( 2,m) - f(21)* ai( 5,m) +
     1        f(22)* ai(11,m)
      vv(12)= f(23)* ai( 4,m) - f(24)* ai( 7,m) + f(25)* ai(13,m)
      vv(13)=vv( 4)
      vv(14)=vv( 8)
      vv(15)=vv(12)
c
      vv(16)=-f(26)*(ai(1,m)-ai(2,m)) + f(27)*(ai(2,m)-ai(5,m)) -
     1        f(28)*(ai(5,m)-ai(11,m))
c
c
      vv(17)= f(29)* ai( 4,m) + f(30)* ai( 7,m) + f(31)* ai(13,m)
      vv(18)= f(32)* ai31p    + f(33)* ai62p    + f(34)* aic5p
      vv(19)= f(35)* ai( 9,m) + f(36)* ai(15,m)
      vv(20)= f(37)* ai(18,m)
      vv(21)= f(44)* ai31m    - f(45)* ai62m    + f(46)* aic5m
      vv(22)= f(38)* ai( 4,m) + f(39)* ai( 7,m) + f(40)* ai(13,m)
      vv(23)= f(41)* ai31p    + f(42)* ai62p    + f(43)* aic5p
      vv(24)=vv(19)
      vv(25)= f(47)* ai(10,m) - f(48)* ai(16,m)
      vv(26)= f(49)* ai31m    - f(50)* ai62m    + f(51)* aic5m
      vv(27)=vv(22)
      vv(28)=vv(18)
      vv(29)= f(52)* ai(19,m)
      vv(30)=vv(25)
      vv(31)=vv(21)
      vv(32)=vv(17)
c
c
c
c
 1000 if (inter.ge.2) go to 1040
c
c
c
c
c        set certain cases to zero in case of inter=1
c
      if (j.ne.0) go to 1021
      vv(2)=0.d0
      vv(4)=0.d0
      vv(5)=0.d0
      vv(6)=0.d0
c
 1021 if (.not.sing) vv(1)=0.d0
      if (.not.trip) vv(2)=0.d0
      if (coup) go to 1030
      do 1025 iv=3,6
 1025 vv(iv)=0.d0
c
 1030 if (heform) go to 1040
c
c
c        transformation into lsj-formalism in case of inter=1
c        (if requested)
      if (j.eq.jj) go to 1035
      jj=j
      aj=float(j)
      aj1=float(j+1)
      d2j1=1.d0/float(2*j+1)
      arjj1=sqrt(aj*aj1)
c
 1035 v3=vv(3)
      v4=vv(4)
      v5=vv(5)
      v6=vv(6)
      v34=-arjj1*(v3-v4)
      v56=arjj1*(v5+v6)
      vv(3)=d2j1*(aj1*v3+aj*v4-v56)
      vv(4)=d2j1*(aj*v3+aj1*v4+v56)
      vv(5)=d2j1*(v34-aj1*v5+aj*v6)
      vv(6)=d2j1*(v34+aj*v5-aj1*v6)
c
c
c        possible different sign depending on the convention used
      vv(5)=-vv(5)
      vv(6)=-vv(6)
c
c
c
c
c        multiply with factors
c        ---------------------
c
c
c
c
 1040 is=mod(j,2)+1
      it=mod(is,2)+1
      indiso=indc(1,im)
      cmc=c(mc,im)
      fc=fff*ff*cmc
      do 1045 iv=1,ive
c
c        multiply with coupling-constant and factors fff and ff
c
      vv(iv)=vv(iv)*fc
      if (inter.ge.2) go to 1045
c
c        multiply with isospin factor
c
      if (.not.indiso) go to 1045
      if (iv.eq.2) go to 1043
      vv(iv)=vv(iv)*tt(is,inter)
      go to 1045
 1043 vv(iv)=vv(iv)*tt(it,inter)
c
c     add up in case of several couplings for one meson-exchange
c     and store
 1045 vj(iv,im)=vj(iv,im)+vv(iv)
c
c
c        if single contributions to one meson are to be printed:
c     write (kwrite,10001) mesong(mg),mc
c     write (kwrite,10002) (vv(iv),iv=1,ive)
10001 format (' contribution from ',a4,' mc =',i2)
10002 format (33x,6d16.6)
c
c
 1095 continue
c
c
      return
      end
      subroutine obai
        use cpot
        use cstate
        use cpoted
        use cob
c
c        obai integrates over theta
c
c
      implicit double precision (a-h,o-z)
c
!      common /cpot/   v(6),xmev,ymev
!      common /cstate/ j,heform,sing,trip,coup,endep,label
!      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,
!     1                znrl,zrel,smev,noced
!      logical heform,sing,trip,coup,endep
!      logical noced
c
!      character*4 label
c
c        common block for all ob-subroutines
c
!      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
!     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
!     2                ez1,ez2,ct(96),wt(96),
!     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
!     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
!     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
!     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
!      logical indc,indxy,indpar
c
c
c        further specifications
      dimension gi(7)
c
      dimension pj(7,96),pjn(96),pjm1(96)
c      real*4 aez,axy2,aomq,am1,am2,am
      logical indj,indint,indee,indepe,indz,indepz
      data ige/7/
      data nnt/-1/,iinter/-1/,jj/-1/
c
c
c        statement functions needed for cray fortran
c
c      float(i)=float(i)
c      sqrt(x)=sqrt(x)
c
c
c
c
      if (inter.eq.iinter) go to 60
      iinter=inter
      indint=.false.
      mino=mint(inter)
      maxo=maxt(inter)
c
      go to (51,51,53),inter
   51 igeint=5
      go to 55
   53 igeint=7
   55 continue
c
      wn=wnn(inter)
      dwn=1.d0/wn
      wnq=wn*wn
c        wd is the mass of the delta
c        wdd(..) is the array of the masses of the delta
      wd=wdd(inter)
      wdq=wd*wd
c
c
c
c
   60 if (j.eq.jj) go to 70
      jj=j
      indj=.false.
c
c
      aj=float(j)
      aj1=float(j+1)
      dj1=1.d0/aj1
      ajdj1=aj*dj1
      aaj=sqrt(ajdj1)
c
c
      aj2=float(j+2)
      aj3=float(j+3)
      ajm1=float(j-1)
      ajm2=float(j-2)
      ajm3=float(j-3)
c
c
      ajj1=aj*aj1
      ajj2=ajm1*aj2
      ajj3=ajm2*aj3
      ajja=aj*ajm3
      ajjb=aj*ajm1
c
      aajj=0.d0
      if (j.gt.1)
     1aajj=aj/sqrt(ajj1*ajj2)
c
      aaj1=aajj*ajm1
      aaj2=aajj*aj1
      aaj3=aajj*2.d0
c
      if (j.gt.1) go to 62
      aajj=0.d0
      go to 63
   62 aajj=1.d0/(aj1*sqrt(ajj2))
c
   63 aaj4=aajj*ajjb
      aaj5=aajj*aj1*2.d0
      aaj6=aajj*(ajj1+2.d0)
      aaj7=aajj*ajj2
c
      if (j.gt.2) go to 64
      aajj=0.d0
      go to 65
   64 aajj=-aj/sqrt(ajj1*ajj2*ajj3)
c
   65 aaj8=aajj*(ajj1+6.d0)
      aaj9=aajj*ajj2
      aaj10=aajj*(ajja+2.d0)
      aaj11=aajj*(ajja-6.d0)
c
      if (j.gt.2) go to 66
      aajj=0.d0
      go to 67
   66 aajj=-1.d0/(aj1*sqrt(ajj2*ajj3))
c
   67 aaj12=aajj*ajjb*ajm2
      aaj13=aajj*(aj*ajjb+4.d0*aj+12.d0)
      aaj14=aajj*(5.d0*ajj1+6.d0)
      aaj15=aajj*3.d0*ajj2
      aaj16=aajj*ajj3*ajm1
      aaj17=aajj*aj1*ajj3
      aaj18=aajj*2.d0*ajj3
c
c
c
c
c
c        find out appropriate number of gauss-points
c        -------------------------------------------
c
c
c
c
   70 c4=c(4,im)
      iprsp=ic(1,im)
c
c
c        prepare starting energy
      if (noced) go to 73
      noced=.true.
      indz=.false.
      indepz=.false.
      if (iprsp.lt.2) go to 74
      if (iprsp.ge.4) go to 72
   71 indz=.true.
      z=2.d0*sqrt(wnq+q0qmev)
      go to 74
   72 if (iprsp.eq.8) go to 74
      indepz=.true.
      eppq=pmev*pmev
      epz=zrel
      go to 74
c
c
   73 if (indxy.and.indint) go to 85
   74 indint=.true.
      indee=.false.
      indepe=.false.
c
      ispm=ic(3,im)
      ez1=0.d0
      if (ispm.eq.1.and.iprsp.lt.2) go to 90
c
c         delta-nucleon mass difference in propagator
c
      ez1=(wd-wn)*dwn
c
c
      if (iprsp.lt.2) go to 90
c
c        prepare for propagators of non-covariant perturbation theory
   75 if (iprsp.ge.4) go to 76
      if (.not.indz) go to 71
      pmevq=0.d0
      ez=z
      go to 77
   76 if (iprsp.eq.8) go to 1800
      if (.not.indepz) go to 72
      pmevq=eppq
      ez=epz
   77 qfmevq=qfmev*qfmev
      qqx=xmev*xmev+pmevq
      ymevq=ymev*ymev
      qqy=ymevq+pmevq
      go to (1100,1200,1300),inter
c
c        ez1 for the nn case
c
 1100 ez1=spe(qqx,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,iprsp,0)
     1  + spe(qqy,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,iprsp,0)
     2  - ez
      ez1=ez1*dwn
      go to 80
c
c        ez1 for the nd case
c
 1200 pq=4.d0*pmevq/(wn+wd)**2
      ispd=ic(2,im)
      ed=spe(ymevq+pq*wdq,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wd,2,
     1 ispd,0)
      en=spe(ymevq+pq*wnq,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,
     1 iprsp,0)
      eza=spe(qqx,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,iprsp,0)
     1 - ez
      ez1=eza+en
      ez2=eza+ed
      ez1=ez1*dwn
      ez2=ez2*dwn
      go to 80
c
c        ez1 for the dd case
c
 1300 ispdd=ic(2,im)
      ez1=spe(qqx,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,iprsp,0)
     1  + spe(qqy,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wd,2,ispdd,0)
     2  - ez
      ez1=ez1*dwn
      go to 80
c
c
c        case iprsp=8
c
 1800 ez1=ex-ey
c
c
c        store ez1 and ez2
c
c
   80 if (iprsp.ge.3) go to 81
      indee=.true.
      ee1=ez1
      ee2=ez2
      go to 89
   81 iiprsp=iprsp
      indepe=.true.
      epe1=ez1
      epe2=ez2
      go to 89
c
c
c        get stored ez1 and ez2
c
c
   85 if (iprsp.lt.2) go to 90
      if (iprsp.ge.3) go to 86
      if (.not.indee) go to 75
      ez1=ee1
      ez2=ee2
      go to 89
   86 if (iprsp.ne.iiprsp) go to 75
      if (.not.indepe) go to 75
      ez1=epe1
      ez2=epe2
c
c
   89 aez=ez1
c
c
c        compute am
c
c
   90 axy2=xy2
      if (iprsp.ne.1) go to 91
      aomq=eem12+c4
      go to 92
   91 aomq=xxpyy+c4
c
   92 am=axy2/aomq
c
      if (iprsp.lt.2) go to 93
      am1=am
      am2=axy2/(aomq+aez*abs(aez))
      if (am2.lt.0.) go to 94
      am=max(am1,am2)
c
c
c        compute number of gausspoints (nt)
c
c
   93 if (am.gt.0.999) go to 94
c
c
      if (am.gt.0.85) am=am**(-log(1.-am)-0.9)
c
c
      nt=float(mino)/(1.-am)+0.9
c
c
      if (nt.gt.maxo) nt=maxo
      go to 95
c
c
   94 nt=maxo
c
c
   95 nt=nt+j
c
c        compute nt, which is suitable for gset
c
      if (nt.le.16) go to 98
      if (nt.gt.24) go to 96
      nt=4*(nt/4)
      go to 98
   96 if (nt.gt.48) go to 97
      nt=8*(nt/8)
      go to 98
   97 nt=16*(nt/16)
      if (nt.gt.96) nt=96
c
   98 if (nt.eq.nnt.and.indj) go to 100
c
c
c
c
c        call gauss-points
c        -----------------
c
c
c
c
      call gset (-1.d0,1.d0,nt,ct,wt)
      nnt=nt
c
c
c
c
c        call legendre-polynoms if necessary
c        -----------------------------------
c
c
c
c
      indxy=.false.
      indj=.true.
c
      call legpc (pjn,pjm1,ct,j,nt)
c
      write(0,*)'nt=',nt
      do 99 i=1,nt
      pj(1,i)=pjn(i)
      pj(3,i)=pjm1(i)
      pj(2,i)=pj(1,i)*ct(i)
      pj(4,i)=pj(2,i)*ct(i)
      pj(6,i)=pj(4,i)*ct(i)
      pj(5,i)=pj(3,i)*ct(i)
   99 pj(7,i)=pj(5,i)*ct(i)
c
c
c
c
c        call integrand
c        --------------
c
c
c
c
  100 call obaa
c
c
c
c
c        prepare for integration
c
c
c
c
      write(0,*)'igeint=',igeint
      do 2001 ig=1,igeint
 2001 gi(ig)=0.d0
c
c
c
c
c        integration-loop of theta
c        -------------------------
c
c
c
c
      do 2005 i=1,nt
      do 2005 ig=1,igeint
 2005 gi(ig)=gi(ig)+pj(ig,i)*aa(i)
c
c
c
      if (j.ne.0) go to 2010
      gi(3)=0.d0
      gi(5)=0.d0
      gi(7)=0.d0
c
c
c
c
c        combinations of integrals
c        -------------------------
c
c
c
c
 2010 ai(1,m)=gi(1)
c
      ai(2,m)=gi(2)
      ai(3,m)= ajdj1*gi(2)+dj1*gi(3)
      gi23m  =gi(2)-gi(3)
      ai(4,m)=aaj*gi23m
c
c
      ai(5,m)=gi(4)
      ai(6,m)= ajdj1*gi(4)+dj1*gi(5)
      gi45m  =gi(4)-gi(5)
      ai(7,m)=aaj*gi45m
c
c
      if (inter.eq.1) go to 3000
c
c
      ai( 8,m)= aaj1*gi(4)-aaj2*gi(1)+aaj3*gi(5)
      aai1    = aaj4*gi(4)+aaj5*gi(1)-aaj6*gi(5)
      aai2    = aaj7*gi23m
      ai( 9,m)= aai2+aai1
      ai(10,m)= aai2-aai1
c
c
      if (inter.ne.3) go to 3000
c
c
      ai(11,m)=gi(6)
      ai(12,m)=ajdj1*gi(6)+dj1*gi(7)
      ai(13,m)=aaj*(gi(6)-gi(7))
c
      ai(14,m)= aaj1*gi(6)-aaj2*gi(2)+aaj3*gi(7)
      aai1    = aaj4*gi(6)+aaj5*gi(2)-aaj6*gi(7)
      aai2    = aaj7*gi45m
      ai(15,m)= aai2+aai1
      ai(16,m)= aai2-aai1
c
      ai(17,m)= aaj8*gi(7)-aaj9*gi(3)-aaj10*gi(6)+aaj11*gi(2)
      aai1    =-aaj12*gi(6)+aaj13*gi(2)-aaj14*gi(7)+aaj15*gi(3)
      aai2    = aaj16*gi(4)-aaj17*gi(1)+aaj18*gi(5)
      ai(18,m)= aai1-aai2
      ai(19,m)= aai1+aai2
c
c
c
 3000 return
      end
      subroutine obaa
        use cstate
        use cpoted
        use cob
c
c        obaa computes the propagators and the cutoffs of ob-exchanges
c
c
      implicit double precision (a-h,o-z)
c
!      common /cstate/ j,heform,sing,trip,coup,endep,label
!      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,
!     1                znrl,zrel,smev,noced
!      logical heform,sing,trip,coup,endep
!      logical noced
c
!      character*4 label
c
c        common block for all ob-subroutines
c
!      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
!     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
!     2                ez1,ez2,ct(96),wt(96),
!     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
!     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
!     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
!     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
!      logical indc,indxy,indpar
c
c
c
c        further specifications
      dimension deltaq(96,3),cut(96),aaa(96)
      logical indp2
      data iinter/-1/
      data ssmev/-1.d0/
c
c
c        statement functions needed for cray fortran
c
c      sqrt(x)=sqrt(x)
c      abs(x)=abs(x)
c      log(x)=log(x)
c      exp(x)=exp(x)
c      sin(x)=sin(x)
c      cos(x)=cos(x)
c      datan(x)=atan(x)
c
c
c
c
      if (inter.eq.iinter) go to 60
      iinter=inter
c
c        dwn is needed for the eikonal cutoff
c
      dwn=1.d0/wnn(inter)
c
c
c
c
c        delta square
c        ------------
c
c
c
c
   60 if (indxy) go to 1000
      indxy=.true.
      do 65 i=1,nt
      xy2t=xy2*ct(i)
c
c        retardation ignored
c
      deltaq(i,1)=xy2t-xxpyy
c     ----------------------
c
c
c        retardation incorporated
c
      deltaq(i,2)=xy2t-eem12
c     ----------------------
c
c
c        for on shell cutoff  ( ferchlaender )
c
   65 deltaq(i,3)=-xy2t-xxpyy
c     -----------------------
c
c
c
c        propagator
c        ----------
c        ----------
c
c
c
c
 1000 c4=c(4,im)
      iprsp=ic(1,im)
      if (iprsp.lt.0) go to 1400
      if (iprsp.ge.2) go to 1050
      iret=iprsp+1
c
      go to (1010,1020,1030), inter
c
c         propagator for the nn case
 1010 do 1011 i=1,nt
 1011 aa(i)=wt(i)/(c4-deltaq(i,iret))
      go to 1500
c         propagator for the nd case
 1020 do 1021 i=1,nt
      omq=c4-deltaq(i,iret)
      om=sqrt(omq)
 1021 aa(i)=wt(i)*(1.d0/omq+1.d0/(om*(om+ez1)))*0.5d0
      go to 1500
c         propagator for the dd case
 1030 do 1031 i=1,nt
      omq=c4-deltaq(i,iret)
      om=sqrt(omq)
 1031 aa(i)=wt(i)/(om*(om+ez1))
      go to 1500
c
c
c        starting energy dependent propagator
c
c
 1050 ispm=ic(3,im)
      go to (1100,1200,1300),inter
c
c        the propagator for the nn case
c
 1100 if (ispm.gt.2) go to 1110
      do 1105 i=1,nt
      omq=c4-deltaq(i,1)
      om=sqrt(omq)
 1105 aa(i)=wt(i)/(om*(om+ez1))
      go to 1500
c
 1110 do 1115 i=1,nt
      omq=c4-deltaq(i,1)
      om=sqrt(omq)
      oms=om
      if (abs(ez1).lt.1.d-12) go to 1115
      oms=oms+smp(-deltaq(i,1),qfmev,ispm)
 1115 aa(i)=wt(i)/(om*(oms+ez1))
      go to 1500
c
c        the propagator for the nd case
c
 1200 if (ispm.gt.2) go to 1210
      do 1205 i=1,nt
      omq=c4-deltaq(i,1)
      om=sqrt(omq)
 1205 aa(i)=wt(i)*(1.d0/(om*(om+ez1))+1.d0/(om*(om+ez2)))*0.5d0
      go to 1500
c
 1210 do 1215 i=1,nt
      omq=c4-deltaq(i,1)
      om=sqrt(omq)
      oms=om
      oms=oms+smp(-deltaq(i,1),qfmev,ispm)
 1215 aa(i)=wt(i)*(1.d0/(om*(oms+ez1))+1.d0/(om*(oms+ez2)))*0.5d0
      go to 1500
c
c        the propagator for the dd case
c
 1300 if (ispm.gt.2) go to 1310
      do 1305 i=1,nt
      omq=c4-deltaq(i,1)
      om=sqrt(omq)
 1305 aa(i)=wt(i)/(om*(om+ez1))
      go to 1500
c
 1310 do 1315 i=1,nt
      omq=c4-deltaq(i,1)
      om=sqrt(omq)
      oms=om
      oms=oms+smp(-deltaq(i,1),qfmev,ispm)
 1315 aa(i)=wt(i)/(om*(oms+ez1))
      go to 1500
c
c
c        "no propagator"
c
 1400 do 1405 i=1,nt
 1405 aa(i)=wt(i)
c
c
 1500 continue
c
c
c
c
c
c        cut-offs
c        --------
c        --------
c
c
c
c
      mi=4
      mm=5
c
c
  999 ityp=ic(mi,im)
      if (ityp.eq.0) go to 2000
      iprspc=ic(mi+1,im)
      iret=iprspc+1
      if (iprspc.eq.3) iret=1
      go to (100,100,300,400,400,400,700,800,900),ityp
c
c
c
c
c        cut-off of dipole type
c        **********************
c
c
c
c
  100 c5=c(mm,im)
      c6=c(mm+1,im)
      nexp=ic(mi+2,im)
c
      do 105 i=1,nt
c
  105 aaa(i)=c5/(c6-deltaq(i,iret))
c     -------------------------
c
      do 106 ii=1,nexp
      do 106 i=1,nt
  106 aa(i)=aa(i)*aaa(i)
c
      if (iprspc.le.2) go to 120
c
      do 110 i=1,nt
c
  110 aaa(i)=c5/(c6-deltaq(i,iprspc))
c     ----------------------------
c
      do 112 ii=1,nexp
      do 112 i=1,nt
  112 aa(i)=aa(i)*aaa(i)
c
c
  120 mi=mi+3
      mm=mm+2
      go to 999
c
c
c
c
c        cut-off of regge type /schierholz/
c        **********************************
c
c
c
c
  300 ax=ex*c(mm,im)
      ay=ey*c(mm,im)
      expo=log((ax+sqrt(ax*ax-1.d0))*(ay+sqrt(ay*ay-1.d0)))
      expo1=c(mm+1,im)*expo
      expo2=c(mm+2,im)*expo
      do 305 i=1,nt
      expon=expo1+expo2*deltaq(i,iret)
      if (expon.lt.-50.d0) go to 302
c
      aa(i)=aa(i)*exp(expon)
c     ---------------------
c
c
c     ---------------------
c
      go to 305
  302 aa(i)=0.d0
  305 continue
      mi=mi+2
      mm=mm+3
      go to 999
c
c
c
c
c        eikonal form factor
c        *******************
c
c
c
c
  400 ieik=ityp-3
c
      eikc5=c(mm,im)
      do 407 i=1,nt
      expon=0.d0
      ieik3=0
      go to (401,402,403),ieik
c
c        t-form
  401 d=-deltaq(i,iret)
      go to 404
c
c        u-form
  402 d=2.d0*xy2*ct(i)-deltaq(i,iret)
      go to 404
c
c        t-form * u-form
  403 ieik3=ieik3+1
      go to (401,402,405),ieik3
c
  404 d1=sqrt(d)
      d2=sqrt(4.d0+d)
c
      expon=eikc5*(2.d0+d)/(d1*d2)*log(0.5d0*(d1+d2))+expon
c     ------------------------------------------------
c
      go to (405,405,403),ieik
  405 if (expon.lt.-50.d0) go to 406
c
      cut(i)=exp(expon)
c     ------------------
c
      go to 407
  406 cut(i)=0.d0
  407 continue
c
c        get or calculate normalization factor
c
      go to (411,412,412),ieik
c
c        normalization of t-form
  411 c6=c(mm+1,im)
      go to 490
c
c        normalization of u-form and t- * u-form
  412 ieik2=ic(mi+2,im)
      go to (421,422,423,421,422,423),ieik2
c
  421 c6=c(mm+2,im)
      go to 490
c
  422 if (smev.eq.ssmev) go to 442
      ssmev=smev
      do 441 il=1,ime
  441 indc(2,il)=.false.
  442 if (indc(2,im)) go to 421
      indc(2,im)=.true.
      ess=smev*dwn
      ess=ess*ess
      go to 450
c
  423 ess=4.d0*ex*ey
c
  450 d=-(4.d0-ess)
      if (ieik2.le.3) d=d+c(4,im)
      if (d.eq.0.d0) go to 460
      d2=sqrt(4.d0+d)
      if (d.lt.0.d0) go to 470
      d1=sqrt(d)
      c6=exp(-eikc5*(2.d0+d)/(d1*d2)*log(0.5d0*(d1+d2)))
      go to 480
c
  460 c6=exp(-eikc5*0.5d0)
      go to 480
c
  470 d1=sqrt(-d)
c     c6=exp(-eikc5*(2.d0+d)/(d1*d2)*datan(d1/d2))
      c6=exp(-eikc5*(2.d0+d)/(d1*d2)*atan2(d1,d2))
      print*,atan(d1/d2),atan2(d1,d2)
c
  480 if (ieik.eq.3) c6=c6*c(mm+1,im)
      if (ieik2.eq.3.or.ieik2.eq.6) go to 490
      c(mm+2,im)=c6
c
c        compute form factor
c
  490 do 495 i=1,nt
  495 aa(i)=aa(i)*cut(i)*c6
c     -------------------
c
      mi=mi+3
      mm=mm+3
      go to 999
c
c
c        exponential form factor
c        ***********************
c
c
  700 c5=c(mm,im)
      do 705 i=1,nt
c
      expo=deltaq(i,iret)*c5
c     ----------------------
c
      if (expo.lt.-50.d0) go to 704
c
      aa(i)=aa(i)*exp(expo)
c     ----------------------
c
      go to 705
  704 aa(i)=0.d0
  705 continue
      mi=mi+2
      mm=mm+1
      go to 999
c
c
c        cloudy bag form factor
c        ***********************
c
c
  800 c5=c(mm,im)
      nexp=ic(mi+2,im)
      do 805 i=1,nt
c
      arg=sqrt(-deltaq(i,iret))*c5
      argc=arg*arg*arg
  805 aaa(i)=3.d0*(sin(arg)-arg*cos(arg))/argc
c
      do 810 ii=1,nexp
      do 810 i=1,nt
  810 aa(i)=aa(i)*aaa(i)
c
      mi=mi+3
      mm=mm+1
      go to 999
c
c
c        propagator of mass-distributed meson
c        ************************************
c
c
  900 c5=c(mm,im)
      c6=c(mm+1,im)
      nspin=ic(mi+2,im)
      indp2=.false.
      if (iprspc.le.1) go to 901
      indp2=.true.
      iret=1
c
  901 do 915 i=1,nt
      d=c6-deltaq(i,iret)
      if (nspin.eq.0) go to 903
      d1=-d*deltaq(i,iret)/c4
  903 d=sqrt(d)
      if (nspin.eq.0) go to 907
      do 905 ii=1,nspin
  905 d=d*d1
c
  907 omq=c4-deltaq(i,iret)+c5*d
      if (indp2) go to 910
      aa(i)=aa(i)/omq
      go to 915
c
  910 om=sqrt(omq)
      aa(i)=aa(i)/(om*(om+ez1))
c
  915 continue
c
      mi=mi+3
      mm=mm+2
      go to 999
c
c
c
c
 2000 return
      end
      subroutine legp (pj,pjm1,x,j)
c
c
c        subroutine legp   computes the legendre polynoms
c
      double precision pj,pjm1,x,a,b
c
c
c
c        compute legendre polynom for j equals zero
c
c
c
      pjm1=1.d0
      if (j.gt.0) go to 1
      pj=1.d0
      return
c
c
c
c        compute legendre polynoms for j equals one
c
c
c
    1 pj=x
      if (j.eq.1) return
c
c
c
c        compute legendre polynom for j greater or equal two
c
c
c
      do 2 i=2,j
      a=x*pj
      b=a-pjm1
      pjm1=pj
    2 pj=-b/float(i)+b+a
c
c
      return
      end
      subroutine legpc (pj,pjm1,x,j,n)
c
c
c        subroutine legp  computes the legendre polynoms
c
      implicit double precision (a-h,o-z)
c
      dimension pj(1),pjm1(1),x(1)
c
c
c        statementfunctions for cray fortran
c
c      float(i)=float(i)
c
c
c        compute legendre ploynom for j equals zero
c
c
c
      do 10 i=1,n
   10 pjm1(i)=1.d0
c
      if(j.gt.0) go to 100
      do 20 i=1,n
   20 pj(i)=1.d0
      go to 9999
c
c
c        compute legendre polynom for j equals one
c
c
  100 do 30 i=1,n
   30 pj(i)=x(i)
      if (j.eq.1) go to 9999
c
c
c        compute legendre polynom for j greater or equal two
c
c
      do 250 ii=2,j
      do 200 i=1,n
      a=x(i)*pj(i)
      b=a-pjm1(i)
      pjm1(i)=pj(i)
  200 pj(i)=-b/float(ii)+b+a
  250 continue
c
c
 9999 return
      end
      double precision function 
     &spe (qq,qfq,ua,ws,uc,ud,wn,iprop,ispex,ipaho)
c
c        single particle energy of one particle below or above
c        the fermi surface of a medium
c
c        spe is called by:
c                          matbnd,
c                          matusf,
c                          matgmt,
c                          obaa,
c                          tbibnn,
c                          tbibd,
c                          tbaann,
c                          tbainn.
c
c        ispe = 0 : same as ispe=2
c        ispe = 1 : same as ispe=2
c        ispe = 2 : continous spectrum with free energies below and
c                   above fermi surface.
c        ispe = 3 : same as ispe=2 plus constant shift of particle
c                   potential above fermi surface by uc.
c        ispe = 4 : conventional choice for sp energy:
c                   bound energy below the fermi surface and
c                   free energy above the fermi surface (i.e. gap).
c        ispe = 5 : continous spectrum with bound energy below the ferm
c                   surface and continous continuation above the fermi
c                   surface until u=0., free energies after
c        ispe = 6 : continuous for ever
c        ispe = 7 : continuous choice with k-dependent parameters
c                   according to fua in case of iprop=2;
c                   in case of iprop=1 same as ispe=6.
c
c
      implicit double precision (a-h,o-z)
c
c
      external fua,fub
c
      data wwn/938.926d0/
c
c
c        statement function needed for cray fortran
c
c      sqrt(x)=sqrt(x)
c      dmin1(x,y)=amin1(x,y)
c
c
      ispe=ispex
      if (ispe.lt.2) ispe=2
c
c
c
c
      go to (1000,2000),iprop
c
c
c
c
 1000 if (ipaho.eq.0) go to 1001
      go to (1100,1200),ipaho
 1001 if (qq.gt.qfq) go to 1200
c
c
 1100 go to (1110,1110,1110,1140,1140,1140,1170),ispe
 1110 spe=0.5d0*qq/wn+wn-wwn
      go to 9000
 
 1140 spe=0.5d0*qq/ws+ua
      if (wn.ne.wwn) spe=spe+0.5d0*qq*(1.d0/wn-1.d0/wwn)+wn-wwn
      go to 9000
 1170 wsk=wn+fua(qq,qfq,wn)
      spe=0.5d0*qq/wn+wsk+fub(qq,qfq,wn)
      go to 9000
c
c
 1200 go to (1210,1210,1210,1210,1250,1140,1170),ispe
 1210 spe=0.5d0*qq/wn+wn-wwn
      go to 8000
 1250 spe1=0.5d0*qq/ws+ua
      if (wn.ne.wwn) spe1=spe1+0.5d0*qq*(1.d0/wn-1.d0/wwn)+wn-wwn
      spe2=0.5d0*qq/wn+wn-wwn
      spe=min(spe1,spe2)
      go to 9000
c
c
c
c
 2000 if (ipaho.eq.0) go to 2001
      go to (2100,2200),ipaho
 2001 if (qq.gt.qfq) go to 2200
c
c
 2100 go to (2110,2110,2110,2140,2140,2140,2270),ispe
 2110 spe=sqrt(qq+wn*wn)
      go to 9000
 2140 spe=sqrt(qq+ws*ws)+wwn-ws+ua
      if (wn.ne.wwn) spe=spe+sqrt(qq+wn*wn)-sqrt(qq+wwn*wwn)
      go to 9000
c
c
 2200 go to (2210,2210,2210,2210,2250,2140,2270),ispe
 2210 spe=sqrt(qq+wn*wn)
      go to 8000
 2250 spe1=sqrt(qq+ws*ws)+wwn-ws+ua
      if (wn.ne.wwn) spe1=spe1+sqrt(qq+wn*wn)-sqrt(qq+wwn*wwn)
      spe2=sqrt(qq+wn*wn)
      spe=min(spe1,spe2)
      go to 9000
 2270 wsk=wn+fua(qq,qfq,wn)
      spe=sqrt(qq+wsk*wsk)+fub(qq,qfq,wn)
c****      if (qq.lt.25.*qfq) go to 9000
c****      spe=sqrt(qq+wn*wn)
      go to 9000
c
c
c
c
 8000 if (mod(ispe,2).eq.1) spe=spe+uc
c
c
c
c
 9000 return
      end
      double precision function smp (qq,qfmev,ismp)
        use crdwrt
!      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c      real*8 qq,qfmev
      logical index
      data index/.false./
      data wn/938.9d0/,uf/197.3217d0/
      data dwn/0.1065d-02/,duf/0.50679d-02/,wnduf/4.758d0/
      if (ismp.le.2) go to 1000
 1000 smp=0.d0
 2000 return
      end
      double precision function fua (qq,qfq,wn)
        use cfuab
c
c        functional form of a;
c        to be applied in case ispe = 7 only.
c
c
      implicit double precision (a-h,o-z)
c
c
!      common /cfuab/ a(5),b(5),c(5)
c
c        statement function
c
c      exp(x)=exp(x)
c
c
c
c
      qqa=qq/qfq
      aa2=(a(2)/a(1))
      fua=a(1)+a(2)*qqa
c
      if (aa2.gt.0.d0) go to 200
c
      if (fua.gt.0.d0) fua=0.d0
      go to 1000
c
c
  200 aa3=aa2/81.d0
      fua=fua*exp(-aa3*qqa*qqa)
c
c
c****      ak=a(1)+a(2)*qqa
c****      ck=c(1)+c(2)*qqa
c
c****      fua=(ak-ck)/(wn+ck)*wn
c
c
c
c
c
 1000 return
      end
      double precision function fub (qq,qfq,wn)
        use cfuab
c
c        functional form of b;
c        to be applied in case ispe = 7 only.
c
c
      implicit double precision (a-h,o-z)
c
!      common /cfuab/ a(5),b(5),c(5)
c
c
c        statement function
c
c      exp(x)=exp(x)
c
c
c
      qqa=qq/qfq
      bb2=(b(2)/b(1))
      fub=b(1)+b(2)*qqa
c
      if (bb2.gt.0.d0) go to 200
c
      if (fub.lt.0.d0) fub=0.d0
      go to 1000
c
c
  200 bb3=bb2/81.d0
      fub=fub*exp(-bb3*qqa*qqa)
c
c
c****      ak=a(1)+a(2)*qqa
c****      bk=b(1)+b(2)*qqa
c****      ck=c(1)+c(2)*qqa
c
c****      cpk=1.d0+ck/wn
c****      wsk=wn+ak
c
c****      fub=(bk*wn+(sqrt(qq*cpk*cpk+wsk*wsk)+bk)*ck)/(wn+ck)
c
c
c
c
c
 1000 return
      end
      subroutine gset(ax,bx,n,z,w)
      implicit double precision (a-h,o-z)
c 
c     n-point gauss zeros and weights for the interval (ax,bx) are
c           stored in  arrays z and w respectively.
c         
      dimension     a(273),x(273),ktab(96)
      dimension z(2),w(2)
c 
c-----table of initial subscripts for n=2(1)16(4)96
      data ktab(2)/1/
      data ktab(3)/2/
      data ktab(4)/4/
      data ktab(5)/6/
      data ktab(6)/9/
      data ktab(7)/12/
      data ktab(8)/16/
      data ktab(9)/20/
      data ktab(10)/25/
      data ktab(11)/30/
      data ktab(12)/36/
      data ktab(13)/42/
      data ktab(14)/49/
      data ktab(15)/56/
      data ktab(16)/64/
      data ktab(20)/72/
      data ktab(24)/82/
      data ktab(28)/82/
      data ktab(32)/94/
      data ktab(36)/94/
      data ktab(40)/110/
      data ktab(44)/110/
      data ktab(48)/130/
      data ktab(52)/130/
      data ktab(56)/130/
      data ktab(60)/130/
      data ktab(64)/154/
      data ktab(68)/154/
      data ktab(72)/154/
      data ktab(76)/154/
      data ktab(80)/186/
      data ktab(84)/186/
      data ktab(88)/186/
      data ktab(92)/186/
      data ktab(96)/226/
c 
c-----table of abscissae (x) and weights (a) for interval (-1,+1).
c 
c**** n=2
      data x(1)/0.577350269189626  e0/, a(1)/1.000000000000000  e0/
c**** n=3
      data x(2)/0.774596669241483  e0/, a(2)/0.555555555555556  e0/
      data x(3)/0.000000000000000  e0/, a(3)/0.888888888888889  e0/
c**** n=4
      data x(4)/0.861136311594053  e0/, a(4)/0.347854845137454  e0/
      data x(5)/0.339981043584856  e0/, a(5)/0.652145154862546  e0/
c**** n=5
      data x(6)/0.906179845938664  e0/, a(6)/0.236926885056189  e0/
      data x(7)/0.538469310105683  e0/, a(7)/0.478628670499366  e0/
      data x(8)/0.000000000000000  e0/, a(8)/0.568888888888889  e0/
c**** n=6
      data x(9)/0.932469514203152  e0/, a(9)/0.171324492379170  e0/
      data x(10)/0.661209386466265 e0/, a(10)/0.360761573048139 e0/
      data x(11)/0.238619186083197 e0/, a(11)/0.467913934572691 e0/
c**** n=7
      data x(12)/0.949107912342759 e0/, a(12)/0.129484966168870 e0/
      data x(13)/0.741531185599394 e0/, a(13)/0.279705391489277 e0/
      data x(14)/0.405845151377397 e0/, a(14)/0.381830050505119 e0/
      data x(15)/0.000000000000000 e0/, a(15)/0.417959183673469 e0/
c**** n=8
      data x(16)/0.960289856497536 e0/, a(16)/0.101228536290376 e0/
      data x(17)/0.796666477413627 e0/, a(17)/0.222381034453374 e0/
      data x(18)/0.525532409916329 e0/, a(18)/0.313706645877887 e0/
      data x(19)/0.183434642495650 e0/, a(19)/0.362683783378362 e0/
c**** n=9
      data x(20)/0.968160239507626 e0/, a(20)/0.081274388361574 e0/
      data x(21)/0.836031107326636 e0/, a(21)/0.180648160694857 e0/
      data x(22)/0.613371432700590 e0/, a(22)/0.260610696402935 e0/
      data x(23)/0.324253423403809 e0/, a(23)/0.312347077040003 e0/
      data x(24)/0.000000000000000 e0/, a(24)/0.330239355001260 e0/
c**** n=10
      data x(25)/0.973906528517172 e0/, a(25)/0.066671344308688 e0/
      data x(26)/0.865063366688985 e0/, a(26)/0.149451349150581 e0/
      data x(27)/0.679409568299024 e0/, a(27)/0.219086362515982 e0/
      data x(28)/0.433395394129247 e0/, a(28)/0.269266719309996 e0/
      data x(29)/0.148874338981631 e0/, a(29)/0.295524224714753 e0/
c**** n=11
      data x(30)/0.978228658146057 e0/, a(30)/0.055668567116174 e0/
      data x(31)/0.887062599768095 e0/, a(31)/0.125580369464905 e0/
      data x(32)/0.730152005574049 e0/, a(32)/0.186290210927734 e0/
      data x(33)/0.519096129206812 e0/, a(33)/0.233193764591990 e0/
      data x(34)/0.269543155952345 e0/, a(34)/0.262804544510247 e0/
      data x(35)/0.000000000000000 e0/, a(35)/0.272925086777901 e0/
c**** n=12
      data x(36)/0.981560634246719 e0/, a(36)/0.047175336386512 e0/
      data x(37)/0.904117256370475 e0/, a(37)/0.106939325995318 e0/
      data x(38)/0.769902674194305 e0/, a(38)/0.160078328543346 e0/
      data x(39)/0.587317954286617 e0/, a(39)/0.203167426723066 e0/
      data x(40)/0.367831498998180 e0/, a(40)/0.233492536538355 e0/
      data x(41)/0.125233408511469 e0/, a(41)/0.249147045813403 e0/
c**** n=13
      data x(42)/0.984183054718588 e0/, a(42)/0.040484004765316 e0/
      data x(43)/0.917598399222978 e0/, a(43)/0.092121499837728 e0/
      data x(44)/0.801578090733310 e0/, a(44)/0.138873510219787 e0/
      data x(45)/0.642349339440340 e0/, a(45)/0.178145980761946 e0/
      data x(46)/0.448492751036447 e0/, a(46)/0.207816047536889 e0/
      data x(47)/0.230458315955135 e0/, a(47)/0.226283180262897 e0/
      data x(48)/0.000000000000000 e0/, a(48)/0.232551553230874 e0/
c**** n=14
      data x(49)/0.986283808696812 e0/, a(49)/0.035119460331752 e0/
      data x(50)/0.928434883663574 e0/, a(50)/0.080158087159760 e0/
      data x(51)/0.827201315069765 e0/, a(51)/0.121518570687903 e0/
      data x(52)/0.687292904811685 e0/, a(52)/0.157203167158194 e0/
      data x(53)/0.515248636358154 e0/, a(53)/0.185538397477938 e0/
      data x(54)/0.319112368927890 e0/, a(54)/0.205198463721296 e0/
      data x(55)/0.108054948707344 e0/, a(55)/0.215263853463158 e0/
c**** n=15
      data x(56)/0.987992518020485 e0/, a(56)/0.030753241996117 e0/
      data x(57)/0.937273392400706 e0/, a(57)/0.070366047488108 e0/
      data x(58)/0.848206583410427 e0/, a(58)/0.107159220467172 e0/
      data x(59)/0.724417731360170 e0/, a(59)/0.139570677926154 e0/
      data x(60)/0.570972172608539 e0/, a(60)/0.166269205816994 e0/
      data x(61)/0.394151347077563 e0/, a(61)/0.186161000015562 e0/
      data x(62)/0.201194093997435 e0/, a(62)/0.198431485327111 e0/
      data x(63)/0.000000000000000 e0/, a(63)/0.202578241925561 e0/
c**** n=16
      data x(64)/0.989400934991650 e0/, a(64)/0.027152459411754 e0/
      data x(65)/0.944575023073233 e0/, a(65)/0.062253523938648 e0/
      data x(66)/0.865631202387832 e0/, a(66)/0.095158511682493 e0/
      data x(67)/0.755404408355003 e0/, a(67)/0.124628971255534 e0/
      data x(68)/0.617876244402644 e0/, a(68)/0.149595988816577 e0/
      data x(69)/0.458016777657227 e0/, a(69)/0.169156519395003 e0/
      data x(70)/0.281603550779259 e0/, a(70)/0.182603415044924 e0/
      data x(71)/0.095012509837637 e0/, a(71)/0.189450610455069 e0/
c**** n=20
      data x(72)/0.993128599185094 e0/, a(72)/0.017614007139152 e0/
      data x(73)/0.963971927277913 e0/, a(73)/0.040601429800386 e0/
      data x(74)/0.912234428251325 e0/, a(74)/0.062672048334109 e0/
      data x(75)/0.839116971822218 e0/, a(75)/0.083276741576704 e0/
      data x(76)/0.746331906460150 e0/, a(76)/0.101930119817240 e0/
      data x(77)/0.636053680726515 e0/, a(77)/0.118194531961518 e0/
      data x(78)/0.510867001950827 e0/, a(78)/0.131688638449176 e0/
      data x(79)/0.373706088715419 e0/, a(79)/0.142096109318382 e0/
      data x(80)/0.227785851141645 e0/, a(80)/0.149172986472603 e0/
      data x(81)/0.076526521133497 e0/, a(81)/0.152753387130725 e0/
c**** n=24
      data x(82)/0.995187219997021 e0/, a(82)/0.012341229799987 e0/
      data x(83)/0.974728555971309 e0/, a(83)/0.028531388628933 e0/
      data x(84)/0.938274552002732 e0/, a(84)/0.044277438817419 e0/
      data x(85)/0.886415527004401 e0/, a(85)/0.059298584915436 e0/
      data x(86)/0.820001985973902 e0/, a(86)/0.073346481411080 e0/
      data x(87)/0.740124191578554 e0/, a(87)/0.086190161531953 e0/
      data x(88)/0.648093651936975 e0/, a(88)/0.097618652104113 e0/
      data x(89)/0.545421471388839 e0/, a(89)/0.107444270115965 e0/
      data x(90)/0.433793507626045 e0/, a(90)/0.115505668053725 e0/
      data x(91)/0.315042679696163 e0/, a(91)/0.121670472927803 e0/
      data x(92)/0.191118867473616 e0/, a(92)/0.125837456346828 e0/
      data x(93)/0.064056892862605 e0/, a(93)/0.127938195346752 e0/
c**** n=32
      data x(94)/0.997263861849481 e0/, a(94)/0.007018610009470 e0/
      data x(95)/0.985611511545268 e0/, a(95)/0.016274394730905 e0/
      data x(96)/0.964762255587506 e0/, a(96)/0.025392065309262 e0/
      data x(97)/0.934906075937739 e0/, a(97)/0.034273862913021 e0/
      data x(98)/0.896321155766052 e0/, a(98)/0.042835898022226 e0/
      data x(99)/0.849367613732569 e0/, a(99)/0.050998059262376 e0/
      data x(100)/0.794483795967942e0/, a(100)/0.058684093478535e0/
      data x(101)/0.732182118740289e0/, a(101)/0.065822222776361e0/
      data x(102)/0.663044266930215e0/, a(102)/0.072345794108848e0/
      data x(103)/0.587715757240762e0/, a(103)/0.078193895787070e0/
      data x(104)/0.506899908932229e0/, a(104)/0.083311924226946e0/
      data x(105)/0.421351276130635e0/, a(105)/0.087652093004403e0/
      data x(106)/0.331868602282127e0/, a(106)/0.091173878695763e0/
      data x(107)/0.239287362252137e0/, a(107)/0.093844399080804e0/
      data x(108)/0.144471961582796e0/, a(108)/0.095638720079274e0/
      data x(109)/0.048307665687738e0/, a(109)/0.096540088514727e0/
c**** n=40
      data x(110)/0.998237709710559e0/, a(110)/0.004521277098533e0/
      data x(111)/0.990726238699457e0/, a(111)/0.010498284531152e0/
      data x(112)/0.977259949983774e0/, a(112)/0.016421058381907e0/
      data x(113)/0.957916819213791e0/, a(113)/0.022245849194166e0/
      data x(114)/0.932812808278676e0/, a(114)/0.027937006980023e0/
      data x(115)/0.902098806968874e0/, a(115)/0.033460195282547e0/
      data x(116)/0.865959503212259e0/, a(116)/0.038782167974472e0/
      data x(117)/0.824612230833311e0/, a(117)/0.043870908185673e0/
      data x(118)/0.778305651426519e0/, a(118)/0.048695807635072e0/
      data x(119)/0.727318255189927e0/, a(119)/0.053227846983936e0/
      data x(120)/0.671956684614179e0/, a(120)/0.057439769099391e0/
      data x(121)/0.612553889667980e0/, a(121)/0.061306242492928e0/
      data x(122)/0.549467125095128e0/, a(122)/0.064804013456601e0/
      data x(123)/0.483075801686178e0/, a(123)/0.067912045815233e0/
      data x(124)/0.413779204371605e0/, a(124)/0.070611647391286e0/
      data x(125)/0.341994090825758e0/, a(125)/0.072886582395804e0/
      data x(126)/0.268152185007253e0/, a(126)/0.074723169057968e0/
      data x(127)/0.192697580701371e0/, a(127)/0.076110361900626e0/
      data x(128)/0.116084070675255e0/, a(128)/0.077039818164247e0/
      data x(129)/0.038772417506050e0/, a(129)/0.077505947978424e0/
c**** n=48
      data x(130)/0.998771007252426e0/, a(130)/0.003153346052305e0/
      data x(131)/0.993530172266350e0/, a(131)/0.007327553901276e0/
      data x(132)/0.984124583722826e0/, a(132)/0.011477234579234e0/
      data x(133)/0.970591592546247e0/, a(133)/0.015579315722943e0/
      data x(134)/0.952987703160430e0/, a(134)/0.019616160457355e0/
      data x(135)/0.931386690706554e0/, a(135)/0.023570760839324e0/
      data x(136)/0.905879136715569e0/, a(136)/0.027426509708356e0/
      data x(137)/0.876572020274247e0/, a(137)/0.031167227832798e0/
      data x(138)/0.843588261624393e0/, a(138)/0.034777222564770e0/
      data x(139)/0.807066204029442e0/, a(139)/0.038241351065830e0/
      data x(140)/0.767159032515740e0/, a(140)/0.041545082943464e0/
      data x(141)/0.724034130923814e0/, a(141)/0.044674560856694e0/
      data x(142)/0.677872379632663e0/, a(142)/0.047616658492490e0/
      data x(143)/0.628867396776513e0/, a(143)/0.050359035553854e0/
      data x(144)/0.577224726083972e0/, a(144)/0.052890189485193e0/
      data x(145)/0.523160974722233e0/, a(145)/0.055199503699984e0/
      data x(146)/0.466902904750958e0/, a(146)/0.057277292100403e0/
      data x(147)/0.408686481990716e0/, a(147)/0.059114839698395e0/
      data x(148)/0.348755886292160e0/, a(148)/0.060704439165893e0/
      data x(149)/0.287362487355455e0/, a(149)/0.062039423159892e0/
      data x(150)/0.224763790394689e0/, a(150)/0.063114192286254e0/
      data x(151)/0.161222356068891e0/, a(151)/0.063924238584648e0/
      data x(152)/0.097004699209462e0/, a(152)/0.064466164435950e0/
      data x(153)/0.032380170962869e0/, a(153)/0.064737696812683e0/
c**** n=64
      data x(154)/0.999305041735772e0/, a(154)/0.001783280721696e0/
      data x(155)/0.996340116771955e0/, a(155)/0.004147033260562e0/
      data x(156)/0.991013371476744e0/, a(156)/0.006504457968978e0/
      data x(157)/0.983336253884625e0/, a(157)/0.008846759826363e0/
      data x(158)/0.973326827789910e0/, a(158)/0.011168139460131e0/
      data x(159)/0.961008799652053e0/, a(159)/0.013463047896718e0/
      data x(160)/0.946411374858402e0/, a(160)/0.015726030476024e0/
      data x(161)/0.929569172131939e0/, a(161)/0.017951715775697e0/
      data x(162)/0.910522137078502e0/, a(162)/0.020134823153530e0/
      data x(163)/0.889315445995114e0/, a(163)/0.022270173808383e0/
      data x(164)/0.865999398154092e0/, a(164)/0.024352702568710e0/
      data x(165)/0.840629296252580e0/, a(165)/0.026377469715054e0/
      data x(166)/0.813265315122797e0/, a(166)/0.028339672614259e0/
      data x(167)/0.783972358943341e0/, a(167)/0.030234657072402e0/
      data x(168)/0.752819907260531e0/, a(168)/0.032057928354851e0/
      data x(169)/0.719881850171610e0/, a(169)/0.033805161837141e0/
      data x(170)/0.685236313054233e0/, a(170)/0.035472213256882e0/
      data x(171)/0.648965471254657e0/, a(171)/0.037055128540240e0/
      data x(172)/0.611155355172393e0/, a(172)/0.038550153178615e0/
      data x(173)/0.571895646202634e0/, a(173)/0.039953741132720e0/
      data x(174)/0.531279464019894e0/, a(174)/0.041262563242623e0/
      data x(175)/0.489403145707052e0/, a(175)/0.042473515123653e0/
      data x(176)/0.446366017253464e0/, a(176)/0.043583724529323e0/
      data x(177)/0.402270157963991e0/, a(177)/0.044590558163756e0/
      data x(178)/0.357220158337668e0/, a(178)/0.045491627927418e0/
      data x(179)/0.311322871990210e0/, a(179)/0.046284796581314e0/
      data x(180)/0.264687162208767e0/, a(180)/0.046968182816210e0/
      data x(181)/0.217423643740007e0/, a(181)/0.047540165714830e0/
      data x(182)/0.169644420423992e0/, a(182)/0.047999388596458e0/
      data x(183)/0.121462819296120e0/, a(183)/0.048344762234802e0/
      data x(184)/0.072993121787799e0/, a(184)/0.048575467441503e0/
      data x(185)/0.024350292663424e0/, a(185)/0.048690957009139e0/
c**** n=80
      data x(186)/0.999553822651630e0/, a(186)/0.001144950003186e0/
      data x(187)/0.997649864398237e0/, a(187)/0.002663533589512e0/
      data x(188)/0.994227540965688e0/, a(188)/0.004180313124694e0/
      data x(189)/0.989291302499755e0/, a(189)/0.005690922451403e0/
      data x(190)/0.982848572738629e0/, a(190)/0.007192904768117e0/
      data x(191)/0.974909140585727e0/, a(191)/0.008683945269260e0/
      data x(192)/0.965485089043799e0/, a(192)/0.010161766041103e0/
      data x(193)/0.954590766343634e0/, a(193)/0.011624114120797e0/
      data x(194)/0.942242761309872e0/, a(194)/0.013068761592401e0/
      data x(195)/0.928459877172445e0/, a(195)/0.014493508040509e0/
      data x(196)/0.913263102571757e0/, a(196)/0.015896183583725e0/
      data x(197)/0.896675579438770e0/, a(197)/0.017274652056269e0/
      data x(198)/0.878722567678213e0/, a(198)/0.018626814208299e0/
      data x(199)/0.859431406663111e0/, a(199)/0.019950610878141e0/
      data x(200)/0.838831473580255e0/, a(200)/0.021244026115782e0/
      data x(201)/0.816954138681463e0/, a(201)/0.022505090246332e0/
      data x(202)/0.793832717504605e0/, a(202)/0.023731882865930e0/
      data x(203)/0.769502420135041e0/, a(203)/0.024922535764115e0/
      data x(204)/0.744000297583597e0/, a(204)/0.026075235767565e0/
      data x(205)/0.717365185362099e0/, a(205)/0.027188227500486e0/
      data x(206)/0.689637644342027e0/, a(206)/0.028259816057276e0/
      data x(207)/0.660859898986119e0/, a(207)/0.029288369583267e0/
      data x(208)/0.631075773046871e0/, a(208)/0.030272321759557e0/
      data x(209)/0.600330622829751e0/, a(209)/0.031210174188114e0/
      data x(210)/0.568671268122709e0/, a(210)/0.032100498673487e0/
      data x(211)/0.536145920897131e0/, a(211)/0.032941939397645e0/
      data x(212)/0.502804111888784e0/, a(212)/0.033733214984611e0/
      data x(213)/0.468696615170544e0/, a(213)/0.034473120451753e0/
      data x(214)/0.433875370831756e0/, a(214)/0.035160529044747e0/
      data x(215)/0.398393405881969e0/, a(215)/0.035794393953416e0/
      data x(216)/0.362304753499487e0/, a(216)/0.036373749905835e0/
      data x(217)/0.325664370747701e0/, a(217)/0.036897714638276e0/
      data x(218)/0.288528054884511e0/, a(218)/0.037365490238730e0/
      data x(219)/0.250952358392272e0/, a(219)/0.037776364362001e0/
      data x(220)/0.212994502857666e0/, a(220)/0.038129711314477e0/
      data x(221)/0.174712291832646e0/, a(221)/0.038424993006959e0/
      data x(222)/0.136164022809143e0/, a(222)/0.038661759774076e0/
      data x(223)/0.097408398441584e0/, a(223)/0.038839651059051e0/
      data x(224)/0.058504437152420e0/, a(224)/0.038958395962769e0/
      data x(225)/0.019511383256793e0/, a(225)/0.039017813656306e0/
c**** n=96
      data x(226)/0.999689503883230e0/, a(226)/0.000796792065552e0/
      data x(227)/0.998364375863181e0/, a(227)/0.001853960788946e0/
      data x(228)/0.995981842987209e0/, a(228)/0.002910731817934e0/
      data x(229)/0.992543900323762e0/, a(229)/0.003964554338444e0/
      data x(230)/0.988054126329623e0/, a(230)/0.005014202742927e0/
      data x(231)/0.982517263563014e0/, a(231)/0.006058545504235e0/
      data x(232)/0.975939174585136e0/, a(232)/0.007096470791153e0/
      data x(233)/0.968326828463264e0/, a(233)/0.008126876925698e0/
      data x(234)/0.959688291448742e0/, a(234)/0.009148671230783e0/
      data x(235)/0.950032717784437e0/, a(235)/0.010160770535008e0/
      data x(236)/0.939370339752755e0/, a(236)/0.011162102099838e0/
      data x(237)/0.927712456722308e0/, a(237)/0.012151604671088e0/
      data x(238)/0.915071423120898e0/, a(238)/0.013128229566961e0/
      data x(239)/0.901460635315852e0/, a(239)/0.014090941772314e0/
      data x(240)/0.886894517402420e0/, a(240)/0.015038721026994e0/
      data x(241)/0.871388505909296e0/, a(241)/0.015970562902562e0/
      data x(242)/0.854959033434601e0/, a(242)/0.016885479864245e0/
      data x(243)/0.837623511228187e0/, a(243)/0.017782502316045e0/
      data x(244)/0.819400310737931e0/, a(244)/0.018660679627411e0/
      data x(245)/0.800308744139140e0/, a(245)/0.019519081140145e0/
      data x(246)/0.780369043867433e0/, a(246)/0.020356797154333e0/
      data x(247)/0.759602341176647e0/, a(247)/0.021172939892191e0/
      data x(248)/0.738030643744400e0/, a(248)/0.021966644438744e0/
      data x(249)/0.715676812348967e0/, a(249)/0.022737069658329e0/
      data x(250)/0.692564536642171e0/, a(250)/0.023483399085926e0/
      data x(251)/0.668718310043916e0/, a(251)/0.024204841792364e0/
      data x(252)/0.644163403784967e0/, a(252)/0.024900633222483e0/
      data x(253)/0.618925840125468e0/, a(253)/0.025570036005349e0/
      data x(254)/0.593032364777572e0/, a(254)/0.026212340735672e0/
      data x(255)/0.566510418561397e0/, a(255)/0.026826866725591e0/
      data x(256)/0.539388108324357e0/, a(256)/0.027412962726029e0/
      data x(257)/0.511694177154667e0/, a(257)/0.027970007616848e0/
      data x(258)/0.483457973920596e0/, a(258)/0.028497411065085e0/
      data x(259)/0.454709422167743e0/, a(259)/0.028994614150555e0/
      data x(260)/0.425478988407300e0/, a(260)/0.029461089958167e0/
      data x(261)/0.395797649828908e0/, a(261)/0.029896344136328e0/
      data x(262)/0.365696861472313e0/, a(262)/0.030299915420827e0/
      data x(263)/0.335208522892625e0/, a(263)/0.030671376123669e0/
      data x(264)/0.304364944354496e0/, a(264)/0.031010332586313e0/
      data x(265)/0.273198812591049e0/, a(265)/0.031316425596861e0/
      data x(266)/0.241743156163840e0/, a(266)/0.031589330770727e0/
      data x(267)/0.210031310460567e0/, a(267)/0.031828758894411e0/
      data x(268)/0.178096882367618e0/, a(268)/0.032034456231992e0/
      data x(269)/0.145973714654896e0/, a(269)/0.032206204794030e0/
      data x(270)/0.113695850110665e0/, a(270)/0.032343822568575e0/
      data x(271)/0.081297495464425e0/, a(271)/0.032447163714064e0/
      data x(272)/0.048812985136049e0/, a(272)/0.032516118713868e0/
      data x(273)/0.016276744849602e0/, a(273)/0.032550614492363e0/
c 
c 
c-----test n
      alpha=0.5e0*(ax+bx)
      beta=0.5e0*(bx-ax)
      if( n.lt.1 .or. n.gt.96 ) go to 100
      if(n.ne.1) go to 1
      z(1)=alpha
      w(1)=bx-ax
      return
c 
    1 if (n.le.16) go to 3
      if (n.gt.24) go to 4
      n=4*(n/4)
      go to 3
    4 if (n.gt.48) go to 5
      n=8*(n/8)
      go to 3
    5 n=16*(n/16)
c 
c----- set k equal to initial subscript and store results
    3 k=ktab(n)
      m=n/2
      do 2 j=1,m
      jtab=k-1+j
      wtemp=beta*a(jtab)
      delta=beta*x(jtab)
      z(j)=alpha-delta
      w(j)=wtemp
      jp=n+1-j
      z(jp)=alpha+delta
      w(jp)=wtemp
    2 continue
      if((n-m-m).eq.0) return
      z(m+1)=alpha
      jmid=k+m
      w(m+1)=beta*a(jmid)
      return
c 
  100 zn=n
      write(6,200) zn
  200 format(1h0/////'0error in gset. n has the non-permissible value',
     1e11.3/'0execution terminated.')
      stop
      end

********************************************************************
      subroutine obnna
      use crdwrt
      use cpot
      use cstate
      use cob
c
c
c        one-boson-exchange nn-nn interaction;
c        version which uses numerical integration
c
c        [changed for Sparc/Sun Feb. 1991, ansi-standard]
c
c
c
      implicit double precision (a-h,o-z)
c
c
!      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c        arguments and values of this subroutine
c
!      common /cpot/   v(6),xmev,ymev
!      common /cstate/ j,heform,sing,trip,coup,endep,label
c
!      character*4 label
c
c        this has been the end of the common-blocks containing
c        the arguments and values of this subroutine in the case of
c        no energy-dependence of the potential;
c        in case of energy-dependence look for the common-block /cped/
c        in obaia and obaaa.
c
c        specifications for these two common blocks
c
!      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
!      common /coba/   vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
!     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
!     2                ez1,ez2,ct(96),wt(96),
!     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
!     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
!     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
!     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
!      logical indc,indxy,indpar
c
c
c        further specifications
c
      logical indmg(12)
      logical index
      character*4 mesong(12)
c
      data pi/3.141592653589793d0/
      data  d3/0.3333333333333333d0/
      data td3/0.6666666666666667d0/
      data fd3/1.3333333333333333d0/
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1                   '1-  ','1-t ','1-tt','1-st','1-ss',
     2                    '1+  ','2+  '/
      data indmg/12*.false./
c
c
c        statement functionc needed for cray fortran
c
c      sqrt(x)=sqrt(x)
c
c4
c
c
      inter=1
c
c
c
c
c        call subroutine obpara once and only once
c
c
      if (index) go to 50
      index=.true.
      if (indpar(inter)) go to 40
c
c
      call obpara
c
c
   40 wdd(inter)=0.d0
      iftgo=ift(inter)+1
      dwn=1.d0/wnn(inter)
      iman=imaa(inter)
      imen=imea(inter)
c
c
c        prepare constant over-all factor
c
      fac=1.d0/(2.d0*pi)*dwn*dwn
c     --------------------------
c
c
c
c
c
c
c
c        prepare expressions depending on x and y
c        ----------------------------------------
c        ----------------------------------------
c
c
c
c
   50 xa=xmev*dwn
      ya=ymev*dwn
      indxy=.false.
      x=xa
      xx=x*x
      y=ya
      yy=y*y
      xy2=x*y*2.d0
      xxpyy=xx+yy
      ex=sqrt(1.d0+xx)
      ey=sqrt(1.d0+yy)
      eem12=(ex*ey-1.d0)*2.d0
c
c
c
c
   55 xy=xy2*0.5d0
c     dxy=1.d0/xy
      ee=ex*ey
      ree=sqrt(ee)
      eem1=ee-1.d0
      eme=ex-ey
      emeh=eme*0.5d0
      emehq=emeh*emeh
      eep1=ee+1.d0
       epe=ex+ey
      xxyy=xx*yy
c
c
c
c
c        prepare over-all factor
c
c
      go to (70,71,72),iftgo
c
c        no additional factor
c
   70 fff=fac
      go to 90
c
c        minimal relativity
c
   71 fff=fac/ree
      go to 90
c
c        factor m/e*m/e
c
   72 fff=fac/ee
c
c
c
c
c
c
   90 do 93 iv=1,6
   93 v(iv)=0.d0
      do 95 il=iman,imen
      do 95 iv=1,32
   95 vj(iv,il)=0.d0
c
c
c
c
c        contributions of mesons
c        -----------------------
c        -----------------------
c
c
c
c
      do 1995 img=1,mge
      mg=mggo(img,inter)
      if (mg.eq.0) go to 2000
      me=mgg(mg,inter)
      go to (100,200,300,400,500,600,600,600,900,1000,1100,1200),mg
c
c
c
c        0-  , pseudo-scalar mesons
c        --------------------------
c
c
c
c
  100 mc=1
c
      ff=1.d0
      f(1)=eem1
      f(2)=-xy
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-eme
      f(8)=-f(7)
c
      call obstra(1,1,me)
      go to 1995
c
c
c
c        0-t , pseudo-vector mesons
c        --------------------------
c
c
c
c
  200 mc=1
c
      ff=1.d0
      f(1)=eem1+emehq*(ee+3.d0)
      f(2)=-xy+emehq*xy
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-eme-eme*(emehq+eem1)
      f(8)=-f(7)
c
      call obstra(1,1,me)
      go to 1995
c
c
c
c        0-st, pseudo-scalar mesons in static limit
c        ------------------------------------------
c
c
c
c
  300 mc=1
c
      ff=1.d0
      f(1)=xxpyy*0.5d0
      f(2)=-xy
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-(xx-yy)*0.5d0
      f(8)=-f(7)
c
      call obstra(1,1,me)
      go to 1995
c
c
c
c
c        0+  , scalar mesons
c        -------------------
c
c
c
c
  400 mc=1
c
      ff=1.d0
      f(1)=-eep1
      f(2)=xy
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=epe
      f(8)=f(7)
c
      call obstra(1,1,me)
      go to 1995
c
c
c
c
c        0+st, scalar mesons in static limit
c        -----------------------------------
c
c
c
c
  500 mc=1
c
c
c
c        central term  '1'  only
c
c
c
      ff=1.d0
      f(1)=-2.d0
      f(2)=0.d0
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-f(1)
      f(8)=f(7)
      f(9)=f(2)
      f(10)=f(2)
      f(11)=f(2)
c
c
c
c        central term  ' k**2 + p**2 '
c
c
c
      f(2)=f(2)+xy
      f(10)=f(10)+xy
      f(11)=f(11)-xy
c
c
c
c        spin orbit
c
c
c
      f(4)=f(4)+xy
      f(5)=f(5)+xy
      f(10)=f(10)-xy
      f(11)=f(11)+xy
c
      call obstra(2,1,me)
c
c
c
c        case of additional sigma-l
c
c
c
c     mc=-1
c     xxyy8=xxyy/8.d0
c     f(1)=-xxyy8
c     f(2)=0.d0
c     f(3)=f(1)
c     f(4)=f(2)
c     f(5)=f(2)
c     f(6)=xxyy8
c     f(7)=f(1)
c     f(8)=f(7)
c     f(9)=f(6)*2.d0
c
c     call obstra(4,1,me)
      go to 1995
c
c
c
c
c        1-  , vector mesons
c        -------------------
c
c
c
c
c        vector coupling
c
c
c
c
  600 mc=1
c
      ff=2.d0
      f(1)=eem1+ee
      f(2)=0.d0
      f(3)=ee
      f(4)=xy
      f(5)=xy2
      f(6)=1.d0
      f(7)=-ey
      f(8)=-ex
c
      call obstra(1,1,me)
      if (mg.eq.7) go to 720
      if (mg.eq.8) go to 810
c
c
c
c
c        tensor coupling
c
c
c
c
      mc=2
c
c
      ff=0.5d0
        e1=-xxpyy+xxyy
      f(1)=eem1*(xxpyy+6.d0)+e1
      f(2)=-(xxpyy+4.d0)*xy
      f(3)=eem1*(xxpyy+2.d0)+e1
        e2=xxpyy+ee
      f(4)=-(1.d0+e2)*xy
      f(5)= (3.d0-e2)*xy
      f(6)=eem1*(xxpyy-2.d0)-xxpyy
        e3=2.d0*eme
        e4=yy*(ey-2.d0*ex)+xx*(ex-2.d0*ey)
      f(7)= e3-e4
      f(8)=-e3-e4
c        factors for additional terms
      f(9)=-xxyy
      f(10)=(1.d0+ee)*xy
      f(11)=-epe*xy
c
      call obstra(2,1,me)
c
c
c
c
c        vector-tensor coupling
c
c
c
c
      mc=3
c
      ff=2.d0
      f(1)=3.d0*eem1-xxpyy
      f(2)=-xy
      f(3)=eem1-xxpyy
      f(4)=-f(2)
      f(5)=3.d0*xy
      f(6)=-(eem1+xxpyy)
        e1=yy*ex+xx*ey
      f(7)= eme+e1
      f(8)=-eme+e1
c
      call obstra(1,1,me)
      go to 1995
c
c
c
c
c        1-t , vector mesons a la gtg
c        ----------------------------
c
c
c
c
c        tensor coupling
c
c
c
c
  720 mc=2
c
      ff=0.25d0
      f(1)=(3.d0*ee+1.d0)*xxpyy
      f(2)=-(6.d0*ee+2.d0-xxpyy)*xy
      f(3)=eem1*xxpyy+4.d0*xxyy
      f(4)=-(4.d0*ee+xxpyy)*xy
      f(5)=(4.d0-3.d0*xxpyy)*xy
      f(6)=6.d0*xxyy-(ee+3.d0)*xxpyy
      f(7)=(ex+3.d0*ey)*xx+eme*yy
      f(8)=(ey+3.d0*ex)*yy-eme*xx
c        factors for additional terms
      f(9)=-2.d0*xxyy
      f(10)=eep1*xy2
      f(11)=-epe*xy2
c
      call obstra(2,1,me)
c
c
c
c
c        vector-tensor coupling
c
c
c
c
      mc=3
c
      ff=1.d0
      f(1)=xxpyy
      f(2)=-xy2
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=6.d0*xy
      f(6)=3.d0*f(3)
      f(7)=(ex*yy+3.d0*ey*xx)
      f(8)=(ey*xx+3.d0*ex*yy)
c
      call obstra(1,1,me)
      go to 1995
c
c
c
c
c        1-tt, vector mesons with consequent ignoration of retardation
c        -------------------------------------------------------------
c
c
c
c
c        vector coupling (additional term)
c
c
c
c
  810 mc=0
c
      f(1)=eep1
      f(2)=xy
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-epe
      f(8)=f(7)
        e1=eme*eme
c
      do 815 mx=1,me
      ff=e1/c(4,ima(mx,mg,inter))
  815 call obstra(1,mx,mx)
c
c
c
c
c
c
c        tensor coupling and vector-tensor coupling
c
c
c
c
      go to 720
c
c
c
c
c        1-st, vector mesons in static limit
c        -----------------------------------
c
c
c
c        vector coupling
c
c
c
c
  900 mc=1
c
c
c
c        central term  '1'  only
c
c
c
      ff=1.d0
      f(1)=2.d0
      f(2)=0.d0
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-f(1)
      f(8)=f(7)
      f(9)=f(2)
      f(10)=f(2)
      f(11)=f(2)
      xxpyyh=xxpyy*0.5d0
c
c
c
c        central term  ' k**2 + p**2 '
c
c
c
      f(1)=f(1)+xxpyyh
      f(2)=f(2)+xy2
      f(3)=f(3)+xxpyyh
      f(6)=f(6)+xxpyyh
      f(7)=f(7)-xxpyyh
      f(8)=f(8)-xxpyyh
      f(10)=f(10)+xy2
      f(11)=f(11)-xy2
c
c
c
c        spin - spin
c
c
c
      f(1)=f(1)+xxpyyh*3.d0
      f(2)=f(2)-xy*3.d0
      f(3)=f(3)-xxpyyh
      f(6)=f(6)-xxpyyh
      f(7)=f(7)+xxpyyh
      f(8)=f(8)+xxpyyh
      f(10)=f(10)+xy
      f(11)=f(11)-xy
c
c
c
c        tensor  with  spin-spin-contribution
c
c
c
      f(1)=f(1)-xxpyyh
      f(2)=f(2)+xy
      f(3)=f(3)+xxpyyh
      f(4)=f(4)-xy
      f(5)=f(5)+xy
      f(6)=f(6)-xxpyyh
      f(7)=f(7)+(xx-yy)*0.5d0
      f(8)=f(8)-(xx-yy)*0.5d0
c
c
c
c        spin - orbit
c
c
c
      xy3=xy*3.d0
      f(4)=f(4)+xy3
      f(5)=f(5)+xy3
      f(10)=f(10)-xy3
      f(11)=f(11)+xy3
c
      call obstra(2,1,me)
c
c
c
c        case of additional sigma-l
c
c
c
c     mc=-1
c     xxyy8=xxyy/8.d0
c     f(1)=xxyy8
c     f(2)=0.d0
c     f(3)=f(1)
c     f(4)=f(2)
c     f(5)=f(2)
c     f(6)=-xxyy8
c     f(7)=f(1)
c     f(8)=f(7)
c     f(9)=f(6)*2.d0
c
c     call obstra(4,1,me)
c
c
c
c
c        tensor coupling
c
c      ( 1-ss , rho-meson in static limit, tensor coupling only )
c
c
 1000 mc=2
c
      if (mg.eq.10) mc=1
c
c
c        spin - spin
c
c
c
      f(1)=xxpyyh*3.d0
      f(2)=-xy*3.d0
      f(3)=-xxpyyh
      f(4)=0.d0
      f(5)=0.d0
      f(6)=-xxpyyh
      f(7)=xxpyyh
      f(8)=xxpyyh
      f(9)=0.d0
      f(10)=xy
      f(11)=-xy
c
c
c
c        tensor  with  spin-spin-contribution
c
c
c
      f(1)=f(1)-xxpyyh
      f(2)=f(2)+xy
      f(3)=f(3)+xxpyyh
      f(4)=f(4)-xy
      f(5)=f(5)+xy
      f(6)=f(6)-xxpyyh
      f(7)=f(7)+(xx-yy)*0.5d0
      f(8)=f(8)-(xx-yy)*0.5d0
c
      call obstra(2,1,me)
      if (mg.eq.10) go to 1995
c
c
c
c        case of additional sigma-l
c
c
c
c     f(1)=xxyy
c     f(2)=0.d0
c     f(3)=f(1)
c     f(4)=f(2)
c     f(5)=f(2)
c     f(6)=-xxyy
c     f(7)=f(1)
c     f(8)=f(7)
c     f(9)=f(6)*2.d0
c
c     call obstra(4,1,me)
c
c
c
c
c        vector-tensor coupling
c
c
c
c
      mc=3
c
c
c
c        central term  ' k**2 + p**2 '
c
c
c
      f(1)=-xxpyy
      f(2)=xy2
      f(3)=f(1)
      f(4)=0.d0
      f(5)=0.d0
      f(6)=f(1)
      f(7)=-f(1)
      f(8)=f(7)
      f(9)=0.d0
      f(10)=f(2)
      f(11)=-f(2)
c
c
c
c        spin - spin
c
c
c
      f(1)=f(1)+xxpyy*3.d0
      f(2)=f(2)-xy2*3.d0
      f(3)=f(3)-xxpyy
      f(6)=f(6)-xxpyy
      f(7)=f(7)+xxpyy
      f(8)=f(8)+xxpyy
      f(10)=f(10)+xy2
      f(11)=f(11)-xy2
c
c
c
c        tensor  with  spin-spin-contribution
c
c
c
      f(1)=f(1)-xxpyy
      f(2)=f(2)+xy2
      f(3)=f(3)+xxpyy
      f(4)=f(4)-xy2
      f(5)=f(5)+xy2
      f(6)=f(6)-xxpyy
      f(7)=f(7)+(xx-yy)
      f(8)=f(8)-(xx-yy)
c
c
c
c        spin - orbit
c
c
c
      xy4=xy*4.d0
      f(4)=f(4)+xy4
      f(5)=f(5)+xy4
      f(10)=f(10)-xy4
      f(11)=f(11)+xy4
c
      call obstra(2,1,me)
c
c
c
c        case of additional sigma-l
c
c
c
c     f(1)=xxyy
c     f(2)=0.d0
c     f(3)=f(1)
c     f(4)=f(2)
c     f(5)=f(2)
c     f(6)=-xxyy
c     f(7)=f(1)
c     f(8)=f(7)
c     f(9)=f(6)*2.d0
c
c     call obstra(4,1,me)
      go to 1995
c
c
c
c
c
c        1+  , a1-meson
c        --------------
c
c
c
c
 1100 mc=1
c
      ff=2.d0
      f(1)=eep1+ee
      f(2)=0.d0
      f(3)=-ee
      f(4)=-xy
      f(5)=xy2
      f(6)=-1.d0
      f(7)=ey
      f(8)=ex
c
      call obstra (1,1,me)
      go to 1995
c
c
c
c
c
c        2+  , f-meson
c        -------------
c
c
c
c
c        g1**2 coupling
c
c
c
c
 1200 mc=1
c
      ff=-8.d0
      ee4=ee*4.d0
      xxyy4=xxyy*4.d0
      eep143=eep1*fd3
        e1=2.d0*(xxpyy+td3)
      f(1)= eep143 +(3.d0*ee+2.d0)*xxpyy+xxyy4
      f(2)=(ee4+td3+xxpyy)*xy
      f(3)=3.d0*xxyy+eep1*e1
      f(4)=(2.d0*xxpyy+3.d0*eep1-d3)*xy
      f(5)=(ee4+11.d0*d3+3.d0*xxpyy)*xy
      f(6)= eep143 +(ee+2.d0)*xxpyy+xxyy4
        e2=-epe*e1
      f(7)=e2+xx*ex
      f(8)=e2+yy*ey
c        factors for additional terms
      f(9)=xxyy
      f(10)=ee*xy
      f(11)=xy
      f(12)=-ey*xy
      f(13)=-ex*xy
c
      call obstra(3,1,me)
      go to 1995
c
c
c
c
c        this has been the end of the contributions of mesons
c        ----------------------------------------------------
c
c
c
c
c        errors and warnings
c        -------------------
c
c
c
c
 9000 if (indmg(mg)) go to 1995
      write (kwrite,19000) mesong(mg)
19000 format(1h0////'0warning in obnna: meson-group  ',a4,'  does not ex
     1ist in this program.'/'0contribution ignored. execution continued.
     2'////)
      indmg(mg)=.true.
c
c
c
c
 1995 continue
c
c
c
c
c        add up contributions of mesons
c        ------------------------------
c
c
c
c
 2000 do 2005 il=iman,imen
      do 2005 iv=1,6
 2005 v(iv)=v(iv)+vj(iv,il)
c
c
c
c
      return
      end
      subroutine obpara
        use crdwrt
        use cstate
        use cob
c
c       obpara reads writes and stores the parameter for all ob-subrout.
c
c
      implicit double precision (a-h,o-z)
c
c
!      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
!      common /cstate/ j,heform,sing,trip,coup,endep,label
!      logical heform,sing,trip,coup,endep
c
!      character*4 label
c
c        common block for all ob-subroutines
c
!      common /coba/   vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
!     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
!     2                ez1,ez2,ct(96),wt(96),
!     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
!     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
!     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
!     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
!      logical indc,indxy,indpar
c
c
c        further specifications
c
      dimension cc(5)
      character*4 name(3),nname(15)
      integer imga(3)
      character*4 cut,cutg,end
      logical index
      logical zerocp,indcut
      character*4 mesong(12)
c
      data cut/'cut '/,cutg/'cutg'/,end/'end '/
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1                   '1-  ','1-t ','1-tt','1-st','1-ss',
     2                    '1+  ','2+  '/
      data index/.false./
      data zerocp/.true./,indcut/.false./
      data uf/197.3286d0/
c
c
c
c
c        statement functions needed for cray fortran
c
c
c      sqrt(x)=sqrt(x)
c      exp(x)=exp(x)
c      atan(x)=atan(x)
c      log(x)=log(x)
c
c
c
c
10000 format (2a4,a2,15a4)
10001 format (1h1)
10002 format (1h0/' jp  name      g**2      f/g       mass    iso-spin
     1iprop/spe'/9x,'cut typ     c u t - o f f   p a r a m e t e r s')
10003 format (2a4,a2,5f10.4)
10004 format (1h0,2a4,a2,2f10.4,f9.2,1x,2(f7.1,3x))
10005 format (1h ,2a4,a2,f3.1,f11.1,f9.4,f14.4,f10.4)
10006 format (2a4,a2,3i3)
10007 format (1h0,2a4,a2,3i3)
10008 format (1h ,61(1h-))
10011 format (1h1  //' obnna:  one-boson-exchange nn-nn interaction (num
     1er. integra.)')
10012 format (1h1  //' obnd:  one-boson-exchange nn-nd interaction (nume
     1r. integra.)')
10013 format (1h1  //' obdd:  one-boson-exchange nn-dd interaction (nume
     1r. integra.)')
10015 format ('0input-parameter-set:'/1h ,20(1h-))
10016 format (1h0,2a4,a2,15a4)
c
c
c
c
      if (index) go to 50
      index=.true.
c
      x=-1.d0
      y=-1.d0
c
c
c
c
c        maxima of certain indices related to the dimension as follows:
c        dimension c(mme,imee),ic(mice,imee),indc(mindce,imee),
c                  mgg(mge,3),mggo(mge,3),mesong(mge),vj(32,imee),
c                  ima(mee,mge,3)
c
      mge=12
      mee=5
      mme=10
      mice=10
      mindce=2
      imb=1
      ime=0
      imee=15
      imec=0
c        mme always ge mice, mindce
c
c        set all meson-parameters and indices to zero or .false.
c
      do 1 int=1,3
      imga(int)=0
      indpar(int)=.false.
      do 1 mgx=1,mge
      mgg(mgx,int)=0
    1 mggo(mgx,int)=0
c
c
      do 2 il=1,imee
      do 2 mm=1,mme
      if (mm.le.mindce) indc(mm,il)=.false.
      if (mm.le.mice) ic(mm,il)=0
    2 c(mm,il)=0.d0
      endep=.false.
c
c
c        indc(2,il) is reserved for information concerning the eikonal
c        form-factor
c
c
c
c
c
c
c        reading and writing of first 4,5 cards
c        --------------------------------------
c        --------------------------------------
c
c
c
c        write headline and read and write name of parameter set
c
   50 go to (51,52,53),inter
   51 write (kwrite,10011)
      go to 55
   52 write (kwrite,10012)
      go to 55
   53 write (kwrite,10013)
   55 write (kwrite,10008)
      write (kwrite,10015)
      read  (kread, 10000) name,nname
      write (kwrite,10016) name,nname
      if (inter.eq.1) label=name(1)
      indpar(inter)=.true.
c
c        read and write index-parameter concerning the factor of the
c        potential
c
      read  (kread, 10006) name,ift(inter)
      write (kwrite,10007) name,ift(inter)
      iftyp=ift(inter)
c**** if (iftyp.lt.0.or.iftyp.gt.4.or.inter.eq.1.and.iftyp.gt.2)goto9003
c
c        read and write parameters for numerical integration
c
      read  (kread, 10006) name,mint(inter),maxt(inter)
      write (kwrite,10007) name,mint(inter),maxt(inter)
c
c        read and write mass of nucleon
c
      read  (kread, 10003) name,wn
      write (kwrite,10004) name,wn
      wnq=wn*wn
      dwn=1.d0/wn
      dwnq=dwn*dwn
      wnn(inter)=wn
      if (inter.lt.2) go to 60
c
c        read and write mass of n-star
c
      read  (kread, 10003) name,wdd(inter)
      write (kwrite,10004) name,wdd(inter)
c
c        write headline for meson parameters
c
   60 write (kwrite,10002)
      write (kwrite,10008)
c
c
c
c
c        read, write and store meson parameters
c        --------------------------------------
c        --------------------------------------
c
c
c
   61 read  (kread, 10003) name,cc
c
c        check if data-card just read contains cut-off parameters
c
      if (name(1).eq.cut.or.name(1).eq.cutg) go to 70
c
c        check if end of mesons
c
      if (name(1).eq.end) go to 2000
c
c
c
c
c        write meson-parameters, which are no cut-off parameters
c        -------------------------------------------------------
c
c
c
c
      indcut=.false.
c
      write (kwrite,10004) name,cc
c
c        check if coupling constants are zero
c
      if (cc(1).ne.0.d0) go to 62
      zerocp=.true.
      go to 61
c
   62 zerocp=.false.
c
c        find out number of meson-group mg
c
      do 63 mg=1,mge
      if (name(1).eq.mesong(mg)) go to 64
   63 continue
      go to 9000
c
c
c
c
c        store meson parameters, which are no cut-off parameters
c        -------------------------------------------------------
c
c
c
c
   64 ime=ime+1
      if (ime.gt.imee) go to 9011
      mgg(mg,inter)=mgg(mg,inter)+1
      m=mgg(mg,inter)
      if (m.gt.mee) go to 9001
      ima(m,mg,inter)=ime
      if (m.ne.1) go to 65
      imga(inter)=imga(inter)+1
      mggo(imga(inter),inter)=mg
   65 continue
c
c        store coupling constant g**2
      c(1,ime)=cc(1)
c        store coupling constant f*g
      c(3,ime)=cc(1)*cc(2)
      if (inter.eq.2.and.mg.eq.10) c(1,ime)=c(1,ime)+c(3,ime)
c        store coupling constant f**2
      c(2,ime)=cc(2)*c(3,ime)
      if (inter.eq.1.and.mg.eq.10)
     1  c(1,ime)=c(1,ime)+c(3,ime)*2.d0+c(2,ime)
c        store meson mass sqare in units of nucleon mass square
      c(4,ime)=cc(3)*cc(3)*dwnq
c
c        test iso-spin
      icc=cc(4)
      if (icc.ne.0.and.icc.ne.1) go to 9004
c         store isospin as logical constant
      if (icc.eq.1) indc(1,ime)=.true.
c        store and test iprsp for meson/delta/nucleon
      icc=cc(5)
      iccc=mod(icc,100)
c        iprsp for nucleons
      ic(1,ime)=mod(iccc,10)
c        ispd for deltas
      ic(2,ime)=iabs(iccc/10)
c        ispm for mesons
      ic(3,ime)=iabs(icc/100)
      if (iabs(ic(1,ime)).gt.8) go to 9005
      if (iabs(ic(1,ime)).ge.2.and.iabs(ic(1,ime)).le.7) endep=.true.
c
c        index values for further storing
      mi=4
      mm=5
      go to 61
c
c
c
c
c        write cut-off parameters
c        ------------------------
c
c
c
c
   70 write (kwrite,10005) name,cc
c
c
c        check if individuel cut or general cut
c
      if (name(1).eq.cut) go to 73
c        case of general cut-off
      if (indcut) go to 90
      if (imec.ge.ime) go to 61
      imac=imec+1
      imec=ime
      if (imac.lt.imb) imac=imb
      go to 90
c        case of individuel cut-off
   73 imac=ime
      imec=ime
      if (zerocp) go to 61
c
c        save present values of indices
c
   90 indcut=.true.
      if (cc(1).eq.0.d0) go to 61
      mix=mi
      mmx=mm
c
c        start loop of mesons, which present cut-off refers to
c
      do 1095 im=imac,imec
      mi=mix
      mm=mmx
c
c
c
c
c        store cut-off parameters
c        ------------------------
c
c
c
c
c        store typ of cut-off
      ic(mi,im)=cc(1)
      ityp=ic(mi,im)
      if (ityp.lt.1.or.ityp.gt.9) go to 9002
c        store and test typ of propagator of cut-off
      ic(mi+1,im)=cc(2)
      if (ic(mi+1,im).lt.0.or.ic(mi+1,im).gt.8) go to 9006
      if (ic(mi+1,im).ge.2.and.ic(mi+1,im).le.7) endep=.true.
      go to (100,100,300,400,400,400,700,800,900),ityp
c
c
c        cut-off of dipole type
c        **********************
c
c
c        store and test exponent of cut-off
  100 ic(mi+2,im)=cc(3)
      if (ic(mi+2,im).lt.0) go to 9009
      if (ic(mi+2,im).gt.0) go to 101
c        exponent is zero, omit cut-off
      ic(mi,im)=0
      ic(mi+1,im)=0
      go to 1000
c        store cut-off mass for denominator
  101 c(mm+1,im)=cc(4)*cc(4)*dwnq
c        store numerator of cut-off
      c(mm,im)=c(mm+1,im)
      if (ityp.eq.2)     c(mm,im)=c(mm,im)-c(4,im)
      mi=mi+3
      mm=mm+2
      go to 1000
c
c
c        cut-off of regge type /schierholz/
c        **********************************
c
c
  300 c(mm,im)=2.d0/sqrt(4.d0-cc(3)*cc(3)*dwnq)
      c(mm+1,im)=cc(4)-1.d0
      c(mm+2,im)=cc(5)*wnq*1.d-6
      mi=mi+2
      mm=mm+3
      go to 1000
c
c
c        eikonal form factor
c        *******************
c
c
c        store gamma as -4.*gamma
  400 c(mm,im)=-cc(3)*4.d0
      ieik=ityp-3
c
c        compute and store normalization factor of t-form
      d=-c(4,im)
      d1=sqrt(-d)
      d2=sqrt(4.d0+d)
c     c(mm+1,im)=exp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*atan(d1/d2))
      c(mm+1,im)=exp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*atan2(d1,d2))
      print*,atan(d1/d2),atan2(d1,d2)
c
      go to (490,412,412),ieik
c
c        compute and store normalization factor of u-form
c        and t- * u-form
  412 ieik2=cc(4)
      if (ieik2.eq.0) ieik2=1
      ic(mi+2,im)=ieik2
      go to (421,422,423,421,422,423),ieik2
c
  421 ess=4.d0*(1.d0+cc(5)*dwn)**2
      go to 450
c
  422 endep=.true.
c
  423 go to 490
c
  450 d=-(4.d0-ess)
      if (ieik2.le.3) d=d+c(4,im)
      if (d.eq.0.d0) go to 460
      d2=sqrt(4.d0+d)
      if (d.lt.0.d0) go to 470
      d1=sqrt(d)
      c(mm+2,im)=exp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*log(0.5d0*(d1+d2)))
      go to 480
c
  460 c(mm+2,im)=exp(2.d0*cc(3))
      go to 480
c
  470 d1=sqrt(-d)
c     c(mm+2,im)=exp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*atan(d1/d2))
      c(mm+2,im)=exp(4.d0*cc(3)*(2.d0+d)/(d1*d2)*atan2(d1,d2))
      print*,atan(d1/d2),atan2(d1,d2)
c
  480 if (ieik.eq.3) c(mm+2,im)=c(mm+2,im)*c(mm+1,im)
c
  490 mi=mi+3
      mm=mm+3
      go to 1000
c
c
c        exponential form factor
c        ***********************
c
c
c        check exponent
  700 if (cc(3).lt.0.d0) go to 9009
      if (cc(3).gt.0.d0) go to 701
c        exponent is zero, omit cutoff
      ic (mi,im)=0
      ic (mi+1,im)=0
      go to 1000
c        compute constant factor for argument of exponential function
  701 c(mm,im)=cc(3)*wnq/(cc(4)*cc(4))
      mi=mi+2
      mm=mm+1
      go to 1000
c
c
c        cloudy bag form factor
c        ***********************
c
c
c        check exponent
  800 if (cc(3).lt.0.d0) go to 9009
      if (cc(3).gt.0.d0) go to 801
c        exponent is zero, omit cutoff
      ic (mi,im)=0
      ic (mi+1,im)=0
      go to 1000
c        store exponent
  801 ic(mi+2,im)=cc(3)
c        store cutoff radius
      c(mm,im)=cc(4)*wn/uf
      mi=mi+3
      mm=mm+1
      go to 1000
c
c
c        propagator of mass-distributed meson
c        ************************************
c
c
c        spin of meson
  900 icc=cc(3)
      ic(mi+2,im)=icc
c        full width
      c(mm,im)=cc(4)*dwn
c        2*pion-mass squared
      c(mm+1,im)=cc(5)*cc(5)*dwnq
c        recalculate width
      d=c(4,im)-c(mm+1,im)
      c(mm,im)=c(mm,im)*sqrt(c(4,im))/sqrt(d)
      if (icc.lt.1) go to 910
      do 905 i=1,icc
  905 c(mm,im)=c(mm,im)/d
  910 mi=mi+3
      mm=mm+2
      go to 1000
c
c
c
c
c        end cut-offs
c        ************
c
c        test dimensions
 1000 if (mi.gt.mice.or.mm-1.gt.mme) go to 9010
c
c
 1095 continue
      go to 61
c
c
c
c
c        last card
c        ---------
c        ---------
c
c
c
c
c        write end mesons
 2000 imaa(inter)=imb
      imea(inter)=ime
      imb=ime+1
      write (kwrite,10004) name
      write (kwrite,10008)
      write (kwrite,10008)
c
c
c
c
      return
c
c
c
c        errors
c        ------
c        ------
c
c
c
c
 9000 write (kwrite,19000) name(1)
19000 format (1h0/////'0error in obpara:  meson-group   ',a4,'   does
     1 not exist in this program.'/'0execution terminated.'////)
      go to 9999
c
c
 9001 write (kwrite,19001)
19001 format (1h0/////'0error in obpara: too many mesons within a meson-
     1group with respect to'/'0the given dimensions. execution terminate
     2.d'////)
      go to 9999
c
c
 9002 write (kwrite,19002) cc(1)
19002 format (1h0/////'0error in obpara: cut-off typ',f10.4,'  does not
     1exist in this program.'/'0execution terminated.'////)
      go to 9999
c
c
 9003 write (kwrite,19003) iftyp
19003 format (1h0/////'0error in obpara: factor typ has the nonp-permiss
     1ible value',i4,' .'/'0execution terminated.'////)
      go to 9999
c
c
 9004 write (kwrite,19004) cc(4)
19004 format (1h0/////'0error in obpara: isospin has the non-permissible
     1value',f10.4,'  .'/'0execution terminated.'////)
      go to 9999
c
c
 9005 write (kwrite,19005) cc(5)
19005 format (1h0/////'0error in obpara: iprop/spe has the non-permissib
     1le value',f10.4,'  .'/'0execution terminated.'////)
      go to 9999
c
c
 9006 write (kwrite,19006) cc(2)
19006 format (1h0/////'0error in obpara: the index for the propagator of
     1 the cut-off has the'/'0non-permissible value',f10.4,'  . executio
     2 n terminated.'////)
      go to 9999
c
c
 9009 write (kwrite,19009)
19009 format (1h0/////'0error in obpara: the exponent of the cut-off is
     1 less than zero.'/'0execution terminated.'////)
      go to 9999
c
c
 9010 write (kwrite,19010)
19010 format (1h0/////'0error in obpara: too many cut-off parameters wit
     1h respect to the given'/'0dimensions. execution terminated.'////)
      go to 9999
c
c
 9011 write (kwrite,19011)
19011 format (1h0/////'0error in obpara:  too many mesons with respect t
     1o the dimensions given'/'0to this program. execution terminated.'
     2////)
      go to 9999
c
c
 9999 stop
      end
      subroutine obstra (icase,max,mex)
        use crdwrt
        use cstate
        use cob
c
c        obstra computes the structure of one-boson-exchanges
c
c
      implicit double precision (a-h,o-z)
c
c
c        common blocks
c
!      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
!      common /cstate/ j,heform,sing,trip,coup,endep,label
!      logical heform,sing,trip,coup,endep
c
!      character*4 label
c
c        common block for all ob-subroutines
c
!      common /coba/   vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
!     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
!     2                ez1,ez2,ct(96),wt(96),
!     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
!     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
!     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
!     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
!      logical indc,indxy,indpar
c
c     further specifications
c
      dimension vv(32)
      dimension tt(2,3)
      logical index
      logical indiso
c
      data jj/-1/
      data index/.false./
c
c
c        statement functions needed for cray fortran
c
c      float(i)=float(i)
c      sqrt(x)=sqrt(x)
c
c
c
c
      if (index) go to 50
      index=.true.
c
c
      tt(1,1)=1.d0
      tt(2,1)=-3.d0
c
      do 1 ii=2,3
      do 1 i=1,2
    1 tt(i,ii)=1.d0
c
c
c
c
c
   50 do 1095 m=max,mex
      im=ima(m,mg,inter)
c
c
      if (mc.ne.1) go to 60
c
c
c
c
c        call integrals
c        --------------
c
c
c
c
      call obaia
c
c
c
c
   60 if (mc.lt.1) mc=1
c
      if (c(mc,im).eq.0.d0) go to 1095
c
c
c
c
      go to (100,200,300),inter
c
c
c
c
c        nn-nn helicity amplitudes /combinations/
c        ----------------------------------------
c
c
c
c
c        ground structure (a factor of 2 is included in v5 and v6)
c
c
  100 ive=6
c
      vv(1)=f(1)*ai(1,m)+f(2)*ai(2,m)
      vv(2)=f(3)*ai(1,m)+f(4)*ai(3,m)
      vv(3)=f(5)*ai(1,m)+f(6)*ai(2,m)
      vv(4)=f(4)*ai(1,m)+f(3)*ai(3,m)
      vv(5)=f(7)*ai(4,m)
      vv(6)=f(8)*ai(4,m)
c
c
      go to (1000,120,130,140),icase
c
c
c        additional terms  in case of tensor coupling
c
c
  120 vv(1)=vv(1)+f(9)*ai(5,m)
      vv(2)=vv(2)+f(10)*ai(2,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(10)*ai(5,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(10)*ai(6,m)
         e1=f(11)*ai(7,m)
      vv(5)=vv(5)+e1
      vv(6)=vv(6)+e1
      go to 1000
c
c
c        additional terms in case of 2+ mesons
c
c
  130 vv(2)=vv(2)+f(10)*ai(2,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(11)*ai(5,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(10)*ai(6,m)
      vv(5)=vv(5)+f(12)*ai(7,m)
      vv(6)=vv(6)+f(13)*ai(7,m)
      go to 1000
c
c
c        additional terms in case of sigma-l in static limit
c
c
  140 vv(1)=vv(1)+f(6)*ai(5,m)
      vv(2)=vv(2)+f(1)*ai(5,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(1)*ai(11,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(1)*ai(12,m)
      vv(5)=vv(5)+f(6)*ai(13,m)
      vv(6)=vv(6)+f(6)*ai(13,m)
      go to 1000
c
c
c
c
c        nn-nd helicity amplitudes
c        -------------------------
c
c
c
c
  200 ive=16
c
      ai3m1=ai(3,m)-ai(1,m)
      ai6m2=ai(6,m)-ai(2,m)
      ai3p1=ai(3,m)+ai(1,m)
      ai6p2=ai(6,m)+ai(2,m)
c
c
      vv( 1)= f( 1)* ai( 4,m) + f( 2)* ai( 7,m)
      vv( 2)= f( 3)* ai( 1,m) + f( 4)* ai( 2,m) + f( 5)* ai( 5,m)
      vv( 3)= f( 6)* ai( 4,m) + f( 7)* ai( 7,m)
      vv( 4)= f( 8)* ai( 8,m)
      vv( 5)= f( 9)* ai( 8,m)
      vv( 6)=-f(10)* ai( 4,m) + f(11)* ai( 7,m)
      vv( 7)= f(12)* ai( 1,m) - f(13)* ai( 2,m) + f(14)* ai( 5,m)
      vv( 8)=-f(15)* ai( 4,m) + f(16)* ai( 7,m)
c
      vv( 9)= f(17)* ai3m1    + f(18)* ai6m2
      vv(10)= f(19)* ai( 4,m) + f(20)* ai( 7,m)
      vv(11)= f(21)* ai3p1    + f(22)* ai6p2
      vv(12)= f(23)* ai( 9,m)
      vv(13)=-f(24)* ai(10,m)
      vv(14)=-f(25)* ai3m1    + f(26)* ai6m2
      vv(15)=-f(27)* ai( 4,m) + f(28)* ai( 7,m)
      vv(16)=-f(29)* ai3p1    + f(30)* ai6p2
      go to 1000
c
c
c
c        nn-dd helicity amplitudes
c        -------------------------
c
c
c
c
  300 ive=32
c
      ai31p=ai( 3,m)+ai( 1,m)
      ai31m=ai( 3,m)-ai( 1,m)
      ai62p=ai( 6,m)+ai( 2,m)
      ai62m=ai( 6,m)-ai( 2,m)
      aic5p=ai(12,m)+ai( 5,m)
      aic5m=ai(12,m)-ai( 5,m)
c
c
      vv( 1)= f( 1)*(ai(1,m)+ai(2,m)) + f( 2)*(ai(2,m)+ai(5,m)) +
     1        f( 3)*(ai(5,m)+ai(11,m))
c
      vv( 2)= f( 4)* ai( 4,m) + f( 5)* ai( 7,m) + f( 6)* ai(13,m)
      vv( 3)= f( 7)* ai( 8,m) + f( 8)* ai(14,m)
      vv( 4)= f( 9)* ai(17,m)
      vv( 5)=vv( 2)
      vv( 6)= f(10)* ai( 1,m) + f(11)* ai( 2,m) + f(12)* ai( 5,m) +
     1        f(13)* ai(11,m)
      vv( 7)= f(14)* ai( 4,m) + f(15)* ai( 7,m) + f(16)* ai(13,m)
      vv( 8)=-f(17)* ai( 8,m) + f(18)* ai(14,m)
      vv( 9)=vv( 3)
      vv(10)=vv( 7)
      vv(11)=-f(19)* ai( 1,m) + f(20)* ai( 2,m) - f(21)* ai( 5,m) +
     1        f(22)* ai(11,m)
      vv(12)= f(23)* ai( 4,m) - f(24)* ai( 7,m) + f(25)* ai(13,m)
      vv(13)=vv( 4)
      vv(14)=vv( 8)
      vv(15)=vv(12)
c
      vv(16)=-f(26)*(ai(1,m)-ai(2,m)) + f(27)*(ai(2,m)-ai(5,m)) -
     1        f(28)*(ai(5,m)-ai(11,m))
c
c
      vv(17)= f(29)* ai( 4,m) + f(30)* ai( 7,m) + f(31)* ai(13,m)
      vv(18)= f(32)* ai31p    + f(33)* ai62p    + f(34)* aic5p
      vv(19)= f(35)* ai( 9,m) + f(36)* ai(15,m)
      vv(20)= f(37)* ai(18,m)
      vv(21)= f(44)* ai31m    - f(45)* ai62m    + f(46)* aic5m
      vv(22)= f(38)* ai( 4,m) + f(39)* ai( 7,m) + f(40)* ai(13,m)
      vv(23)= f(41)* ai31p    + f(42)* ai62p    + f(43)* aic5p
      vv(24)=vv(19)
      vv(25)= f(47)* ai(10,m) - f(48)* ai(16,m)
      vv(26)= f(49)* ai31m    - f(50)* ai62m    + f(51)* aic5m
      vv(27)=vv(22)
      vv(28)=vv(18)
      vv(29)= f(52)* ai(19,m)
      vv(30)=vv(25)
      vv(31)=vv(21)
      vv(32)=vv(17)
c
c
c
c
 1000 if (inter.ge.2) go to 1040
c
c
c
c
c        set certain cases to zero in case of inter=1
c
      if (j.ne.0) go to 1021
      vv(2)=0.d0
      vv(4)=0.d0
      vv(5)=0.d0
      vv(6)=0.d0
c
 1021 if (.not.sing) vv(1)=0.d0
      if (.not.trip) vv(2)=0.d0
      if (coup) go to 1030
      do 1025 iv=3,6
 1025 vv(iv)=0.d0
c
 1030 if (heform) go to 1040
c
c
c        transformation into lsj-formalism in case of inter=1
c        (if requested)
      if (j.eq.jj) go to 1035
      jj=j
      aj=float(j)
      aj1=float(j+1)
      d2j1=1.d0/float(2*j+1)
      arjj1=sqrt(aj*aj1)
c
 1035 v3=vv(3)
      v4=vv(4)
      v5=vv(5)
      v6=vv(6)
      v34=-arjj1*(v3-v4)
      v56=arjj1*(v5+v6)
      vv(3)=d2j1*(aj1*v3+aj*v4-v56)
      vv(4)=d2j1*(aj*v3+aj1*v4+v56)
      vv(5)=d2j1*(v34-aj1*v5+aj*v6)
      vv(6)=d2j1*(v34+aj*v5-aj1*v6)
c
c
c        possible different sign depending on the convention used
      vv(5)=-vv(5)
      vv(6)=-vv(6)
c
c
c
c
c        multiply with factors
c        ---------------------
c
c
c
c
 1040 is=mod(j,2)+1
      it=mod(is,2)+1
      indiso=indc(1,im)
      cmc=c(mc,im)
      fc=fff*ff*cmc
      do 1045 iv=1,ive
c
c        multiply with coupling-constant and factors fff and ff
c
      vv(iv)=vv(iv)*fc
      if (inter.ge.2) go to 1045
c
c        multiply with isospin factor
c
      if (.not.indiso) go to 1045
      if (iv.eq.2) go to 1043
      vv(iv)=vv(iv)*tt(is,inter)
      go to 1045
 1043 vv(iv)=vv(iv)*tt(it,inter)
c
c     add up in case of several couplings for one meson-exchange
c     and store
 1045 vj(iv,im)=vj(iv,im)+vv(iv)
c
c
c        if single contributions to one meson are to be printed:
c     write (kwrite,10001) mesong(mg),mc
c     write (kwrite,10002) (vv(iv),iv=1,ive)
10001 format (' contribution from ',a4,' mc =',i2)
10002 format (33x,6d16.6)
c
c
 1095 continue
c
c
      return
      end
      subroutine obaia
        use cpot
        use cstate
        use cpoted
        use cob
c
c        obaia integrates over theta
c
c
      implicit double precision (a-h,o-z)
c
!      common /cpot/   v(6),xmev,ymev
!      common /cstate/ j,heform,sing,trip,coup,endep,label
!      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,
!     1                znrl,zrel,smev,noced
!      logical heform,sing,trip,coup,endep
!      logical noced
c
!      character*4 label
c
c        common block for all ob-subroutines
c
!      common /coba/   vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
!     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
!     2                ez1,ez2,ct(96),wt(96),
!     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
!     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
!     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
!     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
!      logical indc,indxy,indpar
c
c
c        further specifications
      dimension gi(7)
      dimension pj(7,96)
      logical indj,indint,indee,indepe,indz,indepz
c
      data ige/7/
c
c      real*4 aez,axy2,aomq,am1,am2,am
      data nnt/-1/,iinter/-1/,jj/-1/
c
c
c        statement functions needed for cray fortran
c
c      float(i)=float(i)
c      sqrt(x)=sqrt(x)
c
c
c
c
      if (inter.eq.iinter) go to 60
      iinter=inter
      indint=.false.
      min=mint(inter)
      max=maxt(inter)
c
      go to (51,51,53),inter
   51 igeint=5
      go to 55
   53 igeint=7
   55 continue
c
      wn=wnn(inter)
      dwn=1.d0/wn
      wnq=wn*wn
c        wd is the mass of the delta
c        wdd(..) is the array of the masses of the delta
      wd=wdd(inter)
      wdq=wd*wd
c
c
c
c
   60 if (j.eq.jj) go to 70
      jj=j
      indj=.false.
c
c
      aj=float(j)
      aj1=float(j+1)
      dj1=1.d0/aj1
      ajdj1=aj*dj1
      aaj=sqrt(ajdj1)
c
c
      aj2=float(j+2)
      aj3=float(j+3)
      ajm1=float(j-1)
      ajm2=float(j-2)
      ajm3=float(j-3)
c
c
      ajj1=aj*aj1
      ajj2=ajm1*aj2
      ajj3=ajm2*aj3
      ajja=aj*ajm3
      ajjb=aj*ajm1
c
      aajj=0.d0
      if (j.gt.1)
     1aajj=aj/sqrt(ajj1*ajj2)
c
      aaj1=aajj*ajm1
      aaj2=aajj*aj1
      aaj3=aajj*2.d0
c
      if (j.gt.1) go to 62
      aajj=0.d0
      go to 63
   62 aajj=1.d0/(aj1*sqrt(ajj2))
c
   63 aaj4=aajj*ajjb
      aaj5=aajj*aj1*2.d0
      aaj6=aajj*(ajj1+2.d0)
      aaj7=aajj*ajj2
c
      if (j.gt.2) go to 64
      aajj=0.d0
      go to 65
   64 aajj=-aj/sqrt(ajj1*ajj2*ajj3)
c
   65 aaj8=aajj*(ajj1+6.d0)
      aaj9=aajj*ajj2
      aaj10=aajj*(ajja+2.d0)
      aaj11=aajj*(ajja-6.d0)
c
      if (j.gt.2) go to 66
      aajj=0.d0
      go to 67
   66 aajj=-1.d0/(aj1*sqrt(ajj2*ajj3))
c
   67 aaj12=aajj*ajjb*ajm2
      aaj13=aajj*(aj*ajjb+4.d0*aj+12.d0)
      aaj14=aajj*(5.d0*ajj1+6.d0)
      aaj15=aajj*3.d0*ajj2
      aaj16=aajj*ajj3*ajm1
      aaj17=aajj*aj1*ajj3
      aaj18=aajj*2.d0*ajj3
c
c
c
c
c
c        find out appropriate number of gauss-points
c        -------------------------------------------
c
c
c
c
   70 c4=c(4,im)
      iprsp=ic(1,im)
c
c
c        prepare starting energy
      if (noced) go to 73
      noced=.true.
      indz=.false.
      indepz=.false.
      if (iprsp.lt.2) go to 74
      if (iprsp.ge.4) go to 72
   71 indz=.true.
      z=2.d0*sqrt(wnq+q0qmev)
      go to 74
   72 if (iprsp.eq.8) go to 74
      indepz=.true.
      eppq=pmev*pmev
      epz=zrel
      go to 74
c
c
   73 if (indxy.and.indint) go to 85
   74 indint=.true.
      indee=.false.
      indepe=.false.
c
      ispm=ic(3,im)
      ez1=0.d0
      if (ispm.eq.1.and.iprsp.lt.2) go to 90
c
c         delta-nucleon mass difference in propagator
c
      ez1=(wd-wn)*dwn
c
c
      if (iprsp.lt.2) go to 90
c
c        prepare for propagators of non-covariant perturbation theory
   75 if (iprsp.ge.4) go to 76
      if (.not.indz) go to 71
      pmevq=0.d0
      ez=z
      go to 77
   76 if (iprsp.eq.8) go to 1800
      if (.not.indepz) go to 72
      pmevq=eppq
      ez=epz
   77 qfmevq=qfmev*qfmev
      qqx=xmev*xmev+pmevq
      ymevq=ymev*ymev
      qqy=ymevq+pmevq
      go to (1100,1200,1300),inter
c
c        ez1 for the nn case
c
 1100 ez1=spe(qqx,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,iprsp,0)
     1  + spe(qqy,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,iprsp,0)
     2  - ez
      ez1=ez1*dwn
      go to 80
c
c        ez1 for the nd case
c
 1200 pq=4.d0*pmevq/(wn+wd)**2
      ispd=ic(2,im)
      ed=spe(ymevq+pq*wdq,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wd,2,
     1 ispd,0)
      en=spe(ymevq+pq*wnq,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,
     1 iprsp,0)
      eza=spe(qqx,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,iprsp,0)
     1 - ez
      ez1=eza+en
      ez2=eza+ed
      ez1=ez1*dwn
      ez2=ez2*dwn
      go to 80
c
c        ez1 for the dd case
c
 1300 ispdd=ic(2,im)
      ez1=spe(qqx,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wn,2,iprsp,0)
     1  + spe(qqy,qfmevq,uanspp,wsnspp,ucnspp,udnspp,wd,2,ispdd,0)
     2  - ez
      ez1=ez1*dwn
      go to 80
c
c
c        case iprsp=8
c
 1800 ez1=ex-ey
c
c
c        store ez1 and ez2
c
c
   80 if (iprsp.ge.3) go to 81
      indee=.true.
      ee1=ez1
      ee2=ez2
      go to 89
   81 iiprsp=iprsp
      indepe=.true.
      epe1=ez1
      epe2=ez2
      go to 89
c
c
c        get stored ez1 and ez2
c
c
   85 if (iprsp.lt.2) go to 90
      if (iprsp.ge.3) go to 86
      if (.not.indee) go to 75
      ez1=ee1
      ez2=ee2
      go to 89
   86 if (iprsp.ne.iiprsp) go to 75
      if (.not.indepe) go to 75
      ez1=epe1
      ez2=epe2
c
c
   89 aez=ez1
c
c
c        compute am
c
c
   90 axy2=xy2
      if (iprsp.ne.1) go to 91
      aomq=eem12+c4
      go to 92
   91 aomq=xxpyy+c4
c
   92 am=axy2/aomq
c
      if (iprsp.lt.2) go to 93
      am1=am
      am2=axy2/(aomq+aez*abs(aez))
      if (am2.lt.0.) go to 94
      am=dmax1(am1,am2)
c     am=amax1(am1,am2)     
c
c
c        compute number of gausspoints (nt)
c
c
   93 if (am.gt.0.999) go to 94
c
c
      if (am.gt.0.85) am=am**(-log(1.-am)-0.9)
c
c
      nt=float(min)/(1.-am)+0.9
c
c
      if (nt.gt.max) nt=max
      go to 95
c
c
   94 nt=max
c
c
   95 nt=nt+j
c
c        compute nt, which is suitable for gset
c
      if (nt.le.16) go to 98
      if (nt.gt.24) go to 96
      nt=4*(nt/4)
      go to 98
   96 if (nt.gt.48) go to 97
      nt=8*(nt/8)
      go to 98
   97 nt=16*(nt/16)
      if (nt.gt.96) nt=96
c
   98 if (nt.eq.nnt.and.indj) go to 100
c
c
c
c
c        call gauss-points
c        -----------------
c
c
c
c
      call gset (-1.d0,1.d0,nt,ct,wt)
      nnt=nt
c
c
c
c
c        call legendre-polynoms if necessary
c        -----------------------------------
c
c
c
c
      indxy=.false.
      indj=.true.
      do 99 i=1,nt
      t=ct(i)
      call legp (pj(1,i),pj(3,i),t,j)
      pj(2,i)=pj(1,i)*t
      pj(4,i)=pj(2,i)*t
      pj(6,i)=pj(4,i)*t
      pj(5,i)=pj(3,i)*t
   99 pj(7,i)=pj(5,i)*t
c
c
c
c
c        call integrand
c        --------------
c
c
c
c
  100 call obaaa
c
c
c
c
c        prepare for integration
c
c
c
c
      do 2001 ig=1,igeint
 2001 gi(ig)=0.d0
c
c
c
c
c        integration-loop of theta
c        -------------------------
c
c
c
c
      do 2005 i=1,nt
      do 2005 ig=1,igeint
 2005 gi(ig)=gi(ig)+pj(ig,i)*aa(i)
c
c
c
      if (j.ne.0) go to 2010
      gi(3)=0.d0
      gi(5)=0.d0
      gi(7)=0.d0
c
c
c
c
c        combinations of integrals
c        -------------------------
c
c
c
c
 2010 ai(1,m)=gi(1)
c
      ai(2,m)=gi(2)
      ai(3,m)= ajdj1*gi(2)+dj1*gi(3)
      gi23m  =gi(2)-gi(3)
      ai(4,m)=aaj*gi23m
c
c
      ai(5,m)=gi(4)
      ai(6,m)= ajdj1*gi(4)+dj1*gi(5)
      gi45m  =gi(4)-gi(5)
      ai(7,m)=aaj*gi45m
c
c
      if (inter.eq.1) go to 3000
c
c
      ai( 8,m)= aaj1*gi(4)-aaj2*gi(1)+aaj3*gi(5)
      aai1    = aaj4*gi(4)+aaj5*gi(1)-aaj6*gi(5)
      aai2    = aaj7*gi23m
      ai( 9,m)= aai2+aai1
      ai(10,m)= aai2-aai1
c
c
      if (inter.ne.3) go to 3000
c
c
      ai(11,m)=gi(6)
      ai(12,m)=ajdj1*gi(6)+dj1*gi(7)
      ai(13,m)=aaj*(gi(6)-gi(7))
c
      ai(14,m)= aaj1*gi(6)-aaj2*gi(2)+aaj3*gi(7)
      aai1    = aaj4*gi(6)+aaj5*gi(2)-aaj6*gi(7)
      aai2    = aaj7*gi45m
      ai(15,m)= aai2+aai1
      ai(16,m)= aai2-aai1
c
      ai(17,m)= aaj8*gi(7)-aaj9*gi(3)-aaj10*gi(6)+aaj11*gi(2)
      aai1    =-aaj12*gi(6)+aaj13*gi(2)-aaj14*gi(7)+aaj15*gi(3)
      aai2    = aaj16*gi(4)-aaj17*gi(1)+aaj18*gi(5)
      ai(18,m)= aai1-aai2
      ai(19,m)= aai1+aai2
c
c
c
 3000 return
      end
      subroutine obaaa
        use cstate
        use cpoted
        use cob
        
c
c        obaaa computes the propagators and the cutoffs of ob-exchanges
c
c
      implicit double precision (a-h,o-z)
c
!      common /cstate/ j,heform,sing,trip,coup,endep,label
!      common /cpoted/ q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,
!     1                znrl,zrel,smev,noced
!      logical heform,sing,trip,coup,endep
!      logical noced
c
!      character*4 label
c
c        common block for all ob-subroutines
c
!      common /coba/   vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
!     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
!     2                ez1,ez2,ct(96),wt(96),
!     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
!     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
!     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
!     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
!      logical indc,indxy,indpar
c
c
c
c        further specifications
      dimension deltaq(96,3),cut(96)
      logical indp2
      data iinter/-1/
      data ssmev/-1.d0/
c
c
c        statement functions needed for cray fortran
c
c      sqrt(x)=sqrt(x)
c      abs(x)=abs(x)
c      log(x)=log(x)
c      exp(x)=exp(x)
c      sin(x)=sin(x)
c      cos(x)=cos(x)
c      atan(x)=atan(x)
c
c
c
c
      if (inter.eq.iinter) go to 60
      iinter=inter
c
c        dwn is needed for the eikonal cutoff
c
      dwn=1.d0/wnn(inter)
c
c
c
c
c        delta square
c        ------------
c
c
c
c
   60 if (indxy) go to 1000
      indxy=.true.
      do 65 i=1,nt
      xy2t=xy2*ct(i)
c
c        retardation ignored
c
      deltaq(i,1)=xy2t-xxpyy
c     ----------------------
c
c
c        retardation incorporated
c
      deltaq(i,2)=xy2t-eem12
c     ----------------------
c
c
c        for on shell cutoff  ( ferchlaender )
c
   65 deltaq(i,3)=-xy2t-xxpyy
c     -----------------------
c
c
c
c        propagator
c        ----------
c        ----------
c
c
c
c
 1000 c4=c(4,im)
      iprsp=ic(1,im)
      if (iprsp.lt.0) go to 1400
      if (iprsp.ge.2) go to 1050
      iret=iprsp+1
c
      go to (1010,1020,1030), inter
c
c         propagator for the nn case
 1010 do 1011 i=1,nt
 1011 aa(i)=wt(i)/(c4-deltaq(i,iret))
      go to 1500
c         propagator for the nd case
 1020 do 1021 i=1,nt
      omq=c4-deltaq(i,iret)
      om=sqrt(omq)
 1021 aa(i)=wt(i)*(1.d0/omq+1.d0/(om*(om+ez1)))*0.5d0
      go to 1500
c         propagator for the dd case
 1030 do 1031 i=1,nt
      omq=c4-deltaq(i,iret)
      om=sqrt(omq)
 1031 aa(i)=wt(i)/(om*(om+ez1))
      go to 1500
c
c
c        starting energy dependent propagator
c
c
 1050 ispm=ic(3,im)
      go to (1100,1200,1300),inter
c
c        the propagator for the nn case
c
 1100 do 1105 i=1,nt
      omq=c4-deltaq(i,1)
      om=sqrt(omq)
      oms=om
      if (ispm.le.2) go to 1105
      if (abs(ez1).lt.1.d-12) go to 1105
      oms=oms+smp(-deltaq(i,1),qfmev,ispm)
 1105 aa(i)=wt(i)/(om*(oms+ez1))
      go to 1500
c
c        the propagator for the nd case
c
 1200 do 1205 i=1,nt
      omq=c4-deltaq(i,1)
      om=sqrt(omq)
      oms=om
      if (ispm.le.2) go to 1205
      oms=oms+smp(-deltaq(i,1),qfmev,ispm)
 1205 aa(i)=wt(i)*(1.d0/(om*(oms+ez1))+1.d0/(om*(oms+ez2)))*0.5d0
      go to 1500
c
c        the propagator for the dd case
c
 1300 do 1305 i=1,nt
      omq=c4-deltaq(i,1)
      om=sqrt(omq)
      oms=om
      if (ispm.le.2) go to 1305
      oms=oms+smp(-deltaq(i,1),qfmev,ispm)
 1305 aa(i)=wt(i)/(om*(oms+ez1))
      go to 1500
c
c
c        "no propagator"
c
 1400 do 1405 i=1,nt
 1405 aa(i)=wt(i)
c
c
 1500 continue
c
c
c
c
c
c        cut-offs
c        --------
c        --------
c
c
c
c
      mi=4
      mm=5
c
c
  999 ityp=ic(mi,im)
      if (ityp.eq.0) go to 2000
      iprspc=ic(mi+1,im)
      iret=iprspc+1
      if (iprspc.eq.3) iret=1
      go to (100,100,300,400,400,400,700,800,900),ityp
c
c
c
c
c        cut-off of dipole type
c        **********************
c
c
c
c
  100 c5=c(mm,im)
      c6=c(mm+1,im)
      nexp=ic(mi+2,im)
c
      do 105 i=1,nt
c
      aaa=c5/(c6-deltaq(i,iret))
c     -------------------------
c
      do 105 ii=1,nexp
  105 aa(i)=aa(i)*aaa
c
      if (iprspc.le.2) go to 120
c
      do 110 i=1,nt
c
      aaa=c5/(c6-deltaq(i,iprspc))
c     ----------------------------
c
      do 110 ii=1,nexp
  110 aa(i)=aa(i)*aaa
c
c
  120 mi=mi+3
      mm=mm+2
      go to 999
c
c
c
c
c        cut-off of regge type /schierholz/
c        **********************************
c
c
c
c
  300 ax=ex*c(mm,im)
      ay=ey*c(mm,im)
      expo=log((ax+sqrt(ax*ax-1.d0))*(ay+sqrt(ay*ay-1.d0)))
      expo1=c(mm+1,im)*expo
      expo2=c(mm+2,im)*expo
      do 305 i=1,nt
      expon=expo1+expo2*deltaq(i,iret)
      if (expon.lt.-50.d0) go to 302
c
      aa(i)=aa(i)*exp(expon)
c     ---------------------
c
c
c     ---------------------
c
      go to 305
  302 aa(i)=0.d0
  305 continue
      mi=mi+2
      mm=mm+3
      go to 999
c
c
c
c
c        eikonal form factor
c        *******************
c
c
c
c
  400 ieik=ityp-3
c
      eikc5=c(mm,im)
      do 407 i=1,nt
      expon=0.d0
      ieik3=0
      go to (401,402,403),ieik
c
c        t-form
  401 d=-deltaq(i,iret)
      go to 404
c
c        u-form
  402 d=2.d0*xy2*ct(i)-deltaq(i,iret)
      go to 404
c
c        t-form * u-form
  403 ieik3=ieik3+1
      go to (401,402,405),ieik3
c
  404 d1=sqrt(d)
      d2=sqrt(4.d0+d)
c
      expon=eikc5*(2.d0+d)/(d1*d2)*log(0.5d0*(d1+d2))+expon
c     ------------------------------------------------
c
      go to (405,405,403),ieik
  405 if (expon.lt.-50.d0) go to 406
c
      cut(i)=exp(expon)
c     ------------------
c
      go to 407
  406 cut(i)=0.d0
  407 continue
c
c        get or calculate normalization factor
c
      go to (411,412,412),ieik
c
c        normalization of t-form
  411 c6=c(mm+1,im)
      go to 490
c
c        normalization of u-form and t- * u-form
  412 ieik2=ic(mi+2,im)
      go to (421,422,423,421,422,423),ieik2
c
  421 c6=c(mm+2,im)
      go to 490
c
  422 if (smev.eq.ssmev) go to 442
      ssmev=smev
      do 441 il=1,ime
  441 indc(2,il)=.false.
  442 if (indc(2,im)) go to 421
      indc(2,im)=.true.
      ess=smev*dwn
      ess=ess*ess
      go to 450
c
  423 ess=4.d0*ex*ey
c
  450 d=-(4.d0-ess)
      if (ieik2.le.3) d=d+c(4,im)
      if (d.eq.0.d0) go to 460
      d2=sqrt(4.d0+d)
      if (d.lt.0.d0) go to 470
      d1=sqrt(d)
      c6=exp(-eikc5*(2.d0+d)/(d1*d2)*log(0.5d0*(d1+d2)))
      go to 480
c
  460 c6=exp(-eikc5*0.5d0)
      go to 480
c
  470 d1=sqrt(-d)
c     c6=exp(-eikc5*(2.d0+d)/(d1*d2)*atan(d1/d2))
      c6=exp(-eikc5*(2.d0+d)/(d1*d2)*atan2(d1,d2))
      print*,atan(d1/d2),atan2(d1,d2)
c
  480 if (ieik.eq.3) c6=c6*c(mm+1,im)
      if (ieik2.eq.3.or.ieik2.eq.6) go to 490
      c(mm+2,im)=c6
c
c        compute form factor
c
  490 do 495 i=1,nt
  495 aa(i)=aa(i)*cut(i)*c6
c     -------------------
c
      mi=mi+3
      mm=mm+3
      go to 999
c
c
c        exponential form factor
c        ***********************
c
c
  700 c5=c(mm,im)
      do 705 i=1,nt
c
      expo=deltaq(i,iret)*c5
c     ----------------------
c
      if (expo.lt.-50.d0) go to 704
c
      aa(i)=aa(i)*exp(expo)
c     ----------------------
c
      go to 705
  704 aa(i)=0.d0
  705 continue
      mi=mi+2
      mm=mm+1
      go to 999
c
c
c        cloudy bag form factor
c        ***********************
c
c
  800 c5=c(mm,im)
      nexp=ic(mi+2,im)
      do 805 i=1,nt
c
      arg=sqrt(-deltaq(i,iret))*c5
      argc=arg*arg*arg
      aaa=3.d0*(sin(arg)-arg*cos(arg))/argc
c
      do 805 ii=1,nexp
  805 aa(i)=aa(i)*aaa
c
      mi=mi+3
      mm=mm+1
      go to 999
c
c
c        propagator of mass-distributed meson
c        ************************************
c
c
  900 c5=c(mm,im)
      c6=c(mm+1,im)
      nspin=ic(mi+2,im)
      indp2=.false.
      if (iprspc.le.1) go to 901
      indp2=.true.
      iret=1
c
  901 do 915 i=1,nt
      d=c6-deltaq(i,iret)
      if (nspin.eq.0) go to 903
      d1=-d*deltaq(i,iret)/c4
  903 d=sqrt(d)
      if (nspin.eq.0) go to 907
      do 905 ii=1,nspin
  905 d=d*d1
c
  907 omq=c4-deltaq(i,iret)+c5*d
      if (indp2) go to 910
      aa(i)=aa(i)/omq
      go to 915
c
  910 om=sqrt(omq)
      aa(i)=aa(i)/(om*(om+ez1))
c
  915 continue
c
      mi=mi+3
      mm=mm+2
      go to 999
c
c
c
c
 2000 return
      end

