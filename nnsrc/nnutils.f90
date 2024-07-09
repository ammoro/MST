      subroutine setconstants()
        use constants
        pi = 4d0*atan(1d0)
        amu=931.49432
        hbarc=197.3289
        finec=137.03599d0
      end subroutine setconstants


********************************************************************
* Subroutines shared by different subprograms
* f90 version A.Moro, R. Crespo (2001)
*******************************************************************
      SUBROUTINE LOGFAC(L)
      use factorials
      implicit none
      integer :: i,j,l
c***********************************************************
C  FLOG(I)= LN((I-1)!)
      write(99,*)'logfac: calculating factorials for l<=',l
      if (allocated(fact)) deallocate(fact)
      if (allocated(flog)) deallocate(flog)
      allocate(fact(0:l))
      allocate(flog(l))
C  THE VALUE OF L CORRESPONDS TO THE DIMENSION OF (FLOG(
      flog(1)=0d0
      fact(0)=1d0
      fact(1)=1d0

      if (l> 200) then
         write(99,*)'************************************************'
         write(99,*)'I will calculate factorials for l<=200, not',l
         write(99,*)'*************************************************'
      endif


      DO 100 J=2,L
         if (j<200) fact(j)=fact(j-1)*dble(j)
!          if (j<200) write(40,*)j,fact(j)
  100 FLOG(J)=FLOG(J-1)+LOG(J-1D0)
      RETURN
      END

c****************************************************************
!!$      subroutine fctrls
!!$        use factorials
!!$      implicit double precision (a-h,o-z),integer*4(i-n)
!!$      faclog(1)=0.0
!!$      faclog(2)=0.0
!!$      f(1)=1.0
!!$      f(2)=1.0
!!$      fn=1.0
!!$      do 50 n=3,500
!!$      fn=fn+1.0
!!$      faclog(n)=faclog(n-1)+dlog(fn)
!!$      if(n.gt.52)go to 50
!!$      f(n)=f(n-1)*fn
!!$ 50   continue
!!$      return
!!$      end
c******************************************************************
      function threj0(ia,ib,ic)
        use factorials
      implicit real*8 (a-h,o-z),integer*4(i-n)
      threj0=0.0
      igtw=(ia+ib+ic)/2
      if(mod(igtw,2).ne.0) go to 1000
      ig=igtw/2
      iahf=ia/2
      ibhf=ib/2
      ichf=ic/2
      s1=1-2*mod(ig+ichf,2)
      s1=s1*(1-2*mod(iabs(iahf-ibhf),2))
      iabc=igtw+2
      iabmc=iahf+ibhf-ichf+1
      icamb=ichf+iahf-ibhf+1
      ibcma=ibhf+ichf-iahf+1
      igma=ig-iahf+1
      igmb=ig-ibhf+1
      igmc=ig-ichf+1
!!$      if (.not.allocated(flog)) then
!!$         write(*,*)'threj0:flog not allocated!'
!!$      endif
!!$      write(*,*)'flog(1,2,3)',flog(1),flog(2),flog(3)
!!$      stop
      r1=0.5*(flog(iabmc)+flog(icamb)+flog(ibcma)-flog(iabc))
     1        +flog(ig+1)-flog(igma)-flog(igmb)-flog(igmc)
      threj0=s1*exp(r1)
 1000 return
      end
c****************************************************************
      function racah(j1,j2,j3,j4,j5,j6)
        use factorials
        implicit real*8(a-h,o-z),integer*4(i-n)
        integer dim
!        real*8,pointer :: fact(:)
!        dimension fact(0:51)
!      common /factrl/f(52),faclog(120)
!      equivalence (fact(0),f(1))
!!      fact(0)=f(1)

!!  AMORO (27/11/2003) 
!!
!!       fact(0:51)=f(1:52) !old statement
        
 !      write(*,*)'fact(0,1,2)=',fact(0:2)
!      stop
c      calculates a racah coefficient
      racah = 0.0
      z1 = delr(j1,j2,j5)
      if(z1.eq.0.0) go to 4
      z1 = delr(j3,j4,j5)*z1
      if(z1.eq.0.0) go to 4
      z2 = delr(j1,j3,j6)
      if(z2.eq.0.0) go to 4
      z2 = delr(j2,j4,j6)*z2
      if(z2.eq.0.0) go to 4
      z1 = sqrt(z1/z2)*z2
      jt1 = (j1+j2+j5)/2
      jt2 = (j3+j4+j5)/2
      jt3 = (j1+j3+j6)/2
      jt4 = (j2+j4+j6)/2
      jz1 = (j1+j2+j3+j4)/2
      jz2 = (j1+j4+j5+j6)/2
      jz3 = (j2+j3+j5+j6)/2
      numin = max0(jt1,jt2,jt3,jt4)
      numax = min0(jz1,jz2,jz3)
      if(numax.lt.numin) go to 4
      phase = phasef(numin+jz1)*z1
      do 3 nu=numin,numax
      jy1 = nu-jt1
      jy2 = nu-jt2
      jy3 = nu-jt3
      jy4 = jz1-nu
      jy5 = jz2-nu
      jy6 = jz3-nu
      if(numin.gt.50) go to 1
      fctor = fact(jy1)*fact(jy2)*fact(jy3)*yxfct(nu+1,nu-jt4)*
     >  fact(jy4)*fact(jy5)*fact(jy6)
      go to 2
    1 continue
      fctorl = flog(jy1+1)+flog(jy2+1)+flog(jy3+1)+flog(jy4+1)+
     >  flog(jy5+1)+flog(jy6+1)+flog(nu-jt4+1)-flog(nu+2)
      fctor=exp(fctorl)
    2 continue
      racah = racah+phase/fctor
      phase = -phase
    3 continue
c      debug unit (6),subchk
    4 return
      end


      function delr(j1,j2,j3)
        use factorials
      implicit real*8(a-h,o-z),integer*4(i-n)
!      dimension fact(0:51)
!      common /factrl/f(52),faclog(120)
!      equivalence (fact(0),f(1))
!      fact(0:51)=f(1:52)
      jz1 = (j1+j2-j3)/2
      if(jz1.lt.0) go to 3
      jz2 = (j1-j2+j3)/2
      if(jz2.lt.0) go to 3
      jz3 = (j2+j3-j1)/2
      if(jz3.lt.0) go to 3
      jz4 = (j1+j2+j3)/2+1
      if(jz3.lt.jz2) go to 2
      if(jz3.lt.jz1) go to 1
!      write(40,'(4i4,4g12.6)')jz1,jz2,jz3,jz4,
!     & fact(jz1),fact(jz2),fact(jz3),fact(jz4)
      delr = yxfct(jz4,jz3)*fact(jz1)*fact(jz2)
      return
    1 delr = yxfct(jz4,jz1)*fact(jz2)*fact(jz3)
      return
    2 if(jz2.lt.jz1) go to 1
      delr = yxfct(jz4,jz2)*fact(jz1)*fact(jz3)
      return
    3 delr = 0.0
      return
      end


c******************************************************************
       subroutine pl(x,l,pol,lmx)
c
c         computes legendre polynominals of order zero to order
c         (l) for a given argument (x= cos(theta))
c         also valid for x>1
c
      implicit double precision (a-h,o-z),integer*4(i-n)
      dimension pol(lmx)
      pol(1)=1.
      pol(2)=x
      if(l.le.1)goto20
      do 10 jj=2,l
      s=jj
   10 pol(jj+1)=((2.*s-1.)*x*pol(jj)-(s-1.)*pol(jj-1))/s
   20 return
      end

!****************************************************************
      function cint2db(xtab,ytab,fxytab,xbar,ybar,nnx,nny,
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
      complex*16 :: fxy(10,10),cint2db
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
      cint2db=fxy(nnn,nnn)
c      write(3,*) 'cint2d=',cint2db
c      print*,'Saliendo de cint2db'
      return
      end


c******************************************************
      subroutine newgd(lmu,kmax,xk,xkp,din,dout)
        use jat1
        use parms
        implicit none

        integer :: lmu,kmax,ll,i,j,n,norder
        real*8 :: xk,xkp
        complex*16 :: cint2db,dout(:)
        complex*16, dimension(:,:), intent(in) :: din
        norder = 5
!        write(98,*) 'Entering newgd din=',lmu,kmax,xk,xkp,din
        if (.not.allocated(ttl)) then
           write(*,*) 'ERROR: ttl not allocated. Aborting'
           stop
        endif
        do 50 ll=1,lmu
           do 100 i=1,kmax
              do 100 j=i,kmax
                 n = (i-1)*(2*kmax-i)/2+j
                 ttl(i,j) = -din(ll,n)/ttof
                 ttl(j,i) = ttl(i,j)
!                 write(4,*) i,j,ttl(i,j)
100        continue      
!        write(3,*) 'newgd: xx(1)=',xx(1),xk,xkp,kmax,ttl(1,1)
        dout(ll) = cint2db(xx,xx,ttl,xk,xkp,kmax,kmax,norder,kmax)

  50  continue
      return
      end
 
c******************************************************
      subroutine newgnd(lmu,kmax,xk,xkp,din,dout)
        use jat1
        use parms
c***  Change the sign of the off-diagonal term to be consistent
c     with the i(L-Lp) factor in Goddard subroutine, AMPCAL
c***
      implicit none
      integer :: ll,j,i,lmu,kmax,norder
      real*8 :: xkp,xk
      complex*16, dimension(:,:,:) :: din
      complex*16 ::  dout(:),cint2db
!     complex*16 ::  dout(lmu) !WRONG!!!!!     
      norder = 5
      do 50 ll=1,lmu
        do 100 i=1,kmax
        do 100 j=1,kmax
        ttl(i,j) = -din(ll,i,j)/ttof
 100    continue
      dout(ll) = cint2db(xx,xx,ttl,xk,xkp,kmax,kmax,norder,kmax)
c     dout(ll) = - dout(ll)
  50  continue
      return
      end
 
 
c **********************************************************************
      subroutine ampcal (theta,lmu,a,b,c,d,e,a0)
      use ampnew
      use parms
!      use constants
      implicit none
      integer :: lmu,kr,kl,km,k,ir
      real*8 :: theta, x,s,sq2,c1,c2,c3,pauli,a2l1,
     1        al,al1,al2,c5,c6,c7,c4,alm
      complex*16 :: a0(2),a(2),b(2),c(2),d(2),e(2),u,up
      real*8:: plm(3,3)


c **      kr=1(T=0), kr=2(T=1)
 
      x = cos (theta)
      s = sin (theta)
      plm(2,1) = x
      plm(3,1) = 1.5*x*x - 0.5
      plm(2,2) = -s
      plm(3,2) = -3.0*s*x
      plm(2,3) = 0.0
      plm(3,3) = 3.0*s*s
      sq2 = dsqrt(2.d0)
      c1 = dsqrt(3.d0)
      c2 = sq2*c1
      c3 = -c1/sq2
*----------------------------------------------------------------------
*     sum lowest l parts explicitly
*----------------------------------------------------------------------
      do 20 kr=1,2
      a0(kr) = d0(1)*pauli(1-kr)
     +        + 3.*d0(2)*plm(2,1)*pauli(kr-2)
c     write(22,*)'kr,a0=',kr,a0(kr)
      a(kr) = (dp(1) + dp(1) - sq2*dmm(1))*pauli(1-kr)
     + + (3.0*(dp(2) + dz(2)) - c2*dmm(2))*plm(2,1)*pauli(kr-2)
      b(kr) = (c3*dmm(2) - dp(2) + dm(2))*plm(2,2)*pauli(kr-2)
      c(kr) = 0.0
      d(kr) =(c3*dmm(2) + 1.5*(dp(2) - dz(2)))*plm(2,2)*pauli(kr-2)
      e(kr) = (sq2*dmm(1) + dp(1))*pauli(1-kr)
     + + (c2*dmm(2) + dp(2) + dp(2) + dm(2))*plm(2,1)*pauli(kr-2)
   20 continue


*----------------------------------------------------------------------
*     sum rest of l contributions
*----------------------------------------------------------------------
c **  L+S+T = even
 
      do 100 kl=3, lmu
      a2l1 = kl + kl - 1
      alm = kl - 2
      al = kl - 1
      al1 = kl
      al2 = kl + 1
      c5 = a2l1/(al*al1)
      c6 = al2/al1
      c7 = alm/al
      c1 = dsqrt (alm*al)
      c2 = dsqrt (al2*al1)
      c3 = dsqrt (c7)
      c4 = dsqrt (c6)
      do 50 kr = 1,2
      if (kr.eq.1)then
        if(mod(kl,2).eq.0)go to 50
      else if (kr.eq.2)then
        if(mod(kl,2).ne.0)go to 50
      endif
      u = c1*dpp(kl) + c2*dmm(kl)
      a0(kr) = d0(kl)*a2l1*plm(3,1) + a0(kr)
c     write(22,*)kl,kr,a0(kr)
      a(kr) = a(kr)+
     &       (al2*dp(kl) + a2l1*dz(kl) + alm*dm(kl)-u)*plm(3,1)
      e(kr) = (al1*dp(kl) + al*dm(kl) + u)*plm(3,1)
     + +e(kr)
      u = c3*dpp(kl) - c4*dmm(kl)
      b(kr) = (dm(kl) - dp(kl) + u)*plm(3,2) + b(kr)
      d(kr) = (c6*dp(kl) - c5*dz(kl) - c7*dm(kl)
     + + u)*plm(3,2) + d(kr)
      c(kr) = (dp(kl)/al1 - c5*dz(kl) + dm(kl)/al
     + - dpp(kl)/c1 - dmm(kl)/c2)*plm(3,3) + c(kr)
 50   continue
      do 70 km=1,3
      do 60 k=1,2
      plm(k,km) = plm(k+1,km)
   60 continue
      c1 = kl + km - 2
      c2 = kl - km + 1
      plm(3,km) = (a2l1*x*plm(2,km) - c1*plm(1,km))/c2
   70 continue
  100 continue


      c1 = 2.*pi
      c2 = sq2*c1
      do 300 ir=1,2
      a0(ir) = 2.*c1*a0(ir)
      a(ir) = c1*a(ir)
      b(ir) = c2*b(ir)
      c(ir) = c1*c(ir)
      d(ir) = c2*d(ir)
      e(ir) = 2.*c1*e(ir)



*----------------------------------------------------------------------
  300 continue
      u = a0(1)
      up = a0(2)
      a0(1)=up
      a0(2)=u
      return
      end
c ******************************************************
      function pauli(ix)
      implicit real*8(a-h,o-z)
      integer::ix
      if(ix.eq.0)then
      pauli = 1d0
      else
      pauli = 0d0
      endif
      return
      end
c
c **************************************************************
      subroutine badpt(eps,nqmax,kkmax,ir,ee)
        use j
        use rc
        use jat1
        implicit none
        integer np
        integer :: im,ir,nqmax,kkmax,jq,ik,jqn
        real*8 :: eps,xq,xkk
        complex*16 :: ee,wxxi
        real*8, allocatable :: exx(:)
        complex*16, allocatable :: eeaux(:) 
        dimension ee(2,kkmax,nqmax)
!       dimension eeaux(np),exx(np)

        allocate(eeaux(nqmax))
        allocate(exx(nqmax))
       

        do 100 jq = 1,nqmax
           do 100 ik = 1,kkmax
              if (abs(ctheta(ik,jq)).lt.eps)then
                 xq = xxq(jq)
                 xkk = xxk(ik)
                 im = 0
                 do jqn =1,nqmax
                    if (abs(ctheta(ik,jqn)).gt.eps)then
                       im = im+1
                       eeaux(im) = ee(ir,ik,jqn)
                       exx(im) = xxq(jqn)
                    end if
                 enddo
 !                write(*,*)'ee:nqmax,im,exx=',nqmax,im,exx(1:im)
!! AMORO (26/11/2003)
!!                 ee(ir,ik,jq) = wxxi(xq,exx,eeaux,nqmax,im)
                  ee(ir,ik,jq) = wxxi(xq,exx,eeaux,im,im)
              endif
100           continue
              return
           end 
c

c********************************************************
      subroutine onamp(mat,xkk,xq,kkmax,nqmax,aaon,ir)
        use j
        use ampaux
        use rc
        use jat1
        use toff
        implicit none
        integer :: kkmax,nqmax,ir,norder,jq,ik
        real*8 :: xkk,xq
        complex*16 :: mat(2,kkmax,nqmax),cint2db,aaon
        norder = 5
        do jq=1,nqmax
           do ik=1,kkmax
              ampl(ik,jq) = mat(ir,ik,jq)
           enddo
        enddo
        aaon=cint2db(xxk,xxq,ampl,xkk,xq,kkmax,nqmax,norder,kmax)
       return
      end


c*******************************************************
c      subroutine mesh(xkmax,xqmax,dk,dq,kkmax,nqmax,xxk,xxq,np)
       subroutine mesh(xkmax,xqmax,dk,dq,kkmax,nqmax)
        use block
        use j
        implicit none
c        implicit real*8(a-h,o-z)
        integer :: kkmax,nqmax,np,ik,jq
        real*8 :: xkmax,xqmax,dq,dk,xq,xkk
c
c          sets up evenly spaced k,q mesh
c

      kkmax=xkmax/dk+1.001
      nqmax=xqmax/dq+1.001

      write(99,2030) kkmax,nqmax
 2030 format('+ mesh: allocating',i3,' elements for xxk, xxq')
      if (.not.allocated(xxk)) then 
         allocate(xxk(kkmax))
         allocate(xxq(nqmax))
      endif

      if (.not.allocated(Bq)) then
         allocate (Bq(kkmax))
         allocate (Sq(nqmax))
      endif

      do 10 ik=1,kkmax
      xkk=(ik-1)*dk
      if(ik.eq.1)xkk=dk/100.
      xxk(ik)=xkk
      Bq(ik)=xxk(ik)
 10   continue
      do 20 jq=1,nqmax
      xq=(jq-1)*dq
      if(jq.eq.1)xq=dq/100.
      xxq(jq)=xq
      Sq(jq)=xxq(jq)
 20   continue
      return
 30   print 40
c 30    continue
 40   format(' too many points')
      stop
      end
c
c     library routines for LS programs


c*************************************************************
       subroutine potslj
c*************************************************************
      use toff
      use wk
      use parms

c     Paris-80 subroutine
c
c          generates paris potential vslj(k,k`) for given l, s, j.
c
      implicit real*8 (a-h,o-z)
      integer,save :: ny=12
      real*8:: v0(12,7,2),xmu(12,2),xmt(2),q(30),aux

!          data for paris potential
!      data ny /12/ !AMORO 10/4/04
      data (xmu(i,1),i=1,12)/
     $ 0.69953600d+00,0.16000000d+01,0.23000000d+01,0.30000000d+01,
     $ 0.37000000d+01,0.44000000d+01,0.51000000d+01,0.58000000d+01,
     $ 0.65000000d+01,0.82000000d+01,0.99000000d+01,0.11300000d+02/
      data (v0(i,1,1),i=1,12)/
     $ 0.32290874d+02,-.82465631d+02,0.12329384d+04,-.16859879d+05,
     $ 0.17292683d+06,-.76835277d+06,0.21890475d+07,-.38447287d+07,
     $ 0.27990559d+07,0.50251828d+06,-.26006124d+07,0.15528303d+07/
      data (v0(i,2,1),i=1,12)/
     $ -.10763625d+02,-.42973669d+02,-.71856844d+03,0.42469120d+04,
     $ -.34574024d+05,0.12671169d+06,-.27416841d+06,0.52960724d+06,
     $ -.36606713d+06,-.22303673d+06,0.40683833d+06,-.17678841d+06/
      data (v0(i,3,1),i=1,12)/
     $ -.85980096d-02,0.26814385d-01,-.13280693d+01,0.10324289d+02,
     $ -.11527067d+03,0.69456175d+03,-.23879335d+04,0.42388011d+04,
     $ -.24521604d+04,-.19512821d+04,0.41801160d+04,-.22508692d+04/
      data (v0(i,4,1),i=1,12)/
     $ 0.28660032d-02,-.81798046d-03,-.53314560d+00,0.83162030d+00,
     $ -.31192395d+02,0.30041384d+03,-.12415067d+04,0.24762241d+04,
     $ -.13043030d+04,-.21496577d+04,0.40996917d+04,-.22000183d+04/
      data (v0(i,5,1),i=1,12)/
     $ 0.00000000d+00,-.66176421d+02,0.28903688d+04,-.62592400d+05,
     $ 0.69146141d+06,-.40969146d+07,0.14032093d+08,-.26827468d+08,
     $ 0.23511442d+08,-.14688461d+08,0.12206548d+08,-.47565637d+07/
      data (v0(i,6,1),i=1,12)/
     $ -.10763625d+02,-.46818029d+00,0.60147739d+02,0.35256941d+03,
     $ 0.51432170d+03,0.11637302d+05,-.44595415d+05,0.69211738d+05,
     $ -.48127668d+05,0.70514008d+04,0.30188490d+05,-.29444569d+05/
      data (v0(i,7,1),i=1,12)/
     $ 0.00000000d+00,-.62851020d+00,-.76290197d+02,-.78827581d+03,
     $ -.64904798d+04,0.54734378d+04,-.32941912d+05,0.24949132d+06,
     $ -.16012956d+05,-.86736090d+06,0.55455935d+06,0.18614661d+06/
      data (xmu(i,2),i=1,12)/
     $ 0.68402600d+00,0.16000000d+01,0.23000000d+01,0.30000000d+01,
     $ 0.37000000d+01,0.44000000d+01,0.51000000d+01,0.58000000d+01,
     $ 0.65000000d+01,0.82000000d+01,0.99000000d+01,0.11300000d+02/
      data (v0(i,1,2),i=1,12)/
     $ -.10077427d+02,-.12049564d+03,-.21236460d+03,-.87174198d+04,
     $ 0.54383377d+05,-.21342147d+06,0.49458357d+06,-.66715334d+06,
     $ 0.52957598d+06,-.13703412d+06,-.34697194d+06,0.28509944d+06/
      data (v0(i,2,2),i=1,12)/
     $ 0.33591422d+01,-.86479568d+02,-.46593111d+03,0.18673085d+04,
     $ 0.38509213d+04,-.19674338d+05,0.12323140d+06,-.31449361d+06,
     $ 0.24242440d+06,0.16690404d+06,-.48534364d+06,0.27678611d+06/
      data (v0(i,3,2),i=1,12)/
     $ 0.26851393d-02,0.51092455d-01,-.84264258d+00,0.14736312d+02,
     $ -.14521993d+03,0.84158389d+03,-.27861170d+04,0.50564510d+04,
     $ -.33674205d+04,-.17845529d+04,0.53548266d+04,-.32465460d+04/
      data (v0(i,4,2),i=1,12)/
     $ -.89504644d-03,0.37488481d-01,-.89373089d+00,0.14123475d+02,
     $ -.14660152d+03,0.84191462d+03,-.28394273d+04,0.52653427d+04,
     $ -.35000430d+04,-.24879479d+04,0.73068121d+04,-.45574733d+04/
      data (v0(i,5,2),i=1,12)/
     $ 0.00000000d+00,-.42600359d+03,0.26279517d+05,-.57557033d+06,
     $ 0.60033934d+07,-.34519443d+08,0.11355459d+09,-.20729209d+09,
     $ 0.17131548d+09,-.86418222d+08,0.56271580d+08,-.18345800d+08/
      data (v0(i,6,2),i=1,12)/
     $ 0.33591422d+01,-.85945824d+00,-.10476340d+03,0.12629465d+04,
     $ -.18881061d+05,0.10613246d+06,-.33211910d+06,0.55585762d+06,
     $ -.34916664d+06,-.11945013d+06,0.32952970d+06,-.17079672d+06/
      data (v0(i,7,2),i=1,12)/
     $ 0.00000000d+00,-.52218640d+00,0.18644558d+03,-.37091115d+04,
     $ 0.55913117d+05,-.36998560d+06,0.14537543d+07,-.31352471d+07,
     $ 0.24339081d+07,0.14589455d+07,-.52668478d+07,0.34496510d+07/
!      data (xmt(k),k=1,2)/938.9055d0,5629.d0/
      data (xmt(k),k=1,2)/938.9055d0,938.2592d0/
      k1=kmax+1
      write(99,2020) 2*k1,2*k1
 2020 format('+Allocating memory for ',i3,'x',i3,' elements in  v')
      allocate(vl(2*k1+1,2*k1+1),stat=istat)
      if (istat>0) then
        write(99,*)'Memory allocationg for vl failed. Exiting'
        stop
      endif

c          singlet state lsj matrix elements
      s12=0.0d0
      so2=-dble(ll)*(dble(ll) + 1d0)
      xls=0.d0
      ll2=2*ll
      jj2=2*jj
      if(is.eq.0)go to 50
c         triplet state lsj matrix elements
      c1=-2.d0*dsqrt(30.d0)*(-1.d0)**jj
      x3=threj0(ll2,ll2,4)
      x6=racah(ll2,ll2,2,2,4,jj2)
      y6=racah(ll2,ll2,2,2,4,ll2)
      s12=c1*(ll2+1.d0)*x3*x6
      so2=ll*(ll+1.d0)*(1.d0/3.d0-(1.d0-2.d0*mod(ll+jj,2))*10.d0
     1 *(2.d0*ll+1.d0)*x6*y6)
      xls=jj*(jj+1.d0)/2.d0-ll*(ll+1.d0)/2.d0-1.d0
      if(ic.eq.1)go to 50
c          coupled state lsj matrix elements
c          l=l`=j-1
      s12p=0.0d0
      so2p=0.d0
      llp=ll-2
      llp2=2*llp
      xlsp=jj*(jj+1.d0)/2.d0-llp*(llp+1.d0)/2.d0-1.d0
      if(llp.eq.0)go to 30
      xp3=threj0(llp2,llp2,4)
      xp6=racah(llp2,llp2,2,2,4,jj2)
      yp6=racah(llp2,llp2,2,2,4,llp2)
      s12p=c1*(llp2+1.d0)*xp3*xp6
      so2p=llp*(llp+1.d0)*(1.d0/3.d0-(1.d0-2.d0*mod(llp+jj,2))*10.d0
     1 *(2.d0*llp+1.d0)*xp6*yp6)
c          l-prime=l-2
 30   xod3=threj0(llp2,ll2,4)
      xod6=racah(llp2,ll2,2,2,4,jj2)
      sod12=c1*dsqrt((ll2+1.d0)*(llp2+1.d0))*xod3*xod6

 50   continue
c          set up limits for l for ql loop
      nli=ll-1
      if(ic.eq.2)nli=llp-1
      if(nli.le.0)nli=1
      nlf=ll+3
c          begin k, k-prime loops
      do 10 i=1,k1
      xk=qq(i)
      do 10 j=i,k1
      xkp=qq(j)
      a=xk*xkp*2.d0
      b=xk*xk+xkp*xkp
      vv=0.d0
      vvp=0.d0
      vvod=0.d0
      vvodr=0.d0
c          begin loop for 12 paris potential terms
      do 20 n=1,ny
      xmm=xmu(n,it+1)
      xx=(b+xmm*xmm)/a
c          calculate required ql`s
      do 40 nl=nli,nlf
      l=nl-1
      call ql(xx,l,aux)
      q(nl)=aux
 40   continue
c          set up potential strengths
      vc0=v0(n,is+1,it+1)+v0(n,is+3,it+1)*b*ch*ch/xmt(it+1)
      vls0=v0(n,5,it+1)
      vt0=v0(n,6,it+1)
      vso20=v0(n,7,it+1)
c          add term to potential--uncoupled and larger l coupled
      vv=vv+(vc0-vls0*xls/3.)*all2(xk,xkp,xmm,ll,q)
     1     +vso20*so2*cll(xk,xkp,xmm,ll,q)
      if(is.eq.1)vv=vv+(vls0*xls/3.d0+vt0*s12)*bll(xk,xkp,xmm,ll,q)
      if(ic.eq.1)go to 20
c          lower l coupled
      vvp=vvp+(vc0-vls0*xlsp/3.)*all2(xk,xkp,xmm,llp,q)
     1     +vso20*so2p*cll(xk,xkp,xmm,llp,q)
     2     +(vls0*xlsp/3.d0+vt0*s12p)*bll(xk,xkp,xmm,llp,q)

c          off diagonal coupled
      vvod=vvod+vt0*sod12*bllp(xk,xkp,xmm,ll,q)
      vvodr=vvodr+vt0*sod12*bllp(xkp,xk,xmm,ll,q)
 20   continue
c          end of loop for paris terms
c          fill potential array for this k, k`
 80   vl(i,j)=vv/ch
      vl(j,i)=vl(i,j)
!      write(97,*) xk,' ',xkp,' (i,j)',i,j,' vl(i,j)=',vl(i,j)
      if(ic.eq.1)go to 10
      vl(i+k1,j+k1)=vvp/ch
      vl(j+k1,i+k1)=vvp/ch
      vl(i,j+k1)=vvodr/ch
      vl(j+k1,i)=vvodr/ch
      vl(j,i+k1)=vvod/ch
      vl(i+k1,j)=vvod/ch
 10   continue
c          end of k, k` loops
      write(99,*)'Leaving potlsj...'
      return
      end



!!$
!!$
!!$
!!$
!!$
!!$
!!$c*************************************************************
!!$       subroutine potslj
!!$c*************************************************************
!!$      use toff
!!$      use wk
!!$      use parms
!!$
!!$c     Paris-80 subroutine
!!$c
!!$c          generates paris potential vslj(k,k`) for given l, s, j.
!!$c
!!$      implicit real*8 (a-h,o-z)
!!$      real*8 two,four
!!$
!!$      real*8:: v0(12,7,2),xmu(12,2),xmt(2),q(30),aux
!!$
!!$!          data for paris potential
!!$!!      integer,save :: ny=12
!!$      data ny /12/ !AMORO 10/4/04
!!$      data (xmu(i,1),i=1,12)/
!!$     $ 0.69953600d+00,0.16000000d+01,0.23000000d+01,0.30000000d+01,
!!$     $ 0.37000000d+01,0.44000000d+01,0.51000000d+01,0.58000000d+01,
!!$     $ 0.65000000d+01,0.82000000d+01,0.99000000d+01,0.11300000d+02/
!!$      data (v0(i,1,1),i=1,12)/
!!$     $ 0.32290874d+02,-.82465631d+02,0.12329384d+04,-.16859879d+05,
!!$     $ 0.17292683d+06,-.76835277d+06,0.21890475d+07,-.38447287d+07,
!!$     $ 0.27990559d+07,0.50251828d+06,-.26006124d+07,0.15528303d+07/
!!$      data (v0(i,2,1),i=1,12)/
!!$     $ -.10763625d+02,-.42973669d+02,-.71856844d+03,0.42469120d+04,
!!$     $ -.34574024d+05,0.12671169d+06,-.27416841d+06,0.52960724d+06,
!!$     $ -.36606713d+06,-.22303673d+06,0.40683833d+06,-.17678841d+06/
!!$      data (v0(i,3,1),i=1,12)/
!!$     $ -.85980096d-02,0.26814385d-01,-.13280693d+01,0.10324289d+02,
!!$     $ -.11527067d+03,0.69456175d+03,-.23879335d+04,0.42388011d+04,
!!$     $ -.24521604d+04,-.19512821d+04,0.41801160d+04,-.22508692d+04/
!!$      data (v0(i,4,1),i=1,12)/
!!$     $ 0.28660032d-02,-.81798046d-03,-.53314560d+00,0.83162030d+00,
!!$     $ -.31192395d+02,0.30041384d+03,-.12415067d+04,0.24762241d+04,
!!$     $ -.13043030d+04,-.21496577d+04,0.40996917d+04,-.22000183d+04/
!!$      data (v0(i,5,1),i=1,12)/
!!$     $ 0.00000000d+00,-.66176421d+02,0.28903688d+04,-.62592400d+05,
!!$     $ 0.69146141d+06,-.40969146d+07,0.14032093d+08,-.26827468d+08,
!!$     $ 0.23511442d+08,-.14688461d+08,0.12206548d+08,-.47565637d+07/
!!$      data (v0(i,6,1),i=1,12)/
!!$     $ -.10763625d+02,-.46818029d+00,0.60147739d+02,0.35256941d+03,
!!$     $ 0.51432170d+03,0.11637302d+05,-.44595415d+05,0.69211738d+05,
!!$     $ -.48127668d+05,0.70514008d+04,0.30188490d+05,-.29444569d+05/
!!$      data (v0(i,7,1),i=1,12)/
!!$     $ 0.00000000d+00,-.62851020d+00,-.76290197d+02,-.78827581d+03,
!!$     $ -.64904798d+04,0.54734378d+04,-.32941912d+05,0.24949132d+06,
!!$     $ -.16012956d+05,-.86736090d+06,0.55455935d+06,0.18614661d+06/
!!$      data (xmu(i,2),i=1,12)/
!!$     $ 0.68402600d+00,0.16000000d+01,0.23000000d+01,0.30000000d+01,
!!$     $ 0.37000000d+01,0.44000000d+01,0.51000000d+01,0.58000000d+01,
!!$     $ 0.65000000d+01,0.82000000d+01,0.99000000d+01,0.11300000d+02/
!!$      data (v0(i,1,2),i=1,12)/
!!$     $ -.10077427d+02,-.12049564d+03,-.21236460d+03,-.87174198d+04,
!!$     $ 0.54383377d+05,-.21342147d+06,0.49458357d+06,-.66715334d+06,
!!$     $ 0.52957598d+06,-.13703412d+06,-.34697194d+06,0.28509944d+06/
!!$      data (v0(i,2,2),i=1,12)/
!!$     $ 0.33591422d+01,-.86479568d+02,-.46593111d+03,0.18673085d+04,
!!$     $ 0.38509213d+04,-.19674338d+05,0.12323140d+06,-.31449361d+06,
!!$     $ 0.24242440d+06,0.16690404d+06,-.48534364d+06,0.27678611d+06/
!!$      data (v0(i,3,2),i=1,12)/
!!$     $ 0.26851393d-02,0.51092455d-01,-.84264258d+00,0.14736312d+02,
!!$     $ -.14521993d+03,0.84158389d+03,-.27861170d+04,0.50564510d+04,
!!$     $ -.33674205d+04,-.17845529d+04,0.53548266d+04,-.32465460d+04/
!!$      data (v0(i,4,2),i=1,12)/
!!$     $ -.89504644d-03,0.37488481d-01,-.89373089d+00,0.14123475d+02,
!!$     $ -.14660152d+03,0.84191462d+03,-.28394273d+04,0.52653427d+04,
!!$     $ -.35000430d+04,-.24879479d+04,0.73068121d+04,-.45574733d+04/
!!$      data (v0(i,5,2),i=1,12)/
!!$     $ 0.00000000d+00,-.42600359d+03,0.26279517d+05,-.57557033d+06,
!!$     $ 0.60033934d+07,-.34519443d+08,0.11355459d+09,-.20729209d+09,
!!$     $ 0.17131548d+09,-.86418222d+08,0.56271580d+08,-.18345800d+08/
!!$      data (v0(i,6,2),i=1,12)/
!!$     $ 0.33591422d+01,-.85945824d+00,-.10476340d+03,0.12629465d+04,
!!$     $ -.18881061d+05,0.10613246d+06,-.33211910d+06,0.55585762d+06,
!!$     $ -.34916664d+06,-.11945013d+06,0.32952970d+06,-.17079672d+06/
!!$      data (v0(i,7,2),i=1,12)/
!!$     $ 0.00000000d+00,-.52218640d+00,0.18644558d+03,-.37091115d+04,
!!$     $ 0.55913117d+05,-.36998560d+06,0.14537543d+07,-.31352471d+07,
!!$     $ 0.24339081d+07,0.14589455d+07,-.52668478d+07,0.34496510d+07/
!!$!      data (xmt(k),k=1,2)/938.9055d0,5629.d0/
!!$      data (xmt(k),k=1,2)/938.9055d0,938.2592d0/
!!$      k1=kmax+1
!!$      two=2d0
!!$      four=4d0
!!$      write(99,2020) 2*k1,2*k1
!!$ 2020 format('+Allocating memory for ',i3,'x',i3,' elements in  v')
!!$      allocate(vl(2*k1+1,2*k1+1),stat=istat)
!!$      if (istat>0) then
!!$        write(99,*)'Memory allocationg for vl failed. Exiting'
!!$        stop
!!$      endif
!!$      
!!$      write(97,*)is,ll,jj
!!$     
!!$       
!!$
!!$c          singlet state lsj matrix elements
!!$      s12=0.0d0
!!$      so2=-dble(ll)*(dble(ll) + 1d0)
!!$      xls=0.d0
!!$      ll2=2*ll
!!$      jj2=2*jj
!!$      if(is.eq.0)go to 50
!!$c         triplet state lsj matrix elements
!!$      c1=-2.d0*dsqrt(30.d0)*(-1.d0)**jj
!!$      x3=threj0(ll2,ll2,four)
!!$      
!!$      x6=racah(ll2,ll2,2,2,4,jj2)
!!$      y6=racah(ll2,ll2,2,2,4,ll2)
!!$
!!$!      x6=racah(ll2,ll2,two,two,four,jj2)
!!$!      y6=racah(ll2,ll2,two,two,four,ll2)
!!$
!!$      x6b= rac(dble(ll2),dble(ll2),two,two,four,dble(jj2))
!!$      y6b= rac(dble(ll2),dble(ll2),two,two,four,dble(ll2))
!!$
!!$!      if (is.eq.1) 
!!$      write(98,'(3i4,2x,6g12.5)') is,ll2,jj2,x6,y6,x6b,y6b,x3
!!$      s12=c1*(ll2+1.d0)*x3*x6
!!$      so2=ll*(ll+1.d0)*(1.d0/3.d0-(1.d0-2.d0*mod(ll+jj,2))*10.d0
!!$     1 *(2.d0*ll+1.d0)*x6*y6)
!!$      xls=jj*(jj+1.d0)/2.d0-ll*(ll+1.d0)/2.d0-1.d0
!!$      if(ic.eq.1)go to 50
!!$c          coupled state lsj matrix elements
!!$c          l=l`=j-1
!!$      s12p=0.0d0
!!$      so2p=0.d0
!!$      llp=ll-2
!!$      llp2=2*llp
!!$      xlsp=jj*(jj+1.d0)/2.d0-llp*(llp+1.d0)/2.d0-1.d0
!!$      if(llp.eq.0)go to 30
!!$      xp3=threj0(llp2,llp2,4)
!!$
!!$      xp6=racah(llp2,llp2,2,2,4,jj2)
!!$      yp6=racah(llp2,llp2,2,2,4,llp2)
!!$
!!$!      xp6b=rac(dble(llp2),dble(llp2),two,two,four,dble(jj2))
!!$!      yp6b=rac(dble(llp2),dble(llp2),two,two,four,dble(llp2))
!!$
!!$!     if (is.eq.1) write(98,'(50g12.5)') x6,y6,x6b,y6b
!!$
!!$      s12p=c1*(llp2+1.d0)*xp3*xp6
!!$      so2p=llp*(llp+1.d0)*(1.d0/3.d0-(1.d0-2.d0*mod(llp+jj,2))*10.d0
!!$     1 *(2.d0*llp+1.d0)*xp6*yp6)
!!$c          l-prime=l-2
!!$ 30   xod3=threj0(llp2,ll2,4)
!!$      xod6=racah(llp2,ll2,2,2,4,jj2)
!!$!       xod6=rac(dble(llp2),dble(ll2),two,two,four,dble(jj2))
!!$      sod12=c1*dsqrt((ll2+1.d0)*(llp2+1.d0))*xod3*xod6
!!$
!!$ 50   continue
!!$c          set up limits for l for ql loop
!!$      nli=ll-1
!!$      if(ic.eq.2)nli=llp-1
!!$      if(nli.le.0)nli=1
!!$      nlf=ll+3
!!$c          begin k, k-prime loops
!!$      do 10 i=1,k1
!!$      xk=qq(i)
!!$      do 10 j=i,k1
!!$      xkp=qq(j)
!!$      a=xk*xkp*2.d0
!!$      b=xk*xk+xkp*xkp
!!$      vv=0.d0
!!$      vvp=0.d0
!!$      vvod=0.d0
!!$      vvodr=0.d0
!!$c          begin loop for 12 paris potential terms
!!$      do 20 n=1,ny
!!$      xmm=xmu(n,it+1)
!!$      xx=(b+xmm*xmm)/a
!!$c          calculate required ql`s
!!$      do 40 nl=nli,nlf
!!$      l=nl-1
!!$      call ql(xx,l,aux)
!!$!      write(98,*) 'l,ql=',l,aux
!!$      q(nl)=aux
!!$ 40   continue
!!$c          set up potential strengths
!!$      vc0=v0(n,is+1,it+1)+v0(n,is+3,it+1)*b*ch*ch/xmt(it+1)
!!$      vls0=v0(n,5,it+1)
!!$      vt0=v0(n,6,it+1)
!!$      vso20=v0(n,7,it+1)
!!$c          add term to potential--uncoupled and larger l coupled
!!$      vv=vv+(vc0-vls0*xls/3.)*all2(xk,xkp,xmm,ll,q)
!!$     1     +vso20*so2*cll(xk,xkp,xmm,ll,q)
!!$      if(is.eq.1)vv=vv+(vls0*xls/3.d0+vt0*s12)*bll(xk,xkp,xmm,ll,q)
!!$      if(ic.eq.1)go to 20
!!$c          lower l coupled
!!$      vvp=vvp+(vc0-vls0*xlsp/3.)*all2(xk,xkp,xmm,llp,q)
!!$     1     +vso20*so2p*cll(xk,xkp,xmm,llp,q)
!!$     2     +(vls0*xlsp/3.d0+vt0*s12p)*bll(xk,xkp,xmm,llp,q)
!!$
!!$c          off diagonal coupled
!!$      vvod=vvod+vt0*sod12*bllp(xk,xkp,xmm,ll,q)
!!$      vvodr=vvodr+vt0*sod12*bllp(xkp,xk,xmm,ll,q)
!!$ 20   continue
!!$c          end of loop for paris terms
!!$c          fill potential array for this k, k`
!!$ 80   vl(i,j)=vv/ch
!!$      vl(j,i)=vl(i,j)
!!$      write(97,150) i,j,xk,xkp,vl(i,j)
!!$!xk,' ',xkp,' (i,j)',i,j,' vl(i,j)=',vl(i,j)
!!$150   format("(i,j)=",2i4," xk,xkp=",f8.4,2x,f8.4,3x,"vl(i,j)=",f10.6)
!!$      if(ic.eq.1)go to 10
!!$      vl(i+k1,j+k1)=vvp/ch
!!$      vl(j+k1,i+k1)=vvp/ch
!!$      vl(i,j+k1)=vvodr/ch
!!$      vl(j+k1,i)=vvodr/ch
!!$      vl(j,i+k1)=vvod/ch
!!$      vl(i+k1,j)=vvod/ch
!!$ 10   continue
!!$c          end of k, k` loops
!!$      write(99,*)'Leaving potlsj...'
!!$      return
!!$      end


c**************************************************
c   Not used (replaced by all2): "all" is protected
****************************************************
      function all(x,xp,xm,l,q)
      implicit double precision (a-h,o-z)
      dimension q(30)
      all=q(l+1)/2.d0/x/xp/xm
      return
      end

c******************************************
      function all2(x,xp,xm,l,q)
      implicit real*8 (a-h,o-z)
      integer::l
      real*8::x,xp,xm,q(30)
!      real*8:: q(*) ! AMORO 9/2/04
!      dimension q(30)
      all2=q(l+1)/2.d0/x/xp/xm
      return
      end   
c*******************************************
      function bll(x,xp,xm,l,q)
      implicit real*8 (a-h,o-z)
      real*8:: q(30) ! AMORO 9/2/04
!      dimension q(30)
      bll=0.0d0
      if(l.eq.0)return
      bll=q(l+1)/2.d0/x/xp/xm+3.d0*(q(l)-q(l+2))/2.d0/(2.d0*l+1.d0)
     1 /xm/xm/xm
      return
      end
c*******************************************
      function bllp(x,xp,xm,l,q)
      implicit real*8 (a-h,o-z)
      real*8:: q(30) ! AMORO 9/2/04
!      dimension q(30)
      bllp=(q(l+1)*x/2.d0/xp+q(l-1)*xp/2.d0/x-q(l))/xm/xm/xm
      return
      end
c*******************************************
      function cll(x,xp,xm,l,q)
      implicit real*8 (a-h,o-z)
      integer:: l
      real*8::x,xp,xm,q(30)
!      dimension q(30)
      cll=0.0
      if(l.eq.0)return
      c0=x*x*xp*xp/xm**4/(2.0*l+1.0)
      if(l.gt.1)go to 10
      xl=dlog((xm**2+(x+xp)**2)*(xm**2+(x-xp)**2)/xm**4)
      cll=c0*(-xl/4.0/x/xp/xm+all2(x,xp,xm,3,q)/5.0
     1 -6.0*all2(x,xp,xm,1,q)/5.0)
      return
 10   cll=c0*(all2(x,xp,xm,l-2,q)/(2.0*l-1.0)
     1       +all2(x,xp,xm,l+2,q)/(2.0*l+3.0)
     2       -(4.0*l+2.0)*all2(x,xp,xm,l,q)/(2.0*l-1.0)/(2.0*l+3.0))
      return
      end


c     ****************************************************************
      real*8 function stat(x)
!        use factorials
!        use dmats
      implicit real*8(a-h,o-z)
      stat=sqrt(2.d0*x+1.d0)
      return
      end
c     ******************************************************************
      real*8 function rmat(j,m,n)
        use factorials
        use dmats
      implicit real*8(a-h,o-z)
c     ------------------------------------------------------------------
c     small d rotation matrices for integer j. brink and satchler
c     formula. rotation of beta about y-axis.
c     sbeta2=sin(beta/2), cbeta=cos(beta/2).
c     ------------------------------------------------------------------
!      common/clebma/faclog(500)
!     common/dmats/sc,dc
      integer t
      t=0
      rmat=1.d0
      if(iabs(m).gt.j.or.iabs(n).gt.j) go to 888
      if(j.eq.0) return
      ind1=j+m+1
      ind2=j-n+1
      ind3=n-m+1
      rnum=(flog(ind1)+flog(j-m+1)+flog(ind2)+flog(j+n+1))
      rnum=rnum/2.d0
      rmat=0.d0
      go to 997
   88 den=flog(ind1)+flog(ind2)+flog(ind3)+flog(t+1)
      rmat=rmat+(-1.d0)**t*dexp(rnum-den)*dc**(2*(j-t)+m-n)
     1     *sc**(2*t+n-m)
   99 t=t+1
      ind1=ind1-1
      ind2=ind2-1
      ind3=ind3+1
  997 if((ind1.le.0).or.(ind2.le.0)) return
      if(ind3.le.0) go to 99
      go to 88
  888 rmat=0.d0
      return
      end
c     ****************************************************************
      real*8 function cleb(ria,rid,rib,rie,ric,rif)
        use factorials
        
      implicit real*8(a-h,o-z)
!      common/clebma/faclog(500)
      ia=2.d0*(ria+.0001d0)
      ib=2.d0*(rib+.0001d0)
      ic=2.d0*(ric+.0001d0)
      id=int(sign(1.d0,rid)*2.d0*(abs(rid)+.0001d0))
      ie=int(sign(1.d0,rie)*2.d0*(abs(rie)+.0001d0))
      if=int(sign(1.d0,rif)*2.d0*(abs(rif)+.0001d0))
      wwww=-1.0d0
      cleb=0.0d0
      if(id+ie-if) 7000,105,7000
  105 k1=ia+ib+ic
      if((-1)**k1) 7000,107,107
  107 if(.not.((id.eq.0).and.(ie.eq.0))) go to 110
      k1=k1/2
      if((-1)**k1) 7000,110,110
  110 k1=ia+ib-ic
      k2=ic-iabs(ia-ib)
      k3=min0(k1,k2)
      if(k3) 7000,130,130
  130 if((-1)**(ib+ie)) 7000,7000,140
  140 if((-1)**(ic+if)) 7000,7000,150
  150 if(ia-iabs (id)) 7000,152,152
  152 if(ib-iabs (ie)) 7000,154,154
  154 if(ic-iabs (if)) 7000,160,160
  160 if(ia) 7000,175,165
  165 if(ib) 7000,175,170
  170 if(ic) 7000,180,250
  175 cleb=1.0d0
      go to 7000
  180 fb=float(ib+1)
      cleb=((wwww)**((ia-id)/2))/sqrt(fb)
      go to 7000
  250 fc2=ic+1
      iabcp=(ia+ib+ic)/2+1
      iabc=iabcp-ic
      icab=iabcp-ib
      ibca=iabcp-ia
      iapd=(ia+id)/2+1
      iamd=iapd-id
      ibpe=(ib+ie)/2+1
      ibme=ibpe-ie
      icpf=(ic+if)/2+1
      icmf=icpf-if
      vvv=0.5d0
      sqfclg=vvv*(log(fc2)-flog(iabcp+1)
     1      +flog(iabc)+flog(icab)+flog(ibca)
     2      +flog(iapd)+flog(iamd)+flog(ibpe)
     3      +flog(ibme)+flog(icpf)+flog(icmf))
      nzmic2=(ib-ic-id)/2
      nzmic3=(ia-ic+ie)/2
      nzmi= max0(0,nzmic2,nzmic3)+1
      nzmx= min0(iabc,iamd,ibpe)
      if(nzmx.lt.nzmi) go to 7000
      s1=(wwww)**(nzmi-1)
      do 400 nz=nzmi,nzmx
      nzm1=nz-1
      nzt1=iabc-nzm1
      nzt2=iamd-nzm1
      nzt3=ibpe-nzm1
      nzt4=nz-nzmic2
      nzt5=nz-nzmic3
      termlg=sqfclg-flog(nz)-flog(nzt1)-flog(nzt2)
     1           -flog(nzt3)-flog(nzt4)-flog(nzt5)
      ssterm=s1*exp (termlg)
      cleb=cleb+ssterm
  400 s1=-s1
 7000 return
      end

      real*8 function rac(ria,rib,ric,rid,rie,rif)
        use factorials
c     -------------------------------------------------------------------
c     subroutine calculates the racah coefficient w(abcd;ef) defined
c     according to the convention of brink and satchler.
c     the arguments are real and are the actual values of the angular
c     momenta, ( i.e. they can take half integer values )
c     -------------------------------------------------------------------
      implicit real*8(a-h,o-z)
!      common/clebma/faclog(500)
      dimension lt(6)
      rac=0.0d0
      ia=2.d0*(ria+.0001d0)
      ib=2.d0*(rib+.0001d0)
      ic=2.d0*(ric+.0001d0)
      id=2.d0*(rid+.0001d0)
      ie=2.d0*(rie+.0001d0)
      if=2.d0*(rif+.0001d0)
      k1=ia+ib-ie
      k2=ie-iabs (ia-ib)
      k3=ic+id-ie
      k4=ie-iabs (ic-id)
      k5=ia+ic-if
      k6=if-iabs (ia-ic)
      k7=ib+id-if
      k8=if-iabs(ib-id)
      k9= min0 (k1,k2,k3,k4,k5,k6,k7,k8)
      if(k9) 7000,20,20
   20 k2=k1-2*(k1/2)
      k4=k3-2*(k3/2)
      k6=k5-2*(k5/2)
      k8=k7-2*(k7/2)
      if(max0(k2,k4,k6,k8)) 7000,25,7000
   25 ltmin=min0(ia,ib,ic,id,ie,if)
      if(ltmin) 7000,30,150
   30 lt(1)=ia
      lt(2)=ib
      lt(3)=ic
      lt(4)=id
      lt(5)=ie
      lt(6)=if
      ltmin=lt(1)
      kmin=1
      do 40 n=2,6
      if(lt(n)-ltmin) 35,40,40
   35 ltmin=lt(n)
      kmin=n
   40 continue
      s1=1.0d0
      f1=ie
      f2=if
      go to (55,55,55,55,45,50),kmin
   45 f1=ia
      f2=ic
      s1=(-1.d0)**(k5/2)
      go to 55
   50 f1=ia
      f2=ib
      s1=(-1.d0)**(k1/2)
   55 rac=s1/dsqrt((f1+1.d0)*(f2+1.d0))
      go to 7000
  150 iabep=(ia+ib+ie)/2+1
      icdep=(ic+id+ie)/2+1
      iacfp=(ia+ic+if)/2+1
      ibdfp=(ib+id+if)/2+1
      iabe=iabep-ie
      ieab=iabep-ib
      ibea=iabep-ia
      icde=icdep-ie
      iecd=icdep-id
      idec=icdep-ic
      iacf=iacfp-if
      ifac=iacfp-ic
      icfa=iacfp-ia
      ibdf=ibdfp-if
      ifbd=ibdfp-id
      idfb=ibdfp-ib
      iabcd1=(ia+ib+ic+id+4)/2
      iefmad=(ie+if-ia-id)/2
      iefmbc=(ie+if-ib-ic)/2
      nzmax=min0(iabe,icde,iacf,ibdf)
      nzmi1=-iefmad
      nzmi2=-iefmbc
      nzmin=max0(0,nzmi1,nzmi2)+1
      if(nzmax.lt.nzmin) go to 7000
      sqlog=flog(iabe)+flog(ieab)+flog(ibea)+flog(icde)+flog(i
     1ecd)+flog(idec)+flog(iacf)+flog(ifac)+flog(icfa)+flog(ib
     2df)+flog(ifbd)+flog(idfb)-flog(iabep+1)-flog(icdep+1)-
     & flog(iacfp+1)-flog(ibdfp+1)
      sqlog=0.5d0*sqlog
      do 200 nz=nzmin,nzmax
      nzm1=nz-1
      k1=iabcd1-nzm1
      k2=iabe-nzm1
      k3=icde-nzm1
      k4=iacf-nzm1
      k5=ibdf-nzm1
      k6=nz
      k7=iefmad+nz
      k8=iefmbc+nz
      sslog=sqlog+flog(k1)-flog(k2)-flog(k3)-flog(k4)
     1           -flog(k5)-flog(k6)-flog(k7)-flog(k8)
      ssterm=((-1.d0)**nzm1)*dexp(sslog)
      rac=rac+ssterm
!      write(96,'(20g12.6)') flog(k1),flog(k2),flog(k3),flog(k4),
!     1                      flog(k5),flog(k6),flog(k7),flog(k8) 
  200 continue
 7000 return
      end
c     ********************************************************************
      real*8 function u9(ra,rb,rc,rd,re,rf,rg,rh,ri)
c     -------------------------------------------------------------------
c     nine-j symbol. definition as in brink and satchler.
c     -------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      u9=0.0d0
      k1=idint(2.d0*dabs(ra-ri)+0.01d0)
      k2=idint(2.d0*dabs(rb-rf)+0.01d0)
      k3=idint(2.d0*dabs(rd-rh)+0.01d0)
      minrda=max0(k1,k2,k3)
      k1=idint(2.d0*(ra+ri)+0.01d0)
      k2=idint(2.d0*(rb+rf)+0.01d0)
      k3=idint(2.d0*(rd+rh)+0.01d0)
      maxrda=min0(k1,k2,k3)
      if(minrda-maxrda) 30,30,20
   30 do 50 n1=minrda,maxrda,2
      r1=float(n1)/2.d0
      ramda2=n1
      y9=(ramda2+1.d0)*rac(ra,ri,rd,rh,r1,rg)*rac(rb,rf,rh,rd,r1,re)
      u9=u9+y9*rac(ra,ri,rb,rf,r1,rc)
   50 continue
   20 return
      end
c     ****************************************************************
      real*8 function threej(ria,rid,rib,rie,ric,rif)
c     -------------------------------------------------------------------
c     three-j symbol. definition as in brink and satchler.
c     -------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      threej=cleb(ria,rid,rib,rie,ric,-rif)
      ipas=idint(dabs(ria-rib-rif)+0.01d0)
      threej=dble(float((-1)**ipas))*threej/stat(ric)
      return
      end







      function wxxi(r,rp,wfn,ndim,npts)
c****************************************************************************
c     calculates the interpolated function at point r, given the
c     function (in array wfn) at a set of  npts points (in array rp)
c------------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      integer :: ndim,npts
      complex*16 ::wfn(ndim),y1,y2,y3,y4,y5,y6
      real*8::rp(ndim),r,eps
      complex*16::wxxi
      
      eps=1e-8

!       write(*,*)'In wxxi ndim,npts=',ndim,npts
!       stop
!       write(*,'(7f12.6,1i3)')r,rp(1),rp(ndim),wfn(1),wfn(ndim),ndim
!       write(*,*)r,rp(1),rp(ndim),wfn(1),wfn(ndim),ndim
!       call flush(99)

      nst=0
      do 30 k=1,npts
      if(rp(k)>r) goto 33
!         write(99,*)'r,rp(k-1),rp(k+1)',r,rp(k-1),rp(k+1);goto 33
!      endif
   30 continue
      
33    nst=max0(k-3,1) 
      if ((k<ndim).and.(rp(k+1)<rp(k-1))) then
         write(*,*)'wxxi: x values are not in increasing order!:'
         write(*,*)'k=',k,'rp(k-1)=',rp(k-1),'rp(k+1)=',rp(k+1)
         stop
      endif
!      write(99,*)'nst,r,rp(k-1),rp(k+1)',nst,r,rp(k-1),rp(k+1)
      if((nst.gt.npts-5).or.(nst.eq.0)) nst=npts-5
      x1=rp(nst+0)
      x2=rp(nst+1)
      x3=rp(nst+2)
      x4=rp(nst+3)
      x5=rp(nst+4)
      x6=rp(nst+5)

      y1=wfn(nst+0)
      y2=wfn(nst+1)
      y3=wfn(nst+2)
      y4=wfn(nst+3)
      y5=wfn(nst+4)
      y6=wfn(nst+5)

      pii1=(x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)*(x1-x6)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)*(x2-x5)*(x2-x6)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)*(x3-x5)*(x3-x6)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)*(x4-x5)*(x4-x6)
      pii5=(x5-x1)*(x5-x2)*(x5-x3)*(x5-x4)*(x5-x6)
      pii6=(x6-x1)*(x6-x2)*(x6-x3)*(x6-x4)*(x6-x5)
  777 xd1=r-x1
      xd2=r-x2
      xd3=r-x3
      xd4=r-x4
      xd5=r-x5
      xd6=r-x6

!! AMORO addition (26/11/2003)
!!$      if (abs(xd1)<eps) wxxi=y1;return
!!$      if (abs(xd2)<eps) wxxi=y2;return
!!$      if (abs(xd3)<eps) wxxi=y3;return
!!$      if (abs(xd4)<eps) wxxi=y4;return
!!$      if (abs(xd5)<eps) wxxi=y5;return
!!$      if (abs(xd6)<eps) wxxi=y6;return


      pi1= xd2*xd3*xd4*xd5*xd6
      pi2= xd1*xd3*xd4*xd5*xd6
      pi3= xd1*xd2*xd4*xd5*xd6
      pi4= xd1*xd2*xd3*xd5*xd6
      pi5= xd1*xd2*xd3*xd4*xd6
      pi6= xd1*xd2*xd3*xd4*xd5
!!$      if ((abs(pii1)<1e-15).or.(abs(pii2)<1e-15).or.(abs(pii3)<1e-15)
!!$     & .or.(abs(pii4)<1e-15).or.(abs(pii5)<1e-15).or.(abs(pii6)<1e-15)) 
!!$     & then
!!$         write(99,*)'wxxi: pii1=',pii1
!!$         write(99,*)'wxxi: pii2=',pii2
!!$         write(99,*)'wxxi: pii3=',pii3
!!$         write(99,*)'wxxi: pii4=',pii4
!!$         write(99,*)'wxxi: pii5=',pii5
!!$         write(99,*)'wxxi: pii6=',pii6
!!$         write(99,fmt='"r,x1,x2,x3,x4,x5,x6=",7f10.5') 
!!$     & r,x1,x2,x3,x4,x5,x6 
!!$         wxxi=0.d0
!!$!         write(*,*)'Aborting';stop
!!$      else
         wxxi=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4+
     + y5*pi5/pii5+y6*pi6/pii6
!!$      endif
      return
      end

      subroutine sbesjh(x,lmax,xj,xjp,xh1,xh1p,ifail)
c ***                                                       i.j.thompson
c ***                                                       31 may 1985.
c ***  complex spherical bessel functions from l=0 to l=lmax
c ***    for x in the upper half plane ( im(x) > -3)
c ***
c ***    xj(l)   = j/l(x)          regular solution: xj(0)=sin(x)/x
c ***    xjp(l)  = d/dx j/l(x)
c ***    xh1(l)  = h(1)/l(x)       irregular hankel function:
c ***    xh1p(l) = d/dx h(1)/l(x)            xh1(0) = j0(x) + i. y0(x)
c ***                                               =(sin(x)-i.cos(x))/x
c ***                                               = -i.exp(i.x)/x
c ***  using complex cf1, and trigonometric forms for l=0 solutions.
c ***
      implicit complex*16 (a-h,o-z)
      parameter (limit=20000)
      dimension xj(0:lmax),xjp(0:lmax),xh1(0:lmax),xh1p(0:lmax)
      real*8 zero,one,accur,tm30,absc
      data zero,one/ 0.0d0,1.0d0 /, accur /1.0d-12/, tm30 / 1d-30 /,
     #     ci / (0d0,1d0) /
      absc(w) = abs(real(w)) + abs(dimag(w))
      ifail= -1
      if(absc(x).lt.accur .or. dimag(x).lt.-3.0) go to 5
      xi = one/x
      w  = xi + xi
      pl = lmax*xi
      f = pl + xi
      b  = f + f + xi
      d  = zero
      c  = f
      do 1 l=1,limit
      d  = b - d
      c  = b - one/c
         if(absc(d).lt.tm30) d = tm30
         if(absc(c).lt.tm30) c = tm30
      d = one / d
      del= d * c
      f = f * del
      b = b + w
    1 if(absc(del-one).lt.accur) go to 2
        ifail = -2
        go to 5
c
    2 xj(lmax)   = tm30
      xjp(lmax)  = f * xj(lmax)
c
c *** downward recursion to l=0 (n.b.  coulomb functions)
c
      do 3 l = lmax-1,0,-1
      xj(l) = pl*xj(l+1) + xjp(l+1)
      xjp(l)= pl*xj(l)   - xj(l+1)
    3 pl = pl - xi
c *** calculate the l=0 bessel functions
      xj0  = xi * sin(x)
      xh1(0) = exp(ci*x) * xi * (-ci)
      xh1p(0)= xh1(0) * (ci - xi)
c
c *** rescale xj, xjp,  converting to spherical bessels.
c *** recur   xh1,xh1p             as spherical bessels.
c
        w = one/xj(0)
         pl = xi
      do 4  l = 0,lmax
      xj(l)  =  xj0*(w*xj(l))
      xjp(l) =  xj0*(w*xjp(l)) - xi*xj(l)
         if(l.eq.0) go to 4
      xh1(l) = (pl-xi) * xh1(l-1) - xh1p(l-1)
         pl = pl + xi
      xh1p(l)=- pl     * xh1(l)   + xh1(l-1)
   4  continue
      ifail = 0
      return
    5      write(6,10) ifail
   10      format( 'sbesjh : ifail = ',i4)
      return
      end


c*******************************************************************
      subroutine cpolint(xa,ya,n,x,y,dy)
c ********************************************************************
c given arrays xa and ya, each of length n, and given a value x, this
c routine returns a value y, and an error estimate dy
c numerical recipes - w. h. press et al
c *********************************************************************

      implicit real*8(a-h,o-z)
      complex*16 y,ya,dy,c,d,w,denc
      parameter (nmax=50)
      dimension xa(*),ya(*),c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if(dift.lt.dif)then
         ns =i
         dif=dift
        end if
        c(i) = ya(i)
        d(i) = ya(i)
 11    continue
      y = ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho = xa(i)-x
          hp = xa(i+m)-x
          w = c(i+1) - d(i)
          den = ho - hp
          if(den.eq.0)stop
          denc= w/den
          d(i) = hp*denc
          c(i) = ho*denc
 12      continue
        if (2*ns.lt.n-m)then
          dy = c(ns+1)
          else
          dy = d(ns)
          ns = ns-1
        end if
        y = y+dy
 13    continue
      return
      end




c     ****************************************************************
!!$      subroutine factor
!!$        use factorials
!!$      implicit real*8(a-h,o-z)
!!$      faclog(1)=0.0d0
!!$      faclog(2)=0.0d0
!!$      fn=1.0d0
!!$      do 200 i=3,500
!!$      fn=fn+1.0d0
!!$  200 faclog(i)=faclog(i-1)+dlog(fn)
!!$      return
!!$      end
c     ****************************************************************

      function clebz(rl1,rl2,rl3)
        use factorials
        implicit real*8(a-h,o-z)
        l1cl=2*rl1+0.001
        l2cl=2*rl2+0.001
        l3cl=2*rl3+0.001
        lsum=(l1cl+l2cl+l3cl)/2
        lg=lsum/2
        if(2*lg-lsum) 15,20,20
15      clebz=0.0
        go to 60
20      l1=lsum-l1cl
        l2=lsum-l2cl
        l3=lsum-l3cl
        if(min0(l1,l2,l3)) 15,22,22
22      fact1=sqrt(float(l3cl+1))
        lx=l3cl/2+lg
        if(2*(lx/2)-lx) 30,32,32
30      fact1=-fact1
32      l12=l1/2+1
        l22=l2/2+1
        l32=l3/2+1
        h1=0.5*(flog(l1+1)+flog(l2+1)+flog(l3+1)-flog(lsum+2))+
     &  flog(lg+1)-flog(l12)-flog(l22)-flog(l32)
        clebz=fact1*exp(h1)
   60 return
      end


c     ******************************************************************
      real*8 function dmat(j,m,n,beta)
        use factorials
      implicit real*8(a-h,o-z)
c     ------------------------------------------------------------------
c     small d rotation matrices for integer j. brink and satchler
c     formula. rotation of beta about y-axis.
c     ------------------------------------------------------------------
      integer t
      
      cb2=cos(beta/2.)
      sb2=sin(beta/2.)
      t=0
      dmat=1.d0
      if(iabs(m).gt.j.or.iabs(n).gt.j) go to 888
      if(j.eq.0) return
      ind1=j+m+1
      ind2=j-n+1
      ind3=n-m+1
      rnum=(flog(ind1)+flog(j-m+1)+flog(ind2)+flog(j+n+1))
      rnum=rnum/2.d0
      dmat=0.d0
      go to 997
   88 den=flog(ind1)+flog(ind2)+flog(ind3)+flog(t+1)
      dmat=dmat+(-1.d0)**t*dexp(rnum-den)*cb2**(2*(j-t)+m-n)
     1     *sb2**(2*t+n-m)
   99 t=t+1
      ind1=ind1-1
      ind2=ind2-1
      ind3=ind3+1
  997 if((ind1.le.0).or.(ind2.le.0)) return
      if(ind3.le.0) go to 99
      go to 88
  888 dmat=0.d0
      return
      end
