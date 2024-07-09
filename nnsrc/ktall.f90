c     K to T program
c     E.F. Redish, K. Stricker-Bauer (1987)
c     Formalism described in: Phys. Rev. C36, 513 (1987)
*****************************************************************
c
c          Routine to read the partial wave k matrix from unit 12
c          and write the partial wave t matrix diagonal in l on
c          unit 14.

*****************************************************************
*   f90 version
*   R.Crespo, A.M.Moro (2001)
*
*   Main changes:
*   - kmx, nmx parameters no longer used 
*   - No requires reading "kout.dat" file
*   - common blocks replaced by modules 
*****************************************************************


      subroutine kt !(lmu)
        use toff
        use kl
        use parms
        use jat1
        use ampold
        use constants
        implicit none

      integer :: norder,llp,k1,k1c,nmax,nmx,m,n,j,k,nn,ka,kb,i,
     $  nk,nkp !,lmu
      real*8 :: tcm, alf
c      real*8, allocatable ::  xx(:),wx(:)
      complex, allocatable :: tl(:),tkm(:,:),
     $  vm(:,:),tkh(:,:,:),ttnl(:,:,:,:),ttln(:,:,:,:)
      complex*16 :: zi,d(2,2),cint2d,ak,bk,ck,den,tlonp,
     $  vlon(2,2),tkon(2,2),tlonn(2,2),tlon

      write(99,*) 'Entering kt...'
      norder=5
      zi=cmplx(0.0,1.0)

      k1=kmax+1
      k1c=ic*k1
      nmax=k1c*(k1c+1)/2
      nmx=2*kmax*(2*kmax+1)/2

! Memory allocation
      write(99,*)'+Allocating memory for ttl,tl,tkm,vm,tkh...'
      if (.not.allocated(ttl)) allocate(ttl(kmax,kmax))
      allocate(tl(nmx))
      allocate(tkm(kmax,kmax))
      allocate(vm(kmax,kmax))
      allocate(tkh(2,2,kmax))
      allocate(ttnl(2,2,kmax,kmax))
      allocate(ttln(2,2,kmax,kmax))
      if (.not.allocated(xx)) allocate(xx(kmax))
      if (.not.allocated(wx)) allocate(wx(kmax))


       xx=qq(2:kmax)
       wx=wq(2:kmax)

      write(14,1000)is,ll,jj,it,ic

      if (is.lt.0) return
c      write(16,1000)is,ll,jj,it,ic
 4    format(' p=',f6.4,3x,' s,l,l",j =',4i3)
      do 20 m=1,2
      do 20 n=1,2
      tkon(m,n)=0.
      vlon(m,n)=0.
      do 20 j=1,kmax
      tkh(m,n,j)=0.
 20   continue
 15   continue
      tcm=p*p/2./redm*ch*ch

      alf=2.*redm*p/ch/ch
c         construct half shell k`s and matrix for k-->t
c         calculate on shell v, k and t
      do 10 m=1,ic
      do 10 n=1,ic
      k=nn(m,n,0,0,k1c)
      tkon(m,n)=tkl(k)*ch
      vlon(m,n)=tvl(k)*ch

      do 10 j=1,kmax
      k=nn(m,n,0,j,k1c)
      tkh(m,n,j)=tkl(k)*ch
 10   continue

      ak=1.+tkon(2,2)*zi*alf
      bk=1.+tkon(1,1)*zi*alf
      ck=-tkon(1,2)*zi*alf
      den=ak*bk-ck*ck
      d(1,1)=ak/den
      d(2,2)=bk/den
      d(1,2)=ck/den
      d(2,1)=d(1,2)

c          begin loop for diagonal pieces of coupled states
      do 25 m=1,ic
      do 25 n=1,ic
      tlon=tkon(m,n)
      do 45 ka=1,2
      do 45 kb=1,2
      tlon=tlon-zi*alf*tkon(m,ka)*d(ka,kb)*tkon(kb,n)
 45   continue
      tlonn(m,n)=tlon

c          begin loops to calculate tl on original mesh
      do 30 i=1,kmax
      do 30 j=i,kmax
      k=nn(m,n,i,j,k1c)
      vm(i,j)=tvl(k)*ch
      vm(j,i)=vm(i,j)
      tkm(i,j)=tkl(k)*ch
      tkm(j,i)=tkm(i,j)
      ttl(i,j)=tkm(i,j)
      do 35 ka=1,2
      do 35 kb=1,2
      ttl(i,j)=ttl(i,j)-zi*alf*tkh(ka,m,i)*d(ka,kb)*tkh(kb,n,j)
 35   continue
      ttl(j,i)=ttl(i,j)
      ttln(m,n,i,j)=ttl(i,j)
      ttln(m,n,j,i)=ttl(j,i)
c      write(99,*)'ttl',ttl(i,j)
 30   continue
 25   continue
      
      nmax=kmax*(kmax+1)/2
      do 70 i=1,kmax
      do 70 j=i,kmax
      ttln(1,2,j,i)=ttln(2,1,i,j)
      ttln(2,1,j,i)=ttln(1,2,i,j)
 70   continue
      do 80 m=1,ic
      do 80 n=1,ic
        if(m.eq.n)then
        
        do 40 i=1,kmax
        do 40 j=i,kmax
        nk=(i-1)*(2*kmax-i)/2+j
        tl(nk)=ttln(m,n,i,j)
!        if ((ic.eq.2).and.(m.eq.1)) dm1(jj+2,nk)=tl(nk)
!        else  dp1(jj,:)=tl(nk)
 40     continue
        write(14,1010)(tl(nk),nk=1,nmax)
        else
        write(14,1010)((ttln(m,n,nk,nkp),nk=1,kmax),nkp=1,kmax)        
        endif
 80   continue

      if (is.eq.0) then
         d01(jj+1,:)=tl(:)
      else
         if (ic.eq.1) then
            if (jj.gt.0) then
               dz1(jj+1,:)=tl(:)
            else
               dm1(jj+2,:)=tl(:)
            endif
         else
            do i=1, kmax
               do j=1,kmax
                  nk=(i-1)*(2*kmax-i)/2+j
                  dm1(jj+2,nk)=ttln(1,1,i,j)
                  dp1(jj,nk)=ttln(2,2,i,j)
                  dmm1(jj,i,j)=ttln(1,2,i,j)
                  dpp1(jj+2,i,j)=ttln(2,1,i,j)
               enddo
            enddo
         endif
      end if
      

c *** check program by interpolating to get on shell tl`s
cc      do 120 m=1,ic
cc     do 120 n=1,ic
cc        do 130 i=1,kmax
cc       do 130 j=1,kmax
c       ttl(i,j)=ttln(m,n,i,j)
c 130    continue
c      tlonp=cint2d(xx,xx,ttl,p,p,kmax,kmax,norder,kmax)
cc      write(16,44) m,n
cc      write(16,55) tlonn(m,n)
cc      write(16,66) tlonp
c 120  continue
 44    format(' m=',i2,'  n=',i2/7x,' t',12x)
 55    format(' calc',4e14.5)
 66    format(' int',1x,4e14.5/)
c      go to 100
      write(99,*)'Leaving ktall...'
      return
 90   write(14,1020)
 1000 format(6i5)
 1010 format(8e14.6)
 1020 format('  -1')
      stop
      end

*************************************************************
      function nn(m,n,i,j,k1c)
      k=m+n
      go to (10,20,30,40,10) k
 20   nn=i*(2*k1c-i-1)/2+j+1
      return
 30   if(m.gt.n)go to 35
      nn=i*(2*k1c-i-1)/2+j+k1c/2+1
      return
 35   nn=j*(2*k1c-j-1)/2+i+k1c/2+1
      return
 40   nn=(k1c/2+i)*(3*k1c/2-i-1)/2+k1c/2+j+1
      return
 10   print 9
 9    format(' error in nn')
      return
      end



