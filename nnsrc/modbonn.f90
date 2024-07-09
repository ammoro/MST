c Global modules shared by the subroutines associated
c to the calculation of the Bonn NN potential


!********************************************************
      module crdwrt
!      common /crdwrt/ kread,kwrite,kpunch,kda(9)
        integer :: kread,kwrite,kpunch, kda(9) 
      end module crdwrt

!********************************************************
      module cpot
        real*8 :: v(6),xmev,ymev
      end module cpot


!********************************************************
      module cstate
        character *4 :: label
        integer :: j
        logical :: heform,sing,trip,coup,endep
      end module cstate

!********************************************************
      module cob
      integer ::  ic(10,15),ift(3),mint(3),maxt(3),nt,
     1            mge,mgg(12,3),mggo(12,3),ima(5,12,3),
     2            imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde     
      real*8 :: vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
     1          wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2          ez1,ez2,ct(96),wt(96)

      logical :: indc(2,15),indxy,indpar(3)

      end module cob

!*********************************************************
      module cpoted
      real*8 :: q0qmev,qfmev,pmev,uanspp,wsnspp,ucnspp,udnspp,
     1          znrl,zrel,smev
      logical :: noced
      end module cpoted


!**********************************************************
      module cfuab
      real*8 :: a(5),b(5),c(5)
      end module cfuab

!**********************************************************
!      module clebma
!      real*8 :: faclog(500)
!      end module clebma

!**********************************************************
!     
      module dmats
      real*8::sc,dc
      end module dmats

