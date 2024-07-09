         module bwf
           real*8::grho(0:4,201)
         end module bwf

         module nlspfl
           complex*16, pointer:: u(:,:,:)
           complex*8, pointer:: uc(:,:,:)
           complex*16, allocatable,target:: ul(:,:,:,:)
           complex*8, allocatable,target :: ucoul(:,:,:,:)
         end module nlspfl

         module foptp
           real*8, allocatable,dimension(:,:):: fr,fi
           real*8, allocatable,dimension(:,:):: fcr,fci
         end module foptp

         module kinemt
           real*8 :: k0,el,ecm,muc,mnuc2,s,plab,mnuc,mp2
         end module kinemt

         module params
           integer :: nz,na,nes,nwaves
!           integer :: nes,nwaves
           real*8:: hbarc,pi,mp,mn
         end module params
           
         module switch
           integer :: nifty(20)
         end module switch

         module sizes
           real*8 :: achp,acmp,wsp,achn,acmn,wsn
         end module sizes

         module ranges
           real*8 :: rcoul,rcut
         end module ranges
         
         module inputs
           logical:: kmt,offshell,closure
           integer :: lxmax,kode,nang,ngp,nr,kin,kout,
     &                coulomb,at,zt,zp,second,
     &                ncluster
           real*8:: massp,masst,usnr(10),usni(10),ucnr(10),ucni(10)
           real*8 :: tlab,b,ymin1,ymin2
           real*8 ::ahop,ahon
           real*8:: alfa,beta,gama
           real*8:: bi,br
         end module inputs
           
         module rcomn
           real*8:: rr,ri,rrb,rib
           real*8:: rcr,rci,rcrb,rcib !coulomb (AMoro addition)
         end module rcomn
         
         module spins
           integer:: nspina,nspinb,nspinc,nspind,nspine
         end module spins

         module tcoul
           real*8:: trc , tic , tr3c , ti3c , tr5c , ti5c
         end module tcoul
           
         module tcomn
           real*8,allocatable,dimension(:,:)::tr,ti,tbr,tbi
         end module tcomn

         module tcoulb
           real*8,allocatable::  trcoul(:,:), ticoul(:,:)
         end module tcoulb
         
         module wr
           integer :: iwrite, iread,istrong
         end module wr

         module bgrid
           real*8,allocatable,target ::kk(:)
           real*8,allocatable ::wt(:)
         end module bgrid
           
         module radgr
           real*8 :: radxis(201),radwt(201)
         end module radgr
           
         module bprop
           real*8, allocatable :: deng(:,:)
         end module bprop

         module clebma
           real*8, allocatable :: faclog(:)
         end module clebma

         module bscoul
           real*8,allocatable:: scoul(:)
         end module bscoul
         
         module qnumb
           integer ::lamx,lbmx,l1mx,lfmx
         end module qnumb

         module bkcall 
           integer :: kcall
         end module bkcall

         module blocks
           integer :: option1,option2,option3,option4,option5,option6
         end module blocks

         module radwf
           integer :: igwf,lstore
           real*8:: rmaxwf,stepwf
         end module radwf

         module etapbl
           complex*16, allocatable :: etap(:,:)
         end module etapbl

         module btmatc
           complex*16, allocatable :: tmatc(:,:,:)
         end module btmatc


         
         module bgauss3
         integer :: mquadi,mquado,irmx
         real*8 ::  rmaxr,quin
         end module bgauss3

         module bsdess
         real*8, allocatable :: sbess(:,:,:)
         end module bsdess

         module halo
           integer ::ioptls,ioptcnt
         end module halo

          module ffch
          integer ::ich
         end  module ffch

         module bsflip 
           integer :: isflip
         end module bsflip

         module bion2
           integer :: ion2
         end module bion2

         module bhalodn
           integer:: itydn9,itydn11,itycmdn
         end module bhalodn
         
         module brms
         real*8 :: rms9,rms11
         end module brms

         module bwfnum
           real*8,allocatable :: wfs(:),wfp(:),densr(:)
           real*8 :: drx
           integer:: nramax
         end module bwfnum

         module bandt
           real*8:: imtij
           real*8,allocatable::bnuc(:),bn(:),bnucf(:),bnf(:),retij(:)
         end module bandt
           
         module bdenstg
           integer :: nrmxtg
           real*8,allocatable :: rtg(:),denstg(:,:)
           real*8 :: drxtg
         end module bdenstg
         


         module amphe
         complex*16, allocatable :: an(:,:),cn(:,:)
         end module amphe

         module btmatt
          complex*16, allocatable::tmatt(:,:,:,:) 
         end module btmatt

         module bnoyes  
           complex*16,allocatable::fnoyes(:,:)
         end module bnoyes

         module bcwfn
           real*8,allocatable,dimension(:)::sa,fc,fpc,gc,gpc
         end module bcwfn

         module btaux
           complex*16,allocatable:: tcwf(:,:)
         end module btaux

         module rhos
           real*8,allocatable, target :: ff(:,:,:,:,:)
           real*8, pointer :: ffp(:,:,:,:),ffn(:,:,:,:)
         end module rhos
         

         module tatheta
           integer:: noangs
           real*8,allocatable::theta(:)
           complex*16,allocatable:: tnct(:)
         end module tatheta
