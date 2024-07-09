      MODULE realloc_mod
      CONTAINS
        FUNCTION reallocate(p, n)     
          implicit none
          REAL*8, POINTER, DIMENSION(:) :: p, reallocate
          INTEGER, intent(in) :: n
          INTEGER :: nold, ierr
          ALLOCATE(reallocate(1:n), STAT=ierr)
          IF(ierr /= 0) STOP "allocate error"
          IF(.NOT. ASSOCIATED(p)) RETURN
          nold = MIN(SIZE(p), n)
          reallocate(1:nold) = p(1:nold)
          DEALLOCATE(p) 
        END FUNCTION REALLOCATE
      
        FUNCTION reallocate2(p, n1,n2)               ! reallocate REAL
          implicit none
          REAL*8, POINTER, DIMENSION(:,:):: p,reallocate2
          INTEGER, intent(in) :: n1,n2
          INTEGER :: n1old,n2old, ierr,lb1,lb2
          lb1=lbound(p,1) !lower bound for 1st argument
          lb2=lbound(p,2) !idem 2nd argument
          
          ALLOCATE(reallocate2(lb1:n1,lb2:n2), STAT=ierr)
          IF(ierr /= 0) STOP "allocate error"
          IF(.NOT. ASSOCIATED(p)) RETURN
          n1old = MIN(SIZE(p,1), n1)
          n2old = MIN(SIZE(p,2), n2)
          reallocate2(lb1:n1old,lb2:n2old) = p(lb1:n1old,lb2:n2old)
          DEALLOCATE(p) 
        END FUNCTION REALLOCATE2


c ***   reallocate pointer with 3 arguments
         FUNCTION reallocate3i(p, n1,n2,n3) ! reallocate INT
           implicit none
           integer, POINTER, DIMENSION(:,:,:)::p,reallocate3i
           INTEGER, intent(in) :: n1,n2,n3
           INTEGER :: n1old, n2old,n3old
           INTEGER :: ierr,lb1,lb2,lb3
           lb1=lbound(p,1) !lower bound for 1st argument
           lb2=lbound(p,2) !idem 2nd argument
           lb3=lbound(p,3) !idem 3rd argument          
           ALLOCATE(reallocate3i(lb1:n1,lb2:n2,lb3:n3),STAT=ierr)
           IF(ierr /= 0) STOP "allocate3 error"
           IF(.NOT. ASSOCIATED(p)) RETURN
           n1old = MIN(SIZE(p,1), n1)
           n2old = MIN(SIZE(p,2), n2)
           n3old = MIN(SIZE(p,3), n3)
           reallocate3i(lb1:n1old,lb2:n2old,lb3:n3old) = 
     &       p(lb1:n1old,lb2:n2old,lb3:n3old)
           DEALLOCATE(p) 
        END FUNCTION REALLOCATE3i


c ***   reallocate pointer with 3 arguments
         FUNCTION reallocate3(p, n1,n2,n3) ! reallocate REAL
           implicit none
           REAL*8, POINTER, DIMENSION(:,:,:)::p,reallocate3
           INTEGER, intent(in) :: n1,n2,n3
           INTEGER :: n1old, n2old,n3old
           INTEGER :: ierr,lb1,lb2,lb3
           lb1=lbound(p,1) !lower bound for 1st argument
           lb2=lbound(p,2) !idem 2nd argument
           lb3=lbound(p,3) !idem 3rd argument          
           ALLOCATE(reallocate3(lb1:n1,lb2:n2,lb3:n3),STAT=ierr)
           IF(ierr /= 0) STOP "allocate3 error"
           IF(.NOT. ASSOCIATED(p)) RETURN
           n1old = MIN(SIZE(p,1), n1)
           n2old = MIN(SIZE(p,2), n2)
           n3old = MIN(SIZE(p,3), n3)
           reallocate3(lb1:n1old,lb2:n2old,lb3:n3old) = 
     &       p(lb1:n1old,lb2:n2old,lb3:n3old)
           DEALLOCATE(p) 
        END FUNCTION REALLOCATE3

c ***  
         FUNCTION reallocate5(p, n1,n2,n3,n4,n5) 
          implicit none
          REAL*8, POINTER, DIMENSION(:,:,:,:,:):: p,reallocate5
          INTEGER, intent(in) :: n1,n2,n3,n4,n5
          INTEGER :: n1old, n2old,n3old,n4old,n5old
          INTEGER :: ierr,lb1,lb2,lb3,lb4,lb5
          lb1=lbound(p,1) !lower bound for 1st argument
          lb2=lbound(p,2) !idem 2nd argument
          lb3=lbound(p,3) !idem 3rd argument
          lb4=lbound(p,4) !idem 4st argument
          lb5=lbound(p,5) !idem 5st argument
          
          ALLOCATE(reallocate5(lb1:n1,lb2:n2,lb3:n3,lb4:n4,lb5:n5),
     &      STAT=ierr)
          IF(ierr /= 0) STOP "allocate5 error"
          IF(.NOT. ASSOCIATED(p)) RETURN
          n1old = MIN(SIZE(p,1), n1)
          n2old = MIN(SIZE(p,2), n2)
          n3old = MIN(SIZE(p,3), n3)
          n4old = MIN(SIZE(p,4), n4)
          n5old = MIN(SIZE(p,5), n5)
          reallocate5(lb1:n1old,lb2:n2old,lb3:n3old,
     &       lb4:n4old,lb5:n5old)= 
     &       p(lb1:n1old,lb2:n2old,lb3:n3old,lb4:n4old,lb5:n5old)
          DEALLOCATE(p) 
        END FUNCTION REALLOCATE5

        FUNCTION reallocate5c(p, n1,n2,n3,n4,n5) 
          implicit none
          complex*16, POINTER, DIMENSION(:,:,:,:,:):: p,reallocate5c
          INTEGER, intent(in) :: n1,n2,n3,n4,n5
          INTEGER :: n1old, n2old,n3old,n4old,n5old
          INTEGER :: ierr,lb1,lb2,lb3,lb4,lb5
          lb1=lbound(p,1) !lower bound for 1st argument
          lb2=lbound(p,2) !idem 2nd argument
          lb3=lbound(p,3) !idem 3rd argument
          lb4=lbound(p,4) !idem 4st argument
          lb5=lbound(p,5) !idem 5st argument
          
          ALLOCATE(reallocate5c(lb1:n1,lb2:n2,lb3:n3,lb4:n4,lb5:n5),
     &      STAT=ierr)
          IF(ierr /= 0) STOP "allocate5c error"
          IF(.NOT. ASSOCIATED(p)) RETURN
          n1old = MIN(SIZE(p,1), n1)
          n2old = MIN(SIZE(p,2), n2)
          n3old = MIN(SIZE(p,3), n3)
          n4old = MIN(SIZE(p,4), n4)
          n5old = MIN(SIZE(p,5), n5)
          reallocate5c(lb1:n1old,lb2:n2old,lb3:n3old,
     &       lb4:n4old,lb5:n5old)= 
     &       p(lb1:n1old,lb2:n2old,lb3:n3old,lb4:n4old,lb5:n5old)
          DEALLOCATE(p) 
        END FUNCTION REALLOCATE5c        
      END MODULE realloc_mod


      module rhowf
      implicit none
      integer :: nrho
      integer, allocatable :: qn(:,:)
      real*8 :: rhomax,drho,ca
      real*8, allocatable:: rr(:),cr(:),chnorm(:,:),chknorm(:)
      real*8, allocatable,target:: wfro(:,:,:),wf0(:,:)
      real*8, allocatable ::wfk(:,:,:,:)
      end module rhowf


c *** Modules used in MST program
      module parameters
        parameter (mrxy=300, mstates=285,inel=mstates-1, ma=285)
        parameter(itytrmx=2)
!        parameter(ndel=25,mll=2)
        parameter(ndel=25)
!        parameter(mllt=4)
!       mllt maximum of partial waves = 2X maximum of bound partial waves
        parameter (mrho=200)
        parameter (mthmx=300)
!        parameter(mllmx=4) AMoro 24-10


!! AMORO (30/05/05)        parameter(isnx=1, itnx=1)
         parameter(isnx=1,itnx=1)
        parameter(ikqx=2)
        parameter(lfact=100)
      end module parameters

      module wfns
	use parameters
        implicit real*8(a-h,o-z)
        logical:: iflxly
        integer mllmx,mll,dmax,lxmax,lymax,nchinc
        integer:: parityf,s,nrxy,nwf,elastic,sp
        real*8:: jc,jtot,rstep,smallchan
        real*8,allocatable::wf(:,:,:,:),rv(:),ens(:)
        integer,allocatable::na(:)
        integer,pointer::tnq(:,:,:)
        integer,allocatable:: tqn(:,:,:)
        integer inelastic,itytr
        parameter(ko=6, ipc=3, nfl=7,z=0D0, elastic=1,inelastic=2) 
        character psign(3)
        data psign/'-','?','+'/
      end module wfns

      module trdens
	use parameters
c old dimensions
c      	complex*16::rhoval(ndel,0:mll,0:mll,0:mll,0:inel)
c     	real*8::rhocm(ndel,0:mll,0:inel)
c	integer::bmax(0:inel),cdmax(2,0:mll,0:inel),lvalmax(0:inel)
        complex*16,pointer::rhoval(:,:,:,:)
        real*8,pointer::rhocm(:,:)
!        integer,allocatable::bmax(:),cdmax(:,:,:),lvalmax(:)
! AMORO (July,2005):bmax, cdmax are now in valdens
        integer,allocatable::lvalmax(:)
	real*8:: stepdel
      end module trdens

      module nnamps
        use parameters
c        complex*16:: tnnu(0:isnx,-isnx:isnx,0:1,0:1,mthmx)
c        complex*16:: tauvv(mthmx,0:mllmx,-mllmx:mllmx,0:isnx,0:isnx)  
c        complex*16:: tauvc(mthmx,0:isnx,-isnx:isnx)
        complex*16,allocatable:: tnnu(:,:,:,:,:)
        complex*16,allocatable:: tauvv(:,:,:,:,:)  
        complex*16,allocatable:: tauvc(:,:,:)
      end module nnamps

      module scattering
        use parameters  
        logical:: rel
        integer:: coulmst
        real*8,allocatable::rho(:,:)  !rho(mrho,mrho)
        real*8,allocatable::qdelta(:) !qdelta(mrho)
        real*8:: tlab,k0,muna
        real*8:: mn,zp,jp,masst
        real*8:: kcore, kc0,mucore,kv,kv0,muv,s12,s14,s0
        real*8:: nu,nu12,nu14
        complex*16::zz
        complex*16,allocatable::tint(:,:) !dim=mrho x mrho
        complex*16,allocatable::tvalence(:),tcore(:) !dim=mrho
        complex*16,allocatable::rhoaux(:),rhoauxp(:) !dim=ndel
!AMORO 4/2/04
!        complex*16:: tnorm
        real*8:: tnorm
        complex*16:: tcoreaux, tvalaux, tchkaux
        complex*16:: zrho(mrho,mrho),sum
        complex*16:: rcint ! ,rc1aux(mrho),
        complex*16, allocatable::rc1aux(:)
        real*8, allocatable:: qaux(:) !dim=mrho
        real*8:: imta0,imtb0
        real*8:: alpha,beta,sigtot
        real*8:: m234,m1234,m23,m34,w12,w14
        real*8:: m1,m2,m3,m4
        real*8:: z1,z2,z3,z4
!        real*8:: j2,j3,j4        
        real*8:: qmaxr,quin
        integer:: mquadi,mquado,iqmx
        integer:: mrquadi,mrquado
c       these dimensions are limited by gauss3 subroutine
        real*8:: radxis(201), radwt(201)
        real*8:: radxisd(201), radwtd(201)
        real*8:: radxisr(201),radwtr(201)
        integer:: irho,itydn9,itydn11,itycmdn
        real*8:: s1,s2,s3
        real*8:: thmin,thmax,dth
        integer:: issccm,ihalo,ival,icore,ihvdens
        integer:: itnn,itnnav
        real*8 qmxnn
        integer mq,mdelta
        integer nangles        
!        complex*16 ffval(mrho,0:mll,0:isnx,0:mll,0:inel)
!        complex*16 ffvth(mthmx,0:mll,0:isnx,0:mll,0:inel)
        complex*16,pointer,dimension(:,:,:,:,:):: ffval,ffvth
!        real*8 ffcore(mrho,0:mll,0:inel)
!        real*8 ffcth(mthmx,0:mll,0:inel)
        real*8, pointer,dimension(:,:,:)::ffcore,ffcth
        real*8, allocatable ::dtheta(:,:,:)
!        real*8 dtheta(mthmx,0:mllmx,-mllmx:mllmx)
!        complex*16 :: tnath(mthmx),ssctna(mrho),
        complex*16,allocatable::ssctna(:),tnath(:)
        real*8   xscin(mthmx,0:mstates),  xsconv(mthmx,0:mstates)
        real*8   xstin(mthmx,0:mstates)
!        real*8   doublec(mthmx,mstates),doublev(mthmx,mstates)
!        real*8   doublex(mthmx,mstates), doublet(mthmx,mstates)
        real*8, allocatable:: doublet(:,:)
       
        integer:: jnncb !jnnpin(8),partyin(8),istatcb(8)
        integer, allocatable::jnnpin(:),partyin(:),istatcb(:)
        integer  nwfmax,nustates,quais(100)
        integer k0000,k1100,k0111,k1120,k1121,k1122 
        integer kapa(0:1,0:1,0:2,0:2)
      end module scattering


      
      module bt
        integer:: nq
        real*8,allocatable:: qnnr(:)
        complex*16,allocatable::tnnr(:,:),tsscore(:)
      end module bt

      module bfst
        complex*16,allocatable::fst(:,:)
      end module bfst

      module bq
        real*8:: qmax,dq
      end module bq

      module time
        integer(4) cumtime1,cumtime2
      end module time


