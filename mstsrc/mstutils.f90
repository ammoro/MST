      subroutine  write_time_stat()
        use time
        implicit none
        write(*,*)'Time statistics:'
        if (cumtime1<60)then
           write(*,*)'  rho ->{x,y}:',cumtime1,' secs'
           write(*,*)'  Densities  :',cumtime2,' secs'
        else
            write(*,*)'  rho ->{x,y}:',cumtime1/60.,' min'
           write(*,*)'  Densities  :',cumtime2/60.,' min'
        endif
      end subroutine write_time_stat



****************************************************************
!     Normalization constant appearing in HH
****************************************************************
        FUNCTION znorm(kk,l1,l2)
          use factorials
          IMPLICIT NONE
          
          INTEGER :: kk,l1,l2,n
          REAL*8 :: xn,xk,xlp,xln,znorm,ffactor
          xk=DBLE(kk)
          xlp=DBLE(l1)
          xln=DBLE(l2)
          n=(kk-l1-l2)/2
          xn=DBLE(n)

c          write(*,*)'znorm: k,lx,ly=',kk,l1,l2

          znorm=dsqrt(2.d0*(2.d0 + xk)*
     &       ffactor(xn)*
     &       ffactor(1d0 + xlp + xln+xn)/
     &       ffactor(0.5d0 + xlp + xn)/
     &       ffactor(0.5d0 + xln + xn))

c         write(*,*)'znorm=',znorm
        
        RETURN
        END


**************************************************************
!     Functions C(k,l1,l2,alpha) appearing in the definition of
!     Hyperspherical Harmonics 
!     USES "RUSSIAN" CONVENTION, SO DEFINITION IS CONSISTENT
!     WITH ALPHA=ATAN(X/Y), NOT ALPHA=ATAN(Y/X)
**************************************************************
      function wfalpha(k,lx,ly,ax)
      use factorials
      implicit none
      integer::k,lx,ly
      real*8::wfalpha,normhp,pojac,ax,polj,jacobi
      polj=jacobi(dcos(2*ax),(k-lx-ly)/2,lx+0.5d0,ly+0.5d0)
!!$      wfalpha=normhp((k-lx-ly)/2,dble(lx)+0.5d0,dble(ly)+0.5d0)*
!!$     ccos(ax)**dble(ly+1)*sin(ax)**dble(lx+1)*
!!$     cpojac((k-lx-ly)/2,dble(lx)+0.5d0,dble(ly)+0.5d0,cos(2.d0*ax))

      wfalpha=normhp((k-lx-ly)/2,dble(lx)+0.5d0,dble(ly)+0.5d0)*
     ccos(ax)**dble(ly+1)*sin(ax)**dble(lx+1)*polj
      end function wfalpha


****************************************************************
!     Normalization constant appearing in HH
****************************************************************
      function normhp(n,al,be) 
       use factorials
       implicit none
       real*8::normhp,ffact,al,be,ffactor
       integer::n
!!$       normhp=sqrt((2.d0*ffact(dble(n))*(dble(2*n+al+be)+1.d0)
!!$     c *ffact(dble(n+al+be)))/
!!$     c (ffact(dble(n)+be)*ffact(dble(n)+al)))
         normhp=sqrt((2.d0*ffactor(dble(n))*(dble(2*n+al+be)+1.d0)
     c *ffactor(dble(n+al+be)))/
     c (ffactor(dble(n)+be)*ffactor(dble(n)+al)))
       end function normhp





c================================================================
c     Calculates factorial of znum
c================================================================	      
      real*8 function ffactor(znum)
      real*8 agamma
      real*8 znum
      ffactor=agamma(znum+1.d0,0.d0)
!      print*,'En ffact:x=',znum,' x!=',ffactor
      return
      end


c=================================================================
c    Complex gamma function (borrowed from FRESCO)
c=================================================================
      FUNCTION CLGAMM (Z)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 S,U,V,ZP,Z,CLGAMM

      DIMENSION A(6)
      DATA C,A,D/4.1893853321D-1,-1.3346722722D-2,8.4175084173D-4,     -
     15.9523809524D-4,7.9365079365D-4,-2.777777777D-3,8.33333333D-2,0D0/

      ZRE = Z
      I=9- int(ZRE)
      IF(I.LT.0)I=0
      Z = Z +CMPLX( real(I),D)
      ZP=Z**2
      V=(1.,0.)
      S=(0.,0.)
      DO 1 J=1,6
      S=V*A(J)+S
    1 V=V*ZP
      ZP=S*Z/V
      U=Z-(0.5,0.)
      CLGAMM = LOG(Z)*U-U+ZP+C
      IF(I.EQ.0)RETURN
      DO 2 J=1,I
      Z=Z-(1.,0.)
    2 CLGAMM=CLGAMM - LOG(Z)
      RETURN
      END


c=================================================================
c     GAMMA
C=================================================================
      block data
      IMPLICIT double precision (A-H,O-Z)
      COMMON/GAM/CK(26)
C
C
      DATA CK/1.0D0   ,5.772156649015329D-1,-6.558780715202538D-1,-4.200
     1263503409520D-2, 1.665386113822915D-1,-4.219773455554430D-2,-9.621
     2971527877D-3, 7.218943246663D-3,-1.1651675918591D-3, -2.1524167411
     349D-4,  1.280502823882D-4, -2.01348547807D-5, -1.2504934821D-6, 1.
     4133027232D-6, -2.056338417D-7,  6.116095D-9,  5.0020075D-9, -1.181
     52746D-9,  1.043427D-10,  7.7823D-12, -3.6968D-12,5.1D-13,2.06D-14,
     65.4D-15,1.4D-15,1.0D-16/
      end

      
c*----------------------------------------*
c*     CALCULATES FUNCTION GAMMA(A+IB)    *
c*----------------------------------------*
c
      double precision	function agamma(a,b)
      implicit double precision(a-h,o-z)
      common/gam/ck(26)
      complex z1,z2,cgamma,z3
      PI=DACOS(-1.0D0)
      X=A
      Y=B
      ACM=1.0D0
      ACA=0.0D0
      IF(X.GT.0.0D0)THEN
      IFIN=INT(X+0.5D0)
      DO 100 I=1,IFIN
      IF(DABS(X-1.0D0).LT.1.0D-13)GOTO 100
      X=X-1.0D0
      R0=DSQRT(X*X+Y*Y)
      TET0=DPARG(Y,X)
      ACM=ACM*R0
      ACA=ACA+TET0
 100  CONTINUE
      ELSE
      IFIN=INT(DABS(X)+0.5d0)
      DO 200 I=1,IFIN
      IF(X.EQ.-1.0D0)GOTO 200
      R0=DSQRT(X*X+Y*Y)
      TET0=DPARG(Y,X)
      ACM=ACM/R0
      ACA=ACA-TET0
      X=X+1.0D0
 200  CONTINUE
      ENDIF
      IF(DABS(X-1.0D00).LE.1.0D-13)THEN
      IF(DABS(Y).LE.1.0D-13)THEN
      R1=1.0D0
      TET1=0.0D0
      ELSE
      CALL GAMM2(0.0D0,Y,R,TET)
      R1=Y*R
      TET1=TET+PI/2.0D0
      ENDIF
      ENDIF
C     IF(X.EQ.-1.0D0)THEN
C     CALL GAMM2(0.0D0,Y,R,TET)
C     R1=R/DSQRT(1.0D0+Y*Y)
C     TET1=TET-DPARG(Y,-1.0D0)
C     ENDIF
      IF(DABS(X-1.0D0).GT.1.0D-13)THEN
      CALL GAMM2(X,Y,R1,TET1)
      ENDIF
      TET=ACA+TET1
      R=ACM*R1
	agamma=r
 300  CONTINUE
      RETURN
      END
C
      SUBROUTINE GAMM2(A,B,R,TET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/GAM/CK(26)
      R0=DSQRT(A*A+B*B)
      TET0=DPARG(B,A)
      SUMR=0.0D0
      SUMY=0.0D0
      DO 100 K=1,26
      CN=CK(K)*(R0**K)
      CR=DCOS(TET0*K)
      CY=DSIN(TET0*K)
      SUMR=SUMR+CN*CR
      SUMY=SUMY+CN*CY
 100  CONTINUE
      TET=DPARG(-SUMY,SUMR)
      R=1.0D0/DSQRT(SUMR*SUMR+SUMY*SUMY)
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DPARG2(Y,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      A1=DATAN2(Y,X)
      DPARG2=A1
      RETURN
      END
C
********************* D P A R G 2 *************************
C
      DOUBLE PRECISION FUNCTION DPARG(Y,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PI=DACOS(-1.0D0)
      IF(DABS(R).LE.1.0D-16)THEN
      IF(Y.GT.0.0D0)DPARG=PI/2.0D0
      IF(Y.LT.0.0D0)DPARG=3*PI/2.0D0
      ELSE
      DPARG=DATAN(DABS(Y)/DABS(R))
      ZERO=0.0D0
      IF(Y.GE.ZERO.AND.R.LT.ZERO)THEN
      DPARG=PI-DPARG
      ENDIF
      IF(Y.LT.ZERO.AND.R.LT.ZERO)THEN
      DPARG=PI+DPARG
      ENDIF
      IF(Y.LT.ZERO.AND.R.GT.ZERO)THEN
      DPARG=-DPARG
      ENDIF
      ENDIF
      RETURN
      END



!!$*******************************************************************
!!$C     ALGORITHM 404 COLLECTED ALGORITHMS FROM ACM.
!!$C     ALGORITHM APPEARED IN COMM. ACM, VOL. 14, NO. 01,
!!$C     P. 048.
!!$      FUNCTION CGAMMA(Z)
!!$      implicit real*8 (a-h,o-z)
!!$      COMPLEX*16 Z,ZM,T,TT,SUM,TERM,DEN,CGAMMA,PI,A
!!$      DIMENSION C(12)
!!$      LOGICAL REFLEK
!!$C SET IOUT FOR PROPER OUTPUT CHANNEL OF COMPUTER SYSTEM FOR
!!$C ERROR MESSAGES
!!$      IOUT = 3
!!$      PI = (3.14159265358979311,0.d0)
!!$      zero=0.d0
!!$      X = REAL(Z)
!!$      Y = AIMAG(Z)
!!$C TOL = LIMIT OF PRECISION OF COMPUTER SYSTEM IN SINGLE PRECISI
!!$      TOL = 1.0E-7
!!$      REFLEK = .TRUE.
!!$C DETERMINE WHETHER Z IS TOO CLOSE TO A POLE
!!$C CHECK WHETHER TOO CLOSE TO ORIGIN
!!$      IF(X.GE.TOL) GO TO 20
!!$C FIND THE NEAREST POLE AND COMPUTE DISTANCE TO IT
!!$      XDIST = X-INT(X-.5)
!!$      ZM = CMPLX(XDIST,Y)
!!$      IF(CDABS(ZM).GE.TOL) GO TO 10
!!$C IF Z IS TOO CLOSE TO A POLE, PRINT ERROR MESSAGE AND RETURN
!!$C WITH CGAMMA = (1.E7,0.0E0)
!!$      WRITE(IOUT,900) Z
!!$      CGAMMA = (1.E7,0.E0)
!!$      RETURN
!!$C FOR REAL(Z) NEGATIVE EMPLOY THE REFLECTION FORMULA
!!$C GAMMA(Z) = PI/(SIN(PI*Z)*GAMMA(1-Z))
!!$C AND COMPUTE GAMMA(1-Z).  NOTE REFLEK IS A TAG TO INDICATE THA
!!$C THIS RELATION MUST BE USED LATER.
!!$10    IF(X.GE.0.0) GO TO 20
!!$      REFLEK = .FALSE.
!!$      Z = (1.0,0.0)-Z
!!$      X = 1.0-X
!!$      Y = -Y
!!$C IF Z IS NOT TOO CLOSE TO A POLE, MAKE REAL(Z)>10 AND ARG(Z)<P
!!$20    M = 0
!!$40    IF(X.GE.10.) GO TO 50
!!$      X = X + 1.0
!!$      M = M + 1
!!$      GO TO 40
!!$50    IF(ABS(Y).LT.X) GO TO 60
!!$      X = X + 1.0
!!$      M = M + 1
!!$      GO TO 50
!!$60    T = CMPLX(X,Y)
!!$      TT = T*T
!!$      DEN = T
!!$C COEFFICIENTS IN STIRLING*S APPROXIMATION FOR LN(GAMMA(T))
!!$      C(1) = 1.d0/12.d0
!!$      C(2) = -1.d0/360.d0
!!$      C(3) = 1.d0/1260.d0
!!$      C(4) = -1.d0/1680.d0
!!$      C(5) = 1.d0/1188.d0
!!$      C(6) = -691.d0/360360.d0
!!$      C(7) = 1.d0/156.d0
!!$      C(8) = -3617.d0/122400.d0
!!$      C(9) = 43867.d0/244188.d0
!!$      C(10) = -174611.d0/125400.d0
!!$      C(11) = 77683.d0/5796.d0
!!$! ALOG is archaic form of LOG. AMORO (2002)
!!$!      SUM = (T-(.5d0,0.0d0))*CDLOG(T)-T+
!!$!     &      CMPLX(.5d0*ALOG(2.d0*3.14159265358979311),zero)
!!$      SUM = (T-(.5d0,0.0d0))*CDLOG(T)-T+
!!$     &      CMPLX(.5d0*LOG(2.d0*3.14159265358979311),zero)
!!$      J = 1
!!$70    TERM = C(J)/DEN
!!$C TEST REAL AND IMAGINARY PARTS OF LN(GAMMA(Z)) SEPARATELY FOR
!!$C CONVERGENCE.  IF Z IS REAL SKIP IMAGINARY PART OF CHECK.
!!$      IF(ABS(REAL(TERM)/REAL(SUM)).GE.TOL) GO TO 80
!!$      IF(Y.EQ.0.0) GO TO 100
!!$      IF(ABS(AIMAG(TERM)/AIMAG(SUM)).LT.TOL) GO TO 100
!!$80    SUM = SUM + TERM
!!$      J = J + 1
!!$      DEN = DEN*TT
!!$C TEST FOR NONCONVERGENCE
!!$      IF(J.EQ.12) GO TO 90
!!$      GO TO 70
!!$C STIRLING*S SERIES DID NOT CONVERGE.  PRINT ERROR MESSAGE AND
!!$C PROCEDE.
!!$90    WRITE(IOUT,910) Z
!!$C RECURSION RELATION USED TO OBTAIN LN(GAMMA(Z))
!!$C LN(GAMMA(Z)) = LN(GAMMA(Z+M)/(Z*(Z+1)*...*(Z+M-1)))
!!$C = LN(GAMMA(Z+M)-LN(Z)-LN(Z+1)-...-LN(Z+M
!!$100   IF(M.EQ.0) GO TO 120
!!$      DO 110 I = 1,M
!!$      A = CMPLX(I*1.d0-1.d0,0.0d0)
!!$110   SUM = SUM-CDLOG(Z+A)
!!$C CHECK TO SEE IF REFLECTION FORMULA SHOULD BE USED
!!$120   IF(REFLEK) GO TO 130
!!$      SUM = CDLOG(PI/CDSIN(PI*Z))-SUM
!!$      Z = (1.0d0,0.0d0) -Z
!!$130   CGAMMA = CDEXP(SUM)
!!$      RETURN
!!$900   FORMAT(1X,2E14.7,10X,49HARGUMENT OF GAMMA FUNCTION IS TOO CLOSE TO
!!$     1 A POLE)
!!$910   FORMAT(44H ERROR - STIRLING*S SERIES HAS NOT CONVERGED/14X,4HZ = ,
!!$     12E14.7)
!!$      END





c=====================================================================
c     Polinomio de Jabobi con indices n,alfa,beta
c     Usamos el metodo de recurrencia (formula 22.7.1 Abramovitz)
c=====================================================================      
      real*8 Function jacobi(x,n,al,be)
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 ::x,dk,al,be
!      real*8, allocatable :: AuxJ(:)
      integer ::n
      real*8 :: a1,a2,a3,a4
      real*8 :: aux0,aux1,aux2
c
c     Valores de partida (formula 22.4.1 Abramovitz)
c

      aux0=1d0 
      aux1=0.5D0*(al-be+(al+be+2.d0)*x)

      select case(n)
         case(0) 
            jacobi=aux0
!            return
         case(1)
            jacobi=aux1
!            return
         case default
            Do k=1,n-1
               dk=dble(k)
               a1=2.d0*(dk+1.d0)*(dk+al+be+1.d0)*(2.d0*dk+al+be)
               a2=(2.d0*dk+al+be+1.d0)*(al**2.d0-be**2.d0)
               a3=(2.d0*dk+al+be)*(2.d0*dk+al+be+1.d0)*
     &           (2.d0*dk+al+be+2.d0)
               a4=2.d0*(dk+al)*(dk+be)*(2.d0*dk+al+be+2.d0)
               
               aux2=((a2+a3*x)*aux1-a4*aux0)/a1
               !         AuxJ(k+1)=((a2+a3*x)*AuxJ(k)-a4*AuxJ(k-1))/a1
               aux0=aux1
               aux1=aux2
            End do
            jacobi=aux2
         end select
            
c      print*,'En Jacobi:xk=',x,' n=',n,' alfa=',al,' beta=',be
c      print*,'          Jacobi=',jacobi
      Return
      End
      



      function trian(x,y,z)
        real*8::x,y,z,trian
        trian=(x-y-z)**2 - 4.*y*z
      end function trian

      SUBROUTINE BESSR(L,Q,BSJ,NR,DR)
!        use wfns
***********************************************************
*     Calculates bessel functions with  L=0-4
***********************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      
      complex*16 xj(0:30),xjp(0:30),xh1(0:30),xh1p(0:30),x
      DIMENSION BSJ(NR)
!      lmx2 = 4

!! AMORO July 2005
!!!!!      lmx2=mllmx

      lmx2=l

      R=DR
      DO 100 KR=1,NR
      ZZ=Q*R
      SZ=DSIN(ZZ)
      CZ=DCOS(ZZ)
      IF(L.EQ.0) then
         if (zz.lt.1e-10)then
         BSJ(KR) = 1.d0
         else
         BSJ(KR)=SZ/ZZ
         end if
      end if
      IF(L.EQ.1) BSJ(KR)=(SZ/ZZ-CZ)/ZZ
      IF(L.EQ.2) BSJ(KR)=(SZ*(3.D0/ZZ/ZZ-1.D0)-3.D0/ZZ*CZ)/ZZ

      IF (L.GT.2)THEN
      x = zz
      call sbesjh(x,lmx2,xj,xjp,xh1,xh1p,ifail)
!! AMORO July 2005
!!      IF(L.EQ.3)BSJ(KR)= dreal(xj(3))
!!      IF(L.EQ.4)BSJ(KR)= dreal(xj(4))
      BSJ(KR)=DREAL(XJ(L))
      ENDIF

      R=R+DR

  100 CONTINUE
      RETURN
      END








      FUNCTION WIG3J(A,B,C,AL,BE,GA)
        use factorials
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL FAIL3,FRAC
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      IF(AL+BE+GA) 11,10,11
11    WIG3J = 0.0
      RETURN
10    IF(FAIL3(C,A,B)) GOTO 11
      IF(A-ABS(AL)) 11,14,14
14    IF(B-ABS(BE)) 11,15,15
15    IF(C-ABS(GA)) 11,13,13
13    IA = C-B+AL
      IB = C-A-BE
      IF(IA) 20,21,21
21    IF(IB) 20,23,23
23    MIN = 0
      GO TO 24
20    MIN = -IA
      IF(MIN+IB) 25,24,24
25    MIN = -IB
24    IC = A-AL
      ID = B+BE
      IE = A+B-C
      NIN = MIN
      T=PHASE(MIN)
      S = T
30    MIN = MIN+1
      IZ = IC-MIN+1
      IF(IZ) 29,29,26
26    IY = ID-MIN+1
      IF(IY) 29,29,27
27    IX = IE-MIN+1
      IF(IX) 29,29,28
28    TA = IX*IY*IZ
      TB = MIN*(IA+MIN)*(IB+MIN)
      T = -T*TA/TB
      S = S+T
      GO TO 30
29    IF(S.EQ.0.0) GO TO 11
      I = B-A+GA
      IH = A+AL
      II = B-BE
      IJ = C+GA
      IK = C-GA
      IL = A+C-B
      IM = B+C-A
      IN = A+B+C+1.0
      XDDD = 0.5*(FLOG(IH+1)+FLOG(IC+1)+FLOG(ID+1)
     $           +FLOG(II+1)+FLOG(IJ+1)+FLOG(IK+1)
     $           +FLOG(IE+1)+FLOG(IM+1)+FLOG(IL+1)-FLOG(IN+1))
     $      - ( FLOG(IC-NIN+1)+FLOG(IA+NIN+1)+FLOG(ID-NIN+1)+
     $          FLOG(IB+NIN+1)+FLOG(NIN+1)+FLOG(IE-NIN+1))
      WIG3J = PHASE(I)  *  EXP(XDDD) * S
      RETURN
      END

C     FUNCTION WIG6J
C
C     THIS CALCULATES WIGNER 6-J COEFFICIENTS AS DEFINED IN BRINK AND
C     SATCHLER. IT REQUIRES TWO SUBPROGRAMS, FACT AND COMP... SEE
C     LISTING OF 3-J SYMBOL. IT USES THE CLOSED FORM FOR W-COEFICIENTS
C     DUE TO RACAH ALSO IN BRINK AND SATCHLER. THE SAME RULES APPLY FOR
C     THE INPUT PARAMETERS AS STATED IN THE 3-J PROGRAM.
C
      FUNCTION WIG6J(A,B,C,D,E,F)
        use factorials
c*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, pointer, dimension(:):: FC
      LOGICAL FAIL3,FRAC
     
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)

      fc=>flog
C
      IF(FAIL3(A,B,C)) GO TO 10
      IF(FAIL3(A,E,F)) GO TO 10
      IF(FAIL3(B,D,F)) GO TO 10
      IF(FAIL3(C,D,E)) GO TO 10
      GO TO 14
   10 WIG6J=0.
      RETURN
   14 IC=A+B-C
      ID=E+D-C
      IE=A+E-F
      IG=B+D-F
      IA=C+F-A-D
      IB=C+F-B-E
      IH=A+B+E+D+1.
      M=MIN(IH,IC,ID,IE,IG)
      IF(M)10,17,17
   17 MUP=MIN(IA,IB,0)
      T=PHASE(M)
      N = M
      S=T
   18 M=M-1
      IF(M+MUP)24,16,16
   16 TA=(IA+M+1)*(IB+M+1)*(IH-M)*(M+1)
      TB=(IC-M)*(ID-M)*(IE-M)*(IG-M)
      T=-T*TA/TB
      S=S+T
      GO TO 18
   24 IT=A+B+C+1.
      IU=A+E+F+1.
      IV=B+D+F+1.
      IW=C+D+E+1.
      XD  =  .5*(FC(1+IC)+FC(1+IE+IB)+FC(1+IA+IG)
     1+FC(1+IE)+FC(1+IB+IC)+FC(1+IA+ID)+FC(1+IG)+FC(1+IC+IA)+FC(1+ID+IB)
     2       +FC(1+ID)+FC(1+IA+IE)+FC(1+IB+IG)-FC(1+IT)-FC(1+IU)-
     3      FC(1+IV)-FC(1+IW))
     4             + FC(1+IH-N)-FC(1+N)-FC(1+IA+N)-FC(1+IB+N)-FC(1+IC-N)
     5     -FC(1+ID-N)-FC(1+IE-N)-FC(1+IG-N)
      WIG6J=PHASE(IH-1)  *S*EXP(XD)
      RETURN
      END

C     FUNCTION WIG9J(A,B,C,D,E,F,G,H,Z)
C
C     THIS CALCULATES 9-J SYMBOLS, OR X-COEFICIENTS AS DEFINED IN BRINK
C     AND SATCHLER. IT USES THE FORMULA FOR 9-J SYMBOLS IN TERMS OF 6-J
C     SYMBOLS. IT THEREFORE NEEDS THE 6-J SUBPROGRAM, AND THE CONDITION
C     ON THE INPUT PARAMETERS IS THE SAME AS FOR THE 3-J PROGRAM#      0
C
       FUNCTION WIG9J(A,B,C,D,E,F,G,H,Z)
c********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL FAIL3,FRAC
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      IF(FAIL3(A,B,C)) GO TO 20
      IF(FAIL3(D,E,F)) GO TO 20
      IF(FAIL3(C,F,Z)) GO TO 20
      IF(FAIL3(A,D,G)) GO TO 20
      IF(FAIL3(B,E,H)) GO TO 20
      IF(FAIL3(G,H,Z)) GO TO 20
      GO TO 26
   20 WIG9J=0.
      RETURN
   26 S=0.
      XA=ABS(A-Z)
      XB=ABS(D-H)
      XC=ABS(B-F)
      X=XA
      IF(X-XB)10,11,11
   10 X=XB
   11 IF(X-XC)12,13,13
   12 X=XC
   13 IF(X-A-Z)14,14,15
   14 IF(X-D-H)16,16,15
   16 IF(X-B-F)17,17,15
   17 S=S+(2.*X+1.)*WIG6J(A,Z,X,H,D,G)*WIG6J(B,F,X,D,H,E)*WIG6J(A,Z,X,F,
     1B,C)
      X=X+1.
      GO TO 13
   15 IF(S-0.)18,20,18
   18 K=2.0*(A+B+D+F+H+Z)
      WIG9J=PHASE(K)*S
      RETURN
      END

      FUNCTION CLEB6(A,AL,B,BE,C,M)
        use factorials
c*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL* 8 M !,FACT
      LOGICAL FAIL3,FRAC
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      GA = -M
C
      IF(AL+BE+GA) 11,10,11
C11    WIG3J = 0.0
11    CLEB6 = 0.0
      RETURN
10    IF(FAIL3(C,A,B)) GO TO 11
      IF(A-ABS(AL)) 11,14,14
14    IF(B-ABS(BE)) 11,15,15
15    IF(C-ABS(GA)) 11,13,13
13    IA = C-B+AL
      IB = C-A-BE
      IF(IA) 20,21,21
21    IF(IB) 20,23,23
23    MIN = 0
      GO TO 24
20    MIN = -IA
      IF(MIN+IB) 25,24,24
25    MIN = -IB
24    IC = A-AL
      ID = B+BE
      IE = A+B-C
      NIN = MIN
      T=PHASE(MIN)
      S = T
30    MIN = MIN+1
      IZ = IC-MIN+1
      IF(IZ) 29,29,26
26    IY = ID-MIN+1
      IF(IY) 29,29,27
27    IX = IE-MIN+1
      IF(IX) 29,29,28
28    TA = IX*IY*IZ
      TB = MIN*(IA+MIN)*(IB+MIN)
      T = -T*TA/TB
      S = S+T
      GO TO 30
29    I = B-A+GA
      IF(S.EQ.0.0) GO TO 11
      IH = A+AL
      II = B-BE
      IJ = C+GA
      IK = C-GA
      IL = A+C-B
      IM = B+C-A
      IN = A+B+C+1.0
      XDDD = 0.5*(FLOG(IH+1)+FLOG(IC+1)+FLOG(ID+1)
     $           +FLOG(II+1)+FLOG(IJ+1)+FLOG(IK+1)
     $           +FLOG(IE+1)+FLOG(IM+1)+FLOG(IL+1)-FLOG(IN+1))
     $      - ( FLOG(IC-NIN+1)+FLOG(IA+NIN+1)+FLOG(ID-NIN+1)+
     $          FLOG(IB+NIN+1)+FLOG(NIN+1)+FLOG(IE-NIN+1))
C     WIG3J = (-1.0)**I *  EXP(XDDD) * S
      CLEB6 = SQRT(2*C+1.)*EXP(XDDD) * S
      RETURN
      END

      FUNCTION RACAH2(A,B,C,D,E,F)
      IMPLICIT REAL*8(A-H,O-Z)
      PHASE(I) = (-1)**I
      Z = ABS(A+B+C+D)
      I = Z + 0.5
      RACAH2 = PHASE(I) * WIG6J(A,B,E,D,C,F)
      RETURN
      END



       subroutine drotat
c***************************************************************** 
       use parameters
       use wfns
       use scattering
       use factorials
       use legendre
       use constants
       implicit real*8(a-h,o-z)
       integer :: istat=0

        NJ = mll
        MJ = mll
        if (mj>lfact) then
           write(*,*) 'drotat: insufficient dimension for fact'
           write(*,fmt='("allocated=",1i4,"but required=",1i4)')lfact,mj
           write(*,*) 'Please, increase parameter lfact'
           stop
        endif

        write(99,*)'+ Allocating dtheta with',nangles,mll
       allocate(dtheta(nangles,0:mll,-mll:mll),stat=istat)
!       if (istat>0) then
!          write(*,*)'drotat: error allocating dtheta. Aborting'
!          stop
!       endif
      
       do 200 ith = 1,nangles

       th = thmin + dth*(ith-1)
       thrad = th*pi/180.
       cthna = cos(thrad)
       call PLM(cthna,NJ,MJ)
     
       do 250 nc = 0, mll 
       do 255 ngama = 0, nc          
       plaux = PL(nc+1, ngama+1 )
       dtheta(ith,nc,ngama) = plaux * (-1)**ngama
     *    *sqrt(  exp((flog(nc-ngama)-flog(nc + ngama))) )
       if (nc.gt.0)then
       dtheta(ith,nc,- ngama)=  (-1)**ngama *  dtheta(ith,nc,ngama)      
       endif
 255   continue
 250   continue

 200   continue
       
       return
       end


       subroutine drotatold
c***************************************************************** 
       use parameters
       use wfns
       use scattering
       use factorials
       use legendre
       use constants
       implicit real*8(a-h,o-z)
       integer :: istat=0

        NJ = mllmx + 1
        MJ = 2*NJ + 1
        if (mj>lfact) then
           write(*,*) 'drotat: insufficient dimension for fact'
           write(*,fmt='("allocated=",1i4,"but required=",1i4)')lfact,mj
           write(*,*) 'Please, increase parameter lfact'
           stop
        endif

        write(99,*)'+ Allocating dtheta with',nangles,mllmx
       allocate(dtheta(nangles,0:mllmx,-mllmx:mllmx),stat=istat)
!       if (istat>0) then
!          write(*,*)'drotat: error allocating dtheta. Aborting'
!          stop
!       endif
      
        do 200 ith = 1,nangles
       th = thmin + dth*(ith-1)
       thrad = th*pi/180.
       cthna = cos(thrad)
       call PLM(cthna,NJ,MJ,NAMT)
     
       do 250 nc = 0, mllmx  
       do 255 ngama = 0, nc
          
       plaux = PL(nc+1, ngama+1 )
       dtheta(ith,nc,ngama) = plaux * (-1)**ngama
     *    *sqrt(  exp((flog(nc-ngama+1)-flog(nc + ngama+1))) )
!       write(61,*) ith,nc,ngama,dtheta(ith,nc,ngama)

       if (nc.gt.0)then
       dtheta(ith,nc,- ngama)=(-1)**ngama * (    plaux * (-1)**ngama
     *    *sqrt(  exp((flog(nc-ngama+1)-flog(nc + ngama+1))) )      )
        
       endif
 255   continue
 250   continue
 200   continue
       
       return
       end


       SUBROUTINE PLM(X,N,M,NA)
        use legendre
************************************************************************
*     Associated Legendre Polynomials
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      
      REAL*8 L,X

      if (allocated(pl)) deallocate(pl)
      
      N1 = N+1
      M1 = M+1
      allocate (pl(n1,m1))
      DO 10 J=1,M1
      DO 10 I=1,N1
10    PL(I,J)=0.
      PL(1,1) = 1.
      PL(2,1) = X
      SX = SQRT(1.-X*X)
      PL(2,2) = SX
      DO 20 I=3,N1
      L = I-1
      PL(I,1)=((2.*L-1.)*X*PL(I-1,1) - (L-1.)*PL(I-2,1))/L
      JM = MIN(I,M1)
      DO 20 J=2,JM
      PL(I,J) = (2.*L-1.)*SX*PL(I-1,J-1) + PL(I-2,J)
20    CONTINUE
      RETURN
      END


* Evaluates Rutherford x-sec
      function sigruth(z1,z2,ecm,theta)
        use constants !pi,hbarc,amu,finec
        implicit none
        real*8::q,z1,z2,ecm,theta
        real*8:: sigruth
        
        q=z1*z2*hbarc/(2.*Ecm*finec)
        if (theta<1e-6) theta=1e-6
	sigruth=q*q/(4d0*sin(theta*Pi/180d0/2d0)**4) !En fm^2
        sigruth=sigruth*10.                    !En mb
      end function sigruth

    
	
