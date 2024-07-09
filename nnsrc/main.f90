      program main
        parameter (lfact=100)
        logical::ifmst=.false.
        integer ::fkq
        real*8::tcm,sqmax,bqmax

c     Trace output for debugging
      open(unit=99,file='ampall.log',status='unknown')
        
        tcm=0.
        sqmax=0.
        bqmax=0.
        fkq=-1

        write(99,*)'Calling logfac with lfact=',lfact
        call logfac(lfact)
!        call setconstants()
        write(99,*)'NN main.f90: tcm=',tcm
        call ampnn(tcm,ifmst,sqmax,bqmax,fkq)
      end program main

************************************************************
*  END OF MAIN PROGRAM
************************************************************





