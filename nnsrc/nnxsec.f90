c *** Calculattes x-section for NN scattering
c *** using previously calculated on-shell amplitudes
       subroutine xsec
         use ampon
         use ampnl
         implicit none
         integer:: ith
         real*8::xsect,dth,th
         if (.not.allocated(aon)) then
            write(*,*) 'Cannot calculate xsec'
            write(*,*) 'Set aeon>0'
            stop
         endif
         open(20,file='nnxsec.out',status='unknown')
         dth=180./(nth-1)
         do ith=1,nth
            th=(ith-1)*dth
            xsect=abs(aon(ith))**2 + abs(con(ith))**2
            xsect=xsect*10 ! convert fm^2 -> mb
            write(20,*) th,xsect
         enddo
         close(20)
         return
       end subroutine xsec
