      program msomain_
        implicit real*8 (a-h, o-z)
        real*8:: elab=0.0
        logical :: ifmst=.false.

	open(unit=99,file='mso.log')
        call lptps(elab,ifmst)
      end program msomain_
