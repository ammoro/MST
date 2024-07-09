	parameter (mth=200, mexc=300,mens=mexc)
	real xscin(mth,mexc),xstot(mth),xsexp(mth)
	real ens(mexc),xsens(0:mexc),xsm(0:mens,mth)
	gauss(r,res) = exp(-(r/res)**2)/(res*sqrt(pi))
	
!	open(1,file='doublexs.data')
!	open(4,file='doublexsMC.data',recl=80)
!        open(3,file='ens.test')
	pi = 4d0*atan(1d0)
	rad = pi/180.0
	
        write(*,*)'pause1'
	read(1,*) theta1,theta2,dtheta
	nth = nint((theta2-theta1)/dtheta)+1
	if(nth>mth) stop 'mth'
	        write(*,*)theta1,theta2,dtheta
	read(1,*) nexc
	if(nexc>mexc) stop 'mexc'
                write(*,*)nexc
	
c	do 30 in=1,nexc
c30	read(2,*) ens(in)	
	
	do 60 ith=1,nth
        write(*,*)'ith=',ith
	xstot(ith) = 0.
	do 60 in=1,nexc
	read(1,*) ens(in),    xscin(ith,in)
60	xstot(ith) = xstot(ith) + xscin(ith,in)

	write(*,*) 'read all doublexs.data'
	
	do 90 in=1,nexc
	write(3,70) ens(in)
70	format('#  ens =',f8.3)
	do 80 ith=1,nth
	theta = (ith-1)*dtheta+theta1
80	write(3,*) theta,xscin(ith,in)
90	write(3,*) '&'

	do 100 ith=1,nth
	theta = (ith-1)*dtheta+theta1
	write(20,*) theta,xstot(ith)
100	continue

	thex1 = 24.0
	thex2 = 45.0
	ithex1 = nint((thex1-theta1)/dtheta)+1
	ithex2 = nint((thex2-theta1)/dtheta)+1
	do 160 ith=ithex1,ithex2
	xsexp(ith) = 0.
	do 160 in=1,nexc
160	xsexp(ith) = xsexp(ith) + xscin(ith,in)
	do 200 ith=ithex1,ithex2
	theta = (ith-1)*dtheta+theta1
	write(21,*) theta,xsexp(ith)
200	continue	
	
	xst = 0.0
	do 250 in=1,nexc
	xsens(in) = 0.0
	do 240 ith=ithex1,ithex2
	theta = (ith-1)*dtheta+theta1
	xsens(in) = xsens(in) + sin(theta*rad)*xscin(ith,in)
     x		*dtheta*rad * 2*pi
240	continue
	write(22,*) ens(in),xsens(in)
	xst = xst + xsens(in)
250	continue
	write(*,*) 'Total cross section in range ',thex1,' to ',thex2,
     x		' is ',xst

	dens = 0.1
	ens1 = dens
	ens2 = 10.0
	nens = nint((ens2-ens1)/dens)+1
	if(nens>mens) stop 'mens'
	res = 0.5
	xsm(0,:) = 0.0
	do 300 iens=1,nens
	e = (iens-1)*dens + ens1
	do 290 ith=1,nth
	 t = 0.0
	 do 280 in=1,nexc
280	  t=t + xscin(ith,in)*gauss(e-ens(in),res)
290	xsm(iens,ith) = t
300	continue

	write(4,*) ens1,ens2,dens,nens
	write(4,*) theta1,theta2,dtheta,nth
	do 310 iens=0,nens
	e = (iens-1)*dens + ens1
	write(4,*)  ' e = ',e
310	write(4,*) (xsm(iens,ith),ith=1,nth)
	
	do 350 iens=0,nens
	e = (iens-1)*dens + ens1
	xsens(iens) = 0.0
	do 340 ith=ithex1,ithex2
	theta = (ith-1)*dtheta+theta1
	xsens(iens) = xsens(iens) + sin(theta*rad)*xsm(iens,ith)*
     x		dtheta*rad * 2*pi
340	continue
	write(23,*) e,xsens(iens)
350	continue
	
	stop
	end
