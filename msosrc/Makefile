#Directory containing nnamp source
NNSRC=../$(NNAMP)
MSTSRC=../$(MST)

.SUFFIXES: .f90 .o

MAKEFLAGS=

OBJ=modules.o rc1.o rc2.o rc3.o rc4.o  optpot.o fform.o xsec.o 
OBJ1=modules.o rc1.o rc2.o rc3.o rc4.o   optpot.o msomain.o \
	fform.o xsec.o $(NNSRC)/modnn.o $(NNSRC)/klsp.o \
	$(NNSRC)/ktall.o $(NNSRC)/ampall.o $(NNSRC)/mon.o \
	$(NNSRC)/modbonn.o $(NNSRC)/potbonn.o $(NNSRC)/moff.o  \
	$(NNSRC)/nnmain.o $(NNSRC)/nnutils.o $(NNSRC)/utils.o
all: $(OBJ1)
	$(FC) -o mso $(OBJ1) $(LFLAGS)
	cp *mod $(MSTSRC)
#	cp tatheta.mod amps.mod j2.mod j.mod $(MSTSRC)
msomain.o:msomain.f90 modules.f90
	$(FC) $(FORM) $(FFLAGS) -c msomain.f90
nn:
	cd $(NNSRC); $(MAKE) -f Makefile objs
modules.o:modules.f90
	$(FC) $(FORM) $(FFLAGS) -c modules.f90
rc1.o:modules.o rc1.f90
	$(FC) $(FORM) $(FFLAGS) -c rc1.f90
rc2.o:modules.o rc2.f90
	$(FC) $(FORM) $(FFLAGS) -c rc2.f90
rc3.o:modules.o rc3.f90
	$(FC) $(FORM) $(FFLAGS) -c rc3.f90
rc4.o:modules.o rc4.f90 
	$(FC) $(FORM) $(FFLAGS) -c rc4.f90 
fform.o:modules.o fform.f90 
	$(FC) $(FORM) $(FFLAGS) -c fform.f90
optpot.o:modules.o optpot.f90 
	$(FC) $(FORM) $(FFLAGS) -c optpot.f90
xsec.o:modules.o xsec.f90 
	$(FC) $(FORM) $(FFLAGS) -c xsec.f90
utils.o:modules.o utils.f90
	$(FC) $(FORM) $(FFLAGS) -c utils.f90
deprecated.o:modules.o deprecated.f90 
	$(FC) $(FORM) $(FFLAGS) -c deprecated.f90

clean:
	-rm -f *.o core
	-rm -f $(NNSRC)/*.o $(NNSRC)/core
superclean: clean
	-rm -f *.mod fort.* 
	cd $(NNSRC); $(MAKE) -f Makefile superclean
# DO NOT DELETE







