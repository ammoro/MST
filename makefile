#If compiled with fujitsu f90
#include fujitsu.def

#If compiled with pgf90
#include pgf90.def

#If compiled with Intel ifc
#include ifc.def

#If compiled with Intel ifort
include ifort.def

#If compiled with gfortran
#include gfortran.def

#BINDIR=
BINDIR=$(HOME)/bin

#-------------------------------------------------------------
#Do not edit below this line unless you know what you are doing
#--------------------------------------------------------------


#Directory containing NNAMP source
NNAMP=nnsrc
#Directory containing MSO source
MSO=msosrc
#Directory containing MST source
MST=mstsrc

export FC
export FFLAGS
export FORM
export MSO
export MST
export NNAMP

.SUFFIXES: .f90 .o

MAKEFLAGS=

all: nn msoamp mstamp
#	cd $(NNAMP); $(MAKE) -f Makefile objs
mstamp: msoamp nn
	@echo ------------------------------------
	@echo Compiling MST program:
	@echo -----------------------------------
	cd $(MST); $(MAKE) -f Makefile
msoamp: nn
	@echo ------------------------------------
	@echo Compiling MSO program:
	@echo -----------------------------------
	cd $(MSO); $(MAKE) -f Makefile
nnamp:
	cd $(NNAMP); $(MAKE) -f Makefile
nn:
	@echo ------------------------------------
	@echo Compiling NN program:
	@echo -----------------------------------
	cd $(NNAMP); $(MAKE) -f Makefile objs
	cd $(NNAMP); $(MAKE) -f Makefile
install: all
#if BINDIR
#	echo Please, define BINDIR in the top makefile
#else
#	-mkdir $(BINDIR)
	 @echo 
	 @echo Installing executables on $(BINDIR) 
	 @echo -------------------------------------------------
	-cp $(NNAMP)/nnamp $(BINDIR)/
	-cp $(MSO)/mso     $(BINDIR)/
	-cp $(MST)/mst     $(BINDIR)/
	@echo ---------------------------------------------------
#endif
clean:
	-rm -f *.o core
	cd $(NNAMP); make clean
	cd $(MSO); make clean
	cd $(MST); make clean
superclean: clean
	-rm -f *.mod fort.* 
	cd $(NNAMP); $(MAKE) -f Makefile superclean
	cd $(MSO); $(MAKE) -f Makefile superclean
	cd $(MST); $(MAKE) -f Makefile superclean
# DO NOT DELETE







