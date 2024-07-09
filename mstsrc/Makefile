MSOSRC=../$(MSO)
NNSRC=../$(NNAMP)

.SUFFIXES: .c .o

OBJ1 = modmst.o mstutils.o mstmain.o readsmat.o read3bwf.o
OBJ= $(MSOSRC)/modules.o $(MSOSRC)/rc1.o $(MSOSRC)/rc2.o \
$(MSOSRC)/rc3.o $(MSOSRC)/rc4.o $(NNSRC)/utils.o  \
$(MSOSRC)/optpot.o $(MSOSRC)/fform.o $(MSOSRC)/xsec.o $(NNSRC)/modnn.o \
$(NNSRC)/klsp.o $(NNSRC)/ktall.o $(NNSRC)/ampall.o $(NNSRC)/mon.o \
$(NNSRC)/modbonn.o $(NNSRC)/potbonn.o $(NNSRC)/moff.o  \
$(NNSRC)/nnmain.o $(NNSRC)/nnutils.o
all: $(OBJ1)
	$(FC) -o mst $(OBJ) $(OBJ1) $(LFLAGS)
modmst.o:modmst.f90
	$(FC) $(FORM) $(FFLAGS) -c modmst.f90
mstutils.o:mstutils.f90 modmst.f90
	$(FC) $(FORM) $(FFLAGS) -c mstutils.f90
mstmain.o:mstmain.f90 modmst.f90
	$(FC) $(FORM) $(FFLAGS) -c mstmain.f90
readsmat.o:readsmat.f90 modmst.f90
	$(FC) $(FORM) $(FFLAGS) -c readsmat.f90
read3bwf.o:read3bwf.f90 modmst.f90
	$(FC) $(FORM) $(FFLAGS) -c read3bwf.f90
clean:
	-rm -f *.o core
	-rm -f $(NNSRC)/*.o $(NNSRC)/core
	-rm -f $(MSOSRC)/*.o $(MSOSRC)/core
superclean: clean
	-rm -f *.mod fort.* 
	cd $(NNSRC); $(MAKE) -f Makefile superclean
	cd $(MSOSRC); $(MAKE) -f Makefile superclean
# DO NOT DELETE







