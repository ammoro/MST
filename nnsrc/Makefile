MSOSRC=../$(MSO)
MSTSRC=../$(MST)
.SUFFIXES: .f90 .o

OBJ=modnn.o klsp.o ktall.o ampall.o mon.o modbonn.o potbonn.o \
	moff.o utils.o nnmain.o nnutils.o
OBJ1=modnn.o klsp.o ktall.o ampall.o mon.o modbonn.o potbonn.o \
	moff.o utils.o nnutils.o nnmain.o main.o 
all: $(OBJ1)
	$(FC) -o nnamp $(OBJ1) $(LFLAGS)
objs: $(OBJ)
	cp *mod $(MSOSRC)
	cp *mod $(MSTSRC)
nnutils.o:modnn.f90 nnutils.f90
	$(FC) $(FORM) $(FFLAGS) -c nnutils.f90
utils.o:modnn.f90 utils.f90
	$(FC) $(FORM) $(FFLAGS) -c utils.f90
modnn.o:modnn.f90
	$(FC) $(FORM) $(FFLAGS) -c modnn.f90
	-cp amps.mod AMPS.mod $(MSOSRC)
modcomp.o: modcomp.f90
	$(FC) $(FORM) $(FFLAGS) -c modcomp.f90
modbonn.o:modbonn.f90
	$(FC) $(FORM) $(FFLAGS) -c modbonn.f90
nnmain.o:modnn.o nnmain.f90
	$(FC) $(FORM) $(FFLAGS) -c nnmain.f90
klsp.o:modnn.o klsp.f90
	$(FC) $(FORM) $(FFLAGS) -c klsp.f90
ktall.o:modnn.o ktall.f90 
	$(FC) $(FORM) $(FFLAGS) -c ktall.f90
ampall.o:modnn.o ampall.f90 
	$(FC) $(FORM) $(FFLAGS) -c ampall.f90
mon.o:modnn.o modbonn.o mon.f90
	$(FC) $(FORM) $(FFLAGS) -c mon.f90
ampon.o:modnn.o modbonn.o ampon.f90
	$(FC) $(FORM) $(FFLAGS) -c ampon.f90
potbonn.o:modnn.o modbonn.o potbonn.f90
	$(FC) $(FORM) $(FFLAGS) -c potbonn.f90
moff.o:modnn.o modbonn.o moff.f90
	$(FC) $(FORM) $(FFLAGS) -c moff.f90
main.o:modnn.f90 main.f90
	$(FC) $(FORM) $(FFLAGS) -c main.f90
clean:
	rm -f *.o core
superclean: clean
	rm -f *.mod fort.*
# DO NOT DELETE





