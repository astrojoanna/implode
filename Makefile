F90=gfortran
FFLAGS=-O2 -ffree-form -fopenmp

%.o: %.f90

	$(F90) $(FFLAGS) -c $<

implode: main.o collisions.o maxwell.o parameters.o sort.o types.o eom.o constantsrkf.o constants.o
	$(F90) $(FFLAGS) $(LDFLAGS) $? -o $@

constants.o:    constants.f90
constantsrkf.o: constantsrkf.f90
eom.o:          eom.f90 constantsrkf.o constants.o
types.o:        types.f90
sort.o:         sort.f90 types.o
parameters.o:   parameters.f90 types.o constants.o
maxwell.o:      maxwell.f90 constants.o
collisions.o:   collisions.f90 parameters.o types.o constants.o
main.o:         main.f90 collisions.o maxwell.o parameters.o sort.o types.o eom.o constantsrkf.o constants.o

clean:
	rm -rf implode *.o *.mod energy*.dat output*.dat sizedistr*.dat fort.* 
