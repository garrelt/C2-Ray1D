# F90 compiler
#F90 = mpif90
#F90 = ifort
F90 = gfortran
#F90 = f90
#F90 = pf90
#LDR     = $(F90)

# ifort options
#IFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -DIFORT
#IFORTFLAGS = -g -C -traceback -fpe0 -DIFORT
#IFORTFLAGS = -g -C -DIFORT
GFORTFLAGS = -O

#F90FLAGS1 = $(IFORTFLAGS) 
#F90FLAGS1 = -xW $(IFORTFLAGS) 
#F90FLAGS1 = -xB $(IFORTFLAGS) 
F90FLAGS1 = $(GFORTFLAGS) 

#F90FLAGS = $(F90FLAGS1) -openmp -DMPI
#F90FLAGS = $(F90FLAGS1) -DMPI
#F90FLAGS = $(F90FLAGS1)
#F90FLAGS = $(F90FLAGS1) -openmp
F90FLAGS = $(F90FLAGS1)

# Other compilers
#F90FLAGS = -O3 -u -fpe0 -vec_report -ipo #Lobster 6
#F90FLAGS = -O4 -omp -fpe1# -fpe4 Compaq
#F90FLAGS  = -O3 -openmp

#F90FLAGS = -g -O0 -openmp #seg fault
#F90FLAGS = -g -O0 
#F90FLAGS = -fpe0 -g -u 
#F90FLAGS =  -Wv,-Pprocs,4 -O3 -cpu:opteron

#F90FLAGS = -O3

OPTIONS = $(F90FLAGS)

LDFLAGS = $(OPTIONS)
LIBS = 

UTILS=romberg.o string.o

CONSTANTS = mathconstants.o cgsconstants.o  cgsphotoconstants.o  cgsastroconstants.o c2ray_parameters.o abundances.o atomic.o

C2Ray_1D: precision.o $(CONSTANTS) $(UTILS) file_admin.o sizes.o no_mpi.o clocks.o file_admin.o grid.o tped.o  cosmology.o material.o cooling.o radiation.o thermal.o time.o doric.o photonstatistics.o cosmological_evolution.o evolve.o output.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ precision.o $(UTILS) file_admin.o sizes.o no_mpi.o clocks.o grid.o tped.o  cosmology.o material.o cooling.o radiation.o thermal.o time.o doric.o photonstatistics.o cosmological_evolution.o  evolve.o output.o C2Ray.o

clean : 
	rm -f *.o *.mod *.l *.il

.f.o:
	$(F90) -c $(OPTIONS) $<

.f90.o:
	$(F90) -c $(OPTIONS) $<

.F90.o:
	$(F90) -c $(OPTIONS) $<

f.mod:
	$(F90) -c $(OPTIONS) $<

.SUFFIXES: .f90 .F90 .mod .o


