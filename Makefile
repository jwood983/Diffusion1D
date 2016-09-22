FC=gfortran
FFLAGS=-O3

all: diffusion1d

diffusion1d: solvers.o main.f90 ; $(FC) $(FFLAGS) solvers.o main.f90 -o diffusion1d 
    
solvers.o: solvers.f90 ; $(FC) $(FFLAGS) -c solvers.f90
	
clean: ; rm -rf *.o *.mod *.curve
