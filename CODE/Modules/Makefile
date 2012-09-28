# Makefile

.SUFFIXES : .f90 

export HOSTNAME
# Fortran 90 compiler for MSSMat
#================================
F90=ifort
ifeq ($(HOSTNAME),wayra)
F90 = ifort
endif

ifeq ($(HOSTNAME),shui)
F90 = ifort
endif

ifeq ($(HOSTNAME),pyr)
F90 = ifort
endif

ifeq ($(HOSTNAME),tourab)
F90 = ifort
endif
#===============================

# Fortran 90 compiler for EM2C
ifeq ($(HOSTNAME),madness)
F90 = /usr/voltaire/mpi.pathcc.rsh/bin/mpif90 # Fortran 90 compiler for madness
#F90 = /home/quangt/mpif90
endif




MODULEDIR = .

F90MOD = -module

F90FLAGS = $(F90MOD) $(MODULEDIR)

OBJ := $(patsubst %.f90,%.o,$(wildcard *.f90) )

all : $(OBJ)

clean :
	rm -f *.o *.mod

%.o: %.f90
	$(F90) $(F90FLAGS) $(OPT) -c $<

Domain.o : Element.o Face.o Sources.o Receivers.o logical_input.o Vertex.o TimeP.o Subdomains.o Edge.o Comm.o PlaneW.o Neumann.o Surface.o 

Sources.o : pi.o 
