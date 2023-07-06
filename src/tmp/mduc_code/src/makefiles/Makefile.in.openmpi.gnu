# compile options
#CHEMMECH  = nHeptane.short
CHEMMECH = propane_gri211
MDUC_ACCESS = ACCESS_MPI

# Compiler and archiver
CC  = mpicc
CXX = mpicc
F90 = mpif90
F77 = mpif77
LD  = mpicc
AR  = ar rcv
RL  = ranlib

# External libraries
SDIAL = /opt/sundials/gnu

# Compiler flags
CFLAGS   = 
F90FLAGS =
F77FLAGS =  
LDFLAGS  = -lm
INCFLAGS = 
DBGFLAGS = -g -O0
OPTFLAGS = -O3
SOFLAGS  = -fPIC -shared -g

