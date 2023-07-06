# -*- mode: makefile -*-
# -----------------------------------------------------------------
# Programmer: Radu Serban @ LLNL
# -----------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# -----------------------------------------------------------------
# Makefile for ARKODE fortran parallel examples
#
# This file is generated from a template using various variables
# set at configuration time. It can be used as a template for
# other user Makefiles.
# -----------------------------------------------------------------

SHELL = sh

prefix       = /home/md777659/codes/sundials/sundials-3.2.1/install
exec_prefix  = /home/md777659/codes/sundials/sundials-3.2.1/install
includedir   = /home/md777659/codes/sundials/sundials-3.2.1/install/include
libdir       = /home/md777659/codes/sundials/sundials-3.2.1/install/lib64

FFLAGS       = -O3
F77_LDFLAGS  = -O3
F77_LIBS     =  -lm /usr/lib64/librt.so

LINKFLAGS = -Wl,-rpath,/home/md777659/codes/sundials/sundials-3.2.1/install/lib64

MPIF77       = /opt/intel/impi/2018.4.274/compilers_and_libraries/linux/mpi/bin64/mpiifort
MPI_INC_DIR  = .
MPI_FLAGS    = 
MPIF77_LNKR  = /opt/intel/impi/2018.4.274/compilers_and_libraries/linux/mpi/bin64/mpiifort
MPI_LIB_DIR  = .
MPI_LIBS     = 

TMP_INCS = ${includedir} ${MPI_INC_DIR}
INCLUDES = $(addprefix -I, ${TMP_INCS})
LIBRARIES = -lsundials_farkode -lsundials_arkode -lsundials_fnvecparallel -lsundials_nvecparallel ${LIBS}

EXAMPLES =  fark_diag_non_p fark_diag_kry_bbd_p
OBJECTS  = ${EXAMPLES:=.o}

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .f

.f.o :
	${MPIF77} ${FFLAGS} ${MPI_FLAGS} -I${MPI_INC_DIR} ${INCLUDES} -c $<

# -----------------------------------------------------------------------------------------

all: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${MPIF77_LNKR} -o $${i} $${i}.o ${MPI_FLAGS} ${F77_LDFLAGS} ${F77_LIBS} ${INCLUDES} -L${libdir} ${LIBRARIES} -L${MPI_LIB_DIR} ${MPI_LIBS} ${LINKFLAGS}" ; \
	  ${MPIF77_LNKR} -o $${i} $${i}.o ${MPI_FLAGS} ${F77_LDFLAGS} ${F77_LIBS} ${INCLUDES} -L${libdir} ${LIBRARIES} -L${MPI_LIB_DIR} ${MPI_LIBS} ${LINKFLAGS} ; \
	done

clean:
	rm -f ${OBJECTS}
	rm -f ${EXAMPLES}

# -----------------------------------------------------------------------------------------
