# -*- mode: makefile -*-
# -----------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
#                David J. Gardner @ LLNL
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
# Makefile for ARKODE Fortran serial examples
#
# This file is generated from a template using variables
# set at configuration time. It can be used as a template for
# other user Makefiles.
# -----------------------------------------------------------------

SHELL = sh

prefix       = /home/md777659/codes/sundials/sundials-3.2.1/install
exec_prefix  = /home/md777659/codes/sundials/sundials-3.2.1/install
includedir   = /home/md777659/codes/sundials/sundials-3.2.1/install/include
libdir       = /home/md777659/codes/sundials/sundials-3.2.1/install/lib64

F77         = /opt/intel/Compiler/19.0/1.144/rwthlnk/bin/intel64/ifort
FFLAGS      = -O3
F77_LNKR    = /opt/intel/Compiler/19.0/1.144/rwthlnk/bin/intel64/ifort
F77_LDFLAGS = -O3
F77_LIBS    =  -lm /usr/lib64/librt.so

LINKFLAGS = -Wl,-rpath,/home/md777659/codes/sundials/sundials-3.2.1/install/lib64

# -----------------------------------------------------------------------------------------

LIBRARIES_LAPACK = -lsundials_fsunlinsollapackdense -lsundials_fsunlinsollapackband -lsundials_sunlinsollapackdense -lsundials_sunlinsollapackband  
LINKFLAGS_LAPACK = ${LINKFLAGS}::

INCLUDES_KLU  = 
LIBRARIES_KLU = -lsundials_fsunlinsolklu -lsundials_sunlinsolklu 
LINKFLAGS_KLU = ${LINKFLAGS}:

INCLUDES_SLUMT  = 
LIBRARIES_SLUMT = -lsundials_fsunlinsolsuperlumt -lsundials_sunlinsolsuperlumt  
LINKFLAGS_SLUMT = ${LINKFLAGS}::

TMP_INCS  = ${includedir} ${INCLUDES_SLUMT} ${INCLUDES_KLU}
INCLUDES  = $(addprefix -I, ${TMP_INCS})
LIBRARIES = -lsundials_farkode -lsundials_arkode -lsundials_fnvecserial -lsundials_nvecserial ${LIBS}

# -----------------------------------------------------------------------------------------

EXAMPLES =  fark_diurnal_kry_bp
OBJECTS  = ${EXAMPLES:=.o}

# LAPACK Examples
EXAMPLES_BL = 
OBJECTS_BL  = ${EXAMPLES_BL:=.o}

# KLU Examples
EXAMPLES_KLU = 
OBJECTS_KLU  = ${EXAMPLES_KLU:=.o}

# SuperLU_MT Examples
EXAMPLES_SLUMT = 
OBJECTS_SLUMT  = ${EXAMPLES_SLUMT:=.o}

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .f

.f.o :
	${F77} ${FFLAGS} ${INCLUDES} -c $<

# -----------------------------------------------------------------------------------------

all: examples examples_bl examples_klu examples_slumt

examples: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${F77_LNKR} -o $${i} $${i}.o ${F77_LDFLAGS} ${F77_LIBS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${F77_LNKR} -o $${i} $${i}.o ${F77_LDFLAGS} ${F77_LIBS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done

examples_bl: ${OBJECTS_BL}
	@for i in ${EXAMPLES_BL} ; do \
	  echo "${F77_LNKR} -o $${i} $${i}.o ${F77_LDFLAGS} ${F77_LIBS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_LAPACK} ${LINKFLAGS_LAPACK}" ; \
	  ${F77_LNKR} -o $${i} $${i}.o ${F77_LDFLAGS} ${F77_LIBS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_LAPACK} ${LINKFLAGS_LAPACK} ; \
	done

examples_klu: ${OBJECTS_KLU}
	@for i in ${EXAMPLES_KLU} ; do \
	  echo "${F77_LNKR} -o $${i} $${i}.o ${F77_LDFLAGS} ${F77_LIBS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_KLU} ${LINKFLAGS_KLU}" ; \
	  ${F77_LNKR} -o $${i} $${i}.o ${F77_LDFLAGS} ${F77_LIBS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_KLU} ${LINKFLAGS_KLU} ; \
	done

examples_slumt: ${OBJECTS_SLUMT}
	@for i in ${EXAMPLES_SLUMT} ; do \
	  echo "${F77_LNKR} -o $${i} $${i}.o ${F77_LDFLAGS} ${F77_LIBS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_SLUMT} ${LINKFLAGS_SLUMT}" ; \
	  ${F77_LNKR} -o $${i} $${i}.o ${F77_LDFLAGS} ${F77_LIBS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_SLUMT} ${LINKFLAGS_SLUMT} ; \
	done

clean:
	rm -f ${OBJECTS} ${OBJECTS_BL} ${OBJECTS_KLU} ${OBJECTS_SLUMT}
	rm -f ${EXAMPLES} ${EXAMPLES_BL} ${EXAMPLES_KLU} ${EXAMPLES_SLUMT}

# -----------------------------------------------------------------------------------------
