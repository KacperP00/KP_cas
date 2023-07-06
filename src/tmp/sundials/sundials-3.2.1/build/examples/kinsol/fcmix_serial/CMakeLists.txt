# -*- mode: cmake -*-
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
# CMakeLists.txt for KINSOL serial examples
#
# This file is generated from a template using variables
# set at configuration time. It can be used as a template for
# other user CMakeLists configuration files.
# -----------------------------------------------------------------

cmake_minimum_required(VERSION 2.8)

# Specify project name
PROJECT(FKINSOL_serial_examples Fortran)

MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)

# Examples using SUNDIALS linear solvers
SET(examples  fkinDiagon_kry)

# Examples using LAPACK linear solvers
SET(examples_bl )

# Examples using KLU linear solvers
SET(examples_klu )

# Examples using SuperLU_MT linear solvers
SET(examples_slumt )

# Specify path to SUNDIALS header files
SET(SUNDIALS_INC_DIR
  /home/md777659/codes/sundials/sundials-3.2.1/install/include
  CACHE STRING
  "Location of SUNDIALS header files")

# Add path to SUNDIALS header files
INCLUDE_DIRECTORIES(${SUNDIALS_INC_DIR})

# Set search path for SUNDIALS libraries 
SET(SUNDIALS_LIB_DIR /home/md777659/codes/sundials/sundials-3.2.1/install/lib64)

# Find the SUNDIALS solver's library
FIND_LIBRARY(SUNDIALS_SOLVER_LIB
  sundials_kinsol ${SUNDIALS_LIB_DIR}
  DOC "KINSOL library")

FIND_LIBRARY(SUNDIALS_SOLVER_FLIB
  sundials_fkinsol ${SUNDIALS_LIB_DIR}
  DOC "KINSOL Fortran-C library")

# Find the NVECTOR library
FIND_LIBRARY(SUNDIALS_NVEC_LIB
  sundials_nvecserial ${SUNDIALS_LIB_DIR}
  DOC "NVECTOR library")

FIND_LIBRARY(SUNDIALS_NVEC_FLIB
  sundials_fnvecserial ${SUNDIALS_LIB_DIR}
  DOC "NVECTOR Fortran-C library")

# Set additional libraries
SET(SUNDIALS_EXTRA_LIB  -lm /usr/lib64/librt.so CACHE STRING "Additional libraries")

# Set other libraries libraries
SET(BLAS_LIBRARIES  CACHE STRING "Blas libraries")
SET(LAPACK_LIBRARIES  CACHE STRING "Lapack libraries")

INCLUDE_DIRECTORIES()
SET(SUPERLUMT_LIBRARIES  CACHE STRING "SuperLUMT libraries")

INCLUDE_DIRECTORIES()
SET(KLU_LIBRARIES  CACHE STRING "KLU libraries")

# List of Sundials libraries shared across all examples
SET(SUNDIALS_LIBS
  ${SUNDIALS_SOLVER_FLIB}
  ${SUNDIALS_SOLVER_LIB}
  ${SUNDIALS_NVEC_FLIB}
  ${SUNDIALS_NVEC_LIB}
  ${SUNDIALS_EXTRA_LIB})


# Build examples with SUNDIALS linear solvers one by one
FOREACH(example ${examples})

  # example source files
  ADD_EXECUTABLE(${example} ${example}.f)

  # libraries to link against
  TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})

ENDFOREACH(example ${examples})


# Build each example with LAPACK linear solvers one by one
IF(LAPACK_LIBRARIES)
  FOREACH(example ${examples_bl})

    # Find the Sundials linear solver libraries
    FIND_LIBRARY(SUNDIALS_SUNLINSOLLAPACKBAND_LIB
      sundials_sunlinsollapackband ${SUNDIALS_LIB_DIR}
      DOC "Sundials banded LAPACK linear solver library")

    FIND_LIBRARY(SUNDIALS_SUNLINSOLLAPACKBAND_FLIB
      sundials_fsunlinsollapackband ${SUNDIALS_LIB_DIR}
      DOC "Sundials banded LAPACK linear solver Fortran-C library")

    FIND_LIBRARY(SUNDIALS_SUNLINSOLLAPACKDENSE_LIB
      sundials_sunlinsollapackdense ${SUNDIALS_LIB_DIR}
      DOC "Sundials dense LAPACK linear solver library")

    FIND_LIBRARY(SUNDIALS_SUNLINSOLLAPACKDENSE_FLIB
      sundials_fsunlinsollapackdense ${SUNDIALS_LIB_DIR}
      DOC "Sundials dense LAPACK linear solver Fortran-C library")

    # example source files
    ADD_EXECUTABLE(${example} ${example}.f)
  
    # libraries to link against
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_SUNLINSOLLAPACKBAND_LIB})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_SUNLINSOLLAPACKDENSE_LIB})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_SUNLINSOLLAPACKBAND_FLIB})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_SUNLINSOLLAPACKDENSE_FLIB})
    TARGET_LINK_LIBRARIES(${example} ${LAPACK_LIBRARIES})

    IF(BLAS_LIBRARIES)
      TARGET_LINK_LIBRARIES(${example} ${BLAS_LIBRARIES})
    ENDIF()

  ENDFOREACH(example ${examples_bl})
ENDIF()


# Build each example with KLU linear solvers one by one
IF(KLU_LIBRARIES)
  FOREACH(example ${examples_klu})

    # Find the Sundials linear solver libraries
    FIND_LIBRARY(SUNDIALS_SUNLINSOLKLU_LIB
      sundials_sunlinsolklu ${SUNDIALS_LIB_DIR}
      DOC "Sundials KLU linear solver library")

    FIND_LIBRARY(SUNDIALS_SUNLINSOLKLU_FLIB
      sundials_fsunlinsolklu ${SUNDIALS_LIB_DIR}
      DOC "Sundials KLU linear solver Fortran-C library")

    # example source files
    ADD_EXECUTABLE(${example} ${example}.f)
  
    # libraries to link against
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_SUNLINSOLKLU_LIB})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_SUNLINSOLKLU_FLIB})
    TARGET_LINK_LIBRARIES(${example} ${KLU_LIBRARIES})

  ENDFOREACH(example ${examples_klu})
ENDIF()


# Build each example with SuperLU_MT linear solvers one by one
IF(SUPERLUMT_LIBRARIES)
  FOREACH(example ${examples_slumt})

    # Find the Sundials linear solver libraries
    FIND_LIBRARY(SUNDIALS_SUNLINSOLSLUMT_LIB
      sundials_sunlinsolsuperlumt ${SUNDIALS_LIB_DIR}
      DOC "Sundials KLU linear solver library")

    FIND_LIBRARY(SUNDIALS_SUNLINSOLSLUMT_FLIB
      sundials_fsunlinsolsuperlumt ${SUNDIALS_LIB_DIR}
      DOC "Sundials KLU linear solver Fortran-C library")

    # example source files
    ADD_EXECUTABLE(${example} ${example}.f)
  
    # libraries to link against
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_SUNLINSOLSLUMT_LIB})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_SUNLINSOLSLUMT_FLIB})
    TARGET_LINK_LIBRARIES(${example} ${SUPERLUMT_LIBRARIES})

    IF(BLAS_LIBRARIES)
      TARGET_LINK_LIBRARIES(${example} ${BLAS_LIBRARIES})
    ENDIF()

  ENDFOREACH(example ${examples_slumt})
ENDIF()
