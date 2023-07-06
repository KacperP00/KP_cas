# Install script for directory: /work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/src/sundials

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELEASE")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  MESSAGE("
Install shared components
")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sundials" TYPE FILE FILES
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_band.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_dense.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_direct.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_iterative.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_math.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_nvector.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_fnvector.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_pcg.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_sparse.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_spbcgs.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_spfgmr.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_spgmr.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_sptfqmr.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/sundials/sundials_types.h"
    )
endif()

