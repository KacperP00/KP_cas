# Install script for directory: /work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/src/cvodes

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
Install CVODES
")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/cvodes/libsundials_cvodes.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsundials_cvodes.so.2.9.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsundials_cvodes.so.2"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsundials_cvodes.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/cvodes/libsundials_cvodes.so.2.9.0"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/cvodes/libsundials_cvodes.so.2"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/cvodes/libsundials_cvodes.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsundials_cvodes.so.2.9.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsundials_cvodes.so.2"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsundials_cvodes.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/cvodes" TYPE FILE FILES
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/cvodes/cvodes_band.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/cvodes/cvodes_bandpre.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/cvodes/cvodes_bbdpre.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/cvodes/cvodes_dense.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/cvodes/cvodes_diag.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/cvodes/cvodes_direct.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/cvodes/cvodes.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/cvodes/cvodes_sparse.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/cvodes/cvodes_spbcgs.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/cvodes/cvodes_spgmr.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/cvodes/cvodes_spils.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/cvodes/cvodes_sptfqmr.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/cvodes" TYPE FILE FILES "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/src/cvodes/cvodes_impl.h")
endif()

