# Install script for directory: /work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/src/idas

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
Install IDAS
")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/idas/libsundials_idas.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsundials_idas.so.1.3.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsundials_idas.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsundials_idas.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/idas/libsundials_idas.so.1.3.0"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/idas/libsundials_idas.so.1"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/idas/libsundials_idas.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsundials_idas.so.1.3.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsundials_idas.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsundials_idas.so"
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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/idas" TYPE FILE FILES
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/idas/idas_band.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/idas/idas_bbdpre.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/idas/idas_dense.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/idas/idas_direct.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/idas/idas.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/idas/idas_sparse.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/idas/idas_spbcgs.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/idas/idas_spgmr.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/idas/idas_spils.h"
    "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/include/idas/idas_sptfqmr.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/idas" TYPE FILE FILES "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/src/idas/idas_impl.h")
endif()

