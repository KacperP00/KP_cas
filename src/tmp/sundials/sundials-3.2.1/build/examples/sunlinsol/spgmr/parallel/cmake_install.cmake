# Install script for directory: /home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/md777659/codes/sundials/sundials-3.2.1/install")
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
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/test_sunlinsol_spgmr_parallel.c;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/test_sunlinsol.h;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/test_sunlinsol.c;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/sundials_linearsolver.c;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/sundials_nvector.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel" TYPE FILE FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/test_sunlinsol_spgmr_parallel.c"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../test_sunlinsol.h"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../test_sunlinsol.c"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../../../src/sundials/sundials_linearsolver.c"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../../../src/sundials/sundials_nvector.c"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/test_sunlinsol_spgmr_parallel.c;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/test_sunlinsol.h;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/test_sunlinsol.c;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/sundials_linearsolver.c;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/sundials_nvector.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel" TYPE FILE FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/test_sunlinsol_spgmr_parallel.c"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../test_sunlinsol.h"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../test_sunlinsol.c"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../../../src/sundials/sundials_linearsolver.c"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../../../src/sundials/sundials_nvector.c"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/test_sunlinsol_spgmr_parallel.c;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/test_sunlinsol.h;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/test_sunlinsol.c;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/sundials_linearsolver.c;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/sundials_nvector.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel" TYPE FILE FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/test_sunlinsol_spgmr_parallel.c"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../test_sunlinsol.h"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../test_sunlinsol.c"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../../../src/sundials/sundials_linearsolver.c"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../../../src/sundials/sundials_nvector.c"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/test_sunlinsol_spgmr_parallel.c;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/test_sunlinsol.h;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/test_sunlinsol.c;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/sundials_linearsolver.c;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/sundials_nvector.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel" TYPE FILE FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/test_sunlinsol_spgmr_parallel.c"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../test_sunlinsol.h"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../test_sunlinsol.c"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../../../src/sundials/sundials_linearsolver.c"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/sunlinsol/spgmr/parallel/../../../../src/sundials/sundials_nvector.c"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/CMakeLists.txt")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel" TYPE FILE FILES "/home/md777659/codes/sundials/sundials-3.2.1/build/examples/sunlinsol/spgmr/parallel/CMakeLists.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel/Makefile")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/sunlinsol/spgmr/parallel" TYPE FILE RENAME "Makefile" FILES "/home/md777659/codes/sundials/sundials-3.2.1/build/examples/sunlinsol/spgmr/parallel/Makefile_ex")
endif()

