# Install script for directory: /home/md777659/codes/sundials/sundials-3.2.1/examples/cvode/fcmix_serial

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
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/fcvAdvDiff_bnd.f;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/fcvAdvDiff_bnd.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial" TYPE FILE FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/cvode/fcmix_serial/fcvAdvDiff_bnd.f"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/cvode/fcmix_serial/fcvAdvDiff_bnd.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/fcvDiurnal_kry_bp.f;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/fcvDiurnal_kry_bp.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial" TYPE FILE FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/cvode/fcmix_serial/fcvDiurnal_kry_bp.f"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/cvode/fcmix_serial/fcvDiurnal_kry_bp.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/fcvDiurnal_kry.f;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/fcvDiurnal_kry.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial" TYPE FILE FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/cvode/fcmix_serial/fcvDiurnal_kry.f"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/cvode/fcmix_serial/fcvDiurnal_kry.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/fcvRoberts_dns.f;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/fcvRoberts_dns.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial" TYPE FILE FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/cvode/fcmix_serial/fcvRoberts_dns.f"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/cvode/fcmix_serial/fcvRoberts_dns.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/fcvRoberts_dns_constraints.f;/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/fcvRoberts_dns_constraints.out")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial" TYPE FILE FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/cvode/fcmix_serial/fcvRoberts_dns_constraints.f"
    "/home/md777659/codes/sundials/sundials-3.2.1/examples/cvode/fcmix_serial/fcvRoberts_dns_constraints.out"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/README")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial" TYPE FILE FILES "/home/md777659/codes/sundials/sundials-3.2.1/examples/cvode/fcmix_serial/README")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/CMakeLists.txt")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial" TYPE FILE FILES "/home/md777659/codes/sundials/sundials-3.2.1/build/examples/cvode/fcmix_serial/CMakeLists.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial/Makefile")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/md777659/codes/sundials/sundials-3.2.1/install/examples/cvode/fcmix_serial" TYPE FILE RENAME "Makefile" FILES "/home/md777659/codes/sundials/sundials-3.2.1/build/examples/cvode/fcmix_serial/Makefile_ex")
endif()

