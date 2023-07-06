# Install script for directory: /home/md777659/codes/sundials/sundials-3.2.1/src

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

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/sundials/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/nvec_ser/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunmat_dense/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunmat_band/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunmat_sparse/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunlinsol_band/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunlinsol_dense/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunlinsol_spgmr/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunlinsol_spfgmr/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunlinsol_spbcgs/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunlinsol_sptfqmr/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunlinsol_pcg/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/nvec_par/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/arkode/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/arkode/fcmix/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/cvode/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/cvode/fcmix/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/cvodes/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/ida/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/ida/fcmix/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/idas/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/kinsol/cmake_install.cmake")
  include("/home/md777659/codes/sundials/sundials-3.2.1/build/src/kinsol/fcmix/cmake_install.cmake")

endif()

