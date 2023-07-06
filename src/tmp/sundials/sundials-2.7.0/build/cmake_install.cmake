# Install script for directory: /work/md777659/spray_model_development/codes/sundials/sundials-2.7.0

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sundials" TYPE FILE FILES "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/include/sundials/sundials_config.h")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/sundials/cmake_install.cmake")
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/nvec_ser/cmake_install.cmake")
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/nvec_par/cmake_install.cmake")
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/arkode/cmake_install.cmake")
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/arkode/fcmix/cmake_install.cmake")
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/cvode/cmake_install.cmake")
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/cvode/fcmix/cmake_install.cmake")
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/cvodes/cmake_install.cmake")
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/ida/cmake_install.cmake")
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/ida/fcmix/cmake_install.cmake")
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/idas/cmake_install.cmake")
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/kinsol/cmake_install.cmake")
  include("/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/src/kinsol/fcmix/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/work/md777659/spray_model_development/codes/sundials/sundials-2.7.0/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
