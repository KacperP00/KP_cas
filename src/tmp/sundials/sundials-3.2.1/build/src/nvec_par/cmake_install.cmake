# Install script for directory: /home/md777659/codes/sundials/sundials-3.2.1/src/nvec_par

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
  MESSAGE("
Install NVECTOR_PARALLEL
")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/home/md777659/codes/sundials/sundials-3.2.1/build/src/nvec_par/libsundials_nvecparallel.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_nvecparallel.so.3.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_nvecparallel.so.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_nvecparallel.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/build/src/nvec_par/libsundials_nvecparallel.so.3.2.1"
    "/home/md777659/codes/sundials/sundials-3.2.1/build/src/nvec_par/libsundials_nvecparallel.so.3"
    "/home/md777659/codes/sundials/sundials-3.2.1/build/src/nvec_par/libsundials_nvecparallel.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_nvecparallel.so.3.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_nvecparallel.so.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_nvecparallel.so"
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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nvector" TYPE FILE FILES "/home/md777659/codes/sundials/sundials-3.2.1/include/nvector/nvector_parallel.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/home/md777659/codes/sundials/sundials-3.2.1/build/src/nvec_par/libsundials_fnvecparallel.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_fnvecparallel.so.3.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_fnvecparallel.so.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_fnvecparallel.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/build/src/nvec_par/libsundials_fnvecparallel.so.3.2.1"
    "/home/md777659/codes/sundials/sundials-3.2.1/build/src/nvec_par/libsundials_fnvecparallel.so.3"
    "/home/md777659/codes/sundials/sundials-3.2.1/build/src/nvec_par/libsundials_fnvecparallel.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_fnvecparallel.so.3.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_fnvecparallel.so.3"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_fnvecparallel.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/md777659/codes/sundials/sundials-3.2.1/build/src/nvec_par:"
           NEW_RPATH "")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

