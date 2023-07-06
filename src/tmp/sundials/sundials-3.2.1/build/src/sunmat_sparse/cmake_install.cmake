# Install script for directory: /home/md777659/codes/sundials/sundials-3.2.1/src/sunmat_sparse

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
Install SUNMATRIX_SPARSE
")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunmat_sparse/libsundials_sunmatrixsparse.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_sunmatrixsparse.so.1.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_sunmatrixsparse.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_sunmatrixsparse.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunmat_sparse/libsundials_sunmatrixsparse.so.1.2.1"
    "/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunmat_sparse/libsundials_sunmatrixsparse.so.1"
    "/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunmat_sparse/libsundials_sunmatrixsparse.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_sunmatrixsparse.so.1.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_sunmatrixsparse.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_sunmatrixsparse.so"
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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sunmatrix" TYPE FILE FILES "/home/md777659/codes/sundials/sundials-3.2.1/include/sunmatrix/sunmatrix_sparse.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunmat_sparse/libsundials_fsunmatrixsparse.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_fsunmatrixsparse.so.1.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_fsunmatrixsparse.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_fsunmatrixsparse.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES
    "/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunmat_sparse/libsundials_fsunmatrixsparse.so.1.2.1"
    "/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunmat_sparse/libsundials_fsunmatrixsparse.so.1"
    "/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunmat_sparse/libsundials_fsunmatrixsparse.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_fsunmatrixsparse.so.1.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_fsunmatrixsparse.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libsundials_fsunmatrixsparse.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/md777659/codes/sundials/sundials-3.2.1/build/src/sunmat_sparse:"
           NEW_RPATH "")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

