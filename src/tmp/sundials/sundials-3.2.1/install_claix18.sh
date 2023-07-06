#!/bin/bash
sundials_path=$PWD
build_path="$sundials_path/build"
install_path="$sundials_path/install"

mkdir -p $build_path
mkdir -p $install_path

cd $build_path

cmake .. -DCMAKE_BUILD_TYPE="RELEASE"\
         -DCMAKE_C_FLAGS_RELEASE="-O2 -fp-model precise -DNDEBUG"\
	 -DMPI_MPICC="$MPI_BINDIR/$MPICC"\
         -DEXAMPLES_ENABLE="OFF"\
         -DFCMIX_ENABLE="ON"\
         -DMPI_ENABLE="ON"\
	 -DMPI_RUN_COMMAND="$MPI_BINDIR/$MPIEXEC"\
         -DCMAKE_INSTALL_PREFIX="$install_path"

make install
  
 
