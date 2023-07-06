#!/bin/csh -f
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   Buildscript for FIRE serial (all Fire versions, 64 Bit)
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# $Id: fire_buildscript_serial.csh,v 1.12 2011/03/17 14:16:28 clc2si Exp $
#
# CR/ARH1 Frank, Freikamp, Nikola Jovicic, 14.08.2008
#
# +++ FOR USAGE, PLEASE REFER TO INSTALL.txt +++
#
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   User Subroutines (modify if needed)
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

set fire_mpi = $1

setenv FIRE_USRF "usedef.f,useini.f,useent.f,usepsc.f,useout.f,use_cthtmp.f,use_cthchk.f,cyusesu.f,cyuseso.f"
setenv FIRE_USRO ""

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   specify external libraries for linking flamelet code
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

setenv FIRE_MDUC  "../../lib"
set    SDIAL = `grep SDIAL ../../src/Makefile.in | awk '{print $3}'`
setenv FIRE_SDIAL $SDIAL/lib

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   specify Fire version to be used
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

setenv CDTNG_FIRE_VERSION  '2009.3'
setenv FIRE_VERSION  `which cfdwm_solver | awk -F/ '{print "v"( $(NF-2) )}'`

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   DO NOT MODIFY LINES BELOW
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# define location of include header files needed for the compilation of the code
# used in comp_cfd.csh / com_cfd.mpi.csh
setenv INC_SEARCH_PATH $PWD

# use custom compile settings as defined in comp_cfd.csh
if ( $fire_mpi == 'mpi' ) then
    echo "===== compiling parallel executable ===== "
    setenv FIRE_USRC $PWD/comp_cfd.mpi.csh
    setenv TARGET 'fireMPI.exe'
else
    setenv FIRE_USRC $PWD/comp_cfd.csh
    setenv TARGET 'fireSerial.exe'
endif

# directory of source files
setenv FIRE_SRCDIR $PWD

# add variables from flamelet library if present
if ( -d 'mduc' ) then
    source "mduc/sources.csh"
    echo "===== Building BEST with MDUC ====="
    echo "  sources: $FIRE_USRF"
    echo "  objects: $FIRE_USRO"
    echo "==================================="
endif

# temporarily create a dummy project file
# cfdwm_solver expects it, however it is not really needed
touch dummy.fpr

setenv FIRE_PROJ 'dummy.fpr'

#------------------------------------------------------------------------
# to build a Fire executable the following command is used
set cmd = "cfdwm_solver -exename=$TARGET \
	     -name=fire \
	     -project=$FIRE_PROJ \
	     -uf=$FIRE_USRF \
	     -usr_comp=$FIRE_USRC \
	     -ufo=$FIRE_USRO \
	     -no_start"

#       -files=$FIRE_FILES \


# add linking options if compiling with flamelet library
if ( -d 'mduc' ) then
    set cmd = ($cmd ' -usr_link=$FIRE_USRL')
endif

# add mpi options if parallel
if ( $fire_mpi == 'mpi' ) then
    set cmd = ($cmd ' -mpi')
endif

# run the command
$cmd

#------------------------------------------------------------------------

# remove dummy project file
\rm -rf dummy.fpr

if ( $fire_mpi == 'mpi' ) then
    \mv cfdsolver.*.mpi ../../bin/fireMPI.exe
    \rm start_mpi.csh
else
    \mv cfdsolver.* ../../bin/fireSerial.exe
endif
