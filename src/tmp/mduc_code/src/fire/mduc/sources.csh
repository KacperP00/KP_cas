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

set FIRE_USRF = "$FIRE_USRF,./mduc/mduc_module.f,./mduc/mduc_usedef.f,./mduc/mduc_root.f,./mduc/mduc_sdm.f,./mduc/mduc_useini.f,./mduc/mduc_usepsc.f,./mduc/mduc_useent.f,./mduc/mduc_useout.f,./mduc/mduc_cthchk.f,./mduc/mduc_cthtmp.f,./mduc/mduc_rif.f,./mduc/mduc_scalardiss.f,./mduc/mduc_soot.f,./mduc/mduc_cyuseso.f"

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   specify object files to use based on version and parallelization
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
set FIRE_USRO_VERS = `echo $FIRE_VERSION | sed 's/\.[0-9]//'`
set FIRE_MDUC_USRO = "./mduc/objects/$FIRE_USRO_VERS"
if ( $fire_mpi == 'mpi' ) then
    set FIRE_MDUC_USRO = "$FIRE_MDUC_USRO/mpi"
else
    set FIRE_MDUC_USRO = "$FIRE_MDUC_USRO/serial"
endif

# ===== add objects to overall list
set FIRE_USRO = "$FIRE_USRO,$FIRE_MDUC_USRO/cthdis.o,$FIRE_MDUC_USRO/cthmfr.o"
#set FIRE_USRO = "$FIRE_USRO,$FIRE_MDUC_USRO/cthmfr.o"

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   specify libraries for the linker
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#set AVL_ROOT_DIR = '/software/AVL'
set AVL_ROOT_DIR = `which cfdwm_solver | sed 's|/bin/cfdwm_solver||'`
if ( $FIRE_VERSION == 'v2011.1') then
    set AVL_BLD_DIR  = `echo "$AVL_ROOT_DIR/FIRE/$FIRE_VERSION/bin/bin.x86_64*r5*/bld"`
else
    set AVL_BLD_DIR  = `echo "$AVL_ROOT_DIR/FIRE/$FIRE_VERSION/bin/bin.x86_64*/bld"`
endif

echo $AVL_BLD_DIR

setenv FIRE_USRL "mduc/link_cfd.csh" 

if ( $fire_mpi == 'mpi' ) then
    setenv FIRE_LINK_FILE "$AVL_BLD_DIR/link_cfd.mpi.csh"
else
    setenv FIRE_LINK_FILE "$AVL_BLD_DIR/link_cfd.csh"
endif
