#!/bin/csh -f
#------------------------------------------------------
# link_cfd , version v2009.325 for x86_64-unknown-linux_intel

# do not modify, file is generated automatically
#------------------------------------------------------

echo " ===== Using custom link settings for flamelet library ====="
echo "   Base:     $FIRE_LINK_FILE"
echo "   Flamelet: $FIRE_MDUC"
echo "   Sundials: $FIRE_SDIAL"
echo " ==========================================================="
source $FIRE_LINK_FILE

set cfd_os_libs  = ($FIRE_MDUC/libmduc.a $FIRE_SDIAL/libsundials_nvecserial.a $FIRE_SDIAL/libsundials_cvode.a  $FIRE_SDIAL/libsundials_ida.a  $FIRE_SDIAL/libsundials_kinsol.a $cfd_os_libs )
