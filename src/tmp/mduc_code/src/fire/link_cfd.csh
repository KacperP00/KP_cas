#!/bin/csh -f

# FIRE_LINK_FILE and FIRE_USER_LIBS are exported from makefile if any external
# libraries are present

echo " ===== Using custom link settings for external libraries ====="
echo "   Base:           $FIRE_LINK_FILE"
echo "   User libraries: $FIRE_USER_LIBS"
echo " ==========================================================="
source $FIRE_LINK_FILE

set cfd_os_libs = ( $FIRE_USER_LIBS $cfd_os_libs )
