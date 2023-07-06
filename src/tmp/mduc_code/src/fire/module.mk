# user Subroutines (modify if needed)
FIRE_USRF := usedef.f useent.f useini.f usebnd.f useout.f usepsc.f usetim.f \
  use_cthchk.f use_cthtmp.f cyuseso.f cyusesu.f usecps.f

# user object files required by user subroutines (modify if needed)
FIRE_USRO := 

# user defined libraries to be linked
FIRE_USRLIBS :=

# get some information about the FIRE installation from the environment
# NOTE: this requires that the AVL FIRE directory be found, e.g. from module
FIRE_ROOT    := $(shell which cfdwm_solver | sed 's|/bin/cfdwm_solver||')
FIRE_VERSION := $(notdir $(FIRE_ROOT))
FIRE_YEAR    := $(basename $(FIRE_VERSION))

# exit if no FIRE present
ifeq ($(strip $(FIRE_ROOT)),)
$(error AVL FIRE must be present in environment to determine version)
else
$(warning FIRE installation -- $(FIRE_ROOT))
$(warning FIRE version      -- $(FIRE_VERSION))
endif

# specify anything version specific
ifeq ($(findstring 2011,$(FIRE_YEAR)),2011)
FIRE_BLD := $(FIRE_ROOT)/FIRE/v$(FIRE_VERSION)/bin/bin.x86_64*r5*/bld
else
FIRE_BLD := $(FIRE_ROOT)/FIRE/v$(FIRE_VERSION)/bin/bin.x86_64*/bld
endif
FIRE_LINK := $(FIRE_BLD)/link_cfd.csh

# default values for serial compilation
FIRE_INC    := -I$(PWD)/fire
FIRE_USRC   := -usr_comp=fire/comp_cfd.csh
FIRE_USRL   :=

# include the stuff from the interface
include fire/mduc/module.mk

# set some environment variables necessary for AVL compilation
export FIRE_INC_PATH:= $(FIRE_INC)

# some dummy functions for replacing commas with spaces later
comma := ,
empty :=
space := $(empty) $(empty)
FIRE_LIST   = $(subst $(space),$(comma),$(1))
ifneq ($(strip $(FIRE_USRLIBS)),)
export FIRE_LINK_FILE:= $(FIRE_LINK)
export FIRE_USER_LIBS:= $(FIRE_USRLIBS)
FIRE_USRL := -usr_link=fire/link_cfd.csh
endif

# set the target names
fire:     $(BINDIR)/fireSerial.exe
fire_mpi: $(BINDIR)/fireMPI.exe

# override global variables for parallel compilation
$(BINDIR)/fireMPI.exe: FIRE_MPI := -mpi
$(BINDIR)/fireMPI.exe: FIRE_USRC := -usr_comp=fire/comp_cfd.mpi.csh
$(BINDIR)/fireMPI.exe: export FIRE_LINK_FILE:= $(FIRE_BLD)/link_cfd.mpi.csh

# default rule for creating a fire executable
$(BINDIR)/%.exe: %.fpr
	@cfdwm_solver -name=fire -project=$< -no_start -exename=$@ \
		-uf=$(call FIRE_LIST,$(addprefix fire/,$(FIRE_USRF))) \
	  -ufo=$(call FIRE_LIST,$(addprefix fire/$(MDUC_OBJDIR)/,$(MDUC_OBJS))) \
		$(FIRE_USRC) $(FIRE_USRL) $(FIRE_MPI)
	@mv $(@F) $@
	@mv comp_out.log link_out.log $(OBJDIR)
	@rm -f start_mpi.csh


#	@mv $@ $(BINDIR)/$@

# use this empty rule to track whether files have changed
%.fpr: $(FIRE_USRF) $(FIRE_USRLIBS)
	@echo "Files changed: $?"
	@touch $@
