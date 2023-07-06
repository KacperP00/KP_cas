# Directories
HOMEDIR := $(shell pwd | sed -e 's/\/src.*//')
BASEDIR := $(basename $(HOMEDIR))
LIBDIR  := $(HOMEDIR)/lib
OBJDIR  := $(HOMEDIR)/obj
DEPDIR  := $(HOMEDIR)/dep
BINDIR  := $(HOMEDIR)/bin
INCDIR  := $(HOMEDIR)/include
VPATH   := $(LIBDIR) $(OBJDIR) $(DEPDIR)

# Default Options 
LD  = $(CC)
AR  = ar rcv
RL  = ranlib
DBGFLAGS = -g -O0
OPTFLAGS = -O3

# take options from input file
include Makefile.in

# process any info from Makefile.in
LIBSUFFIX = a
ifeq ($(strip $(CC)), gcc)
CFLAGS += -std=c99
endif
ifndef CHEMMECH
$(warning no compiled mechanism present)
MECHDIR  = mechanism
CHEMMECH = dummy_mech
else
MECHDIR = $(dir $(shell find mechanism -name $(CHEMMECH).c))
$(info compiling mechanism from directory $(MECHDIR))
VPATH += $(MECHDIR)
endif
MDUC_DEFINE = $(MDUC_ACCESS:%=-D%)

# seems to be no way to do this with target
ifneq ($(filter shared shared.% %.shared fluent,$(MAKECMDGOALS)),)
LIBSUFFIX = so
endif

# Files
C_DIRS  := 
TOOLS   := main plot rxnrates
MODULES := mechanism solve io chem stats flamelet eos engine
CPINCS  := solve/accessMDUC.h

ifneq ($(findstring fire,$(MAKECMDGOALS)),)
TOOLS += fire
endif

# for backup files
DATE		:= $(shell date +%y%m%d_%H%M)
TARFILE := $(DATE)_$(notdir $(BASEDIR))

# look for include files in each of the modules
#INCFLAGS += $(patsubst %,-I%,$(MODULES)) 
VPATH += $(MODULES) $(TOOLS)
#vpath %.c $(HOMEDIR)/$(MODULES) $(HOMEDIR)/$(TOOLS)
#vpath %.h $(HOMEDIR)/$(MODULES) $(HOMEDIR)/$(TOOLS)
#$(warning VPATH is $(VPATH))

# binaries
BINS :=

# extra libraries if required
LIBS := $(LIBDIR)/libmduc.$(LIBSUFFIX)

# set up Sundials library
#ifeq ($(strip $(SDIAL)),)
#SDIAL = $(dir $(shell echo $$LD_LIBRARY_PATH | awk \
#'BEGIN{FS=":"};{for(i=1;i<=NF;i++){if(match($$i,"sundials")){print $$i;exit;}}}'))
#ifeq ($(strip $(SDIAL)),)
#$(error no path to Sundials library detected)
#else
#$(warning using Sundials from $(SDIAL))
#endif
#endif
SDIAL_INC  := -I$(SDIAL)/include
SDIAL_LINK := -lsundials_cvodes -lsundials_ida -lsundials_kinsol \
	           -lsundials_nvecserial
#ifneq ($(findstring mpicc,$(strip $(CC))),)
SDIAL_LINK += -lsundials_nvecparallel
#endif
SDIAL_LIB := -L$(SDIAL)/lib $(SDIAL_LINK) -Wl,-rpath,$(SDIAL)/lib
#SDIAL_LIB := -L$(SDIAL)/lib64 $(SDIAL_LINK) -Wl,-rpath,$(SDIAL)/lib64

# C source files will be added to this
CFILES :=

# function for substituting OFILES
OFILE = $(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(1)))

# include the source files from for each module/tool
include $(patsubst %,%/module.mk,$(MODULES))

# Default target is to just build the libraries
.DEFAULT_GOAL:
default: $(LIBS) Makefile Makefile.in
	@cp $(CPINCS) $(INCDIR)

archive:
	@tar -cf $(BASEDIR)/../$(TARFILE).tar $(BASEDIR)
	@gzip $(BASEDIR)/../$(TARFILE).tar

# include the source for the binaries
# NOTE: currently, this must be done after default goal, but I'm unsure why
include $(patsubst %,%/module.mk,$(TOOLS))

# function for converting case
toupper = $(shell echo $(1) | tr a-z A-Z)
empty  :=
space  := $(empty) $(empty)

# combine all sources to main libarary
TMP_SRC  := $(patsubst %,MDUC_%_SRC,$(call toupper, $(MODULES)))
MDUC_SRC := $(foreach mod, $(TMP_SRC), $($(mod)))

# convert all the C files to object files
OFILES := $(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(CFILES))))

# determine the dependency files
DFILES := $(addprefix $(DEPDIR)/,$(patsubst %.o,%.d,$(notdir $(OFILES))))

# add the necessary flags to compile
export C_INCLUDE_PATH=.:$(subst $(space),:,$(MODULES))
CFLAGS  += $(MDUC_DEFINE) $(SDIAL_INC)
# $(patsubst %,-I%,$(MODULES))
LDFLAGS += $(SDIAL_LIB)

# ----- additional targets -----
.PHONY: clean all shared debug opt fluent
dist: $(BINS)
all:  $(BINS)

# for fluent compilation (assumes loaded in path)
fluent: FLUENT_DIR = $(patsubst %/fluent/bin/,%,$(dir $(shell which fluent)))
fluent: CC = $(shell find $(FLUENT_DIR)/commonfiles/MPI -name mpicc\
	 | grep -i platform)
fluent: LDFLAGS += -L$(shell find $(FLUENT_DIR)/commonfiles/MPI -name libmpio.so\
	 | grep linux_amd64/lib) -lmpi
fluent: INCFLAGS += -I$(shell find $(FLUENT_DIR)/commonfiles/MPI -name mpi.h\
	 | grep include/mpi.h)
fluent: MDUC_DEFINE += -DACCESS_MPI
fluent: shared
	$(info ==== compiling for Ansys Fluent $(FLUENT_DIR))

# parallel
mpi %.mpi: CC = mpicc
mpi %.mpi: LD = mpicc
mpi %.mpi: LDFLAGS += -lsundials_nvecparallel
mpi %.mpi: MDUC_DEFINE += -DACCESS_MPI
mpi: default
%.mpi: %
	$(info ==== building in parallel)
# shared flags
shared %.shared: CFLAGS += -fPIC -shared
shared: default
%.shared: %
	$(info ==== making shared libraries)
# debug flags
debug %.debug: CFLAGS += $(DBGFLAGS)
debug: default
%.debug: %
	$(info === making in debug mode)
# opt flags
opt %.opt: CFLAGS += $(OPTFLAGS)
opt: default
%.opt: %
	$(info === making in opt mode)

clean:
	rm -f $(OBJDIR)/*.o;
	rm -f $(DEPDIR)/*.d;
	rm -f $(INCDIR)/$(notdir $(CPINCS));
	rm -f $(strip $(LIBS));
	rm -f $(LIBDIR)/*.so;
	rm -f $(addprefix $(BINDIR)/,$(strip $(BINS)))

# ----- patterns -----
# general rule for creating object files from C source
$(OBJDIR)/%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

# calculate the C include dependencies
$(DEPDIR)/%.d: %.c
	gcc -MM -MG $< | sed -e "s|^\(.*\)\.o:|$@ $(OBJDIR)/$*.o:|" > $@

# general rule for creating libraries
.SECONDEXPANSION:
$(LIBDIR)/lib%.$(LIBSUFFIX): $$(call OFILE,$$($$(call toupper,$$(*F))_SRC))
	$(AR) $@ $(call OFILE,$($(call toupper,$(*F))_SRC));
	$(RL) $@

# include all the dependency files
ifeq ($(filter-out clean archive,$(MAKECMDGOALS)),$(MAKECMDGOALS))
-include $(DFILES)
endif

