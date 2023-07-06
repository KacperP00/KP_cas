# Directories
RUNDIR  := $(shell pwd | sed -e 's/\/src.*//')
HOMEDIR := $(shell pwd | sed -e 's/\/src.*//')
LIBDIR  := $(HOMEDIR)/lib
OBJDIR  := $(HOMEDIR)/obj
DEPDIR  := $(HOMEDIR)/dep
BINDIR  := $(HOMEDIR)/bin
INCDIR  := $(HOMEDIR)/include
#VPATH   := $(LIBDIR) $(BINDIR) $(OBJDIR)
VPATH   := $(LIBDIR) $(OBJDIR) $(DEPDIR)

# Default Options 

# take options from input file
include Makefile.in

# Files
C_DIRS  := 
TOOLS   := main plot stats
MODULES := mechanism solve

# look for include files in each of the modules
#INCFLAGS += $(patsubst %,-I%,$(MODULES)) 
VPATH += $(MODULES) $(TOOLS)
#vpath %.c $(HOMEDIR)/$(MODULES) $(HOMEDIR)/$(TOOLS)
#vpath %.h $(HOMEDIR)/$(MODULES) $(HOMEDIR)/$(TOOLS)
#$(warning VPATH is $(VPATH))

# binaries
BINS :=

# extra libraries if required
LIBS :=

SDIAL_INC = -I$(SDIAL)/include
SDIAL_LIB = -L$(SDIAL)/lib -lsundials_cvode -lsundials_ida -lsundials_kinsol \
													 -lsundials_nvecserial

# C source files will be added to this
CFILES :=

# function for substituting OFILES
OFILE = $(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(1)))

# include the source files from for each module/tool
include $(patsubst %,%/module.mk,$(MODULES))
include $(patsubst %,%/module.mk,$(TOOLS))

# function for converting case
toupper = $(shell echo $(1) | tr a-z A-Z)

# convert all the C files to object files
OFILES := $(patsubst %.c,%.o,$(filter %.c,$(CFILES)))
OFILED := $(addprefix $(OBJDIR)/%, $(notdir $(OFILES)))

# determine the dependency files
DFILES := $(addprefix $(DEPDIR)/,$(patsubst %.o,%.d,$(notdir $(OFILES))))

# add the necessary flags to compile
CFLAGS  += $(MDUC_DEFINE) $(SDIAL_INC) $(patsubst %,-I%,$(MODULES))
LDFLAGS += $(SDIAL_LIB)

# Targets
.DEFAULT_GOAL:
default: $(LIBS) Makefile Makefile.in

all: $(warning binaries $(BINS)) $(BINS)

.PHONY: clean
clean:
	rm -f $(OFILED);
	rm -f $(DFILES);
	rm -f $(LIBDIR)/$(strip $(LIBS));
	rm -f $(BINDIR)/$(strip $(BINS))

# general rule for creating object files from C source
$(OBJDIR)/%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

# calculate the C include dependencies
$(DEPDIR)/%.d: %.c
	$(CC) -MM -MG $< | sed -e "s|^\(.*\)\.o:|$@ $(OBJDIR)/$*.o:|" > $@

# general rule for creating libraries
.SECONDEXPANSION:
$(LIBDIR)/lib%.$(LIBSUFFIX): $$(call OFILE,$$($$(call toupper,$$(*F))_SRC))
	$(AR) $@ $(call OFILE,$($(call toupper,$(*F))_SRC));
	$(RL) $@

# include all the dependency files
ifneq ($(MAKECMDGOALS),clean)
-include $(patsubst %.o, $(DEPDIR)/%.d, $(notdir $(OFILES)))
endif

