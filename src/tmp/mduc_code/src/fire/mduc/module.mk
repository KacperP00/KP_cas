# try to get the sundials directory from a module if not specified 
ifeq ($(strip $(FIRE_SDIAL)),)
SDIAL_PATH = $(shell echo $$LD_LIBRARY_PATH | awk 'BEGIN{FS=":"};{for(i=1;i<=NF;i++){if(match($$i,"sundials")){print $$i}}}')
# the above will have the /lib directory
FIRE_SDIAL = $(dir $(SDIAL_PATH))
endif

# for compiling without BEST, get library locations from main Makefile (MDUC)
ifdef HOMEDIR
FIRE_MDUC := $(HOMEDIR)
endif
ifdef SDIAL
FIRE_SDIAL := $(SDIAL)
endif
$(warning Using flamelet installation -- $(FIRE_MDUC))
$(warning Using Sundials installation -- $(FIRE_SDIAL))

# define sources/objects for flamelet linking
MDUC_SRC += mduc_module.f mduc_usedef.f mduc_root.f mduc_sdm.f mduc_useini.f \
	mduc_usebnd.f mduc_usepsc.f mduc_useent.f mduc_useout.f mduc_usetim.f \
	mduc_cthchk.f mduc_cthtmp.f mduc_cyuseso.f \
	mduc_rif.f mduc_scalardiss.f mduc_soot.f 

# FIRE specific objects for coupling flamelet
MDUC_OBJS   := cthdis.o cthmfr.o
MDUC_OBJDIR := mduc/objects/v$(FIRE_YEAR)/serial

# necessary sundials libraries to be linked
SDIAL_LIBS := libsundials_nvecserial.a libsundials_cvode.a libsundials_ida.a \
	libsundials_kinsol.a libsundials_nvecparallel.a

# add to global user sources
FIRE_USRF    += $(addprefix mduc/,$(MDUC_SRC))
FIRE_USRLIBS += $(FIRE_MDUC)/lib/libmduc.a
FIRE_USRLIBS += $(addprefix $(FIRE_SDIAL)/lib/,$(SDIAL_LIBS))
FIRE_INC     += -I$(PWD)/mduc

# override global variables for parallel compilation
fireMPI.exe: MDUC_OBJDIR := mduc/objects/v$(FIRE_YEAR)/mpi

