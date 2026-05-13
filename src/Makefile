# Define locations
BINDIR = ../bin
MODDIR = ../mod
OBJDIR = ../obj

#withMDUC = YES

#Only for ndodecane
#withSOOT = YES

ifeq ($(withMDUC),YES)

# The compiler
FC = mpiifx

# LINK MDUC LIBRARY
ifeq ($(withSOOT),YES)
MDUC_DIR = .
MDUC_INC = -I$(MDUC_DIR)/include
MDUC_LIB = -L$(MDUC_DIR)/lib -lmduc -lmduc_chem -lmduc_mechanism
else
MDUC_DIR = .
MDUC_INC = -I$(MDUC_DIR)/include
MDUC_LIB = -L$(MDUC_DIR)/lib_ndodecaneCai -lmduc -lmduc_chem -lmduc_mechanism
endif

# LINK SUNDIALS
ifeq ($(withSOOT),YES)
SDIAL = . 
else
SDIAL = .
endif

# set up Sundials library
ifeq ($(strip $(SDIAL)),)
SDIAL = $(dir $(shell echo $$LD_LIBRARY_PATH | awk \
'BEGIN{FS=":"};{for(i=1;i<=NF;i++){if(match($$i,"sundials")){print $$i;exit;}}}'))
ifeq ($(strip $(SDIAL)),)
$(error no path to Sundials library detected)
else
$(warning using Sundials from $(SDIAL))
endif
endif
SDIAL_INC  := -I$(SDIAL)/include
SDIAL_LINK := -lsundials_cvode -lsundials_ida -lsundials_kinsol \
                   -lsundials_nvecserial
ifneq ($(findstring mpiifx,$(strip $(FC))),)
SDIAL_LINK += -lsundials_nvecparallel
endif
SDIAL_LIB := -L$(SDIAL)/lib $(SDIAL_LINK) -Wl,-rpath,$(SDIAL)/lib

INCFLAGS = $(MODDIR) $(MDUC_INC) $(SDIAL_INC)
LDFLAGS = -mkl -limf -lm $(MDUC_LIB) $(SDIAL_LIB)

ifeq ($(withSOOT),YES)
MDUCFLAGS = -DMDUC_MPI -DMDUC -DSOOT -DSUNDIALS_FROM_3_2 -DPREITER -DAPPROXEXP3 #-DINTEL
else
MDUCFLAGS = -DMDUC_MPI -DMDUC -DSUNDIALS_FROM_2_7 -DPREITER -DAPPROXEXP3 #-DINTEL
endif

else

FC = ifx

INCFLAGS = $(MODDIR)
LDFLAGS = #-mkl -limf -lm

MDUCFLAGS = 

endif
# ======================================================================
# Let's start with the declarations
# ======================================================================

.SUFFIXES:
.SUFFIXES: .o .f90 .mod

# flags for debugging or for maximum performance, comment as necessary
DBGFLAGS = -g -fpp $(MDUCFLAGS) #-heap-arrays #-fbounds-check

OPTFLAGS = -O3 -fpp $(MDUCFLAGS) #-heap-arrays

DEBUG_OPTIONS = "FLAGS_ADD = $(DBGFLAGS)"
OPT_OPTIONS   = "FLAGS_ADD = $(OPTFLAGS)"

# flags forall (e.g. look for system .mod files, required in gfortran)
ICFLAGS  = -I $(MODDIR)

MODFLAGS = -module $(MODDIR)

OBJFLAGS = -J $(OBJDIR)

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
CASDBG = cas
CASOPT = cas_opt


# objects
ifeq ($(withMDUC),YES)
F90FILES = precision.f90 math.f90 string.f90 \
      pc_defs.f90 pc_database.f90 pc_func.f90 fpt_defs.f90 fpt_func.f90 \
      rk_defs.f90 solver_defs.f90 spray_defs.f90 parser.f90 rk_func.f90 solver_func.f90 combust_rif_mduc_defs.f90 \
      combust_mduc.f90 spray_combust.f90 spray_func.f90 cas.f90

else
F90FILES = precision.f90 math.f90 string.f90 \
      pc_defs.f90 pc_database.f90 pc_func.f90 fpt_defs.f90 fpt_func.f90 \
      rk_defs.f90 solver_defs.f90 spray_defs.f90 parser.f90 rk_func.f90 solver_func.f90 spray_func.f90 \
      cas.f90
endif

# object files and module files 
OFILES = $(F90FILES:.f90=.o) $(CFILES:.c=.o)
MODFILES = $(F90FILES:.f90=.mod)

##$(info $$OFILES is [${OFILES}])

default: alldebug

all_debug: alldebug

all_opt: allopt

alldebug:
	@$(MAKE) $(CASDBG) $(DEBUG_OPTIONS)

allopt:
	@$(MAKE) $(CASOPT) $(OPT_OPTIONS)

$(CASOPT): $(OFILES)
	$(FC) $(FLAGS_ADD) -I$(INCFLAGS) $(OBJDIR)/*o $(LDFLAGS) -o $(BINDIR)/$(CASOPT) $(MODFLAGS)

$(CASDBG): $(OFILES)
	$(FC) $(FLAGS_ADD) -I$(INCFLAGS) $(OBJDIR)/*o $(LDFLAGS) -o $(BINDIR)/$(CASDBG) $(MODFLAGS)



$(OFILES): $(F90FILES)
	${FC} $(FLAGS_ADD) -c $*.f90 -o $(OBJDIR)/$@ $(MODFLAGS)

# Utility targets
.PHONY: clean veryclean

clean:
	cd $(MODDIR); rm -f *.mod
	cd $(OBJDIR); rm -f *.o
	rm -rf *.o *.mod

veryclean:
	cd $(MODDIR); rm -f *.mod
	cd $(OBJDIR); rm -f *.o
	cd $(BINDIR); rm -f $(CASOPT) $(CASDBG)
	rm -rf *.o *.mod
