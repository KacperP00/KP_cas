# Define locations
BINDIR = ../bin
MODDIR = ../mod
OBJDIR = ../obj

# ======================================================================
# Let's start with the declarations
# ======================================================================

# The compiler
FC = ifort

.SUFFIXES:
.SUFFIXES: .o .f90 .mod

# flags for debugging or for maximum performance, comment as necessary
DBGFLAGS = -g#-fbounds-check

OPTFLAGS = -O3

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
F90FILES = precision.f90 math.f90 \
      pc_defs.f90 pc_database.f90 pc_func.f90 \
      rk_defs.f90 solver_defs.f90 spray_defs.f90 rk_func.f90 solver_func.f90 spray_func.f90 parser.f90 \
      cas.f90

# object files and module files 
OFILES = $(F90FILES:.f90=.o) $(CFILES:.c=.o)
MODFILES = $(F90FILES:.f90=.mod)

default: alldebug

all_debug: alldebug

all_opt: allopt

alldebug:
	@$(MAKE) $(CASDBG) $(DEBUG_OPTIONS)

allopt:
	@$(MAKE) $(CASOPT) $(OPT_OPTIONS)

$(CASOPT): $(OFILES)
	$(FC) $(FLAGS_ADD) -I$(MODDIR) $(OBJDIR)/*o -o $(BINDIR)/$(CASOPT) $(MODFLAGS)

$(CASDBG): $(OFILES)
	$(FC) $(FLAGS_ADD) -I$(MODDIR) $(OBJDIR)/*o -o $(BINDIR)/$(CASDBG) $(MODFLAGS)



$(OFILES): $(F90FILES)
	${FC} $(FLAGS_ADD) -c $*.f90 -o $(OBJDIR)/$@ $(MODFLAGS)

# Utility targets
.PHONY: clean veryclean

clean:
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(OBJDIR); rm -f $(OFILES)
	rm -rf *.o *.mod

veryclean:
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(BINDIR); rm -f $(CASOPT) $(CASDBG)
	rm -rf *.o *.mod
