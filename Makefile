HOMEDIR = $(shell pwd | sed -e 's/\/src.*//')
LIBDIR  = $(HOMEDIR)/lib
MODDIR  = $(HOMEDIR)/mod
OBJDIR  = $(HOMEDIR)/obj
BINDIR  = $(HOMEDIR)/bin
VPATH   = $(LIBDIR) $(BINDIR) $(OBJDIR)

# Compiler and archiver
CC  = mpicc
CXX = mpicxx
F90 = mpif90
F77 = mpif90
LD = mpif90
AR  = xiar rcsv
RL  = echo

# Compiler flags
CFLAGS   =
F90FLAGS = -fpp
F77FLAGS =
INCFLAGS =
LDFLAGS  = -mkl -limf -lm
MODFLAGS = -module $(MODDIR)
DBGFLAGS = -g -CA -CB -CS -CV -traceback -debug all -WB
#DBGFLAGS = -g -CA -CB -CS -CV -traceback -debug all -WB -warn all -check all
OPTFLAGS = -O3 -ip -prec-div -xhost -axSSE4.2,SSE3,SSE2 -static-intel -heap-array

DEBUG_OPTIONS = "FLAGS_ADD = $(DBGFLAGS)" "LDFLAGS_ADD = $(DBGLDFLAGS)"
OPT_OPTIONS = "FLAGS_ADD = $(OPTFLAGS)" "LDFLAGS_ADD = $(OPTLDFLAGS)"

LIBDIR  = ./lib
MODDIR  = ./mod
BINDIR  = ./bin

.SUFFIXES: 
.SUFFIXES: .o .f90 .c

F90FILES = precision.f90 spray_defs.d90 parser.f90 cas.f90

CFILES = 
LIBFILE =

OFILES = $(F90FILES:.f90=.o) $(CFILES:.c=.o)
MODFILES = $(F90FILES:.f90=.mod)

default: $(LIBDIR)/$(LIBFILE)

$(LIBDIR)/$(LIBFILE): $(OFILES)
	cd $(OBJDIR); $(AR) $@ $(OFILES); $(RL) $@

# $(MODFLAGS) is replaced, since it is unavailable in some compilers
.f90.o:
	$(F90) $(F90FLAGS) $(FLAGS_ADD) $(INCFLAGS) -I$(MODDIR) -c $*.f90  -o $(OBJDIR)/$*.o $(MODFLAGS)
.c.o:
	$(CC) $(CFLAGS) $(INCFLAGS) -c $*.c -o $(OBJDIR)/$*.o

opt:
	@$(MAKE) $(OPT_OPTIONS)

debug:
	@$(MAKE) $(DEBUG_OPTIONS)

clean:
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(LIBDIR); rm -f $(LIBFILE)
	cd $(OBJDIR); rm -f $(OFILES)
