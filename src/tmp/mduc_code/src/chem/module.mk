CHEM_SRC      := mducChem.c
MDUC_CHEM_SRC := chemInit.c chemMech.c transport.c chemTable.c chemUtils.c \
	               soot.c sootHMOM.c sootHMOM_f.c wiebe.c

LIBS   += $(LIBDIR)/libmduc_chem.$(LIBSUFFIX)
BINS   += mducChem
CFILES += $(patsubst %.c,chem/%.c,$(MDUC_CHEM_SRC))
CFILES += $(patsubst %.c,chem/%.c,$(CHEM_SRC))

mducChem: $(call OFILE,$(CHEM_SRC)) $(LIBS) Makefile.in
	$(LD) $(call OFILE,$(CHEM_SRC)) -o $(BINDIR)/mducChem \
	-L$(LIBDIR) -lmduc $(LDFLAGS)
