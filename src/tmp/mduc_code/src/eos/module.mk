EOS_SRC      := mducEOS.c
MDUC_EOS_SRC := eosFluid.c pengRobinson.c helmholtzEOS.c

LIBS   += $(LIBDIR)/libmduc_eos.$(LIBSUFFIX)
BINS   += mducEOS
CFILES += $(patsubst %.c,eos/%.c,$(MDUC_EOS_SRC))
CFILES += $(patsubst %.c,eos/%.c,$(EOS_SRC))

mducEOS: $(call OFILE,$(EOS_SRC)) $(LIBS) Makefile.in
	$(LD) $(call OFILE,$(EOS_SRC)) -o $(BINDIR)/mducEOS \
	-L$(LIBDIR) -lmduc $(LDFLAGS)
