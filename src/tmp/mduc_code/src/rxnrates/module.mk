RXN_SRC := mducRxnRates.c

# add to global
BINS   += mducRxnRates
CFILES += $(patsubst %.c,plot/%.c,$(RXN_SRC))

mducRxnRates: $(call OFILE,$(RXN_SRC)) $(LIBS) Makefile.in
	$(LD) $(call OFILE,$(RXN_SRC)) -o $(BINDIR)/mducRxnRates \
	-L$(LIBDIR) -lmduc $(LDFLAGS)
