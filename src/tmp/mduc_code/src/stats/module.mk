STATS_SRC      := mducStats.c statsTools.c
MDUC_STATS_SRC := statsUtils.c betaPDF.c smlPDF.c gaussPDF.c

# add to global
LIBS   += $(LIBDIR)/libmduc_stats.$(LIBSUFFIX)
BINS   += mducStats
CFILES += $(patsubst %.c,stats/%.c,$(MDUC_STATS_SRC))
CFILES += $(patsubst %.c,stats/%.c,$(STATS_SRC))

mducStats: $(call OFILE,$(STATS_SRC)) $(LIBS) Makefile.in
	$(LD) $(call OFILE,$(STATS_SRC)) -o $(BINDIR)/mducStats \
	-L$(LIBDIR) -lmduc $(LDFLAGS)
