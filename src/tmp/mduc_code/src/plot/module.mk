PLOT_SRC := mducPlot.c gnuplotScript.c

# add to global
BINS   += mducPlot
CFILES += $(patsubst %.c,plot/%.c,$(PLOT_SRC))

mducPlot: $(call OFILE,$(PLOT_SRC)) $(LIBS) Makefile.in
	$(LD) $(call OFILE,$(PLOT_SRC)) -o $(BINDIR)/mducPlot \
	-L$(LIBDIR) -lmduc $(LDFLAGS)
