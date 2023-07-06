MAIN_SRC := mducMain.c mducTools.c postProcess.c

# add to global
BINS   += mduc
CFILES += $(patsubst %.c,main/%.c,$(MAIN_SRC))

mduc: $(call OFILE,$(MAIN_SRC)) $(LIBS) Makefile.in
	$(LD) $(call OFILE,$(MAIN_SRC)) -o $(BINDIR)/mduc \
	-L$(LIBDIR) -lmduc $(LDFLAGS)
