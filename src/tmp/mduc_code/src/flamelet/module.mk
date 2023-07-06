MDUC_FLAMELET_SRC := flameletInit.c flameletUtils.c flamelet1D.c flamelet2D.c \
	steadyFlamelet.c

LIBS   += $(LIBDIR)/libmduc_flamelet.$(LIBSUFFIX)
CFILES += $(patsubst %.c,flamelet/%.c,$(MDUC_FLAMELET_SRC))
