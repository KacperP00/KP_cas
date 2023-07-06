MDUC_MECHANISM_SRC := $(CHEMMECH).c

LIBS   += $(LIBDIR)/libmduc_mechanism.$(LIBSUFFIX)
CFILES += $(patsubst %.c,$(MECH_DIR)/%.c,$(MDUC_MECHANISM_SRC))
