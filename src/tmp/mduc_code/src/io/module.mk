MDUC_IO_SRC := chemKinIO.c flameMasterIO.c ceaIO.c vtkIO.c iFileIO.c

LIBS   += $(LIBDIR)/libmduc_io.$(LIBSUFFIX)
CFILES += $(patsubst %.c,io/%.c,$(MDUC_IO_SRC))
