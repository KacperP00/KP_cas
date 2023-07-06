# simple functions to make libraries from this level
SUB_LIBS = mduc spray export
SUB_LIBS_CLEAN = $(patsubst %,%.clean,$(SUB_LIBS))

.PHONY: $(SUB_LIBS) clean
all: $(SUB_LIBS)
all.clean: $(SUB_LIBS_CLEAN)

$(SUB_LIBS):
	$(MAKE) -C $@

$(SUB_LIBS_CLEAN):
	$(MAKE) -C $(@:.clean=) clean
