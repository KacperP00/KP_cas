MDUC_ENGINE_SRC := engine.c engineUtils.c

CFILES += $(patsubst %.c,engine/%.c,$(MDUC_ENGINE_SRC))
