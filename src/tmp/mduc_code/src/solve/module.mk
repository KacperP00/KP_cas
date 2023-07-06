MDUC_SOLVE_SRC := init.c destroy.c fileIO.c parser.c display.c utils.c \
                  rk4solver.c solve1D.c solve2D.c metric.c solveUtils.c scalarDiss.c \
	 	  accessMDUC.c fortranMDUC.c \
                  series.c errFunc.c hypergeometric.c polyInt.c

#LIBS   += $(LIBDIR)/libmduc_solve.$(LIBSUFFIX)
CFILES += $(patsubst %.c,solve/%.c,$(MDUC_SOLVE_SRC))
