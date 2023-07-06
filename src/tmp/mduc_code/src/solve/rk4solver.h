#ifndef RK_SOLVER_H
#define RK_SOLVER_H
#include "mducStructs.h"

int rk4init(void *mem);
int rk4solve(void *mem);
int bCommFlameletRK(double *Y, void *data);

int rk4getStep(void *mem);
double rk4getStepSize(void *mem);
int rk4setStopTime(void *mem, double tEnd);
int rk4setMaxStep(void *mem, double dt_max);
#endif /* RK_SOLVER_H */
