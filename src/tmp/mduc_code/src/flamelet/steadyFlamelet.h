#ifndef MDUC_STEADY_FLAMELET_H
#define MDUC_STEADY_FLAMELET_H

#include "mducStructs.h"

void initSteadyFlamelet(flameStruct *flame, sdialStruct *kin);

int (*kinFlamelet)(N_Vector, N_Vector, void *);
int kinFlamelet1D(N_Vector u, N_Vector fval, void *f_data);
#ifdef SUNDIALS_FROM_2_7
int kinFlamelet1D_Loc(long int Nlocal, N_Vector u, N_Vector fval, void *data);
#else
#ifdef SUNDIALS_FROM_3_0
int kinFlamelet1D_Loc(sunindextype Nlocal, N_Vector u, N_Vector fval, void *data);
#else /* SUNDIALS LOWER THAN 2.7 */
int kinFlamelet1D_Loc(int Nlocal, N_Vector u, N_Vector fval, void *data);
#endif /* SUNDIALS_FROM_3_0 */
#endif /* SUNDIALS_FROM_2_7 */
int kinFlamelet2D(N_Vector u, N_Vector fval, void *f_data);
void newtFlamelet1D(int nF, double *a, double *f, void *args);

void solveSteadyFlamelet1D(void *mem);
void solveSteadyFlamelet2D(void *mem);

void steadyFlameletMonitor(void *mem);

#endif /* MDUC_STEADY_FLAMELET_H */
