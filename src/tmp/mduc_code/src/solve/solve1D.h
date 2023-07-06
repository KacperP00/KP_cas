#ifndef MDUC_SOLVE_FLAMELET_1D_H
#define MDUC_SOLVE_FLAMELET_1D_H

#include "mducStructs.h"

/* Function called by CVODE to compute the RHS of the system of equations */
int  cvodeRHS_0D(realtype t, N_Vector u, N_Vector udot, void *data);

/* function wrappers to compute residual for IDA solver */
int idaRHS_0D(realtype t, N_Vector y, N_Vector yp, N_Vector r, void *data);

/* General function to compute RHS of system of equations */
void calcY_RHS_0D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                  flowStruct *flow, int iRow, int iCol, double T, double *Y,
                  double *RHS);
void calcT_RHS_0D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                  flowStruct *flow, int iRow, int iCol, double *Y, double *RHS);

#endif /* MDUC_SOLVE_FLAMELET_1D_H */
