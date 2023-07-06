#ifndef MDUC_SOLVE_FLAMELET_2D_H
#define MDUC_SOLVE_FLAMELET_2D_H

#include "mducStructs.h"

/* functions to be called from SUNDIALS */
int cvodeRHS_2D(realtype t, N_Vector u, N_Vector udot, void *data);
int idaRHS_2D(realtype t, N_Vector y, N_Vector yp, N_Vector r, void *data);

/* compute the off-direction explicity terms */
void calcExplicit_2D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                     flowStruct *flow, int iRow, double *Cp, double *Y);

/* solve a row/column implicitly */
void calcY_RHS_imp(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                   flowStruct *flow, int iRow, double *T, double *Y,double *RHS);
void calcT_RHS_imp(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                   flowStruct *flow, int iRow, double *Y, double *RHS);

/* computes the parameters that are invariant over the time step */
void solvePrestep2D(double **chi1, double **chi12, double **chi2,
                    driveStruct *driver, flameStruct *flame, flowStruct *flow);
void updateFlow2D(double t, int iRow, driveStruct *driver, timeStruct *time,
                  flameStruct *flame, flowStruct *flow);
void updateEnergy2D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                    flowStruct *flow, double dt, int k, double *Y, double *Yp);

void updateYpre(flameStruct *flame);

/* enforce maximum time step for accuracy of explicit terms in ADI */
double calcMaxTimeStep(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                       flowStruct *flow, double CFL, double *cpmix, double *Y);

/* computes the 2D mixture specific heat */
void calcCpmix2D(chemStruct *chem, flameStruct *flame, double *Y, double *cpmix);

/* functions to switch pointers on solution structs for column solve of ADI */
void setFlameCol(flameStruct *flame, flameStruct *flameCol);
void setFlowCol(flowStruct *flow, flowStruct *flowCol);

#endif /* MDUC_SOLVE_FLAMELET_2D_H */
