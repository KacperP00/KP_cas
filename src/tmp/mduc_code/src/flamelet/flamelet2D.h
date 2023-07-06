#ifndef MDUC_FLAMELET_2D_H
#define MDUC_FLAMELET_2D_H

#include "mducStructs.h"


/* ----- drivers for the various solvers ----- */
void solveFlamelet2D(double dt,double P,double phi1min,double phi1max,
                     double phi2min,double phi2max,double **chi1,double **chi12,
                     double **chi2,int inert,void *mem);

/* TODO: placeholders for implicit solves  */
int cvodeFlamelet2D(realtype t, N_Vector u, N_Vector udot, void *data);
int idaFlamelet2D(realtype t, N_Vector y, N_Vector yp, N_Vector r, void *data);

#ifdef ACCESS_MPI
#ifdef SUNDIALS_FROM_2_7
int cvodeFlamelet2D_local(long int Nlocal, realtype t, N_Vector u, N_Vector udot,
                          void *data);
#else
#ifdef SUNDIALS_FROM_3_0
int cvodeFlamelet2D_local(sunindextype Nlocal, realtype t, N_Vector u, N_Vector udot,
                          void *data);
#else /* SUNDIALS LOWER THAN 2.7 */
int cvodeFlamelet2D_local(int Nlocal, realtype t, N_Vector u, N_Vector udot,
                          void *data);
#endif /* SUNDIALS_FROM_3_0 */
#endif /* SUNDIALS_FROM_2_7 */
/* parallel communicators */
int bCommFlamelet2D(realtype *Y, void *data);
#endif /* ACCESS_MPI */

/* ----- provide the RHS of the flamelet equations ----- */
void flameletRHS_2D(double t, double *Y, double *RHS, void *args);
/* species equations */
void flameletY_RHS_2D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                      flowStruct *flow, double *T, double *Y, double *RHS);
/* temperature form of energy equation */
void flameletT_RHS_2D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                      flowStruct *flow, double *Y, double *RHS);

#endif /* MDUC_FLAMELET_2D_H */
