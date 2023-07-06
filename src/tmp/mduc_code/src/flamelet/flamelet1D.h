#ifndef MDUC_FLAMELET_1D_H
#define MDUC_FLAMELET_1D_H

#include "mducStructs.h"

/* ----- main solver functions ----- */
void solveFlamelet1D(double dt,double P,double phimin,double phimax,double *chi,
                     int inert,void *mem);
void flameletPrestep1D(driveStruct *driver,flameStruct *flame,chemStruct *chem,flowStruct *flow, ioStruct *io, timeStruct *time);

void flameletPreiteration1D(driveStruct *driver,flameStruct *flame, chemStruct *chem, flowStruct *flow,timeStruct *time, mpiStruct *mpi, double t);
/* ----- drivers for the various solvers ----- */
/* CVODE: differential solver */
int cvodeFlamelet1D(realtype t, N_Vector u, N_Vector udot, void *data);
/* IDA: differential-algebraic solver.  Required when solving enthalpy as
   energy equation */
int idaFlamelet1D(realtype t, N_Vector y, N_Vector yp, N_Vector r, void *data);

#ifdef ACCESS_MPI
#ifdef SUNDIALS_FROM_2_7
int cvodeFlamelet1D_local(long int Nlocal, realtype t, N_Vector u, N_Vector udot,
                          void *data);
#else
#ifdef SUNDIALS_FROM_3_0
int cvodeFlamelet1D_local(sunindextype Nlocal, realtype t, N_Vector u, N_Vector udot,
                          void *data);
#else
int cvodeFlamelet1D_local(int Nlocal, realtype t, N_Vector u, N_Vector udot,
                          void *data);
#endif /* SUNDIALS_FROM_3_0 */
#endif /* SUNDIALS_FROM_2_7 */
/* parallel communicators */
int bCommFlamelet1D(realtype *Y, void *data);
#endif /* ACCESS_MPI */

/* ----- provide the RHS of the flamelet equations ----- */
int flameletRHS_1D(double t, double *Y, double *RHS, void *args);
/* species equations */
void flameletY_RHS_1D(sdialStruct *sdial, driveStruct *driver, chemStruct *chem, flameStruct *flame,flowStruct *flow, double *T, double *Y, double *RHS);
void flameletConvVel1D(driveStruct *driver, chemStruct *chem, flameStruct *flame, flowStruct *flow,double *Y, double *T);
/* reaction rate quadrature */
void flameletQ_RHS_1D(sdialStruct *sdial, driveStruct *driver, chemStruct *chem, double *RHS);
/* reaction rate quadrature */
int cvodeFlameletQ1D(realtype t, N_Vector u, N_Vector udot, void *data);
/* temperature form of energy equation */
void flameletT_RHS_1D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                      flowStruct *flow, double *Y, double *RHS);
/* for solving the scalar dissipation rate as part of the solution */
void flameletChi_RHS_1D(driveStruct *driver,chemStruct *chem,flameStruct *flame,
                        flowStruct *flow, double *T, double *Y, double *RHS);

/* ----- functions to update properties during the time step ----- */
void updateEnergy1D(driveStruct *driver, chemStruct *chem, flameStruct *flame,
                    flowStruct *flow, double dt, double *Y, double *Yp);
void updateFlow1D(double t, driveStruct *driver, timeStruct *time,
                  flameStruct *flame, flowStruct *flow);

/* ----- monitoring functions ----- */
void flamelet1D_monitorHeader(ioStruct *io, sdialStruct *sdial);
void flamelet1D_monitor(int nst, double t, double dtPre, sdialStruct *sdial,
                        ioStruct *io, flowStruct *flow, flameStruct *flame);

int getTimeIndex(ioStruct *io, double tIn);
void flameletUpdateSoot1D(chemStruct *chem,flameStruct *flame,driveStruct *driver,ioStruct *io,
		          flowStruct *flow,int inert);
void flameletAllgatherv(chemStruct *chem,flameStruct *flame,driveStruct *driver,ioStruct *io,mpiStruct *mpi);
#endif /* MDUC_FLAMELET_1D_H */

