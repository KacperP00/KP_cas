#ifndef MDUC_FLAMELET_INIT_H
#define MDUC_FLAMELET_INIT_H

#include "mducStructs.h"

flameStruct *flameletInit(ioStruct *io, driveStruct *driver, chemStruct *chem,
                          flowStruct *flow);

/* ----- drivers for the various solvers ----- */
int sdialFlameletInit(double *Y0, driveStruct *driver, flameStruct *flame,
                      chemStruct *chem, sdialStruct *sdial);
/* CVODE: differential solver */
int cvodeInitFlamelet(double t0, double *Y0, sdialStruct *cvode);
/* IDA: differential-algebraic solver.  Required when solving enthalpy as
   energy equation */
int idaInitFlamelet(double t0, double *Y0, sdialStruct *cvode);
/* RK: just use this one for debugging! */
int rk4initFlamelet(double t0, double *Y0, void *mem);

#endif /* MDUC_FLAMELET_INIT_H */
