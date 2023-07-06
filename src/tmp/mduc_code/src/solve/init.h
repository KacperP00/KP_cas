#ifndef MDUC_INIT_H
#define MDUC_INIT_H

#include "mducStructs.h"
#include "pdf_parameters.h"

void setDefaults(ioStruct *io, driveStruct *driver, timeStruct *time,
                 flameStruct *flame, flowStruct *flow, engineStruct *engine);

void initIO(ioStruct *io, driveStruct *driver, chemStruct *chem);

void initGrid(int n, double xmin, double xmax, int nP, double *pivot,
              double *perc, double **x, double **dx);

void initSolver(ioStruct *io, driveStruct *driver, chemStruct *chem,
                flameStruct *flame, flowStruct *flow);

void initTime(driveStruct *driver, timeStruct *time);

void initSUNDIALS(ioStruct *io, sdialStruct *sdial);

void initCVODE(int dim, int nMesh, int nVars, double tInit, double *YInit,
               driveStruct *driver, sdialStruct *cvode);

void initIDA(int dim, int nMesh, int nVars, double tInit, double reltol,
             double abstol, double *YInit, sdialStruct *ida);

void initFlow(driveStruct *driver, flameStruct *flame, flowStruct *flow);

void initX2dimension(driveStruct *driver, ioStruct *io, chemStruct *chem,
                     flameStruct *flame, flowStruct *flow, double *Y2);

double MDUCcalcSpaldingFunc(double Ts, void *args);

int cmpfunc(const void *x, const void *y);
#endif /* MDUC_INIT_H */
