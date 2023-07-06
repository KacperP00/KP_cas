#ifndef MDUC_ENSIGHT_IO_H
#define MDUC_ENSIGHT_IO_H

#include "mducStructs.h"

//typedef struct ensightStruct {

//};
/* -------------------------------------------------------------------------- */

void initEnsight(ioStruct *io);

void writeEnsightGeom(ioStruct *io, driveStruct *driver, timeStruct *time,
                      flameStruct *flame, flowStruct *flow);

void writeEnsightCaseHead(ioStruct *io, driveStruct *driver,
                          timeStruct *time, flameStruct *flame, flowStruct *flow);

void writeEnsightData(char *fileName, char *varName, int m, int n,
                      int offset, int index, double *data);

void updateEnsightCaseTime(double newTime, ioStruct *io, timeStruct *time);

void writeBinaryArrayReal(FILE *pFile, int n, double *val);

#endif /* MDUC_ENSIGHT_IO_H */
