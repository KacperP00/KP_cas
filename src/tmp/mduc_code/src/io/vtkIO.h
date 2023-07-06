#ifndef MDUC_VTK_IO_H
#define MDUC_VTK_IO_H

#include "mducStructs.h"

void VTKwriteFlamelet(char *fileName, ioStruct *io, driveStruct *driver,
                      timeStruct *time, flameStruct *flame, flowStruct *flow);

void readVTK(char *fileName);
void readFlameletVTK(char *fileName, ioStruct *io, driveStruct *driver,
                     timeStruct *time, flameStruct *flame, flowStruct *flow);

void VTKrealFieldData(FILE *pF,int tab,char *type,char *name,int n,char *format,
                      long long int offset,double *val);
void VTKmetaDataArray(FILE *pF,int tab,char *type,char *name,int n,char *format,
                      long long int offset);
void VTKwriteDataArray(FILE *pF,int m,int n,int offset,int index,double *data);

#endif /* MDUC_VTK_IO_H */
