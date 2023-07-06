#ifndef MDUC_UDF_VTKIO_H
#define MDUC_UDF_VTKIO_H

#include "exportUDF.h"

/* functions for writing to Ensight*/
void VTKexport();
void VTKwriteTimeInfo(char *dir, char *prefix);
void VTKwriteMultiBlock(int nZone, char **zoneNames);

int  VTKsetZones(int nZ, char **zNames, Thread ***zT);
int  VTKsetSurfs(int nZ, char **zNames, Thread ***zT);

void VTKwriteUnstructData(char *zone, int nZ, char **zNames);
void VTKwriteParticleData(char *inj_id, int nNoz, Injection **streams);
void VTKwriteWallSurf(char *zone, int nZ, char **zNames);
void VTKwriteCxSurf(char *zone, int nZ, char **zNames);

/* write data */
void VTKwriteDummyData(int n, FILE *pFile);

/* functions to write VTK XML lines */
void VTKinfoDataSet(FILE *p,int tab,double t,char *group,char *part,char *file);
void VTKstrDataSet(char *p,int tab,double t,char *group,char *part,char *file);
void VTKinfoDataArray(FILE *p,int tab,char *type,char *name,int n,char *format,
                      long long int offset);
void VTKstrDataArray(char *p,int tab,char *type,char *name,int n,char *format,
                     long long int offset);

int VTKcellType(int verts);
int XMLparseLine(char *xmlIn, char **type, char ***tags, char ***vals);

#endif /* MDUC_UDF_VTKIO_H */
