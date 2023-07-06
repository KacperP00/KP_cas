#ifndef MDUC_UDF_ENSIGHTIO_H
#define MDUC_UDF_ENSIGHTIO_H

#include "exportUDF.h"

/* functions for writing to Ensight*/
void UDFexportEnsight(char *dir, char *prefix);
void UDFwriteEnsightCase(char *dir, char *prefix);
void UDFwriteEnsightGeom(char *dir, char *prefix);
void UDFwriteEnsightNodes(Thread *t, FILE *pFile, int id);
void UDFwriteEnsightUnstructGeom(Thread *t, FILE *pFile, int verts, char *cell);

void UDFwriteEnsightVelocity(char *dir, char *prefix);
void UDFwriteEnsightInjection(char *dir, char *prefix, char *ident,
                              int nNozzles, Injection **stream);




/* functions for writing scalars in parts */
void UDFwriteEnsightScalar(char *dir,char *prefix,char *name,char *abbrv,
                           int ind, int type);

void writeEnsightVector(Thread *t, FILE *pFile, int verts, char *cell);
void writeEnsightScalar(Thread *t, FILE *pFile, int verts, char *cell, int ind,
                        int type);

/* functions for writing geometry in single zone */
void UDFwriteEnsightGeom_SP(char *dir, char *prefix, char *ident);
void UDFwriteEnsightVelocity_SP(char *dir, char *prefix);
void UDFwriteEnsightScalar_SP(char *dir, char *prefix, char *name, char *abbrv,
                              int ind, int type);

#endif /* MDUC_UDF_ENSIGHTIO_H */
