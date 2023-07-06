#ifndef MDUC_FLAMELET_UTILS_H
#define MDUC_FLAMELET_UTILS_H

#include "mducStructs.h"

flameStruct *paramGetFlamelet(char *name, driveStruct *driver, chemStruct *chem,
                              flowStruct *flow, parseStruct *p);

/* char paramGetBC(char **name,int *m,int *nSp,char ***spNames,double **spY, */
/*                 double *T,double *H, char **inp); */
void writeFlameletBC(FILE *pF,bndStruct *bnd);
bndStruct * paramGetBC(char *name,int n,char **inp);

#endif /* MDUC_FLAMELET_UTILS_H */
