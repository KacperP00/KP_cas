#ifndef MDUC_SOOT_INIT_H
#define MDUC_SOOT_INIT_H

#include "chemStruct.h"

void initSoot(chemStruct *chem);

/* function for initial seeding of soot moments */
void initSootMoments(int nMom, double *moments,double tol);
void reinitSootMoments(int nMom, double *moments, double tol);

#endif /* MDUC_SOOT_INIT_H */
