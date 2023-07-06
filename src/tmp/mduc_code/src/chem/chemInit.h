#ifndef MDUC_CHEMINIT_H
#define MDUC_CHEMINIT_H

#include "chemStruct.h"
#include "ioStruct.h"

int setChemDefaults(chemStruct *chem);
int chemInitPeriodicTable(chemStruct *chem);
chemStruct *chemInit(char *mech, parseStruct *p);

#endif /* MDUC_CHEMINIT_H */
