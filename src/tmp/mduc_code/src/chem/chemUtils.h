#ifndef MDUC_CHEM_UTILS_H
#define MDUC_CHEM_UTILS_H

#include "mducStructs.h"

void initChemFuel(chemStruct *chem);

double calcStoichMixRatio(chemStruct *chem);

rxnStruct * paramGetReaction(char *name, parseStruct *p);

int getSpeciesIndices(int nCheck, char **spCheck, int *iSp, int nM, char **spM);

/* functions to operate on chemistry structs */
void copyRxn(rxnStruct *rc, rxnStruct *r);
void destroySpecies(molStruct *sp);
void setRxnDefault(rxnStruct *r);
void destroyRxn(rxnStruct *r);
void destroyChem(chemStruct *chem);

#endif /* MDUC_CHEM_UTILS_H */
