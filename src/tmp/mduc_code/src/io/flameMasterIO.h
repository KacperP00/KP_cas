#ifndef MDUC_FLAMEMASTER_IO_H
#define MDUC_FLAMEMASTER_IO_H

#include "chemStruct.h"

/* TODO: write functions to parse ChemKin mechanism and thermo files  */
int readFlameMasterMech(char *fileName, int *nAtoms, char ***atoms, int *nSpecs,
                        char ***species, rxnStruct **rxns, chemStruct *chem);

#endif  /* MDUC_FLAMEMASTER_IO_H */
