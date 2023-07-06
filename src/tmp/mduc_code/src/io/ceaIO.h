#ifndef MDUC_CHEM_CEA_IO_H
#define MDUC_CHEM_CEA_IO_H

#include "chemStruct.h"

int CEAreadThermo(char *fileName, chemStruct *chem);
int CEAwriteThermo(char *fileName, chemStruct *chem);

#endif /* MDUC_CHEM_CEA_IO_H */
