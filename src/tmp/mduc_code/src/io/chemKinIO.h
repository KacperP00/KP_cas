#ifndef MDUC_CHEMKIN_IO_H
#define MDUC_CHEMKIN_IO_H

#include "chemStruct.h"

/* get the number of species from the mechanism file */
int readChemKinSpecies(char *fileName, char ***species);
int parseReactionSpecies(char *rxn, char *delim, double **nu, char ***sp);

/* read in chemical mechanism from ChemKin file format */
int readChemKinMech(char *fileName, int *nAtoms, char ***atoms, int *nSpecs,
                    char ***species, rxnStruct **rxns, chemStruct *chem);

void writeChemKinMech(char *fileName, chemStruct *chem);
void writeChemKinRxn(FILE *pF, rxnStruct r);
int  writeChemKinRxnSpecies(FILE *pF, int nSp, double *nu, char **sp, int fo);

void parseMolecule(char *molID, molStruct *mol, int nElem, char **elems,
                   chemStruct *chem);
int getElementIndex(char *element, chemStruct *chem);

/* read in thermodynamic coefficients from ChemKin file format according to
   http://www2.galcit.caltech.edu/EDL/public/formats/chemkin.html  */
int readChemKinThermo(char *fileName, double *Trange, char ***species,
                      char ***info, char **phase, double ***spT, int **nAtoms,
                      int ***nAt, char ****atoms, double ***coeff);
int writeChemKinThermo(char *fileName, chemStruct *chem);

/* read/write transport properties from ChemKin formatted file, where columns
   are as follows:
     1 Species Name
     2 Geometry: 0 single atom 1 = linear molecule,  2 = non-linear molecule
     3 eps/k:    Lennard-Jones potential well depth / Boltzmann constant     [K]
     4 sigma:    Lennard-Jones collision diameter                     [Angstrom]
     5 mu:       Dipole moment                [Debye = 10^-18cm^(3/2)/erg^(1/2)]
     6 Zrot:     rotational relaxation collision number at 298K              [-]
 */
int readChemKinTrans(char *fileName, char ***species, double **eps_k,
                     double **sigma);
int writeChemKinTrans(char *fileName, chemStruct *chem);


void writeMechCompiledC(chemStruct *chem);
void writeMechRateCoeff_C(FILE *pF, int i, int f, double A, double b, double E);
void writeFallOffRateFunc(FILE *pF);

#endif  /* MDUC_CHEMKIN_IO_H */
