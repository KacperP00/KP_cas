#ifndef UTILS_H
#define UTILS_H

// functions for transforming between (Z1,Z2) and (x,y) space
void getZfromGrid(int trans, double x1, double x2, double *Z1, double *Z2);
void getGridfromZ(int trans, double Z1, double Z2, double *x1, double *x2);
void transformZ(double Z1, double Z2, double *x1, double *x2);
void retrieveZ(double x1, double x2, double *Z1, double *Z2);

// returns the location of the species (name) from the list of species (species)
int getSpeciesIndex(int nSpecies, char **species, char *name);

// functions to convert between mole/mass fractions
void convertMassToMoleFrac(int nSpecies, double *W, double *Y, double *X);
void convertMoleToMassFrac(int nSpecies, double *W, double *X, double *Y);

// returns the average mixture molecular weight
double calcMixtureW(int Ngi, double *Y, double *W);

// returns the average mixture specific heat
double calcMixtureCp(int numSpecies, double *cp, double *Y);

// returns the location of the maximum value in an array
int getIndexMax(int n, int offset, double *vec);

// returns the maximum value of an array (vec)
double getMinValue(int n1, int n2, int index, int offset, double *vec);
double getMaxValue(int n1, int n2, int index, int offset, double *vec);
double getArrayValueStr(char *val, int n1, int n2, int n3, double ***vec);

// find the adjacent grid point in x1 that is just less than x
int findAdjacentGridPoint(int n1, double *x1, double x);
int gridLocate(int n1, double *x1, double x);

/* string functions */
int  emptyString(char *a);
int  countChars(char *str, char a);
void stripString(char *str);
void prefixString(char *pre, char *str);

double clip(double x, double min, double max);
double percError(double val, double base);
int    factorial(int x);

/*
 return the value of y at a location x based on linear interpolation
 between mesh points x0 and xN with values of y0 and yN
 */
double linInterp(double x, double x0, double y0, double xN, double yN);
double biLinInterp(double x, double y, double x0, double x1, double y0,
                   double y1, double Y00, double Y01, double Y10, double Y11);
double lookupProfile(double xv, int n, double *x, double *y);
void interpSolution1D(int n1, double *x, int nI, double *xI, int nVars,
                      double *Y, double *YI);
void interpSolution2D(int n1, int n2, double *x1, double *x2, int n1I, int n2I,
                      double *x1I, double *x2I, int nVars, double *Y,double *YI);
void interpDNS_1D(int n, double *x, int norm, double xmax, double *var);

double calcMean(int x0, int xN, int index, int offset, double *pdf, double *Y);

double calcMeanX(int x0, int xN, int index, int offset, double *pdf, double *Y,
                 int nSpecies, double *W);

/* calculate grid growth rate given number of points, initial spacing, length*/
double calcGridGrowth(int n, double a0, double L);

/* ----- operations on arrays ----- */
void copyArray(int n, double *x, double *y);
void convertArrayToFloat(int n, double *x, float *y);
// extract variable from array multi unstruct. array with given index
void extractVarArray(int n1, int n2, int n3, int index, double *Y, double *var);
// set variable to multi unstruct. array at given index to value
void setVarArray(int n1, int n2, int n3, int index, double *Y, double *var);
// convert unstructured array to structured
void convertArray1D_MD(int n1, int n2, int n3, double *YU, double ***YS);
// convert structured array to unstructured
void convertArrayMD_1D(int n1, int n2, int n3, double ***YS, double *YU);
// transpose unstructured array from rows/columns to columns/rows
void transposeRowCol(int n1, int n2, int nVars, double *Y);

void interpSolution1D(int n1, double *x, int nI, double *xI, int nVars,
                      double *Y, double *YI);
void interpSolution2D(int n1, int n2, double *x1, double *x2, int n1I, int n2I,
                      double *x1I,double *x2I, int nVars, double *Y, double *YI);

// simple functions
int IMIN(int a, int b);
int IMAX(int a, int b);
#endif // UTILS_H
