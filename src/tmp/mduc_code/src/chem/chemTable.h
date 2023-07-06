#ifndef MDUC_CHEMTABLE_H
#define MDUC_CHEMTABLE_H

#include "mducStructs.h"

typedef struct TABLE_STRUCT {
  int    nDim;
  int    *i,*n;
  double *dat;
} tabStruct;

typedef struct CHEMTABLE_STRUCT {
  /* ----- table properties ----- */
  char format[16];     /* fortran or C order */
  char model[64];      /* type of model for this table */

  /* ----- table co-ordinates ----- */
  int    nDim;            /* number of independent co-ordinates */
  int    *n;              /* length of each co-ordinate mapping array */
  double **x;             /* co-ordinate maps */
  double *xmin,*xmax;     /* co-ordinate extrema */
  char   **xNames;        /* co-ordinate names */

  /* ----- table variables ----- */
  int    nVar;            /* number of variables stored in table */
  char   **names;         /* variable names */
  int    nData;           /* length of each variable data array */
  double *data;           /* storage for all variables */
  double *datMin,*datMax; /* extrema of each data point */
  int    *mask;           /* mask for valid table region */
  int    iRho;            /* storage for density location */
  int    iProg;           /* storage for progress variable */
  int    iHR;             /* storage for heat release rate */
  int    nSpecies;        /* number of species */
  char   **species;       /* species */
  int    *iSp;            /* species indices */
  double *progMin;        /* min of progress variable */
  double *progMax;        /* max of progress variable*/

  /* ----- model specific variables ----- */
  double Pref;           /* reference pressure*/
  double Tref;           /* reference temperature*/
} chemTab;

/* ----- function for manipulating chemtables ----- */
void  *loadChemTable(char *chemTabFile, int nDim, char *type);
double lookupChemTable(char *name, double *x, double *vars, void *tab);
double lookupChemTableIndex(int *index, double *x, double *vars, void *tab);
double lookupChemTablePT(char *name, double *x, double *vars, void *tab,
                         void *tabEx);
int chemTableCheckLookup(char *name,int *nV,int *iV0,int *offset,chemTab *cT);
void   printChemTableInfo(void *tab);
double chemTableVarMin(char *name, void *tab);
double chemTableVarMax(char *name, void *tab);
double chemTableProgMax(double *x, void *tab);
double chemTableProgMin(double *x, void *tab);
double chemTableProgScaleP(double P, double *x, void *tab);

/* ----- retrieve table info ----- */
int    getChemTabNdim(void *tab);
int    getChemTabNvar(void *tab);
int    getChemTabNspecies(void *tab);
double getChemTabPref(void *tab);
double getChemTabTref(void *tab);
void   getChemTabVarNames(char **names, void *tab);
void   getChemTabSpecies(char **species, void *tab);

/* ----- interpolators ----- */
double chemTableInterp2D(char *name, double *vals, double *vars, chemTab *cT);
double chemTableInterp3D(char *name, double *vals, double *vars, chemTab *cT);
double chemTableInterp3DIndex(int *index,double *vals,double *vars,chemTab *cT);
double chemTableInterp4D(char *name, double *vals, double *vars, chemTab *cT);
double chemTableInterp4DIndex(int *index,double *vals,double *vars,chemTab *cT);

/* ----- index functions ----- */
/*  return correct index for 4D table in column major order */
int (*indx)(int, int*, int*);
int indxC(int nx, int *i, int *n);
int indxF(int nx, int *i, int *n);
int cTab_indx(int i, int j, int k, int p, int q, int *n);
int cti2(int i, int j, int k, int *n);
int cti3(int i, int j, int k, int p, int *n);
int cti4(int i, int j, int k, int l,int p, int *n);

#endif /* MDUC_CHEMTABLE_H */
