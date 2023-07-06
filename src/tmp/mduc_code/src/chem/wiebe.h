#ifndef MDUC_WIEBE_H
#define MDUC_WIEBE_H

/* declare necessary structs forward */
#include "chemStruct.h"

/* declare struct for storing Wiebe parameters */
typedef struct wiebeStruct_ {
  /* function parameters */
  double a;       /* leading coefficient of Wiebe function */
  double m;       /* exponent of Wiebe function */
  double mol0;    /* initial number of moles */
  double molF;    /* initial number of fuel moles */
  double V;       /* current volume */
  double dVdt;    /* volume rate of change */

  /* duration parameters */
  double alpha;   /* duration constant */
  double beta;    /* duration linear coefficient */
  double t;       /* storage for current t */
  double t0;      /* start of combustion */
  double dt;      /* combustion duration */
  double tConv;   /* need to store for converting CAD to seconds */
  
  /* for I/O */
  double src;
} wiebeStruct;

/* ============================== Functions ================================= */

int initWiebeSpecies(int *nA, char ***atoms, int *nSpMech, char ***species, 
                     chemStruct *chem);


/* function to return the "chemical source term" based on a Wiebe function 
   in the same format as the detailed chemistry.  That is, it will return
   cdot in kmol/(m^3-s) */
void ComputeProductionRatesWiebe(double *cdot, double *w, double *k, double *c,
                                 double *M, double temp, double pressure,
                                 chemStruct *chem);

/* Wiebe function and its derivative */
double chemCalcWiebe(double t, wiebeStruct w);
double chemCalcWiebe_dt(double t, wiebeStruct w);

#endif /* MDUC_WIEBE_H */
