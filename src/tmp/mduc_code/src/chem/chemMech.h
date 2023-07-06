#ifndef MDUC_CHEM_MECH_H
#define MDUC_CHEM_MECH_H

#include "chemStruct.h"

/* function to compute the species enthalpy and specific heat for a given
   temperature */
void (*chemComputeThermoData)(double *, double *, double, chemStruct *);
void ComputeThermoDataComp(double *h, double *cp, double T, chemStruct *chem);
void ComputeThermoDataInterp(double *h, double *cp, double T, chemStruct *chem);
void calcSpeciesThermo(molStruct *sp,double *h,double *cp,double *S,double T);

/* function to compute the production rates */
void (*chemComputeProductionRates)(double *, double *, double *, double *, 
                                   double *, double, double, chemStruct *);
void ComputeProductionRatesComp(double *cdot, double *w, double *k, double *c, 
                                double *M, double temp, double pressure,
                                chemStruct *chem);
void ComputeProductionRatesInterp(double *cdot, double *w, double *k, double *c,
                                  double *M, double temp, double pressure,
                                  chemStruct *chem);

double calcReverseRxnRate(double T, double logT, double kf, rxnStruct *r,
                          chemStruct *chem);

double calcFallOffRateCoeff(double T, double logT, double RT, double P,
                            double k0, double kInf, double conc, rxnStruct *r);

void   chemFitReverseRxnRate(rxnStruct *r, chemStruct *chem);
void   fitArrhenius(double x, double *p, int np);
double calcGibbsFreeEnergy(double T, molStruct *sp);

/* species list manipulation */
int addSpecies(double newNu, char *newSp, int n, double **nu, char ***sp);
int removeSpecies(char *remSp, int n, double **nu, char ***sp);
int removeDuplicateSpecies(int n, double **nu, char ***sp);

/* Newton iteration to compute the temperature for a given total enthalpy
   and species mass fractions.  The guess should be set to the temperature
   at the neighbouring location if known */
double calcTfromH(int nSpecies, int nSolve, double H, double *Y, double guess,
                  chemStruct *chem);

void calcTotalEnthalpy(int n1, int n2, int nVars, int sT, double *Y, double *H,
                       chemStruct *chem);

#endif /* MDUC_CHEM_MECH_H */
