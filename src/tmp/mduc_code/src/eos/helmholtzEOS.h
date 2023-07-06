#ifndef MDUC_HELMHOLTZ_EOS_H
#define MDUC_HELMHOLTZ_EOS_H

#include "eosFluid.h"

typedef struct hhGrad_ {
  double O;   /* value evaluated at delta and tau */
  double t;   /* 1st derivative w.r.t. tau */
  double tt;  /* 2nd derivative w.r.t. tau */
  double d;   /* 1st derivative w.r.t. delta */
  double dd;  /* 2nd derivative w.r.t. delta */
  double ddd; /* 3rd derivative w.r.t. delta */
  double dt;  /* cross derivative w.r.t. delta and tau */
} hhGrad;

typedef struct eosHH_ {
  int nT;   /* number of terms in the equation */
  int nPol; /* number of polynomial terms */
  int nExp; /* number of polynomial terms */
  double u[4],v[4]; /* Einstein coefficients */

  /* required coefficients and exponents for the functional form
     NOTE: functional form is different for polar and non polar, will be set
     based on flag defined above */
  int    dExp[24], powExp[24];
  double tExp[24], n[24];
} eosHH;


/* standard functions */
double HH_calcP_rhoT(double rho, double T, eosFluid *fluid);
double HH_calcT_rhoP(double rho, double P, eosFluid *fluid);
double HH_calcRho_PT(double P,   double T, eosFluid *fluid);

/* functions to calculate the */
void HH_calcAlphaIdeal(double tau,double delta,hhGrad *alpha0, eosFluid *fluid);
void HH_calcAlphaReal(double tau,double delta,hhGrad *alphaR,eosFluid *fluid);
void HH_calcAlphaD(double tau,double delta,hhGrad *alphaR,eosFluid *fluid);

void HH_setStateRhoT(double rho, double T, eosFluid *fluid, eosState *state);
void HH_setStateRhoP(double rho, double P, eosFluid *fluid, eosState *state);
void HH_setStateTP(double T, double P, eosFluid *fluid, eosState *state);


/* ----- functions for calculating saturation properties ----- */
void HH_calcSatT(double T,eosFluid *fluid,eosState *stateL,eosState *stateV);
void HH_calcSatP(double P,eosFluid *fluid,eosState *stateL,eosState *stateV);
void HH_calcSatRho(double rho,eosFluid *fluid,eosState *state);

/* characteristic functions/Jacobians */
void   HH_resSatT(int n,double *delta,double *res, void *args);
void   HH_jacSatT(int n,double *delta,double **jac,void *args);
double HH_brentSatP(double T, void *args);
double HH_brentSatRho(double T, void *args);
/* void HH_calcCharFunc(int n, double *delta, double *res, void *args); */
/* void HH_calcJacob(int n,double *res,double **jac,void *args); */
/* void HH_calcCharFunc(double *delta, double *res, int **ARGS); */
/* void HH_calcJacob(double *delta,double *res,double **jac,int **ARGS, */
/*                   void (*func)(double*, double*, int**)); */

#endif /* MDUC_HELMHOLTZ_EOS_H */
