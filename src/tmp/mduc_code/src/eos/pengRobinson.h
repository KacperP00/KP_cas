#ifndef MDUC_PENG_ROBINSON_H
#define MDUC_PENG_ROBINSON_H

#include "eosFluid.h"

typedef struct eosPR_ {
  int      nFl;
  eosFluid *f1,*f2;

  double kappa;  /* constant characteristic for alpha */
} eosPR;

eosPR * PR_init(eosFluid *f);

double PR_calcP_rhoT(double rho, double T, eosFluid *fluid);
double PR_calcT_rhoP(double rho, double P, eosFluid *fluid);
double PR_calcRho_PT(double P,   double T, eosFluid *fluid);

void PR_calcParams(double T, double *a, double *b, eosFluid *fluid);

/* double PR_calcPfromRhoT(double rho, double T, eosPR *pr); */
/* double PR_calcTfromRhoP(double rho, double P, eosPR *pr); */
/* double PR_calcRhofromPT(double P,   double T, eosPR *pr); */

#endif /* MDUC_PENG_ROBINSON_H */
