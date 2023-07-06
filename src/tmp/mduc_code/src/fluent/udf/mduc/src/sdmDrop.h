#ifndef MDUC_SDM_DROPLET_H
#define MDUC_SDM_DROPLET_H

#include "mducUDF.h"

typedef struct SDMSTRUCT {
  char monitor[256];

  /* for conditional mean source term */
  int    nZ;             /* number of points on mixture fraction grid for SDM */
  double *Z;             /* mixture fraction grid */
  double *pdf;           /* memory for pdf storage */
  double W;              /* fuel source term*/
  double *Wv;            /* mean fuel source conditioned on mixture fraction*/
  double alpha, Zconst;  /* SDM model parameters for W|Z = alpha*Z^xi */
  double Zsat;           /* saturation mixture fraction */
  double Re_t;           /* turbulent cell Reynolds number */
  double C_x;            /* constant representing tau/tau_z */

  /* store some constanst to avoid computing in each cell */
  double one_third, six_over_pi;

  /* for Spalding number computation */
  char   NIST_File[200];
  double P, Yf_inf, T_inf, Q, T0;
  int    nNIST;
  double *T, *Psat, *cpG, *cpL, *hfg;
  double Mfuel, Mair;
  double xiMax;
} sdmStruct;

/* initialisation functions */
void SDMinit(void);
int  checkSDM(Domain *d);
void getInjection(char *name, int nNozzle, Injection ***I);

/* function definitions for root solve */
void   SDMcalcSrc(Domain *d, int nInj, int nNoz, Injection ***stream);
double SDMcalcVarsrc(double fVol, double fSrc, double mass, double tau,
                     double B, double mu_L, double Zmean, double Zvar,
                     double *fdt, double *vSrc, double *vEps);
double calcSpaldingFunc(double Ts, void *args);
double SDMmeanSrcFunc(double xi, void *args);

/* functions for computing droplet properties */
double calcHvapHeptane(double T);

/* files for reading data files generated from NIST Chemistry Webbook */
void getThermoNIST(char *fileName, double *hfg, double *cpgas, double *cpliq,
                   double *Psat, double T);
int readThermoNIST(char *fileName, double **hfg, double **cpgas, double **cpliq,
                   double **Psat, double **T);

void interpThermoNIST(double *hfg, double *cpG, double *cpL, double *Psat,
                      double T, int nNIST, double *NIST_hfg, double *NIST_cpG,
                      double *NIST_cpL, double *NIST_Psat, double *NIST_T);

void printThermoNIST(int nNIST, double *hfg, double *cpG, double *cpL,
                     double *Psat, double *T);



#endif /* MDUC_SDM_DROPLET_H */
