#ifndef MDUC_SCALAR_DISSIPATION_H
#define MDUC_SCALAR_DISSIPATION_H

/* ----- functional forms of the scalar dissipation rate  ----- */
double (*calcChiForm)(double, double, double);
double calcLogForm(double Z, double Zref, double Zmax);
double calcErfcForm(double Z, double Zref, double Zmax);
double calcQuadForm(double Z, double Zref, double Zmax);

/* functions for computing 1D profiles for the analytical solutions of a
 counterflow diffustion flame and unsteady mixing layer */
void calcChi1DcounterflowDiffFlame(int nZ, double *Z, int j, double a,
                                   double *chi);
void calcChi1DUnsteadyLamMixLayer(int nZ, double *Z, double t, double *chi);
void calcChi1D(char *form, int n, double *x, double Zref, double Zmax, 
               double chiRef, double *chi1);

/* compute 2D scalar dissipation rate profile */
void calcChi2D(char *form, int *n, double **x, double *Zref, double *chiRef,
               double **chi1, double **chi12, double **chi2, int trans,
               int scaled);

/* ----- functions for calculating 1D profiles ----- */
#ifdef ACCESS_CFD
void calcChi1DPitsch(int nCFD1, int nCFD2, int CFDn3, double ***C_Z1MEAN,
                     double ***C_Z1VAR, double ***C_K, double ***C_D,
                     double ***C_MASS, int nZ, double *Z, double Z1ref,
                     double *chi);

void binChi1D(int nCFD1, int nCFD2, int nCFD3, double ***C_Z1MEAN,
              double ***C_Z1VAR, double ***C_K, double ***C_D, double ***C_MASS,
              int nZ, double *Z, double Z1ref, double *chi);
#endif // ACCESS_CFD

double calcBetaPDFStoich(double stoich, double mean, double var);

#endif /* MDUC_SCALAR_DISSIPATION_H */
