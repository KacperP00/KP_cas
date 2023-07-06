#ifndef MDUC_FLOW_STRUCT_H
#define MDUC_FLOW_STRUCT_H

typedef struct {
  char *func;     /* window function */
  int    nW;      /* number of weights (filter width) */
  double *w;      /* filter weights */
} filterStruct;

typedef struct {
  double min;
  double max;
  double del;
} flowBound;

struct FLOWSTRUCT
{
  // flags
  int constLewisNum;
  
  // Background/Operating/Average pressure in domain
  double P, Dt_P, P_m1, P_m2, P_pre;
  double phi1min, phi2min, phi1max, phi2max, D_phi1, D_phi2;
  double D_phi1_m1, D_phi2_m1;
  double phi1min_m1, phi2min_m1, phi1max_m1, phi2max_m1;
  double Dt_phi1min, Dt_phi2min, Dt_phi1max, Dt_phi2max;
/*   double phiMin[3][5], phiMax[3][5], D_phi[3]; */
/*   double **phiMin, **phiMax, **D_phi; */
  filterStruct fhat;
  flowBound **phi, **Dt_phi, *phiStar;
  flowBound **phiRaw, **phiHat, **Dt_phiRaw, **Dt_phiHat, **DtHat_phi;
/*   double Dt_phiMin[3], Dt_phiMax[3]; */
  // Extents without the dPdt source for enthalpy direction
  double phi1min_corr, phi2min_corr, phi1max_corr, phi2max_corr;

  /* ----- reference properties ----- */
  double chiRef[3]; // reference scalar dissipation rate
  double strain[3]; // strain rate for scalar dissipation
  double ZRef[3];   // reference mixture fraction

  /* ----- scalar dissipation rates ----- */
  double *chi1, *chi2, *chi12;          // current
  double *chi1_m1, *chi2_m1, *chi12_m1; // previous 
  double chiZ1max;

  // conditional mean of enthalpy source term on mixture fraction
  double *avgSdot;    // current
  double *avgSdot_m1; // previous

  // scalar dissipation rate equation variables
  double *strainCondZ; // strain rate conditioned on mixture fraction
  double *avgSdot_M;   // source term from evaporation for chi equation

  // conditional mean of enthalpy on mixture fraction
  double *avgH;

  /* ----- statistical properties ----- */
  double H_mean;     // Favre average enthalpy
  double phiMean[3]; // Favre average of axis variable
  double phiVar[3];  // variance of axis variable

  double *pdf; // probability density function for each direction

  /* ----- general system properties ----- */
  double mass;       /* mass in the system */
  double volume;     /* volume of the system */
  double meanRho;    /* mean density of the flow */
};
typedef struct FLOWSTRUCT flowStruct;

#endif /* MDUC_FLOW_STRUCT_H */
