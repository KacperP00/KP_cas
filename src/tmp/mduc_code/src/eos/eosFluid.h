#ifndef MDUC_EOS_FLUID_H
#define MDUC_EOS_FLUID_H

#include "parser.h"

/* struct to store all state variables of the fluid */
typedef struct eosState_ {
  double P;   /* pressure (Pa) */
  double T;   /* temperature (K) */
  double rho; /* density (kg/m^3) */
  double a;   /* Helmholtz energy (J/kg) */
  double u;   /* specific internal energy (J/kg) */
  double h;   /* specific enthalpy (J/kg) */
  double s;   /* specific entropy (J/kg-K) */
  double g;   /* Gibbs energy (J/kg) */
  double cp;  /* isobaric specific heat (J/kg-K) */
  double cv;  /* isochoric specific heat (J/kg-K) */
  double Z;   /* compressibility factor */
  double w;   /* speed of sound (m/s) */
  double mu;  /* viscosity */
  double phi; /* ???? */
  
  double dP_drho, dP_dT;
} eosState;

typedef struct eosFluid_ {
  char name[64];   /* full name of fluid */
  char ident[16];  /* short standard identifier */
  char polar;      /* p: polar, n: non-polar */
  double Ru;       /* universal gas constant [J/(mol-K)] (TODO: move global) */
  
  /* ----- critical properties ----- */
  double Pc;    /* critical pressure    [Pa] */
  double Tc;    /* critical temperature [K] */
  double rhoc;  /* critical density     [kg/m^3] */
  double rhoMc; /* critical density     [mol/L] */
  double Zc;    /* critical compressibility [-] */
  double Tbp;   /* boiling temperature [K] */
  double Ttp;   /* ?? [K] */

  /* limiting ranges */
  double Pmin,Pmax,Tmin,Tmax,rhoMin,rhoMax;

  /* properties */
  double Mw;    /* molecular weight [g/mol] */
  double R;     /* gas constant [J/(kg-K)] */
  double omega; /* acentric factor [-] */

  /* equation of state point */
  char  eosName[32];
  void *eos;

  /* saturation properties */
  double *antoine;

  /* Coefficients for Ideal Gas Heat Capacity and Helmholtz Energy Equations */
  double a[2], c[3];

  /* properties for Helmholtz equation (TODO: here for now) */
  int nT;   /* number of terms in the equation */
  int nPol; /* number of polynomial terms */
  int nExp; /* number of polynomial terms */
  double u[4],v[4]; /* Einstein coefficients */

  /* required coefficients and exponents for the functional form
     NOTE: functional form is different for polar and non polar, will be set
     based on flag defined above */
  int    dExp[24], powExp[24];
  double tExp[24], n[24];

  /* numerical parameters */
  double satTOL; /* tolerance used when finding saturation point */
} eosFluid;

typedef struct eosSat_ {
  int       phase;
  double    P,T,rho;
  eosFluid *fluid;
  eosState *liq,*vap;
} eosSat;

typedef struct eosMixture_ {
  char name[64];
  int       nFl;
  char    **flNames;
  eosFluid *fluids;
  double   *X,*Y,*vol;   /* mole/mass fractions of fluid mixture */
} eosMixture;

eosFluid * fluidInit(char *name,char *eos,parseStruct *p);
eosFluid * paramGetFluid(char *name, parseStruct *p);
void       eosPrintFluid(eosFluid *f);

/* standard functions for computing parameters */
double (*eosCalcP_rhoT)(double, double, eosFluid *);
double (*eosCalcT_rhoP)(double, double, eosFluid *);
double (*eosCalcRho_PT)(double, double, eosFluid *);

double idealCalcP_rhoT(double rho, double T, eosFluid *fluid);
double idealCalcT_rhoP(double rho, double P, eosFluid *fluid);
double idealCalcRho_PT(double P,   double T, eosFluid *fluid);

double antoineP(double T, double *ant);
double antoineT(double P, double *ant);

double eosCalcAcentricFac(double pSat, eosFluid *f);

void eosCopyState(eosState *state1, eosState *state2);

#endif /* MDUC_EOS_FLUID_H */
