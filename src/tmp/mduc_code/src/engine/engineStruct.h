#ifndef MDUC_ENGINE_STRUCT_H
#define MDUC_ENGINE_STRUCT_H

#include "wiebe.h"

struct VALVESTRUCT {
  /* ----- general information ----- */
  int  type;
  char ident[25];

  /* ----- valve geometry ----- */
  double L_max;   /* maximum lift */
  double D_port;  /* port diameter (D_p) */
  double D_seat;  /* seat diameter (D) */
  double D_head;  /* valve head diameter (D_v) */
  double D_stem;  /* valve stem diameter (D_s) */
  double W;       /* seat width */
  double theta;   /* seat angle */

  /* ----- lift information ----- */
  char   liftFile[128];
  double L;       /* current lift */
  double *L_prof; /* valve profile */
  
  /* ----- timing ----- */
  double open, close;
  int    nP, nL;
  double *CA_P, *CA_L;

  /* ----- manifold conditions ----- */
  char   P_file[128];
  double *P_prof; /* pressure profile */
  double P;       /* current pressure */
  double Pr;      /* pressure ratio across the valve */
  double T;       /* inlet/outlet temperature */
  double H;       /* inlet/outlet enthalpy */
  double T_wall;  /* wall temperature */
  double rho;     /* density */
  double gamma;   /* ratio of specific heats */
  double nY;      /* number of species */
  int    *iY;     /* species indexing */
  double *Y;      /* species mass fractions */
  int    *maskY;  /* mask of species */
  double wmix,R;  /* mixture properties */

  /* ----- flow properties ----- */
  double A_f;     /* flow area */
  double mdot;    /* mass flow */
  char   mdotFile[128];
  double nMdot;
  double *CA_mdot;
  double *mdotProf;
};
typedef struct VALVESTRUCT valveStruct;
/* -------------------------------------------------------------------------- */

struct WOSCHNISTRUCT {
  double a,b,c,d;          /* exponents */
  double C1,C2,C3;         /* coefficients */
  double alpha;            /* scaling coefficient */
  double Pref, Tref, Vref; /* reference values */
  double Pscale;           /* pressure scaling for fit */
};
typedef struct WOSCHNISTRUCT woschniStruct;
/* -------------------------------------------------------------------------- */

struct INJECTIONSTRUCT {
  char   injFile[128];
  double CA0, CA1, dCA;    /* timing info */
  double T;                /* fuel temperature */
  int    nMdot;            /* number of points in profile */
  double *t, *mdot;        /* mass flow rate profile */
  int    nY;               /* number of fuel species */
  int    *iY;              /* indices of fuel species in mechanism */
  double *Y;               /* boundary condition fuel species */
};
typedef struct INJECTIONSTRUCT injStruct;
/* -------------------------------------------------------------------------- */

struct ENGINESTRUCT {
  char   type[32];  /* HCCI, spark, SACI */
  char   monFile[256]; /* monitor file */
  FILE   *monitor;

  /* ----- cylinder geometry ----- */
  double bore;        /* cylinder diameter (m) */
  double stroke;      /* piston stroke (m) */
  double rod;         /* connecting rod length (m) */
  double vol;         /* volume of cylinder (m^3) */
  double zp;          /* piston position (m) */
  double zTDC;        /* piston location at TDC (m) */
  double *vol_prof;   /* volume profile */
  double *dVdt;       /* volume rate of change */
  double V_TDC;       /* volume at top-dead-centre (m^3) */
  double V_BDC;       /* volume at bottom-dead-centre (m^3) */
  double V_SWP;       /* swept volume of the cylinder (m^3) */
  double Acyl;        /* cylinder cross-sectional area (m^2) */
  double mass, m0;    /* mass of gas in cylinder (kg) */
  double mfuel,mf0;   /* mass of fuel (kg) */
  double mair,mair0;  /* mass of air (kg) */
  double megr,megr0;  /* mass of EGR (kg) */
  double CR;          /* compression ratio */
  double phi,lambda;  /* stoiciometric */


  /* ----- speed and position ----- */
  double RPM;      /* engine speed */
  double CAD;      /* crank-angle-degrees */
  double CAD_0;    /* starting crank angle */
  double CAD_end;  /* stopping crank angle */
  double Dt_CAD;   /* time step in CAD */
  int    nCAD;     /* number of points in crank-angle profile */
  double *CAD_vec; /* vector of CAD for cycle */

  /* ----- global performance parameters ----- */
  double CA50;     /* combustion phasing */
  
  /* ----- valve information ----- */
  valveStruct intake, exhaust;
  char manifoldFile[128];
  double IVO, IVC;             /* intake valve timing */
  double EVO, EVC;             /* exhaust valve timing */

  /* ----- combustion models ----- */
  wiebeStruct *wiebe, wiebeSI, wiebeHCCI;
  
  /* ----- spark ignition  ----- */
  char   sparkModel[32];
  int    spark;     /* 0: off, 1: wiebe */
  int    ignited;   /* flag for checking if already ignited */
  int    fullBurn;  /* */
  double sparkAdv;  /* spark advance */
  double sparkCA;   /* spark timing in crank angle */
  double sparkD;    /* initial spark diameter */
  double sparkX;
  double sparkY;
  double Dflame;    /* diameter of premixed flame */
  double Aflame;    /* area of premixed flame */
  double m_prem;    /* mass of premixed zone */
  double nPremSp;   /* number of burnt product species to consider */
  char   **premSp;  /* names burnt species in premixed zone */
  int    *iPremSp;  /* indices of premixed species in mechanism */
  double *Yburnt;   /* burnt species for premixed zone */
  double Tprem;     /* temperature of premixed zone */
  double rhoPrem;   /* density of the premixed zone */
  double cpPrem;    /* specific heat of the premixed zone */
  double hB,hU;     /* specific enthalpy of burned/unburned zones */
  double MwPrem;    /* premixed zone molecular weight */
  double Rprem;     /* gas constant of premixed zone */
  double *Yprem;    /* species for premixed zone (all species) */
  double Vprem;     /* volume burned region */
  double dVbdt;     /* rate of burnt volume change */
  double Y_fuel;
  double qdotPrem;  /* premixed heat release rate [J/(m^3-s)] */
  double qPrem, qPrem_m1;
  double burnRateScale;
  char  burnRateFile[256];
  int    nBurnRate;
  double *burnRateCA;
  double *burnFrac, *burnRate;

  /* ----- injection ----- */
  int       nInj;   /* number of injections */
  char **injProf;   /* storage for reading injection profiles */
  double *inj_t;    /* storage for reading injection timings */
  double *injDur;   /* storage for reading injection timings */
  injStruct *inj;   /* injection structs  */

  /* ----- temperature and heat transfer info ----- */
  double T_wall;
  int    htModel;
  woschniStruct woschni,hohen;

  /* ----- I/O ----- */
  FILE *pMan;
};
typedef struct ENGINESTRUCT engineStruct;
/* -------------------------------------------------------------------------- */

#endif /* MDUC_ENGINE_STRUCT_H */
