#ifndef MDUC_FLUENT_UDF_H
#define MDUC_FLUENT_UDF_H

/* general C header files */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* general FLUENT header files */
#include "udf.h"
#include "sg.h"
#include "version.h"

/* FLUENT header files for parallel run */
#ifdef MDUC_MPI
#include "prf.h"
#include "mpi.h"
#endif /* MDUC_MPI */

/* udf library header files */
#include "udfScalars.h"
#include "udfMemory.h"
#include "udfTiming.h"
#include "sdmDrop.h"

typedef struct {
  char *func;     /* window function */
  int    nW;      /* number of weights (filter width) */
  double *w;      /* filter weights */
} filterStruct;  /* 28 Byte */

/* struct for storing UDF library options */
typedef struct MDUC_UDF_STRUCT {
  /* ansys stuff */
  int    ansysRel[3];
#ifdef PARALLEL
  MPI_Comm comm;
#endif PARALLEL
  
  /* memory for flamelet solver */
  void   *mduc;
  int    transform;

  /* flags for activating models */
  int    rif;          /* use Representative Interactive Flamelet */
  int    premixed;     /* use premixed chem-table */
  int    rifPremixed;  /* use combined rif-premixed */
  int    restart;      /* initialize or restart? */
  int    sdm_model;    /* use single droplet evaporation model? */
  int    wiebe;        /* use Wiebe combustion model */
  int    monitor;      /* use internal monitor functions */
  int    rifInert;     /* use inert solution for RIF update */
  int    rifEnthalpy;  /* enthalpy is a RIF co-ordinate */
  int    rifDim;       /* number of dimensions for RIF */

  /* flow variables to be stored */
  int    injection;   /* current number injections that have occured */
  int    sdm_active;
  double tPilot;      /* time to initialise pilot injection */
  double tMain;       /* time to initialise main injection */
  double PAvg;        /* volume average pressure */
  double PAvg_M1;     /* volume averaged pressure at t(n-1) */
  double dPdt;        /* pressure rate */

  /* store the indices of the scalars to be solved */
  int    nUDS, nUDM;
  char   **udsNames, **udmNames;
  int    iPhi[3], iPhiV[3];
  int    iZ1, iZ1V, iZ2, iZ2V, iZ0V;
  int    iH, iHV, iHMOD, iHRrate, iHRprem;
  int    iEGR, iEGRV, iProg, iProgV, iProgSrc;
  int    iRIFSrc, iProgTotal,iProgMax,iProgMin,iDiffSc,initHfromTable;
  double Tmin, Tmax, delT;
  double Hmin, Hmax, delH;
  double HminHat, HmaxHat;
  double dHminHat_dt, dHmaxHat_dt, dHmin_dtHat, dHmax_dtHat;
  double phiMin[3], phiMax[3], deltaPhi[3], deltaPhi2[3], phiZero[3];
  double phiMin_m1[3], phiMax_m1[3], dPhiMin_dt[3], dPhiMax_dt[3];
  double phiMinV[3], phiMaxV[3];
  double **phiMinHat, **phiMaxHat;
  double dPhiMinHat_dt[3], dPhiMaxHat_dt[3];
  double dPhiMin_dt_hat[3], dPhiMax_dt_hat[3];
  double phiBuf[3];

  /* UDM memory */
  int    iFuelV[3], iCHIR[3];
  int    iZV_SRC[3],  iZV_EPS[3], iFDT[3];
  /* for debugging SDM */
  int    iXI[3], iAlphaB[3], iSPALD, iReTURB;
  int    iZV_PRD[3], iZV_D2[3], iZV_SRC_R[3], iZV_SRC_L[3], iZV_EPS_R[3];
  int    iSL, iST, iCrif;
  int    iRegime; /* 0: no reaction, 1: rif, -1: premixed */

/*   filterStruct phiMinHat[3], phiMaxhat[3]; */
/*   filterStruct dPhiMin_dt_hat[3], dPhiMax_dt_hat[3]; */

  /* variables for the filtered */
  filterStruct fhat;

  /* premixed/partially premixed variables */
  void   *chemTabPrem; /* premixed chemistry table */
  void   *chemTabPres; /* premixed pressure table */
  int    chemTabPdep;
  int    nDimPrem;     /* number of dimensions in premixed table */
  int    nTabSp;       /* number of species in table */
  int    *iTabSp;      /* indexing of table species to mechanism */
  double *chemTabX;    /* array for storing values for chemtable */
  double Y_O2;         /* threshold value for regmie indicator */
  int    nProgSp;      /* number of species defining progress variable */
  int    *iProgSp;     /* index of progress variable species */
  char   **progSp;     /* progress variable species */
  double progMax;      /* maximum value of progress variable */
  double progSrcScale; 
  double a_thresh;     /* threshold for RIF/Premixed switching */
  double gradC_prem;   /* maximum gradient along flame front from flamelet */
  int    ignited;      /* flag for ignited */
  double ignTime;      /* time to ignite flame */
  double ignScale;     /* scaling of max progress variable in flame kernel */
  char   *ignGeom;     /* geometry of ignition region */
  double *ignCoord;    /* ignition co-ordinates */
  double ignRadius;    /* flame kernel radius for sphereical initialization */
  double *ignSize;     /* size of flame kernel rectangle: L x W x H */
  char   *sTmodel;     /* turbulent flame speed model */
  int    iAC, iBC, iPscale;
  int    avgMix;

  /* spray variables */
  int    nInj;         /* number of injections */
  int    nNoz;         /* number of injector nozzles */
  double T_drop_ref;   /* reference temperature of the droplet */
  double H_vap_ref;    /* reference enthalpy of vapour for source term */
  char   **inject_id;  /* names of the injections */
  Injection ***stream; /* Fluent pointer to particle streams for injections */

  /* info for mapping species tracked in Fluent */
  int   nY;         /* number of species in Fluent */
  int  *indexY;     /* index of FLUENT species in mechanism */
  int  *indexTab;   /* index of Table species in mechanism -> -1
                       if species is not in table, tableindex if stored in table
                       (Ymean[i]=Ytab[indexTab]) */
  char **species;   /* names of fluent species according to mech. */
  int  *Yinit;      /* initial species mass fractions (optional) */
  double *Wmol;     /* molecular weight of Fluent species */
  int    nEGR;         /* grid points in manifold */
  char **EGRspecies;   /* main initial EGR species */
  double *Y_EGR;       /* initial EGR mass fractions */
  int     nManZ;
  double *Zman;
  double **Yin, **Yex; /* intake/exhaust manifold species */
  double ***Yinert;
  int    iFuel, iO2, iCO2, iH2O, iN2;
  double Mw_fuel, Mw_CO2, Mw_H2O, Mw_O2, Mw_N2;

  /* data required for Wiebe function */
  double LHV;      /* lower heating value */
  double m_fuel;   /* mass of fuel */
  double m_cyl;    /* total mass in the system */

  /* data for initialising a second injection */
  int    nZ2, nFuel2;
  double Tfuel2, *Yfuel2;
  char   **fuel2;

  /* options for pdf and scalar dissipation rate treatment */
  char *chiForm1D, *chiForm2D, *chiMode1D, *chiMode2D;
  char *chiNorm1D, *chiNorm2D, *pdf1D, *pdf2D;
  int  chiForm[3], chiMode[3], chiNorm[3], pdfMode[3];
  double varScale;             /* scales the variance by a factor */
  double chiScale;             /* scales the scalar dissipation */

  /* variables required for 2D enthalyp/mixture flamelet */
  int    nCondZ;
  double *condZ, *src_condZ;

  /* ----- information about threads to solve in the system ----- */
  int      nThreadSolve;         /* cylinder */
  int     *solve_thread_ids;
  char   **solve_thread_names;
  Thread **solve_threads;

  int      nThreadIntake;        /* intake */
  int     *intake_thread_ids;
  char   **intake_thread_names;
  Thread **intake_threads;

  int      nThreadExhaust;       /* exhaust */
  int     *exhaust_thread_ids;
  char   **exhaust_thread_names;
  Thread **exhaust_threads;

  /* I/O variables */
  char *filePrefix, *caseDir;

  /* timing information */
  int      nTimer;
  char     timingFile[1000];
  udfTimer *timer;

  /* ignition model */
  double rk,mk,Tk,Tad,t,Diff;
};
typedef struct MDUC_UDF_STRUCT udfStruct;

/* define the species that will be used for the fuel source */
#define C_FUEL1_SRC(c,t) C_DPMS_YI(c,t,0)
#define C_FUEL2_SRC(c,t) C_DPMS_YI(c,t,1)

#endif /* MDUC_FLUENT_UDF_H */
