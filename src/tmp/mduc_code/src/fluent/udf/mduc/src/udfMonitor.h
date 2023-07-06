#ifndef MDUC_UDF_MONITOR_H
#define MDUC_UDF_MONITOR_H

#include "mducUDF.h"
#include "dx.h"
#include "surf.h"
#include "cxsurf.h"
#include "cxiface.h"


typedef struct MDUC_MONFILE_STRUCT {
  /* number of zones to monitor and corresponding names */
  char   id[256];
  FILE  *p;
  int    n;
  double **vars;
  char   **names;
};
typedef struct MDUC_MONFILE_STRUCT monFileStruct;

typedef struct MDUC_MONZONE_STRUCT {

  /* number of zones to monitor and corresponding names */
  char     name[256];
  int      nThreads;
  char   **threadNames;
  Thread **threads;
  int     *SID;

  /* struct to store information about the files in this zone */
  monFileStruct *out;

  /* ----- variables that will be monitored by default ----- */
  /* states */
  double P;            /* volume averaged pressure in each zone */
  double T;            /* mass averaged temperature */
  double H;            /* mass averaged enthalpy */
  double HRrate;       /* integrated heat release */
  double HRprem;       /* integrated heat release from premixed combustion */
  double rho;          /* mass averaged density */
  double mass;         /* total mass in each zone */
  double massU;        /* Unburned mass */
  double massB;        /* Burned mass */
  double massBprem;    /* Burned mass premixed*/
  double massFl;       /* mass in flame*/
  double WMixU;        /* mass weighted molecular weight in unburned zone */
  double WMixB;        /* mass weighted molecular weight in burned zone */
  double vol;          /* system volume */
  double area;         /* surface area */
  double *Yavg, Ysum;  /* mass of species in each zone */
  double Rgas, gamma;  /* gas properties */
  double cp;
  double q, qdot;      /* average/integrated heat flux */
  double mdot;         /* integrated mass flow */

  /* flow properties */
  double k, eps;      /* turbulent parameters */

  /* keep track of min/max states */
  double Pm[2], Pm_c[2][ND_ND];
  double Tm[2], Tm_c[2][ND_ND];
  double Hm[2], Hm_c[2][ND_ND];

  /* user defined variables */
  double *uds, *udsMin, *udsMax;  /* user-define scalars */
  double *udm, *udmMin, *udmMax;  /* user-defined memory */
  double fuel1, fuel2;      /* fuel vapour evaporated from each injection */
};
typedef struct MDUC_MONZONE_STRUCT monZoneStruct;

typedef struct MDUC_MON_STRUCT {

  /* number of zones to monitor and corresponding names */
  int    nM;        /* number of different monitor files */
  int    nW, nS;
  int    nZones, nWalls, nSurf;
  char **zoneNames, **wallNames, **surfNames;
  char  *suffix;

  /* */
  int nMassAvg, nMassInt;
  int nVolAvg, nVolInt;
  char *massAvgVars, *massIntVars;
  char *volAvgVars, *volIntVars;

  /* a separate struct for all the variables */
  monZoneStruct *zone, *wall, *surf;

  /* additional monitors (not zone based) */
  monFileStruct udfVars;   /* RIF model variables */ 
  monFileStruct residual;  /* equation residuals */
};
typedef struct MDUC_MON_STRUCT monStruct;

/* initialization functions */
void UDFinitMonitor();
void UDFinitResMonitor();
int  UDFmonitorGetThreads(monZoneStruct *mZ);
int  UDFinitThreadNames(int nList, char **list, monZoneStruct *m);

/* functions to write the monitors to a file */
void UDFwriteMonitors();
void UDFcheckMonitorZones();
void UDFinitMonitorFile(int r, char *zoneName, char *suffix, monFileStruct *mF);
void UDFwriteMonitorLine(monFileStruct *mF);
void UDFwriteResidualLine();

/* functions to compute the values to monitor */
double UDFmonitorStandardFluid(monZoneStruct *mZ);
double UDFmonitorStandardWall(monZoneStruct *mZ);
double UDFmonitorStandardSurface(monZoneStruct *mZ);

int  UDFgetSurfaces(int n, char **names, int *sid);
void facet_area_3D(Surface *s, int i, double *area);

#endif /* MDUC_UDF_MONITOR_H */
