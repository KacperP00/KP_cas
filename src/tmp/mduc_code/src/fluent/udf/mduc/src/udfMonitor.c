#include "udfMonitor.h"
#include "udfTools.h"
#include "udfRIF.h"
#include "udfPremixed.h"
#include "accessMDUC.h"

extern udfStruct udf;
static monStruct *mon;

DEFINE_ON_DEMAND(initMonitor)
{
  UDFinitMonitor();
}
/* -------------------------------------------------------------------------- */

void UDFinitMonitor()
{
  Message0(" ----- Initializing custom monitors -----\n");
  int i, j, k;
  mon = (monStruct*) calloc(1, sizeof(monStruct));
  monZoneStruct *mZ;
  monFileStruct *mF;

  /* get the number of zones and their threads from the environment */
  if (RP_is_Defined("mduc/monitor/zone-names"))
    mon->nZones = parseFluentString(RP_Get_String("mduc/monitor/zone-names"),
                                    &mon->zoneNames);

  if (RP_is_Defined("mduc/case-id"))
    mon->suffix = RP_Get_String("mduc/case-id");
  else {
    mon->suffix = (char*) calloc(8, sizeof(char));
    strcpy(mon->suffix, "01");
  }
  
  int nThreadList = 0;
  char **threadList;
  if (RP_is_Defined("mduc/monitor/zone-threads"))
    nThreadList = parseFluentString(RP_Get_String("mduc/monitor/zone-threads"),
                                    &threadList);

  /* set initial memory */
  mon->zone = (monZoneStruct*) calloc(mon->nZones, sizeof(monZoneStruct));
  for (i=0; i<mon->nZones; i++) {
    strcpy(mon->zone[i].name, mon->zoneNames[i]);
    mon->zone[i].nThreads    = 0;
    mon->zone[i].threadNames = (char**)   calloc(1, sizeof(char*));
  }
  
  /* initialize the thread zone names */
  char cStart, cEnd;
  int iZ, iT, nT, nZ, len;
  iZ = nZ = 0;
  for (i=0; i<nThreadList; i++) {
    mZ = &mon->zone[iZ];
    mZ->nThreads = mZ->nThreads + 1;
    iT = mZ->nThreads-1;
    
    mZ->threadNames = (char**) realloc(mZ->threadNames,
                                       mZ->nThreads*sizeof(char*));
    mZ->threadNames[iT] = (char*) calloc(128, sizeof(char));
    
    /* get the current string characteristics */
    len    = strlen(threadList[i]);
    cStart = threadList[i][0];
    cEnd   = threadList[i][len-1];
    
    if (cStart == '(') {
      if (cEnd == ')') { /* only one thread in zone */
        strncpy(mZ->threadNames[iT], &threadList[i][1], len-2);
        mZ->threadNames[iT][len-2] = '\0';
        iZ++;
        nZ++;
      }
      else {
        strcpy(mZ->threadNames[iT], &threadList[i][1]);
      }
    }
    else if (cEnd == ')') {
      strncpy(mZ->threadNames[iT], threadList[i], len-1);
      iZ++;
      nZ++;
    }
    else {
      strcpy(mZ->threadNames[iT], threadList[i]);
    }
  } /* end loop over threads */
  
  if (nZ != mon->nZones)
    Message0("Warning: number of thread zones doesn't match number of zones");

  /* get the threads and write an interim summary */
  int nFound;
  Message0("");
  Message0("  Zones initialized: %i\n", mon->nZones);
  for (i=0; i<mon->nZones; i++) {
    mZ = &mon->zone[i];
    Message0("    Threads for zone %i: %s ( ", i, mon->zone[i].name);
    for (j=0; j<mZ->nThreads; j++)
      Message0("%s ", mZ->threadNames[j]);
    Message0(")\n");
    
    /* allocate and get the threads */
    mZ->threads = (Thread**) calloc(mZ->nThreads, sizeof(Thread*));
    nFound = UDFgetThreads(mZ->nThreads, mZ->threadNames, mZ->threads);
    Message0("         threads found: %i\n", nFound);
  } /* end loop over zones */
  
  /* initialize the different monitor files */
  mon->nM = 4;
  for (i=0; i<mon->nZones; i++) {
    mZ = &mon->zone[i]; /* current zone */
    mZ->out = (monFileStruct*) calloc(mon->nM, sizeof(monFileStruct));
    
    mZ->Yavg   = (double*) calloc(udf.nY, sizeof(double));
    mZ->uds    = (double*) calloc(N_UDS, sizeof(double));
    mZ->udsMin = (double*) calloc(N_UDS, sizeof(double));
    mZ->udsMax = (double*) calloc(N_UDS, sizeof(double));
    mZ->udm    = (double*) calloc(N_UDM, sizeof(double));
    mZ->udmMin = (double*) calloc(N_UDM, sizeof(double));
    mZ->udmMax = (double*) calloc(N_UDM, sizeof(double));
    
    /* ----- states file ----- */
    mF = &mon->zone[i].out[0]; /* point to the current monitor file setup */
    strcpy(mF->id, "state");
    mF->n = 12;
    if (udf.rifPremixed)
      mF->n += 7;
    mF->vars  = (double**) calloc(mF->n, sizeof(double*));
    mF->names = (char**)   calloc(mF->n, sizeof(char*));
    for (j=0; j<mF->n; j++)
      mF->names[j] = (char*) calloc(32, sizeof(char));
    k = 0;
    strcpy(mF->names[k], "pressure[Pa]");    mF->vars[k] = &mZ->P;      k++;
    strcpy(mF->names[k], "temperature[K]");  mF->vars[k] = &mZ->T;      k++;
    strcpy(mF->names[k], "density[kg/m^3]"); mF->vars[k] = &mZ->rho;    k++;
    strcpy(mF->names[k], "mass[kg]");        mF->vars[k] = &mZ->mass;   k++;
    strcpy(mF->names[k], "volume[m^3]");     mF->vars[k] = &mZ->vol;    k++;
    strcpy(mF->names[k], "Rgas[J/(kg-K)]");  mF->vars[k] = &mZ->Rgas;   k++;
    strcpy(mF->names[k], "gamma[-]");        mF->vars[k] = &mZ->gamma;  k++;
    strcpy(mF->names[k], "cp[J/(kg-K)]");    mF->vars[k] = &mZ->cp;     k++;
    strcpy(mF->names[k], "TKE[m^2/s^2]");    mF->vars[k] = &mZ->k;      k++;
    strcpy(mF->names[k], "EPS[m^2/s^3]");    mF->vars[k] = &mZ->eps;    k++;
    strcpy(mF->names[k], "H[J/kg]");         mF->vars[k] = &mZ->H;      k++;
    strcpy(mF->names[k], "HR-rate[W]");      mF->vars[k] = &mZ->HRrate; k++;
    if (udf.rifPremixed) {
      strcpy(mF->names[k], "HR-prem[W]");    mF->vars[k] = &mZ->HRprem; k++;
      strcpy(mF->names[k], "massU[kg]");      mF->vars[k] = &mZ->massU;  k++;
      strcpy(mF->names[k], "massB[kg]");      mF->vars[k] = &mZ->massB;  k++;
      strcpy(mF->names[k], "massFl[kg]");     mF->vars[k] = &mZ->massFl; k++;
      strcpy(mF->names[k], "WmixU[g/mol]");   mF->vars[k] = &mZ->WMixU;  k++;
      strcpy(mF->names[k], "WmixB[g/mol]");   mF->vars[k] = &mZ->WMixB;  k++;
      strcpy(mF->names[k], "massB-prem[kg]"); mF->vars[k] = &mZ->massBprem; k++;
    }
    for (j=0;j<mF->n;j++)
      Message0("chamber-state-names %s\n", mF->names[j]);

    /* ----- mass/species file ----- */
    mF = &mon->zone[i].out[1]; /* point to the current monitor file setup */
    strcpy(mF->id, "species");
    mF->n = 1 + udf.nY + 1 + 2;
    mF->vars  = (double**) calloc(mF->n, sizeof(double*));
    mF->names = (char**)   calloc(mF->n, sizeof(char*));
    for (j=0; j<mF->n; j++)
      mF->names[j] = (char*) calloc(32, sizeof(char));
    k = 0;
    strcpy(mF->names[k], "mass[kg]");    mF->vars[k] = &mZ->mass; k++;
    for (j=0; j<udf.nY; j++, k++) { /* loop over species */
      sprintf(mF->names[k], "Y_%s[-]", udf.species[j]);
      mF->vars[k] = &mZ->Yavg[j];
    }
    strcpy(mF->names[k], "Ysum[-]");        mF->vars[k] = &mZ->Ysum;  k++;
    strcpy(mF->names[k], "fuel-vap-1[kg]"); mF->vars[k] = &mZ->fuel1; k++;
    strcpy(mF->names[k], "fuel-vap-2[kg]"); mF->vars[k] = &mZ->fuel2; k++;
    
    /* ----- uds file ----- */
    mF = &mon->zone[i].out[2]; /* point to the current monitor file setup */
    strcpy(mF->id, "mixture");
    mF->n = (N_UDS + N_UDM)*3;
    mF->vars  = (double**) calloc(mF->n, sizeof(double*));
    mF->names = (char**)   calloc(mF->n, sizeof(char*));
    for (j=0; j<mF->n; j++)
      mF->names[j] = (char*) calloc(32, sizeof(char));
    k = 0;
    for (j=0; j<N_UDS; j++, k++) { /* loop over UDS */
      sprintf(mF->names[k], "%s[-]", user_scalar_vars[j].name);
      mF->vars[k] = &mZ->uds[j]; k++;
      sprintf(mF->names[k], "%s-min[-]", user_scalar_vars[j].name);
      mF->vars[k] = &mZ->udsMin[j]; k++;
      sprintf(mF->names[k], "%s-max[-]", user_scalar_vars[j].name);
      mF->vars[k] = &mZ->udsMax[j];
    }
    for (j=0; j<N_UDM; j++, k++) { /* loop over UDM */
      sprintf(mF->names[k], "%s[-]", user_memory_vars[j].name);
      mF->vars[k] = &mZ->udm[j]; k++;
      sprintf(mF->names[k], "%s-min[-]", user_memory_vars[j].name);
      mF->vars[k] = &mZ->udmMin[j]; k++;
      sprintf(mF->names[k], "%s-max[-]", user_memory_vars[j].name);
      mF->vars[k] = &mZ->udmMax[j];
    }
    
    /* ----- bounds ----- */
    mF = &mon->zone[i].out[3]; /* point to the current monitor file setup */
    strcpy(mF->id, "bounds");
    mF->n = 3*(2 + 2*ND_ND);
    mF->vars  = (double**) calloc(mF->n, sizeof(double*));
    mF->names = (char**)   calloc(mF->n, sizeof(char*));
    for (j=0; j<mF->n; j++)
      mF->names[j] = (char*) calloc(32, sizeof(char));
    k = 0;
    strcpy(mF->names[k], "p-min[Pa]"); mF->vars[k] = &mZ->Pm[0]; k++;
    for (j=0; j<ND_ND; j++, k++) {
      sprintf(mF->names[k], "p-min-x%i[m]", j+1);
      mF->vars[k] = &mZ->Pm_c[0][j];
    }
    strcpy(mF->names[k], "p-max[Pa]"); mF->vars[k] = &mZ->Pm[1]; k++;
    for (j=0; j<ND_ND; j++, k++) {
      sprintf(mF->names[k], "p-max-x%i[m]", j+1);
      mF->vars[k] = &mZ->Pm_c[1][j];
    }
    strcpy(mF->names[k], "T-min[K]"); mF->vars[k] = &mZ->Tm[0]; k++;
    for (j=0; j<ND_ND; j++, k++) {
      sprintf(mF->names[k], "T-min-x%i[m]", j+1);
      mF->vars[k] = &mZ->Tm_c[0][j];
    }
    strcpy(mF->names[k], "T-max[K]"); mF->vars[k] = &mZ->Tm[1]; k++;
    for (j=0; j<ND_ND; j++, k++) {
      sprintf(mF->names[k], "T-max-x%i[m]", j+1);
      mF->vars[k] = &mZ->Tm_c[1][j];
    }
    strcpy(mF->names[k], "H-min[K]"); mF->vars[k] = &mZ->Hm[0]; k++;
    for (j=0; j<ND_ND; j++, k++) {
      sprintf(mF->names[k], "H-min-x%i[m]", j+1);
      mF->vars[k] = &mZ->Hm_c[0][j];
    }
    strcpy(mF->names[k], "H-max[K]"); mF->vars[k] = &mZ->Hm[1]; k++;
    for (j=0; j<ND_ND; j++, k++) {
      sprintf(mF->names[k], "H-max-x%i[m]", j+1);
      mF->vars[k] = &mZ->Hm_c[1][j];
    }
  } /* end loop over zones */
  Message0("Finished looping over zones\n");
  
  /* set the UDF variable file */
  mF = &mon->udfVars;
  strcpy(mF->id, "vars");
  mF->n = 4+1;
  mF->names = (char**)   calloc(mF->n, sizeof(char*));
  mF->vars  = (double**) calloc(mF->n, sizeof(double*));
  for (i=0; i<mF->n; i++) mF->names[i] = (char*) calloc(32, sizeof(char));

  k = 0;
  strcpy(mF->names[k], "dPdt[Pa/s]"); mF->vars[k] = &udf.dPdt; k++;
  for (j=1; j<3; j++) {
    strcpy(mF->names[k], "dPhiMin_%i");
    mF->vars[k] = &udf.dPhiMin_dt[j]; k++;
    strcpy(mF->names[k], "dPhiMax_%i");
    mF->vars[k] = &udf.dPhiMax_dt[j]; k++;
  }
#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    UDFinitMonitorFile(0, "udf", mon->suffix, &mon->udfVars);

  /* ----- set up the wall monitors ----- */
  /* get the number of walls and their threads */
  if (RP_is_Defined("mduc/monitor/wall-names"))
    mon->nWalls = parseFluentString(RP_Get_String("mduc/monitor/wall-names"),
                                    &mon->wallNames);
  int nWallList = 0;
  char **wallList;
  if (RP_is_Defined("mduc/monitor/wall-threads"))
    nWallList = parseFluentString(RP_Get_String("mduc/monitor/wall-threads"),
                                  &wallList);
  /* set initial memory */
  mon->wall = (monZoneStruct*) calloc(mon->nWalls, sizeof(monZoneStruct));
  for (i=0; i<mon->nWalls; i++) {
    strcpy(mon->wall[i].name, mon->wallNames[i]);
    mon->wall[i].nThreads    = 0;
    mon->wall[i].threadNames = (char**) calloc(1, sizeof(char*));
  }

  /* parse the thread names */
  UDFinitThreadNames(nWallList, wallList, mon->wall);

  /* get the threads and write an interim summary */
  Message0("");
  Message0("  Walls initialized: %i\n", mon->nWalls);
  for (i=0; i<mon->nWalls; i++) {
    mZ = &mon->wall[i];
    Message0("    Threads for wall %i: %s ( ", i, mon->wall[i].name);
    for (j=0; j<mZ->nThreads; j++)
      Message0("%s ", mZ->threadNames[j]);
    Message0(")\n");
    
    /* allocate and get the threads */
    mZ->threads = (Thread**) calloc(mZ->nThreads, sizeof(Thread*));
    nFound = UDFgetThreads(mZ->nThreads, mZ->threadNames, mZ->threads);
    Message0("         threads found: %i\n", nFound);
  } /* end loop over zones */

  /* initialize the different monitor files */
  mon->nW = 1;
  for (i=0; i<mon->nWalls; i++) {
    mZ = &mon->wall[i]; /* current zone */
    mZ->out = (monFileStruct*) calloc(mon->nW, sizeof(monFileStruct));
    
    mZ->Yavg = (double*) calloc(udf.nY, sizeof(double));
    mZ->uds  = (double*) calloc(N_UDS, sizeof(double));
    mZ->udm  = (double*) calloc(N_UDM, sizeof(double));
    
    /* ----- states file ----- */
    mF = &mon->wall[i].out[0]; /* point to the current monitor file setup */
    strcpy(mF->id, "state");
    mF->n = 8;
    mF->vars  = (double**) calloc(mF->n, sizeof(double*));
    mF->names = (char**)   calloc(mF->n, sizeof(char*));
    for (j=0; j<mF->n; j++)
      mF->names[j] = (char*) calloc(32, sizeof(char));
    k = 0;
    strcpy(mF->names[k], "pressure[Pa]");      mF->vars[k] = &mZ->P;     k++;
    strcpy(mF->names[k], "temperature[K]");    mF->vars[k] = &mZ->T;     k++;
    strcpy(mF->names[k], "density[kg/m^3]");   mF->vars[k] = &mZ->rho;   k++;
    strcpy(mF->names[k], "area[m^2]");         mF->vars[k] = &mZ->area;  k++;
    strcpy(mF->names[k], "mass_flux[kg/m^2]"); mF->vars[k] = &mZ->mass;  k++;
    strcpy(mF->names[k], "mass_flow[kg/s]");   mF->vars[k] = &mZ->mdot;  k++;
    strcpy(mF->names[k], "heat_flux[J/m^2]");  mF->vars[k] = &mZ->q;     k++;
    strcpy(mF->names[k], "heat_flow[J/s]");    mF->vars[k] = &mZ->qdot;  k++;
  }

  /* ----- set up the surface monitors ----- */
  /* get the number of walls and their threads */
  if (RP_is_Defined("mduc/monitor/surf-names"))
    mon->nSurf = parseFluentString(RP_Get_String("mduc/monitor/surf-names"),
                                    &mon->surfNames);
  int nSurfList = 0;
  char **surfList;
  if (RP_is_Defined("mduc/monitor/surf-threads"))
    nSurfList = parseFluentString(RP_Get_String("mduc/monitor/surf-threads"),
                                  &surfList);
  /* set initial memory */
  mon->surf = (monZoneStruct*) calloc(mon->nSurf, sizeof(monZoneStruct));
  for (i=0; i<mon->nSurf; i++) {
    strcpy(mon->surf[i].name, mon->surfNames[i]);
    mon->surf[i].nThreads    = 0;
    mon->surf[i].threadNames = (char**) calloc(1, sizeof(char*));
  }

  /* get the reference names from the rp variables */
  int nCx, nCx2, *cxRefId;
  char **cxRefName, **cxRefIdS;
  if (RP_is_Defined("mduc/monitor/surf-ref-names")) {
    nCx = parseFluentString(RP_Get_String("mduc/monitor/surf-ref-names"),
                            &cxRefName);
  }
  else {
    Message0("Warning: couldn't find reference surface names\n");
  }
  if (RP_is_Defined("mduc/monitor/surf-ref-ids")) {
    nCx2 = parseFluentString(RP_Get_String("mduc/monitor/surf-ref-ids"),
                             &cxRefIdS);
    if (nCx2 != nCx)
      Message0("Warning: no. of reference surf ids and names don't match\n");
    cxRefId = (int*) calloc(nCx, sizeof(int));
    for (i=0; i<nCx; i++) {
      cxRefId[i] = atoi(cxRefIdS[i]);
    }
  }
  else {
    Message0("Warning: couldn't find reference surface ids\n");
  }
  if (mon->nSurf > 0) {
    Message0("Reference surfaces: \n");
    for (i=0; i<nCx; i++) {
      Message0("   %s --> %i\n",cxRefName[i],cxRefId[i]);
    }
  }

  /* parse the thread names */
  UDFinitThreadNames(nSurfList, surfList, mon->surf);

  for (i=0; i<mon->nSurf; i++) {
    Message0(" Threads for %s --> ( ", mon->surf[i].name);
    for (j=0; j<mon->surf[i].nThreads; j++)
      Message0(" %s ", mon->surf[i].threadNames[j]);
    Message0(" )\n");
  }

  /* get the threads and write an interim summary */
  Message0("");
  Message0("  Surfaces initialized: %i\n", mon->nSurf);
  for (i=0; i<mon->nSurf; i++) {
    mZ = &mon->surf[i];
    Message0("    Threads for surface %i: %s ( ", i, mon->surf[i].name);
    for (j=0; j<mZ->nThreads; j++)
      Message0("%s ", mZ->threadNames[j]);
    Message0(")\n");
    
    /* allocate and get the threads */
    mZ->SID = (int*) calloc(mZ->nThreads, sizeof(int));

    nFound = 0;
    for (j=0; j<mZ->nThreads; j++) {
      for (k=0; k<nCx; k++) {
        Message0("ref. Surf: %s %i\n",cxRefName[k],cxRefId[k]);
        if (!strcmp(mZ->threadNames[j], cxRefName[k])) {
          mZ->SID[j] = cxRefId[k];
          Message0("found: %i\n",mZ->SID[j]);
          nFound++;
        }
      }
    }
    Message0("         threads found: %i\n", nFound);
  } /* end loop over zones */
  if (mon->nSurf > 0) {
  	for (i=0; i<nCx; i++)  free(cxRefName[i]); free(cxRefName);
  	for (i=0; i<nCx2; i++) free(cxRefIdS[i]);  free(cxRefIdS);
  	free(cxRefId);
  }

  /* initialize the different monitor files */
  mon->nS = 1;
  for (i=0; i<mon->nSurf; i++) {
    mZ = &mon->surf[i]; /* current zone */
    mZ->out = (monFileStruct*) calloc(mon->nS, sizeof(monFileStruct));
    
    mZ->Yavg = (double*) calloc(udf.nY, sizeof(double));
    mZ->uds  = (double*) calloc(N_UDS, sizeof(double));
    mZ->udm  = (double*) calloc(N_UDM, sizeof(double));
    
    /* ----- states file ----- */
    mF = &mon->surf[i].out[0]; /* point to the current monitor file setup */
    strcpy(mF->id, "state");
    mF->n = 8;
    mF->vars  = (double**) calloc(mF->n, sizeof(double*));
    mF->names = (char**)   calloc(mF->n, sizeof(char*));
    for (j=0; j<mF->n; j++)
      mF->names[j] = (char*) calloc(32, sizeof(char));
    k = 0;
    strcpy(mF->names[k], "pressure[Pa]");      mF->vars[k] = &mZ->P;     k++;
    strcpy(mF->names[k], "temperature[K]");    mF->vars[k] = &mZ->T;     k++;
    strcpy(mF->names[k], "density[kg/m^3]");   mF->vars[k] = &mZ->rho;   k++;
    strcpy(mF->names[k], "area[m^2]");         mF->vars[k] = &mZ->area;  k++;
    strcpy(mF->names[k], "mass_flux[kg/m^2]"); mF->vars[k] = &mZ->mass;  k++;
    strcpy(mF->names[k], "mass_flow[kg/s]");   mF->vars[k] = &mZ->mdot;  k++;
    strcpy(mF->names[k], "heat_flux[J/m^2]");  mF->vars[k] = &mZ->q;     k++;
    strcpy(mF->names[k], "heat_flow[J/s]");    mF->vars[k] = &mZ->qdot;  k++;
  }

  /* set residual file to NULL pointer to check for initialization after
     first iteration */
  mon->residual.p = NULL;

  /* initialize all the files */
#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      for (i=0; i<mon->nZones; i++)
        for (j=0; j<mon->nM; j++)
          UDFinitMonitorFile(0, mon->zone[i].name, mon->suffix,
                             &mon->zone[i].out[j]);

      for (i=0; i<mon->nWalls; i++)
        for (j=0; j<mon->nW; j++)
          UDFinitMonitorFile(0, mon->wall[i].name, mon->suffix,
                             &mon->wall[i].out[j]);

      for (i=0; i<mon->nSurf; i++)
        for (j=0; j<mon->nS; j++)
          UDFinitMonitorFile(0, mon->surf[i].name, mon->suffix,
                             &mon->surf[i].out[j]);
    }

#if !RP_HOST
  /* write the current state*/
  UDFwriteMonitors();
#endif /* !RP_HOST */
  Message0("\n ------------------------------------\n");
}
/* -------------------------------------------------------------------------- */

void UDFwriteMonitors()
{
  int i, j;

  for (i=0; i<mon->nZones; i++) {
    Message0("    writing in zone %s: ", mon->zone[i].name);
    /* compute the standard variables */
    UDFmonitorStandardFluid(&mon->zone[i]);
    
    /* write the initialized monitors */
    for (j=0; j<mon->nM; j++) {
      Message0("%s ", mon->zone[i].out[j].id);
      UDFwriteMonitorLine(&mon->zone[i].out[j]);
    }
    Message0("\n");
  } /* end loop over zones */
  
  for (i=0; i<mon->nWalls; i++) {
    Message0("    writing on wall %s: ", mon->wall[i].name);
    /* compute the standard variables */
    UDFmonitorStandardWall(&mon->wall[i]);
    
    /* write the initialized monitors */
    for (j=0; j<mon->nW; j++) {
      Message0("%s ", mon->wall[i].out[j].id);
      UDFwriteMonitorLine(&mon->wall[i].out[j]);
    }
    Message0("\n");
  } /* end loop over walls */

  for (i=0; i<mon->nSurf; i++) {
    Message0("    writing on surface %s: ", mon->surf[i].name);
    /* compute the standard variables */
    UDFmonitorStandardSurface(&mon->surf[i]);
    
    /* write the initialized monitors */
    for (j=0; j<mon->nS; j++) {
      Message0("%s ", mon->surf[i].out[j].id);
      UDFwriteMonitorLine(&mon->surf[i].out[j]);
    }
    Message0("\n");
  } /* end loop over walls */
}
/* -------------------------------------------------------------------------- */

void UDFcheckMonitorZones()
{
  int i, nFound;
  monZoneStruct *mZ;
  Message0("    checking monitor zones: \n");
  for (i=0; i<mon->nZones; i++) {
    mZ = &mon->zone[i];
    nFound = UDFgetThreads(mZ->nThreads, mZ->threadNames, mZ->threads);
    Message0("    %s -> threads found: %i\n", mon->zone[i].name, nFound);
  } /* end loop over zones */

  for (i=0; i<mon->nWalls; i++) {
    mZ = &mon->wall[i];
    nFound = UDFgetThreads(mZ->nThreads, mZ->threadNames, mZ->threads);
    Message0("    %s -> threads found: %i\n", mon->wall[i].name, nFound);
  } /* end loop over zones */
}
/* -------------------------------------------------------------------------- */

void UDFmonitorCountThreads(int *nCells, int *nWalls)
{
  Domain *d = Get_Domain(1);
  Thread *t;
  int iT, nT = 0;
  *nCells = *nWalls = 0;
  thread_loop(t,d,iT) {
    if (FLUID_CELL_THREAD_P(t)) {
      *nCells = *nCells + 1;
    }
    else if (THREAD_TYPE(t) == THREAD_F_WALL) {
      *nWalls = *nWalls + 1;
    }
  } /* end loop over all threads */
}
/* -------------------------------------------------------------------------- */

void UDFinitMonitorFile(int r, char *zoneName, char *suffix, monFileStruct *mF)
{
  /* generate the file name */
  int i, nT, cur_n, n_comp;
  char fileName[512], *id, fileLine[2056], dummy[2056], *prev;
  long nPrev;

  /* open state file and write header */
  sprintf(fileName, "%s/monitor/%s-%s.%s", udf.caseDir,zoneName,mF->id,suffix);

  /* check to see if the file exists */
  mF->p = NULL;
  mF->p = fopen(fileName, "r");
  if (mF->p) {
    if (r == 1)
      n_comp = N_ITER;
    else
      n_comp = N_TIME;
    
    Message0("  monitor %s exists\n", fileName);
    
    /* read in the file and search for the time step */
    cur_n = -1;
    while (!feof(mF->p)) {
      nPrev = ftell(mF->p);
      fgets(fileLine, 2056, mF->p);
      
      if (fileLine[0] != '#' && fileLine[0] != '\n')
        sscanf(fileLine, "%i %*s", &cur_n);
      
      if (cur_n >= n_comp)
        break;
    } /* end loop over file */
    
      /* go to beginning of file */
    fseek (mF->p, 0, SEEK_SET);
    
    /* read in everything up until the current time step */
    prev = (char*) calloc(nPrev, sizeof(char));
    fread(prev,1,nPrev,mF->p);
    
    /* reopen the file and write the previous information */
    freopen(fileName, "w", mF->p);
    fwrite(prev,1,nPrev,mF->p);
    
    free(prev);
  }
  else { /* create the file new and write the header */
    mF->p = fopen(fileName, "w");
    
    if (r == 1) {
      nT = 2;
      fprintf(mF->p, "niter ");
    }
    else {
      nT = 3;
      fprintf(mF->p, "# nt<1> time<2> ");
      if (in_cylinder_p) {
        fprintf(mF->p, "CAD<3> ");
        nT++;
      }
    }
    for (i=0; i<mF->n; i++)
      fprintf(mF->p, "%s<%i> ", mF->names[i], nT+i);
    fprintf(mF->p, "\n");
  }
  fflush(mF->p);
}
/* -------------------------------------------------------------------------- */

void UDFwriteMonitorLine(monFileStruct *mF)
{
#if PARALLEL
  if (myid != node_zero) return;
#endif /* PARALLEL */

  int i;
  double CAD;
  fprintf(mF->p, "%6i %15.6e ", N_TIME, CURRENT_TIME);
  if (in_cylinder_p) {
    CAD = RP_Get_Real("dynamesh/in-cyn/crank-angle");
    fprintf(mF->p, "%15.6e ", CAD);
  }
  
  for (i=0; i<mF->n; i++)
    fprintf(mF->p, "%15.6e ", mF->vars[i][0]);

  fprintf(mF->p, "\n");
  fflush(mF->p);
}
/* -------------------------------------------------------------------------- */

void UDFwriteResidualLine()
{
  monFileStruct *mF = &mon->residual;

  if (N_ITER < 1)
    return;

  if (mF->p == NULL)
#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    UDFinitResMonitor();

  Domain *d = Get_Domain(1);
  Residual *res_str = d->residual_top;

#if PARALLEL
  if (myid != node_zero) return;
#endif /* PARALLEL */

  int i;
  fprintf(mF->p, "%8i ", N_ITER);

  for (i=0; i<mF->n; i++)
    fprintf(mF->p, "%15.6e ", res_str->res[i]/res_str->res_scale[i]);

  fprintf(mF->p, "\n");
  fflush(mF->p);
}
/* -------------------------------------------------------------------------- */

void UDFinitResMonitor()
{
  Message0(" ##### initializing residuals file -> ");
  int i;
  Domain *d = Get_Domain(1);
  Residual *res_str = d->residual_top;
  strcpy(mon->residual.id, "all");
  mon->residual.n = res_str->noofeqns; /* number of equations solved */
  mon->residual.names = (char**)   calloc(mon->residual.n, sizeof(char*));
  for (i=0; i<mon->residual.n; i++) {
    /* store the names */
    mon->residual.names[i] = (char*) calloc(32, sizeof(char));
    strcpy(mon->residual.names[i], d->equation_labels[res_str->Mapeqns[i]]);
  } /* end loop over equations */

  /* initialize the file */
  UDFinitMonitorFile(1, "residual", mon->suffix, &mon->residual);
  Message0(" finished #####\n");
}
/* -------------------------------------------------------------------------- */

int UDFinitThreadNames(int nList, char **list, monZoneStruct *m)
{
  int i, iZ, nZ, iT, len;
  char cStart, cEnd;
  monZoneStruct *mZ;

  if (nList == 0) return;

  iZ = nZ = 0;
  for (i=0; i<nList; i++) {
    mZ = &m[iZ]; /* current zone */
    mZ->nThreads = mZ->nThreads + 1;
    iT = mZ->nThreads-1;
    
    mZ->threadNames = (char**) realloc(mZ->threadNames,
                                       mZ->nThreads*sizeof(char*));
    mZ->threadNames[iT] = (char*) calloc(128, sizeof(char));
    
    /* get the current string characteristics */
    len    = strlen(list[i]);
    cStart = list[i][0];
    cEnd   = list[i][len-1];

    if (cStart == '(') {
      if (cEnd == ')') { /* only one thread in zone */
        strncpy(mZ->threadNames[iT], &list[i][1], len-2);
        mZ->threadNames[iT][len-2] = '\0';
        iZ++;
        nZ++;
      }
      else {
        strcpy(mZ->threadNames[iT], &list[i][1]);
      }
    }
    else if (cEnd == ')') {
      strncpy(mZ->threadNames[iT], list[i], len-1);
      iZ++;
      nZ++;
    }
    else {
      strcpy(mZ->threadNames[iT], list[i]);
    }
  } /* end loop over threads */

  return mZ->nThreads;
}
/* -------------------------------------------------------------------------- */

/* ##### start functions that only work for serial/compute node ##### */
#if !RP_HOST
int UDFmonitorGetThreads(monZoneStruct *mZ)
{
  Domain *d = Get_Domain(1);
  Thread *t;
  int i, iT, nFound;

  /* set all threads to null */
  for (i=0; i<mZ->nThreads; i++) mZ->threads[i] = NULL;

  /* loop over all threads */
  i = 0;
  nFound = 0;
  thread_loop(t,d,iT) {
    for (i=0; i<mZ->nThreads; i++) {
      if (!strcmp(mZ->threadNames[i], THREAD_NAME(t))) {
        mZ->threads[i] = t;
        nFound++;
      }
    }
  } /* end loop over all threads */
  return nFound;
}
/* -------------------------------------------------------------------------- */

double UDFmonitorStandardFluid(monZoneStruct *mZ)
{
  int i,j,iT,iMin,iMax;
  Thread *t;
  cell_t  c;

  double mass, totMass, totVol, *rwork;
  double massU,massB,massBprem,massFl;
  double *Pmin, *Pmax, **Pmin_c, **Pmax_c;
  double *Tmin, *Tmax, **Tmin_c, **Tmax_c;
  double *Hmin, *Hmax, **Hmin_c, **Hmax_c;
  double **udsMin, **udsMax, ***udsMin_c, ***udsMax_c;
  double **udmMin, **udmMax, ***udmMin_c, ***udmMax_c;
  double WMix,WMixU,WMixB,WMixInv;
#if !PARALLEL
  int myid = 0;
  compute_node_count = 1;
#endif /* !PARALLEL */

  rwork  = (double*)  calloc(1024, sizeof(double));
  Pmin   = (double*)  calloc(compute_node_count, sizeof(double));
  Pmax   = (double*)  calloc(compute_node_count, sizeof(double));
  Pmin_c = (double**) calloc(compute_node_count, sizeof(double*));
  Pmax_c = (double**) calloc(compute_node_count, sizeof(double*));
  Tmin   = (double*)  calloc(compute_node_count, sizeof(double));
  Tmax   = (double*)  calloc(compute_node_count, sizeof(double));
  Tmin_c = (double**) calloc(compute_node_count, sizeof(double*));
  Tmax_c = (double**) calloc(compute_node_count, sizeof(double*));
  Hmin   = (double*)  calloc(compute_node_count, sizeof(double));
  Hmax   = (double*)  calloc(compute_node_count, sizeof(double));
  Hmin_c = (double**) calloc(compute_node_count, sizeof(double*));
  Hmax_c = (double**) calloc(compute_node_count, sizeof(double*));
  udsMin = (double**) calloc(compute_node_count, sizeof(double*));
  udsMax = (double**) calloc(compute_node_count, sizeof(double*));
  udmMin = (double**) calloc(compute_node_count, sizeof(double*));
  udmMax = (double**) calloc(compute_node_count, sizeof(double*));
  for (i=0; i<compute_node_count; i++) {
    Pmin[i] = Pmax[i] = 0.0;
    Tmin[i] = Tmax[i] = 0.0;
    Pmin_c[i] = (double*) calloc(ND_ND, sizeof(double));
    Pmax_c[i] = (double*) calloc(ND_ND, sizeof(double));
    Tmin_c[i] = (double*) calloc(ND_ND, sizeof(double));
    Tmax_c[i] = (double*) calloc(ND_ND, sizeof(double));
    Hmin_c[i] = (double*) calloc(ND_ND, sizeof(double));
    Hmax_c[i] = (double*) calloc(ND_ND, sizeof(double));
    for (j=0; j<ND_ND; j++) {
      Pmin_c[i][j] = Pmax_c[i][j] = Tmin_c[i][j] = Tmax_c[i][j] = 0.0;
      Hmin_c[i][j] = Hmax_c[i][j] = 0.0;
    }

    udsMin[i] = (double*) calloc(N_UDS, sizeof(double*));
    udsMax[i] = (double*) calloc(N_UDS, sizeof(double*));
    udmMin[i] = (double*) calloc(N_UDM, sizeof(double*));
    udmMax[i] = (double*) calloc(N_UDM, sizeof(double*));
  }
  
  /* initialize all variables */
  mZ->P = mZ->T = mZ->H = mZ->rho = mZ->mass = mZ->vol = 0.0;
  mZ->WMixU = mZ->WMixB = 0.0;
  mZ->HRrate = mZ->HRprem = 0.0;
  mZ->cp = mZ->Rgas = mZ->gamma = 0.0;
  mZ->k = mZ->eps = 0.0;
  mZ->fuel1 = mZ->fuel2 = 0.0;
  for (i=0; i<udf.nY; i++) mZ->Yavg[i] = 0.0;
  for (i=0; i<N_UDS; i++) {
    mZ->uds[i]  = 0.0;
    udsMin[myid][i] =  1e20;
    udsMax[myid][i] = -1e20;
  }
  for (i=0; i<N_UDM; i++) {
    mZ->udm[i]  = 0.0;
    udmMin[myid][i] =  1e20;
    udmMax[myid][i] = -1e20;
  }
  Pmin[myid] = Tmin[myid] = Hmin[myid] =  1e20;
  Pmax[myid] = Tmax[myid] = Hmax[myid] = -1e20;
  totMass = totVol = WMixU = WMixB = 0.0;
  massU   = massB  = massFl = massBprem = 0.0; 

  for (iT=0; iT<mZ->nThreads; iT++) {
    t = mZ->threads[iT]; /* current thread */
    if (t == NULL) continue;

    begin_c_loop_int(c,t) {

      WMix = 0.0;
      for (i=0; i<udf.nY; i++) {
        if (udf.indexY[i]>= 0) {
          WMix += C_YI(c,t,i)/udf.Wmol[udf.indexY[i]]/(C_R(c,t)*C_VOLUME(c,t));
        }
      }
      /* integrate the total mass */
      mass = C_R(c,t)*C_VOLUME(c,t);
      mZ->vol += C_VOLUME(c,t);
      
      /* volume average the pressure */
      mZ->P += C_VOLUME(c,t)*ABS_P(C_P(c,t), op_pres);
      
      /* mass average the other quantities */
      mZ->T     += mass*C_T(c,t);
      if (udf.iH >= 0)
        mZ->H     += mass*C_HTOT(c,t);
      mZ->rho   += mass*C_R(c,t);
      mZ->k     += mass*C_K(c,t);
      mZ->eps   += mass*C_D(c,t);
      mZ->gamma += mass*C_GAMMA(c,t);
      mZ->Rgas  += mass*C_RGAS(c,t);
      mZ->cp    += mass*C_CP(c,t);

      /* integrate the total heat release rate */
      if (udf.iHRrate >= 0) {
        if (udf.rifPremixed) {
          if (C_REGIME(c,t) >= 0.0)
            mZ->HRrate += C_VOLUME(c,t)*C_HR_RATE(c,t);
        }
        else
          mZ->HRrate += C_VOLUME(c,t)*C_HR_RATE(c,t);
      }

      if (udf.rifPremixed && udf.iHRprem >= 0) {
        if (C_REGIME(c,t) < 0.0)
          mZ->HRprem += C_VOLUME(c,t)*C_HR_PREM(c,t);
      }
      
      /* compute the mass average of each species */
      for (i=0; i<udf.nY; i++)
        mZ->Yavg[i] += mass*C_YI(c,t,i);
      
      /* compute the user defined variables  */
      for (i=0; i<N_UDS; i++) {
        mZ->uds[i] += mass*C_UDSI(c,t,udsOffset+i);
        if (C_UDSI(c,t,udsOffset+i) < udsMin[myid][i])
          udsMin[myid][i] = C_UDSI(c,t,udsOffset+i);
        if (C_UDSI(c,t,udsOffset+i) > udsMax[myid][i])
          udsMax[myid][i] = C_UDSI(c,t,udsOffset+i);
      }
      for (i=0; i<N_UDM; i++) {
        mZ->udm[i] += mass*C_UDMI(c,t,udmOffset+i);
        if (C_UDMI(c,t,udmOffset+i) < udmMin[myid][i])
          udmMin[myid][i] = C_UDMI(c,t,udmOffset+i);
        if (C_UDMI(c,t,udmOffset+i) > udmMax[myid][i])
          udmMax[myid][i] = C_UDMI(c,t,udmOffset+i);
      }
      
      /* integrate the amount of fuel evaporated */
      if (udf.iZ1 >= 0) mZ->fuel1 += mass*C_Z1MEAN(c,t);
      if (udf.iZ2 >= 0) mZ->fuel2 += mass*C_Z2MEAN(c,t);
      
      /* keep track of min/max of some states */
      if (ABS_P(C_P(c,t), op_pres) < Pmin[myid]) {
        Pmin[myid] = ABS_P(C_P(c,t), op_pres);
        C_CENTROID(Pmin_c[myid],c,t);
      }
      if (ABS_P(C_P(c,t), op_pres) > Pmax[myid]) {
        Pmax[myid] = ABS_P(C_P(c,t), op_pres);
        C_CENTROID(Pmax_c[myid],c,t);
      }
      if (C_T(c,t) < Tmin[myid]) {
        Tmin[myid] = C_T(c,t);
        C_CENTROID(Tmin_c[myid],c,t);
      }
      if (C_T(c,t) > Tmax[myid]) {
        Tmax[myid] = C_T(c,t);
        C_CENTROID(Tmax_c[myid],c,t);
      }
      if (udf.iH >= 0) {
        if (C_HTOT(c,t) < Hmin[myid]) {
          Hmin[myid] = C_HTOT(c,t);
          C_CENTROID(Hmin_c[myid],c,t);
        }
        if (C_HTOT(c,t) > Hmax[myid]) {
          Hmax[myid] = C_HTOT(c,t);
          C_CENTROID(Hmax_c[myid],c,t);
        }
      }
      
      /* accumulate total quantities */
      totMass += mass;
      if (udf.rifPremixed > 0 ) {
        if (C_PROG(c,t) < 0.01) {
          massU += mass;
          WMixU += 1/WMix;
        }
        else if (C_PROG(c,t) > 0.24) {
          massB += mass;
          WMixB += 1/WMix;
        }
        else {
          massFl += mass;
        }
        massBprem += mass*C_PROG(c,t);
      }
      totVol  += C_VOLUME(c,t);
    } end_c_loop_int(c,t);
  } /* end loop over threads */
  
  /* accumulate across processors */
#if RP_NODE /* parallel */
  mZ->P      = PRF_GRSUM1(mZ->P);
  mZ->T      = PRF_GRSUM1(mZ->T);
  mZ->H      = PRF_GRSUM1(mZ->H);
  mZ->HRrate = PRF_GRSUM1(mZ->HRrate);
  mZ->HRprem = PRF_GRSUM1(mZ->HRprem);
  mZ->rho    = PRF_GRSUM1(mZ->rho);
  mZ->k      = PRF_GRSUM1(mZ->k);
  mZ->eps    = PRF_GRSUM1(mZ->eps);
  mZ->gamma  = PRF_GRSUM1(mZ->gamma);
  mZ->Rgas   = PRF_GRSUM1(mZ->Rgas);
  mZ->cp     = PRF_GRSUM1(mZ->cp);
  mZ->fuel1  = PRF_GRSUM1(mZ->fuel1);
  mZ->fuel2  = PRF_GRSUM1(mZ->fuel2);
  totMass    = PRF_GRSUM1(totMass);
  totVol     = PRF_GRSUM1(totVol);
  if (udf.rifPremixed > 0) {
    massU      = PRF_GRSUM1(massU);
    massB      = PRF_GRSUM1(massB);
    massBprem  = PRF_GRSUM1(massBprem);
    massFl     = PRF_GRSUM1(massFl);
    WMixU      = PRF_GRSUM1(WMixU);
    WMixB      = PRF_GRSUM1(WMixB);
  }
  if (udf.nY > 0) PRF_GRSUM(mZ->Yavg, udf.nY, rwork);
  if (N_UDS > 0)  PRF_GRSUM(mZ->uds, N_UDS, rwork);
  if (N_UDM > 0)  PRF_GRSUM(mZ->udm, N_UDM, rwork);
  PRF_GRSUM(Pmin, compute_node_count, rwork);
  PRF_GRSUM(Pmax, compute_node_count, rwork);
  PRF_GRSUM(Tmin, compute_node_count, rwork);
  PRF_GRSUM(Tmax, compute_node_count, rwork);
  PRF_GRSUM(Hmin, compute_node_count, rwork);
  PRF_GRSUM(Hmax, compute_node_count, rwork);
  for (i=0; i<compute_node_count; i++) {
    PRF_GRSUM(Pmin_c[i], ND_ND, rwork);
    PRF_GRSUM(Pmax_c[i], ND_ND, rwork);
    PRF_GRSUM(Tmin_c[i], ND_ND, rwork);
    PRF_GRSUM(Tmax_c[i], ND_ND, rwork);
    PRF_GRSUM(Hmin_c[i], ND_ND, rwork);
    PRF_GRSUM(Hmax_c[i], ND_ND, rwork);
    if (N_UDS > 0) {
      PRF_GRSUM(udsMin[i], N_UDS, rwork);
      PRF_GRSUM(udsMax[i], N_UDS, rwork);
    }
    if (N_UDM > 0) {
      PRF_GRSUM(udmMin[i], N_UDM, rwork);
      PRF_GRSUM(udmMax[i], N_UDM, rwork);
    }
  }
#endif /* RP_NODE */
  mZ->mass = totMass;
  mZ->vol  = totVol;
  if (udf.iProg > 0) {
    mZ->massU     = massU;
    mZ->massB     = massB;
    mZ->massBprem = massBprem;
    mZ->massFl    = massFl;
    mZ->WMixU     = WMixU/massU;
    mZ->WMixB     = WMixB/massB;
  }
  /* complete normalization */
  mZ->P     /= totVol;
  mZ->T     /= totMass;
  mZ->H     /= totMass;
  mZ->rho   /= totMass;
  mZ->k     /= totMass;
  mZ->eps   /= totMass;
  mZ->gamma /= totMass;
  mZ->Rgas  /= totMass;
  mZ->cp    /= totMass;
  for (i=0; i<udf.nY; i++) mZ->Yavg[i] /= totMass;
  mZ->Ysum = 0.0;
  for (i=0; i<udf.nY; i++) mZ->Ysum    += mZ->Yavg[i];
  for (i=0; i<N_UDS; i++)  mZ->uds[i]  /= totMass;
  for (i=0; i<N_UDM; i++)  mZ->udm[i]  /= totMass;

  /* find the minimum state across the zone */
  iMin = iMax = 0;
  for (i=1; i<compute_node_count; i++) {
    if (Pmin[i] < Pmin[iMin]) iMin = i;
    if (Pmax[i] > Pmin[iMax]) iMax = i;
  }
  mZ->Pm[0] = Pmin[iMin]; mZ->Pm[1] = Pmax[iMax];
  for (i=0; i<ND_ND; i++) {
    mZ->Pm_c[0][i] = Pmin_c[iMin][i];
    mZ->Pm_c[1][i] = Pmax_c[iMax][i];
  }

  iMin = iMax = 0;
  for (i=1; i<compute_node_count; i++) {
    if (Tmin[i] < Tmin[iMin]) iMin = i;
    if (Tmax[i] > Tmin[iMax]) iMax = i;
  }
  mZ->Tm[0] = Tmin[iMin]; mZ->Tm[1] = Tmax[iMax];
  for (i=0; i<ND_ND; i++) {
    mZ->Tm_c[0][i] = Tmin_c[iMin][i];
    mZ->Tm_c[1][i] = Tmax_c[iMax][i];
  }

  iMin = iMax = 0;
  for (i=1; i<compute_node_count; i++) {
    if (Hmin[i] < Hmin[iMin]) iMin = i;
    if (Hmax[i] > Hmin[iMax]) iMax = i;
  }
  mZ->Hm[0] = Hmin[iMin]; mZ->Hm[1] = Hmax[iMax];
  for (i=0; i<ND_ND; i++) {
    mZ->Hm_c[0][i] = Hmin_c[iMin][i];
    mZ->Hm_c[1][i] = Hmax_c[iMax][i];
  }

  for (j=0; j<N_UDS; j++) {
    iMin = iMax = 0;
    for (i=1; i<compute_node_count; i++) {
      if (udsMin[i][j] < udsMin[iMin][j]) iMin = i;
      if (udsMax[i][j] > udsMax[iMax][j]) iMax = i;
    } /* end loop over nodes */
    mZ->udsMin[j] = udsMin[iMin][j];
    mZ->udsMax[j] = udsMax[iMax][j];
  } /* end loop over scalars */

  for (j=0; j<N_UDM; j++) {
    iMin = iMax = 0;
    for (i=1; i<compute_node_count; i++) {
      if (udmMin[i][j] < udmMin[iMin][j]) iMin = i;
      if (udmMax[i][j] > udmMax[iMax][j]) iMax = i;
    } /* end loop over nodes */
    mZ->udmMin[j] = udmMin[iMin][j];
    mZ->udmMax[j] = udmMax[iMax][j];
  } /* end loop over scalars */

  /* free local memory */
  free(rwork);
  for (i=0; i<compute_node_count; i++) {
    free(Pmin_c[i]); free(Pmax_c[i]);
    free(Tmin_c[i]); free(Tmax_c[i]);
    free(Hmin_c[i]); free(Hmax_c[i]);
    free(udsMin[i]); free(udsMax[i]);
    free(udmMin[i]); free(udmMax[i]);
  }
  free(Pmin); free(Pmax); free(Pmin_c); free(Pmax_c);
  free(Tmin); free(Tmax); free(Tmin_c); free(Tmax_c);
  free(Hmin); free(Hmax); free(Hmin_c); free(Hmax_c);
  free(udsMin); free(udsMax); free(udmMin); free(udmMax);
}
/* -------------------------------------------------------------------------- */
 
double UDFmonitorStandardWall(monZoneStruct *mZ)
{
  int i, iT;
  Thread *t;
  face_t  f;

  double A[ND_ND], area, totArea;

  mZ->P = mZ->T = mZ->rho = mZ->area = 0.0;
  mZ->k = mZ->eps = 0.0;
  mZ->qdot = mZ->mdot = 0.0;
  totArea = 0.0;
  for (iT=0; iT<mZ->nThreads; iT++) {
    t = mZ->threads[iT]; /* current thread */
    if (t == NULL) continue;

    /* loop over all faces */
    begin_f_loop(f,t) {
        /* ensure faces not counted twice at partition boundaries, always true
           in serial */
        if (PRINCIPAL_FACE_P(f,t)) {
          F_AREA(A, f, t);
          area = NV_MAG(A);

          /* area weighted values */
          mZ->P   += area*ABS_P(F_P(f,t), op_pres);
          mZ->T   += area*F_T(f,t);
          mZ->k   += area*F_K(f,t);
          mZ->eps += area*F_D(f,t);

          /* compute fluxes */
          mZ->qdot += BOUNDARY_HEAT_FLUX(f,t);  /* W */
          mZ->mdot += F_FLUX(f,t);              /* kg/s */

          totArea += area;
        } /* end principal face */
      } end_f_loop(f,t);
  } /* end loop over wall threads */

  /* accumulate across processors */
#if RP_NODE /* parallel */
  mZ->P     = PRF_GRSUM1(mZ->P);
  mZ->T     = PRF_GRSUM1(mZ->T);
  mZ->k     = PRF_GRSUM1(mZ->k);
  mZ->eps   = PRF_GRSUM1(mZ->eps);
  mZ->qdot  = PRF_GRSUM1(mZ->qdot);
  mZ->mdot  = PRF_GRSUM1(mZ->mdot);
  totArea   = PRF_GRSUM1(totArea);
#endif /* RP_NODE */
  mZ->area = totArea;
  mZ->P   /= totArea;
  mZ->T   /= totArea;
  mZ->k   /= totArea;
  mZ->eps /= totArea;
}
/* -------------------------------------------------------------------------- */

double UDFmonitorStandardSurface(monZoneStruct *mZ)
{
  int i, k, iT;
  Thread *t;
  cell_t c;
  Surface *s;
  face_t  f;

  double A[ND_ND], vel[ND_ND], area, totArea, mass, totMass, totQ;

  mZ->P = mZ->T = mZ->rho = mZ->area = 0.0;
  mZ->k = mZ->eps = 0.0;
  mZ->qdot = mZ->mdot = 0.0;
  totArea = totMass = totQ = 0.0;
  
  /* loop over all surface desired surface ideas */
  for (iT=0; iT<mZ->nThreads; iT++) {
    s = SurfaceList + mZ->SID[iT]; /* current thread */
    if (s == NULL) continue;

    for (i=0, k=0; k<s->nf; i+=s->facets[i]+1, k++) {
      if (s->facets[i] >= 3) {
        /* get the facet area */
        facet_area_3D(s, i, A);
        area = NV_MAG(A);

        /* get the cell associated with the */
        c = RP_CELL(s->cells + k);
        t = RP_THREAD(s->cells + k);

        /* get the velocity vector and "mass" (really momentum term) */
        NV_D(vel,=,C_U(c,t),C_V(c,t),C_W(c,t)) ;
        mass = C_R(c,t)*area;

        /* area weighted values */
        mZ->rho += area*C_R(c,t);
        mZ->P   += area*ABS_P(C_P(c,t), op_pres);
        mZ->T   += mass*C_T(c,t);
        mZ->k   += mass*C_K(c,t);
        mZ->eps += mass*C_D(c,t);
        
        /* compute fluxes */
        mZ->qdot += C_R(c,t)*C_H(c,t)*NV_DOT(A, vel);  /* W */
        mZ->mdot += C_R(c,t)*NV_DOT(A, vel);  /* kg/s */
        
        totArea += area;
        totMass += mass;
      }
    }
  } /* end loop over surfaces */

  /* accumulate across processors */
#if RP_NODE /* parallel */
  mZ->P     = PRF_GRSUM1(mZ->P);
  mZ->T     = PRF_GRSUM1(mZ->T);
  mZ->k     = PRF_GRSUM1(mZ->k);
  mZ->eps   = PRF_GRSUM1(mZ->eps);
  mZ->qdot  = PRF_GRSUM1(mZ->qdot);
  mZ->mdot  = PRF_GRSUM1(mZ->mdot);
  totArea   = PRF_GRSUM1(totArea);
  totMass   = PRF_GRSUM1(totMass);
#endif /* RP_NODE */
  mZ->area = totArea;
  mZ->rho /= totArea;
  mZ->P   /= totArea;
  mZ->T   /= totMass;
  mZ->k   /= totMass;
  mZ->eps /= totMass;
  mZ->mass = totMass/totArea;
}
/* -------------------------------------------------------------------------- */

int UDFgetSurfaces(int n, char **names, int *sid)
{
  Domain *d = Get_Domain(1);
  Thread *t;
  Surface *s, *surfaces;
  int i, iT, nFound;

  Message0("Checking %i surfaces\n", n);
  /* Message0("Looking at %i surfaces\n", maxSurfaceId); */

  /* set all threads to null */
  for (i=0; i<n; i++) {
    Message0(" Surface 1: %s\n", names[i]);
    sid[i] = 0;
  }


/*   for (i=0; i<maxSurfaceId; i++) { */
/*     s = SurfaceList + i; /\* current thread *\/ */
/*     Message0("Surface name: %s\n", Surface_Name(i, s)); */
/*   } */

  /* loop over all threads */
/*   i = 0; */
/*   nFound = 0; */
/*   begin_surface_loop(s,iT,SurfaceList); { */
/*    /\*  Message0("Surface_name: %s\n", Surface_Name(iT, s)); *\/ */
/*     /\* for (i=0; i<n; i++) { *\/ */
/* /\*       if (!strcmp(names[i], Surface_Name(iT, s))) { *\/ */
/* /\*         sid[i] = iT; *\/ */
/* /\*         nFound++; *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */
/*   } end_surface_loop(s,iT,surfaces); /\* end loop over all surfaces *\/ */

  return nFound;
}
/* -------------------------------------------------------------------------- */

void facet_area_3D(Surface *s, int i, double *area)
{
   double v[MAX_FACE_NODES][3];
   int k, m;
   area[0] = area[1] = area[2]=0.0;

   for (k=0; k < s->facets[i]; k++) {
     v[k][0] = Get_Surface_Point_Coord(&(s->points[s->facets[i+k+1]]), X_DIM); 
     v[k][1] = Get_Surface_Point_Coord(&(s->points[s->facets[i+k+1]]), Y_DIM); 
     v[k][2] = Get_Surface_Point_Coord(&(s->points[s->facets[i+k+1]]), Z_DIM);
   }
   
   for (k=0; k < s->facets[i]; k++) {
     m = (k+1)%(s->facets[i]);
     area[0] += 0.5*((v[m][1] + v[k][1]) * (v[m][2] - v[k][2]));
     area[1] += 0.5*((v[m][2] + v[k][2]) * (v[m][0] - v[k][0]));
     area[2] += 0.5*((v[m][0] + v[k][0]) * (v[m][1] - v[k][1]));
   }
}
/* -------------------------------------------------------------------------- */

void UDFcalcMassAverage(Domain *d)
{

}
/* -------------------------------------------------------------------------- */

void UDFcalcMassIntegral()
{
  
}
/* -------------------------------------------------------------------------- */

void UDFcalcVolumeAverage()
{

}
/* -------------------------------------------------------------------------- */


void UDFcalcVolumeIntegral()
{
  
}
/* -------------------------------------------------------------------------- */
#endif /* RPHOST */
/* ##### end functions that only work for serial/compute node ##### */
