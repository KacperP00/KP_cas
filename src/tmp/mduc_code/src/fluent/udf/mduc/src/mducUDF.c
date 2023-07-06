#include "mducUDF.h"
#include "udfRIF.h"
#include "udfPremixed.h"
#include "udfScalarDiss.h"
#include "udfTools.h"
#include "udfMonitor.h"
#include "accessMDUC.h"

/* global struct to store UDF options */
udfStruct udf;
static int last_ts = -1;

DEFINE_ON_DEMAND(init)
{
  Message0("\n ==================== Initialising MDUC ====================\n");
  int  i, j, k;
  FILE *pFile;

  /* get the Fluent release */
  Inquire_Release(&udf.ansysRel[0],&udf.ansysRel[1],&udf.ansysRel[2]);
  Message0("UDF library with Fluent Release %i.%i.%i\n",udf.ansysRel[0],
           udf.ansysRel[1],udf.ansysRel[2]);

  /* read all the fluent options */
  initFluentOptions(&udf);

  /* Reserve necessary memory for the library and set the names of the
     user-defined variables for post-processing */
  if (udf.nUDS > N_UDS) Error("Allocate additional UDS memory");
  for (i=0; i<udf.nUDS; i++)
    Set_User_Scalar_Name(udsOffset+i, udf.udsNames[i]);

  if (udf.nUDM > N_UDM) Error("Allocate additional UDM memory");
  for (i=0; i<udf.nUDM; i++)
    Set_User_Memory_Name(udmOffset+i, udf.udmNames[i]);

  /* initialize the timers */
  udf.nTimer = UDFtimerInit(&udf.timer);

  /* initialise the scalar dissipation rate options */
  UDFinitScalarDiss(&udf);

#if !RP_HOST
  Domain *d;
  Thread *t;
  cell_t  c;
  Material *mat, *sp;
  d = Get_Domain(1); /* get the domain from Fluent */

  /* initialize the threads that should be solved */
  initThreadSolve(d);

  /* compute the initial pressure average values */
  int iT;
  double PAvg, totVolume;
  PAvg = totVolume = 0.0;
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;

    begin_c_loop_int(c,t) {
      PAvg      += C_VOLUME(c,t)*ABS_P(C_P(c,t), op_pres);
      totVolume += C_VOLUME(c,t);
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
#if RP_NODE /* parallel */
  totVolume = PRF_GRSUM1(totVolume);
  PAvg      = PRF_GRSUM1(PAvg);
#endif /* RP_NODE */
  PAvg /= totVolume;
  udf.PAvg = PAvg;
  udf.PAvg_M1 = PAvg;
  udf.dPdt = 0.0;

  /* ensure diffusion coefficients are at least 1 */
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      for (i=0; i<3; i++) {
        if (udf.iFDT[i] >= 0)
          C_UDMI(c,t,udmOffset+udf.iFDT[i]) = 1.0;
      }
    } end_c_loop_int(c,t);
  } /* end loop over threads */

  /* ----- allocate the indices for species in CFD mixture ----- */
  mat    = mixture_material(d);
  if (mat != NULL) {
    udf.nY = MIXTURE_NSPECIES(mat);
    Message0("\n Indexing %i species in %s mixture to mechanism:\n", udf.nY, 
             DOMAIN_MATERIAL_NAME(d));
    /* store species names from material */
    udf.species = (char**) calloc(udf.nY, sizeof(char*));
    mixture_species_loop(mat,sp,i) {
      udf.species[i] = (char*) calloc(25, sizeof(char));
      strcpy(udf.species[i], MIXTURE_SPECIE_NAME(mat,i));
      /* convert species to upper case */
      if (strcmp(udf.species[i], "dummy"))
        for (j=0; j<25; j++) udf.species[i][j] = toupper(udf.species[i][j]);
      /* store indices of major species */
      if (!strcmp(udf.species[i], "O2")) {
        udf.iO2  = i;
        udf.Mw_O2 = MATERIAL_PROP(sp, PROP_mwi);
      }
      else if (!strcmp(udf.species[i], "N2")) {
        udf.iN2 = i;
        udf.Mw_N2 = MATERIAL_PROP(sp, PROP_mwi);
      }
      else if (!strcmp(udf.species[i], "H2O")) {
        udf.iH2O = i;
        udf.Mw_H2O = MATERIAL_PROP(sp, PROP_mwi);
      }
      else if (!strcmp(udf.species[i], "CO2")) {
        udf.iCO2 = i;
        udf.Mw_CO2 = MATERIAL_PROP(sp, PROP_mwi);
      }
      else if (!strcmp(udf.species[i], "I-C8H18")) {
        udf.iFuel = i;
        udf.Mw_fuel = MATERIAL_PROP(sp, PROP_mwi);
      }
    } /* end loop over species */

    Message0("Material properties read from Fluent:\n");
    Message0("    O2: %8.3f [g/mol]\n", udf.Mw_O2);
    Message0("   H2O: %8.3f [g/mol]\n", udf.Mw_H2O);
    Message0("   CO2: %8.3f [g/mol]\n", udf.Mw_CO2);
    Message0("  Fuel: %8.3f [g/mol]\n", udf.Mw_fuel);
    Message0("\n");
  }

  /* check whether mixture fraction is present for SDM model */
  if (udf.sdm_model)
    udf.sdm_active = checkSDM(d);

  /* write the binned pdf of the solver variables */
  if (!udf.rif)
    getPhiBounds(d);
  if (udf.iPhi[1] >= 0 && udf.iPhi[2] >= 0)
    UDFwritePDF_2D(d, 0);
#endif /* !RP_HOST */

  /* ----- initialize the activated models ----- */
  double curT;
  if (udf.rif) {
    if (udf.rifPremixed)
      UDFinitPremixed(); /*Initialize table before calling UDFinitRIF() or add table initialization to UDFinitRIF*/
    UDFinitRIF();

    if (in_cylinder_p)
      curT = RP_Get_Real("dynamesh/in-cyn/crank-angle");
    else
      curT = CURRENT_TIME;

    if (udf.restart && curT >= udf.ignTime) {
      udf.ignited = 1;
    }
  }
  else if (udf.premixed) { /* premixed chemtable */
    UDFinitPremixed();
  }

  if (udf.sdm_model) { 
    /* set the model flag to inactive and set a timer */
    udf.nTimer = UDFtimerCreate("SDM", &udf.timer);
    SDMinit();
  }
  if (udf.wiebe) {
    udf.nTimer = UDFtimerCreate("Wiebe", &udf.timer);
    UDFinitWiebe();
  }

  if (udf.rif && udf.wiebe)
    Message0("Warning! RIF and Wiebe models both active.");

  if (udf.monitor) {
    udf.nTimer = UDFtimerCreate("monitor", &udf.timer);
    UDFinitMonitor();
  }

  /* create the file name and write the header */
  strcpy(udf.timingFile, udf.caseDir);
  strcat(udf.timingFile, "/timing");
#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      Message0("intializing the timing file: %s\n", udf.timingFile);
      pFile = fopen(udf.timingFile, "w");
      fprintf(pFile, "iter  Total[s]  ");
      for (i=1; i<udf.nTimer; i++)
        fprintf(pFile, "%s[s] %s[%%] ", udf.timer[i].name, udf.timer[i].name);
      fprintf(pFile, "rest[s] rest[%%]\n");
      fclose(pFile);
      Message0("done\n");
    }

  Message0("\n ===========================================================\n");
}
/* -------------------------------------------------------------------------- */

DEFINE_ADJUST(adjust, d)
{
#if !RP_HOST  
  int status;
  int curr_ts;
  int it_start_ts;
  curr_ts = N_TIME;
  if (last_ts != curr_ts) {
    last_ts = curr_ts;
    it_start_ts = N_ITER;
  }
  /* start the timer for the overall time step here
     should be stopped at the end of and execute-at-end function */

  if (UDFtimerStatus("Total", udf.timer) == 0)
    UDFtimerStart("Total", udf.timer);

  /* set the default iteration values */
  UDFsetScalarDefaultSources(d);

  /* compute terms for modified enthalpy if necessary */
  if (udf.rifEnthalpy && udf.rifDim == 2) {
    computeModH_src(d);
  }

  /* start the adjust timers */
  if (udf.sdm_model) {
    if (UDFtimerStatus("SDM", udf.timer) == 0)
      UDFtimerStart("SDM", udf.timer);
    else
      UDFtimerSet("SDM", udf.timer);
    
    /* compute the single droplet source term */
    if (udf.sdm_active > 0)
      SDMcalcSrc(d, udf.nInj, udf.nNoz, udf.stream);
    
    /* update the timer */
    UDFtimerAdd("SDM", udf.timer);
  }
  
  if (udf.wiebe) {
    if (UDFtimerStatus("Wiebe", udf.timer) == 0)
      UDFtimerStart("Wiebe", udf.timer);
    else
      UDFtimerSet("Wiebe", udf.timer);
    
    /* compute the Wiebe src term */
    WiebeCalcSrc(d);
    
    /* update the timer */
    UDFtimerAdd("Wiebe", udf.timer);
  }

  if (udf.monitor)
    UDFwriteResidualLine();
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

DEFINE_EXECUTE_AT_END(solve)
{
#if !RP_HOST
  Message0("\n ######################## MDUC solver ###################### \n");

  /* error checking for correct initialisation */
  if (udf.rif && udf.mduc == NULL)
    Error("ERROR: flamelet workspace memory is missing!\n");

  if (udf.sdm_model)
    UDFtimerStop("SDM", udf.timer);

  if (udf.wiebe)
    UDFtimerStop("Wiebe", udf.timer);

  int    i, j, k, iT;
  double PAvg, totVolume, totMass;

  /* Fluent variables */
  Domain *d;
  Thread *t;
  cell_t c;

  /* get the domain from FLUENT */
  d = Get_Domain(1);

  /* check if the threads have changed */
  checkThreadSolve(d);

  if (udf.rif)
    getPhiBounds(d);

  /* check whether at what stage of injection it is */
  udf.injection = checkInjection(d);
  if (udf.sdm_model)
    udf.sdm_active = checkSDM(d);

  /* enforce proper bounds on the mean(s) and variance(s) */
  clipVariance(d);

  if ( udf.iHV >= 0 ) {
    Message0("  Clipping enthalpy variance\n");
    clipEnthalpyVariance(d);
  }

  /* compute average values */
  PAvg = totVolume = 0.0;
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;
    
    begin_c_loop_int(c,t) {
      PAvg      += C_VOLUME(c,t)*ABS_P(C_P(c,t), op_pres);
      totVolume += C_VOLUME(c,t);
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
#if RP_NODE /* parallel */
  totVolume = PRF_GRSUM1(totVolume);
  PAvg      = PRF_GRSUM1(PAvg);
#endif /* RP_NODE */
  PAvg /= totVolume;
  udf.PAvg = PAvg;

  Message0(" injection = %i\n", udf.injection);
  Message0(" SDM       = %i\n", udf.sdm_active);
  udf.dPdt = (PAvg - udf.PAvg_M1)/(CURRENT_TIMESTEP);
  Message0(" Current cylinder properties:\n");
  Message0("   Pressure = %15.6e Pa\n", PAvg);
  Message0("   dP/dt    = %15.6e Pa/s\n", udf.dPdt);
  udf.PAvg_M1 = PAvg; /* update previous pressure for next time step */

  /* check ignition if necessary */
  double curT;
  if (udf.premixed || udf.rifPremixed) {
    if (in_cylinder_p)
      curT = RP_Get_Real("dynamesh/in-cyn/crank-angle");
    else
      curT = CURRENT_TIME;

    if (udf.ignited == 0 && curT >= udf.ignTime) {
      Message0(" !!!!! igniting premixed !!!!!\n");
      igniteFlameKernel();
      udf.ignited = 1;
    }
  }

  /* call the RIF solver if turned on */
  if (udf.rif) {
    UDFsolveRIF(d, udf.injection, PAvg);
  }
  else if (udf.premixed) {
    if (udf.ignited)
      UDFpremixedUpdateCFD(d);
  }

  if (udf.wiebe)
    writeWiebe(d);

  /* write the monitors */
  if (udf.monitor) {
    Message0(" writing monitor files\n");
    UDFtimerStart("monitor", udf.timer);
    UDFcheckMonitorZones();
    UDFwriteMonitors();
    UDFtimerStop("monitor", udf.timer);
  }
  
  /* write the binned pdf of the solver variables */
  if (!udf.rif)
    getPhiBounds(d);
  if (udf.iPhi[1] >= 0 && udf.iPhi[2] >= 0)
    UDFwritePDF_2D(d, 0);

  UDFtimerStop("Total", udf.timer);
#if PARALLEL
  PRF_GSYNC();
  if (myid == node_zero)
#endif /* PARALLEL */
    UDFtimerMonitor(udf.timingFile, udf.timer);
#endif /* !RP_HOST */
  Message0(" ########################################################### \n\n");
}
/* ------------------------------------------------------------------------- */

/* ##### start test functions ##### */
DEFINE_ON_DEMAND(testChi1D)
{
  Message0("testing 1D scalar dissipation rate...\n");
  int j;
  double *chi;

  /* error checking for correct initialisation */
  if (udf.mduc == NULL)
    Error("ERROR: flamelet workspace memory is missing!\n");

  Message0("Getting flamelet variables...\n");
  int    nGridMDUC; /* number of grid points in mixture fraction space */
  double Z1ref;     /* reference mixture fraction (usually stoichiometric) */
  double *MDUCgrid; /* mixture fraction grid */

  /* initialise the grid/species information of the flamelet */
  nGridMDUC = MDUCgetGrid(&MDUCgrid, udf.mduc);
  Z1ref     = MDUCgetZref(1, udf.mduc);

  /* allocate memory */
  chi = (double*) calloc(nGridMDUC, sizeof(double));

#if !RP_HOST
  Domain *d;

  d = Get_Domain(1); /* get the domain from FLUENT */

  /* compute scalar dissipation rate based on method from the fluent script */
  if (!strcmp(RP_Get_String("mduc/chi1-mode"), "distribution"))
    UDFcalcChi1D(d, nGridMDUC, MDUCgrid, chi);
  else if (!strcmp(RP_Get_String("mduc/chi1-mode"), "bin"))
    UDFbinChi1D(d, nGridMDUC, MDUCgrid, chi);
  else if (!strcmp(RP_Get_String("mduc/chi1-mode"), "global-stoich"))
    UDFcalcChiPitsch(d, nGridMDUC, MDUCgrid, Z1ref, chi);

#if RP_NODE
  /* if parallel, send result to host */
  /* send using node zero (all nodes should have the same solution) */
  if (myid == 0)
    PRF_CSEND_REAL(node_host, chi, nGridMDUC, myid);
#endif /* RP_NODE */

#endif /* !RP_HOST */

#if !RP_NODE
  /* ---- print chi to a file (on host or serial) ----- */
  Message("Writing scalar dissipation rate to file...\n");
  char fileName[] = "chi1D_";
  strcat(fileName, RP_Get_String("mduc/chi-mode"));
  strcat(fileName, ".dat");

#if RP_HOST
  /* if parallel, receive result from node zero */
  PRF_CRECV_REAL(node_zero, chi, nGridMDUC, node_zero);
#endif /* RP_HOST */

  FILE *pFile = fopen(fileName, "w");
  fprintf(pFile, "# Z chi\n");

  for (j=0; j<nGridMDUC; j++)
    fprintf(pFile, "%8.6f %15.6e\n", MDUCgrid[j], chi[j]);

  fclose(pFile);
#endif /* !RP_NODE */

  free(chi);
  free(MDUCgrid);
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(testChi2D)
{
  int    j, k;
  double **chi1, **chi12, **chi2;

  /* error checking for correct initialisation */
  if (udf.mduc == NULL)
    Error("ERROR: flamelet workspace memory is missing!\n");

  if (MDUCgetNdim(udf.mduc) == 1)
    Error("ERROR: Need to allocate 2D flamelet!\n");

  /* get flamelet variables */
  int    *nGridMDUC; /* number of grid points in mixture fraction space */
  double **MDUCgrid; /* mixture fraction grid */

  /* initialise the grid/species information of the flamelet */
  nGridMDUC = MDUCgetGrid2D(&MDUCgrid, udf.mduc);

  /* allocate memory */
  chi1  = (double**) calloc(nGridMDUC[1], sizeof(double*));
  chi12 = (double**) calloc(nGridMDUC[1], sizeof(double*));
  chi2  = (double**) calloc(nGridMDUC[1], sizeof(double*));
  for (k=0; k<nGridMDUC[1]; k++)
    {
      chi1[k]  = (double*) calloc(nGridMDUC[0], sizeof(double));
      chi12[k] = (double*) calloc(nGridMDUC[0], sizeof(double));
      chi2[k]  = (double*) calloc(nGridMDUC[0], sizeof(double));
    }

#if !RP_HOST
  Domain *d;

  d = Get_Domain(1); /* get the domain from FLUENT */

  /* compute scalar dissipation rate based on method from the fluent script */
  if (!strncmp(RP_Get_String("mduc/chi-2d-mode"), "distribution", 4))
    UDFcalcChi2D(d, nGridMDUC, MDUCgrid, chi1, chi12, chi2);
  else if (!strcmp(RP_Get_String("mduc/chi-2d-mode"), "bin"))
    UDFbinChi2D(d,nGridMDUC[0],nGridMDUC[1],MDUCgrid[0],MDUCgrid[1],chi1,chi2);

#if RP_NODE
  /* if parallel, send result to host */
  /* send using node zero (all nodes should have the same solution) */
  if (myid == 0)
    {
      for (k=0; k<nGridMDUC[1]; k++)
        {
          PRF_CSEND_REAL(node_host, chi1[k], nGridMDUC[0], myid);
          PRF_CSEND_REAL(node_host, chi12[k], nGridMDUC[0], myid);
          PRF_CSEND_REAL(node_host, chi2[k], nGridMDUC[0], myid);
        }
    }
#endif /* RP_NODE */

#endif /* !RP_HOST */

#if !RP_NODE
  /* ---- print chi to a file ----- */
  char fileName[] = "chi2D_";
  strcat(fileName, RP_Get_String("mduc/chi-2d-form"));
  strcat(fileName, "_");
  strcat(fileName, RP_Get_String("mduc/chi-2d-mode"));
  strcat(fileName, ".dat");
  Message0("Writing results to file: %s\n", fileName);

  FILE *pFile = fopen(fileName, "w");
  fprintf(pFile, "# x1 x2 chi1 chi12 chi2\n");

#if RP_HOST
  /* if parallel, receive result from node zero */
  for (k=0; k<nGridMDUC[1]; k++)
    {
      PRF_CRECV_REAL(node_zero, chi1[k], nGridMDUC[0], node_zero);
      PRF_CRECV_REAL(node_zero, chi12[k], nGridMDUC[0], node_zero);
      PRF_CRECV_REAL(node_zero, chi2[k], nGridMDUC[0], node_zero);
    }
#endif /* RP_HOST */

  for (k=0; k<nGridMDUC[1]; k++)
    {
      for (j=0; j<nGridMDUC[0]; j++)
        {
          fprintf(pFile, "%8.6f %8.6f %15.6e %15.6e %15.6e\n", MDUCgrid[0][j],
                  MDUCgrid[1][k], chi1[k][j], chi12[k][j], chi2[k][j]);
        } /* end loop over columns */
      fprintf(pFile, "\n");
    } /* end loop over rows */

  fclose(pFile);
#endif /* !RP_NODE */

  for (k=0; k<nGridMDUC[1]; k++)
    {
      free(chi1[k]);
      free(chi2[k]);
    }
  free(chi1); free(chi12); free(chi2);
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(writeUDS)
{
  Message("writing mixture fraction data to file...\n");
  Domain *d;
  Thread *t;
  cell_t c;

  d = Get_Domain(1); /* get the domain from FLUENT */

  FILE *pFile = fopen("UDSfield.dat", "w");

  int nCells = 0;
  double TINY = 1.0e-9;

  thread_loop_c(t,d)
    {
      begin_c_loop(c,t)
        {
          if ( (C_Z1MEAN(c,t) > TINY) || (C_Z2MEAN(c,t) > TINY) )
            nCells++;
        }
      end_c_loop(c,t) /* end loop over all cells */
    } /* end loop over cell threads */

  fprintf(pFile, "# nCells\n");
  fprintf(pFile, "%i\n", nCells);
  fprintf(pFile, "# Z1Mean Z1Var Z2Mean Z2Var Z12Var\n");

  thread_loop_c(t,d)
    {
      begin_c_loop(c,t)
        {
          if ( (C_Z1MEAN(c,t) > TINY) || (C_Z2MEAN(c,t) > TINY) )
            fprintf(pFile,"%15.6e %15.6e %15.6e %15.6e %15.6e\n",C_Z1MEAN(c,t),
                    C_Z1VAR(c,t), C_Z2MEAN(c,t), C_Z2VAR(c,t), C_Z0VAR(c,t));
        }
      end_c_loop(c,t) /* end loop over all cells */
    } /* end loop over cell threads */

  fclose(pFile);
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(testCFDupdate)
{
  Message("updating domain from flamelet solution...\n");

#if !RP_HOST
  Domain *d;
  Thread *t;
  cell_t c;
  int    i, j, k;

  d = Get_Domain(1); /* get the domain from FLUENT */

  /* check whether at what stage of injection it is */
  udf.injection = checkInjection(d);

#if RP_NODE
  udf.injection = PRF_GIHIGH1(udf.injection); /* get global injection from all nodes */
  node_to_host_int_1(udf.injection);      /* make sure host has injection */
#endif /* RP_NODE */

  /* error checking for correct initialisation */
  if (udf.mduc == NULL)
    Error("ERROR: flamelet workspace memory is missing!\n");

  int nSpecies = MDUCgetNspecies(udf.mduc);

  if (udf.injection == 1)
    {
      int    nGridMDUC; /* number of grid points in mixture fraction space */
      double *MDUCgrid; /* mixture fraction grid */
      nGridMDUC = MDUCgetGrid(&MDUCgrid, udf.mduc);

      double *pdf;
      double **Y, *Temp;

      /* allocate local memory */
      Y    = (double**) calloc(nGridMDUC, sizeof(double*));
      Temp = (double*)  calloc(nGridMDUC, sizeof(double));
      pdf  = (double*)  calloc(nGridMDUC, sizeof(double));
      for (j=0; j<nGridMDUC; j++)
        Y[j] = (double*) calloc(nSpecies, sizeof(double));

      MDUCgetSolution(&Y, &Temp, udf.mduc);

      updateCFD_1D(d,udf.nY,udf.indexY,nGridMDUC,MDUCgrid,Y,Temp,pdf);

      Message0("Successfully updated CFD domain\n");

      for (j=0; j<nGridMDUC; j++) free(Y[j]);
      free(Y); free(Temp); free(MDUCgrid);
    }
  else if (udf.injection == 2)
    {
      int    *nGridMDUC; /* number of grid points in mixture fraction space */
      double **MDUCgrid; /* mixture fraction grid */
      nGridMDUC = MDUCgetGrid2D(&MDUCgrid, udf.mduc);

      double **pdf;
      double ***Y, **Temp;

      /* allocate local memory */
      Y    = (double***) calloc(nGridMDUC[1], sizeof(double**));
      Temp = (double**)  calloc(nGridMDUC[1], sizeof(double*));
      pdf  = (double**)  calloc(nGridMDUC[1], sizeof(double*));
      for (k = 0; k < nGridMDUC[1]; k++)
        {
          Y[k]    = (double**) calloc(nGridMDUC[0], sizeof(double*));
          Temp[k] = (double*)  calloc(nGridMDUC[0], sizeof(double));
          pdf[k]  = (double*)  calloc(nGridMDUC[0], sizeof(double));

          for (j = 0; j < nGridMDUC[0]; j++)
            Y[k][j] = (double*) calloc(nSpecies, sizeof(double));
        }

      MDUCgetSolution(Y, Temp, udf.mduc);

      updateCFD_2D(d,udf.nY,udf.indexY,nGridMDUC,MDUCgrid,Y,Temp,pdf);

      for (k = 0; k < nGridMDUC[1]; k++)
        {
          for (j = 0; j < nGridMDUC[0]; j++)
            free(Y[k][j]);

          free(Y[k]); free(Temp[k]);
        }
      free(Y); free(Temp);
      free(nGridMDUC); free(MDUCgrid[0]); free(MDUCgrid[1]); free(MDUCgrid);
    }
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(testBetaPDF)
{
  Message0("testing beta PDF...\n");

  int i, j, k;
  double ***Y, **Temp;   /* storage for mass fractions and temperature */

  /* ----- get the flamlet variables ---- */
  int    *nGridMDUC;
  double **MDUCgrid;
  int nSpecies = MDUCgetNspecies(udf.mduc);
  /* initialise the grid/species information of the flamelet */
  nGridMDUC = MDUCgetGrid2D(&MDUCgrid, udf.mduc);

  int n1 = nGridMDUC[0];
  int n2 = nGridMDUC[1];
  double *x1 = MDUCgrid[0];
  double *x2 = MDUCgrid[1];
  double dx, dy;

  int    nZ1, nZ2, nZ, counter, nCount = 500, lineCount = 200;
  double *Ymean, **pdf, *rwork;
  int    **pdfBound, pdfFlag;
  double **pdfGlobal, **pdfDen; /* for computing the global PDF */

  double TINY = 1e-6;

  /* debugging variables */
  char fileName[50], prefix[25], header[200], buffer[200];
  char pdf1Mode[50], pdf2Mode[50];

  /* allocate local memory */
  Y         = (double***) calloc(n2, sizeof(double**));
  Temp      = (double**)  calloc(n2, sizeof(double*));
  Ymean     = (double*)  calloc(nSpecies, sizeof(double));
  pdf       = (double**) calloc(n2, sizeof(double*));
  pdfGlobal = (double**) calloc(n2, sizeof(double*));
  pdfDen    = (double**) calloc(n2, sizeof(double*));
  pdfBound  = (int**)    calloc(n2, sizeof(int*));
  rwork     = (double*)  calloc(n1, sizeof(double));

  for (k = 0; k < n2; k++)
    {
      Y[k]         = (double**) calloc(n1, sizeof(double*));
      Temp[k]      = (double*)  calloc(n1, sizeof(double));
      pdf[k]       = (double*) calloc(n1, sizeof(double));
      pdfGlobal[k] = (double*) calloc(n1, sizeof(double));
      pdfDen[k]    = (double*) calloc(n1, sizeof(double));
      pdfBound[k]  = (int*) calloc(2, sizeof(int));

      for (j = 0; j < n1; j++)
        Y[k][j] = (double*) calloc(nSpecies, sizeof(double));
    }

  double sumPDF, meanPDF, meanPDF2, tmpPDF, mass;
  int negPDF;

  strcpy(pdf1Mode, RP_Get_String("mduc/pdf-1d"));
  strcpy(pdf2Mode, RP_Get_String("mduc/pdf-2d"));

  double vFac = 1.0; /* factor to scale the variance by */
  vFac = RP_Get_Real("mduc/variance-scale");

#if !RP_HOST
  Domain *d;
  Thread *t;
  cell_t c;

  d = Get_Domain(1); /* get the domain from FLUENT */
        
  counter = nZ = nZ1 = nZ2 = 0;
  thread_loop_c(t,d)
    {
      begin_c_loop_int(c,t)
        {
          negPDF = 0;
          mass = C_R(c,t)*C_VOLUME(c,t); /* mass of this cell */

          /* ----- compute the mean mixture fraction of all species ----- */
          if ( (C_Z1MEAN(c,t) > TINY) || (C_Z2MEAN(c,t) > TINY) )
            {
              /* compute the mean of all species for this cell
                 this function will also return the pdf */
              pdfFlag =  MDUCcalcMeanY_2D(pdf2Mode, nGridMDUC,MDUCgrid, nSpecies,
                                          C_Z1MEAN(c,t), vFac*C_Z1VAR(c,t),
                                          C_Z2MEAN(c,t), vFac*C_Z2VAR(c,t),
                                          vFac*C_Z0VAR(c,t), Y, Ymean,
                                          pdf, pdfBound, 1);

              /* add the contribution of this cell to the global PDF */
              for (k=0; k<n2; k++)
                {
                  for (j=0; j<n1; j++)
                    {
                      pdfGlobal[k][j] += pdf[k][j]*mass;
                      pdfDen[k][j]    += mass;
                    }
                }

              /* count how many cells had each mixture fraction */              
              if (pdfFlag == 1 || pdfFlag == 11)
                nZ1++;
              else if (pdfFlag == 2 || pdfFlag == 21)
                nZ2++;
              else
                nZ++;
            } /* end 1D for Z2 */
          else
            {
              /* no mixture fraction, update based on the oxidizer species */
              /* essentially delta function at pure oxidizer */
              for (i=0; i < nSpecies; i++)
                Ymean[i] = Y[0][0][i];

              /* add mass of cell to oxidizer */
              pdfGlobal[0][0] += mass;
              pdfDen[0][0]    += mass;
            }
                                        
          counter++;
          if (counter%nCount == 0)
            Message(".");
          if (counter%(lineCount*nCount) == 0)
            Message0("\n");
        }
      end_c_loop_int(c,t) /* end loop over all cells */
    } /* end loop over cell threads */

#if RP_NODE
  /* if parallel, get global pdf contributions from all nodes */
  for (k=0; k<n2; k++)
    {
      PRF_GRSUM(pdfGlobal[k], n1, rwork);
      PRF_GRSUM(pdfDen[k], n1, rwork);
    }
  nZ1     = PRF_GISUM1(nZ1);
  nZ2     = PRF_GISUM1(nZ2);
  nZ      = PRF_GISUM1(nZ);
  counter = PRF_GISUM1(counter);
#endif /* RP_NODE */

  /* normalise the global pdf */
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++)
        {
          if (pdfDen[k][j] > 0.0)
            pdfGlobal[k][j] = pdfGlobal[k][j]/pdfDen[k][j];
          else
            pdfGlobal[k][j] = 0.0;
        }
    }

  Message0("\n%i cells with mixture updated (percent = %6.2f)\n", nZ1+nZ2+nZ,
          (double)(nZ1+nZ2+nZ)/counter*100.0);
  Message0("%i with Z1 (percent = %6.2f)\n", nZ1, (double)(nZ1)/counter*100.0);
  Message0("%i with Z2 (percent = %6.2f)\n", nZ2, (double)(nZ2)/counter*100.0);
  Message0("%i with Z  (percent = %6.2f)\n", nZ,  (double)(nZ)/counter*100.0);

#if RP_NODE
  if (myid == node_zero)
    PRF_SEND_REAL(node_host, pdfGlobal[k], n1, myid);
#endif /* RP_NODE */

#endif /* !RP_HOST */

#if !RP_NODE
  /* write the global PDF for this time step */
  sprintf(fileName, "globalPDF_2D_TEST.dat");
  sprintf(header, "#2 %i %i\n# x1 x2 pdf", n1, n2);

#if RP_HOST
  PRF_RECV_REAL(node_zero, pdfGlobal[k], n1, node_zero);
#endif /* RP_HOST */

  write2DArray(fileName, header, n1, n2, x1, x2, pdfGlobal);

#endif /* !RP_NODE */

  for (k=0; k < n2; k++)
    {
      for (j=0; j < n1; j++) free(Y[k][j]);

      free(pdf[k]); free(pdfGlobal[k]); free(pdfDen[k]); free(pdfBound[k]);
      free(Y[k]); free(Temp[k]);
    }
  free(Y); free(Temp); free(Ymean);
  free(pdf); free(pdfGlobal); free(pdfDen); free(pdfBound);
  free(rwork);
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(calcTotMass)
{
#if !RP_HOST
  Message("Calculating total mass...\n");
  Domain *d;
  Thread *t;
  cell_t c;

  d = Get_Domain(1); /* get the domain from FLUENT */

  int i, nY;
  double cellMass, mass = 0.0;

  char **species;
  /* N-C7H16, dummy,O2, CO2, CO, H2O, N2  */
  nY = RP_Get_Integer("mduc/num-species");
  readFluentString(nY, RP_Get_String("mduc/species-names"), &species);
  double *mSpecies = calloc(udf.nY, sizeof(double));
  double *rwork    = calloc(udf.nY, sizeof(double));
  double mCheck = 0;

  thread_loop_c(t,d)
    {
      begin_c_loop_int(c,t)
        {
          cellMass = C_R(c,t)*C_VOLUME(c,t);
          mass += cellMass;
          
          for (i=0; i<nY; i++)
            mSpecies[i] += C_YI(c,t,i)*cellMass;
        }
      end_c_loop(c,t)
    }

#if RP_NODE
  PRF_GRSUM(mSpecies, nY, rwork);
#endif RP_NODE

  Message0("Total mass in cylinder: %15.6e kg\n", mass);
  Message0("Mass of species:\n");
  for (i=0; i<nY; i++)
    {
      Message0("%10s: %15.6e kg\n", species[i], mSpecies[i]);
      mCheck += mSpecies[i];
    }
  Message0("%15.6e\n", mCheck);

  free(mSpecies);
  free(rwork);
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(adjMass)
{
#if !RP_HOST
  Domain *d = Get_Domain(1);
  adjustMass(RP_Get_Real("mduc/cylinder-mass"), d);
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(adjTemp)
{
#if !RP_HOST
  Domain *d = Get_Domain(1);
  adjustTemp(RP_Get_Real("mduc/initial-temperature"), d);
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(initSwirlVelocity)
{
#if !RP_HOST
  Message0("Initialising swirl velocity with ");
  Domain *d;
  Thread *t;
  cell_t c;

  d = Get_Domain(1); /* get the domain from FLUENT */

  double PI = 3.1415926535897931;
  double r;          /* current radial location in x-y plane */
  double RPM;        /* swirl in rev/min */
  double phi, theta; /* angle or radial */
  double xc[ND_ND];  /* stores centroid of cell */

  double omega;      /* angular velocity (rads/s) */
  double cw;         /* 1.0 for clockwise, -1.0 for counter-clockwise */

  /* determine the direction of rotation */
  cw = RP_Get_Integer("swirl/clockwise");

  /* compute the angular velocity */
  RPM   = RP_Get_Real("swirl/rpm");
  omega = 2*PI/60.0*RPM;
  Message0("omega = %8.6f rad/s\n", omega);

  /* loop over all cells to set the velocity */
  thread_loop_c(t,d)
    {
      begin_c_loop_int(c,t)
        {
          /* get position and angles */
          C_CENTROID(xc,c,t);
          r     = sqrt(xc[0]*xc[0] + xc[1]*xc[1]);
          theta = atan(ABS_D(xc[1]/xc[0]));
          phi   = 0.5*PI - theta; 
          
          /* set velocity */
          C_U(c,t) = cw*r*omega*cos(phi)*SGN_D(xc[1]);
          C_V(c,t) = -1.0*cw*r*omega*sin(phi)*SGN_D(xc[0]);
          C_W(c,t) = 0.0;
        }
      end_c_loop(c,t)
    }
  Message0("Finished swirl initialisation\n");
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(updateSolverThreadNames)
{
#if !RP_HOST
  Message0("\n ===== Updating solver thread names ===== \n");
  Domain *d = Get_Domain(1);
  UDFupdateSolverThreadNames(d);
  Message0("\n ======================================= \n");
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(writePDF_2D)
{
  Message0("\n ===== Binning 2D pdf ===== \n");

  Message0(" Initializing directory name --> ");
  if (RP_is_Defined("mduc/directory"))
    udf.caseDir = RP_Get_String("mduc/directory");
  else { /* write to local directory */
    udf.caseDir = (char*) calloc(32, sizeof(char));
    strcpy(udf.caseDir, "./");
  }
  Message0(" done\n");

  /* set defaults for scalars */
  Message0(" Initializing scalars --> ");
  udf.iPhi[1]  = 0; /* Z1 */
  udf.iPhiV[1] = 1; /* Z1-variance */
  udf.iPhi[2]  = 2; /* Z2 */
  udf.iPhiV[2] = 3; /* Z2-variance */
  udf.iPhiV[0] = 4; /* Z2-variance */
  udf.rifEnthalpy = 0;
  if (RP_is_Defined("mduc/rif/x-coord"))
    checkRIF_coord(RP_Get_String("mduc/rif/x-coord"),
                   &udf.iPhi[1],&udf.iPhiV[1]);
  if (RP_is_Defined("mduc/rif/y-coord"))
    checkRIF_coord(RP_Get_String("mduc/rif/y-coord"),
                   &udf.iPhi[2],&udf.iPhiV[2]);
  Message0(" done\n");

  /* initialize the threads to condition over */
  Message0(" Initializing threads --> ");
  udf.nThreadSolve = -1;
  if (RP_is_Defined("mduc/solve-zones"))
    udf.nThreadSolve = parseFluentString(RP_Get_String("mduc/solve-zones"),
                                         &udf.solve_thread_names);

#if !RP_HOST
  Domain *d = Get_Domain(1);
  initThreadSolve(d);
  Message0(" done\n");

  /* get the domain bounds */
  Message0(" Getting domain bounds --> ");
  getPhiBounds(d);
  Message0(" done\n");

  /* write out the pdf */
  UDFwritePDF_2D(d, 1);
  Message0("\n ========================== \n");
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */
