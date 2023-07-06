#include "udfRIF.h"
#include "udfPremixed.h"
#include "udfTools.h"
#include "udfScalarDiss.h"
#include "accessMDUC.h"
#include "math.h"

extern udfStruct udf;

void UDFinitRIF()
{
  int i,j,k,iT;
  Message0("===== Initializing Representative Interactive Flamelet ===== \n\n");
  /* get the input file name */
  char *fileName, *inName, *exName;
  if (RP_is_Defined("mduc/inputfile"))
    fileName = RP_Get_String("mduc/inputfile");
  else
    Error("Missing MDUC input file specification (mduc/inputfile)");

  /* check to make sure the input file is present     */
  FILE *pFile = NULL;
  pFile = fopen(fileName, "r");
  if (pFile == NULL)
    Error("Flamelet initialisation failed! %s not found\n", fileName);
  else
    fclose(pFile);
  
#if !RP_HOST && PARALLEL && MDUC_MPI
  /* Need to duplicate the communicator in Fluent for some reason */
  MPI_Comm_dup(MPI_COMM_WORLD, &udf.comm);
  MDUCsetPrintProc(udf.comm);
#else
  MDUCsetPrintProc(NULL);
#endif /* PARALLEL */

  /* initialise the working memory */
  Message0("  Initializing MDUC from file: %s\n", fileName);
  udf.mduc = MDUCinit(fileName);

  /* override the file/directory output from the flamelet file */
  char mducDir[1000];
  strcpy(mducDir, udf.caseDir);
  strcat(mducDir, "/mduc");
  Message0("  Setting MDUC directory to: %s\n", mducDir);
  MDUCsetDataDir(mducDir, udf.mduc);
  MDUCsetFilePrefix(udf.filePrefix, udf.mduc);

  if (MDUCgetNdim(udf.mduc) == 1) {
    double *Zgrid;
    MDUCgetGrid(&Zgrid, udf.mduc);
    udf.phiZero[1] = 0.5*(Zgrid[1]-Zgrid[0]);
    udf.phiZero[2] = udf.phiZero[1];
    free(Zgrid);
  }
  else if (MDUCgetNdim(udf.mduc) == 2) {
    int *nZgrid, iHdir;
    double **Zgrid;
    nZgrid = MDUCgetGrid2D(&Zgrid, udf.mduc);
    udf.phiZero[1] = 0.5*(Zgrid[0][1]-Zgrid[0][0]);
    udf.phiZero[2] = 0.5*(Zgrid[1][1]-Zgrid[1][0]);

    iHdir = MDUCgetHdim(udf.mduc);
    Message0("  Enthalpy direction: %i\n", iHdir);

    if (MDUCgetHdim(udf.mduc) != 0) {
      Message0("  Initializing memory for conditional source term\n");
      udf.nCondZ = nZgrid[iHdir-1];
      udf.condZ = (double*) calloc(udf.nCondZ, sizeof(double));
      udf.src_condZ = (double*) calloc(udf.nCondZ, sizeof(double));
      for (j=0; j<udf.nCondZ; j++)
        udf.condZ[j] = Zgrid[iHdir-1][j];
    }
    free(Zgrid[0]); free(Zgrid[1]); free(Zgrid);
  }

  /* determine whether CFD should be initialised from flamelet, or is this
     a restart */
  udf.restart = 0;
  if (RP_is_Defined("mduc/restart"))
    udf.restart = RP_Get_Boolean("mduc/restart");
  else
    Message0("MDUC Warning: restart is not defined\n");
  /* Mixture fraction scaling factor for table lookup */
  double factorZ = 1;
  if (RP_is_Defined("mduc/scale-z")) {
    factorZ = RP_Get_Real("mduc/scale-z");
    Message0("Z scaled for table lookup by %4.3f\n",factorZ);
  }
  /* get information about whether the mass or temperature should be adjusted */
  double adjMass = -1.0, adjT = -1.0;
  if (RP_is_Defined("mduc/cylinder-mass"))
    adjMass = RP_Get_Real("mduc/cylinder-mass");
  if (RP_is_Defined("mduc/initial-temperature"))
    adjT = RP_Get_Real("mduc/initial-temperature");
  if (RP_is_Defined("mduc/init-h-from-table"))
    udf.initHfromTable = RP_Get_Integer("mduc/init-h-from-table");
  if (udf.initHfromTable > 0)
    Message0("Enthalpy will be initialized with tabulated species\n");
  else
    Message0("Enthalpy will be initialized with RIF species\n");
  /* create the timers */
  udf.nTimer = UDFtimerCreate("chi", &udf.timer);
  udf.nTimer = UDFtimerCreate("flamelet", &udf.timer);
  udf.nTimer = UDFtimerCreate("update", &udf.timer);

#if !RP_HOST
 
#if PARALLEL && MDUC_MPI
  /* set info necessary for running parallel */
  MDUCsetParallel(udf.comm, udf.mduc);
#endif /* PARALLEL && MDUC_MPI */

  /* ensure that the time information is up-to-date with FLUENT */
  MDUCupdateTime(CURRENT_TIME, CURRENT_TIMESTEP, N_TIME, udf.mduc);
  Message0("\nSolution time: %8.6f ms\n", CURRENT_TIME*1e3);

  MDUCsetSpeciesIndex(udf.nY, udf.species, &udf.indexY, udf.mduc);
  Message0("\n Fluent species indices in mechanism:\n");
  int specWarn = 0;
  for (i=0; i<udf.nY; i++) {
    Message0("  %2i %10s %4i\n",i,udf.species[i],udf.indexY[i]);
    if (udf.indexY[i] < 0) specWarn = 1;
  }
  if (specWarn)
    Message0("\n  Warning: some species not found in mechanism!\n");
  /* get the molecular weights of all species and store */
  udf.Wmol = (double*) calloc(MDUCgetNspecies(udf.mduc), sizeof(double));
  MDUCgetMolWeightAll(udf.Wmol, udf.mduc);

  int  nSpecies = MDUCgetNspecies(udf.mduc);
  double *cp    = (double*) calloc(nSpecies, sizeof(double));
  double *h     = (double*) calloc(nSpecies, sizeof(double));
  double *YOx   = (double*) calloc(nSpecies, sizeof(double));
  double *YFuel = (double*) calloc(nSpecies, sizeof(double));

  /* get the species at the oxidizer boundary from the flamelet solution */
  MDUCgetPointY(0.0, 0.0, YOx, udf.mduc);
  MDUCgetPointY(1.0, 0.0, YFuel, udf.mduc);

  /* get enthalpy of the fuel vapour at the reference temperature */
  MDUCcomputeThermoData(h, cp, udf.T_drop_ref,udf.mduc);
  for (i=0, udf.H_vap_ref=0.0; i<nSpecies; i++)
    udf.H_vap_ref += YFuel[i]*h[i];
  Message0("\n Fuel properties:\n");
  Message0("   reference temperature: %15.2f K\n", udf.T_drop_ref);
  Message0("   reference enthalpy:    %15.6e J/kg\n", udf.H_vap_ref);

  Domain *d;
  d = Get_Domain(1); /* get the domain from Fluent */
  int nDim = MDUCgetNdim(udf.mduc);

  /* get the min/max bounds in the domain */
  getPhiBounds(d);

  /* TODO: get rid of this for mixture fractions, or put somewhere else  */
  if (udf.iPhi[1] != udf.iH) {
    udf.phiMin[1]   = 0.0;
    udf.phiMax[1]   = 1.0;
    udf.deltaPhi[1] = 1.0;
  }
  if (udf.iPhi[2] >= 0 && udf.iPhi[2] != udf.iH) {
    udf.phiMin[2]   = 0.0;
    udf.phiMax[2]   = 1.0;
    udf.deltaPhi[2] = 1.0;
  }

  if (udf.restart) {
    Message0("\nRestarting times:\n");
    Message0("    Simulation: %15.6e\n", CURRENT_TIME);
    Message0("    Flamelet:   %15.6e\n", MDUCgetTime(0,udf.mduc));
    
    /* get the pressure data from the previous time step read from file */
    udf.PAvg_M1 = MDUCgetPressure(-1, udf.mduc);
    udf.dPdt    = MDUCgetDPDT(udf.mduc);
    Message0("\nProperties from t(n-1):\n");
    Message0("  Pressure = %15.6e\n", udf.PAvg_M1);
    Message0("  dP/dt    = %15.6e\n", udf.dPdt);
    
    /* check whether at what stage of injection it is */
    udf.injection = checkInjection(d);
#if RP_NODE /* parallel */
    udf.injection = PRF_GIHIGH1(udf.injection); /* get global injection */
    if (myid == node_zero)
#endif /* RP_NODE */
      MDUCwriteFlamelet(strcat(mducDir,"restart.rif"), udf.mduc);
    Message0("\n----- %iD flamelet restart successful -----\n",
             MDUCgetNdim(udf.mduc));
  }
  else {
    Thread *t, *t_f;
    cell_t c;
    face_t f;
    
    double PAvg, HFav, TFav;   /* total average quantities */
    double Hmin, Hmax;
    double mass, WmixInv, totVolume, totMass;
    double Tox;
    
    /* adjust the total mass in the cylinder if specified */
    if (adjMass > 0.0)
      adjustMass(adjMass, d);
    
    /* adjust the average cylinder temperature if specified */
    if (adjT > 0.0)
      adjustTemp(adjT, d);
    
    /* allocate memory for computing mean species */
    double *Ymean = (double*) calloc(nSpecies, sizeof(double));
    double *Ytab  = (double*) calloc(udf.nTabSp, sizeof(double));
    int i;
    
    int     n1, n2, *n, pdfFlag;
    double *x1, *x2, **grid2D;
    
    /* get dimension of initialization */
    n1 = MDUCgetGridSize(1, udf.mduc);
    n2 = MDUCgetGridSize(2, udf.mduc);
    
    /* allocate 1D arrays */
    double *pdf1D   = (double*)  calloc(n1, sizeof(double));
    double *Temp1D  = (double*)  calloc(n1, sizeof(double));
    double **Y1D    = (double**) calloc(n1, sizeof(double*));
    int *pdfBound1D = (int*)     calloc(2, sizeof(int));
    for (i=0; i<n1; i++)
      Y1D[i] = (double*) calloc(nSpecies, sizeof(double));
    
    /* allocate 2D arrays */
    double **pdf   = (double**)  calloc(n2, sizeof(double*));
    double **Temp  = (double**)  calloc(n2, sizeof(double*));
    double ***Y    = (double***) calloc(n2, sizeof(double**));
    int **pdfBound = (int**)     calloc(n2, sizeof(int*));
    for (k=0; k<n2; k++) {
      pdf[k]      = (double*)  calloc(n1, sizeof(double));
      Temp[k]     = (double*)  calloc(n1, sizeof(double));
      Y[k]        = (double**) calloc(n1, sizeof(double*));
      pdfBound[k] = (int*)     calloc(2, sizeof(int));
      
      for (j=0; j<n1; j++)
        Y[k][j] = (double*) calloc(nSpecies, sizeof(double));
    }

    if (nDim == 1) { /* 1D initialization */
      Message0("Initializing from 1D solution\n");
      /* get the grid size */
      n1 = MDUCgetGrid(&x1, udf.mduc);
      n2 = MDUCgetGridSize(2, udf.mduc);
      
      /* get the flamelet solution */
      MDUCgetSolution(&Y1D, &Temp1D, udf.mduc);
      Message0("pdf type: %s\n", udf.pdf1D);
    }
    else if (nDim == 2) { /* 2D initialisation */
      Message0("Initializing from 2D solution\n");
      /* get the grid size */
      n = MDUCgetGrid2D(&grid2D, udf.mduc);
      n1 = n[0];      n2 = n[1];
      x1 = grid2D[0]; x2 = grid2D[1];
      
      /* get the flamelet solution */
      MDUCgetSolution(Y, Temp, udf.mduc);
      
      Message0("pdf type: %s\n", udf.pdf2D);
    }
    Message0("nspecies: %i, n1: %i, n2: %i\n", nSpecies, n1, n2);
    
    double Z1mean, Z1var;
    double phiMean[3], phiVar[3], ct_x[4];
    /* if recondition temperature call, function */
    /* adjust_T_cond_dTdZ(d, n1, x1); */
    
    int    counter, nCount = 500, lineCount = 200;

    /* compute the average flow field quantities required for flamelet */
    HFav = TFav = PAvg = totMass = totVolume = 0.0;
    Hmin = 1e20; Hmax = -1e20;
    counter = 0;
    for (iT=0; iT<udf.nThreadSolve; iT++) {
      t = udf.solve_threads[iT]; /* current thread */
      if (t == NULL) continue;
      
      begin_c_loop_int(c,t) {
        /* normalize the co-ordinate */
        phiMean[1] = phiVar[1] = 0.0;
        if (udf.iPhi[1] >= 0) {
          phiMean[1] = (C_PHI1(c,t) - udf.phiMin[1])/udf.deltaPhi[1];
          if (udf.iPhiV[1] >= 0)
            phiVar[1]  = C_PHIV1(c,t)/udf.deltaPhi2[1];
        }
        phiMean[2] = phiVar[2] = 0.0;
        if (udf.iPhi[2] >= 0) {
          phiMean[2] = (C_PHI2(c,t) - udf.phiMin[2])/udf.deltaPhi[2];
          if (udf.iPhiV[2] >= 0)
            phiVar[2]  = C_PHIV2(c,t)/udf.deltaPhi2[2];
        }
        phiVar[0] = 0.0;
        if (udf.iPhiV[0] >= 0)
          phiVar[0]  = C_PHIV0(c,t)/(udf.deltaPhi[1]*udf.deltaPhi[2]);

        if (udf.rifEnthalpy && nDim == 1) {
          /* TODO DEBUG: read species from boundary */
          phiMean[1] = phiVar[1] = 0.0;
        }
        
        /* initialise the species from the flamelet input file */
        for (i=0; i<nSpecies; i++) Ymean[i] = 0.0;
        if (nDim == 1)
          pdfFlag = MDUCcalcMeanY_1D(udf.pdf1D,n1,x1,nSpecies,
                                     phiMean[1],phiVar[1],Y1D,Ymean,
                                     pdf1D, pdfBound1D);
        else if (nDim == 2)
          pdfFlag =  MDUCcalcMeanY_2D(udf.pdf2D, n, grid2D, nSpecies,
                                      phiMean[1], phiVar[1],
                                      phiMean[2], phiVar[2], phiVar[0],Y,Ymean,
                                      pdf, pdfBound, udf.transform);
        
        if (pdfFlag == -10) {
          Message0("pdf not converged cell: %i\n", c);
        }
        
        WmixInv = 0.0;
        for (i=0; i<udf.nY; i++) {
          C_YI(c,t,i) = 0.0;
          if (udf.indexY[i] >= 0) {
            C_YI(c,t,i) = Ymean[udf.indexY[i]];
            WmixInv += C_YI(c,t,i)/udf.Wmol[udf.indexY[i]];
          }
        }

        /* ensure that the density is correct */
        /* C_R(c,t) = ABS_P(C_P(c,t),op_pres)/(C_RGAS(c,t)*C_T(c,t)); */
        if (adjMass < 0.0 || adjT < 0.0) {
          C_R(c,t) = ABS_P(C_P(c,t),op_pres)/
            (UNIVERSAL_GAS_CONSTANT*WmixInv*C_T(c,t));
        }

       	/* compute the mass local cell mass */
       	mass = C_R(c,t)*C_VOLUME(c,t);

        /* INSERT OPTION FROM WHICH MODEL (RIF OR FPVA) ENTHALPY IS INITIALIZED */
        if (udf.rifPremixed && udf.initHfromTable > 0) {
           ct_x[0] = C_Z1MEAN(c,t)/(1-C_EGR(c,t))*factorZ; 
           ct_x[1] = C_EGR(c,t); 
           ct_x[2] = C_Z1VAR(c,t)/(1-C_EGR(c,t))/(1-C_EGR(c,t))*factorZ*factorZ; 
           ct_x[3] = 0.0; 
           MDUClookupChemTable("species", ct_x, Ytab, udf.chemTabPrem);
           for (i=0; i<nSpecies; i++) Ymean[i] = 0.0;
           for (i=0; i<udf.nTabSp; i++) {
              if (udf.iTabSp[i] >= 0)
                 Ymean[udf.iTabSp[i]] = Ytab[i];
           }
        } 
        
       	/* compute the total enthalpy for this cell */
       	C_HTOT(c,t) = 0.0;
       	MDUCcomputeThermoData(h, cp, C_T(c,t), udf.mduc);
       	for (i=0; i<nSpecies; i++)
       	  C_HTOT(c,t) += Ymean[i]*h[i];
       
       	/* during initialization assign total enthalpy to modified enthalpy */
       	if (udf.iHMOD >= 0)
       	  C_HMOD(c,t) += C_HTOT(c,t);
        	
       	/* compute the Favre average enthalpy of the domain */
       	HFav += mass*C_HTOT(c,t);
       	TFav += mass*C_T(c,t);

       	if (C_HTOT(c,t) < Hmin) Hmin = C_HTOT(c,t); 
       	if (C_HTOT(c,t) > Hmax) Hmax = C_HTOT(c,t);

        /* Initialize progress variable */
        if (udf.rifPremixed) { 
          ct_x[0] = C_Z1MEAN(c,t)/(1-C_EGR(c,t))*factorZ; 
          ct_x[1] = C_EGR(c,t); 
          ct_x[2] = C_Z1VAR(c,t)/(1-C_EGR(c,t))/(1-C_EGR(c,t))*factorZ*factorZ; 
          ct_x[3] = 0.0; 
          C_PROG(c,t) = MDUCchemTableProgMin(ct_x, udf.chemTabPrem);
          /* TODO MARCO: Rename C_PROGTOTAL -> Rather C_RIF */
          C_PROGTOTAL(c,t) = MDUCchemTableProgMin(ct_x, udf.chemTabPrem);
        }
        /* END CHANGE MARCO */

        /* volume average the pressure */
        PAvg += ABS_P(C_P(c,t), op_pres)*C_VOLUME(c,t);
        
        /* update total mass/volume */
        totMass   += mass;
        totVolume += C_VOLUME(c,t);
        counter++;
        if (counter%nCount == 0)             Message0(".");
        if (counter%(lineCount*nCount) == 0) Message0("\n");
      } end_c_loop_int(c,t);
    } /* end loop over cell threads */
    
#if RP_NODE
    /* Get the total average quantities from all compute nodes */
    counter   = PRF_GISUM1(counter);
    totMass   = PRF_GRSUM1(totMass);
    totVolume = PRF_GRSUM1(totVolume);
    PAvg      = PRF_GRSUM1(PAvg);
    TFav      = PRF_GRSUM1(TFav);
    HFav      = PRF_GRSUM1(HFav);
    Hmin      = PRF_GRLOW1(Hmin);
    Hmax      = PRF_GRHIGH1(Hmax);
#endif /* RP_NODE */
    PAvg /= totVolume;
    TFav /= totMass;
    HFav /= totMass;
    udf.Hmin = Hmin;
    udf.Hmax = Hmax;

    double tmp = 100.0/counter;
    Message0("\n Initialization summary:\n");
    Message0(" ----- cells updated -----  number ----- percent -----\n");
    Message0("\t cylinder   \t %8i \t %6.2f\n", counter, (double)(counter)*tmp);
    Message0(" -----------------------------------------------------\n");

    /* free memory no longer needed */
    free(pdf1D); free(pdfBound1D); free(x1);
    for (i=0; i<n1; i++) free(Y1D[i]);
    free(Y1D); free(Temp1D); free(Ymean);

    /* set the previous pressure and DPDT */
    udf.PAvg_M1 = PAvg;
    
    /* set the previous min/max bounds in the domain */
    for (i=0; i<3; i++) {
      udf.phiMin_m1[i] = udf.phiMin[i];
      udf.phiMax_m1[i] = udf.phiMax[i];
    }

    if (udf.rifEnthalpy && nDim == 1 && udf.iH == udf.iPhi[1]) {
      /* TODO DEBUG: this should be fixed, shouldn't have to compute enthalpy
         above*/
      udf.phiMin[1]   = Hmin;
      udf.phiMax[1]   = Hmax;
      udf.deltaPhi[1] = Hmax - Hmin;
    }

    /* get the enthalpy bounds if necessary */
    if (udf.iH == udf.iPhi[1]) {
      udf.Hmin = udf.phiMin[1];
      udf.Hmax = udf.phiMax[1];
      udf.delH = udf.deltaPhi[1];
      /*           udf.phiMin[1] -= udf.phiBuf[1]; */
      /*           udf.phiMax[1] += udf.phiBuf[1]; */
    }
    else if (udf.iH == udf.iPhi[2]) {
      udf.Hmin = Hmin;
      udf.Hmax = Hmax;
      udf.delH = Hmax - Hmin;
    }
    else
      getEnthalpyBounds(d);

      /* initialize the filter for the bounds */
/*       double tmp, sumW, sigma = 4.0; */
/*       phiMinHat = (double**) calloc(3, sizeof(double*)); */
/*       phiMaxHat = (double**) calloc(3, sizeof(double*)); */
/*       for (i=0; i<3; i++) { */
/*         phiMinHat[i] = (double*) calloc(udf.fhat.nW, sizeof(double)); */
/*         phiMaxHat[i] = (double*) calloc(udf.fhat.nW, sizeof(double)); */

/*         /\* compute the filter weights *\/ */
/*         udf.fhat.w = (double*) calloc(udf.fhat.nW, sizeof(double)); */
/*         sumW = 0.0; */
/*         for (j=0; j<udf.fhat.nW; j++) { */
/*           tmp = (double)(j*j)/(udf.fhat.nW*udf.fhat.nW)/(sigma*sigma); */
/*           udf.w[j] = 0.39894228*pow(sigma,0.5)*exp(-pow(tmp,2)); */
/*           sumW += udf.w[j]; */
/*         } */
/*         /\* normalize the filter *\/ */
/*         for (j=0; j<udf.fhat.nW; j++) */
/*           udf.w[j] /= sumW; */
        
/*       } */

    Message0("\nCylinder Properties:\n");
    Message0("  Volume      = %15.3f cm^3\n", totVolume*1e6);
    Message0("  Mass        = %15.3f mg\n",  totMass*1e6);
    Message0("  Pressure    = %15.2f bar\n", PAvg*1e-5);
    Message0("  dP/dt       = %15.2f MPa/s\n", udf.dPdt*1e-6);
    Message0("  Temperature = %15.2f K\n", TFav);
    Message0("  Enthalpy    = %15.2f kJ/kg\n", HFav*1e-3);
    Message0("     min      = %15.2f kJ/kg\n", udf.Hmin*1e-3);
    Message0("     max      = %15.2f kJ/kg\n", udf.Hmax*1e-3);
    Message0("  phi1        = %15.2f\n", phiMean[1]);
    Message0("     min      = %15.2f\n", udf.phiMin[1]);
    Message0("     max      = %15.2f\n", udf.phiMax[1]);

    /* this needs to be generalized for enthalpy */
    MDUCinitFlow(CURRENT_TIME, PAvg, udf.phiMin[1], udf.phiMax[1],
                 udf.phiMin[2], udf.phiMax[2], udf.mduc);
    
    Message0("\n----- %iD successfully flamelet initialised -----\n",
             MDUCgetNdim(udf.mduc));
  } /* finish initialisation */
  free(cp); free(h); free(YOx); free(YFuel);
  
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void UDFsolveRIF(Domain *d, int injection, double PAvg)
{
  Message0("\n ***** Solving RIF Combustion ***** \n");
  int i;

  /* get the species at the oxidizer boundary from the flamelet solution */
  int nSpecies = MDUCgetNspecies(udf.mduc);
  double T0, T1;
  double rho0, rho1, wmix0, wmix1, ysum0, ysum1, w0, w1, y0, y1;
  double *Y0 = (double*) calloc(nSpecies, sizeof(double));
  double *Y1 = (double*) calloc(nSpecies, sizeof(double));
  MDUCgetPointY(0.0, 0.0, Y0, udf.mduc);
  MDUCgetPointY(1.0, 0.0, Y1, udf.mduc);
  T0 = MDUCgetPointT(0.0, 0.0, udf.mduc);
  T1 = MDUCgetPointT(1.0, 0.0, udf.mduc);

  ysum0 = ysum1 = w0 = w1 = 0.0;
  for (i=0; i<nSpecies; i++) {
    y0 = MIN_D(MAX_D(Y0[i],0.0),1.0);
    y1 = MIN_D(MAX_D(Y1[i],0.0),1.0);
    w0 += y0/udf.Wmol[i];
    w1 += y1/udf.Wmol[i];
    ysum0 += y0;
    ysum1 += y1;
  }
  wmix0 = ysum0/w0;
  wmix1 = ysum1/w1;

  rho0 = wmix0*udf.PAvg/(8314.34*T0);
  rho1 = wmix1*udf.PAvg/(8314.34*T1);

  Message0("    Conditions at flamelet bounds:\n");
  Message0("        %15s %15s %15s\n", "T[K]", "rho[kg/m^3]","mixtureW[g/mol]");
  Message0("      0 %15.6e %15.6e %15.6e\n", T0, rho0, wmix0);
  Message0("      1 %15.6e %15.6e %15.6e\n", T1, rho1, wmix1);

  /* get the min/max bounds in the domain */
  getPhiBounds(d);

  /* TODO: get rid of this for mixture fractions, or put somewhere else  */
  if (udf.iPhi[1] != udf.iH) {
    udf.phiMin[1]   = 0.0;
    udf.phiMax[1]   = 1.0;
    udf.deltaPhi[1] = udf.deltaPhi2[1] = 1.0;
  }
  if (udf.iPhi[2] != udf.iH) {
    udf.phiMin[2]   = 0.0;
    udf.phiMax[2]   = 1.0;
    udf.deltaPhi[2] = udf.deltaPhi2[2] = 1.0;
  }

  for (i=0; i<3; i++) {
    udf.dPhiMin_dt[i] = (udf.phiMin[i] - udf.phiMin_m1[i])/CURRENT_TIMESTEP;
    udf.dPhiMax_dt[i] = (udf.phiMax[i] - udf.phiMax_m1[i])/CURRENT_TIMESTEP;
  }

  Message0("    Co-ordinate Bounds:\n");
  Message0("      phi1 (%10s)  min:  %15.6e\n", udf.udsNames[udf.iPhi[1]],
           udf.phiMin[1]);
  Message0("             %10s  rate: %15.6e\n", "", udf.dPhiMin_dt[1]);
  Message0("             %10s  max:  %15.6e\n", "", udf.phiMax[1]);
  Message0("             %10s  rate: %15.6e\n", "", udf.dPhiMax_dt[1]);
  Message0("             %10s  min:  %15.6e\n", "variance", udf.phiMinV[1]);
  Message0("             %10s  max:  %15.6e\n", "", udf.phiMaxV[1]);
  if (udf.iPhi[2] >= 0) {
    Message0("      phi2 (%10s)  min:  %15.6e\n",  udf.udsNames[udf.iPhi[2]],
             udf.phiMin[2]);
    Message0("             %10s  rate: %15.6e\n", "", udf.dPhiMin_dt[2]);
    Message0("             %10s  max:  %15.6e\n", "", udf.phiMax[2]);
    Message0("             %10s  rate: %15.6e\n", "", udf.dPhiMax_dt[2]);
    Message0("             %10s  min:  %15.6e\n", "variance", udf.phiMinV[2]);
    Message0("             %10s  max:  %15.6e\n", "", udf.phiMaxV[2]);
    Message0("\n");
  }
  Message0("   CURRENT_TIMESTEP: %15.6e\n", CURRENT_TIMESTEP);

/*   Message0("    Computing dHmax/dt based on adiabatic rate:\n"); */
/*   udf.phiMax[1]     = udf.phiMax_m1[1] + 1.0/rho1*udf.dPdt; */
/*   udf.dPhiMax_dt[1] = (udf.phiMax[1] - udf.phiMax_m1[1])/CURRENT_TIMESTEP; */
/*   Message0("       new rate: %15.6e\n", udf.dPhiMax_dt[1]); */

  /* only for enthalpy */
/*   if (udf.iPhi[1] == 5) { */
/*     if (udf.dPdt >= 0.0) { */
/*       if (udf.dPhiMin_dt[1] < 1.0/rho0*udf.dPdt) { */
/*         fprintf(stderr, "  Clipping dPhi_min/dt: %15.6e %15.6e\n", */
/*                 udf.dPhiMin_dt[1], 1.0/rho0*udf.dPdt); */
/*         fprintf(stderr, "udf.dPhiMin_dt[1]: %15.6e\n", udf.dPhiMin_dt[1]); */
/*         fprintf(stderr, "rho0:    %15.6e\n", rho0); */
/*         fprintf(stderr, "udf.dPdt:       %15.6e\n", udf.dPdt); */
/*         udf.dPhiMin_dt[1] = 1.0/rho0*udf.dPdt; */
/*       } */
/*       if (udf.dPhiMax_dt[1] > 1.0/rho1*udf.dPdt) { */
/*         fprintf(stderr, "  Clipping dPhi_max/dt: %15.6e %15.6e\n", */
/*                 udf.dPhiMax_dt[1], 1.0/rho1*udf.dPdt); */
/*         udf.dPhiMax_dt[1] = 1.0/rho1*udf.dPdt; */
/*       } */
/*     } */
/*     else { /\* negative dP/dt *\/ */
/*       if (udf.dPhiMin_dt[1] > 1.0/rho0*udf.dPdt) { */
/*         fprintf(stderr, "  Clipping dPhi_min/dt (expansion): %15.6e %15.6e\n", */
/*                 udf.dPhiMin_dt[1], 1.0/rho0*udf.dPdt); */
/*         fprintf(stderr, "udf.dPhiMin_dt[1]: %15.6e\n", udf.dPhiMin_dt[1]); */
/*         fprintf(stderr, "rho0:    %15.6e\n", rho0); */
/*         fprintf(stderr, "udf.dPdt:       %15.6e\n", udf.dPdt); */
/*         udf.dPhiMin_dt[1] = 1.0/rho0*udf.dPdt; */
/*       } */
/*       if (udf.dPhiMax_dt[1] < 1.0/rho1*udf.dPdt) { */
/*         fprintf(stderr, "  Clipping dPhi_max/dt (expansion): %15.6e %15.6e\n", */
/*                 udf.dPhiMax_dt[1], 1.0/rho1*udf.dPdt); */
/*         udf.dPhiMax_dt[1] = 1.0/rho1*udf.dPdt; */
/*       } */
/*     } */
/*   } */


  /* compute chi and solve flamelet based on injection */
  if (injection == 10) { /* don't solve reaction pre-injection */
    UDF_rif_1D(d, 1, PAvg, udf.mduc);
  } /* end non-reactive */
  else if (injection == 11) { /* solve pilot injection */
    UDF_rif_1D(d, 0, PAvg, udf.mduc);
  } /* end first injection */
  else if (injection >= 20) { /* solve main injection */
    /* initialize a 2D flamelet if necessary */
    if (MDUCgetNdim(udf.mduc) == 1 && injection != 21) {
      Message0("Initialising a 2D flamelet from 1D solution with:\n");
      Message0("at time step %i\n", N_TIME);
      
      if (udf.nZ2 < 0) udf.nZ2 = MDUCgetGridSize(1, udf.mduc);
      
      if (udf.Tfuel2 < 0)
        Error("MDUC: main injection temperature not set");
      if (udf.fuel2 == NULL || udf.Yfuel2 == NULL)
        Error("MDUC: main injection fuel species not set");
      
      Message0("  n = %i\n", udf.nZ2);
      Message0("  T = %8.2f K\n", udf.Tfuel2);
      Message0("  species   mass-fraction\n");
      for (i=0; i<udf.nFuel2; i++)
        Message0("  %s   %8.6f\n", udf.fuel2[i], udf.Yfuel2[i]);
      
      /* call initialisation function */
      MDUCinit2Dfrom_BC_and_1D(udf.nZ2, udf.nFuel2, udf.fuel2, udf.Yfuel2,
                               udf.Tfuel2,udf.mduc);
          
      /* update the grid tiny value */
      int     *nGridMDUC;
      double **MDUCgrid;
      nGridMDUC = MDUCgetGrid2D(&MDUCgrid, udf.mduc);
      udf.phiZero[2] = 0.5*(MDUCgrid[1][1]-MDUCgrid[1][0]);
      free(nGridMDUC); free(MDUCgrid[0]); free(MDUCgrid[1]); free(MDUCgrid);
#if RP_NODE
#if MDUC_MPI
      MDUCsetParallel(MPI_COMM_WORLD, udf.mduc);
#endif /* MDUC_MPI */
      if (myid == node_zero)
#endif /* RP_NODE */
        MDUCwriteFlamelet(strcat(RP_Get_String("mduc/directory"),
                                 "mduc_init2D.rif"), udf.mduc);
    }

    /* set enthalpy variables if necessary */
    if (udf.rifEnthalpy) {
      MDUCsetAvgSdot(0, udf.nCondZ, udf.src_condZ, udf.mduc);
    }

    /* call appropriate solver */
    if (injection == 20)
      UDF_inert_2D(d, 1, PAvg, udf.mduc);
    else if (injection == 21) /* TODO: not implemented, solve 1D for now*/
      UDF_rif_1D(d, 0, PAvg, udf.mduc);
    else if (injection == 21)
      UDF_rif_2D(d, 0, PAvg, udf.mduc);
  } /* end second injection */
  
  for (i=0; i<3; i++) {
    udf.phiMin_m1[i] = udf.phiMin[i];
    udf.phiMax_m1[i] = udf.phiMax[i];
  }

  /* ----- update the intake/exhaust manifolds ----- */
  Message0("\n ***** RIF Combustion Completed ***** \n");
}
/* -------------------------------------------------------------------------- */

/* ##### start functions that only work for serial/compute node ##### */
#if !RP_HOST
int checkInjection(Domain *d)
{
  Thread *t;
  cell_t c;
  int iT, inj, mix[3], solution;
  double phiMean[3];
  
  /* set mixture to zero */
  mix[1] = mix[2] = 0;

  /* set default to inert 1D */
  inj = 10;

  /* ----- check which mixtures are present ----- */
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;

    begin_c_loop_int(c,t) {
      /* normalize the co-ordinate */
      phiMean[1] = phiMean[2] = 0.0;
      if (udf.iPhi[1] >= 0)
        phiMean[1] = (C_PHI1(c,t)  - udf.phiMin[1])/udf.deltaPhi[1];
      if (udf.iPhi[2] >= 0)
        phiMean[2] = (C_PHI2(c,t)  - udf.phiMin[2])/udf.deltaPhi[2];
      
      if ( phiMean[1] > udf.phiZero[1] ) mix[1] = 1;
      if ( phiMean[2] > udf.phiZero[2] ) mix[2] = 1;
      if (mix[1] + mix[2] == 2) break;
    } end_c_loop_int(c,t);
  } /* end loop over threads  */

#if RP_NODE /* parallel */
   mix[1] = PRF_GIHIGH1(mix[1]);
   mix[2] = PRF_GIHIGH1(mix[2]);
#endif /* RP_NODE */
   mix[0] = mix[1] + mix[2];

  /* ----- check status for solving flamelet ----- */
  solution = 0;
  /*  pilot injection  */
  if (udf.tPilot > 0.0) { /* check if explicitly set */
    if (CURRENT_TIME > udf.tPilot)
      solution = 1;
  }
  else if (udf.tPilot >= -1.0) {
    if (mix[1] >= 1)
      solution = 1;
  }

  /* main injection  */
  if (udf.tMain > 0.0) { /* check if explicitly set */
    if (CURRENT_TIME > udf.tMain)
      solution = 2;
  }
  else if (udf.tMain >= -1.0) {
    if (mix[2] == 2)
      solution = 2;
  } /* end main */
  
  /* TODO: generalize this, currently for DEBUG of enthalpy */
  if (udf.iPhi[1] == udf.iH) {
    solution = 1;
    mix[0]   = 1;
  }
  else if (udf.rifDim == 2 && udf.iPhi[2] == udf.iH) {
    solution = 2;
  }

  /* TODO: do this in cleaner way */
  if (solution == 2 && udf.rifDim == 1) {
    Message0("Warning: using 1D flamelet although two mixture fractions\n");
    solution = 1;
    mix[0]   = 1;
  }

  switch (solution) {
  case 0: /* inert */
    if (mix[0] < 2)
      inj = 10;
    else if (mix[0] == 2)
      inj = 20;
    break;
  case 1: /* reactive 1D */
    if (mix[0] <= 1)
      inj = 11;
    else if (mix[0] == 2)
      inj = 21;
      break;
  case 2: /* reactive 2D */
    inj = 22;
    break;
  }

  return inj;
}
/* ------------------------------------------------------------------------- */

void UDF_rif_1D(Domain *d, int inert, double PAvg, void *mduc)
{
  Message0("\n ----- pilot injection -----\n");
  /* ----- declare local variables ----- */
  int i, j, iT;
  /* Fluent variables */
  Thread *t;
  cell_t c;

  /* Flamelet variables */
  double *chi;         /* 1D scalar dissipation rate */
  int    *pdfBoundG;   /* global pdf non-zero bounds */
  double *pdf;         /* global pdf returned from updated after solve */
  double **Y, *Temp;   /* storage for mass fractions and temperature */
  double Tox, Tfuel;

  /* ----- get the flamlet variables ---- */
  int    n1, nSpecies;
  double Z1ref;
  double *x1;

  /* initialise the grid/species information of the flamelet */
  nSpecies = MDUCgetNspecies(mduc);
  n1       = MDUCgetGrid(&x1, mduc);
  Z1ref    = MDUCgetZref(1, mduc);
  pdfBoundG    = (int*) calloc(2, sizeof(int));
  pdfBoundG[0] = 0;
  pdfBoundG[1] = n1-1;

  /* allocate local memory */
  Y    = (double**) calloc(n1, sizeof(double*));
  Temp = (double*)  calloc(n1, sizeof(double));
  chi  = (double*)  calloc(n1, sizeof(double));
  pdf  = (double*)  calloc(n1, sizeof(double));
  for (j=0; j<n1; j++)
    Y[j] = (double*) calloc(nSpecies, sizeof(double));

  /* ----- get the scalar dissipation rate ----- */
  UDFtimerStart("chi", udf.timer);
  if (udf.chiMode[1] == 1) /* compute a distribution */
    UDFcalcChi1D(d, n1, x1, chi);
  else if (udf.chiMode[1] == 2) /* simple binning */
    UDFbinChi1D(d, n1, x1, chi);
  else if (udf.chiMode[1] == 3) /* Pitsch method */
    UDFcalcChiPitsch(d, n1, x1, Z1ref, chi);
  else
    Message0("Warning: chi1-mode not defined!");
  UDFtimerStop("chi", udf.timer);

  /* ----- solve the 1D flamelet ----- */
  Message0("Solving a 1D flamelet...\n");
  Message0("    Co-ordinate Bounds (RIF):\n");
  Message0("      phi1 (%10s)  min: %15.6e\n", udf.udsNames[udf.iPhi[1]],
           udf.phiMin[1]);
  Message0("             %10s  max: %15.6e\n", "", udf.phiMax[1]);
  Message0("\n");
  UDFtimerStart("flamelet", udf.timer);

  if (inert) { /* Get the cylinder average temperature if inert solution */
    Tox = calcAvgTemp(d);
    MDUCsetTempBounds(Tox,-1.0,-1.0,mduc);
  }

  Message0("**  min: %15.6e, max: %15.6e\n **",udf.phiMin[1],udf.phiMax[1]);
  /* solve flamelet (NOTE: communication now taken care of in library  */
  MDUCresetODE(mduc);
  MDUCflamelet1D(CURRENT_TIMESTEP,PAvg,udf.phiMin[1],udf.phiMax[1],chi,inert,
                 Y,Temp,mduc);

  UDFtimerStop("flamelet", udf.timer);

  /* ----- update the CFD solution ----- */
  UDFtimerStart("update", udf.timer);
  /* CHANGE MARCO - Ensuring same bounds in RIF and FLUENT */
  /* This can cause unphyiscal backreaction close to the wall though, since no heat loss is considered in RIF */
  /* TODO: Should be optional to enforce same enthalpy bounds -> Input flag */
#if RP_HOST
  udf.phiMin[1] = MDUCgetHbounds(-1, udf.mduc);
  udf.phiMax[1] = MDUCgetHbounds( 1, udf.mduc); 
#endif /* RP_HOST */
  /* END CHANGE */
  updateCFD_1D(d, udf.nY, udf.indexY, n1, x1, Y, Temp, pdf);
  UDFtimerStop("update", udf.timer);

  /* for now, write file from node 0.
     According to Fluent, this should be done from the host only.
     However, this is a problem for me because mduc is only up-to-date
     on node. */
#if RP_NODE
  if (myid == node_zero)
#endif /* RP_NODE */
    {
      MDUCwriteFlowData(&chi, NULL, &pdf, NULL, mduc);
      MDUCwriteScalars(CURRENT_TIME,NULL,NULL,&pdf,&pdfBoundG,-1,mduc);
    }

  for (j=0; j<n1; j++) free(Y[j]);
  free(Y); free(Temp); free(chi); free(x1);
  free(pdf); free(pdfBoundG);
  Message0("\n ----- exiting pilot injection ----- \n");
}
/* -------------------------------------------------------------------------- */

/* NOTE: don't think I even need to pass temperature, leave for now
   This will be done on every node */
void updateCFD_1D(Domain *d, int nYFav, int *indexYFav, int n1, double *x1,
                  double **Y, double *Temp, double *pdfGlobal)
{
  Message0("\n ----- Averaging flamelet solution and updating cells ----- \n");
  Thread *t;
  cell_t c;
  face_t f;   
  Thread *t_f;
  d = Get_Domain(1); /* get the domain from Fluent */
  
  int    i, j, k, iT, nA, nZ1, nStoich, nPremix,found;
  /* counters for tracking updated  cells */
  int    counter, nCount = 500, lineCount = 200;
  double *pdf, *Ymean,*YmeanU, *Ytab, *Yun, *rwork;
  double tmp, *fp,*fs;
  double P, p_tab_ref, p_ratio, facP; /*Pressure expansion variables */
  double H, h_tab_ref, h_ratio, facH; /*Temperature expansion variables TODO: Change to facT...*/
  double facCa, facCb, facCM, cRef, cMaxRef, sum, minRef; /* Temperature shift correction variables */
  double progMax,progMin,L_f, gradC_ref, progSRC; /* gradC_ref is not needed anymore I think */
  int    pdfFlag, *pdfBound;
  double mass, totMass;
  double phiMean, phiVar, ct_x[4],ct_Unburned[4], ct_xRef[4], ct_noC[3]; /* TODO: Table dimensions hardcoded! */
  int nSpecies = MDUCgetNspecies(udf.mduc); /* consider passing as argument */
  double sumY,sumYTab;
  int debug,debug2,iTab,npar=11; /*Set Number of Expansion Parameters */
  double GradZ[3], taui;
  double V_b,rk,lt;

  /* for error reporting */
  char errFileName[1024];
  char tabSpI[64];
  FILE *errFile;
 
  /* allocate all local memory */
  pdfBound = (int*) calloc(2,          sizeof(int));
  pdf   = (double*) calloc(n1,         sizeof(double));
  rwork = (double*) calloc(n1,         sizeof(double));
  Ymean = (double*) calloc(nSpecies,   sizeof(double));
  YmeanU= (double*) calloc(nSpecies,   sizeof(double));
  Ytab  = (double*) calloc(udf.nTabSp, sizeof(double));
  Yun   = (double*) calloc(udf.nTabSp, sizeof(double));
  fp    = (double*) calloc(udf.nTabSp*npar, sizeof(double));
  fs    = (double*) calloc(npar, sizeof(double));
  
  /* zero out the global pdf */
  for (j=0; j<n1; j++) pdfGlobal[j] = 0.0;

  double HFav = 0.0;
  double Hmin =  1e20;
  double Hmax = -1e20;
  int jM, kM, Atot, A00, A01, A10, A11;
  int *indexTab;

  double z1Mean,z1VarMean,egrMean;
  double C_starRIF;
  double ZRIF;
  double EGRRIF;
  double PROGRIFmin;
  double PROGRIFmax;
  double RIF_PROG, RIF_SRC;
  double alpha;
  double T_start;
  double factorZ = 1.0;

  /* Store index in udf struct */
  /* if species is stored in table udf.indexTab gives table index j to
     mechanism index i else udf.indexTab = -1 */
  indexTab = (int*) calloc(nSpecies,sizeof(int));
  for (i=0; i<nSpecies; i++) {
    indexTab[i] = -1;
    for (j=0; j<udf.nTabSp;j++) {
      if (udf.iTabSp[j] == i) { 
        indexTab[i] = j;
      }
   }
  }
  mass = totMass = 0;
  if (udf.rifPremixed) {
  	z1Mean  =  z1VarMean = egrMean = V_b = 0;
  	factorZ = RP_Get_Real("mduc/scale-z");
  	if (factorZ>0.0) 
  	  Message0("Z scaled by %4.3f\n for table lookup" ,factorZ);
  	
  	if (RP_Get_Boolean("mduc/avg-mixture")) 
		Message0("Mixture distribution not considered -> taking average values\n");
  }
  /* Comment Marco: Determining Hmin & Hmax here is not necessary, 
 	since it is not used (apart from messages)! */
  for (iT=0; iT<udf.nThreadSolve; iT++) { 
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;
    
    begin_c_loop_int(c,t) {
      mass       = C_R(c,t)*C_VOLUME(c,t);
      HFav      += mass*C_HTOT(c,t);
      totMass   += mass;
      if (udf.rifPremixed) {
      	z1Mean    += mass*C_Z1MEAN(c,t)*factorZ;
      	z1VarMean += mass*C_Z1VAR(c,t)*factorZ*factorZ;
      	if (udf.iEGR > 0)
        	egrMean += mass*C_EGR(c,t);
        /* TODO MARCO: */
        if (C_PROG(c,t) < 0.5) {
          if (C_HTOT(c,t) < Hmin) Hmin = C_HTOT(c,t);
          if (C_HTOT(c,t) > Hmax) Hmax = C_HTOT(c,t);
        }
      	if (C_PROG(c,t) > 0.9) V_b += C_VOLUME(c,t);
      }
      else {
      	if (C_HTOT(c,t) < Hmin) Hmin = C_HTOT(c,t);
      	if (C_HTOT(c,t) > Hmax) Hmax = C_HTOT(c,t);
      }
      /* END CHANGE MARCO */
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
  
#if RP_NODE
  /* if parallel, synchronize all node solutions here by summing the
     contributions from each node */
  HFav    = PRF_GRSUM1(HFav);
  totMass = PRF_GRSUM1(totMass); /* total mass in the cylinder */
  Hmin    = PRF_GRLOW1(Hmin);
  Hmax    = PRF_GRHIGH1(Hmax);
  if (udf.rifPremixed) {
  	V_b     = PRF_GRSUM1(V_b);
  	z1Mean  = PRF_GRSUM1(z1Mean);
  	egrMean = PRF_GRSUM1(egrMean);
  	z1VarMean = PRF_GRSUM1(z1VarMean);
  }
#endif /* RP_NODE */
  HFav   /= totMass;
  /* TODO: fix, this is for symmetric geometries*/
  if (udf.rifPremixed)
  	V_b    *= 2;

  MDUCupdateTime(CURRENT_TIME, CURRENT_TIMESTEP, N_TIME, udf.mduc);
  Message0("Enthalpy Summary\n:");
  Message0("  avg. = %15.3f kJ/kg\n", HFav*1e-3);
  Message0(" min_m1= %15.3f kJ/kg\n", udf.phiMin_m1[1]*1e-3);
  Message0("  min. = %15.3f (%15.3f) kJ/kg\n", Hmin*1e-3, udf.phiMin[1]*1e-3);
  Message0(" max_m1= %15.3f kJ/kg\n", udf.phiMax_m1[1]*1e-3);
  Message0("  max. = %15.3f (%15.3f) kJ/kg\n", Hmax*1e-3, udf.phiMax[1]*1e-3);
  Message0("max(max) . = %15.3f (%15.3f) kJ/kg\n",Hmax*1e-3,udf.phiMax[1]*1e-3);
  if (udf.rifPremixed) {
    Message0("  Burned volume   = %15.3f mm^3\n", V_b*1e9);
    rk = pow(3.0/(4.0*3.14159265)*V_b,1.0/3.0);
    Message0("  Flamekernel-size= %15.3f mm\n", rk*1e3);
    z1Mean    /= totMass;
    egrMean   /= totMass;
    z1VarMean /= totMass;
    Message0("Mixture Summary:\n");
    Message0(" Zavg    = %15.6f\n", z1Mean);
    Message0(" EGRavg  = %15.6f\n", egrMean);
    Message0(" ZVaravg = %15.6f\n", z1VarMean);
  
    ct_x[0]  = z1Mean/(1-egrMean);
    ct_x[1]  = egrMean;
    ct_x[2]  = z1VarMean/(1-egrMean)/(1-egrMean);
    ct_x[3]  = 0;
    progMin  = MDUCchemTableProgMin(ct_x, udf.chemTabPrem);
    progMax  = MDUCchemTableProgMax(ct_x, udf.chemTabPrem);
    ct_x[3]  = (progMin+progMax)/2;
  }
  totMass = 0;
  debug = debug2 = 0;

  char **tabSp   = (char**) calloc(udf.nTabSp, sizeof(char*));
  if (udf.rifPremixed) {
    for (i=0; i<udf.nTabSp; i++)
      tabSp[i] = (char*) calloc(64, sizeof(char));
    MDUCgetChemTabSpecies(tabSp, udf.chemTabPrem);
  }

  counter = nZ1 = nStoich = nPremix = mass = totMass = 0;
  if (udf.rifPremixed)
  	Message0("Warning: Mixture fraction variance considered for table lookup (not for expansion though)!!!\n");

  for (iT=0; iT<udf.nThreadSolve; iT++) {    
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;
    
    begin_c_loop_all(c,t) {
      mass = C_R(c,t)*C_VOLUME(c,t); /* cell mass */
      if (udf.rifPremixed) {
        C_PROGSRC(c,t)   = 0.0;
        C_RIFSRC(c,t)    = 0.0;
        C_HR_PREM(c,t)   = 0.0;
         
        /* normalize the co-ordinate */
        /* CHANGE MARCO */
        /* Get max RIF progress variable to detect auto-ignition */
        for (i=0; i<nSpecies; i++) Ymean[i] = 0.0;
        if ( (RP_Get_Real("mduc/avg-z")   >= 0) &&
             (RP_Get_Real("mduc/avg-egr") >= 0) ) {
          EGRRIF  = RP_Get_Real("mduc/avg-egr");
          ZRIF    = RP_Get_Real("mduc/avg-z");
          ct_x[0] = ZRIF/(1-EGRRIF);
          ct_x[1] = EGRRIF; 
          ct_x[2] = C_Z1VAR(c,t)/(1-EGRRIF)/(1-EGRRIF);
          ct_x[3] = 0;
        }
        else {
          ct_x[0] = z1Mean/(1-egrMean);
          ct_x[1] = egrMean;
          ct_x[2] = z1VarMean/(1-ct_x[1])/(1-ct_x[1]);
          ct_x[3] = 0;
        } 

        MDUClookupChemTable("species", ct_x, Yun, udf.chemTabPrem);
        for (i=0; i<nSpecies; i++) YmeanU[i] = 0.0;
        
        for (i=0; i<udf.nTabSp; i++) { 
          if (udf.iTabSp[i] >= 0)  
            YmeanU[udf.iTabSp[i]] = Yun[i];
        }
        /* calculate minimum progress variable rif based on average mixture
           -> used for scaling */
        PROGRIFmin = 0.0; 
        for (i=0; i<udf.nProgSp; i++)
          PROGRIFmin += YmeanU[udf.iProgSp[i]]; 

        /* look up burned species */
        ct_x[3] = 1.0;
        MDUClookupChemTable("species", ct_x, Yun, udf.chemTabPrem);
        for (i=0; i<nSpecies; i++) YmeanU[i] = 0.0;
        
        for (i=0; i<udf.nTabSp; i++) { 
          if (udf.iTabSp[i] >= 0)  
            YmeanU[udf.iTabSp[i]] = Yun[i];
        }
        /* calculate maximum progress variable rif based on average mixture
           -> used for scaling*/
        PROGRIFmax = 0.0;
        for (i=0; i<udf.nProgSp; i++)
          PROGRIFmax += YmeanU[udf.iProgSp[i]]; 
      
        phiMean = (MIN_D(udf.phiMax[1],MAX_D(udf.phiMin[1],C_PHI1(c,t))) - udf.phiMin[1])/(udf.phiMax[1]-udf.phiMin[1]);
        phiVar  = C_PHIV1(c,t)/udf.deltaPhi2[1];
      }
      else {
      	phiMean = (C_PHI1(c,t) - udf.phiMin[1])/udf.deltaPhi[1];
      	phiVar  = C_PHIV1(c,t)/udf.deltaPhi2[1];
        phiMean = MAX_D(udf.phiMin[1],MIN_D(udf.phiMax[1],phiMean));
        phiVar  = MAX_D(0.0,phiVar);
      }

      if (phiMean < 0.0 || phiMean > 1.0) {
        Message0("Warning: mean is out of bounds\n");
        Message0("  phiMean = %8.6f --> phiVar = %15.6e (%8.6f/%8.6f)\n",
                 phiMean, phiVar, udf.phiMin[1], udf.phiMax[1]);
      }

      C_HR_RATE(c,t) = 0.0;
      for (i=0; i<nSpecies; i++) Ymean[i] = 0.0;
      /* ----- compute the mean mixture fraction of all species  -----  */
      if ( phiMean > udf.phiZero[1] || udf.rifEnthalpy) {
        if ( phiMean >= 0.045 && phiMean <= 0.065 ) nStoich++;

        /* compute mean mass fractions for this cell, also return PDF */
        pdfFlag = MDUCcalcMeanY_1D(udf.pdf1D,n1,x1,nSpecies,phiMean,phiVar,
                                   Y,Ymean,pdf,pdfBound);
       	C_HR_RATE(c,t) = MDUCcalcLocalMeanVar("heat-release",&phiMean,&phiVar,
                                              &pdf,&pdfBound,pdfFlag,udf.mduc);

        if (udf.rifPremixed) {
          /* TODO: maybe put this elsewhere */
          /* compute turbulent diffusivity scaling */
          lt = 0.37*pow(C_K(c,t),3.0/2.0)/C_D(c,t);
          C_DIFFSC(c,t) = MIN_D(pow(2*rk/lt,4.0/3.0),1.0);

          /* look up the maximum progress variable */
          ct_x[3] = 0;
	  C_PROGMAX(c,t) = MDUCchemTableProgMax(ct_x, udf.chemTabPrem);
          C_PROGMIN(c,t) = MDUCchemTableProgMin(ct_x, udf.chemTabPrem);

          /* Scaling coefficient */
          alpha = (PROGRIFmax-PROGRIFmin);

          /* Scale the RIF progress variable -> consistency */
          /* Compute RIF source term */
          C_CRIF(c,t) = 0.0;
          for (i=0; i<udf.nProgSp; i++)
            C_CRIF(c,t) += Ymean[udf.iProgSp[i]];
          C_CRIF(c,t) = MAX_D(0.0,MIN_D(1.0,(C_CRIF(c,t)-PROGRIFmin)/alpha));
          /* source term enforces C_RIFTOTAL reaches C_RIF 
 	     ensures non negative source term (RIF decreases close to walls)*/
          C_RIFSRC(c,t) = MAX_D(0.0,C_R(c,t)*(C_CRIF(c,t)-C_PROGTOTAL(c,t))/
                                CURRENT_TIMESTEP);

          if (RP_Get_Boolean("mduc/avg-mixture")) {
            if ( RP_Get_Real("mduc/avg-z")   >= 0 &&
                 RP_Get_Real("mduc/avg-egr") >= 0 ) {
              EGRRIF  = RP_Get_Real("mduc/avg-egr");
              ZRIF    =  RP_Get_Real("mduc/avg-z");
              ct_x[0] = ZRIF/(1-EGRRIF);
              ct_x[1] = EGRRIF; 
              ct_x[2] = C_Z1VAR(c,t)/(1-EGRRIF)/(1-EGRRIF);
              ct_x[3] = 0;
            }
            else {
              ct_x[0] = z1Mean/(1-egrMean);
              ct_x[1] = egrMean;
              ct_x[2] = z1VarMean/(1-ct_x[1])/(1-ct_x[1]);
              ct_x[3] = 0;
            } 
          }
          else {
            ct_x[0] = C_Z1MEAN(c,t)/(1-C_EGR(c,t))*factorZ;
            ct_x[1] = C_EGR(c,t);
            ct_x[2] = C_Z1VAR(c,t)/(1-ct_x[1])/(1-ct_x[1])*factorZ*factorZ;
            ct_x[3] = 0;
          }
          ct_x[3] = C_PROG(c,t);

	  /* clipping progress variable for table look up - don't need to
             do since should be included in table */
          if (ct_x[3] > C_PROGMAX(c,t)) ct_x[3] = C_PROGMAX(c,t);
          if (ct_x[3] < C_PROGMIN(c,t)) ct_x[3] = C_PROGMIN(c,t);

          
          /* compute the thresholds */
          /* heat release term should only be counted if not premixed */
          if (C_PROG(c,t) > 0.5)
            C_HR_RATE(c,t) = 0.0;
          
          /* store the flame speed */
          /* no temp or pressure expansion for laminar buring velocity yet */
          C_SL(c,t) = MDUClookupChemTable("sL", ct_x, NULL, udf.chemTabPrem); 
          
          if (udf.chemTabPdep && C_PROG(c,t)>1e-3+C_PROGMIN(c,t)) {
            /* compute the pressure ratio for the scaling */
            P         = ABS_P(C_P(c,t),op_pres);
            /* MARCO - limiting pressure expansion to 30 bar */
            P         = MAX_D(5.0e5,MIN_D(P,30.0e5));
            /* H         = C_HTOT(c,t);*/
            p_tab_ref = MDUCgetChemTabPref(udf.chemTabPres);
            p_ratio   = log(P/p_tab_ref);
     	    
            ct_Unburned[0] = ct_x[0];
            ct_Unburned[1] = ct_x[1];
            ct_Unburned[2] = ct_x[2];
            ct_Unburned[3] = 0;
            
            h_tab_ref = MDUCgetChemTabTref(udf.chemTabPres);
            /* get unburned species */
            MDUClookupChemTable("species", ct_Unburned, Yun, udf.chemTabPrem);
            for (i=0; i<nSpecies; i++) YmeanU[i] = 0.0;
            
            for (i=0; i<udf.nTabSp; i++) { 
              if (udf.iTabSp[i] >= 0)  
                YmeanU[udf.iTabSp[i]] = Yun[i];
            }
            /* Calc unburned Temperature */
            H = MDUCcalcTfromH(nSpecies,C_HTOT(c,t),YmeanU,1000.0,udf.mduc);
            /* MARCO limiting temperature expansion to 1200K */
            H = MAX_D(MIN_D(H,1200.0),450.0);
            h_ratio   = H/h_tab_ref;
            
            ct_noC[0] = ct_x[0];
            ct_noC[1] = ct_x[1];
            /* no variance for coefficients */
            ct_noC[2] = ct_x[2]*0; 

            /* SOURCE TERM */
            
            fs[0] = MDUClookupChemTable("SRC_PROG_Par1",ct_noC,NULL,udf.chemTabPres);
            fs[1] = MDUClookupChemTable("SRC_PROG_Par2",ct_noC,NULL,udf.chemTabPres);
            fs[2] = MDUClookupChemTable("SRC_PROG_Par3",ct_noC,NULL,udf.chemTabPres);
            fs[3] = MDUClookupChemTable("SRC_PROG_Par4",ct_noC,NULL,udf.chemTabPres);
            fs[4] = MDUClookupChemTable("SRC_PROG_Par5",ct_noC,NULL,udf.chemTabPres);
            fs[5] = MDUClookupChemTable("SRC_PROG_Par6",ct_noC,NULL,udf.chemTabPres);
            fs[6] = MDUClookupChemTable("SRC_PROG_Par7",ct_noC,NULL,udf.chemTabPres);
            fs[7] = MDUClookupChemTable("SRC_PROG_Par8",ct_noC,NULL,udf.chemTabPres);
            fs[8] = MDUClookupChemTable("SRC_PROG_Par9",ct_noC,NULL,udf.chemTabPres);
            fs[9] = MDUClookupChemTable("SRC_PROG_Par10",ct_noC,NULL,udf.chemTabPres);
            fs[10] = MDUClookupChemTable("SRC_PROG_Par11",ct_noC,NULL,udf.chemTabPres);
            
            facP    = exp(fs[0]*p_ratio*p_ratio + fs[1]*p_ratio);
            facCM   =    MAX_D(fs[3]*(h_ratio-1.0)+fs[4]*h_ratio*h_ratio+1.0-fs[4],1e-10);
            facCa   =    MAX_D(fs[5]*(h_ratio-1.0)+1.0,1e-10);
            facCb   =    MAX_D(fs[6]*(h_ratio-1.0)+1.0,1e-10);
            
            facH    = exp(fs[8]*log(h_ratio)*log(h_ratio) + fs[7]*log(h_ratio));
            cMaxRef = fs[9];
            minRef  = fs[10];
            /* transform the progress variable to the reference profile */
            if (ct_x[3]-cMaxRef*facCM <  0.0)   
              cRef = (ct_x[3] - cMaxRef*facCM) / facCa + cMaxRef;
            else 
              cRef = (ct_x[3] - cMaxRef*facCM) / facCb + cMaxRef;   
            /* look up the reference profile */
            for (j=0;j<3;j++) ct_xRef[j] = ct_x[j];
            ct_xRef[3] = cRef;
            /* scale the source term */
            if (C_PROG(c,t)<C_PROGMAX(c,t) && C_PROG(c,t)>1e-3+C_PROGMIN(c,t)) {
              progSRC = MDUClookupChemTable("SRC_PROG", ct_xRef, NULL,udf.chemTabPrem);
              C_PROGSRC(c,t) = MAX_D(progSRC,0.0)*facP*facH; 
              
            }
            else {
              C_PROGSRC(c,t) = 0;
            }
            /* check if already auto-ignited. If so, set source term to 0 to conserve C_FPVA -> regime tracking */
            if (C_PROG(c,t)<C_PROGTOTAL(c,t)) C_PROGSRC(c,t) = 0;
            /* Marking regime */
            if (C_PROGSRC(c,t)> C_RIFSRC(c,t)+1.0) {
              C_REGIME(c,t) = 1.0; /* flame */
            }
            else {
              if (C_PROG(c,t)>C_PROGTOTAL(c,t)) C_REGIME(c,t) = 0.0;
              else C_REGIME(c,t) = 2.0;
            }
            /* check if cell is in flame front using both source terms. If yes,
               read in species mass fraction with C_PROG (only source term
               from FPVA). Otherwise use C_PROGRIF  */
            /* this is done to ensure right flame speed during flame propagation
               and reading fully burned products after auto-ignition (C_PROG
               mixes after auto-ignition --> taking this progress variable
               would */
            /* lead to species from the reaction zone */
            /* if (C_PROGSRC(c,t)<C_RIFSRC(c,t)+1e-3) ct_x[3]=MAX_D(C_PROG(c,t),C_CRIF(c,t));*/
            ct_x[3] = MAX_D(C_PROG(c,t),C_PROGTOTAL(c,t));
            /* SPECIES */
            MDUClookupChemTable("species", ct_noC, fp, udf.chemTabPres);
     	    
            for (i=0; i<udf.nTabSp; i++) {
              facP    = exp(fp[i*npar+0]*p_ratio*p_ratio + fp[i*npar+1]*p_ratio);
              facCM   =  MAX_D(fp[i*npar+3]*(h_ratio-1.0)+fp[i*npar+4]*h_ratio*h_ratio+1.0-fp[i*npar+4],1e-10);
              facCa   =    MAX_D(fp[i*npar+5]*(h_ratio-1.0)+1.0,1e-10);
              facCb   =    MAX_D(fp[i*npar+6]*(h_ratio-1.0)+1.0,1e-10);
              facH    = exp(fp[i*npar+8]*log(h_ratio)*log(h_ratio) + fp[i*npar+7]*log(h_ratio));
              cMaxRef =     fp[i*npar+9];
              minRef  =     fp[i*npar+10];
     	      
              /* transform the progress variable to the reference profile */
              if (ct_x[3]-cMaxRef*facCM <  0.0 )   
                cRef = (ct_x[3] - cMaxRef*facCM) / facCa + cMaxRef;
              else 
                cRef = (ct_x[3] - cMaxRef*facCM) / facCb + cMaxRef;   
              /* look up the reference profile */
              for (j=0;j<3;j++) ct_xRef[j] = ct_x[j];
              ct_xRef[3] = cRef;
              
              strcpy(tabSpI, "Y_");
              strcat(tabSpI,tabSp[i]);
     	      
              Ytab[i] = MDUClookupChemTable(tabSpI, ct_xRef, NULL, udf.chemTabPrem); 
              
              /* scale the profile */
              Ytab[i] = facP*((Ytab[i]-minRef)*facH+minRef);    
              if ((Ytab[i] > 1.0) || (Ytab[i] < 0.0))       {
                sprintf(errFileName, "%s/mduc/Ymean_err_node%i_cell%i-%05i",
                        udf.caseDir, myid, c, N_TIME);
                Message0("Problem with species %i:\n", i);
                Message0("Massfraction : %1.3f\n", Ytab[i]);
                Message0("Warning: writing Ymean_err: %s\n", errFileName);
                errFile = fopen(errFileName, "w");
                fprintf(errFile, "Problem with species: %s\n",tabSp[i]);
                fprintf(errFile, "Ytab     = %8.6f\n",  MDUClookupChemTable(tabSpI, ct_xRef, NULL, udf.chemTabPrem));
                fprintf(errFile, "phiMean  = %15.6e\n", phiMean);
                fprintf(errFile, "phiVar   = %15.6e\n", phiVar);
                fprintf(errFile, "phiMin   = %15.6e\n", udf.phiMin[1]);
                fprintf(errFile, "phiMax   = %15.6e\n", udf.phiMax[1]);
                fprintf(errFile, "deltaPhi = %15.6e\n", udf.deltaPhi[1]);
                if (udf.rifPremixed) {
                  fprintf(errFile, "Z1       = %15.6e\n", C_Z1MEAN(c,t));
                  fprintf(errFile, "EGR      = %15.6e\n", C_EGR(c,t));
                  fprintf(errFile, "Z1-var   = %15.6e\n", C_Z1VAR(c,t));
                  fprintf(errFile, "prog     = %15.6e\n", C_PROG(c,t));
                  fprintf(errFile, "prog-src = %15.6e\n", C_PROGSRC(c,t));
                  fprintf(errFile, "HTOT     = %15.6e\n", H);
                  fprintf(errFile, "Href     = %15.6e\n", h_tab_ref);
                  fprintf(errFile, "pratio   = %15.6e\n", p_ratio);
                  fprintf(errFile, "hratio   = %15.6e\n", h_ratio);
                  fprintf(errFile, "facP     = %15.6e\n", facP);
                  fprintf(errFile, "facH     = %15.6e\n", facH);
                  fprintf(errFile, "cRef     = %15.6e\n", cRef);
                  fprintf(errFile, "facCa    = %15.6e\n", facCa);
                  fprintf(errFile, "facCb    = %15.6e\n", facCb);
                  for (j=0;j<npar;j++) 
                    fprintf(errFile, "fp[%i]   = %15.6e\n", j, fp[i*npar+j]);
                }
                fclose(errFile);
              }
            }
     	    
          }
          else {
            C_PROGSRC(c,t) = MDUClookupChemTable("SRC_PROG",ct_x,NULL,
                                                 udf.chemTabPrem);
            
            if (C_PROG(c,t) < C_PROGTOTAL(c,t)) C_PROGSRC(c,t) = 0;
            /*if (C_PROGSRC(c,t)<C_RIFSRC(c,t)+1e-3) ct_x[3]=MAX_D(C_PROG(c,t),C_CRIF(c,t));*/
            /* Marking regime */
            if (C_PROGSRC(c,t)> C_RIFSRC(c,t)+1.0) {
              C_REGIME(c,t) = 1.0; /* flame */
            }
            else {
              if (C_PROG(c,t)>C_PROGTOTAL(c,t)) C_REGIME(c,t) = 0.0;
              else C_REGIME(c,t) = 2.0;
            }
            
            ct_x[3] = MAX_D(C_PROG(c,t),C_PROGTOTAL(c,t));
            for (i=0; i<udf.nTabSp; i++) { 
              strcpy(tabSpI, "Y_"); 
              strcat(tabSpI,tabSp[i]); 
              Ytab[i] = MDUClookupChemTable(tabSpI,ct_x,NULL,udf.chemTabPrem); 
            }  
          }
          /* Scaled progress variable from FPVA model as "regime fraction"*/
          /* C_PROG = 1 is already completely burned and has not mixed yet */
          /* C_PROG = 0 was never touched by premixed flame */
          /* rest is in the flame front or mixed after auto-igntion */
          
          /* IN FRONT AND BEHIND FLAME FRONT */
          /* usage of "regime-fraction" */
          /* --> unburned region --> RIF */
          /* burned region depends on which regime was present during combution
             and mixing after combustion */
          if (C_PROGSRC(c,t) < C_RIFSRC(c,t)+1e-2) {
            for (i=0; i<nSpecies; i++) {
              if (indexTab[i] >= 0)
                Ymean[i] = Ytab[indexTab[i]]*MIN_D(1.0,MAX_D(0.0,C_PROG(c,t)))+Ymean[i]*(1.0-MIN_D(1.0,MAX_D(0.0,C_PROG(c,t))));
              else
                Ymean[i] = Ymean[i]*(1.0-MIN_D(1.0,MAX_D(0.0,C_PROG(c,t))));
            }
          }
          
          /* IN THE FLAME FRONT */
          /* taking fraction based species impact the flamespeed!!!
             (wrong species --> wrong T and W --> wrong density -->
             d(rho*C)/d(t)=C_SRC */
          /* --> only tabled species considered here */
          else {
            for (i=0; i<nSpecies; i++) Ymean[i]=0.0;
            for (i=0; i<nSpecies; i++) {
              if (indexTab[i] >= 0)
                Ymean[i] = Ytab[indexTab[i]];
            }
          }
          /* for (i=0; i<nSpecies; i++) Ymean[i]=0.0;
             for (i=0; i<nSpecies; i++) {
             if (indexTab[i] >= 0)
             Ymean[i] = Ytab[indexTab[i]];*/
          /* Clipping progress Variable */
          if (C_PROG(c,t) > C_PROGMAX(c,t)+1e-3)
            C_PROG(c,t)=C_PROGMAX(c,t)+1e-3;
          if (C_PROGTOTAL(c,t) > C_PROGMAX(c,t)+1e-3) {
            C_PROGTOTAL(c,t)=C_PROGMAX(c,t)+1e-3;
            /* dump warning only once */
            if (debug2 == 0 || debug2 == 2) {
              Message0("Exceeding maximum progress variable --> Clipping\n");
              debug2 += 1;
            }
          }    
          if (C_PROG(c,t)>0.9*C_PROGMAX(c,t)) nPremix++;
          
          if (C_PROG(c,t)<0.5) {
            for (j=pdfBound[0]; j<=pdfBound[1]; j++)
              pdfGlobal[j] += pdf[j]*mass;
          }
        } /* end rifPremixed */
        else {
          /* store local contribution to global PDF (mass averaged) */
          for (j=pdfBound[0]; j<=pdfBound[1]; j++)
            pdfGlobal[j] += pdf[j]*mass;
        }
        
        nZ1++; /* count the number of cells computed */
      } /* end cell computation for non-zero Z1 */
      else {
        /* no mixture fraction, update based on the oxidizer species */
        /* essentially delta function at pure oxidizer */
        for (i=0; i < nSpecies; i++)
          Ymean[i] = Y[0][i];

        /* store the heat release rate */
        phiMean     = 0.0;
        phiVar      = 0.0;
        pdf[0]      = 1.0;
        pdf[1]      = 0.0;
        pdfBound[0] = 0;
        pdfBound[1] = 1;
        C_HR_RATE(c,t) = MDUCcalcLocalMeanVar("heat-release",&phiMean,&phiVar,
                                              &pdf,&pdfBound,-1,udf.mduc);

        pdfGlobal[0] += mass;
      }
      
      /* ----- update the cell properties ----- */
      /* update the mixture fractions tracked by Fluent */
      sumYTab=0;
      for (i=0;i<udf.nTabSp; i++) {
        sumYTab += Ytab[i];
      }
      
      sumY = 0;
      for (i=0; i<nYFav; i++) {
        if (indexYFav[i] >= 0 ) { 
          C_YI(c,t,i)= Ymean[indexYFav[i]];
          sumY += Ymean[indexYFav[i]];
        }
      }

      /* update the cell temperature by iterating from the total
         enthalpy and mass fractions */
      if (((sumY>1.05)||(sumY<0.95))&&debug==0) {      
        sprintf(errFileName, "%s/mduc/Ymean_err_cell%i-%05i",
                udf.caseDir, c, N_TIME);
        Message0("Sum of mass fraction not equal to 1!!!");
        Message0("Sum(Y)           : %3.3f\n", sumY);
        Message0("Z1       = %3.6f\n", C_Z1MEAN(c,t));
        Message0("HTOT             : %1.3f\n", H);
        if (udf.rifPremixed) {
          Message0("H_tab_ref        : %1.3f\n", h_tab_ref);
          Message0("h_ratio          : %1.3f\n", h_ratio);
          Message0("P                : %1.3f\n", P);
          Message0("P_tab_ref        : %1.3f\n", p_tab_ref);
          Message0("p_ratio          : %1.3f\n", p_ratio);
          Message0("PROG     = %3.6f\n", C_PROG(c,t));
          Message0("PROG-SRC = %3.6f\n", C_PROGSRC(c,t));
          Message0("FactorP          : %1.3f\n", facP);
          Message0("FactorH          : %1.3f\n", facH);
          Message0("Tabulated species   : %i\n", udf.nTabSp);
        }
        
        for (i=0; i<nYFav; i++)
          Message0("%7s = %8.6f\n", udf.species[i], C_YI(c,t,i));
        Message0("Warning: writing : %s\n", errFileName);
        errFile = fopen(errFileName, "w");
        if (errFile == NULL) {
          Message0("File not created");
        }
        fprintf(errFile, "Sum(YMean)= %8.6f\n", sumY);
        fprintf(errFile, "phiMean  = %15.6e\n", phiMean);
        fprintf(errFile, "phiVar   = %15.6e\n", phiVar);
        fprintf(errFile, "phiMin   = %15.6e\n", udf.phiMin[1]);
        fprintf(errFile, "phiMax   = %15.6e\n", udf.phiMax[1]);
        fprintf(errFile, "deltaPhi = %15.6e\n", udf.deltaPhi[1]);
        fprintf(errFile, "Z1       = %15.6e\n", C_Z1MEAN(c,t));
        /* fprintf(errFile, "Z2       = %15.6e\n", C_Z2MEAN(c,t));*/
        fprintf(errFile, "Z1-var   = %15.6e\n", C_Z1VAR(c,t));
        if (udf.rifPremixed) {
          fprintf(errFile, "prog     = %15.6e\n", C_PROG(c,t));
          fprintf(errFile, "prog-src = %15.6e\n", C_PROGSRC(c,t));
        }
        fclose(errFile);
        debug = 1;
      }

      T_start  = C_T(c,t);
      C_T(c,t) = MDUCcalcTfromH(nSpecies,C_HTOT(c,t),Ymean,T_start,udf.mduc);
      if (C_T(c,t)<100 || C_T(c,t)>5000) {
        Message("Warning: Temperature iteration failed ->"
                "Clipping temperature from %5.1f to %5.1f\n",C_T(c,t),
                (T_start+C_T(c,t))/2);
        C_T(c,t) = (T_start+C_T(c,t))/2;
      }

      /* Error checking on Fluent temperature, reasonably bounded */
      if ( (C_T(c,t) > 5000) || (C_T(c,t) < 100) || isnan(C_T(c,t)) ) {
        double cX[ND_ND];
        C_CENTROID(cX,c,t);
        sprintf(errFileName, "%s/mduc/T_err_node%i_cell%i_%05i", udf.caseDir,
                myid, c, N_TIME);
        Message("Warning: writing T_err file: %s\n", errFileName);
        errFile = fopen(errFileName, "w");
        fprintf(errFile, "# zone: %s\n", THREAD_NAME(t));
        fprintf(errFile, "# co-ords:    (%15.6e %15.6e %15.6e)\n",
                cX[0], cX[1], cX[2]);
        fprintf(errFile, "# phiMean     = %15.6e\n", phiMean);
        fprintf(errFile, "# phiVar      = %15.6e\n", phiVar);
        fprintf(errFile, "# phiMin      = %15.6e\n", udf.phiMin[1]);
        fprintf(errFile, "# phiMax      = %15.6e\n", udf.phiMax[1]);
        fprintf(errFile, "# deltaPhi    = %15.6e\n", udf.deltaPhi[1]);
        fprintf(errFile, "# C_T(c,t)    = %15.6e K\n", C_T(c,t));
        fprintf(errFile, "# C_HTOT(c,t) = %15.6e J/kg\n", C_HTOT(c,t));
        if (udf.rifPremixed) {
          fprintf(errFile, "# Z1          = %15.6e\n", C_Z1MEAN(c,t));
          fprintf(errFile, "# Z1-var      = %15.6e\n", C_Z1VAR(c,t));
          fprintf(errFile, "# regime      = %15.6e\n", C_REGIME(c,t));
          fprintf(errFile, "# prog        = %15.6e\n", C_PROG(c,t));
          fprintf(errFile, "# prog-max    = %15.6e\n", progMax);
          fprintf(errFile, "# prog-rif    = %15.6e\n", C_CRIF(c,t));
          fprintf(errFile, "# prog-src    = %15.6e\n", C_PROGSRC(c,t));
        }
        fflush(errFile);
        fprintf(errFile, "# index/x1 Yi/pdf\n");
        for (i=0; i<nSpecies; i++)
          fprintf(errFile, "%8i %15.6e\n",i,Ymean[i]);

        fprintf(errFile, "\n\n");
        for (j=0; j<n1; j++)
          fprintf(errFile, "%8.6f %15.6e\n", x1[j], pdf[j]);
        fclose(errFile);

#if RP_NODE
        PRF_GSYNC();
#endif /* RP_NODE */
        Error("Temperature out of bounds!\n");
      }
      
      /* update counter and print progress to screen */
      counter++;
      if (counter%nCount == 0)             Message0(".");
      if (counter%(lineCount*nCount) == 0) Message0("\n");

#if DEBUG_SDM
      taui = C_D(c,t)/C_K(c,t); /* inverse turb. timescale */
      
      GradZ[0] = (C_Z1MEAN_G(c,t)[0] + C_Z2MEAN_G(c,t)[0]);
      GradZ[1] = (C_Z1MEAN_G(c,t)[1] + C_Z2MEAN_G(c,t)[1]);
      GradZ[2] = (C_Z1MEAN_G(c,t)[2] + C_Z2MEAN_G(c,t)[2]);
      
      C_Z0VAR_PRD(c,t)   = 2.0*C_MU_T(c,t)/0.7*NV_MAG2(GradZ);
      C_Z1VAR_PRD(c,t)   = 2.0*C_MU_T(c,t)/0.7*NV_MAG2(C_Z1MEAN_G(c,t));
      C_Z2VAR_PRD(c,t)   = 2.0*C_MU_T(c,t)/0.7*NV_MAG2(C_Z2MEAN_G(c,t));
      C_Z0VAR_D2(c,t)    = 2.0*taui*C_R(c,t)*C_Z0VAR(c,t);
      C_Z1VAR_D2(c,t)    = 2.0*taui*C_R(c,t)*Z1var;
      C_Z2VAR_D2(c,t)    = 2.0*taui*C_R(c,t)*C_Z2VAR(c,t);
      C_Z0VAR_SRC_R(c,t) = C_R(c,t)*C_Z0VAR_SRC(c,t);
      C_Z1VAR_SRC_R(c,t) = C_R(c,t)*C_Z1VAR_SRC(c,t);
      C_Z2VAR_SRC_R(c,t) = C_R(c,t)*C_Z2VAR_SRC(c,t);
      C_Z0VAR_EPS_R(c,t) = C_R(c,t)*C_Z0VAR_EPS(c,t)*C_Z0VAR(c,t);
      C_Z1VAR_EPS_R(c,t) = C_R(c,t)*C_Z1VAR_EPS(c,t)*Z1var;
      C_Z2VAR_EPS_R(c,t) = C_R(c,t)*C_Z2VAR_EPS(c,t)*C_Z2VAR(c,t);
      
      C_Z0VAR_SRC_L(c,t) *= C_R(c,t)*C_Z0VAR(c,t);
      C_Z1VAR_SRC_L(c,t) *= C_R(c,t)*Z1var;
      C_Z2VAR_SRC_L(c,t) *= C_R(c,t)*C_Z2VAR(c,t);
#endif /* DEBUG_SDM */
      
      /* sum total mass */
      totMass += mass;
    } end_c_loop_all(c,t);
  } /* end loop over cell threads */

#if RP_NODE
 /* if parallel, synchronize all node solutions here by summing the
    contributions from each node */
  PRF_GRSUM(pdfGlobal, n1, rwork);
  totMass = PRF_GRSUM1(totMass); /* total mass in the cylinder */
  nZ1     = PRF_GISUM1(nZ1);     /* total cells updated on all nodes */
  nStoich = PRF_GISUM1(nStoich); /* total cells near stoichiometric */
  nPremix = PRF_GISUM1(nPremix); /* total cells near stoichiometric */
  counter = PRF_GISUM1(counter); /* total cells */
#endif /* RP_NODE */

  /* normalise the global pdf */
  for (j=0; j<n1; j++)
    pdfGlobal[j] /= totMass;

  tmp = 100.0/counter;
  Message0("\n Update summary:\n");
  Message0(" ----- cells updated -----  number ----- percent -----\n");
  Message0("\t Z1        \t %8i \t %6.2f\n", nZ1, (double)(nZ1)*tmp);
  Message0("\t stoich.   \t %8i \t %6.2f\n", nStoich, (double)(nStoich)*tmp);
  Message0("\t premixed  \t %8i \t %6.2f\n", nPremix, (double)(nPremix)*tmp);
  Message0(" -----------------------------------------------------\n");

  free(pdf); free(pdfBound);
  free(Ymean);free(YmeanU); free(Ytab);free(Yun); free(rwork);
  for (i=0; i<udf.nTabSp; i++) free(tabSp[i]);
  free(tabSp);
  free(fp);
  free(fs);
  free(indexTab);
}
/* ------------------------------------------------------------------------- */

/* ===== start 2D functions ===== */
void UDF_rif_2D(Domain *d, int inert, double PAvg, void *mduc)
{
  Message0("\n ----- main injection -----\n");
  /* ----- declare local variables ----- */
  int i, j, k;
  /* Fluent variables */

  /* flamelet variables */
  double **chi1, **chi12, **chi2;    /* 2D scalar dissipation rate */
  int    **pdfBoundG; double **pdf;  /* global pdf returned from cell update */
  double ***Y, **Temp;               /* mass fractions and temperature */
  double Tox, Tfuel1, Tfuel2;

  /* ----- get the flamlet variables ---- */
  int     nSpecies;
  int    *nGridMDUC, n1, n2;
  double **MDUCgrid, *x1, *x2;

  /* initialise the grid/species information of the flamelet */
  nSpecies  = MDUCgetNspecies(mduc);
  nGridMDUC = MDUCgetGrid2D(&MDUCgrid, mduc);
  n1 = nGridMDUC[0]; n2 = nGridMDUC[1];
  x1 = MDUCgrid[0];  x2 = MDUCgrid[1];

  /* allocate local memory */
  Y         = (double***) calloc(n2, sizeof(double**));
  Temp      = (double**)  calloc(n2, sizeof(double*));
  chi1      = (double**)  calloc(n2, sizeof(double*));
  chi12     = (double**)  calloc(n2, sizeof(double*));
  chi2      = (double**)  calloc(n2, sizeof(double*));
  pdf       = (double**)  calloc(n2, sizeof(double*));
  pdfBoundG = (int**)     calloc(n2, sizeof(int*));

  for (k=0; k<n2; k++) {
    Y[k]         = (double**) calloc(n1, sizeof(double*));
    Temp[k]      = (double*)  calloc(n1, sizeof(double));
    chi1[k]      = (double*)  calloc(n1, sizeof(double));
    chi12[k]     = (double*)  calloc(n1, sizeof(double));
    chi2[k]      = (double*)  calloc(n1, sizeof(double));
    pdf[k]       = (double*)  calloc(n1, sizeof(double));
    pdfBoundG[k] = (int*)     calloc(2,  sizeof(int));
    pdfBoundG[k][0] = 0;
    pdfBoundG[k][1] = n1-1;
    
    for (j=0; j<n1; j++)
      Y[k][j] = (double*) calloc(nSpecies, sizeof(double));
  }

  /* ----- compute 2D scalar dissipation rate ----- */
  UDFtimerStart("chi", udf.timer);
  if (udf.chiMode[2] == 1) /* compute a distribution */
    UDFcalcChi2D(d, nGridMDUC, MDUCgrid, chi1, chi12, chi2);
  else if (udf.chiMode[2] == 2) /* simple binning */
    UDFbinChi2D(d,n1,n2,x1,x2,chi1,chi2);
  else
    Message0("Warning: chi2-mode not defined!\n");
  UDFtimerStop("chi", udf.timer);

  /* ----- call the 2D flamelet solver ----- */
  if (inert) {
    /* compute average temperatures */
    Tox    = calcAvgTemp(d);
    Tfuel1 = MDUCgetPointT(1.0, 0.0, mduc);
    Tfuel2 = MDUCgetPointT(0.0, 1.0, mduc);
#if RP_NODE
    if (myid == node_zero)
#endif /* RP_NODE */
      MDUCupdateInert2D(CURRENT_TIMESTEP,PAvg,udf.phiMin[1],udf.phiMax[1],
                        udf.phiMin[2],udf.phiMax[2],Tox,Tfuel1,Tfuel2,
                        chi1,chi12,chi2,Y, Temp, mduc);
  }
  else {
    UDFtimerStart("flamelet", udf.timer);
    MDUCsolve2D(CURRENT_TIMESTEP, PAvg,udf.phiMin[1],udf.phiMax[1],
                udf.phiMin[2],udf.phiMax[2],chi1,chi12,chi2,Y,Temp,mduc);
  }

#if RP_NODE
  PRF_GSYNC(); /* make sure all nodes are finished before proceeding */

  if (myid == node_zero) {
    /* send the solution to all other nodes for updating */
    for (i=1; i<compute_node_count; i++) {
      for (k=0; k<n2; k++) {
        for (j=0; j<n1; j++)
          PRF_CSEND_REAL(i, Y[k][j], nSpecies, myid);
        PRF_CSEND_REAL(i, Temp[k], n1, myid);
      } /* end loop over rows */
    } /* end loop over compute nodes */
  }
  else { /* receive solution on nodes, needed for cell update */
    for (k=0; k<n2; k++) {
      for (j=0; j<n1; j++)
        PRF_CRECV_REAL(node_zero, Y[k][j], nSpecies, node_zero);
      PRF_CRECV_REAL(node_zero, Temp[k], n1, node_zero);
    } /* end loop over rows */
  }
#endif /* RP_NODE */
  UDFtimerStop("flamelet", udf.timer);

  /* ----- update the CFD solution ----- */
  UDFtimerStart("update", udf.timer);
  if (!inert)
    updateCFD_2D(d, udf.nY, udf.indexY, nGridMDUC, MDUCgrid, Y, Temp, pdf);
  UDFtimerStop("update", udf.timer);

  /* ----- write the flow data to a file ----- */
#if RP_NODE
  if (myid == node_zero)
#endif /* RP_NODE */
    {
      MDUCwriteFlowData(chi1, chi2, pdf, NULL, mduc);
      MDUCwriteScalars(CURRENT_TIME,NULL,NULL,pdf,pdfBoundG,-1,mduc);
    }

  for (k = 0; k<n2; k++) {
    for (j = 0; j<n1; j++) free(Y[k][j]);
    free(Y[k]); free(Temp[k]);
    free(chi1[k]); free(chi12[k]); free(chi2[k]);
    free(pdf[k]); free(pdfBoundG[k]);
  }
  free(Y); free(Temp); free(chi1); free(chi12); free(chi2);
  free(nGridMDUC); free(MDUCgrid[0]); free(MDUCgrid[1]); free(MDUCgrid);
  free(pdf); free(pdfBoundG);
  Message0("\n ----- exiting main injection -----\n");
}
/* -------------------------------------------------------------------------- */

void updateCFD_2D(Domain *d, int nYFav, int *indexYFav, int *n, double **x,
                  double ***Y, double **Temp, double **pdfGlobal)
{
  Message0("\n ----- Averaging flamelet solution and updating cells -----\n");
  Thread *t;
  cell_t c;

  int    i, j, k, iT;
  int    nZ1, nZ2, nZ, counter, nCount = 500, lineCount = 200;
  int    nZ1L, nZ2L;
  int    nStoich, nDiric;
  char   pdf2Mode[50];
  int    pdfTreatment = 0;
  double phiMean[3], phiVar[3];
  double *Ymean, **pdf, *pdfX2, *pdfX1, **Y2, *rwork;
  int    **pdfBound, pdfBound1D[2], pdfBoundX1[2], pdfBoundX2[2], pdfFlag;
  double mass, totMass, pdfTmp, tmp;
  double mean2Loc, mean1Loc, meantot;
  int    nNotCvg, nNotCvgL, nNegRho, nRhoP;
  double **meanPDF, **varPDF, **meanC, **varC;
  double **X;
  double covar, rhoCoeff;
  int ii, row;
  double vFac = 1.0; /* factor to scale the variance by */
  vFac = udf.varScale;

  Message0("variance scale: %g\n", vFac);

  /* flamelet grid parameters */
  int    nSpecies = MDUCgetNspecies(udf.mduc);
  int    n1       = n[0];
  int    n2       = n[1];
  double *x1      = x[0];
  double *x2      = x[1];
  double dx1_0    = x1[1] - x1[0];
  double dx2_0    = x2[1] - x2[0];

  /* zero values */
  double TINY    = 1e-4;
  double TINYVAR = 1e-4;

  /* debugging */
  double GradZ[3], taui;

  /* debugging variables */
  FILE *pFile;
  char fileName[100], prefix[25], header[200], buffer[200];
  int yOutofBounds, Ysum;

  /* allocate all local memory */
  Ymean    = (double*)  calloc(nSpecies, sizeof(double));
  Y2       = (double**) calloc(n2, sizeof(double*));
  pdf      = (double**) calloc(n2, sizeof(double*));
  pdfX1    = (double*)  calloc(n1, sizeof(double));
  pdfX2    = (double*)  calloc(n2, sizeof(double));
  pdfBound = (int**)    calloc(n2, sizeof(int*));
  rwork    = (double*)  calloc(n1, sizeof(double));
  meanPDF  = (double**) calloc(1, sizeof(double*));
  varPDF   = (double**) calloc(1, sizeof(double*));
  meanC    = (double**) calloc(1, sizeof(double*));
  varC     = (double**) calloc(1, sizeof(double*));

  for (k=0; k < n2; k++) {
    Y2[k]       = (double*) calloc(nSpecies, sizeof(double));
    pdf[k]      = (double*) calloc(n1, sizeof(double));
    pdfBound[k] = (int*) calloc(2, sizeof(int));
  }

  /* construct the grid for computing moments */
  X = (double**) calloc(3, sizeof(double*));
  for (i=0; i<3; i++) X[i] = (double*) calloc(n2*n1, sizeof(double));
  for (k=0; k<n2; k++) {
    row = k*n1;
    for (j=0; j<n1; j++) {
      if (udf.transform == 0) {
        X[0][row+j] = x1[j] + x2[k];
        X[1][row+j] = x1[j];
        X[2][row+j] = x2[k];
      }
      else if (udf.transform == 1) {
        X[0][row+j] = x1[j] + x2[k]*(1.0-x1[j]);
        X[1][row+j] = x1[j];
        X[2][row+j] = x2[k]*(1.0-x1[j]);
      }
      else if (udf.transform == 2) {
        X[0][row+j] = x1[j];
        X[1][row+j] = x1[j]*(1 - x2[k]);
        X[2][row+j] = x2[k]*x1[j];
      }
    } /* end loop over columns */
  } /* end loop over rows */

  /* store the first column of the solution for use when computing a 1D pdf
     if Z1 = 0 */
  for (k=0; k<n2; k++)
    for (i=0; i<nSpecies; i++)
      Y2[k][i] = Y[k][0][i];

  /* make sure the global PDF is set to zero */
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      pdfGlobal[k][j] = 0.0;

  totMass = 0.0;
  counter = nZ = nZ1 = nZ2 = nDiric = nNotCvg = nNegRho = nRhoP = nStoich = 0;
  nZ1L = nZ2L = 0;
  nNotCvgL = 0;
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL)
      continue;
    
    begin_c_loop_all(c,t) {
      yOutofBounds = 0;
      mass     = C_R(c,t)*C_VOLUME(c,t); /* mass of this cell */
      mean1Loc = C_PHI1(c,t);
      mean2Loc = C_PHI2(c,t);
      meantot = mean1Loc + mean2Loc;
      if ( (meantot >= 0.045) && (meantot <= 0.065) )
        nStoich++;

      /* normalize the co-ordinate */
      phiMean[1] = (C_PHI1(c,t) - udf.phiMin[1])/udf.deltaPhi[1];
      phiVar[1]  = C_PHIV1(c,t)/udf.deltaPhi2[1];
      phiMean[2] = (C_PHI2(c,t) - udf.phiMin[2])/udf.deltaPhi[2];
      phiVar[2]  = C_PHIV2(c,t)/udf.deltaPhi2[2];

      for (i=0; i<nSpecies; i++) Ymean[i] = 0.0;
      /* ----- compute the mean mixture fraction of all species ----- */
      if ( phiMean[1] > udf.phiZero[1] ) {
        /* 2D distribution */
        if ( phiMean[2] > udf.phiZero[2] ) {
          /* set the default joint pdf treatment */
          strcpy(pdf2Mode, udf.pdf2D);
          
          covar    = 0.5*(C_PHIV0(c,t) - C_PHIV1(c,t) - C_PHIV2(c,t));
          rhoCoeff = covar/sqrt(C_PHIV1(c,t)*C_PHIV2(c,t));
          if (rhoCoeff < 0.0)
            nNegRho++;
          if (fabs(1.0-rhoCoeff) < 0.02)
            nRhoP++;
          
          /* for means less than the first grid size, use a Dirichlet */
          /* if ( C_Z1MEAN(c,t) < dx1_0 || C_PHI2(c,t) < dx2_0) */
          /* 		    { */
          /* 		      strcpy(pdf2Mode, "dirichlet"); */
          /* 		      nDiric++; */
          /* 		    } */
          
          /* this cell has both mixture fractions, compute 2D
             distribution */
          pdfFlag =  MDUCcalcMeanY_2D(pdf2Mode, n, x, nSpecies,
                                      C_PHI1(c,t), vFac*C_PHIV1(c,t),
                                      C_PHI2(c,t), vFac*C_PHIV2(c,t),
                                      vFac*C_PHIV0(c,t), Y, Ymean,
                                      pdf, pdfBound, udf.transform);
          
          /* keep track of errors in PDF */
          if (pdfFlag == -10) {
              nNotCvg++;
              nNotCvgL++;
              ii = nNotCvgL-1;
              
              meanPDF = (double**) realloc(meanPDF, nNotCvgL*sizeof(double*));
              varPDF  = (double**) realloc(varPDF,  nNotCvgL*sizeof(double*));
              meanC   = (double**) realloc(meanC,   nNotCvgL*sizeof(double*));
              varC    = (double**) realloc(varC,    nNotCvgL*sizeof(double*));
              meanPDF[ii] = (double*) calloc(3, sizeof(double));
              varPDF[ii]  = (double*) calloc(4, sizeof(double));
              meanC[ii]   = (double*) calloc(3, sizeof(double));
              varC[ii]    = (double*) calloc(4, sizeof(double));
              
              meanC[ii][0] = C_PHI1(c,t) + C_PHI2(c,t);
              meanC[ii][1] = C_PHI1(c,t);
              meanC[ii][2] = C_PHI2(c,t);
              varC[ii][0]  = C_PHIV0(c,t);
              varC[ii][1]  = C_PHIV1(c,t);
              varC[ii][2]  = C_PHIV2(c,t);
              varC[ii][3]  = 0.5*(varC[ii][0]-varC[ii][1]-varC[ii][2]);
              
              /* compute the moments from the PDF */
              MDUCcalcMomentsPDF_2D(n,x,X,meanPDF[ii],varPDF[ii],
                                    pdf,pdfBound,pdfFlag,udf.transform);
          }
          else if (pdfFlag == 102) {
            nZ1L++;
          }
          else if (pdfFlag == 202) {
            nZ2L++;
          }
          
          /*  add contribution to global pdf */
          for (k=0; k<n2; k++)
            for(j=0; j<n1; j++)
              pdfGlobal[k][j] += pdf[k][j]*mass;
          
          nZ++;
        } /* end 2D contribution */
        else { /* this cell has just Z1, compute 1D distribution */
          pdfFlag = MDUCcalcMeanY_1D(udf.pdf1D, n1, x1, nSpecies,
                                     C_PHI1(c,t), vFac*C_PHIV1(c,t),
                                     Y[0], Ymean, pdfX1, pdfBound1D);
          
          /* add to global pdf */
          for (j=pdfBound1D[0]; j<=pdfBound1D[1]; j++)
            pdfGlobal[0][j] += pdfX1[j]*mass;
          
          nZ1++;
        } /* end 1D for Z1 */
      } /* end computation for cells with Z1 */
      else if ( phiMean[2] > udf.phiZero[2] ) {
        /* this cell has just Z2, compute 1D distribution */
        pdfFlag = MDUCcalcMeanY_1D(udf.pdf1D, n2, x2, nSpecies,
                                   C_PHI2(c,t), vFac*C_PHIV2(c,t),
                                   Y2, Ymean, pdfX2, pdfBound1D);
        
        /* add to global pdf */
        for (k=pdfBound1D[0]; k<=pdfBound1D[1]; k++)
          pdfGlobal[k][0] += pdfX2[k]*mass;
        
        nZ2++;
      } /* end 1D for Z2 */
      else { /* no mixture fraction, update based on the oxidizer species */
             /* essentially delta function at pure oxidizer */
        for (i=0; i < nSpecies; i++)
          Ymean[i] = Y[0][0][i];
        
        /* add mass of cell to oxidizer */
        pdfGlobal[0][0] += mass;
      }
      
      /* ----- update the cell properties ----- */
      /* update the mixture fractions tracked by Fluent */
      Ysum = 0.0;
      for (i=0; i<nYFav; i++) {
        if (indexYFav[i] >= 0)
          C_YI(c,t,i) = Ymean[indexYFav[i]];
        
        if (isnan(C_YI(c,t,i)) || ( (C_YI(c,t,i)>1.0 || C_YI(c,t,i)<0.0) &&
                                    ABS_D(C_YI(c,t,i)) > 1e-2) )
          yOutofBounds = 1;
        
        Ysum += C_YI(c,t,i);
      }
          
      /* update the cell temperature by iterating from the total
         enthalpy and mean mass fractions */
      C_T(c,t) = MDUCcalcTfromH(nSpecies,C_HTOT(c,t),Ymean,C_T(c,t),udf.mduc);
      
      /* error checking on mixture fraction and temperature */
      if (yOutofBounds == 1 || isnan(C_T(c,t))) {
        Message0("species mass fractions for cell %i:\n", c);
        Message0("Z1mean = %15.6e, Z1var = %15.6e\n",
                 C_PHI1(c,t), C_PHIV1(c,t));
        Message0("Z2mean = %15.6e, Z2var = %15.6e\n",
                 C_PHI2(c,t), C_PHIV2(c,t));
        Message0("Zvar = %15.6e\n", C_PHIV0(c,t));
        Message0("Htot = %15.6e\n", C_HTOT(c,t));
        Message0(" --- species --- Y ---\n");
        for (i=0; i<nYFav; i++)
          Message0("%s = %8.6f\n", udf.species[i], C_YI(c,t,i));
        
        pFile = fopen("YmeanError.dat", "w");
        for (i=0; i<nSpecies; i++) fprintf(pFile, "%15.6e\n", Ymean[i]);
        fclose(pFile);
        
        pFile = fopen("pdfError.dat", "w");
        fprintf(pFile, "# %15.6e %15.6e %15.6e %15.6e %15.6e\n",
                C_PHI1(c,t), C_PHI2(c,t), C_PHIV1(c,t), C_PHIV2(c,t),
                C_PHIV0(c,t));
        for (k=0; k<n2; k++) {
          for (j=0; j<n1; j++) {
            fprintf(pFile, "%8.6f %8.6f %15.6e\n", x1[j], x2[k], pdf[k][j]);
          }
          fprintf(pFile, "\n");
        }
        fclose(pFile);
        
        if ( yOutofBounds == 1 )
          Error("Mass fraction has gone out of bounds! (cell %i)\n", c);
        else if (isnan(C_T(c,t)))
          Error("Temperature has become undefined! (cell %i)\n", c);
      }
      
      counter++;
      if (counter%nCount == 0)             Message0(".");
      if (counter%(lineCount*nCount) == 0) Message0("\n");
      
#if DEBUG_SDM
      taui = C_D(c,t)/C_K(c,t); /* inverse turb. timescale */
      
      GradZ[0] = (C_Z1MEAN_G(c,t)[0] + C_Z2MEAN_G(c,t)[0]);
      GradZ[1] = (C_Z1MEAN_G(c,t)[1] + C_Z2MEAN_G(c,t)[1]);
      GradZ[2] = (C_Z1MEAN_G(c,t)[2] + C_Z2MEAN_G(c,t)[2]);
      
      C_Z0VAR_PRD(c,t)   = 2.0*C_MU_T(c,t)/0.7*NV_MAG2(GradZ);
      C_Z1VAR_PRD(c,t)   = 2.0*C_MU_T(c,t)/0.7*NV_MAG2(C_Z1MEAN_G(c,t));
      C_Z2VAR_PRD(c,t)   = 2.0*C_MU_T(c,t)/0.7*NV_MAG2(C_Z2MEAN_G(c,t));
      
      /* dissipation terms with constant coefficient */
      C_Z0VAR_D2(c,t)    = 2.0*taui*C_R(c,t)*C_Z0VAR(c,t);
      C_Z1VAR_D2(c,t)    = 2.0*taui*C_R(c,t)*C_Z1VAR(c,t);
      C_Z2VAR_D2(c,t)    = 2.0*taui*C_R(c,t)*C_Z2VAR(c,t);
      
      C_Z0VAR_SRC_R(c,t) = C_R(c,t)*C_Z0VAR_SRC(c,t);
      C_Z1VAR_SRC_R(c,t) = C_R(c,t)*C_Z1VAR_SRC(c,t);
      C_Z2VAR_SRC_R(c,t) = C_R(c,t)*C_Z2VAR_SRC(c,t);
      C_Z0VAR_EPS_R(c,t) = C_R(c,t)*C_Z0VAR_EPS(c,t)*C_Z0VAR(c,t);
      C_Z1VAR_EPS_R(c,t) = C_R(c,t)*C_Z1VAR_EPS(c,t)*C_Z1VAR(c,t);
      C_Z2VAR_EPS_R(c,t) = C_R(c,t)*C_Z2VAR_EPS(c,t)*C_Z2VAR(c,t);
      
      C_Z0VAR_SRC_L(c,t) *= C_R(c,t)*C_Z0VAR(c,t);
      C_Z1VAR_SRC_L(c,t) *= C_R(c,t)*C_Z1VAR(c,t);
      C_Z2VAR_SRC_L(c,t) *= C_R(c,t)*C_Z2VAR(c,t);
#endif /* DEBUG_SDM */
      
      totMass += mass;
    } end_c_loop_all(c,t);
  } /* end loop over cell threads */
  Message0("\n");

#if RP_NODE
  /* if parallel, get global pdf contributions from all nodes */
  for (k=0; k<n2; k++)
    PRF_GRSUM(pdfGlobal[k], n1, rwork);
  totMass = PRF_GRSUM1(totMass);
  nZ1     = PRF_GISUM1(nZ1);
  nZ2     = PRF_GISUM1(nZ2);
  nZ      = PRF_GISUM1(nZ);
  counter = PRF_GISUM1(counter);
  nStoich = PRF_GISUM1(nStoich);
  nDiric  = PRF_GISUM1(nDiric);
  nNotCvg = PRF_GISUM1(nNotCvgL);
  nNegRho = PRF_GISUM1(nNegRho);
  nRhoP   = PRF_GISUM1(nRhoP);
  nZ1L    = PRF_GISUM1(nZ1L);
  nZ2L    = PRF_GISUM1(nZ2L);
#endif /* RP_NODE */

  /* normalise the global pdf */
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      pdfGlobal[k][j] /= totMass;

  tmp = 100.0/counter;
  Message0("\n Update summary:\n");
  Message0(" ----- cells updated -----  number ----- percent -----\n");
  Message0("\t total     \t %8i \t %6.2f\n", nZ1+nZ2+nZ,
           (double)(nZ1+nZ2+nZ)*tmp);
  Message0("\t Z1        \t %8i \t %6.2f\n", nZ1, (double)(nZ1)*tmp);
  Message0("\t Z2        \t %8i \t %6.2f\n", nZ2, (double)(nZ2)*tmp);
  Message0("\t Z         \t %8i \t %6.2f\n", nZ,  (double)(nZ)*tmp);
  Message0("\t stoich.   \t %8i \t %6.2f\n", nStoich,(double)(nStoich)*tmp);
  if (strcmp(udf.pdf2D, "dirichlet")) {
    tmp = 100.0/nZ; /* want percentages relative to cells with Z */
    Message0("\t Dirichlet \t %8i \t %6.2f\n", nDiric, (double)(nDiric)*tmp);
    Message0("\t Not cvgd. \t %8i \t %6.2f\n", nNotCvg,(double)(nNotCvg)*tmp);
    Message0("\t Neg. rho  \t %8i \t %6.2f\n", nNegRho,(double)(nNegRho)*tmp);
    Message0("\t rho 1     \t %8i \t %6.2f\n", nRhoP,(double)(nRhoP)*tmp);
    Message0("\t nZ1L      \t %8i \t %6.2f\n", nZ1L,(double)(nZ1L)*tmp);
    Message0("\t nZ2L      \t %8i \t %6.2f\n", nZ2L,(double)(nZ2L)*tmp);
  }
  Message0(" -----------------------------------------------------\n");

  int *tmpInt = calloc(1, sizeof(int));
  char errFileName[1000];
  if (nNotCvg > 0) {
#if RP_NODE
    if (myid == node_zero)
#endif /* PARALLEL */
      {    
        strcpy(errFileName, udf.caseDir);
        sprintf(buffer, "/pdf_error_%05i", N_TIME);
        strcat(errFileName, buffer);
        
        pFile = fopen(errFileName, "w");
        
        /* write the header */
        fprintf(pFile, "# id mean1 mean2 var1 var2 covar rho ");
        fprintf(pFile, "meanPDF1 errM1[%%] meanPDF2 errM2[%%] ");
        fprintf(pFile, "varPDF1 errV1[%%] varPDF2  errV2[%%] ");
        fprintf(pFile, "covarPDF errCV[%%]\n");
        
        /* write the non converged cells on the root */
        for (k=0; k<nNotCvgL; k++) {
          for (i=1; i<3; i++) fprintf(pFile, "%15.6e ", meanC[k][i]);
          for (i=1; i<4; i++) fprintf(pFile, "%15.6e ", varC[k][i]);
          fprintf(pFile, "%15.6e ", varC[k][3]/sqrt(varC[k][1]*varC[k][2]));
          for (i=1; i<3; i++)
            fprintf(pFile, "%15.6e %15.6e ", meanPDF[k][i],
                    fabs(meanPDF[k][i]-meanC[k][i])/meanC[k][i]*100);
          for (i=1; i<4; i++)
            fprintf(pFile, "%15.6e %15.6e ", varPDF[k][i],
                    fabs(varPDF[k][i]-varC[k][i])/varC[k][i]*100);
              fprintf(pFile, "\n");
        }

#if RP_NODE
        /* write the info from the other nodes */
        for (k=0; k<nNotCvgL; k++) {
          free(meanPDF[k]); free(varPDF[k]);
          free(meanC[k]); free(varC[k]);
        }
        free(meanPDF); free(varPDF); free(meanC); free(varC);
        
        for (i=1; i<compute_node_count; i++) {
          /* get the number of cells to receive */
          PRF_CRECV_INT(i, tmpInt, 1, i);
          /* reallocate the memory */
          nNotCvgL = tmpInt[0];
          meanPDF = (double**) calloc(nNotCvgL, sizeof(double*));
          varPDF  = (double**) calloc(nNotCvgL, sizeof(double*));
          meanC   = (double**) calloc(nNotCvgL, sizeof(double*));
          varC    = (double**) calloc(nNotCvgL, sizeof(double*));
          
          for (k=0; k<nNotCvgL; k++) {
            meanPDF[k] = (double*) calloc(3, sizeof(double));
            varPDF[k]  = (double*) calloc(4, sizeof(double));
            meanC[k]   = (double*) calloc(3, sizeof(double));
            varC[k]    = (double*) calloc(4, sizeof(double));
            PRF_CRECV_REAL(i, meanPDF[k], 3, i);
            PRF_CRECV_REAL(i, varPDF[k], 4, i);
            PRF_CRECV_REAL(i, meanC[k], 3, i);
            PRF_CRECV_REAL(i, varC[k], 4, i);
            
            for (j=1; j<3; j++) fprintf(pFile, "%15.6e ", meanC[k][j]);
            for (j=1; j<4; j++) fprintf(pFile, "%15.6e ", varC[k][j]);
            fprintf(pFile, "%15.6e ", varC[k][3]/sqrt(varC[k][1]*varC[k][2]));
            for (j=1; j<3; j++)
              fprintf(pFile, "%15.6e %15.6e ", meanPDF[k][j],
                      fabs(meanPDF[k][j]-meanC[k][j])/meanC[k][j]*100);
            for (j=1; j<4; j++)
              fprintf(pFile, "%15.6e %15.6e ", varPDF[k][j],
                      fabs(varPDF[k][j]-varC[k][j])/varC[k][j]*100);
            fprintf(pFile, "\n");
          } /* end loop over cells */
          /* free up the memory again */
          for (k=0; k<nNotCvgL; k++) {
            free(meanPDF[k]); free(varPDF[k]);
            free(meanC[k]); free(varC[k]);
          }
          free(meanPDF); free(varPDF); free(meanC); free(varC);
        } /* end loop over compute nodes */
#endif /* PARALLEL */
        
        fclose(pFile);
      } /* end node zero */
#if RP_NODE
    else {
      tmpInt[0] = nNotCvgL;
      PRF_CSEND_INT(node_zero, tmpInt, 1, myid);
      for (k=0; k<nNotCvgL; k++) {
        PRF_CSEND_REAL(node_zero, meanPDF[k], 3, myid);
        PRF_CSEND_REAL(node_zero, varPDF[k], 4, myid);
        PRF_CSEND_REAL(node_zero, meanC[k], 3, myid);
        PRF_CSEND_REAL(node_zero, varC[k], 4, myid);
        free(meanPDF[k]); free(varPDF[k]);
        free(meanC[k]); free(varC[k]);
      }
      free(meanPDF); free(varPDF); free(meanC); free(varC);
    }
#endif /* PARALLEL */
  } /* end not converged output*/

  for (k=0; k < n2; k++) {
    free(Y2[k]);
    free(pdf[k]); free(pdfBound[k]);
  }
  free(Y2); free(Ymean);
  free(pdf); free(pdfX1); free(pdfX2); free(pdfBound);
  free(rwork);
  for (i=0; i<3; i++) free(X[i]); free(X);
  free(tmpInt);
}
/* ------------------------------------------------------------------------- */

void adjust_T_cond_dTdZ(Domain *d, int n1, double *x1)
{
  int i,j,k, iT;
  Thread *t;
  cell_t  c;
  double phiMean[3];
  double mass, totMass;
  FILE *pFile;

  int  nSpecies = MDUCgetNspecies(udf.mduc);
  double *Ymean = (double*) calloc(nSpecies, sizeof(double));
  double *cp    = (double*) calloc(nSpecies, sizeof(double));
  double *h     = (double*) calloc(nSpecies, sizeof(double));

  /* get the min/max of mixture fraction */
  double EGRmean;
  double EGRmin = 1.0e20, EGRmax = -1.0e20, delEGR;
  totMass = 0.0;
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;

    begin_c_loop_int(c,t) {
      mass = C_R(c,t)*C_VOLUME(c,t); /* local cell mass */
      if (C_Z2MEAN(c,t) < EGRmin) EGRmin = C_Z2MEAN(c,t);
      if (C_Z2MEAN(c,t) > EGRmax) EGRmax = C_Z2MEAN(c,t);
      EGRmean += mass*C_Z2MEAN(c,t);
      totMass += mass;
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
#if RP_NODE
  totMass = PRF_GRSUM1(totMass);
  EGRmin  = PRF_GRLOW1(EGRmin);
  EGRmax  = PRF_GRHIGH1(EGRmax);
#endif /* RP_NODE */
  delEGR   = EGRmax-EGRmin;
  EGRmean /= totMass;
  Message0("EGR mean: %15.6e\n", EGRmean);
  Message0("EGR min/max: %15.6e/%15.6e (%15.6e)\n", EGRmin, EGRmax, delEGR);

  /* compute the average dT/dH */
  double dTdZ;
  double *T_condZ = (double*) calloc(n1, sizeof(double));
  double *T_denZ  = (double*) calloc(n1, sizeof(double));
  double *rwork   = (double*) calloc(n1, sizeof(double));
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;

    begin_c_loop_int(c,t) {
      mass = C_R(c,t)*C_VOLUME(c,t); /* local cell mass */

      phiMean[1] = (C_Z2MEAN(c,t) - EGRmin)/delEGR;

      /* get the location of the current mixture fraction in the grid */
      j = UDFgridLocate(n1, x1, phiMean[1]);

      if (phiMean[1] < 0.5*(x1[j]+x1[j+1])) {
        T_condZ[j] += mass*C_T(c,t);
        T_denZ[j]  += mass;
      }
      else {
        T_condZ[j+1] += mass*C_T(c,t);
        T_denZ[j+1]  += mass;
      }
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
#if RP_NODE
  /* if parallel, get the contribution to chi from each node */
  PRF_GRSUM(T_condZ, n1, rwork);
  PRF_GRSUM(T_denZ, n1, rwork);
#endif /* RP_NODE */  
  /* normalize conditioned temperature over the whole grid */
  for (j=0; j<n1; j++) {
    if (T_denZ[j] > 0.0)
      T_condZ[j] = T_condZ[j]/T_denZ[j];
  } /* end loop over conditional grid */

    /* write to file */
#if RP_NODE
  if (myid == node_zero)
#endif /* RP_NODE */  
    {
      pFile = fopen("T_cond_Z_init.dat", "w");
      fprintf(pFile, "# phiMean T[K]\n");
      for (j=0; j<n1; j++)
        fprintf(pFile, "%15.6e %15.6e\n", x1[j], T_condZ[j]);
      fclose(pFile);
    }

  double T_adjusted, deltaZ, T_EGRmean;
  /* get the temperature at the mean mixture location */
  j = UDFgridLocate(n1,x1,(EGRmean - EGRmin)/delEGR);
  T_EGRmean = T_condZ[j]
    + (EGRmean-x1[j])*(T_condZ[j+1]-T_condZ[j])/(x1[j+1]-x1[j]);

  /* get the slope of the temperature */
  /* neglect ends */
  j = UDFgridLocate(n1, x1, 0.2);
  k = UDFgridLocate(n1, x1, 0.6);
  dTdZ = (T_condZ[k] - T_condZ[j])/((x1[k]-x1[j])*delEGR);
  Message0("dTdZ = %15.6e\n", dTdZ);

  /* adjust the enthalpy based on a linear interpolation to mean mixture */
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;

    begin_c_loop_int(c,t) {

      /* get the current species distribution */
      for (i=0; i<udf.nY; i++) {
        if (udf.indexY[i] >= 0) {
          Ymean[udf.indexY[i]] = C_YI(c,t,i);
        }
      }

      /* find the adjusted temperature */
      T_adjusted = C_T(c,t) - (C_Z2MEAN(c,t) - EGRmean)*dTdZ;

      /* recompute the enthalpy field based on this temperature */
      C_HTOT(c,t) = 0.0;
      MDUCcomputeThermoData(h, cp, T_adjusted, udf.mduc);
      for (i=0; i<nSpecies; i++)
        C_HTOT(c,t) += Ymean[i]*h[i];

    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

  free(T_condZ); free(T_denZ); free(rwork);
  free(Ymean); free(h); free(cp);
}
/* ------------------------------------------------------------------------- */

void computeModH_src(Domain *d)
{
  int i,j,iT;
  Thread *t;
  cell_t  c;

  int     n1 = udf.nCondZ; /* size of conditional grid */
  double *x1 = udf.condZ;  /* grid for computing conditional source */
  double src;              /* source term */
  double mass, Zmean;
  double *denZ, *rwork;

  /* get the min/max of enthalpy */
  /* getEnthalpyBounds(d); */
  
  denZ  = (double*) calloc(n1, sizeof(double));
  rwork = (double*) calloc(n1, sizeof(double));

  /* condition the modified enthalpy source term on mixture fraction */
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;
    
    begin_c_loop_int(c,t) {
      mass = C_R(c,t)*C_VOLUME(c,t); /* local cell mass */

      /* compute the source term for this cell */
      src = udf.dPdt + (C_DPMS_ENERGY(c,t) +
                        udf.H_vap_ref*C_SOURCE_MASS_DPM(c,t))/C_VOLUME(c,t);

      /* divide the source with the local density before averaging.
         Density will be multiplied again in the actual source term */
      src *= 1.0/C_R(c,t);

      /* TODO: figure out how to find if this cell is a boundary cell with a
         head flux that needs to be considered */

      /* TODO: generalize this calculation of Zmean */
      /* interpolate mixture based on Air+EGR as Z=0 and fuel as Z=1 */
      Zmean = C_Z1MEAN(c,t);
      
      /* get the location of the current mixture fraction in the grid */
      j = UDFgridLocate(n1, x1, Zmean);
      
      if (Zmean < 0.5*(x1[j]+x1[j+1])) {
        udf.src_condZ[j] += mass*src;
        denZ[j]          += mass;
      }
      else {
        udf.src_condZ[j+1] += mass*src;
        denZ[j+1]          += mass;
      }
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
#if RP_NODE
  /* if parallel, get the contribution to chi from each node */
  PRF_GRSUM(udf.src_condZ, n1, rwork);
  PRF_GRSUM(denZ, n1, rwork);
#endif /* RP_NODE */  
  /* normalize conditioned temperature over the whole grid */
  for (j=0; j<n1; j++) {
    if (denZ[j] > 0.0)
      udf.src_condZ[j] /= denZ[j];
  } /* end loop over conditional grid */

  free(denZ); free(rwork);
}
/* ------------------------------------------------------------------------- */
#endif /* RPHOST */
/* ##### end functions that only work for serial/compute node ##### */

void checkRIF_coord(char *name, int *iMean, int *iVar)
{
  if  (!strcmp(name, "default")) {
    Message0("Using default RIF co-ordinates\n");
  }
  else if (!strcmp(name, "Z1")) {
    iMean[0] = udf.iZ1;
    iVar[0]  = udf.iZ1V;
  }
  else if (!strcmp(name, "Z2")) {
    iMean[0] = udf.iZ2;
    iVar[0]  = udf.iZ2V;
  }
  else if (!strcmp(name, "H") || !strcmp(name, "enthalpy")) {
    iMean[0] = udf.iH;
    iVar[0]  = udf.iHV;
    udf.rifEnthalpy = 1;
  }
  else if (!strcmp(name, "EGR")) {
    iMean[0] = udf.iEGR;
    iVar[0]  = udf.iEGRV;
  }
  else {
    Message0("Warning: %s co-ordinate undefined for RIF\n", name);
  }
}
/* ------------------------------------------------------------------------- */
/* TODO: This function needs to be checked - Currently not used */
void UDFupdateProgSRC(Domain *d)
{
  Thread *t;
  cell_t c;
  double progMax,progMin, progSRC,*fs, P, p_tab_ref, p_ratio, H, h_tab_ref, h_ratio,facP,facH,facCa,facCb,facCM, cRef, cMaxRef,sum,minRef,*YmeanU,*Yun,ct_x[4],ct_Unburned[4],ct_noC[3],ct_xRef[4];
  int i,j,iT,npar = 11;
  int nSpecies = MDUCgetNspecies(udf.mduc);
  YmeanU= (double*) calloc(nSpecies,   sizeof(double));
  Yun   = (double*) calloc(udf.nTabSp, sizeof(double));
  fs    = (double*) calloc(npar, sizeof(double));
  Message0("Updating Progress Variable Source Term\n"); 
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;
    begin_c_loop_int(c,t) {
             ct_x[0] = C_Z1MEAN(c,t)/(1-C_EGR(c,t));
             ct_x[1] = C_EGR(c,t);
             ct_x[2] = C_Z1VAR(c,t)/(1-C_EGR(c,t))/(1-C_EGR(c,t));
             ct_x[3] = 0;
	     progMax = MDUCchemTableProgMax(ct_x, udf.chemTabPrem);
             progMin = MDUCchemTableProgMin(ct_x, udf.chemTabPrem);
	     ct_x[3] = C_PROG(c,t);

             if (ct_x[3]>progMax) ct_x[3]=progMax;
             if (ct_x[3]<progMin) ct_x[3]=progMin;

     	      if (udf.chemTabPdep && C_PROG(c,t)>1e-2+C_PROGMIN(c,t)) {
     	        /* compute the pressure ratio for the scaling */
     	  
     	            P         = ABS_P(C_P(c,t),op_pres);
                    /* MARCO - limiting pressure expansion to 30 bar */
                    P         = MIN_D(P,30.0e5);
     	           /* H         = C_HTOT(c,t);*/
     	            p_tab_ref = MDUCgetChemTabPref(udf.chemTabPres);
     	            p_ratio   = log(P/p_tab_ref);
     	            
     	            ct_Unburned[0]=ct_x[0];
     	            ct_Unburned[1]=ct_x[1];
     	            ct_Unburned[2]=ct_x[2];
     	            ct_Unburned[3]=0;

                    h_tab_ref = MDUCgetChemTabTref(udf.chemTabPres);
                    /* get unburned species */
                    MDUClookupChemTable("species", ct_Unburned, Yun, udf.chemTabPrem);
      	            for (i=0; i<nSpecies; i++) YmeanU[i] = 0.0;

                    for (i=0; i<udf.nTabSp; i++) { 
                 		if (udf.iTabSp[i] >= 0)  
                    	   YmeanU[udf.iTabSp[i]] = Yun[i];
                    }
                   /* Calc unburned Temperature */
      	            H = MDUCcalcTfromH(nSpecies,C_HTOT(c,t),YmeanU,1000.0,udf.mduc);
                    /* MARCO limiting temperature expansion to 1200K */
                    H = MIN_D(H,1200.0);
     	            h_ratio   = H/h_tab_ref;

     	            ct_noC[0] = ct_x[0];
     	            ct_noC[1] = ct_x[1];
     	            ct_noC[2] = ct_x[2]*0; /*<---------------------------------------------------- no variance for coefficients */
     	            /* SOURCE TERM */
                    
     	            fs[0] = MDUClookupChemTable("SRC_PROG_Par1",ct_noC,NULL,udf.chemTabPres);
     	            fs[1] = MDUClookupChemTable("SRC_PROG_Par2",ct_noC,NULL,udf.chemTabPres);
     	            fs[2] = MDUClookupChemTable("SRC_PROG_Par3",ct_noC,NULL,udf.chemTabPres);
     	            fs[3] = MDUClookupChemTable("SRC_PROG_Par4",ct_noC,NULL,udf.chemTabPres);
     	            fs[4] = MDUClookupChemTable("SRC_PROG_Par5",ct_noC,NULL,udf.chemTabPres);
     	            fs[5] = MDUClookupChemTable("SRC_PROG_Par6",ct_noC,NULL,udf.chemTabPres);
     	            fs[6] = MDUClookupChemTable("SRC_PROG_Par7",ct_noC,NULL,udf.chemTabPres);
     	            fs[7] = MDUClookupChemTable("SRC_PROG_Par8",ct_noC,NULL,udf.chemTabPres);
     	            fs[8] = MDUClookupChemTable("SRC_PROG_Par9",ct_noC,NULL,udf.chemTabPres);
     	            fs[9] = MDUClookupChemTable("SRC_PROG_Par10",ct_noC,NULL,udf.chemTabPres);
     	            fs[10] = MDUClookupChemTable("SRC_PROG_Par11",ct_noC,NULL,udf.chemTabPres);

	            facP    = exp(fs[0]*p_ratio*p_ratio + fs[1]*p_ratio);
                    facCM   =    MAX_D(fs[3]*(h_ratio-1.0)+fs[4]*h_ratio*h_ratio+1.0-fs[4],1e-10);
                    facCa   =    MAX_D(fs[5]*(h_ratio-1.0)+1.0,1e-10);
                    facCb   =    MAX_D(fs[6]*(h_ratio-1.0)+1.0,1e-10);
                      
     	            facH    = exp(fs[8]*log(h_ratio)*log(h_ratio) + fs[7]*log(h_ratio));
     	            cMaxRef = fs[9];
     	            minRef  = fs[10];
     	            /* transform the progress variable to the reference profile */
     	            if (ct_x[3]-cMaxRef*facCM <  0.0)   
     	               cRef = (ct_x[3] - cMaxRef*facCM) / facCa + cMaxRef;
     	            else 
     	               cRef = (ct_x[3] - cMaxRef*facCM) / facCb + cMaxRef;   
     	            /* look up the reference profile */
     	            for (j=0;j<3;j++) ct_xRef[j] = ct_x[j];
     	            ct_xRef[3] = cRef;
     	            /* scale the source term */
     	            if (C_PROG(c,t)<C_PROGMAX(c,t) && C_PROG(c,t)>1e-2+C_PROGMIN(c,t)) {
     	               progSRC = MDUClookupChemTable("SRC_PROG", ct_xRef, NULL,udf.chemTabPrem);
     	               C_PROGSRC(c,t) = MAX_D(progSRC,0.0)*facP*facH; 
                       
     	            }
     	            else {
     	               C_PROGSRC(c,t) = 0;
     	            }
                    /* check if already auto-ignited. If so, set source term to 0 to conserve C_FPVA -> regime tracking */
                    if (C_PROG(c,t)<=C_CRIF(c,t)) C_PROGSRC(c,t) = 0;
               }
               else {
     	       	 C_PROGSRC(c,t) = MDUClookupChemTable("SRC_PROG", ct_x, NULL,udf.chemTabPrem);
                 if (C_PROG(c,t)<=C_CRIF(c,t)) C_PROGSRC(c,t) = 0;
     	       }

    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
  free(YmeanU); free(Yun);free(fs);
}
/*-------------------------------------------------------------------------- */
void UDFupdate(Domain *d)
{
  int    j,n1, nSpecies;
  double *x1;
  double *pdf;         /* global pdf returned from updated after solve */
  double **Y, *Temp;   /* storage for mass fractions and temperature */
  nSpecies = MDUCgetNspecies(udf.mduc);
  n1       = MDUCgetGrid(&x1, udf.mduc);
  Y    = (double**) calloc(n1, sizeof(double*));
  Temp = (double*)  calloc(n1, sizeof(double));
  pdf  = (double*)  calloc(n1, sizeof(double));
  for (j=0; j<n1; j++)
    Y[j] = (double*) calloc(nSpecies, sizeof(double));
  MDUCcopyYandTemp(Y, Temp, udf.mduc);
  updateCFD_1D(d, udf.nY, udf.indexY, n1, x1, Y, Temp, pdf);
}
