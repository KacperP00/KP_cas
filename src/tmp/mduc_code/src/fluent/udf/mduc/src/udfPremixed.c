#include "udfPremixed.h"
#include "mducUDF.h"
#include "udfScalars.h"
#include "udfTools.h"
#include "accessMDUC.h"
#include "math.h"

extern udfStruct udf;
static double Sc_t = 0.7; /* turbulent Schmidt number */

/*
 Source term for progress variable equation:
     rho*Sdot = rho_u*s_t*mag(C)  :  [kg/(m^3-s)]
*/
DEFINE_SOURCE(prog_src, c, t, dS, eqn)
{
  double w, C;      /* source term to return */
  double gradC;     /* magnitude of gradient of progress variable */
  double sL;        /* laminar flame speed (m/s) */
  double sT;        /* turbulent flame speed (m/s) */
  double rho_unb;   /* unburned density (kg/m^3) */
  double c_x[3];

  /* return the progress variable from the table looked-up in update routine */
  return C_PROGSRC(c,t);

  gradC = NV_MAG(C_PROG_G(c,t));

  if (gradC <= 1.0e-9) {
    dS[eqn] = 0.0;
    return 0.0;
  }

  /* TODO: generalize laminar flame speed */
  udf.chemTabX[0] = C_Z1MEAN(c,t);
  udf.chemTabX[1] = C_Z1VAR(c,t);

  /* lookup the unburned density for this mixture */
  udf.chemTabX[2] = 0.0;
  rho_unb = MDUClookupChemTable("RHO", udf.chemTabX, NULL, udf.chemTabPrem);

  /* look-up the laminar flame speed for current progress variable */
  udf.chemTabX[2] = C_PROG(c,t);
  sL = MDUClookupChemTable("sL", udf.chemTabX, NULL, udf.chemTabPrem);
  C_SL(c,t) = sL;

  /* compute the flame speed */
  sT = calcTurbFlameSpeed(sL, c, t);
  C_ST(c,t) = sT;

  /* compute the source term */
  C = rho_unb*sT;
  w = C*NV_MAG(C_PROG_G(c,t));

  /* TODO: should add the gradient 
           d|C|/dC = (C_x*C_xx + C_y*C_yy + C_z*C_zz)/|C| */
  /* dS[eqn] = 0; */

  return w;
}
/* -------------------------------------------------------------------------- */

DEFINE_DIFFUSIVITY(Diff_prog, c, t, i)
{
  return C_MU_T(c,t)/Sc_t;
}
/* ------------------------------------------------------------------------ */

DEFINE_DIFFUSIVITY(Diff_prog_ign, c, t, i)
{
  return C_MU_T(c,t)/Sc_t*C_DIFFSC(c,t);
}
/* -------------------------------------------------------------------------- */

void UDFinitPremixed()
{
  int i,j,k,n,iT;
  double ct_x[4];
  
  Message0(" ===== Initializing Premixed properties ===== \n\n");

  /* get the location of the tables */
  char tabFile[1024], presFile[1024];
  if (RP_is_Defined("mduc/premixed/table")) {
    sprintf(tabFile, "%s/../init/%s", udf.caseDir,
            RP_Get_String("mduc/premixed/table"));
  }
  if (RP_is_Defined("mduc/premixed/pressure-table")) {
    sprintf(presFile, "%s/../init/%s", udf.caseDir,
            RP_Get_String("mduc/premixed/pressure-table"));
  }

  /* error check to make sure tables are present */
  FILE *pFile = NULL;
  udf.chemTabPrem = udf.chemTabPres = NULL;
  if (tabFile) {
    pFile = fopen(tabFile, "r");
    if (pFile == NULL)
      Error("Premixed table %s not found\n", tabFile);
    else
      fclose(pFile);
    
    /* check whether dimensions are defined */
    udf.nDimPrem = 4; /* default */
    if (RP_is_Defined("mduc/premixed/table-dim"))
      n = RP_Get_Integer("mduc/premixed/table-dim");
    if (n > 0)
      udf.nDimPrem = n;
        
    /* load the table into memory */
    udf.chemTabPrem = MDUCloadChemTable(tabFile, udf.nDimPrem, "premixed");

    /* check if there is a pressure table associated.  No need to check outside
       of tabFile check since not valid without tabFile */
    udf.chemTabPdep = 0;
    if (strlen(presFile) > 0) {
      pFile = NULL;
      pFile = fopen(presFile, "r");
      if (pFile == NULL)
        Error("Premixed table %s not found\n", tabFile);
      else
      fclose(pFile);

      /* load the table into memory */
      /* TODO: check whether it always needs the same dimensions */
      udf.chemTabPres = MDUCloadChemTable(presFile,udf.nDimPrem-1, "pressure"); 
      udf.chemTabPdep = 1;
    }
  }

  /* ----- flame speed parameters ----- */
  if (RP_is_Defined("mduc/premixed/prog-threshold"))
    udf.a_thresh = RP_Get_Real("mduc/premixed/prog-threshold");
  if (RP_is_Defined("mduc/premixed/prog-src-scale"))
    udf.progSrcScale = RP_Get_Real("mduc/premixed/prog-src-scale");
  if (RP_is_Defined("mduc/premixed/turb-flame-speed"))
    udf.sTmodel = RP_Get_String("mduc/premixed/turb-flame-speed");
  if (!strncmp(udf.sTmodel, "zim", 3))
    calcTurbFlameSpeed = &calcZimontTurbFlameSpeed;
  else if (!strncmp(udf.sTmodel, "pet", 3))
    calcTurbFlameSpeed = &calcPetersTurbFlameSpeed;
  else if (!strcmp(udf.sTmodel, "off") || !strcmp(udf.sTmodel, "none"))
    calcTurbFlameSpeed = &calcNoTurbFlameSpeed;
  else {
    Message0("Warning: turbulent flame speed model %s not recognized!\n",
             udf.sTmodel);
  }

  /* ----- get any ignition parameters ----- */
  udf.ignited = 0;
  /* type and location */
  if (RP_is_Defined("mduc/premixed/ignition-time"))
    udf.ignTime = RP_Get_Real("mduc/premixed/ignition-time");
  if (RP_is_Defined("mduc/premixed/ignition-geom"))
    udf.ignGeom = RP_Get_String("mduc/premixed/ignition-geom");
  if (RP_is_Defined("mduc/premixed/ignition-loc"))
    n = parseFluentReal(RP_Get_String("mduc/premixed/ignition-loc"),
                        &udf.ignCoord);
  if (n < 3)
    Message0("Warning: only %i co-ordinates for ignition location given!\n");

  /* initial size */
  if (RP_is_Defined("mduc/premixed/ignition-radius"))
    udf.ignRadius = RP_Get_Real("mduc/premixed/ignition-radius");
  if (RP_is_Defined("mduc/premixed/ignition-size"))
    n = parseFluentReal(RP_Get_String("mduc/premixed/ignition-size"),
                        &udf.ignSize);

  /* error checking */
  if (!strncmp(udf.ignGeom, "spher", 5)) {
    if (udf.ignRadius <= 0.0)
      Message0("Warning: ignition radius not specified\n");
  }
  else if (!strncmp(udf.ignGeom, "rect", 4)) {
    if (n < 3) {
      Message0("Warning: not enough sides specified for ignition region\n");
    }
  }
  else {
    Error("Ignition geometry %s not recognized!", udf.ignGeom);
  }
  if (RP_is_Defined("mduc/premixed/ignition-c-scale"))
    udf.ignScale = RP_Get_Real("mduc/premixed/ignition-c-scale");

  /* ----- initialize any variables ----- */
  udf.chemTabX = (double*) calloc(udf.nDimPrem, sizeof(double));
  /* get the maximum progress variable */
  udf.progMax = MDUCchemTableVarMax("PROG", udf.chemTabPrem);

  /* compute the representative gradient of the progress variable */
  double L_f, cmax;
  ct_x[0] = 0.0;
  ct_x[1] = 0.0;
  ct_x[2] = 0.0;
  ct_x[3] = 0.0;

  /* get the species from the table */
  /* TODO: only works with compiled mechanism for now */
  if (udf.mduc == NULL) {
    Message0(" Initializing chemistry from mechanism\n");
    /* TODO: this should only be called if not using RIF as well,
       only need chemistry for now */
   /* udf.mduc = MDUCinitChem(NULL);
    if (udf.mduc == NULL) {
      Message0("why is this still NULL?\n");
      return;
    }*/
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
    
    /* initialise the working memory */
    Message0("  Initializing MDUC from file: %s in UDFPremixed.c\n", fileName);
    udf.mduc = MDUCinit(fileName);
    Message0("  MDUC initialized");
    MDUCsetSpeciesIndex(udf.nY, udf.species, &udf.indexY, udf.mduc);
    udf.Wmol = (double*) calloc(MDUCgetNspecies(udf.mduc), sizeof(double));
    MDUCgetMolWeightAll(udf.Wmol, udf.mduc);
  }

  /* get the species from the table */
  udf.nTabSp = MDUCgetChemTabNspecies(udf.chemTabPrem);
  char **tabSp   = (char**) calloc(udf.nTabSp, sizeof(char*));
  for (i=0; i<udf.nTabSp; i++)
    tabSp[i] = (char*) calloc(64, sizeof(char));
  MDUCgetChemTabSpecies(tabSp, udf.chemTabPrem);
  MDUCsetSpeciesIndex(udf.nTabSp, tabSp, &udf.iTabSp, udf.mduc);

  /* get the definition of progress variable */
  if (RP_is_Defined("mduc/prog-var-species"))
    udf.nProgSp = parseFluentString(RP_Get_String("mduc/prog-var-species"),
                                    &udf.progSp);
  if (udf.nProgSp < 0) {
    Message0("Warning: using default progress variable species\n");
    udf.nProgSp = 4;
    udf.progSp = (char**) calloc(udf.nProgSp, sizeof(char*));
    for (i=0; i<udf.nProgSp; i++)
      udf.progSp[i] = (char*) calloc(64, sizeof(char));
    strcpy(udf.progSp[0], "CO2");
    strcpy(udf.progSp[1], "H2O");
    strcpy(udf.progSp[2], "CO");
    strcpy(udf.progSp[3], "H2");
  }
  MDUCsetSpeciesIndex(udf.nProgSp, udf.progSp, &udf.iProgSp, udf.mduc);
  
  Message0(" ========== Species indices in mechanism ==========\n");
  Message0("   Tracked species:\n");
  int specWarn = 0;
  for (i=0; i<udf.nY; i++) {
    Message0("   %2i %20s %4i\n",i,udf.species[i],udf.indexY[i]);
    if (udf.indexY[i] < 0) specWarn = 1;
  }
  if (specWarn)
    Message0("\n  Warning: some species not found in mechanism!\n");

  Message0("   Table species:\n");
  for (i=0, specWarn=0; i<udf.nTabSp; i++) {
    Message0("   %2i %20s %4i\n",i,tabSp[i],udf.iTabSp[i]);
    if (udf.iTabSp[i] < 0) specWarn = 1;
  }
  if (specWarn)
    Message0("\n  Warning: some species not found in mechanism!\n");

  Message0("   Progress variable species:\n");
  for (i=0, specWarn=0; i<udf.nProgSp; i++) {
    Message0("   %2i %20s %4i\n",i,udf.progSp[i],udf.iProgSp[i]);
    if (udf.iTabSp[i] < 0) specWarn = 1;
  }
  if (specWarn)
    Message0("\n  Warning: some species not found in mechanism!\n");
  Message0(" ===================================================\n");
  L_f = MDUClookupChemTable("lF", ct_x, NULL, udf.chemTabPrem);
  cmax = MDUCchemTableProgMax(ct_x, udf.chemTabPrem);
  udf.gradC_prem = cmax/L_f;

#if !RP_HOST
  /* look up the species if necesary */
  if (!udf.restart && !udf.rif) {
    Domain *d = Get_Domain(1);
    /* compute the initial enthalpy from the species/temperature */
    UDFupdateEnthalpy(d, udf.mduc);

    /* update the solution from the table */
    UDFpremixedUpdateCFD(d);
  }
#endif /* !RP_HOST */

  /* ----- output summary ----- */
  Message0(" \n");
  Message0("  Maximum progress variable: %.3f\n", udf.progMax);
  Message0("  Gradient progress var:     %.6e\n", udf.gradC_prem);
  Message0("  Gradient threshold:        %.6e\n", udf.a_thresh);
  Message0("  Progress var scaling:      %.6e\n", udf.progSrcScale);
  Message0("  Progress ignition scaling: %.6e\n", udf.ignScale);
  Message0(" \n");
  Message0("  Ignition characteristics:\n");
  if (in_cylinder_p)
    Message0("    CAD:      %.2f CAD\n", udf.ignTime);
  else
    Message0("    time:     %.6f s\n", udf.ignTime);
  Message0("    shape:    %s\n", udf.ignGeom);
  Message0("    location: (%.6f, %.6f, %.6f)\n", udf.ignCoord[0],
           udf.ignCoord[1], udf.ignCoord[2]);
  if (!strncmp(udf.ignGeom, "spher", 5))
    Message0("    radius:   %.6f\n", udf.ignRadius);
  else if (!strncmp(udf.ignGeom, "rect", 4))
    Message0("    LxWxH:    (%.6f, %.6f, %.6f)\n",udf.ignSize[0],udf.ignSize[1],
             udf.ignSize[2]);

  Message0("\n ============================================\n");
  if (udf.chemTabPdep) {
    Message0("  Table reference Pressure: %.6e\n",MDUCgetChemTabPref(udf.chemTabPres));
    Message0("  Table reference Temperature: %.6e\n",MDUCgetChemTabTref(udf.chemTabPres));
  }
  for (i=0; i<udf.nTabSp; i++) free(tabSp[i]);
  free(tabSp);
}
/* -------------------------------------------------------------------------- */

void igniteFlameKernel()
{
  int iT;
  Thread *t;
  cell_t c;

  double r, x[ND_ND], tmp2[3], scale, WmixInv;
  double L_f     = 3e-5;      /* flame thickness */
  double progMax;  /* maximum progress variable */
  double progMin, mass_c;
  
  int    i,j,k;
  int    nY, nSp  = MDUCgetNspecies(udf.mduc);
  int    nVar     = MDUCgetChemTabNvar(udf.chemTabPrem);
  double *Ymean   = (double*) calloc(nSp,  sizeof(double));
  double *Ytab    = (double*) calloc(udf.nTabSp, sizeof(double));
  double *h       = (double*) calloc(nSp,  sizeof(double));
  double *cp      = (double*) calloc(nSp,  sizeof(double));

  int ignCell = 0;
  int ignType = -1;
  if (!strncmp(udf.ignGeom, "spher", 5))
    ignType = 1;
  else if (!strncmp(udf.ignGeom, "rect", 4))
    ignType = 2;

  /* only loop over solution threads (don't track flame in non-reactive region*/
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;

    begin_c_loop_int(c,t) {
      /* get the location of current cell */
      C_CENTROID(x,c,t);

      /* look up the flame thickness */
      udf.chemTabX[0] = C_Z1MEAN(c,t)/(1-C_EGR(c,t));
      udf.chemTabX[1] = C_EGR(c,t);
      udf.chemTabX[2] = C_Z1VAR(c,t)/(1-C_EGR(c,t))/(1-C_EGR(c,t));
      udf.chemTabX[3] = C_PROG(c,t);
      L_f = MDUClookupChemTable("lF", udf.chemTabX, NULL, udf.chemTabPrem); 
      progMin = MDUCchemTableProgMin(udf.chemTabX, udf.chemTabPrem);
      progMax = MDUCchemTableProgMax(udf.chemTabX, udf.chemTabPrem);
      
      /* check to see if cell is in ignition region */
      ignCell = 0;
      if (ignType == 1) { /* spherical */
        tmp2[0] = x[0] - udf.ignCoord[0]; 
        tmp2[1] = x[1] - udf.ignCoord[1]; 
        tmp2[2] = x[2] - udf.ignCoord[2]; 
        r = pow(tmp2[0]*tmp2[0]+tmp2[1]*tmp2[1]+tmp2[2]*tmp2[2], 0.5);

        scale = 0.5*(1.0+tanh((r - udf.ignRadius)*8.0/L_f));
        C_PROG(c,t) = udf.ignScale*MAX_D((1.0-scale)*(progMax-progMin)+progMin, progMin);
        ignCell = 1;
      }
      else if (ignType == 2) { /* rectangle */
        /* TODO: put a profile here.  For now just set everything directly */
        if ( (x[0] > udf.ignCoord[0] - 0.5*udf.ignSize[0]) &&
             (x[0] < udf.ignCoord[0] + 0.5*udf.ignSize[0]) &&
             (x[1] > udf.ignCoord[1] - 0.5*udf.ignSize[1]) &&
             (x[1] < udf.ignCoord[1] + 0.5*udf.ignSize[1]) &&
             (x[2] > udf.ignCoord[2] - 0.5*udf.ignSize[2]) &&
             (x[2] < udf.ignCoord[2] + 0.5*udf.ignSize[2]) )
          C_PROG(c,t) = udf.ignScale*progMax;
        ignCell = 1;
      }

      /* update properties in ignited cell */
      if (ignCell) {
        udf.chemTabX[3] = C_PROG(c,t);

        /* lookup the mean species and temperature from the table */
        MDUClookupChemTable("species", udf.chemTabX, Ytab, udf.chemTabPrem);

        /* lookup the species from the chemtable variables */
        for (i=0; i<nSp; i++) Ymean[i] = 0.0;
        for (i=0; i<udf.nTabSp; i++) {
          if (udf.iTabSp[i] > 0)
            Ymean[udf.iTabSp[i]] = Ytab[i];
        }

        /* update the tracked species */
        for (i=0; i<udf.nY; i++) {
          C_YI(c,t,i) = 0.0;
          if (udf.indexY[i] >= 0)
            C_YI(c,t,i) = Ymean[udf.indexY[i]];
        }

        WmixInv = 0.0;
        for (i=0; i<nSp; i++)
          WmixInv += Ymean[i]/udf.Wmol[i];

        /* compute the temperature from the enthlapy */
        C_T(c,t) = MDUCcalcTfromH(nSp,C_HTOT(c,t),Ymean,1000.0,udf.mduc);
        /* C_T(c,t) = MDUClookupChemTable("T",udf.chemTabX,NULL,udf.chemTabPrem); */

        /* ensure that the density is consistent */
/*         C_R(c,t) = ABS_P(C_P(c,t),op_pres)/ */
/*           (UNIVERSAL_GAS_CONSTANT*WmixInv*C_T(c,t)); */
        /* C_R(c,t) = MDUClookupChemTable("RHO",udf.chemTabX,NULL,udf.chemTabPrem); */
      }
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
  free(Ymean); free(Ytab); free(h); free(cp);
}
/* -------------------------------------------------------------------------- */

void flameKernelDev()
{
  Thread *t;
  cell_t c;
  double rk, mtot,mass_c,z1,z2,z1var,tmp2[3],x[ND_ND],htot,r;
  double ct_x[4],ct_xRef[4],c_xu[4],c_xb[4],c_xnoC[3];
  double rho_u,rho_b,cp,progMin,progMax,L_f,dt;
  double Tb,dmkdt,dTkdt,sT,eta_eff,dQdt,curT;
  double pi,scale,sL,fp1,fp2,fT1,fT2;
  double P,P_tab_ref,p_ratio,facP;
  double T,T_tab_ref,T_ratio,facT;
  double facCM,facCa,facCb,cMaxRef,minRef,cRef, Pavg;
  int i,iT,ignCell,j,npar = 11;
  double  C_star_q, lt, C_mu = 0.09; 
  double Cd = 0.37; /* default constant for turb. length scale */
  double C_s = 2.0;
  double Dt,D = 0;  /* diffusivity neglegted against turbulent diffusivity */
  double l_star_q,uprime,k_spk;
  double *Ymean,*Ytab,*fp,*fpS;
  char tabSpI[64];
  char **tabSp   = (char**) calloc(udf.nTabSp, sizeof(char*));
  for (i=0; i<udf.nTabSp; i++)
    tabSp[i] = (char*) calloc(64, sizeof(char));
  MDUCgetChemTabSpecies(tabSp, udf.chemTabPrem);
  int nSpecies= MDUCgetNspecies(udf.mduc); 
  Ymean = (double*) calloc(nSpecies,   sizeof(double));
  Ytab  = (double*) calloc(udf.nTabSp, sizeof(double));
  fp    = (double*) calloc(udf.nTabSp*npar, sizeof(double));
  fpS = (double*) calloc(npar, sizeof(double));
  /* ---------------CONSTANTS ------------*/
  pi = 3.14159265359;

  if (udf.ignited == 0) udf.rk = udf.ignRadius; 

  mass_c= 0;z1=0;z2=0;z1var=0;mtot=0;htot=0,sT=0;Pavg=0;Dt=0;
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;

    begin_c_loop_int(c,t) {
      /* get the location of current cell */
      C_CENTROID(x,c,t);
      tmp2[0] = x[0] - udf.ignCoord[0]; 
      tmp2[1] = x[1] - udf.ignCoord[1]; 
      tmp2[2] = x[2] - udf.ignCoord[2]; 
      r = pow(tmp2[0]*tmp2[0]+tmp2[1]*tmp2[1]+tmp2[2]*tmp2[2], 0.5);
      
      /* Averaging table inputs */
      if (r<=udf.rk) {
        mass_c = C_R(c,t)*C_VOLUME(c,t);
        ct_x[0]= C_Z1MEAN(c,t)/(1-C_EGR(c,t));
        ct_x[1]= C_EGR(c,t);
        ct_x[2]= C_Z1MEAN(c,t)/(1-C_EGR(c,t));
        ct_x[3]= C_PROG(c,t);
        sL    = MDUClookupChemTable("sL", ct_x, NULL, udf.chemTabPrem);
        if (udf.chemTabPdep) {
          
          ct_x[0]= C_Z1MEAN(c,t)/(1-C_EGR(c,t));
          ct_x[1]= C_EGR(c,t);
          ct_x[2]= C_Z1MEAN(c,t)/(1-C_EGR(c,t))/(1-C_EGR(c,t));
          ct_x[3]= 0;
          
          for (i=0;i<3;i++) c_xnoC[i] = ct_x[i];
          /*Expansion Parameter*/
          fp1 = MDUClookupChemTable("sL_Par1", c_xnoC, NULL, udf.chemTabPres);
          fp2 = MDUClookupChemTable("sL_Par2", c_xnoC, NULL, udf.chemTabPres);
          fT1 = MDUClookupChemTable("sL_Par4", c_xnoC, NULL, udf.chemTabPres);
          fT2 = MDUClookupChemTable("sL_Par5", c_xnoC, NULL, udf.chemTabPres);
          
          /*Pressure*/
          P         = ABS_P(C_P(c,t),op_pres);
          P_tab_ref = MDUCgetChemTabPref(udf.chemTabPres);
          p_ratio   = log(P/P_tab_ref);
          facP      = exp(fp1*p_ratio*p_ratio + fp2*p_ratio);
          /*Temperature*/
          T_tab_ref = MDUCgetChemTabTref(udf.chemTabPres);
          /* get unburned species */
          MDUClookupChemTable("species", ct_x, Ytab, udf.chemTabPrem);
          for (i=0; i<nSpecies; i++) Ymean[i] = 0.0;
          
          for (i=0; i<udf.nTabSp; i++) {
            if (udf.iTabSp[i] >= 0)
              Ymean[udf.iTabSp[i]] = Ytab[i];
          }
          T         = MDUCcalcTfromH(nSpecies,C_HTOT(c,t),Ymean,1000,udf.mduc);
          T_ratio   = log(T/T_tab_ref);
          facT      = exp(fT2*T_ratio*T_ratio + fT1*T_ratio);
          sL       *= facT*facP;  
        }
        uprime   = pow(2.0*C_K(c, t)/3.0,0.5);
        lt    += Cd*uprime*uprime*uprime/(C_D(c, t)+1.e-6) * mass_c;
        sT    += calcPetersTurbFlameSpeedIgn(sL, c, t)*mass_c;
        z1    += C_Z1MEAN(c,t)*mass_c;
        z2    += C_EGR(c,t)*mass_c;
        z1var += C_Z1VAR(c,t)*mass_c;
        mtot  += mass_c;
        htot  += C_HTOT(c,t)*mass_c;
        Pavg  += ABS_P(C_P(c,t),op_pres)*mass_c; 
        Dt    += C_MU_T(c,t)/Sc_t*C_DIFFSC(c,t)*mass_c; 
      }
    } end_c_loop_int(c,t);
  }
  
#if RP_NODE
 /* if parallel, synchronize all node solutions here by summing the
    contributions from each node */
  mtot = PRF_GRSUM1(mtot); 
  z1   = PRF_GRSUM1(z1);
  z2   = PRF_GRSUM1(z2);
  z1var= PRF_GRSUM1(z1var);
  htot = PRF_GRSUM1(htot);
  sT   = PRF_GRSUM1(sT);
  Pavg = PRF_GRSUM1(Pavg);
  Dt   = PRF_GRSUM1(Dt);
  lt   = PRF_GRSUM1(lt);
#endif /* RP_NODE */
  z1   /= mtot;
  z2   /= mtot;
  z1var/= mtot;
  htot /= mtot;
  sT   /= mtot;
  Pavg /= mtot;
  Dt   /= mtot;
  lt   /= mtot;
 
  Message0("--------- Flame Kernel --------\n");
  Message0("	z1   :%15.10e\n", z1    ); 
  Message0("	z2   :%15.10e\n", z2    );
  Message0("	z1var:%15.10e\n", z1var );
  Message0("	htot :%15.10e\n", htot  );
  Message0("	sT   :%15.10e\n", sT    );
  Message0("	Pavg :%15.10e\n", Pavg  );
  Message0("	Dt   :%15.10e\n", Dt    );
  Message0("	lt   :%15.10e\n", lt    );
  ct_x[0] = z1/(1-z2);
  ct_x[1] = z2;
  ct_x[2] = z1var/(1-z2)/(1-z2);
  ct_x[3] = 0;

  progMax = MDUCchemTableProgMax(ct_x, udf.chemTabPrem);
  progMin = MDUCchemTableProgMin(ct_x, udf.chemTabPrem);
  /* lookup the burned density for this mixture */

  for (i=0;i<3;i++) c_xb[i] = ct_x[i];
  for (i=0;i<3;i++) c_xu[i] = ct_x[i];
  for (i=0;i<3;i++) c_xnoC[i] = ct_x[i];
  c_xb[3] = progMax;
  c_xu[3] = progMin;

  /* TODO: Expansion for cp and diffusivity */
  
/*
  if (udf.chemTabPdep) {
       P         = Pavg;
       P_tab_ref = MDUCgetChemTabPref(udf.chemTabPres);
       p_ratio   = log(P/P_tab_ref);
       
       T_tab_ref = MDUCgetChemTabTref(udf.chemTabPres);
*/
       /* get unburned species */
/*
       MDUClookupChemTable("species", c_xu, Ytab, udf.chemTabPrem);
       for (i=0; i<nSpecies; i++) Ymean[i] = 0.0;
  
       for (i=0; i<udf.nTabSp; i++) { 
              if (udf.iTabSp[i] >= 0) Ymean[udf.iTabSp[i]] = Ytab[i];
       }
*/  
       /* Calc unburned Temperature */
/*       T = MDUCcalcTfromH(nSpecies,htot,Ymean,1000,udf.mduc); */ /* <-------------------------- Adjust to kernel temperature */
/*
       Message0("Unburned temperature %10.5e\n",T);
       T_ratio   = T/T_tab_ref;
       for (i=0; i<3;i++) c_xnoC[i] = ct_x[i];
       MDUClookupChemTable("species", c_xnoC, fp, udf.chemTabPres);
        
       for (i=0; i<udf.nTabSp; i++) {
                   facP    = exp(fp[i*npar+0]*p_ratio*p_ratio + fp[i*npar+1]*p_ratio);
                   facCM   =  MAX_D(fp[i*npar+3]*(T_ratio-1.0)+fp[i*npar+4]*T_ratio*T_ratio+1.0-fp[i*npar+4],1e-10);
                   facCa   =    MAX_D(fp[i*npar+5]*(T_ratio-1.0)+1.0,1e-10);
                   facCb   =    MAX_D(fp[i*npar+6]*(T_ratio-1.0)+1.0,1e-10);
                   facT    = exp(fp[i*npar+8]*log(T_ratio)*log(T_ratio) + fp[i*npar+7]*log(T_ratio));
                   cMaxRef =     fp[i*npar+9];
                   minRef  =     fp[i*npar+10];
*/                
                 /* transform the progress variable to the reference profile */
/*                    if (ct_x[3]-cMaxRef*facCM <  0.0 )   
                       cRef = (ct_x[3] - cMaxRef*facCM) / facCa + cMaxRef;
                    else 
                       cRef = (ct_x[3] - cMaxRef*facCM) / facCb + cMaxRef;   
*/
                 /* look up the reference profile */
/*                    for (j=0;j<3;j++) ct_xRef[j] = ct_x[j];
                    ct_xRef[3] = cRef;

                    strcpy(tabSpI, "Y_");
                    strcat(tabSpI,tabSp[i]);
                 
                    Ytab[i] = MDUClookupChemTable(tabSpI, ct_xRef, NULL, udf.chemTabPrem);
                    Ytab[i] = facP*((Ytab[i]-minRef)*facT+minRef);
  
  	for (i=0;i<npar;i++) fpS[i] = MDUClookupChemTable(sprintf("RHOPar","%i",i), c_xnoC, NULL, udf.chemTabPres);
        facP    = exp(fpS[0]*p_ratio*p_ratio + fpS[1]*p_ratio);
        facCM   = MAX_D(fpS[3]*(T_ratio-1.0)+fpS[4]*T_ratio*T_ratio+1.0-fpS[4],1e-10);
        facCa   = MAX_D(fpS[5]*(T_ratio-1.0)+1.0,1e-10);
        facCb   = MAX_D(fpS[6]*(T_ratio-1.0)+1.0,1e-10);
        facT    = exp(fpS[8]*log(T_ratio)*log(T_ratio) + fpS[7]*log(T_ratio));
        cMaxRef = fpS[9];
        minRef  = fpS[10];
*/
     /* Burned density*/
        /* transform the progress variable to the reference profile */
/*        if (c_xb[3]-cMaxRef*facCM <  0.0 )
           cRef = (c_xb[3] - cMaxRef*facCM) / facCa + cMaxRef;
        else 
           cRef = (c_xb[3] - cMaxRef*facCM) / facCb + cMaxRef;
*/
        /* look up the reference profile */
/*        for (j=0;j<3;j++) ct_xRef[j] = c_xb[j];
        ct_xRef[3] = cRef;  
  	rho_b  = MDUClookupChemTable("RHO", c_xb, NULL, udf.chemTabPrem);
        rho_b  = facP*((rho_b-minRef)*facT+minRef);
*/     /* Unburned density*/  
        /* transform the progress variable to the reference profile */
/*        if (c_xu[3]-cMaxRef*facCM <  0.0 )
           cRef = (c_xu[3] - cMaxRef*facCM) / facCa + cMaxRef;
        else 
           cRef = (c_xu[3] - cMaxRef*facCM) / facCb + cMaxRef;
*/        /* look up the reference profile */
/*        for (j=0;j<3;j++) ct_xRef[j] = c_xu[j];
        ct_xRef[3] = cRef;  
  	rho_u  = MDUClookupChemTable("RHO", c_xb, NULL, udf.chemTabPrem);
        rho_u  = facP*((rho_u-minRef)*facT+minRef);
*/     /* Diffusivity */
/*        D    = MDUClookupChemTable("diffusivity", ct_x, NULL, udf.chemTabPrem);
*/        /*Expansion Parameter*/
/*        fp1 = MDUClookupChemTable("diffusivity_Par1", c_xnoC, NULL, udf.chemTabPres);
        fp2 = MDUClookupChemTable("diffusivity_Par2", c_xnoC, NULL, udf.chemTabPres);
        fT1 = MDUClookupChemTable("diffusivity_Par4", c_xnoC, NULL, udf.chemTabPres);
        fT2 = MDUClookupChemTable("diffusivity_Par5", c_xnoC, NULL, udf.chemTabPres);
        
*/        /*Pressure*/
/*        facP      = exp(fp1*p_ratio*p_ratio + fp2*p_ratio);
*/        /*Temperature*/
/*	facT      = exp(fT2*log(T_ratio)*log(T_ratio) + fT1*log(T_ratio));
*/        /* Scaling */
/*        D        *= facT*facP;  
        
*/        /* Cp */
/*  	cp     = MDUClookupChemTable("CP", c_xb, NULL, udf.chemTabPrem);
*/        /*Expansion Parameter*/
/*        fp1 = MDUClookupChemTable("CP_Par1", c_xnoC, NULL, udf.chemTabPres);
        fp2 = MDUClookupChemTable("CP_Par2", c_xnoC, NULL, udf.chemTabPres);
        fT1 = MDUClookupChemTable("CP_Par4", c_xnoC, NULL, udf.chemTabPres);
        fT2 = MDUClookupChemTable("CP_Par5", c_xnoC, NULL, udf.chemTabPres);
        
*/        /*Pressure*/
/*        facP      = exp(fp1*p_ratio*p_ratio + fp2*p_ratio);
*/        /*Temperature*/
/*	facT      = exp(fT2*log(T_ratio)*log(T_ratio) + fT1*log(T_ratio));
*/        /*Scaling*/
/*        D       *= facT*facP;  
       }
  }

  else {
*/
  	MDUClookupChemTable("species", c_xb, Ytab, udf.chemTabPrem);
  	rho_b  = MDUClookupChemTable("RHO", c_xb, NULL, udf.chemTabPrem);
  	rho_u  = MDUClookupChemTable("RHO", c_xu, NULL, udf.chemTabPrem);
  	cp     = MDUClookupChemTable("CP", c_xb, NULL, udf.chemTabPrem);
  	/*D      = MDUClookupChemTable("diffusivity", c_xb, NULL, udf.chemTabPrem);*/
/*  }*/
  /* lookup the species from the chemtable variables */
  for (i=0; i<nSpecies; i++) Ymean[i] = 0.0;
  for (i=0; i<udf.nTabSp; i++) {
    if (udf.iTabSp[i] > 0) Ymean[udf.iTabSp[i]] = Ytab[i];
  }
  Tb      = MDUCcalcTfromH(nSpecies,htot,Ymean,1000.0,udf.mduc);
  udf.Tad = Tb;

  /* Update cells */
  /* INITIALIZATION */
  if (udf.ignited == 0) {
    udf.t  = 0.0;
    udf.mk = rho_b*4/3*pi*pow(udf.rk,3);
    udf.Tk = Tb;
    L_f = MDUClookupChemTable("lF", ct_x, NULL, udf.chemTabPrem); 
    for (iT=0; iT<udf.nThreadSolve; iT++) { 
      t = udf.solve_threads[iT]; /* current thread */
      if (t == NULL) continue;
      
      begin_c_loop_int(c,t) {
        /* get the location of current cell */
        C_CENTROID(x,c,t);
	
        /* check to see if cell is in ignition region */
        tmp2[0] = x[0] - udf.ignCoord[0]; 
        tmp2[1] = x[1] - udf.ignCoord[1]; 
        tmp2[2] = x[2] - udf.ignCoord[2]; 
        r = pow(tmp2[0]*tmp2[0]+tmp2[1]*tmp2[1]+tmp2[2]*tmp2[2], 0.5);
        
        ct_x[0]= C_Z1MEAN(c,t)/(1-C_EGR(c,t));
        ct_x[1]= C_EGR(c,t);
        ct_x[2]= C_Z1MEAN(c,t)/(1-C_EGR(c,t))/(1-C_EGR(c,t));
        ct_x[3]= 0;
        
        progMax = MDUCchemTableProgMax(ct_x, udf.chemTabPrem);
        progMin = MDUCchemTableProgMin(ct_x, udf.chemTabPrem);
	
        scale = 0.5*(1.0+tanh((r - udf.ignRadius)*8.0/L_f));
        C_PROG(c,t) = udf.ignScale*MAX_D((1.0-scale)*(progMax-progMin)+progMin, progMin);
        /* C_ALPHA(c,t) = udf.ignScale*MAX_D((1.0-scale),0); */
        C_DIFFSC(c,t) = 0;              
      } end_c_loop_int(c,t);
    } /* end loop over cell threads */
    /* INITIALIZATION DONE */ 
    udf.ignited = 1;
    Dt          = 0;
  } 
  else {
    for (iT=0; iT<udf.nThreadSolve; iT++) { 
      t = udf.solve_threads[iT]; /* current thread */
      if (t == NULL) continue;
      
      begin_c_loop_int(c,t) {
        uprime   = pow(2.0*C_K(c, t)/3.0,0.5);
        lt       = Cd*uprime*uprime*uprime/(C_D(c, t)+1.e-6);
        k_spk    = C_K(c,t) * MIN_D(1,pow(2*udf.rk/lt,2.0/3.0));
        l_star_q = 1-exp(-C_s*C_D(c,t)/C_K(c,t)*udf.t);
        C_DIFFSC(c,t) = k_spk*k_spk/C_K(c,t)/C_K(c,t)*l_star_q;
      } end_c_loop_int(c,t);
    } /* end loop over cell threads */
  }
  /* UPDATE KERNEL*/
  sT  = MAX_D(sT-2/udf.rk*(D + Dt),0);
  Message0("	sT_spk   %15.10e\n", sT);
  Message0("	mk (old) %15.10e\n", udf.mk);
  udf.rk = pow((3.0*udf.mk/(4.0*pi*rho_b)),1.0/3.0);
  Message0("	rk       %15.10e\n", udf.rk);
  dt     = CURRENT_TIMESTEP;
  dmkdt  = 4.0*pi*pow(udf.rk,2.0)*rho_u*sT;
  /* TODO: Temperature rise because of heat flux */
  if (udf.ignited==0) {
    dQdt        = 2e-3/CURRENT_TIMESTEP;
    eta_eff     = (0.94+0.5)/2.0;
    udf.ignited = 1;
  }
  else {
    dQdt        = 30.0;
    eta_eff     = 0.3;
  }
  dTkdt  = dQdt*eta_eff/(udf.mk*cp)+1.0/(rho_b*cp)*udf.dPdt-dmkdt/udf.mk*(udf.Tk-Tb);
  udf.Tk = udf.Tk + dTkdt*dt;
  
  udf.mk += dmkdt*dt;
  Message0("	mk (new) %15.10e\n", udf.mk);
  
  udf.t += dt;
  free(Ymean); free(Ytab);
  for (i=0; i<udf.nTabSp; i++) free(tabSp[i]);
  free(tabSp);
  free(fp);
  free(fpS);
}
/* -------------------------------------------------------------------------- */

void UDFpremixedUpdateProgSrc(Domain *d)
{
  Thread *t;
  cell_t  c;
  int iT;
  double c_x[3], progMax, p_tab_ref, P, p_ratio, f1, f2, tmp;

  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;
    
    begin_c_loop_int(c,t) {
      /* set up the chemtable co-ordinates */
      c_x[0] = C_Z1MEAN(c,t);
      c_x[1] = C_Z1VAR(c,t);
      c_x[2] = C_PROG(c,t);
      
      /* rescale progress variable if pressure dependent */
      if (udf.chemTabPdep) {
        P      = ABS_P(C_P(c,t),op_pres);
        c_x[2] = MDUCchemTableProgScaleP(P, c_x, udf.chemTabPres);
      }

      /* lookup the progress variable source term */
      C_PROGSRC(c,t) = 0.0;
      progMax = MDUCchemTableProgMax(c_x, udf.chemTabPrem);
      if (C_PROG(c,t) < progMax)
        C_PROGSRC(c,t) = MDUClookupChemTable("SRC_PROG", c_x, NULL,
                                             udf.chemTabPrem);

      /* rescale source term if pressure dependent */
      if (udf.chemTabPdep) {
        /* get the pressure scaling ratio */
        p_tab_ref = MDUCgetChemTabPref(udf.chemTabPres);
        p_ratio   = log(P/p_tab_ref);

        /* look up the polynomial coefficients */ 
        c_x[2] = 1.0;
        f1 = MDUClookupChemTable("SRC_PROG", c_x, NULL, udf.chemTabPres);
        c_x[2] = 2.0;
        f2 = MDUClookupChemTable("SRC_PROG", c_x, NULL, udf.chemTabPres);

        if (ABS_D(f1) > 1.0e-12 || ABS_D(f2) > 1.0e-12) 
          C_PROGSRC(c,t) *= exp(f1*p_ratio*p_ratio + f2*p_ratio);
      }

    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
}
/* -------------------------------------------------------------------------- */

void UDFpremixedUpdateCFD(Domain *d)
{
  Thread *t;
  cell_t  c;

  int    i,j,k,iT,counter,nCount=500,lineCount=100,nZ1;
  int    nSp = MDUCgetNspecies(udf.mduc);
  double *Ymean, *Ytab;
  double ct_x[3], P, tmp;
  double phiMean, phiVar, progMax;

  /* for chemtable pressure dependence */
  double *Y_f1, *Y_f2;
  double cdot_f1, cdot_f2;
  double p_ratio, p_tab_ref;

  /* allocate local memory */
  Ymean = (double*) calloc(nSp,        sizeof(double));
  Ytab  = (double*) calloc(udf.nTabSp, sizeof(double)); 
  if (udf.chemTabPdep) {
    Y_f1 = (double*) calloc(udf.nTabSp, sizeof(double));
    Y_f2 = (double*) calloc(udf.nTabSp, sizeof(double));
  }

  /* for error reporting */
  char errFileName[1024];
  FILE *errFile;

  /* loop over reactive threads */
  counter = nZ1 = 0;
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;
    
    begin_c_loop_int(c,t) {
      phiMean = C_PHI1(c,t);
      phiVar  = C_PHIV1(c,t);
      /* set up the chemtable co-ordinates */
      /* TODO: generalize to enable specifying these,
         currently only 3D (Z,Zvar,C) */
      ct_x[0] = C_Z1MEAN(c,t);
      ct_x[1] = C_Z1VAR(c,t);
      ct_x[2] = C_PROG(c,t);

      if (udf.chemTabPdep) {
        P       = ABS_P(C_P(c,t),op_pres);
        ct_x[2] = MDUCchemTableProgScaleP(P, ct_x, udf.chemTabPres);
      }

      /* store the new source term */
      C_PROGSRC(c,t) = 0.0;
      progMax = MDUCchemTableProgMax(ct_x, udf.chemTabPrem);
      if (ct_x[2] < progMax) 
        C_PROGSRC(c,t)=MDUClookupChemTable("SRC_PROG", ct_x, NULL,
                                           udf.chemTabPrem);

      /* lookup the mean species from the table */
      MDUClookupChemTable("species", ct_x, Ytab, udf.chemTabPrem);

      /* lookup the species from the chemtable variables */
      for (i=0; i<nSp; i++) Ymean[i] = 0.0;
      for (i=0; i<udf.nTabSp; i++) {
        if (udf.iTabSp[i] >= 0)
          Ymean[udf.iTabSp[i]] = Ytab[i];
      }

      /* correct for pressure dependence if necessary */
      if (udf.chemTabPdep) {
        /* pressure ratio for scaling */
        p_tab_ref = MDUCgetChemTabPref(udf.chemTabPres);
        p_ratio   = log(P/p_tab_ref);

        /* look up the polynomial coefficients from the pressure chemtable
             log(phi/phi_ref) = f1*log(P/Pref)^2 + f2*log(P/Pref) + f3 */
        ct_x[2] = 1.0;
        MDUClookupChemTable("species", ct_x, Y_f1, udf.chemTabPres);
        cdot_f1 = MDUClookupChemTable("SRC_PROG", ct_x, NULL, udf.chemTabPres);
        ct_x[2] = 2.0;
        MDUClookupChemTable("species", ct_x, Y_f2, udf.chemTabPres);
        cdot_f2 = MDUClookupChemTable("SRC_PROG", ct_x, NULL, udf.chemTabPres);

        /* scale the species */
        for (i=0; i<udf.nTabSp; i++) {
          tmp = exp(Y_f1[i]*p_ratio*p_ratio + Y_f2[i]*p_ratio);
          if (ABS_D(Y_f1[i]) > 1.0e-12 || ABS_D(Y_f2[i]) > 1.0e-12)
            Ytab[i] *= tmp;
        }

        /* scale the progress variable source term */
        C_PSCALE(c,t) = tmp = exp(cdot_f1*p_ratio*p_ratio + cdot_f2*p_ratio);
        if (tmp > 0)
          C_PROGSRC(c,t) *= tmp;
      }

      /* ----- update the cell properties ----- */
      /* update the mixture fractions tracked by Fluent */
      for (i=0; i<udf.nY; i++) {
        if (udf.indexY[i] >= 0)
          C_YI(c,t,i) = Ymean[udf.indexY[i]];
        
        /* Error checking on Fluent species, make sure properly bounded
           NOTE: should consider enforcing sum(Y) = 1.0 */
        if ( isnan(C_YI(c,t,i)) || C_YI(c,t,i) > 1.0 ||
             ( (C_YI(c,t,i) < 0.0) && (ABS_D(C_YI(c,t,i)) > 1e-2) ) ) {
          sprintf(errFileName, "%s/mduc/Ymean_err_node%i_cell%i-%05i",
                  udf.caseDir, myid, c, N_TIME);
          Message0("Warning: writing Ymean_err: %s\n", errFileName);
          errFile = fopen(errFileName, "w");
          fprintf(errFile, "Problem with species %i:\n", i);
          fprintf(errFile, "Ymean    = %8.6f\n",  C_YI(c,t,i));
          fprintf(errFile, "phiMean  = %15.6e\n", phiMean);
          fprintf(errFile, "phiVar   = %15.6e\n", phiVar);
          fprintf(errFile, "phiMin   = %15.6e\n", udf.phiMin[1]);
          fprintf(errFile, "phiMax   = %15.6e\n", udf.phiMax[1]);
          fprintf(errFile, "deltaPhi = %15.6e\n", udf.deltaPhi[1]);
          fclose(errFile);
        }
      }
      
      /* update the cell temperature by iterating from the total
         enthalpy and mass fractions */
      /* TODO: option for looking up from table? */
      C_T(c,t) = MDUCcalcTfromH(nSp,C_HTOT(c,t),Ymean,1000.0,udf.mduc);
      /* C_T(c,t) = MDUClookupChemTable("T",ct_x,NULL,udf.chemTabPrem); */
      
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
        fprintf(errFile, "# progress    = %15.6e\n", C_PROG(c,t));
        fprintf(errFile, "# C_T(c,t)    = %15.6e K\n", C_T(c,t));
        fprintf(errFile, "# C_HTOT(c,t) = %15.6e J/kg\n", C_HTOT(c,t));
        fflush(errFile);
        fprintf(errFile, "# index Yi\n");
        for (i=0; i<nSp; i++)
          fprintf(errFile, "%8i %15.6e\n",i,Ymean[i]);

        fprintf(errFile, "\n\n");
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
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
  Message0("\n");

  tmp = 100.0/counter;
  Message0("\n Premixed Chemtable Update summary:\n");
  Message0(" ----- cells updated -----  number ----- percent -----\n");
  Message0("\t Z1        \t %8i \t %6.2f\n", nZ1, (double)(nZ1)*tmp);
  Message0(" -----------------------------------------------------\n");
  free(Ymean); free(Ytab);
  if (udf.chemTabPdep) {
    free(Y_f1); free(Y_f2);
  }
}
/* -------------------------------------------------------------------------- */

double calcKeckLaminarFlameSpeed(double phi, double Tu, double Pu)
{
  double C1,C2,C3;        /* constants */
  double gamma,beta;      /* exponents */
  double Tu_ref = 298;
  double Pu_ref = 101325;
  double sL_ref;          /* reference laminar flame speed */

  gamma =  2.18 - 0.80*(phi - 1.0);
  beta  = -0.16 + 0.22*(phi - 1.0);

  sL_ref = C1 + C2*(phi - C3)*(phi - C3);

  return sL_ref*pow(Tu/Tu_ref,gamma)*pow(Pu/Pu_ref,beta);
}
/* -------------------------------------------------------------------------- */

double calcNoTurbFlameSpeed(double sL, cell_t c, Thread *t)
{
  return sL;
}
/* -------------------------------------------------------------------------- */

double calcZimontTurbFlameSpeed(double sL, cell_t c, Thread *t)
{
  double A  = 0.52; /* default flame speed constant */
  double Cd = 0.37; /* default constant for turb. length scale */
  double uprime;    /* r.m.s. velocity (m/s) */
  double alpha;     /* thermal diffusivity (m^2/s) */
  double L_t;       /* turbulent length scale (m) */
  double tau_t;     /* turbulent time scale (s) */
  double tau_c;     /* chemical time scale (s) */

  /* compute the prereqs */
  uprime = pow(2.0*C_K(c,t)/3.0,0.5);
  alpha  = C_K_L(c,t)/(C_R(c,t)*C_CP(c,t));
  L_t    = Cd*uprime*uprime*uprime/(C_D(c,t)+1.e-6);

  /* compute the time scales */
  tau_t  = L_t/uprime;
  tau_c  = alpha/(sL*sL);

  return A*pow(tau_t/tau_c,0.25);
}
/* -------------------------------------------------------------------------- */

double calcPetersTurbFlameSpeed(double sL, cell_t c, Thread *t)
{
  double b1,b3,q,c_mu,C_s; /* constants */
  double Cd = 0.37;        /* default constant for turb. length scale */
  double a4;
  double uprime;           /* velocity scale (m/s) */
  double delta;            /* laminar flame thickness (m) */
  double L_t;              /* turbulent length scale (m) */
  double L_f;              /* flame brush thickness (m) */
  double L_alg;            /* algebraic flame brush thickness (m) */
  double L_star_q = 1.0;   /* Ewald's corrector */
  double tmp,sigma_t;

  /* set the default constants */
  b1   = 2.0;
  b3   = 1.0;
  q    = 0.66;
  c_mu = 0.09; /* from k-epsilon model */
  C_s  = 2.0;

  /* compute the prereqs */
  uprime = pow(2.0*C_K(c, t)/3.0,0.5);
  L_t    = Cd*uprime*uprime*uprime/(C_D(c, t)+1.e-6);
  delta  = (C_K_L(c,t)/C_CP(c,t))/(C_R(c,t)*sL);
  L_alg  = pow(2*c_mu/(C_s*Sc_t),0.5)*L_t;
  a4     = pow(3*c_mu*C_s/Sc_t,0.5);
  L_f    = pow(L_t*C_MU_T(c,t)/(uprime*C_R(c,t)*Sc_t),0.5);

  /* compute the function */
  tmp     = 0.25*a4*b3*b3*L_f/(b1*delta)*L_star_q;
  sigma_t = tmp+pow(tmp*tmp+0.5*C_s*b3*b3*uprime*L_f*L_f/(sL*delta*L_t),0.5);

  return sL*(1+sigma_t);
}
/* -------------------------------------------------------------------------- */

double calcPetersTurbFlameSpeedIgn(double sL, cell_t c, Thread *t)
{
  double b1,b3,q,c_mu,C_s; /* constants */
  double Cd = 0.37;        /* default constant for turb. length scale */
  double a4;
  double uprime;           /* velocity scale (m/s) */
  double delta;            /* laminar flame thickness (m) */
  double L_t;              /* turbulent length scale (m) */
  double L_f;              /* flame brush thickness (m) */
  double L_alg;            /* algebraic flame brush thickness (m) */
  double l_star_q = 1.0;   /* Ewald's corrector */
  double tmp,sigma_t;

  /* set the default constants */
  b1   = 2.0;
  b3   = 1.0;
  q    = 0.66;
  c_mu = 0.09; /* from k-epsilon model */
  C_s  = 2.0;

  /* compute the prereqs */
  uprime = pow(2.0*C_K(c, t)/3.0,0.5);
  L_t    = Cd*uprime*uprime*uprime/(C_D(c, t)+1.e-6);
  delta  = (C_K_L(c,t)/C_CP(c,t))/(C_R(c,t)*sL);
  L_alg  = pow(2*c_mu/(C_s*Sc_t),0.5)*L_t;
  a4     = pow(3*c_mu*C_s/Sc_t,0.5);
  L_f    = pow(L_t*C_MU_T(c,t)/(uprime*C_R(c,t)*Sc_t),0.5);

  /* compute the function */
  tmp     = 0.25*a4*b3*b3*L_f/(b1*delta)*l_star_q;
  
  l_star_q = (1-exp(-C_s*C_D(c,t)/C_K(c,t)*udf.t));
  sigma_t  = tmp*pow(l_star_q,1+q)
    + pow(tmp*tmp*pow(l_star_q,q+1)*pow(l_star_q,q+1)
          + l_star_q * 0.5*C_s*b3*b3*uprime*L_f*L_f/(sL*delta*L_t),0.5);

  return sL*(1+sigma_t);
}
/* -------------------------------------------------------------------------- */
