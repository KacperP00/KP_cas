#include "sdmDrop.h"
#include "udfTools.h"
#include "accessMDUC.h"

extern udfStruct udf;
/* quantities for computing the single droplet model (SDM) */
static sdmStruct *sdm = NULL;

void SDMinit(void)
{
  Message0("\n ----- Initialising SDM droplet -----\n");
  int    i, j;
  double dZ;

  /* initialize the data structure */
  sdm = (sdmStruct*) calloc(1, sizeof(sdmStruct));
  sdmStruct *D = sdm;

  /* set the output file name */
  if (RP_is_Defined("mduc/directory")) {
    strcpy(sdm->monitor, RP_Get_String("mduc/directory"));
    strcat(sdm->monitor, "/sdm_monitor.dat");
  }
  else
    strcpy(sdm->monitor, "sdm_monitor.dat");

  /* get the data from the specified file */
  if (RP_is_Defined("mduc/sdm/nist-sat-data"))
    strcpy(D->NIST_File, RP_Get_String("mduc/sdm/nist-sat-data"));
  else Error("Need NIST saturation data for SDM model");

  D->xiMax = 8;
  if (RP_is_Defined("mduc/sdm/xi-max"))
    D->xiMax = RP_Get_Real("mduc/sdm/xi-max");

  Message0("Reading NIST thermo file: %s\n", D->NIST_File);
  D->nNIST = readThermoNIST(D->NIST_File, &D->hfg, &D->cpG, &D->cpL, &D->Psat,
                            &D->T);

  /* printThermoNIST(D->nNIST, D->hfg, D->cpG, D->cpL, D->Psat, D->T); */

  D->T0    = 300;     /* initial droplet temperature [K] */
  D->Mfuel = 100.204; /* molecular weight of the fuel [g/mol] */
  D->Mair  = 28.9;    /* molecular weight of the air [g/mol] */

  /* initialise the SDM mixture fraction grid */
  D->nZ  = 201;
  D->Z   = (double*) calloc(D->nZ, sizeof(double));
  D->pdf = (double*) calloc(D->nZ, sizeof(double));
  D->Wv  = (double*) calloc(D->nZ, sizeof(double));

  /* initialise constants */
  D->C_x         = 2.0;
  D->one_third   = 1.0/3.0;
  D->six_over_pi = 6.0/M_PI;

  dZ = 1.0/(D->nZ-1);
  D->Z[0] = 0.0;
  for (j=1; j<D->nZ; j++) D->Z[j] = D->Z[j-1] + dZ;

  /* make sure the diffusion coefficient is specified to 1.0 */
#if !RP_HOST
  Domain *d = Get_Domain(1);
  Thread *t;
  cell_t c;

  /* ensure diffusion coefficients are at least 1 */
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      for (i=0; i<3; i++) {
        if (udf.iFDT[i] >= 0)
          C_UDMI(c,t,udmOffset+udf.iFDT[i]) = 1.0;
      }
    } end_c_loop_int(c,t);
  } /* end loop over threads */
#endif /* !RP_HOST */

  Message0("\n ------------------------------------\n");
}
/* ------------------------------------------------------------------------- */

int checkSDM(Domain *d)
{
#if !RP_HOST
  Thread *t;
  cell_t c;

  int    mix[3] = {0, 0, 0};
  double TINY = 1e-6;

  /* ----- check which mixture fractions are present ----- */
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      if (udf.iZ1 >= 0)
        if ( C_Z1MEAN(c,t) > TINY ) mix[1] = 1;
      if (udf.iZ2 >= 0)
        if ( C_Z2MEAN(c,t) > TINY ) mix[2] = 1;
      if (mix[1] + mix[2] == 2) break;
    } end_c_loop_int(c,t); /* end loop over all cells */
  } /* end loop over cell threads */

#if RP_NODE /* parallel */
  mix[1] = PRF_GIHIGH1(mix[1]);
  mix[2] = PRF_GIHIGH1(mix[2]);
#endif /* RP_NODE */
  
  return mix[1] + mix[2];
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

void SDMcalcSrc(Domain *d, int nInj, int nNoz, Injection ***stream)
{
#if !RP_HOST
  Thread   *t;
  Particle *p;
  cell_t    c;

  int    i, j;
  int    calcSpald;
  double mass;         /* cell mass */
  double tau;          /* turbulent eddy turnover time [s] */
  double B;            /* Spalding number (mass transfer of droplet) */
  double Psat;         /* saturation pressure [Pa] */
  double hfg;          /* enthalpy of vaporisation for droplet */
  double cpG, cpL;     /* specific heat of gas/liquid */
  double Xfs, Yfs;     /* fuel mole/mass fraction at droplet surface */
  double Ts;           /* temperature at droplet surface */
  double Re_t;         /* turbulent cell Reynolds number */
  double Lo, Hi;       /* for root solve */
  double Tmin, Tmax;
  double tmp, TINY = 1.0e-6;

  /* TODO: !!!!! this needs to be generalized later !!!!! */
  double nu_f = 3.6723703644565089;

  /* loop over all injections/particles and compute the volume of each fuel in
     each cell */
  for (i=0; i<nInj; i++) {
    for (j=0; j<nNoz; j++) {
      loop(p, stream[i][j]->p) {
        /* Get the cell and thread that the particle is currently in */
        c = P_CELL(p);   
        t = P_CELL_THREAD(p);
        
        if (i == 0)
          C_FUEL1_VOL(c,t) += P_MASS(p)/P_RHO(p)*p->number_in_parcel;
        else if (i == 1)
          C_FUEL2_VOL(c,t) += P_MASS(p)/P_RHO(p)*p->number_in_parcel;
      } /* end loop over particles */
    } /* end loop over nozzles */
  } /* end loop over injections */
  
  /* compute the source term in all cells */
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      /* get the cell properties */
      mass = C_R(c,t)*C_VOLUME(c,t);
      tau  = C_K(c,t)/(C_D(c,t)+1.0e-6);
      
      /* set the dissipation coefficient to the standard just in case */
      for (i=0; i<3; i++) {
        if (udf.iZV_EPS[i] >= 0) /* variance dissipation */
          C_UDMI(c,t,udmOffset+udf.iZV_EPS[i]) = sdm->C_x/tau;
      }
  
      /*   Re_turb = rho*k^2/(epsilon*mu) */
      Re_t = C_R(c,t)*C_K(c,t)*tau/C_MU_L(c,t);
#if DEBUG_SDM
      C_RE_TURB(c,t) = Re_t;
#endif DEBUG_SDM

      /* check whether Spalding number needs to be computed */
      calcSpald = 0;
      if (udf.iFuelV[1] >= 0) {
        if (C_FUEL1_VOL(c,t) > 0.0 && C_FUEL1_SRC(c,t) > 0.0)
          calcSpald = 1;
      }
      if (udf.iFuelV[2] >= 0) {
        if (C_FUEL2_VOL(c,t) > 0.0 && C_FUEL2_SRC(c,t) > 0.0)
          calcSpald = 1;
      }

      /* Spalding number is the same if using the same fuel in each
         injection, only need to compute once */
      if (calcSpald == 1) {
        /* assign ambient conditions */
        sdm->Yf_inf = C_YI(c,t,0);
        sdm->T_inf  = C_T(c,t);
        sdm->P      = ABS_P(C_P(c,t), op_pres);
        
        /* compute brackets to see whether solution is possible */
        Tmin = 300; /* assume minimum liquid temperature */
        Tmax = sdm->T[sdm->nNIST - 1];
        Lo   = calcSpaldingFunc(Tmin, sdm);
        Hi   = calcSpaldingFunc(Tmax, sdm);
        
        if (Lo*Hi < 0.0) /* root exists */
          Ts = rootBrent(calcSpaldingFunc, Tmin, Tmax, 1.0e-4, sdm);
        else /* assume at boiling limit */
          Ts = Tmax;
        
        /* get the properties at the surface temperature and compute
           the Spalding number */
        interpThermoNIST(&hfg, &cpG, &cpL, &Psat, Ts, sdm->nNIST,
                         sdm->hfg, sdm->cpG, sdm->cpL, sdm->Psat, sdm->T);
        sdm->Q = hfg + cpL*(Ts - sdm->T0);
        B = cpG*(sdm->T_inf - Ts)/sdm->Q;

        /* compute the saturation mixture fraction */
        /* TODO: !!!!! needs to be generalised !!!!! */
        Yfs       = (B + C_YI(c,t,0))/(1.0 + B);
        sdm->Zsat = (nu_f*Yfs - C_YI(c,t,2) + 0.2153)/(nu_f + 0.2153);
        sdm->Zsat = MIN_D(1.0, MAX_D(0.0, sdm->Zsat));
#if DEBUG_SDM
        C_SPALD_NO(c,t) = B;
#endif /* DEBUG_SDM */
      } /* end Spalding number */

      /* compute the source terms */
      if (udf.iZ1 >= 0) {
        if (C_FUEL1_VOL(c,t) > 0.0 && C_FUEL1_SRC(c,t) > 0.0)
          SDMcalcVarsrc(C_FUEL1_VOL(c,t), C_FUEL1_SRC(c,t), mass, tau, B,
                        C_MU_L(c,t), C_Z1MEAN(c,t), C_Z1VAR(c,t),
                        &C_FDT_1(c,t), &C_Z1VAR_SRC(c,t), &C_Z1VAR_EPS(c,t));
      }
       
      if (udf.iZ2 >= 0) {
        if (C_FUEL2_VOL(c,t) > 0.0 && C_FUEL2_SRC(c,t) > 0.0)
          SDMcalcVarsrc(C_FUEL2_VOL(c,t), C_FUEL2_SRC(c,t), mass, tau, B,
                        C_MU_L(c,t), C_Z2MEAN(c,t), C_Z2VAR(c,t),
                        &C_FDT_2(c,t), &C_Z2VAR_SRC(c,t), &C_Z2VAR_EPS(c,t));
      }
      
      if (udf.iZ1 >= 0 && udf.iZ2 >= 0 && udf.iZ0V >= 0) {
        if ((C_FUEL1_VOL(c,t) > 0.0 && C_FUEL1_SRC(c,t) > 0.0) ||
            (C_FUEL2_VOL(c,t) > 0.0 && C_FUEL2_SRC(c,t) > 0.0) )
          SDMcalcVarsrc(C_FUEL1_VOL(c,t) + C_FUEL2_VOL(c,t),
                        C_FUEL1_SRC(c,t) + C_FUEL2_SRC(c,t), mass, tau, B,
                        C_MU_L(c,t), C_Z1MEAN(c,t)+C_Z2MEAN(c,t), C_Z0VAR(c,t),
                        &C_FDT_0(c,t), &C_Z0VAR_SRC(c,t), &C_Z0VAR_EPS(c,t));
      }

    } end_c_loop_int(c,t);
  } /* end loop over thread */
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

double SDMcalcVarsrc(double fVol, double fSrc, double mass, double tau,
                     double B, double mu_L, double Zmean, double Zvar,
                     double *fdt, double *vSrc, double *vEps)
{
  int    i,j, brac;

  double tau_tauZ;         /* time scale ratio (turbulent/scalar) */
  double C;                /* time scale constant for Hasse formulation */
  double Sc_t = 0.7;       /* turbulent Schmidt number */
  double dropD;            /* effective droplet diameter for SDM model */
  double *Wv  = sdm->Wv;   /* mean fuel source conditioned on mixture fraction*/
  double *pdf = sdm->pdf;  /* local cell PDF */
  double Zsat = sdm->Zsat; /* saturation mixture fraction */
  int pdfBound[2];
  
  double srcP, srcM;
  double ab, xi;
  double Lo, Hi; /* for bracketing */
  double tmp;

  /* get the mixture fraction grid (initialised in struct) */
  int nZ      = sdm->nZ;
  double *eta = sdm->Z;
  double C_x  = sdm->C_x;

  /* compute the effective diameter of a single droplet */
  dropD = pow(sdm->six_over_pi*fVol, sdm->one_third);

  /* the change in mass fraction from this time step */
  sdm->Zconst = fSrc/mass*CURRENT_TIMESTEP;

  /* compute the constant part of alpha */
  sdm->alpha = 2.0*M_PI*mu_L/(mass*sdm->Re_t*Sc_t)*log(1.+B)*dropD;

  /* set the source term constraint */
  sdm->W = fSrc/mass; /* [-/s] */

  /* compute the scaled pdf to be used for integration */
  MDUCcalcPDF_1D("beta", nZ, eta, Zmean/Zsat, Zvar/(Zsat*Zsat), pdf, pdfBound);

  /* solve for the exponent */
  Lo = 1.0; Hi = 10.0;
  brac = rootBrac(SDMmeanSrcFunc, &Lo, &Hi, sdm);
  xi = 0.0;
  if (brac == 1) /* solution exists */
    xi = rootBisect(SDMmeanSrcFunc, Lo, Hi, 1.0e-4, sdm);

  ab = 0.0;
  if (xi >= 1.0 && xi < sdm->xiMax)
    ab = sdm->alpha*pow(sdm->Zconst, -xi);

  /* compute the source terms for the variance equation based on 
     the linearised SDM model from Hasse's thesis */
  tau_tauZ = C_x;
  C = 0.0;
  if (Zmean > 1.0e-6 && ab > 0.0) {
    /* compute the constant for the time scale */
    C = ab*( xi*pow(Zmean, xi-1.0) - (xi+1.0)*pow(Zmean, xi) );
    if (C*tau > 0.0) {
      tau_tauZ = 3.*C_x*C*tau/(1.-exp(-3.*C*tau));
      (*fdt)   = tau_tauZ/(tau_tauZ - 2.*C*tau);
    }
  }
  
#ifdef MDUC_SPRAY_LINEAR
  if (C > 0) {
    (*fdt)  = tau_tauZ/(tau_tauZ - 2.*C*tau);
    (*vSrc) = 2*C;
  }
#else
  /* compute the conditional mean source of fuel */
  for (j=0; j<nZ; j++) Wv[j] = 0;
  if (sdm->Zconst > 0.0 && xi >= 1.0 && xi < sdm->xiMax) {
    for (j=1; j<nZ; j++)
      Wv[j] = sdm->alpha*pow(eta[j]*Zsat/sdm->Zconst, xi);
  }

  /* compute the source terms for the variance equation based on 
     the SDM model */
  srcP = srcM = 0.0;
  for (j=pdfBound[0]; j<=pdfBound[1]; j++) {
    tmp   = (eta[j]*Zsat - Zmean)*Wv[j]*pdf[j];
    srcP += tmp;
    srcM += (eta[j]*Zsat - Zmean)*tmp;
  }
  srcP *= Zsat;
  srcM *= Zsat;
          
  (*vSrc) = 2*(1.0 - Zmean)*srcP - srcM;
#endif /* MDUC_SPRAY_LINEAR */

  /* compute the dissipation term and ensure fdt is bounded */
  (*vEps) = MAX_D(tau_tauZ/tau, C_x/tau);
  (*fdt)  = MIN_D(1.5, MAX_D(1.0, (*fdt)));
}
/* -------------------------------------------------------------------------- */

void getInjection(char *name, int nNozzle, Injection ***I)
{
  Message0("Getting injection: %s\n", name);
  int i;
  char ident[80];
  I[0] = (Injection**) calloc(nNozzle, sizeof(Injection*));

  for (i=0; i<nNozzle; i++) {
    sprintf(ident, "%s-%i", name, i+1);
    I[0][i] = Pick_Injection(ident);
  }
}
/* ------------------------------------------------------------------------- */

double calcSpaldingFunc(double Ts, void *args)
{
  sdmStruct *vars = (sdmStruct*) args;

  double BT, Bf, Psat, Q, cpG, cpL, hfg, Xfs, Yfs;

  /* get the thermo tabulated from NIST */
  interpThermoNIST(&hfg, &cpG, &cpL, &Psat, Ts, vars->nNIST, vars->hfg,
                   vars->cpG, vars->cpL, vars->Psat, vars->T);

  /* compute mass fraction at the surface based on saturation pressure*/
  Xfs = Psat/vars->P;
  Yfs = vars->Mfuel*Xfs/(vars->Mfuel*Xfs+(1-Xfs)*vars->Mair);

  Q = hfg + cpL*(Ts - vars->T0);

  BT = cpG*(vars->T_inf - Ts)/Q;
  Bf = (Yfs - vars->Yf_inf)/(1.0 - Yfs);

  return BT - Bf;
}
/* -------------------------------------------------------------------------- */

double SDMmeanSrcFunc(double xi, void *args)
{
  sdmStruct *vars = (sdmStruct*) args;

  int j;
  double W = 0.0;

  for (j=1; j<vars->nZ; j++)
    W += vars->alpha*pow(vars->Z[j]*vars->Zsat/vars->Zconst, xi)*vars->pdf[j];
  
  return W - vars->W;
}
/* -------------------------------------------------------------------------- */

/* from NIST Chemistry Webbook for n-Heptane
   referenced from Majer and Svoboda, 1985
   valid in temperature range from 298 - 363 Kelvin
*/
double calcHvapHeptane(double T)
{
  double Tc   = 540.2;   /* critical temperature [K] */
  double A    = 53.66e3; /* constant [J/mol] */
  double beta = 0.2831;
  double Tr; /* reduced temperature */

  Tr = T/Tc;
  return A*exp(-beta*Tr)*pow(1-Tr,beta);
}
/* ------------------------------------------------------------------------- */

void getThermoNIST(char *fileName, double *hfg, double *cpgas, double *cpliq,
                   double *Psat, double T)
{
  int i;
  double HL, HG, dT;
  int nVars, sT, sPsat, sCpL, sCpG, sHG, sHL;

  nVars = 25;
  sT    = 0;
  sPsat = 1;
  sHL   = 5;
  sCpL  = 8;
  sHG   = 17;
  sCpG  = 20;

  FILE *pFile = fopen(fileName, "r");

  if (!pFile) {
    fprintf(stderr, "Error: thermo file %s not found in getThermoNIST\n",
            fileName);
    exit(1);
  }

  char comment[800];
  float rTmp;
  double line[25], line_m1[25];
  fpos_t pos;

  /* get the header line */
  fgets(comment, 800, pFile);

  /* read the first line */
  for (i=0; i<nVars; i++) {
    fscanf(pFile, "%f", &rTmp);
    line_m1[i] = (double)(rTmp);
  }
  if (T == line_m1[0]) {
    *Psat  = line_m1[sPsat]*1e6;
    *cpliq = line_m1[sCpL];
    *cpgas = line_m1[sCpG];
    *hfg   = (line_m1[sHG]-line_m1[sHL])*1e3;
    fclose(pFile);
    return;
  }

  while (!feof(pFile)) {
    fgetpos(pFile, &pos);
    for (i=0; i<25; i++) {
      fscanf(pFile, "%f", &rTmp);
      line[i] = (double)(rTmp);
    }
    
    if (T == line[0]) {
      *Psat  = line[sPsat]*1e6;
      *cpliq = line[sCpL];
      *cpgas = line[sCpG];
      *hfg   = (line[sHG]-line[sHL])*1e3;
      fclose(pFile);
      return;
    }
    else if (T < line[0] && T > line_m1[0]) {
      dT     = line[sT]-line_m1[sT];
      *Psat  = line_m1[sPsat] + (line[sPsat]-line_m1[sPsat])/dT;
      *cpliq = line_m1[sCpL]  + (line[sCpL]-line_m1[sCpL])/dT;
      *cpgas = line_m1[sCpG]  + (line[sCpG]-line_m1[sCpG])/dT;
      HG     = line_m1[sHG]   + (line[sHG]-line_m1[sHG])/dT;
      HL     = line_m1[sHL]   + (line[sHL]-line_m1[sHL])/dT;
      *Psat  = *Psat*1e6;
      *hfg   = (HG - HL)*1e3;
      fclose(pFile);
      return;
    }
    
    for (i=0; i<25; i++) line_m1[i] = line[i];
  }
  
  fprintf(stderr, "ERROR: temperature range not sufficient for T = %8.2f\n", T);
  fclose(pFile);
}
/* ------------------------------------------------------------------------- */

void interpThermoNIST(double *hfg, double *cpG, double *cpL, double *Psat,
                      double T, int nNIST, double *NIST_hfg, double *NIST_cpG,
                      double *NIST_cpL, double *NIST_Psat, double *NIST_T)
{
  int jR, jL;
  double C;

  /* ensure that the solution is bounded */
  T = MAX_D(NIST_T[0], MIN_D(NIST_T[nNIST-1],T));

  /* find the coarse mesh points either side of current fine mesh */
  jR = 0;
  do { jR++; } while ( NIST_T[jR] < T && jR < nNIST );
  jL = jR - 1; /* assume for within 1 grid point */

  /* scaling constant */
  C     = (T - NIST_T[jL])/(NIST_T[jR] - NIST_T[jL]);

  *hfg  = NIST_hfg[jL]  + C*(NIST_hfg[jR]  - NIST_hfg[jL]);
  *cpG  = NIST_cpG[jL]  + C*(NIST_cpG[jR]  - NIST_cpG[jL]);
  *cpL  = NIST_cpL[jL]  + C*(NIST_cpL[jR]  - NIST_cpL[jL]);
  *Psat = NIST_Psat[jL] + C*(NIST_Psat[jR] - NIST_Psat[jL]);
}
/* ------------------------------------------------------------------------- */

int readThermoNIST(char *fileName, double **hfg, double **cpgas, double **cpliq,
                   double **Psat, double **T)
{
  int i, j, len;
  double HL, HG, dT;
  int nVars, sT, sPsat, sCpL, sCpG, sHG, sHL;

  char comment[800];
  float rTmp;
  double line[25];

  nVars = 25;
  sT    = 0;
  sPsat = 1;
  sHL   = 5;
  sCpL  = 8;
  sHG   = 17;
  sCpG  = 20;

  FILE *pFile = fopen(fileName, "r");

  if (!pFile) {
    fprintf(stderr, "Error: thermo file %s not found in getThermoNIST\n",
            fileName);
    exit(1);
  }

  len = 0;
  /* get the number of data lines */
  while (!feof(pFile)) {
    strcpy(comment, "\n");
    fgets(comment, 1000, pFile);
    if (strcmp(comment, "\n"))
      len++;
  }
  len -= 1;
  rewind(pFile);

  /* fprintf(stderr, "file length = %i\n", len); */

  /* allocate memory */
  T[0]     = (double*) calloc(len, sizeof(double));
  Psat[0]  = (double*) calloc(len, sizeof(double));
  cpgas[0] = (double*) calloc(len, sizeof(double));
  cpliq[0] = (double*) calloc(len, sizeof(double));
  hfg[0]   = (double*) calloc(len, sizeof(double));

  /* get the header line */
  fgets(comment, 800, pFile);

  for (j=0; j<len; j++) {
    /* scan in the line */
    for (i=0; i<25; i++) {
      fscanf(pFile, "%f", &rTmp);
      line[i] = (double)(rTmp);
    }
    
    /* store relevant data */
    T[0][j]     = line[sT];
    Psat[0][j]  = line[sPsat]*1.0e6;
    cpgas[0][j] = line[sCpG];
    cpliq[0][j] = line[sCpL];
    hfg[0][j]   = (line[sHG]-line[sHL])*1.0e3;
  }
  
  fclose(pFile);
  return len;
}
/* ------------------------------------------------------------------------- */

void printThermoNIST(int nNIST, double *hfg, double *cpG, double *cpL,
                     double *Psat, double *T)
{
  int i;

  fprintf(stderr,"# T[K] Psat[Pa] cp-vapour[J/mol-K] cp-liquid[J/mol-K] hfg\n");
  fprintf(stderr, "n = %i\n", nNIST);
  for (i=0; i<nNIST; i++)
    fprintf(stderr, "%6.2f %15.6e %15.6e %15.6e %15.6e\n", T[i], Psat[i],
              cpG[i], cpL[i], hfg[i]);
}
/* ------------------------------------------------------------------------- */
