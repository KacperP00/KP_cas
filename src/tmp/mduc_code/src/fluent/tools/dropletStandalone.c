#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
/* #include "sdmDrop.h" */

double calcHvapHeptane(double T);
double antoinePvap(double T);
void getThermo(char *fileName, double *hfg, double *cpgas, double *cpliq, double *Psat, double T);

/* int readThermoNIST(char *fileName, double **hfg, double **cpgas, double **cpliq, */
/*                    double **Psat, double **T); */

/* void interpThermoNIST(double *hfg, double *cpG, double *cpL, double *Psat, */
/*                       double T, int nNIST, double *NIST_hfg, double *NIST_cpG, */
/*                       double *NIST_cpL, double *NIST_Psat, double *NIST_T); */

/* void printThermoNIST(int nNIST, double *hfg, double *cpG, double *cpL, */
/*                      double *Psat, double *T); */

/* void initSDMDrop(sdmStruct **drop) */
/* { */
/*   int j; */
/*   double dZ; */

/*   drop[0] = (sdmStruct*) calloc(1, sizeof(sdmStruct)); */

/*   sdmStruct *d = drop[0]; */

/*   /\* get the data from the specified file *\/ */
/* /\*   strcpy(d->NIST_File, "/home/emdoran/M47/libudf_mduc/src/heptaneSat.dat"); *\/ */
/*   strcpy(d->NIST_File, "heptaneSat.dat"); */

/*   fprintf(stderr, "d->NIST_file = %s\n", d->NIST_File); */
/*   drop[0]->nNIST = readThermoNIST(drop[0]->NIST_File, &drop[0]->hfg, &drop[0]->cpG, &drop[0]->cpL, &drop[0]->Psat, &drop[0]->T); */
 
/*   printThermoNIST(d->nNIST, d->hfg, d->cpG, d->cpL, d->Psat, d->T); */

/*   d->T0    = 300;      /\* initial droplet temperature [K] *\/ */
/*   d->Mfuel = 100.204; /\* molecular weight of the fuel [g/mol] *\/ */
/*   d->Mair  = 28.9;    /\* molecular weight of the air [g/mol] *\/ */

/*   /\* initialise the SDM mixture fraction grid *\/ */
/*   d->nZ = 101; */
/*   d->Z   = (double*) calloc(d->nZ, sizeof(double)); */
/*   d->pdf = (double*) calloc(d->nZ, sizeof(double)); */

/*   dZ = 1.0/(d->nZ-1); */
/*   d->Z[0] = 0.0; */
/*   for (j=1; j<d->nZ; j++) d->Z[j] = d->Z[j-1] + dZ; */
/* } */
/* /\* ------------------------------------------------------------------------- *\/ */

int main()
{
  int i, j, k;

  double P_inf, Psat, PsatCC, hgCC;
  double Xfs, *Yfs, Yf_inf;
  double Mfuel, Mair;
  double Ts, T0, T1, Tinf, Tbp, Tavg, dT;
  double alpha;

  alpha = 2/3;

  double hfg, cpG, cpL, cpAvg;

  double Q, *BT, *Bf;

  double *Taxis;

  int nT = 100;

  double R, Ru = 8.314;

  /* sdmStruct *sdmDrop = NULL; */

  /* initSDMDrop(&sdmDrop); */

  /* interpThermoNIST(&hfg, &cpG, &cpL, &Psat, 438.0, sdmDrop->nNIST, sdmDrop->hfg, */
/*                    sdmDrop->cpG, sdmDrop->cpL, sdmDrop->Psat, sdmDrop->T); */

  // for heptane
  Tbp = 371.53; // 447.1;
  cpL = 225.14; // [J/mol-K]

  // set range of surface temperatures
  T0 = 300;
  T1 = 500;

  // set ambient properties
  Tinf = 500;
  P_inf = 1.e5;
  Yf_inf = 0.0;

  Mfuel = 100;
  Mair = 28.9;

  R = Ru/(Mfuel/1000);

  Taxis = calloc(nT, sizeof(double));
  Yfs = calloc(nT, sizeof(double));
  BT = calloc(nT, sizeof(double));
  Bf = calloc(nT, sizeof(double));

  dT = (535 - T0)/(nT-1);

  Taxis[0] = T0;
  for (i=1; i<nT; i++) Taxis[i] = Taxis[i-1] + dT;

  for (i=0; i<nT; i++)
    {
      //Tavg = alpha*Taxis[i] + (1-alpha)*(Tinf-Taxis[i]);
      Tavg = Taxis[i];
/*       getThermo("heptaneSat.dat", &hfg, &cpG, &cpL, &Psat, Tavg); */
      getThermo("nOctane.dat", &hfg, &cpG, &cpL, &Psat, Tavg);

      //hgCC = calcHvapHeptane(Tavg);
      //PsatCC = 1.0e5*exp(hfg/Ru*(1.0/Tavg - 1.0/Tbp));
      //PsatCC = antoinePvap(Tavg);
                        
      Q = hfg + cpL*(Tavg-T0);

      Xfs = Psat/P_inf;
      Yfs[i] = Mfuel*Xfs/(Mfuel*Xfs+(1-Xfs)*Mair);

      BT[i] = cpG*(Tinf-Tavg)/Q;
      Bf[i] = (Yfs[i]-Yf_inf)/(1-Yfs[i]);

      //printf("T = %8.2f, hfg = %15.6e, cpG = %15.6e, cpL = %15.6e\n", Taxis[i], hfg, cpG, cpL);
      printf("T = %6.1f, Psat = %15.6e, CC = %15.6e, diff = %8.6f\n",
             Tavg, Psat, PsatCC, fabs(Psat-PsatCC)/Psat*100);
      //printf("T = %6.1f, hfg = %15.6e, hgCC = %15.6e, diff = %8.6f\n",
      // Tavg, hfg, hgCC, fabs(hgCC-hfg)/hfg*100);
    }

  FILE *pFile = fopen("transferNumber.dat", "w");
        
  fprintf(pFile, "# Ts Yfs BT Bf\n");

  for (i=0; i<nT; i++)
    fprintf(pFile, "%8.2f %15.6e %15.6e %15.6e\n", Taxis[i], Yfs[i], BT[i], Bf[i]);

  fclose(pFile);

  return 0;
}

/* from NIST Chemistry Webbook for n-Heptane
   referenced from Majer and Svoboda, 1985
   valid in temperature range from 298 - 363 Kelvin
*/
/* double calcHvapHeptane(double T) */
/* { */
/*   double Tc   = 540.2;   // critical temperature [K] */
/*   double A    = 53.66e3; // constant [J/mol] */
/*   double beta = 0.2831; */
/*   double Tr; // reduced temperature */

/*   Tr = T/Tc; */
/*   return A*exp(-beta*Tr)*pow(1-Tr,beta); */
/* } */

/* double antoinePvap(double T) */
/* { */
/*   double A, B, C; */
/*   A = 4.02832; */
/*   B = 1268.636; */
/*   C = -56.199; */

/*   return pow(A - (B/(T+C)), 10)*1.0e5; */
/* } */

void getThermo(char *fileName, double *hfg, double *cpgas, double *cpliq, double *Psat, double T)
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
  char comment[800];
  float rTmp;
  double line[25], line_m1[25];
  fpos_t pos;

  // get the header line
  fgets(comment, 800, pFile);

  // read the first line
  for (i=0; i<nVars; i++)
    {
      fscanf(pFile, "%f", &rTmp);
      line_m1[i] = (double)(rTmp);
    }
  if (T == line_m1[0])
    {
      *Psat  = line_m1[sPsat]*1e6;
      *cpliq = line_m1[sCpL];
      *cpgas = line_m1[sCpG];
      *hfg   = (line_m1[sHG]-line_m1[sHL])*1e3;
      return;
    }

  while (!feof(pFile))
    {
      fgetpos(pFile, &pos);
      for (i=0; i<25; i++)
        {
          fscanf(pFile, "%f", &rTmp);
          line[i] = (double)(rTmp);
        }

      if (T == line[0])
        {
          *Psat  = line[sPsat]*1e6;
          *cpliq = line[sCpL];
          *cpgas = line[sCpG];
          *hfg   = (line[sHG]-line[sHL])*1e3;
          return;
        }
      else if (T < line[0] && T > line_m1[0])
        {
          dT     = line[sT]-line_m1[sT];
          *Psat  = line_m1[sPsat] + (line[sPsat]-line_m1[sPsat])/dT;
          *cpliq = line_m1[sCpL]  + (line[sCpL]-line_m1[sCpL])/dT;
          *cpgas = line_m1[sCpG]  + (line[sCpG]-line_m1[sCpG])/dT;
          HG     = line_m1[sHG]   + (line[sHG]-line_m1[sHG])/dT;
          HL     = line_m1[sHL]   + (line[sHL]-line_m1[sHL])/dT;
          *Psat  = *Psat*1e6;
          *hfg   = (HG - HL)*1e3;
          return;
        }

      for (i=0; i<25; i++) line_m1[i] = line[i];
    }

  printf("ERROR: temperature range not sufficient for T = %8.2f\n", T);
  fclose(pFile);
}

void interpThermoNIST(double *hfg, double *cpG, double *cpL, double *Psat,
                      double T, int nNIST, double *NIST_hfg, double *NIST_cpG,
                      double *NIST_cpL, double *NIST_Psat, double *NIST_T)
{
  int jR, jL;
  double C;

  if (T < 300.0 || T > 538.80)
    {
      fprintf(stderr, "Temperature out of range\n");
      T = 447.1;
    }

  // find the coarse mesh points either side of current fine mesh
  jR = 0;
  do { jR++; } while ( NIST_T[jR] < T && jR < nNIST );
  jL = jR - 1; // assume for within 1 grid point

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
  fprintf(stderr, "Reading NIST thermo file: %s\n", fileName);
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

  if (!pFile)
    {
      fprintf(stderr, "Error: thermo file %s not found in getThermoNIST\n",
              fileName);
      exit(1);
    }

  len = 0;
  // get the number of data lines
  while (!feof(pFile))
    {
      strcpy(comment, "\n");
      fgets(comment, 1000, pFile);
      if (strcmp(comment, "\n"))
        len++;
    }
  len -= 1;
  rewind(pFile);

  fprintf(stderr, "file length = %i\n", len);

  /* allocate memory */
  T[0]     = (double*) calloc(len, sizeof(double));
  Psat[0]  = (double*) calloc(len, sizeof(double));
  cpgas[0] = (double*) calloc(len, sizeof(double));
  cpliq[0] = (double*) calloc(len, sizeof(double));
  hfg[0]   = (double*) calloc(len, sizeof(double));

  // get the header line
  fgets(comment, 800, pFile);

  for (j=0; j<len; j++)
    {
      /* scan in the line */
      for (i=0; i<25; i++)
        {
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

  fprintf(stderr, "# T[K] Psat[Pa] cp-vapour[J/mol-K] cp-liquid[J/mol-K] hfg\n");
  fprintf(stderr, "n = %i\n", nNIST);
  for (i=0; i<nNIST; i++)
    fprintf(stderr, "%6.2f %15.6e %15.6e %15.6e %15.6e\n", T[i], Psat[i],
              cpG[i], cpL[i], hfg[i]);
}
/* ------------------------------------------------------------------------- */

