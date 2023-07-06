#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "chemMech.h"
#include "chemUtils.h"
#include "accessChem.h"
#include "solveUtils.h"

void ComputeThermoDataComp(double *h, double *cp, double T, chemStruct *chem)
{
  /* just call the compiled mechanism function */
  ComputeThermoData(h, cp, T);
}
/* -------------------------------------------------------------------------- */

void ComputeThermoDataInterp(double *h, double *cp, double T, chemStruct *chem)
{
  int i,j,iT;
  double *a, *b, R = 8314.34, R_W;
  double T_low, T_mid, T_high;

  for (i=0; i<chem->nSpecies; i++) {
    a = b = NULL;
    /* ----- determine the correct temperature range ----- */
    T_low  = chem->specs[i].Trange[0];
    T_mid  = chem->specs[i].Trange[1];
    T_high = chem->specs[i].Trange[2];
    if (T < T_high && T > T_mid)
      iT = 1;
    else if (T >= T_low)
      iT = 0;
    else {
      fprintf(stderr,"Error: temperature (%f) is out of range for species %s\n",
              T, chem->specs[i].name);
      exit(2);
    }

    /* assign the appropriate coefficients */
    a = chem->specs[i].aH[iT];
    b = chem->specs[i].aCp[iT];

    /* compute leading coefficient */
    R_W = R/chem->W[i];

    /* ----- compute the enthalpy ----- */
    h[i]  = R_W*(a[0] + T*(a[1] + T*(a[2] + T*(a[3] + T*(a[4] + T*a[5])))));

    /* ----- compute the specific heat ----- */
    cp[i] = R_W*(b[0] + T*(b[1] + T*(b[2] + T*(b[3] + T*b[4]))));
  }
}
/* -------------------------------------------------------------------------- */

void calcSpeciesThermo(molStruct *sp,double *h,double *cp,double *S,double T)
{
  double *a, *b, *s, R = 8314.34, R_W;
  if (T < sp->Trange[2] && T > sp->Trange[1]) {
    a = sp->aH[1];
    b = sp->aCp[1];
    s = sp->aS[1];
  }
  else if (T >= sp->Trange[0]) {
    a = sp->aH[0];
    b = sp->aCp[0];
    s = sp->aS[0];
  }
  else {
    fprintf(stderr, "Error: temperature (%f) is out of range for species %s\n",
            T, sp->name);
    exit(2);
  }

  R_W = R/sp->W;

  /* ----- compute the enthalpy ----- */
  *h  = R_W*(a[0] + T*(a[1] + T*(a[2] + T*(a[3] + T*(a[4] + T*a[5])))));
  
  /* ----- compute the specific heat ----- */
  *cp = R_W*(b[0] + T*(b[1] + T*(b[2] + T*(b[3] + T*b[4]))));

  /* ----- compute the entropy ----- */
  *S  = R_W*(s[0]*log(T) + s[1] + T*(s[2] + T*(s[3] + T*(s[4] + T*s[5]))));
}
/* -------------------------------------------------------------------------- */

void ComputeProductionRatesComp(double *cdot, double *w, double *k, double *c, 
                                double *M, double temp, double pressure,
                                chemStruct *chem)
{
  ComputeProductionRates(cdot,w,k,c,M,temp,pressure);
}
/* -------------------------------------------------------------------------- */

void ComputeProductionRatesInterp(double *cdot, double *w, double *k, double *c,
                                  double *M, double temp, double pressure,
                                  chemStruct *chem)
{
  int i,j,iR,iSp;
  double R    = 8314.34; /* J/(kmol-K) */
  double logT = log(temp);
  double RT   = R*temp;
  double cM;
  double kF, kR;
  double wF, wR;
  double k0, kInf;
  molStruct *sp;
  rxnStruct *r;

  /* ----- compute the rate coefficients ----- */
  for (iR=0; iR<chem->nRxn; iR++) {
    r = &chem->rxns[iR];
    kF = kR = 0.0;

    /* compute the forward rate */
    cM = 1.0;
    /* check if the reaction is pressure dependent */
    if (r->pressure > 0) {
      for (i=0, cM=0.0; i<chem->nSpecies; i++)
        cM += r->concF[i]*c[i];
    }

    kF = r->A*exp(r->b*logT - r->E/RT);    
    if (r->pressure > 1) { /* reaction rate is in the fall-off region */
      /* compute the high and low pressure limit rate parameters */
      k0   = r->A_low*exp(r->b_low*logT - r->E_low/RT);
      kInf = kF;

      /* apply the fall-off function */
      kF = calcFallOffRateCoeff(temp,logT,RT,pressure,k0,kInf,cM,r);
    }

    /* compute the reverse rate */
    if (r->reverse) {
      if (r->reverse == 3) { /* backward coefficients exist */
        kR = r->Ar*exp(r->br*logT - r->Er/RT);
        if (r->pressure > 1) {
          /* compute the high and low pressure limit rate parameters */
          k0   = r->Ar_low*exp(r->br_low*logT - r->Er_low/RT);
          kInf = kR;
          
          /* apply the fall-off function */
          kR = calcFallOffRateCoeff(temp,logT,RT,pressure,k0,kInf,cM,r);
        }
      }
      else { /* compute exact rate based on equilibrium */
        kR = calcReverseRxnRate(temp, logT, kF, r, chem);
      }
    }

    /* compute the rate of progress */
    wF = kF;
    wR = 0.0; 
    for (i=0; i<r->nR; i++) { /* forward */
      for (j=0; j<r->nuR[i]; j++)
        if (r->iR[i] >= 0)
          wF *= c[r->iR[i]];
    } /* end loop over reactants */
    if (r->reverse) { /* reverse */
      /* compute the reverse rate */
      wR = kR;
      for (i=0; i<r->nP; i++) {
        for (j=0; j<r->nuP[i]; j++)
          if (r->iP[i] >= 0)
            wR *= c[r->iP[i]];
      } /* end loop over products */
    } /* end reverse */
    w[iR] = wF - wR;
    if (r->pressure == 1)
      w[iR] *= cM;
  } /* end loop over reactions */

  /* ----- compute the production rates for each species ----- */
  for (iSp=0; iSp<chem->nSpecies; iSp++) {
    sp = &chem->specs[iSp];
    cdot[iSp] = 0.0;
    /* loop over the number of reactions */
    for (i=0; i<sp->nRxn; i++)
      cdot[iSp] += sp->nuRxn[i]*w[sp->iRxn[i]];
  } /* end loop over species */
}
/* -------------------------------------------------------------------------- */

/* Reverse reaction rate computed according to:
      k_r = k_f/K_c
   where K_c is the equilibrium constant, which is computed according to:
      K_c = (P_atm/RT)^sum(nu_i)exp(delS/R - delH/RT)
 */
double calcReverseRxnRate(double T, double logT, double kf, rxnStruct *r,
                          chemStruct *chem)
{
  int i;
  double R  = 8314.34;
  double RT = R*T;
  double mu, sumNu, logKp, logKc, *a;
  molStruct *sp;

  /* loop over reactants */
  sumNu = logKp = 0.0;
  for (i=0; i<r->nS; i++) {
    sp = &chem->specs[r->iS[i]];

    mu = calcGibbsFreeEnergy(T, sp);
    
    sumNu += r->nu[i];
    logKp += r->nu[i]*mu;
  } /* end loop over species */

  logKc = -logKp/RT - sumNu*log(RT/101325);

  return kf/exp(logKc);
}
/* -------------------------------------------------------------------------- */

double calcGibbsFreeEnergy(double T, molStruct *sp)
{
  double R  = 8314.34;
	double mu = 0.0;
  double *a;

  /* get the correct coefficients for the temperature */
	if (T < sp->Trange[2] && T > sp->Trange[1]) {
    a = sp->therm[1];
  }
  else if (T > sp->Trange[0]) {
    a = sp->therm[0];
  }

	mu  = a[0]*(1.0-log(T)) + a[5]/T - a[6];
	mu -= 0.5*T*(a[1] + T*(a[2]/3.0 + T*(a[3]/6.0 + 0.1*T*a[4])));

	return mu*R*T;
}
/* -------------------------------------------------------------------------- */

double calcFallOffRateCoeff(double T, double logT, double RT, double P,
                            double k0, double kInf, double conc, rxnStruct *r)
{
  double Pr,logPr;  /* reduced pressure*/
  double F = 1.0;   /* Troe function (default unity value is for Lindemann */
  double c,d,n,Fc;  /* Troe or SRI coefficients */
  double logFc;
  double tmp;

  if (conc <= 0.0) {
    conc = P/(RT);
  }

  /* compute the reduced pressure */
  Pr = k0*conc/kInf;

  if (r->troe) { /* use Troe form */
    logPr = log10(Pr);
    Fc = r->troeA[0]*exp(-T/r->troeT[0]) + r->troeA[1]*exp(-T/r->troeT[1])
      + r->troeA[2]*exp(-r->troeT[2]/T);
    logFc = log10(Fc);
    c = - 0.4  - 0.67*logFc;
    n =   0.75 - 1.27*logFc;
    d =   0.14;
    tmp = (logPr+c)/(n-d*(logPr+c));
    F = pow(Fc, 1.0/(1.0+tmp*tmp));
  }
  else if (r->sri) { /* use SRI form */
    logPr = log10(Pr);
    tmp   = 1.0/(1.0+logPr*logPr);
    F     = pow(r->sriA[0]*exp(-r->sriA[1]/T)+exp(-T/r->sriA[2]), tmp);
    F    *= r->sriA[3]*pow(T, r->sriA[4]);
  }

  return kInf*(Pr/(1.0+Pr))*F;
}
/* -------------------------------------------------------------------------- */

void chemFitReverseRxnRate(rxnStruct *r, chemStruct *chem)
{
  int    i;
  double RGAS = 8314.34;

  /* grid and data  parameters */
  double P_low  = 1e5;   /* TODO: should this be at different pressure */
  double P_high = 20e5;
  int    nT = 100;
  double dTi;
  double Tmin = 500;
  double Tmax = 2000;
  double *T, *Tinv, logT, RT;
  double cM;
  double kf, kf_low;
  double *kr, *kr_low;
 
  T      = (double*) calloc(nT, sizeof(double));
  Tinv   = (double*) calloc(nT, sizeof(double));
  kr     = (double*) calloc(nT, sizeof(double));
  kr_low = (double*) calloc(nT, sizeof(double));

  /* polynomial fit parameters */
  int     nFit   = 3; /* linear */
  double  chisq;
  double *yDat   = (double*)  calloc(nT,   sizeof(double));
  double *sig    = (double*)  calloc(nT,   sizeof(double));
  double *aCoeff = (double*)  calloc(nFit, sizeof(double));
  double **covar = (double**) calloc(nFit, sizeof(double*));
  double **uFit  = (double**) calloc(nT,   sizeof(double*));
  double **vFit  = (double**) calloc(nFit, sizeof(double*));
  double *wFit   = (double*)  calloc(nFit, sizeof(double));
  int    *ia     = (int*)     calloc(nFit, sizeof(int));
  for (i=0; i<nT; i++) 
    uFit[i]  = (double*) calloc(nFit, sizeof(double));
  for (i=0; i<nFit; i++) {
    covar[i] = (double*) calloc(nFit, sizeof(double));
    vFit[i]  = (double*) calloc(nFit, sizeof(double));
    ia[i]    = 1;
  }

  /* use an equidistant grid in 1/T */
  dTi = (1.0/Tmax - 1.0/Tmin)/(nT-1);
  for (i=0; i<nT; i++) {
    /* build the grid */
    Tinv[i] = 1.0/Tmax - dTi*i;
    T[i] = 1.0/Tinv[i];
    logT = log(T[i]);
    RT   = RGAS*T[i];
    /* compute the forward rate coefficient */
    kf = r->A*exp(r->b*logT - r->E/RT);
    /* compute and store the reverse rate */
    kr[i] = calcReverseRxnRate(T[i], logT, kf, r, chem);
    if (r->pressure > 1) {
      kf_low    = r->A_low*exp(r->b_low*logT - r->E_low/RT);
      kr_low[i] = calcReverseRxnRate(T[i], logT, kf_low, r, chem);
    }

    sig[i] = 1.0;
  } /* end loop over Tinv grid */

  /* use a least-squares fit */
  for (i=0; i<nFit; i++) aCoeff[i] = 0.0;
  for (i=0; i<nT; i++)
    yDat[i] = log(kr[i]);
  lfit(T,yDat,sig,nT,aCoeff,ia,nFit,covar,&chisq,fitArrhenius);
/*   svdfit(T,yDat,sig,nT,aCoeff,nFit,uFit,vFit,wFit,&chisq,fitArrhenius); */

  /* assign the parameters */
  r->Ar = exp(aCoeff[0]);
  r->br = aCoeff[1];
  r->Er = aCoeff[2]*RGAS;
  
  if (r->pressure > 1) {
    for (i=0; i<nFit; i++) aCoeff[i] = 0.0;
    for (i=0; i<nT; i++)
      yDat[i] = log(kr_low[i]);
    lfit(T,yDat,sig,nT,aCoeff,ia,nFit,covar,&chisq,fitArrhenius);
/*     svdfit(T,yDat,sig,nT,aCoeff,nFit,uFit,vFit,wFit,&chisq,fitArrhenius); */
    
    /* assign the parameters */
    r->Ar_low = exp(aCoeff[0]);
    r->br_low = aCoeff[1];
    r->Er_low = aCoeff[2]*RGAS;
  }

  /* write the data to a file for checking */
  if (chem->revRxnDBG) {
    char fileName[256];
    sprintf(fileName, "reaction_fit_%03i.dat", r->no);
    FILE *pF = fopen(fileName, "w");
    double *kr_fit     = (double*) calloc(nT, sizeof(double));
    double *kr_fit_low = (double*) calloc(nT, sizeof(double));
    double k_ref;
    
    fprintf(pF, "# A  = %15.6e, b  = %.6f, E  = %15.6e\n", r->A, r->b, r->E);
    fprintf(pF, "# Ar = %15.6e, br = %.6f, Er = %15.6e\n", r->Ar, r->br, r->Er);
    if (r->pressure > 1) {
      fprintf(pF, "# A_low  = %15.6e, b_low  = %.6f, E_low  = %15.6e\n",
              r->A_low, r->b_low, r->E_low);
      fprintf(pF, "# Ar_low = %15.6e, br_low = %.6f, Er_low = %15.6e\n",
              r->Ar_low, r->br_low, r->Er_low);
    }
    fprintf(pF, "# eConverter = %g\n", 1e6);
    fprintf(pF, "# aConverter = %g\n", pow(1.0e-3, r->orderF[1]-1));
    fprintf(pF, "# aRConvert  = %g\n", pow(1.0e-3, r->orderR[1]-1));
    fprintf(pF, "# T[K]<1> 1000/T[1/K]<2> kr[]<3> ln(kr)[]<4> kr_fit[]<5> ");
    fprintf(pF, "ln(kr_fit)[]<6> ");
    if (r->pressure > 1) {
      fprintf(pF, "kr_low<7> ln(kr_low)<8> kr_low_fit<9> ln(kr_low_fit)<10>");
    }
    fprintf(pF, "\n");
    
    for (i=0; i<nT; i++) {
      logT = log(T[i]);
      RT   = RGAS*T[i];
      kr_fit[i] = r->Ar*exp(r->br*logT - r->Er/RT);
      if (r->pressure > 1) {
        kr_fit_low[i] = r->Ar_low*exp(r->br_low*logT - r->Er_low/RT);
      }
      
      fprintf(pF, "%15.6e %15.6e %15.6e %15.6e %15.6e  %15.6e ",T[i],1000/T[i],
              kr[i], log(kr[i]), kr_fit[i], log(kr_fit[i]));
      if (r->pressure > 1) {
        fprintf(pF, "%15.6e %15.6e %15.6e %15.6e ", kr_low[i], log(kr_low[i]),
                kr_fit_low[i], log(kr_fit_low[i]));
      }
      fprintf(pF, "\n");
    }
    fclose(pF);
    free(kr_fit); free(kr_fit_low);
  }
    
  free(T); free(Tinv); free(kr); free(kr_low);
  free(yDat); free(sig); free(aCoeff); free(ia);
  for (i=0; i<nT; i++) free(uFit[i]);
  free(uFit); free(wFit);
  for (i=0; i<nFit; i++) {
    free(covar[i]); free(vFit[i]);
  }
  free(covar); free(vFit);
}
/* -------------------------------------------------------------------------- */

/* Fits the function 
      ln(k) = ln(A) + b*log(x) - E/(R*x)
   Expects x = 1/T
 */
void fitArrhenius(double x, double *p, int np)
{
  p[0] = 1.0;
  p[1] = log(x);
  p[2] = -1.0/x; /* -1.0/(R*T) */
}
/* -------------------------------------------------------------------------- */

int addSpecies(double newNu, char *newSp, int n, double **nu, char ***sp)
{
  /* add species */
  n++;
  sp[0]      = (char**)  realloc(sp[0], n*sizeof(char*));
  sp[0][n-1] = (char*) calloc(32, sizeof(char));
  strcpy(sp[0][n-1], newSp);

  /* only add coefficient if reaction */
  if (newNu > 0.0 && nu != NULL) {
    nu[0]      = (double*) realloc(nu[0], n*sizeof(double));
    nu[0][n-1] = newNu;
  }
  return n;
}
/* -------------------------------------------------------------------------- */

int removeSpecies(char *remSp, int n, double **nu, char ***sp)
{
  int i,rem = -1;

  /* find index of species to remove */
  for (i=0; i<n; i++) {
    if (!strcmp(remSp, sp[0][i])) {
      rem = i;
      break;
    }
  } /* end loop over species */

  if (rem < 0) {
    fprintf(stderr,"Warning: couldn't remove species %s from reaction\n",remSp);
    return n;
  }

  /* copy all the species back */
  for (i=rem; i<n-1; i++) {
    nu[0][i] = nu[0][i+1];
    strcpy(sp[0][i], sp[0][i+1]);
  }

  /* reallocate memory */
  free(sp[0][n-1]);
  n--;
  sp[0] = (char**)  realloc(sp[0], n*sizeof(char*));
  nu[0] = (double*) realloc(nu[0], n*sizeof(double));
  return n;
}
/* -------------------------------------------------------------------------- */

int removeDuplicateSpecies(int n, double **nu, char ***sp)
{
  int i,j,k;

  for (i=0; i<n-1; i++) {
    for (j=i+1; j<n; j++) {
      if (!strcmp(sp[0][i], sp[0][j])) {
        /* remove the species by shifting all species forward */
        for (k=j+1; k<n; k++) {
          strcpy(sp[0][k-1], sp[0][k]);
          nu[0][k-1] = nu[0][k];
        }
        /* update the coefficients and number of reactants */
        nu[0][i] += nu[0][j];
        n--;
        /* update memory allocation */
        free(sp[0][n]);
        sp[0] = (char**)  realloc(sp[0], n*sizeof(char*));
        nu[0] = (double*) realloc(nu[0], n*sizeof(double));
      }
    } /* end loop over rest of species */
  } /* end loop over species */
  return n;
}
/* -------------------------------------------------------------------------- */

void calcTotalEnthalpy(int n1, int n2, int nVars, int sT, double *Y, double *H,
                       chemStruct *chem)
{
  int i, j, k, r, row, col;

  double *cp = (double*) calloc(chem->nSpecies, sizeof(double));
  double *h  = (double*) calloc(chem->nSpecies, sizeof(double));

  for (k=0; k<n2; k++) {
    row = k*n1*nVars;
    r   = k*n1;
    for (j=0; j<n1; j++) {
      col = j*nVars;
      chemComputeThermoData(h, cp, Y[row+col+sT], chem);
      H[r+j] = 0.0;
      for (i=0; i<chem->nSpecSolve; i++)
        H[r+j] += h[i]*Y[row+col+i];
    } /* end loop over columns */
  } /* end loop over rows */

  free(cp); free(h);
}
/* -------------------------------------------------------------------------- */

double calcTfromH(int nSpecies, int nSolve, double H, double *Y, double guess,
                  chemStruct *chem)
{
  int i, k, count = 0;
  double temp = (guess > 0.0) ? guess : 1000;
  double deltaT = 0.0;
  double hSum, cpSum;

  double *h  = (double*) calloc(nSpecies, sizeof(double));
  double *cp = (double*) calloc(nSpecies, sizeof(double));

  do {
    chemComputeThermoData(h, cp, temp, chem);
    hSum = cpSum = 0.0;
    for (i=0; i<nSolve; i++) {
      cpSum += cp[i]*Y[i];
      hSum  +=  h[i]*Y[i];
    }
    deltaT = -(hSum - H)/cpSum;
    temp += deltaT;
    if (++count > 1000) {
      fprintf(stderr,"ERROR: temperature iteration for h = %g not converged\n",H);
      exit(2);
    }
  } while (fabs(deltaT/temp) > 1.0e-8);
  free(h); free(cp);
  return temp;
}
/* -------------------------------------------------------------------------- */

