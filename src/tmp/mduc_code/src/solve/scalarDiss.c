#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "scalarDiss.h"
#include "accessMDUC.h"
#include "errFunc.h"
#include "statsUtils.h"
#include "utils.h"

void calcChi1DcounterflowDiffFlame(int nZ, double *Z, int j, double a,
                                   double *chi)
{
  int i;

  double tmp, PI = 3.1415926535897931;
  double inverfc_sq;

  for (i=0; i<nZ; i++)
    {
      // store [erfc^-1(2Z)]^2
      tmp        = inverfc(2*Z[i]);
      inverfc_sq = tmp*tmp;

      // compute the scalar dissipation rate for the local mixture fraction
      chi[i] = a*(j+1)/PI*exp(-2*inverfc_sq);
    }
}
/* ------------------------------------------------------------------------- */

void calcChi1DunsteadyLamMixLayer(int nZ, double *Z, double t, double *chi)
{
  int i;

  double tmp, PI = 3.1415926535897931;
  double inverfc_sq;

  for (i=0; i<nZ; i++)
    {
      // store [erfc^-1(2Z)]^2
      tmp        = inverfc(2*Z[i]);
      inverfc_sq = tmp*tmp;

      // compute the scalar dissipation rate for the local mixture fraction
      chi[i] = 1/(2*PI*t)*exp(-2*inverfc_sq);
    }
}
/* ------------------------------------------------------------------------- */

void calcChi1D(char *form, int n, double *x, double Zref, double Zmax, 
               double chiRef, double *chi)
{
  int i;
  // set the form
  if (!strncmp(form, "log", 3))
    calcChiForm = &calcLogForm;
  else if (!strncmp(form, "erfc", 4) || !strncmp(form, "inverfc",5))
    calcChiForm = &calcErfcForm;

  for (i=0; i<n; i++)
    chi[i] = chiRef*calcChiForm(x[i], Zref, Zmax);
}
/* ------------------------------------------------------------------------- */

void calcChi2D(char *form, int *n, double **x, double *Zref, double *chiRef,
               double **chi1, double **chi12, double **chi2, int trans,
               int scaled)
{
  int j, k;
  double Z1, Z2, Zmax[3] = {1.0, 1.0, 1.0};
  int     n1 = n[0];
  int     n2 = n[1];
  double *x1 = x[0];
  double *x2 = x[1];

  // set the form
  if (!strncmp(form, "log", 3))
    calcChiForm = &calcLogForm;
  else if (!strncmp(form, "erfc", 4) || !strncmp(form, "inverfc",5))
    calcChiForm = &calcErfcForm;

  // scalar dissipation rate of Z0 = 1 - Z1 - Z2 (alternatively Z1+Z2)
  double **chi0 = (double**) calloc(n2, sizeof(double**));
  for (k=0; k<n2; k++)
    chi0[k] = (double*) calloc(n1, sizeof(double));

  // compute dissipation according to specified form
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++)
        {
          // get the actual mixture fractions for this grid point
          getZfromGrid(trans, x1[j], x2[k], &Z1, &Z2);

          if (scaled)
            {
              Zmax[0] = 1.0 - Z1 - Z2;
              Zmax[1] = 1.0 - Z2;
              Zmax[2] = 1.0 - Z1;
            }

          // compute scaled scalar dissipation
          chi0[k][j] = chiRef[0]*calcChiForm(Z1+Z2, Zref[0], Zmax[0]);
          chi1[k][j] = chiRef[1]*calcChiForm(Z1,    Zref[1], Zmax[1]);
          chi2[k][j] = chiRef[2]*calcChiForm(Z2,    Zref[2], Zmax[2]);

          // enforce positive values
          chi0[k][j] = fmax(chi0[k][j], 0.0);
          chi1[k][j] = fmax(chi1[k][j], 0.0);
          chi2[k][j] = fmax(chi2[k][j], 0.0);

          // compute the cross-dissipation term
          chi12[k][j] = 0.5*(chi0[k][j] - chi1[k][j] - chi2[k][j]);
        } /* end loop over columns */
    } /* end loop over rows */

  for (k=0; k<n2; k++) free(chi0[k]);
  free(chi0);
}
/* ------------------------------------------------------------------------- */

double calcLogForm(double Z, double Zref, double Zmax)
{
  if ( Z < 1.0e-10 || Z >= 1.0 || Z/Zmax >= 1.0)
    return 0.0;
  else if ( Zref < 1.0e-10 || Zref >= 1.0 || Zref/Zmax >= 1.0 )
    return 0.0;
  else
    return Z*Z/(Zref*Zref)*log(Z/Zmax)/log(Zref/Zmax);
}
/* ------------------------------------------------------------------------- */

double calcErfcForm(double Z, double Zref, double Zmax)
{
  double tmp, erfc1, erfc2;

  if ( Z < 1.0e-10 || Z >= 1.0)
    return 0.0;
  else if ( Z/Zmax < 1.0e-10 || Z/Zmax >= 1.0 )
    return 0.0;
  else
    erfc1 = inverfc(2.0*Z/Zmax);

  if ( (Zref < 1.0e-10) || (Zref >= 1.0) )
    erfc2 = 0.0;
  else if ( Zref/Zmax < 1.0e-10 || Zref/Zmax >= 1.0)
    erfc2 = 0.0;
  else
    erfc2 = inverfc(2.0*Zref/Zmax);

  return Zmax*Zmax*exp(-2.0*erfc1*erfc1 + 2.0*erfc2*erfc2);
}
/* ------------------------------------------------------------------------- */

double calcQuadForm(double Z, double Zref, double Zmax)
{
  if (Zref > 0.0 && Z < Zmax && Zref < Zmax)
    return Z/Zref*(Z-Zmax)/(Zref-Zmax);
  else
    return 0.0;
}
/* ------------------------------------------------------------------------- */

double calcBetaPDFStoich(double stoich, double mean, double var)
{
  mean = fmin(1.0, fmax(0.0, mean));
  var  = fmin(fabs(var), mean*(1.0 - mean));

  if ( (var == mean*(1.0-mean)) || (mean*(1.0-mean)/var > 2000) )
    return 0.0;
  else
    return betaDist(stoich, mean, var);
}
/* ------------------------------------------------------------------------- */

#ifdef ACCESS_CFD
void calcChi1DPitsch(int nCFD1, int nCFD2, int nCFD3, double ***C_Z1MEAN,
                     double ***C_Z1VAR, double ***C_K, double ***C_D,
                     double ***C_MASS, int n1, double *x1, double Z1ref,
                     double *chi)
{
  // define CFD variables
  double ***Z1Mean;  // mean mixture fraction
  double ***Z1Var;   // mixture fraction variance
  double ***K;       // turbulent kinetic energy
  double ***D;       // turbulent dissipation
  double ***Mass;    // cell mass
  double ***chi1_st; // storage for intermediate chi in each cell

  // set the pointers appropriately for the dimension of the array
  // NOTE: 1D and 2D passed as address of pointer
  if (nCFD3 == 1 && nCFD2 == 1)
    {
      Z1Mean = &C_Z1MEAN;
      Z1Var  = &C_Z1VAR;
      K      = &C_K;
      D      = &C_D;
      Mass   = &C_MASS;
    }
  else
    {
      Z1Mean = C_Z1MEAN;
      Z1Var  = C_Z1VAR;
      K      = C_K;
      D      = C_D;
      Mass   = C_MASS;
    }

  // other local variables for flow solver
  int    i, j, k, p;
  double TINY  = 1.0e-9; // for use in potential divide-by-zero
  double tmp;
  double chiDen;         // temporary storage for integral of denominator
  double pdf1_st, *pdf;  // storage for PDFs
  int    pdfBound[2];    // stores bounds of non-zero PDF values
  double tauI;           // turbulent eddy turnover time
  double c_chi = 2.0;    // const. for modelling mean scalar dissipation
  double chiMean;        // unconditional mean scalar dissipation
  double chiStAvg;       // domain-averaged conditional mean chi
  double sigma;          // standard deviation
  double Zmax;           // maximum mixture fraction in the domain

  // allocate all local memory
  pdf = (double*) calloc(n1, sizeof(double));
  chi1_st = (double***) calloc(nCFD1, sizeof(double**));
  for (i=0; i<nCFD1; i++)
    {
      chi1_st[i] = (double**) calloc(nCFD2, sizeof(double*));

      for(j=0; j<nCFD2; j++)
        chi1_st[i][j] = (double*) calloc(nCFD3, sizeof(double));
    }

  /* ----- get the max in the domain ----- */
  /* this is defined as the mean + 2*sigma */
  Zmax = 0.0;
  for (k=0; k<nCFD3; k++)
    {
      for (j=0; j<nCFD2; j++)
        {
          for (i=0; i<nCFD1; i++)
            {                         
              if ( (Z1Mean[k][j][i] > 0.0) && (Z1Var[k][j][i] > TINY) )
                {
                  sigma = pow(fmax(Z1Var[k][j][i],0), 0.5);
                                                        
                  if (Z1Mean[k][j][i] + 2*sigma > Zmax)
                    Zmax = Z1Mean[k][j][i] + 2*sigma;
                }
            }
        }
    } // end loop over CFD cells

  /* ----- compute mean scalar dissipation for each cell ----- */
  printf("Computing mean scalar dissipation for cells...\n");
  for (k=0; k<nCFD3; k++)
    {
      for (j=0; j<nCFD2; j++)
        {
          for (i=0; i<nCFD1; i++)
            {
              tauI = D[k][j][i]/(K[k][j][i]+TINY); // eddy turnover time

              // compute the unconditional mean scalar dissipation rate
              chiMean = c_chi*tauI*Z1Var[k][j][i];

              /* compute the PDF based on beta function*/
              /* NOTE: this function returns the PDF already scaled with dx
                 i.e., which integrating equivalent to pdf(x)*dx */
              MDUCcalcPDF_1D("beta",n1,x1,Z1Mean[k][j][i],Z1Var[k][j][i],pdf,
                              pdfBound);

              // compute the standar deviation
              sigma = pow(fmax(Z1Var[k][j][i],0), 0.5);

              /* integrate the denominator using the logarithmic functional
                 form of Pitsch */
              chiDen = 0.0;
              for (p=pdfBound[0]; p<=pdfBound[1]; p++)
                chiDen += calcLogForm(x1[p], Z1ref, Zmax)*pdf[p];

              /* compute the conditional mean chi at stoichiometric */
              if (chiDen > TINY)
                chi1_st[k][j][i] = chiMean/chiDen;

            }
        }
    } // end loop over CFD cells

  /* ----- compute the reference scalar dissipation rate ----- */
  /* domain average the conditional mean scalar dissipation rate weighted
     by the stoichiometric surface area */
  printf("Computing reference scalar dissipation rate...\n");
  chiStAvg  = 0.0;
  chiDen    = 0.0;
  for (k=0; k<nCFD3; k++)
    {
      for (j=0; j<nCFD2; j++)
        {
          for (i=0; i<nCFD1; i++)
            {
              /* compute the PDF(Z1_st, 0) */
              pdf1_st = calcBetaPDFStoich(Z1ref,Z1Mean[k][j][i],Z1Var[k][j][i]);
                                        
              /* add to integral for chi1 */
              if (chi1_st[k][j][i] > 0.0)
                {
                  /* sum the contributions to numerator and denominator */
                  tmp = pow(chi1_st[k][j][i],0.5)*Mass[k][j][i]*pdf1_st;
                  chiStAvg += chi1_st[k][j][i]*tmp;
                  chiDen   += tmp;

                  /* Error checking on average chi */
                  if (isnan(chiStAvg))
                    {
                      printf("in cell [%i][%i][%i]\n", k, j, i);
                      printf("C_Z1MEAN = %15.6e, C_Z1VAR = %15.6e\n",
                              Z1Mean[k][j][i], Z1Var[k][j][i]);
                      printf("chi1_st = %15.6e, pdf1_st = %15.6e, tmp = %15.6e\n",
                              chi1_st[k][j][i], pdf1_st, tmp);
                      printf("<chiStoich> is undefined\n");
                      exit(2);
                    } /* end error checking */
                } /* end chi computation for non-zero cell */
            }
        }
    } // end loop over CFD cells
  chiStAvg = chiStAvg/chiDen; // complete the averaging
  printf("<chi1_stoich> = %15.6e\n", chiStAvg);

  /* ---- compute the scalar dissipation rate profiles ----- */
  // again, based on the log form of Pitsch
  for (p=0; p<n1; p++)
    chi[j] = chiStAvg*calcLogForm(x1[p], Z1ref, Zmax);

  
  // free all local memory
  for (i=0; i<nCFD1; i++)
    {
      for (j=0; j<nCFD2; j++)
        free(chi1_st[i][j]);

      free(chi1_st[i]);
    }
  free(chi1_st);
  free(pdf);
}
/* -------------------------------------------------------------------------- */

void binChi1D(int nCFD1, int nCFD2, int nCFD3, double ***C_Z1MEAN,
              double ***C_Z1VAR, double ***C_K, double ***C_D, double ***C_MASS,
              int n1, double *x1, double Z1ref, double *chi)
{
  printf("Binning 1D scalar dissipation rate...\n");

  // define CFD variables
  double ***Z1Mean;  // mean mixture fraction
  double ***Z1Var;   // mixture fraction variance
  double ***K;       // turbulent kinetic energy
  double ***D;       // turbulent dissipation
  double ***Mass;    // cell mass
  double ***chi1_st; // storage for intermediate chi in each cell

  // set the pointers appropriately for the dimension of the array
  // NOTE: 1D and 2D passed as address of pointer
  if (nCFD3 == 1 && nCFD2 == 1)
    {
      Z1Mean = &C_Z1MEAN;
      Z1Var  = &C_Z1VAR;
      K      = &C_K;
      D      = &C_D;
      Mass   = &C_MASS;
    }
  else
    {
      Z1Mean = C_Z1MEAN;
      Z1Var  = C_Z1VAR;
      K      = C_K;
      D      = C_D;
      Mass   = C_MASS;
    }

  /* local variables for flow solver */
  int    i, j, k, p;
  double TINY  = 1.0e-9; /* for use in potential divide-by-zero */
  double tmp;
  double chiMean;
  double *chiDen;        /* temporary storage for integral of denominator */
  double tau;            /* turbulent eddy turnover time */
  double c_chi = 2.0;    /* const. for modelling mean scalar dissipation */

  double dx;
  double w1, w2;         /* weights for distributing chi */

  /* allocate all local memory */
  chiDen = (double*) calloc(n1, sizeof(double));

  for (j=0; j<n1; j++)
    {
      chi[j]    = 0.0;
      chiDen[j] = 0.0;
    }

  for (k=0; k<nCFD3; k++)
    {
      for (j=0; j<nCFD2; j++)
        {
          for (i=0; i<nCFD1; i++)
            {
              if ( (Z1Mean[k][j][i] > 0.0) && (Z1Var[k][j][i] > TINY) )
                {
                  tau  = D[k][j][i]/(K[k][j][i]+TINY); // eddy turnover time

                  // compute the unconditional mean scalar dissipation rate
                  chiMean = c_chi*tau*Z1Var[k][j][i];

                  /* search for the bin that it fits in */
                  p = 0;
                  while (x1[p] < Z1Mean[k][j][i]) p++;
                  p--; /* this gives the grid point immediately left of mean */

                  dx = x1[p+1] - x1[p];

                  /* take care of the boundary points */
                  if (p == 0)
                    {
                      w1 = 0.0;
                      w2 = 1.0;
                    }
                  else if (p == n1-2)
                    {
                      w1 = 1.0;
                      w2 = 0.0;
                    }
                  else
                    {
                      w1 = (x1[p+1] - Z1Mean[k][j][i])/dx;
                      w2 = (Z1Mean[k][j][i] - x1[p])/dx;
                    }

                  chi[p]    += Mass[k][j][i]*w1*chiMean;
                  chiDen[p] += Mass[k][j][i];

                  chi[p+1]    += Mass[k][j][i]*w2*chiMean;
                  chiDen[p+1] += Mass[k][j][i];
                }
            }
        }
    }

  /* compute chi(Z) over the whole grid */
  for (p=0; p<n1; p++)
    {
      if (chiDen[p] > 0.0)
        chi[p] = chi[p]/chiDen[p];
    } /* end loop over flamelet grid */

  free(chiDen);
}
/* ------------------------------------------------------------------------- */
#endif // ACCESS_CFD
