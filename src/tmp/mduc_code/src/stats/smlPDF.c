#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "smlPDF.h"
#include "solveUtils.h"

int calcSML_PDF(int n, double *x, double mean, double var, double *pdf,
                int *pdfBound)
{
  int    i, j, k, check, pdfFlag;
  double dx;          // grid spacing for scaling
  double sigma;       // standard deviation
  double sumPDF;      // for scaling the pdf
  double A[3];        // vector containing parameters for root-finding
  statStruct sml;     // struct for storing info needed to evaluation function

  // zero out the pdf
  for (j=0; j<n; j++) pdf[j] = 0.0;

  // ensure that the mean and variance are properly bounded
  clipMoments1D(&mean, &var);

  // check for special cases of pdf, e.g. delta function
  pdfFlag = checkPDF_1D(n, x, mean, var, pdf, pdfBound);

  // if pdfFlag is greater than zero, special case was found, so exit
  if (pdfFlag > 0)
    return pdfFlag;

  // allocate the grid properties to be passed to the root finder
  sml.nMom = 3;
  sml.n1 = n;
  sml.n2 = 1;
  sml.x1 = x;
  sml.x2 = NULL;

  // allocate a priori pdf to be constant (1)
  sml.q    = (double**) calloc(1, sizeof(double*));
  sml.q[0] = (double*)  calloc(n, sizeof(double));
  for (j=0; j<n; j++) sml.q[0][j] = 1.0;

  // set the mean and variance
  sml.mean[0] = mean;
  sml.var[0]  = var;

  // initial guess of lagrange multipliers based on Gaussian
  A[0] = -0.5*mean*mean/var;
  A[1] = mean/var;
  A[2] = -0.5/var;

  // use globally convergent Newton method to find roots
  newt(A, 3, &check, &sml, smlRes1D, smlJac1D);

  // compute the pdf using the values returned
  sumPDF = 0.0;
  for (j=0; j<n; j++)
    {
      if (j == 0)
        dx = 0.5*(x[1]-x[0]);
      else if (j == n-1)
        dx = 0.5*(x[j] - x[j-1]);
      else
        dx = 0.5*(x[j+1] - x[j-1]);

      pdf[j] = smlDist1D(x[j], A)*dx;

      sumPDF += pdf[j];
    }

  // normalize the pdf
  for (j=0; j<n; j++)
    pdf[j] /= sumPDF;

  pdfBound[0] = 0;
  pdfBound[1] = n-1;

  free(sml.q[0]);
  free(sml.q);
  return 1;
}
/* -------------------------------------------------------------------------- */

int calcSML_PDF2D(int *n, double **x, double mean1, double var1, double mean2,
                  double var2, double var, double **pdf, int **pdfBound,
                  int transform)
{
  double PI = 3.1415926535897931;
  int    i, j, k, check, flag, pdfFlag;
  // grid variables
  int    n1, n2;   
  double *x1, *x2, Z1, Z2, x1cur, x2cur;
  double dx1, dx2, dA;
  double var1Max, var2Max, covar;
  double sigma1, sigma2;

  double A[6];               // array containing parameters for root find
  double sumPDF, P;
  double rho_n, rho_sq, tmp; // tmp vars for computing 2D Gaussian
  int initGuess = 1;

  statStruct  sml;  // allocate struct for passing info

  // init grid parameters
  sml.nMom = 6;
  n1 = sml.n1 = n[0];
  n2 = sml.n2 = n[1];
  x1 = sml.x1 = x[0];
  x2 = sml.x2 = x[1];
  sml.transform = transform;

  // zero out the pdf
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      pdf[k][j] = 0.0;

  // ensure that the mean and variance are properly bounded
  clipMoments2D(&mean1, &var1, &mean2, &var2, &var, &covar);

  // treat special cases for the pdf
  pdfFlag = checkPDF_2D(n, x, mean1, var1, mean2, var2, var, pdf, pdfBound,
                        transform, 1);

  // exit if special case was found
  if (pdfFlag > 0)
    return pdfFlag;

  // allocate a priori pdf (constant)
  sml.q  = (double**) calloc(n2, sizeof(double*));
  for (k=0; k<n2; k++)
    {
      sml.q[k] = (double*) calloc(n1, sizeof(double));
      for (j=0; j<n1; j++)
        {
          sml.q[k][j]   = 1.0;
        } // end loop over columns
    } // end loop over rows

  // set the mean and variance vectors for passing
  sml.moments[0] = 1.0;
  sml.moments[1] = sml.mean[0] = mean1;
  sml.moments[2] = sml.mean[1] = mean2;
  sml.moments[3] = sml.var[0]  = var1;
  sml.moments[4] = sml.var[1]  = var2;
  sml.moments[5] = sml.var[2]  = covar;

  // initial guess of lagrange multipliers
  if (initGuess == 0)
    {
      // uniform distribution
      A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = 0.0;
    }
  else if (initGuess == 1)
    {
      // 2D Gaussian
      rho_sq = covar*covar/(var1*var2);
      rho_n  = covar/(var1*var2);
      if (rho_sq != 1.0)
        tmp = -0.5/(1.0-rho_sq);
      else
        tmp = 0.0;

      A[0] = tmp*(mean1*mean1/var1 - 2.0*mean1*mean2*rho_n + mean2*mean2/var2);
      A[1] = -2.0*tmp*(mean1/var1 - mean2*rho_n);
      A[2] = -2.0*tmp*(mean2/var2 - mean1*rho_n);
      A[3] = tmp/var1;
      A[4] = tmp/var2;
      A[5] = -2.0*tmp*rho_n;
    }

#ifdef MDUC_PDF_KINSOL
  // init KINSOL if using
  sdialStruct kin;  // allocate struct for holding KINSOL memory
  /* set up the solver */
  initKinSML(&kin, &sml, A);

  /* set the residual and Jacobian functions */
  sml.fcn = &smlRes2D;
  sml.jac = &smlJac2D;
  
  /* call the KINSOL solver */
  flag = KINSol(kin.mem, kin.Y, KIN_LINESEARCH, kin.Yscale, kin.Fscale);
  if (flag < 0)
    fprintf(stderr, "KINSOL Error: %i\n", flag);

  /* set the solution into the parameter array */
  for (i=0; i<6; i++)
    A[i] = kin.Ydata[i];

  destroyKinSML(&kin);
#elif MDUC_PDF_BROYDEN
  /* use Broyden's method to get parameters */
  broyden(A, 6, &check, &sml, smlFunc2D, smlJac2D);
#elif MDUC_PDF_LM
  double lmInfo[LM_INFO_SZ];

  /* call the external Levenberg-Marquardt library */
  dlevmar_der(smlFunc2D_LM, smlJac2D_LM, A, sml.moments, 6, 6, 200, NULL,
              lmInfo, NULL, NULL, &sml);

  if (lmInfo[6] == 3) check = 1;
#else
  /* use Newton solver to get parameters */
  newt(A, 6, &check, &sml, smlRes2D, smlJac2D);
#endif /* MDUC_PDF_KINSOL */

  // sum over the entire grid
  sumPDF = 0.0;
  for (k=0; k<n2; k++)
    {
      if (k == 0)
        {
          dx2 = 0.5*(x2[1]-x2[0]);
          x2cur = 0.5*dx2;
        }
      else if (k == n2-1)
        {
          dx2 = 0.5*(x2[k] - x2[k-1]);
          x2cur = x2[k] - 0.5*dx2;
        }
      else
        {
          dx2 = 0.5*(x2[k+1] - x2[k-1]);
          x2cur = x2[k];
        }

      // start loop over columns
      for (j=0; j<n1; j++)
        {
          if (j == 0)
            {
              dx1 = 0.5*(x1[1]-x1[0]);
              x1cur = 0.5*dx1;
            }
          else if (j == n1-1)
            {
              dx1 = 0.5*(x1[j] - x1[j-1]);
              x1cur = x1[j] - 0.5*dx1;
            }
          else
            {
              dx1 = 0.5*(x1[j+1] - x1[j-1]);
              x1cur = x1[j];
            }

          // get current local values of Z1 and Z2
          if (transform == 0)
            {
              if (j == n1-1-k)
                {
                  Z1 = x1cur - dx1/6.0;
                  Z2 = x2cur - dx2/6.0;
                  dA = 0.5*dx1*dx2;
                }
              else if (j > n1-1-k)
                {
                  Z1 = 0.2;
                  Z2 = 0.2;
                  dA = 0.0;
                }
              else
                {
                  Z1 = x1cur;
                  Z2 = x2cur;
                  dA = dx1*dx2;
                }
            }
          else if (transform == 1)
            {
              Z1 = x1cur;
              Z2 = x2cur*(1-x1cur);
              dA = dx1*dx2*(1-x1cur);
            }
          else if (transform == 2)
            {
              Z1 = x1cur*(1.0 - x2cur);
              Z2 = x1cur*x2cur;
              dA = dx1*dx2*x1cur;
            }

          P = smlDist2D(Z1, Z2, A);
          
          pdf[k][j] = P*dA;

          sumPDF += pdf[k][j];
        } // end loop over columns
    } // end loop over rows

  FILE *pFile = fopen("smlParams.dat", "r");

  if (!pFile)
    {
      // write header
      pFile = fopen("smlParams.dat", "w");
      fprintf(pFile, "# A0 A1 A2 A3 A4 A5\n");
    }
  fclose(pFile);

  pFile = fopen("smlParams.dat", "a");
  for (i=0; i<6; i++)
    fprintf(pFile, "%15.6e ", A[i]);
  fprintf(pFile, "\n");

  fclose(pFile);

  // normalise the pdf
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      pdf[k][j] /= sumPDF;
  
  // set bounds of pdf
  for (k=0; k<n2; k++)
    {
      pdfBound[k][0] = 0;
      pdfBound[k][1] = n1-1;
      free(sml.q[k]);
    }
  free(sml.q);

  return -1;
}
/* -------------------------------------------------------------------------- */

void smlRes1D(int nF, double *a, double *f, void *args)
{
  /* dereference storage variables */
  statStruct *sml = (statStruct*) args;

  /* compute the moments from the parameters */
  smlFunc1D(nF, a, f, args);

  /* compute the residual */
  f[0] -= 1.0;
  f[1] -= sml->mean[0];
  f[2] -= sml->var[0];
}
/* -------------------------------------------------------------------------- */

void smlFunc1D(int nF, double *a, double *f, void *args)
{
  int i, j, k, m;
  int n1;
  double *x1, dx1;
  double mean, var, *q;
  double P;

  // dereference storage variables
  statStruct *sml = (statStruct*) args;
  n1   = sml->n1;
  x1   = sml->x1;
  q    = sml->q[0];
  mean = sml->mean[0];
  var  = sml->var[0];

  // set all values to zero
  for (i=0; i<nF; i++) f[i] = 0.0;

  for (j=0; j<n1; j++)
    {
      if (j == 0)
        dx1 = 0.5*(x1[1]-x1[0]);
      else if (j == n1-1)
        dx1 = 0.5*(x1[j] - x1[j-1]);
      else
        dx1 = 0.5*(x1[j+1] - x1[j-1]);

      P = smlDist1D(x1[j], a);

      // normalisation coefficient
      f[0] += q[j]*P*dx1;

      // mean coefficient
      f[1] += q[j]*P*x1[j]*dx1;

      // variance coefficient
      f[2] += q[j]*P*(x1[j]-mean)*(x1[j]-mean)*dx1;

    } // end loop over columns
} 
/* -------------------------------------------------------------------------- */

void smlJac1D(int nF, double *a, double **r, void *args)
{
  int i, j, k, m;
  int n1;
  double *x1, dx1;
  double mean, var, *q;
  double P, dTerm[3];


  // dereference storage variables
  statStruct *sml = (statStruct*) args;
  n1   = sml->n1;
  x1   = sml->x1;
  q    = sml->q[0];
  mean = sml->mean[0];
  var  = sml->var[0];

  // set all values to zero
  for (k=0; k<nF; k++)
    for (i=0; i<nF; i++)
      r[k][i] = 0.0;

  for (j=0; j<n1; j++)
    {
      if (j == 0)
        dx1 = 0.5*(x1[1]-x1[0]);
      else if (j == n1-1)
        dx1 = 0.5*(x1[j] - x1[j-1]);
      else
        dx1 = 0.5*(x1[j+1] - x1[j-1]);

      P = smlDist1D(x1[j], a);
      dTerm[0] = 1.0;
      dTerm[1] = x1[j];
      dTerm[2] = x1[j]*x1[j];

      // loop over number or roots
      for (i=0; i<nF; i++)
        {
          // normalisation coefficient
          r[0][i] += q[j]*P*dTerm[i]*dx1;

          // mean coefficient
          r[1][i] += q[j]*P*dTerm[i]*x1[j]*dx1;

          // variance coefficient
          r[2][i] += q[j]*P*dTerm[i]*(x1[j]-mean)*(x1[j]-mean)*dx1;
        }
    } // end loop over columns
} 
/* -------------------------------------------------------------------------- */

void smlRes2D(int nF, double *a, double *f, void *args)
{
  /* dereference storage variables */
  statStruct *sml = (statStruct*) args;

  /* compute the moments from the parameters */
  smlFunc2D(nF, a, f, args);

  /* compute the residual */
  f[0] -= 1.0;
  f[1] -= sml->mean[0];
  f[2] -= sml->mean[1];
  f[3] -= sml->var[0];
  f[4] -= sml->var[1];
  f[5] -= sml->var[2];
}
/* -------------------------------------------------------------------------- */

void smlFunc2D(int nF, double *a, double *f, void *args)
{
  int i, j, k, n, m;
  int nDim, *nx, n1, n2;
  double **grid, *x1, *x2, Z1, Z2;
  double *mean, *var;
  double **q;                   // a priori pdf
  double dx1, dx2, dA;
  double P;
  double x1cur, x2cur;
  
  // dereference variables passed from solver
  statStruct *sml = (statStruct*) args;

  // reallocate to local variables for readability
  n1 = sml->n1;
  n2 = sml->n2;
  x1 = sml->x1;
  x2 = sml->x2;
  q = sml->q;
  mean = sml->mean;
  var = sml->var;

  // set all values to zero
  for (i=0; i<nF; i++) f[i] = 0.0;

  // sum over the entire grid
  for (k=0; k<n2; k++)
    {
      if (k == 0)
        {
          dx2 = 0.5*(x2[1]-x2[0]);
          x2cur = 0.5*dx2;
        }
      else if (k == n2-1)
        {
          dx2 = 0.5*(x2[k] - x2[k-1]);
          x2cur = x2[k] - 0.5*dx2;
        }
      else
        {
          dx2 = 0.5*(x2[k+1] - x2[k-1]);
          x2cur = x2[k];
        }

      for (j=0; j<n1; j++)
        {
          if (j == 0)
            {
              dx1 = 0.5*(x1[1]-x1[0]);
              x1cur = 0.5*dx1;
            }
          else if (j == n1-1)
            {
              dx1 = 0.5*(x1[j] - x1[j-1]);
              x1cur = x1[j] - 0.5*dx1;
            }
          else
            {
              dx1 = 0.5*(x1[j+1] - x1[j-1]);
              x1cur = x1[j];
            }

          if (sml->transform == 0)
            {
              if (j == n1-1-k)
                {
                  Z1 = x1cur - dx1/6.0;
                  Z2 = x2cur - dx2/6.0;
                  dA = 0.5*dx1*dx2;
                }
              else if (j > n1-1-k)
                {
                  Z1 = 0.2;
                  Z2 = 0.2;
                  dA = 0.0;
                }
              else
                {
                  Z1 = x1cur;
                  Z2 = x2cur;
                  dA = dx1*dx2;
                }
            }
          else if (sml->transform == 1)
            {
              Z1 = x1cur;
              Z2 = x2cur*(1-x1cur);
              dA = dx1*dx2*(1-x1cur);
            }
          else if (sml->transform == 2)
            {
              Z1 = x1cur*(1.0 - x2cur);
              Z2 = x1cur*x2cur;
              dA = dx1*dx2*x1cur;
            }

          P = smlDist2D(Z1, Z2, a);

          // normalisation condition
          f[0] += q[k][j]*P*dA;

          // compute the means
          f[1] += q[k][j]*P*Z1*dA; // mean1
          f[2] += q[k][j]*P*Z2*dA; // mean2

          // compute the variances
          f[3] += q[k][j]*P*(Z1-mean[0])*(Z1-mean[0])*dA;
          f[4] += q[k][j]*P*(Z2-mean[1])*(Z2-mean[1])*dA;
          f[5] += q[k][j]*P*(Z1-mean[0])*(Z2-mean[1])*dA;
        } // end loop over columns
    } // end loop over rows
} 
/* -------------------------------------------------------------------------- */

void smlJac2D(int nF, double *a, double **r, void *args)
{
  int i, j, k, n, m;
  int nDim, *nx, n1, n2;
  double **grid, *x1, *x2, Z1, Z2;
  double *mean, *var;
  double **q;                   // a priori pdf
  double dx1, dx2, dA;
  double P, dTerm[6];
  double x1cur, x2cur;
  
  // dereference storage variables
  statStruct *sml = (statStruct*) args;

  // reallocate to local variables for readability
  n1 = sml->n1;
  n2 = sml->n2;
  x1 = sml->x1;
  x2 = sml->x2;
  q = sml->q;
  mean = sml->mean;
  var = sml->var;

  // set all values to zero
  for (k=0; k<nF; k++)
    for (i=0; i<nF; i++)
      r[k][i] = 0.0;

  // sum over the entire grid
  for (k=0; k<n2; k++)
    {
      if (k == 0)
        {
          dx2 = 0.5*(x2[1]-x2[0]);
          x2cur = 0.5*dx2;
        }
      else if (k == n2-1)
        {
          dx2 = 0.5*(x2[k] - x2[k-1]);
          x2cur = x2[k] - 0.5*dx2;
        }
      else
        {
          dx2 = 0.5*(x2[k+1] - x2[k-1]);
          x2cur = x2[k];
        }

      for (j=0; j<n1; j++)
        {
          if (j == 0)
            {
              dx1 = 0.5*(x1[1]-x1[0]);
              x1cur = 0.5*dx1;
            }
          else if (j == n1-1)
            {
              dx1 = 0.5*(x1[j] - x1[j-1]);
              x1cur = x1[j] - 0.5*dx1;
            }
          else
            {
              dx1 = 0.5*(x1[j+1] - x1[j-1]);
              x1cur = x1[j];
            }

          if (sml->transform == 0)
            {
              if (j == n1-1-k)
                {
                  Z1 = x1cur - dx1/6.0;
                  Z2 = x2cur - dx2/6.0;
                  dA = 0.5*dx1*dx2;
                }
              else if (j > n1-1-k)
                {
                  Z1 = 0.2;
                  Z2 = 0.2;
                  dA = 0.0;
                }
              else
                {
                  Z1 = x1cur;
                  Z2 = x2cur;
                  dA = dx1*dx2;
                }
            }
          else if (sml->transform == 1)
            {
              Z1 = x1cur;
              Z2 = x2cur*(1-x1cur);
              dA = dx1*dx2*(1-x1cur);
            }
          else if (sml->transform == 2)
            {
              Z1 = x1cur*(1.0 - x2cur);
              Z2 = x1cur*x2cur;
              dA = dx1*dx2*x1cur;
            }

          P = smlDist2D(Z1, Z2, a);

          dTerm[0] = 1.0;
          dTerm[1] = Z1;
          dTerm[2] = Z2;
          dTerm[3] = Z1*Z1;
          dTerm[4] = Z2*Z2;
          dTerm[5] = Z1*Z2;

          for (i=0; i<nF; i++)
            {
              // normalisation condition
              r[0][i] += q[k][j]*P*dTerm[i]*dA;

              // compute the means
              r[1][i] += q[k][j]*P*dTerm[i]*Z1*dA; // mean1
              r[2][i] += q[k][j]*P*dTerm[i]*Z2*dA; // mean2

              // compute the variances
              r[3][i] += q[k][j]*P*dTerm[i]*(Z1-mean[0])*(Z1-mean[0])*dA;
              r[4][i] += q[k][j]*P*dTerm[i]*(Z2-mean[1])*(Z2-mean[1])*dA;
              r[5][i] += q[k][j]*P*dTerm[i]*(Z1-mean[0])*(Z2-mean[1])*dA;
            }
        } // end loop over columns
    } // end loop over rows
} 
/* -------------------------------------------------------------------------- */

/* wrapper functions for external Levenberg-Marquardt library */
void smlFunc1D_LM(double *p, double *hx, int m, int n, void *args)
{
  smlFunc1D(m, p, hx, args);
}
/* -------------------------------------------------------------------------- */

void smlJac1D_LM(double *p, double *j, int m, int n, void *args)
{
  int i,k;
  double **jac;
  jac = (double**) calloc(n, sizeof(double*));
  for (k=0; k<n; k++)
    jac[k] = (double*) calloc(m, sizeof(double));

  smlJac1D(m, p, jac, args);

  for (k=0; k<n; k++)
    for (i=0; i<m; i++)
      j[k*m+i] = jac[k][i];

  for (k=0; k<n; k++) free(jac[k]); free(jac);
}
/* -------------------------------------------------------------------------- */

void smlFunc2D_LM(double *p, double *hx, int m, int n, void *args)
{
  smlFunc2D(m, p, hx, args);
}
/* -------------------------------------------------------------------------- */

void smlJac2D_LM(double *p, double *j, int m, int n, void *args)
{
  int i,k;
  double **jac;
  jac = (double**) calloc(n, sizeof(double*));
  for (k=0; k<n; k++)
    jac[k] = (double*) calloc(m, sizeof(double));

  smlJac2D(m, p, jac, args);

  for (k=0; k<n; k++)
    for (i=0; i<m; i++)
      j[k*m+i] = jac[k][i];

  for (k=0; k<n; k++) free(jac[k]); free(jac);
}
/* -------------------------------------------------------------------------- */

#ifdef MDUC_PDF_KINSOL
/* functions for KINSOL Newton solver */
void initKinSML(sdialStruct *kin, statStruct *sml, double *a)
{
  int i, flag;

  kin->nEq = sml->nMom;

  // allocate the N_Vector space needed
  kin->Y      = N_VNew_Serial(kin->nEq);
  kin->tmp    = N_VNew_Serial(kin->nEq);
  kin->Yscale = N_VNew_Serial(kin->nEq);
  kin->Fscale = N_VNew_Serial(kin->nEq);

  // set pointers to data in N_Vector elements
  kin->Ydata      = NV_DATA_S(kin->Y);
  kin->YscaleData = NV_DATA_S(kin->Yscale);
  kin->FscaleData = NV_DATA_S(kin->Fscale);

  // create and allocate the KINSOL memory object
  kin->mem = KINCreate();
#if SUNDIAL_2_3
  flag = KINMalloc(kin->mem, kinSML_Res, kin->tmp);
  if (flag < 0)
    fprintf(stderr, "KINMalloc error: %i\n", flag);

  // attach a dense solver
  flag = KINDense(kin->mem, kin->nEq);
  if (flag < 0)
    fprintf(stderr, "KINDense error: %i\n", flag);

  // set pointers to data and function calls
  flag = KINSetFdata(kin->mem, sml);
  if (flag < 0)
    fprintf(stderr, "KINSetFdata error: %i\n", flag);
  flag = KINDenseSetJacFn(kin->mem, kinSML_Jac, sml); 
  if (flag < 0)
    fprintf(stderr, "KINDenseSetJacFn error: %i\n", flag);
#else

  flag = KINInit(kin->mem, kinSML_Res, kin->tmp);
  flag = KINDense(kin->mem, kin->nEq);
  flag = KINSetUserData(kin->mem, sml);
  flag = KINDlsSetDenseJacFn(kin->mem, kinSML_Jac); 

#endif /* SUNDIAL_2_3 */

  // set the solution scaling vectors (not at all sure what to do here)
  kin->YscaleData[0] = 1.0;
  kin->YscaleData[1] = 1.0;
  kin->YscaleData[2] = 1.0;
  kin->YscaleData[3] = 1.0;
  kin->YscaleData[4] = 1.0;
  kin->YscaleData[5] = 1.0;

  // set the function scaling vectors
  kin->FscaleData[0] = 1.0;
  kin->FscaleData[1] = fabs(1.0/(sml->mean[0]+1e-9));
  kin->FscaleData[2] = fabs(1.0/(sml->mean[1]+1e-9));
  kin->FscaleData[3] = fabs(1.0/(sml->var[0]+1e-9));
  kin->FscaleData[4] = fabs(1.0/(sml->var[1]+1e-9));
  kin->FscaleData[5] = fabs(1.0/(sml->var[2]+1e-9));

  // copy the initial guess into the solution vector
  for (i=0; i<kin->nEq; i++)
    kin->Ydata[i] = a[i];
}
/* -------------------------------------------------------------------------- */

void destroyKinSML(sdialStruct *kin)
{
  N_VDestroy_Serial(kin->Y);
  N_VDestroy_Serial(kin->tmp);
  N_VDestroy_Serial(kin->Yscale);
  N_VDestroy_Serial(kin->Fscale);
  KINFree(&kin->mem);
}
/* -------------------------------------------------------------------------- */

int kinSML_Res(N_Vector u, N_Vector fval, void *f_data)
{
  // dereference storage variables
  statStruct *sml = (statStruct*) f_data;
  
  realtype *a, *f;

  a = NV_DATA_S(u);
  f = NV_DATA_S(fval);

  sml->fcn(sml->nMom, a, f, sml);
  return 0;
}
/* -------------------------------------------------------------------------- */

#if SUNDIALS_2_3
int kinSML_Jac(long int N, DenseMat J, N_Vector u, N_Vector fu,
               void *jac_data, N_Vector tmp1, N_Vector tmp2)
#else
int kinSML_Jac(int N, N_Vector u, N_Vector fu, DlsMat J, void *jac_data,
               N_Vector tmp1, N_Vector tmp2)
#endif /* SUNDIAL_2_3 */
{
  int i, j;
  // dereference storage variables
  statStruct *sml = (statStruct*) jac_data;
  realtype *a, **Jtmp;

  // allocate memory for temporary Jacobian
  Jtmp = (double**) calloc(N, sizeof(double*));
  for (i=0; i<N; i++) Jtmp[i] = (double*) calloc(N, sizeof(double));

  /* pointers to data */
  a = NV_DATA_S(u);

  /* call the Jacobian function */
  sml->jac(sml->nMom, a, Jtmp, sml);

  // copy the computed Jacobian to J
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      DENSE_ELEM(J,i,j) = Jtmp[i][j];

  for (i=0; i<N; i++) free(Jtmp[i]); free(Jtmp);
  return 0;
}
/* -------------------------------------------------------------------------- */
#endif /* MDUC_PDF_KINSOL */

double smlDist1D(double x, double *a)
{
  return exp(a[0] + a[1]*x + a[2]*x*x);
}
/* -------------------------------------------------------------------------- */

double smlDist2D(double x1, double x2, double *a)
{
  return exp(a[0] + a[1]*x1 + a[2]*x2 + a[3]*x1*x1 + a[4]*x2*x2 + a[5]*x1*x2);
}
/* -------------------------------------------------------------------------- */
