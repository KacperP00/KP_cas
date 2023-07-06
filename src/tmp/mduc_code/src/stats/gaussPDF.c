#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "gaussPDF.h"
#include "solveUtils.h"

int calcGaussPDF(int n, double *x, double mean, double var, double *pdf,
                 int *pdfBound)
{
  int    j, pdfFlag;
  double P, xD;             // for computing delta functions
  double deltaA, deltaB;    // delta functions at boundaries for clipping
  double dx, dxD = 0.01;    // grid spacing for scaling
  double sumPDF;            // for scaling the pdf
  double truncate = 1.0e-6; // stop computing pdf below this limit
  double a0, a1, a2, a3;

  // zero out the pdf
  for (j=0; j<n; j++) pdf[j] = 0.0;

  // ensure that the mean and variance are properly bounded
  clipMoments1D(&mean, &var);

  // check for special cases of pdf, e.g. delta function
  pdfFlag = checkPDF_1D(n, x, mean, var, pdf, pdfBound);

  // if pdfFlag is greater than zero, special case was found, so exit
  if (pdfFlag > 0)
    return pdfFlag;

  // initial guess of lagrange multipliers based on Gaussian
  a0 = 0.398942280401/sqrt(var);   // 1/(sigma*sqrt(2*pi))
  a1 = -0.5*mean*mean/var;
  a2 = mean/var;
  a3 = -0.5/var;

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

      pdf[j] = gaussDist1D(x[j], a0, a1, a2, a3)*dx;

      sumPDF += pdf[j];
    }

  // compute a delta function at x = 0 to represent probability from -inf to 0
  xD = 0.0;
  deltaA = 0.0;
  do {
    xD -= dxD;
    P = gaussDist1D(xD, a0, a1, a2, a3)*dxD;
    deltaA += P;
  } while (P > truncate);

  // compute a delta function at x = 1 to represent probability from 1 to inf
  xD = 1.0;
  deltaB = 0.0;
  do {
    xD += dxD;
    P = gaussDist1D(xD, a0, a1, a2, a3)*dxD;
    deltaB += P;
  } while (P > truncate);

  // add delta functions to the boundaries
  pdf[0]   += deltaA;
  pdf[n-1] += deltaB;
  sumPDF = sumPDF + deltaA + deltaB;

  // normalize the pdf
  for (j=0; j<n; j++)
    pdf[j] /= sumPDF;

  pdfBound[0] = 0;
  pdfBound[1] = n-1;

  return 1;
}
/* -------------------------------------------------------------------------- */

int calcLogNormPDF(int n, double *x, double mean, double var, double *pdf,
                   int *pdfBound)
{
  int    j, pdfFlag;
  double P, xD, xc;         // for computing delta functions
  double deltaB;            // delta functions at boundary for clipping
  double dx, dxD = 0.01;    // grid spacing for scaling
  double sumPDF;            // for scaling the pdf
  double truncate = 1.0e-6; // stop computing pdf below this limit
  double lMean, lVar;       // mean and variance variables for log normal
  double a0, a1, a2, a3;

  // zero out the pdf
  for (j=0; j<n; j++) pdf[j] = 0.0;

  // ensure that the mean and variance are properly bounded
  clipMoments1D(&mean, &var);

  // check for special cases of pdf, e.g. delta function
  pdfFlag = checkPDF_1D(n, x, mean, var, pdf, pdfBound);

  // if pdfFlag is greater than zero, special case was found, so exit
  if (pdfFlag > 0)
    return pdfFlag;

  // compute the log normal parameters
  lVar  = log(var/(mean*mean) + 1.0);
  lMean = log(mean) - 0.5*lVar;

  // initial guess of lagrange multipliers based on Gaussian
  a0 = 0.398942280401/sqrt(lVar);   // 1/(sigma*sqrt(2*pi))
  a1 = -0.5*lMean*lMean/lVar;
  a2 = lMean/lVar;
  a3 = -0.5/lVar;

  // compute the pdf using the values returned
  sumPDF = 0.0;
  for (j=0; j<n; j++)
    {
      if (j == 0)
        {
          dx = 0.5*(x[1]-x[0]);
          xc = 0.5*dx;
        }
      else if (j == n-1)
        {
          dx = 0.5*(x[j] - x[j-1]);
          xc = x[j] - 0.5*dx;
        }
      else
        {
          dx = 0.5*(x[j+1] - x[j-1]);
          xc = x[j];
        }

      pdf[j] = 1.0/xc*gaussDist1D(log(xc), a0, a1, a2, a3)*dx;

      sumPDF += pdf[j];
    }

  // compute a delta function at x = 1 to represent probability from 1 to inf
  xD = 1.0;
  deltaB = 0.0;
  do {
    xD += dxD;
    P = 1.0/xD*gaussDist1D(log(xD), a0, a1, a2, a3)*dxD;
    deltaB += P;
  } while (P > truncate);

  // add delta functions to the boundaries
  pdf[n-1] += deltaB;
  sumPDF   += deltaB;

  // normalize the pdf
  for (j=0; j<n; j++)
    pdf[j] /= sumPDF;

  pdfBound[0] = 0;
  pdfBound[1] = n-1;

  return 1;
}
/* -------------------------------------------------------------------------- */

int calcGaussPDF2D(int *n, double **x, double mean1, double var1, double mean2, 
                   double var2, double var, double **pdf, int **pdfBound,
                   int transform)
{
  double PI = 3.1415926535897931;
  int    i, j, k, check, flag, pdfFlag;
  double Z1, Z2, x1c, x2c;   // grid variables
  double dx1, dx2, dA;       // grid spacing
  double covar;              // covariance
  double a[7];               // array containing parameters distribution
  double sumPDF, P;          // sum of PDF, local probability
  double rho_n, rho_sq, tmp; // tmp vars for computing 2D Gaussian
  double *deltaA, *deltaB, *deltaC, *deltaD;
  double truncate = 1.0e-6; // stop computing pdf below this limit

  // init grid parameters
  int     n1 = n[0];
  int     n2 = n[1];
  double *x1 = x[0];
  double *x2 = x[1];

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

  // allocate local memory
  deltaA = (double*) calloc(n1, sizeof(double));
  deltaB = (double*) calloc(n2, sizeof(double));
  deltaC = (double*) calloc(IMAX(n1,n2), sizeof(double));

  rho_sq = covar*covar/(var1*var2);
  rho_n  = covar/(var1*var2);
  if (rho_sq != 1.0)
    tmp = -0.5/(1.0-rho_sq);
  else
    tmp = 0.0;

  a[0] = 1.0/(2.0*PI*sqrt(var1*var2*(1.0-rho_sq)));
  a[1] = tmp*(mean1*mean1/var1 - 2.0*mean1*mean2*rho_n + mean2*mean2/var2);
  a[2] = -2.0*tmp*(mean1/var1 - mean2*rho_n);
  a[3] = -2.0*tmp*(mean2/var2 - mean1*rho_n);
  a[4] = tmp/var1;
  a[5] = tmp/var2;
  a[6] = -2.0*tmp*rho_n;

  // sum over the entire grid
  sumPDF = 0.0;
  for (k=0; k<n2; k++)
    {
      if (k == 0)
        {
          dx2 = 0.5*(x2[1]-x2[0]);
          x2c = 0.5*dx2;
        }
      else if (k == n2-1)
        {
          dx2 = 0.5*(x2[k] - x2[k-1]);
          x2c = x2[k] - 0.5*dx2;
        }
      else
        {
          dx2 = 0.5*(x2[k+1] - x2[k-1]);
          x2c = x2[k];
        }

      // start loop over columns
      for (j=0; j<n1; j++)
        {
          if (j == 0)
            {
              dx1 = 0.5*(x1[1]-x1[0]);
              x1c = 0.5*dx1;
            }
          else if (j == n1-1)
            {
              dx1 = 0.5*(x1[j] - x1[j-1]);
              x1c = x1[j] - 0.5*dx1;
            }
          else if ((j == n1 - 1 - k) && (transform == 0))
            {
              // if not transformed space, take into account diagonal boundary
              dx1 = 0.5*(x1[j]-x1[j-1]);
              dx2 = 0.5*(x2[k]-x2[k-1]);
              x1c = x1[j] - 0.5*dx1;
              x2c = x2[k] - 0.5*dx2;
            }
          else
            {
              dx1 = 0.5*(x1[j+1] - x1[j-1]);
              x1c = x1[j];
            }

          // get current local values of Z1 and Z2
          if (transform == 0)
            {
              Z1 = x1c;
              Z2 = x2c;
              if (Z1 + Z2 <= 1.0)
                dA = dx1*dx2;
              else
                dA = 0.0;
            }
          else if (transform == 1)
            {
              Z1 = x1c;
              Z2 = x2c*(1-x1c);
              dA = dx1*dx2*(1-x1c);
            }
          else if (transform == 2)
            {
              Z1 = x1c*(1.0 - x2c);
              Z2 = x1c*x2c;
              dA = dx1*dx2*x1c;
            }

          P = gaussDist2D(Z1, Z2, a[0], a[1], a[2], a[3], a[4], a[5], a[6]);

          pdf[k][j] = P*dA;

          sumPDF += pdf[k][j];
        } // end loop over columns
    } // end loop over rows

  // get the delta function at lower bounds
  double dZ = 0.01;

  // clipping for everything below Z1 axis (Z2 = )
  for (j=0; j<n1; j++)
    {
      if (j == 0)
        {
          dx1 = 0.5*(x1[1]-x1[0]);
          x1c = 0.5*dx1;
        }
      else if (j == n1-1)
        {
          dx1 = 0.5*(x1[j] - x1[j-1]);
          x1c = x1[j] - 0.5*dx1;
        }
      else
        {
          dx1 = 0.5*(x1[j+1] - x1[j-1]);
          x1c = x1[j];
        }

      Z1 = x1c;
      Z2 = 0.0;
      dA = dZ*dx1;
      deltaA[j] = 0.0;
      do {
        Z2 -= dZ;
        P   = gaussDist2D(Z1, Z2, a[0], a[1], a[2], a[3], a[4], a[5], a[6])*dA;
        deltaA[j] += P;
      } while (P > truncate);

      // add the delta function to the boundary
      pdf[0][j] += deltaA[j];
      sumPDF    += deltaA[j];
    }

  // clipping for everything to the left of the Z2 axis 
  for (k=0; k<n2; k++)
    {
      if (k == 0)
        {
          dx2 = 0.5*(x2[1]-x2[0]);
          x2c = 0.5*dx2;
        }
      else if (k == n2-1)
        {
          dx2 = 0.5*(x2[k] - x2[k-1]);
          x2c = x2[k] - 0.5*dx2;
        }
      else
        {
          dx2 = 0.5*(x2[k+1] - x2[k-1]);
          x2c = x2[k];
        }

      Z1 = 0.0;
      Z2 = x2c;
      dA = dZ*dx2;
      deltaB[k] = 0.0;
      do {
        Z1 -= dZ;
        P   = gaussDist2D(Z1, Z2, a[0], a[1], a[2], a[3], a[4], a[5], a[6])*dA;
        deltaB[k] += P;
      } while (P > truncate);

      // add the delta function to the boundary
      pdf[k][0] += deltaB[k];
      sumPDF    += deltaB[k];
    }

  // normalise the pdf
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      pdf[k][j] /= sumPDF;
  
  // set bounds of pdf
  for (k=0; k<n2; k++)
    {
      pdfBound[k][0] = 0;
      pdfBound[k][1] = n1-1;
    }

  free(deltaA);
  free(deltaB);
  free(deltaC);

  return -1;
}
/* -------------------------------------------------------------------------- */

int calcLogNormPDF2D(int *n,double **x, double mean1, double var1, double mean2, 
                     double var2, double var, double **pdf, int **pdfBound,
                     int transform)
{
  double PI = 3.1415926535897931;
  int    i, j, k, check, flag, pdfFlag;
  double Z1, Z2, x1c, x2c;   // grid variables
  double dx1, dx2, dA;       // grid spacing
  double covar;              // covariance
  double a[7];               // array containing parameters distribution
  double sumPDF, P;          // sum of PDF, local probability
  double rho_n, rho_sq, tmp; // tmp vars for computing 2D Gaussian
  double *deltaA, *deltaB, *deltaC, *deltaD;
  double truncate = 1.0e-6; // stop computing pdf below this limit
  double Lmean1, Lmean2, Lvar1, Lvar2, Lvar12;

  // init grid parameters
  int     n1 = n[0];
  int     n2 = n[1];
  double *x1 = x[0];
  double *x2 = x[1];

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

  // allocate local memory
  deltaA = (double*) calloc(n1, sizeof(double));
  deltaB = (double*) calloc(n2, sizeof(double));
  deltaC = (double*) calloc(IMAX(n1,n2), sizeof(double));

  // compute the adjusted means and variances for the log normal
  Lvar1  = log(var1/(mean1*mean1) + 1.0);
  Lvar2  = log(var2/(mean2*mean2) + 1.0);
  Lvar12 = log(covar/(mean1*mean2) + 1.0);
  Lmean1 = log(mean1) - 0.5*Lvar1;
  Lmean2 = log(mean2) - 0.5*Lvar2;

  rho_sq = Lvar12*Lvar12/(Lvar1*Lvar2);
  rho_n  = Lvar12/(Lvar1*Lvar2);
  if (rho_sq != 1.0)
    tmp = -0.5/(1.0-rho_sq);
  else
    tmp = 0.0;

  a[0] = 1.0/(2.0*PI*sqrt(Lvar1*Lvar2*(1.0-rho_sq)));
  a[1] = tmp*(Lmean1*Lmean1/Lvar1-2.0*Lmean1*Lmean2*rho_n+Lmean2*Lmean2/Lvar2);
  a[2] = -2.0*tmp*(Lmean1/Lvar1 - Lmean2*rho_n);
  a[3] = -2.0*tmp*(Lmean2/Lvar2 - Lmean1*rho_n);
  a[4] = tmp/Lvar1;
  a[5] = tmp/Lvar2;
  a[6] = -2.0*tmp*rho_n;

  // sum over the entire grid
  sumPDF = 0.0;
  for (k=0; k<n2; k++)
    {
      if (k == 0)
        {
          dx2 = 0.5*(x2[1]-x2[0]);
          x2c = 0.5*dx2;
        }
      else if (k == n2-1)
        {
          dx2 = 0.5*(x2[k] - x2[k-1]);
          x2c = x2[k] - 0.5*dx2;
        }
      else
        {
          dx2 = 0.5*(x2[k+1] - x2[k-1]);
          x2c = x2[k];
        }

      // start loop over columns
      for (j=0; j<n1; j++)
        {
          if (j == 0)
            {
              dx1 = 0.5*(x1[1]-x1[0]);
              x1c = 0.5*dx1;
            }
          else if (j == n1-1)
            {
              dx1 = 0.5*(x1[j] - x1[j-1]);
              x1c = x1[j] - 0.5*dx1;
            }
          else if ((j == n1 - 1 - k) && (transform == 0))
            {
              // if not transformed space, take into account diagonal boundary
              dx1 = 0.5*(x1[j]-x1[j-1]);
              dx2 = 0.5*(x2[k]-x2[k-1]);
              x1c = x1[j] - 0.5*dx1;
              x2c = x2[k] - 0.5*dx2;
            }
          else
            {
              dx1 = 0.5*(x1[j+1] - x1[j-1]);
              x1c = x1[j];
            }

          // get current local values of Z1 and Z2
          if (transform == 0)
            {
              Z1 = x1c;
              Z2 = x2c;
              if (Z1 + Z2 <= 1.0)
                dA = dx1*dx2;
              else
                dA = 0.0;
            }
          else if (transform == 1)
            {
              Z1 = x1c;
              Z2 = x2c*(1-x1c);
              dA = dx1*dx2*(1-x1c);
            }
          else if (transform == 2)
            {
              Z1 = x1c*(1.0 - x2c);
              Z2 = x1c*x2c;
              dA = dx1*dx2*x1c;
            }

          P = gaussDist2D(log(Z1),log(Z2),a[0],a[1],a[2],a[3],a[4],a[5],a[6]);

          pdf[k][j] =  1.0/(Z1*Z2)*P*dA;

          sumPDF += pdf[k][j];
        } // end loop over columns
    } // end loop over rows

  // normalise the pdf
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      pdf[k][j] /= sumPDF;
  
  // set bounds of pdf
  for (k=0; k<n2; k++)
    {
      pdfBound[k][0] = 0;
      pdfBound[k][1] = n1-1;
    }

  free(deltaA);
  free(deltaB);
  free(deltaC);

  return -1;
}
/* -------------------------------------------------------------------------- */

double gaussDist1D(double x, double a0, double a1, double a2, double a3)
{	
  return a0*exp(a1 + a2*x + a3*x*x);
}
/* ------------------------------------------------------------------------- */

double gaussDist2D(double x1, double x2, double a0, double a1, double a2,
                   double a3, double a4, double a5, double a6)
{	
  return a0*exp(a1 + a2*x1 + a3*x2 + a4*x1*x1 + a5*x2*x2 + a6*x1*x2);
}
/* ------------------------------------------------------------------------- */
