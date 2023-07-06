#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "statsUtils.h"
#include "utils.h"
#include "polyInt.h"

void clipMoments1D(double *mean, double *var)
{
  *mean  = fmin(1.0, fmax(0.0, *mean));
  *var   = fmax(0.0, fmin(*var, *mean*(1.0 - *mean)));
}
/* -------------------------------------------------------------------------- */

void clipMoments2D(double *mean1, double *var1, double *mean2, double *var2,
                   double *var, double *covar)
{
  double mean, var1Max, var2Max, covarMax;

  // clip means
  *mean1 = fmin(1.0, fmax(0.0, *mean1));
  *mean2 = fmin(1.0, fmax(0.0, *mean2));
  mean   = *mean1 + *mean2;

  // clip variances
  var1Max  = *mean1*(1.0 - *mean1);
  var2Max  = *mean2*(1.0 - *mean2);
  *var1    = fmax(0.0, fmin(*var1, var1Max));
  *var2    = fmax(0.0, fmin(*var2, var2Max));
  *var     = fmax(0.0, fmin(*var,  var1Max + var2Max));
  *covar   = 0.5*(*var - *var1 - *var2);
  covarMax = 0.98*sqrt((*var1)*(*var2));
  *covar  = fmax(-covarMax, fmin(*covar, covarMax));
}
/* -------------------------------------------------------------------------- */

void getPDFrowBounds(int n, double *x, double *bounds, double tol)
{
  int j = 0;
  while (x[j] < tol) { j++; };
  bounds[0] = j;

  j = n-1;
  while (x[j] < tol) { j--; };
  bounds[1] = j; 
}
/* -------------------------------------------------------------------------- */

int  checkPDF_1D(int n, double *x, double mean, double var, double *pdf,
                 int *pdfBound)
{
  int    j0, j1, jMean;
  double dx, sigma, sumPDF;

  // delta at x = 0
  if (mean <= 1.0e-10)
    {
      pdf[0] = 1.0;
      pdfBound[0] = pdfBound[1] = 0;
      return 1;
    }

  // delta at x = 1
  if (mean >= (1.0 - 1.0e-10) )
    {
      pdf[n-1] = 1.0;
      pdfBound[0] = pdfBound[1] = n-1;
      return 1;
    }

  // double delta function
  if (var >= mean*(1.0-mean))
    {
      pdf[0]    = 1.0 - mean;
      pdf[n-1] = mean;
      pdfBound[0] = 0;
      pdfBound[1] = n-1;
      return 1;
    }

  // find the location of the mesh point before the mean
  jMean = -1;
  do { jMean++; } while (x[jMean] < mean);
  jMean--;

  sigma = sqrt(var);    // compute the standard deviation

  /* if the mesh size about the mean is less than 1.5 standard deviations
     compute at two points about the mean such that the mean is conserved */
  if (x[jMean+1] - x[jMean] > 1.5*sigma)
    {
      j0 = jMean;
      j1 = jMean + 1;
               
      dx = (x[j1]-x[j0]);
      pdf[j0] = fabs(x[j1] - mean)/dx;
      pdf[j1] = fabs(mean - x[j0])/dx;
                        
      // normalise the PDF
      sumPDF = pdf[j0]+pdf[j1];
      pdf[j0] = pdf[j0]/sumPDF;
      pdf[j1] = pdf[j1]/sumPDF;
                        
      // store bounds of non-zero PDF values
      pdfBound[0] = j0;
      pdfBound[1] = j1;
      return 1;
    }

  return 0; // not a special case
}
/* -------------------------------------------------------------------------- */

int  checkPDF_2D(int *n, double **x, double mean1, double var1, double mean2, 
                 double var2, double var, double **pdf, int **pdfBound,
                 int transform, int method)
{
  int    i, j, k, jM, kM, pdfFlag;
  double sumPDF;
  double x1M, x2M;        // the values of x1 and x2 for given means
  double dx1, dx2, x1min, x1max, x2min, x2max;
  double sigma1, sigma2;
  double sigma;           // standard deviation
  double *pdfTmp;
  int    pdfBoundTmp[2];

  // set up the grid parameters
  int    n1  = n[0];
  int    n2  = n[1];
  double *x1 = x[0];
  double *x2 = x[1];

  // set the function pointer to the appropriate 1D pdf
  if (method == 0)
    calcPDF_1D = calcBetaPDF;
  else if (method == 1)
    calcPDF_1D = calcSML_PDF;

  // 1D profile along abscissa (x1)
  if (mean2 <= 1.0e-10)
    {
      // compute the 1D pdf
      // NOTE: special 1D cases are handled in function
      pdfFlag = (*calcPDF_1D)(n1, x1, mean1, var1, pdf[0], pdfBound[0]);

      // set the pdf bounds for the rest of the domain
      for (k=1; k<n2-1; k++)
        pdfBound[k][0] = pdfBound[k][1] = 0;

      if (pdfFlag == 0)
        return 1;
      else if (pdfFlag == 1)
        return 11;
    }

  // 1D profile along ordinate (x2)
  if (mean1 <= 1.0e-10)
    {
      // copy the solution into the first column and set bounds
      if (transform == 1)
        {
          pdfTmp = (double*) calloc(n2, sizeof(double));
          pdfFlag = (*calcPDF_1D)(n2, x2, mean2, var2, pdfTmp, pdfBoundTmp);

          for (k=0; k<n2; k++)
            {
              pdf[k][0] = pdfTmp[k];
              pdfBound[k][0] = pdfBound[k][1] = 0;
            }
          free (pdfTmp);
        }
      else if (transform == 2)
        {
          pdfFlag = (*calcPDF_1D)(n1,x1,mean2,var2,pdf[n2-1],pdfBound[n2-1]);

          for (k=0; k<n2-1; k++)
            pdfBound[k][0] = pdfBound[k][1] = 0;
        }

      if (pdfFlag == 0)
        return 2;
      else if (pdfFlag == 1)
        return 21;
    }

  // 1D profile along Z1+Z2 = 1
  if ((1.0-mean1-mean2) <= 1.0e-10)
    {
      if (transform == 1)
        {
          double meanFuel2;

          meanFuel2 = sqrt(mean1*mean1+(1.0-mean2)*(1.0-mean2));
          // compute the 1D pdf
          pdfFlag = (*calcPDF_1D)(n1, x1, mean1, 0.5*(var1+var2+var),
                                  pdf[n2-1],pdfBound[n2-1]);

          // set the pdf bounds for the rest of the domain
          for (k=0; k<n2-1; k++)
            pdfBound[k][0] = pdfBound[k][1] = 0;
        }
      else if (transform == 2)
        {
          pdfTmp = (double*) calloc(n2, sizeof(double));

          // not sure this is right
          pdfFlag = (*calcPDF_1D)(n2, x2, mean2, 0.5*(var1+var2+var),
                                  pdfTmp, pdfBoundTmp);

          for (k=0; k<n2-1; k++)
            {
              pdf[k][n1-1] = pdfTmp[k];
              pdfBound[k][0] = pdfBound[k][1] = n1-1;
            }
          free(pdfTmp);
        }

      if (pdfFlag == 0)
        return 4;
      else if (pdfFlag == 1)
        return 41;
    }

  // set the appropriate location based on transform
  if (transform == 0)
    {
      x1M = mean1;
      x2M = mean2;
      sigma1 = sqrt(var1);
      sigma2 = sqrt(var2);
    }
  else if (transform == 1)
    {
      x1M = mean1;
      x2M = mean2/(1.0 - x1M);
      sigma1 = sqrt(var1);
      sigma2 = sqrt(var2);
    }
  else if (transform == 2)
    {
      x1M = mean1 + mean2;
      x2M = mean2/x1M;
      sigma1 = sqrt(var);
      sigma2 = sqrt(var2);
    }

  // find the location of the mesh point before the mean
  jM  = gridLocate(n1, x1, x1M);
  kM  = gridLocate(n2, x2, x2M);
  dx1 = x1[jM+1] - x1[jM];
  dx2 = x2[kM+1] - x2[kM];

  double delSig = 2.0;
  double Atot, A00, A01, A10, A11; /* for area weighting */
  double *x1_c, *x2_c;
  gPoint aP[4], bP[4];

  // determine whether grid resolution is sufficient
  if ( dx1 > delSig*sigma1 )      /* line solution along x2 */
    {
      if ( dx2 > delSig*sigma2 )  /* delta function*/
        {
          /* add weighted distribution to each of the four points
             surrounding the value */
          /* compute the relative areas */
          Atot = (x1[jM+1] - x1[jM])*(x2[kM+1] - x2[kM]); /* dx*dy */
          A00  = (x1[jM+1] - x1M   )*(x2[kM+1] - x2M   );
          A01  = (x1[jM+1] - x1M   )*(x2M      - x2[kM]);
          A10  = (x1M      - x1[jM])*(x2[kM+1] - x2M   );
          A11  = (x1M      - x1[jM])*(x2M      - x2[kM]);

          // distribute pdf to four surrounding points based on weighted areas
          pdf[kM][jM]     = A00/Atot; // bottom left corner 
          pdf[kM+1][jM]   = A01/Atot; // top left corner
          pdf[kM][jM+1]   = A10/Atot; // bottom right corner
          pdf[kM+1][jM+1] = A11/Atot; // top right corner

          // normalise the distribution
          sumPDF = pdf[kM][jM]+pdf[kM+1][jM]+pdf[kM][jM+1]+pdf[kM+1][jM+1];
          pdf[kM][jM]     /= sumPDF;
          pdf[kM+1][jM]   /= sumPDF;
          pdf[kM][jM+1]   /= sumPDF;
          pdf[kM+1][jM+1] /= sumPDF;

          // set bounds for all rows to zero
          for (k=0; k<n2; k++)
            pdfBound[k][0] = pdfBound[k][1] = 0;

          // set the bounds to just compute 4 points
          pdfBound[kM][0] = pdfBound[kM+1][0] = jM;
          pdfBound[kM][1] = pdfBound[kM+1][1] = jM+1;

          return 5;
        }
      else
        {
          x2max   = 1.0-mean1;
          pdfTmp  = (double*) calloc(n2, sizeof(double));
          pdfFlag = (*calcPDF_1D)(n2, x2, mean2/x2max, var2/(x2max*x2max),
                                  pdfTmp, pdfBoundTmp);

          for (k=0; k<n2; k++)
            for (j=0; j<n1; j++)
              pdf[k][j] = 0.0;

          /* set the thickness to be the rms */
          x1min = mean1 - 0.5*delSig*(sigma1 + 1.0e-6);
          x1max = mean1 + 0.5*delSig*(sigma1 + 1.0e-6);

          /* compute the midpoints and the grid spacing for the grid */
          x1_c   = (double*) calloc(n1+1, sizeof(double));
          x1_c[0] = x1[0];
          for (j=0; j<n1-1; j++) x1_c[j+1] = 0.5*(x1[j+1]+x1[j]);
          x1_c[n1] = x1[n1-1];

          x2_c = (double*) calloc(n2+1, sizeof(double));
          x2_c[0] = x2[0];
          for (k=0; k<n2-1; k++) x2_c[k+1] = 0.5*(x2[k+1]+x2[k]);
          x2_c[n2] = x2[n2-1];

          for (k=0; k<n2; k++)
            {
              /* get the grid point from the mean mixture fractions */
              /* use clockwise order */
              getGridfromZ(transform,x1min,x2max*x2_c[k],  &aP[0].x, &aP[0].y);
              getGridfromZ(transform,x1min,x2max*x2_c[k+1],&aP[1].x, &aP[1].y);
              getGridfromZ(transform,x1max,x2max*x2_c[k+1],&aP[2].x, &aP[2].y);
              getGridfromZ(transform,x1max,x2max*x2_c[k],  &aP[3].x, &aP[3].y);

              /* for testing, loop over entire other grid */
              Atot = polyIntersect(4, aP, 4, aP);
              for (i=0; i<n2; i++)
                {
                  bP[0].y = bP[3].y = x2_c[i];
                  bP[1].y = bP[2].y = x2_c[i+1];

                  for (j=0; j<n1; j++)
                    {
                      bP[0].x = bP[1].x = x1_c[j];
                      bP[2].x = bP[3].x = x1_c[j+1];

                      if ( x1_c[j+1] > fmin(aP[0].x, aP[1].x) &&
                           x2_c[i+1] > fmin(aP[0].y, aP[3].y) )
                        {
                          if (x1_c[j] > fmax(aP[2].x, aP[3].x))
                            break;

                          /* compute the intersection area */
                          A00 = polyIntersect(4, aP, 4, bP);

                          /* add the probability from 1 cell to this */
                          pdf[i][j] += A00/Atot*pdfTmp[k];
                        }
                    }
                  if (x2_c[i] > fmax(aP[1].y, aP[2].y))
                    break;
                }
            } /* end loop over rows */

          sumPDF = 0.0;
          for (k=0; k<n2; k++)
            {
              for (j=0; j<n1; j++)
                {
                  sumPDF += pdf[k][j];
                }
              pdfBound[k][0] = 0;
              pdfBound[k][1] = n1-1;
            }

          for (k=0; k<n2; k++)
            for (j=0; j<n1; j++)
              pdf[k][j] /= sumPDF;

          free(pdfTmp); free(x1_c); free(x2_c);
          return 202;
        }
    }
  else if ( dx2 > 1.5*sigma2 ) /* line solution along x1 */
    {
      pdfTmp = (double*) calloc(n1, sizeof(double));
      x1max = 1.0-mean2;
      pdfFlag = (*calcPDF_1D)(n1, x1, mean1/x1max, var1/(x1max*x1max),
                              pdfTmp, pdfBoundTmp);

      for (k=0; k<n2; k++)
        for (j=0; j<n1; j++)
          pdf[k][j] = 0.0;

      /* set the thickness to be the rms */
      x2min = mean2 - 0.5*delSig*(sigma2 + 1.0e-6);
      x2max = mean2 + 0.5*delSig*(sigma2 + 1.0e-6);

      /* compute the midpoints and the grid spacing for the grid */
      x1_c   = (double*) calloc(n1+1, sizeof(double));
      x1_c[0] = x1[0];
      for (j=0; j<n1-1; j++) x1_c[j+1] = 0.5*(x1[j+1]+x1[j]);
      x1_c[n1] = x1[n1-1];

      x2_c = (double*) calloc(n2+1, sizeof(double));
      x2_c[0] = x2[0];
      for (k=0; k<n2-1; k++) x2_c[k+1] = 0.5*(x2[k+1]+x2[k]);
      x2_c[n2] = x2[n2-1];

      for (j=0; j<n1; j++)
        {
          /* get the grid point from the mean mixture fractions */
          /* use clockwise order */
          getGridfromZ(transform, x1max*x1_c[j],   x2min, &aP[0].x, &aP[0].y);
          getGridfromZ(transform, x1max*x1_c[j],   x2max, &aP[1].x, &aP[1].y);
          getGridfromZ(transform, x1max*x1_c[j+1], x2max, &aP[2].x, &aP[2].y);
          getGridfromZ(transform, x1max*x1_c[j+1], x2min, &aP[3].x, &aP[3].y);

          /* for testing, loop over entire other grid */
          Atot = polyIntersect(4, aP, 4, aP);
          for (k=0; k<n2; k++)
            {
              bP[0].y = bP[3].y = x2_c[k];
              bP[1].y = bP[2].y = x2_c[k+1];

              for (i=0; i<n1; i++)
                {
                  bP[0].x = bP[1].x = x1_c[i];
                  bP[2].x = bP[3].x = x1_c[i+1];

                  if ( x1_c[i+1] > fmin(aP[0].x, aP[1].x) &&
                       x2_c[k+1] > fmin(aP[0].y, aP[3].y) )
                    {
                      if (x1_c[i] > fmax(aP[2].x, aP[3].x))
                        break;

                      /* compute the intersection area */
                      A00 = polyIntersect(4, aP, 4, bP);

                      /* add the probability from 1 cell to this */
                      pdf[k][i] += A00/Atot*pdfTmp[j];
                    }
                }
              if (x2_c[k] > fmax(aP[1].y, aP[2].y))
                break;
            }
        } /* end loop over rows */

      sumPDF = 0.0;
      for (k=0; k<n2; k++)
        {
          for (j=0; j<n1; j++)
            {
              sumPDF += pdf[k][j];
            }
          pdfBound[k][0] = 0;
          pdfBound[k][1] = n1-1;
        }

      for (k=0; k<n2; k++)
        for (j=0; j<n1; j++)
          pdf[k][j] /= sumPDF;

      free(pdfTmp); free(x1_c); free(x2_c);
      return 102;
    }
  return 0;
}
/* -------------------------------------------------------------------------- */

double calcLinMeanY(int index, int offset, double alpha, double beta, int n, 
                    double *x, double *Y, double *pdf, int *bounds, int pdfFlag)
{
  int j, k;
  double Ymean = 0.0;
  double dx0 = 0.5*(x[1]-x[0]);
  double dxN = 0.5*(x[n-1]-x[n-2]);
  double fgam;

  if (pdfFlag == 0)
    fgam = lnGamma(alpha + beta) - lnGamma(alpha) - lnGamma(beta);

  for (j=bounds[0]; j<=bounds[1]; j++)
    {
      Ymean += Y[j*offset+index]*pdf[j];

      if (pdfFlag == 0)
        {
          if (j == 0)
            Ymean += (Y[offset+index]-Y[index])/(x[1]-x[0])*
              exp((alpha+1)*log(dx0) + fgam)/(alpha+1);
          else if (j == n-1)
            Ymean += (Y[(j-1)*offset+index]-Y[j*offset+index])/(2*dxN)*
              exp((beta+1)*log(dxN) + fgam)/(beta+1);
        }
    }
  return Ymean;
}

/* -------------------------------------------------------------------------- */
void calcLinAllMeanY(int nSpecies, double alpha, double beta, int n, 
                    double *x, double *Y, double *Ymean, double *pdf, int *bounds, int pdfFlag)
{
  int i, j, k;
  double dx0 = 0.5*(x[1]-x[0]);
  double dxN = 0.5*(x[n-1]-x[n-2]);
  double fgam;

  for (i=0;i<nSpecies;i++) Ymean[i] = 0.0;

  if (pdfFlag == 0)
    fgam = lnGamma(alpha + beta) - lnGamma(alpha) - lnGamma(beta);

  for (j=bounds[0]; j<=bounds[1]; j++)
    {
      for (i=0;i<nSpecies;i++) Ymean[i] += Y[j*nSpecies+i]*pdf[j];

      if (pdfFlag == 0)
        {
          if (j == 0) {
            for (i=0;i<nSpecies;i++) Ymean[i] += (Y[nSpecies+i]-Y[i])/((x[1]-x[0])*(alpha+1))*
              exp((alpha+1)*log(dx0) + fgam);
	  }
          else if (j == n-1) {
            for (i=0;i<nSpecies;i++) Ymean[i] += (Y[(j-1)*nSpecies+i]-Y[j*nSpecies+i])/((2*dxN)*(beta+1))*
              exp((beta+1)*log(dxN) + fgam);
 	  }
        }
    }
  return;
}
/* -------------------------------------------------------------------------- */

double calcLinMeanX(int index, int offset, double alpha, double beta, int n, 
                    double *x, double *Y, double *pdf, int *bounds, int pdfFlag,
                    int nSpecies, double *W)
{
  int j, k;
  double Ymean = 0.0;
  double dx0 = 0.5*(x[1]-x[0]);
  double dxN = 0.5*(x[n-1]-x[n-2]);
  double fgam, WMix1 = 0.0, WMix2 = 0.0;

  if (pdfFlag == 0)
    fgam = lnGamma(alpha + beta) - lnGamma(alpha) - lnGamma(beta);

  for (j=bounds[0]; j<=bounds[1]; j++)
    {
      k = j*offset;
      WMix1 = calcMixtureW(nSpecies, &Y[k], W);
      Ymean += pdf[j]*Y[k+index]*WMix1/W[index];

      if (pdfFlag == 0)
        {
          if (j == 0)
            {
              WMix2 = calcMixtureW(nSpecies, &Y[(j+1)*offset], W);
              Ymean += (Y[offset+index]*WMix2-Y[index]*WMix1)/(x[1]-x[0])*
                exp((alpha+1)*log(dx0) + fgam)/(alpha+1)/W[index];
            }
          else if (j == n-1)
            {
              WMix2 = calcMixtureW(nSpecies, &Y[(j-1)*offset], W);
              Ymean += (Y[(j-1)*offset+index]*WMix2-Y[j*offset+index]*WMix1)/
                (2*dxN)*exp((beta+1)*log(dxN) + fgam)/(beta+1)/W[index];
            }
        }
    }
  return Ymean;
}
/* -------------------------------------------------------------------------- */

double calcLinMeanY_2D(int index, int offset, int *n, double **x, double *Y,
                       double beta1, double beta2, double beta3, double **pdf,
                       int **pdfBound, int pdfFlag, int transform)
{
  int    j, k, row, col, iX1, iX2;
  int    n1, n2;    // number of points in each grid
  double *x1, *x2;  // grids in each direction
  double x1cur, x2cur;
  double dx1, dx2, tmp, mean;

  // temporary variables used when 1D solution along x2 axis
  int    pdfBoundTmp[2];
  double *pdf1, *pdf2, *Y2;

  // temporary variables for boundary computations
  double fgam;                // beta function (from gammas)
  double F_Y, G_Y, Q_Y;       // for functions of Y
  double B1, B2, B1_1, B2_1;  // for weird partial beta distributions
  double eps0, epsN, eps0_1, epsN_1;
  double b1_0, b1_N, b2_0, b2_N;
  double Y00, Y01, Y10, Y11;

  // assign local grid parameters
  n1 = n[0];
  n2 = n[1];
  x1 = x[0];
  x2 = x[1];

  // half grid sizes at the boundaries
  double dx1_0  = 0.5*(x1[1]-x1[0]);
  double dx1_N  = 0.5*(x1[n1-1]-x1[n1-2]);
  double dx2_0  = 0.5*(x2[1]-x2[0]);
  double dx2_N  = 0.5*(x2[n2-1]-x2[n2-2]);

  // set the indices for boundaries
  iX1 = (n1-1)*offset;      // start of x1 = 1, x2 = 0
  iX2 = (n2-1)*(n1*offset); // start of x1 = 0, x2 = 1

  /* ----- solution along Z1 (Z2 = 0) ----- */
  if (pdfFlag == 1)
    return calcLinMeanY(index, offset, beta1, beta3, n1, x1, &Y[0], pdf[0],
                        pdfBound[0], 0);

  if (pdfFlag == 11)
    return calcLinMeanY(index, offset, beta1, beta3, n1, x1, &Y[0], pdf[0],
                        pdfBound[0], 1);

  /* ----- solution along Z2 (Z1 = 0) ----- */
  if ( (pdfFlag == 2) || (pdfFlag == 21) )
    {
      Y2   = (double*) calloc(n2, sizeof(double));
      pdf2 = (double*) calloc(n2, sizeof(double));
      pdfBoundTmp[0] = 0;
      pdfBoundTmp[1] = n2-1;

      for (k=0; k<n2; k++)
        {
          row = k*n1*offset;
          pdf2[k] = pdf[k][0];
          Y2[k] = Y[row+index];
        }

      if (pdfFlag == 2)
        mean = calcLinMeanY(0, 1, beta2, beta3, n2, x2, Y2, pdf2, pdfBoundTmp,0);
      else if (pdfFlag == 21)
        mean = calcLinMeanY(0, 1, beta2, beta3, n2, x2, Y2, pdf2, pdfBoundTmp,1);

      free(Y2);
      free(pdf2);
      return mean;
    }

  /* ----- solution along Z1+Z2 = 1 ----- */
  if (pdfFlag == 4)
    return calcLinMeanY(index, offset, beta1, beta2, n1, x1, &Y[iX2], pdf[n2-1],
                        pdfBound[n2-1], 0);

  if (pdfFlag == 41)
    return calcLinMeanY(index, offset, beta1, beta2, n1, x1, &Y[iX2], pdf[n2-1],
                        pdfBound[n2-1], 1);

  /* ----- triple delta function ----- */
  if (pdfFlag == 3)
    return Y[0+index]*pdf[0][0] + Y[iX1+index]*pdf[0][n1-1]
      + Y[iX2+index]*pdf[n2-1][0];

  /* ----- delta function ----- */
  if (pdfFlag == 5)
    {
      mean = 0.0;
      for(k=0; k<n2; k++)
        {
          row = k*n1*offset;
          for (j=0; j<n1; j++)
            {
              col = j*offset;
              mean += Y[row+col+index]*pdf[k][j];
            }
        }
      return mean;
    }

  /* ----- extract the specific value from the unstructured array ----- */
  double **Yi = (double**) calloc(n2, sizeof(double*));
  for (k=0; k<n2; k++)
    Yi[k] = (double*) calloc(n1, sizeof(double));

  for(k=0; k<n2; k++)
    {
      row = k*n1*offset;
      for (j=0; j<n1; j++)
        {
          col = j*offset;
          Yi[k][j] = Y[row + col + index];
        } // end loop over columns
    } // end loop over rows

  /* ----- compute the contribution from the center points ----- */
  mean = 0.0;
  for (k=0; k < n2; k++)
    for (j=0; j < n1; j++)
      mean += Yi[k][j]*pdf[k][j];

  // error check; can't use boundary expansion with negative beta
  if (pdfFlag < 0 || pdfFlag > 100 || transform == 0)
    {
      for (k=0; k<n2; k++) free(Yi[k]);
      free(Yi);
      return mean;
    }

  /* ----- add additional terms on boundary assuming linear ----- */
  fgam = lnGamma(beta1+beta2+beta3)-lnGamma(beta1)-lnGamma(beta2)-lnGamma(beta3);

  if (transform == 1)
    {
      b1_0 = beta1;
      b1_N = beta2+beta3;
      b2_0 = beta2;
      b2_N = beta3;
    }
  else if (transform == 2)
    {
      b1_0 = beta1+beta2;
      b1_N = beta3;
      b2_0 = beta2;
      b2_N = beta1;
    }

  // add extra terms to oxidizer corner (0,0)
  double T1, T2, T3;
  dx1 = x1[1]-x1[0];
  dx2 = x2[1]-x2[0];
  F_Y = (Yi[0][1] - Yi[0][0])/dx1;
  G_Y = (Yi[1][0] - Yi[0][0])/dx2;
  Q_Y = (Yi[1][1] - Yi[1][0] - Yi[0][1] + Yi[0][0])/(dx1*dx2);
  T1 = F_Y*exp((b1_0+1)*log(dx1_0) + b2_0*log(dx2_0) + fgam)/((b1_0+1)*b2_0);
  T2 = G_Y*exp(b1_0*log(dx1_0) + (b2_0+1)*log(dx2_0) + fgam)/(b1_0*(b2_0+1));
  T3 = Q_Y*exp((b1_0+1)*log(dx1_0) + (b2_0+1)*log(dx2_0) + fgam)/
    ((b1_0+1)*(b2_0+1));
  mean = mean + T1 + T2 + T3;

  // add extra terms to corner (1,0)
  dx1 = x1[n1-1]-x1[n1-2];
  F_Y = (Yi[0][n1-2] - Yi[0][n1-1])/dx1;
  G_Y = (Yi[1][n1-1] - Yi[0][n1-1])/dx2;
  Q_Y = (Yi[1][n1-2] - Yi[1][n1-1] - Yi[0][n2-2] + Yi[0][n2-1])/(dx1*dx2);
  T1 = F_Y*exp((b1_N+1)*log(dx1_N) + b2_0*log(dx2_0) + fgam)/((b1_N+1)*b2_0);
  T2 = G_Y*exp(b1_N*log(dx1_N) + (b2_0+1)*log(dx2_0) + fgam)/(b1_N*(b2_0+1));
  T3 = Q_Y*exp((b1_N+1)*log(dx1_N) + (b2_0+1)*log(dx2_0) + fgam)/
    ((b1_N+1)*(b2_0+1));
  mean = mean + T1 + T2 + T3;

  // add extra terms to corner (0,1)
  dx1 = x1[1]-x1[0];
  dx2 = x2[n2-1]-x2[n2-2];
  F_Y = (Yi[n2-1][1] - Yi[n2-1][0])/dx1;
  G_Y = (Yi[n2-2][0] - Yi[n2-1][0])/dx2;
  Q_Y = (Yi[n2-2][1] - Yi[n2-2][0] - Yi[n2-1][1] + Yi[n2-1][0])/(dx1*dx2);
  T1 = F_Y*exp((b1_0+1)*log(dx1_0) + b2_N*log(dx1_N) + fgam)/((b1_0+1)*b2_N);
  T2 = G_Y*exp(b1_0*log(dx1_0) + (b2_N+1)*log(dx2_N) + fgam)/(b1_0*(b2_N+1));
  T3 = Q_Y*exp((b1_0+1)*log(dx1_0) + (b2_N+1)*log(dx2_N) + fgam)
    /((b1_0+1)*(b2_N+1));
  mean = mean + T1 + T2 + T3;

  // add extra terms to corner (1,1)
  /*dx1 = x1[n1-1]-x1[n1-2];
  dx2 = x2[n2-1]-x2[n2-2];
  F_Y = (Yi[n2-1][n1-2] - Yi[n2-1][n1-1])/dx1;
  G_Y = (Yi[n2-2][n1-1] - Yi[n2-1][n1-1])/dx2;
  Q_Y = (Yi[n2-2][n1-2] - Yi[n2-2][n1-1] - Yi[n2-1][n1-2] + Yi[n2-1][n1-1])/
    (dx1*dx2);
  T1 = F_Y*exp((b1_N+1)*log(dx1_N) + b2_N*log(dx1_N) + fgam)/((b1_N+1)*b2_N);
  T2 = G_Y*exp(b1_N*log(dx1_N) + (b2_N+1)*log(dx2_N) + fgam)/(b1_N*(b2_N+1));
  T3 = Q_Y*exp((b1_N+1)*log(dx1_N) + (b2_N+1)*log(dx2_N) + fgam)
    /((b1_N+1)*(b2_N+1));
    mean = mean + T1 + T2 + T3; */
 
  /*
  // add additional terms to left boundary
  eps0   = b1_0*log(dx1_0);
  eps0_1 = (b1_0+1.0)*log(dx1_0);
  dx1 = x1[1]-x1[0];
  for (k=1; k<n2-1; k++)
    {
      dx2 = 0.5*(x2[k+1]-x2[k-1]);
      B2   = (b2_0-1.0)*log(x2[k]) + (b2_N-1.0)*log(1.0-x2[k]);
      B2_1 = b2_0*log(x2[k]) + (b2_N-1.0)*log(1.0-x2[k]);

      Y00 = 0.5*(Yi[k][0] + Yi[k-1][0]);
      Y01 = 0.5*(Yi[k][1] + Yi[k-1][1]);
      Y10 = 0.5*(Yi[k][0] + Yi[k+1][0]);
      Y11 = 0.5*(Yi[k][1] + Yi[k+1][1]);

      F_Y = (Y01 - Y00)/dx1;
      G_Y = (Y10 - Y00)/dx2;
      Q_Y = (Y11 - Y10 - Y01 + Y00)/(dx1*dx2);

      T1 = F_Y*exp(eps0_1 + B2   + fgam)*dx2/(b1_0+1.0);
      T2 = G_Y*exp(eps0   + B2_1 + fgam)*dx2/(b1_0);
      T3 = Q_Y*exp(eps0_1 + B2_1 + fgam)*dx2/(b1_0+1.0);

      mean = mean + T1 + T2 + T3;
    }

  // add additional terms to right boundary
  epsN   = b1_N*log(dx1_N);
  epsN_1 = (b1_N+1.0)*log(dx1_N);
  dx1 = x1[n1-1]-x1[n1-2];
  for (k=1; k<n2-1; k++)
    {
      dx2 = 0.5*(x2[k+1]-x2[k-1]);
      B2   = (b2_0-1.0)*log(x2[k]) + (b2_N-1.0)*log(1.0-x2[k]);
      B2_1 = b2_0*log(x2[k]) + (b2_N-1.0)*log(1.0-x2[k]);

      Y00 = 0.5*(Yi[k][n1-1] + Yi[k-1][n1-1]);
      Y01 = 0.5*(Yi[k][n1-2] + Yi[k-1][n1-2]);
      Y10 = 0.5*(Yi[k][n1-1] + Yi[k+1][n1-1]);
      Y11 = 0.5*(Yi[k][n1-2] + Yi[k+1][n1-2]);

      F_Y = (Y01 - Y00)/dx1;
      G_Y = (Y10 - Y00)/dx2;
      Q_Y = (Y11 - Y10 - Y01 + Y00)/(dx1*dx2);

      T1 = F_Y*exp(epsN_1 + B2   + fgam)*dx2/(b1_N+1.0);
      T2 = G_Y*exp(epsN   + B2_1 + fgam)*dx2/(b1_N);
      T3 = Q_Y*exp(epsN_1 + B2_1 + fgam)*dx2/(b1_N+1.0);

      mean = mean + T1 + T2 + T3;
    }

  // bottom boundary
  dx2 = x2[k]-x2[0];
  eps0   = b2_0*log(dx2_0);
  eps0_1 = (b2_0+1.0)*log(dx2_0);
  for (j=1; j<n1-1; j++)
    {
      dx1 = 0.5*(x1[j+1]-x1[j-1]);
      B1   = (b1_0-1.0)*log(x1[j])+(b1_N-1.0)*log(1.0-x1[j]);
      B1_1 = b1_0*log(x1[j]) + (b2_N-1.0)*log(1.0-x1[j]);

      Y00 = 0.5*(Yi[0][j] + Yi[0][j-1]);
      Y01 = 0.5*(Yi[0][j] + Yi[0][j+1]);
      Y10 = 0.5*(Yi[1][j] + Yi[1][j-1]);
      Y11 = 0.5*(Yi[1][j] + Yi[1][j+1]);

      F_Y = (Y01 - Y00)/dx1;
      G_Y = (Y10 - Y00)/dx2;
      Q_Y = (Y11 - Y10 - Y01 + Y00)/(dx1*dx2);

      T1 = F_Y*exp(eps0   + B1_1 + fgam)*dx1/(b2_0);
      T2 = G_Y*exp(eps0_1 + B1   + fgam)*dx1/(b2_0+1.0);
      T3 = Q_Y*exp(eps0_1 + B1_1 + fgam)*dx1/(b2_0+1.0);

      mean = mean + T1 + T2 + T3; 
    }
  
  // top boundary
  dx2    = x2[n2-1]-x2[n2-2];
  epsN   = b2_N*log(dx2_N);
  epsN_1 = (b2_N+1.0)*log(dx2_N);
  for (j=1; j<n1-1; j++)
    {
      dx1 = 0.5*(x1[j+1]-x1[j-1]);
      B1   = (b1_0-1.0)*log(x1[j])+(b1_N-1.0)*log(1.0-x1[j]);
      B1_1 = b1_0*log(x1[j]) + (b1_N-1.0)*log(1.0-x1[j]);

      Y00 = 0.5*(Yi[n2-1][j] + Yi[n2-1][j-1]);
      Y01 = 0.5*(Yi[n2-1][j] + Yi[n2-1][j+1]);
      Y10 = 0.5*(Yi[n2-2][j] + Yi[n2-2][j-1]);
      Y11 = 0.5*(Yi[n2-2][j] + Yi[n2-2][j+1]);

      F_Y = (Y01 - Y00)/dx1;
      G_Y = (Y10 - Y00)/dx2;
      Q_Y = (Y11 - Y10 - Y01 + Y00)/(dx1*dx2);

      T1 = F_Y*exp(epsN   + B1_1 + fgam)*dx1/(b2_N);
      T2 = G_Y*exp(epsN_1 + B1   + fgam)*dx1/(b2_N+1.0);
      T3 = Q_Y*exp(epsN_1 + B1_1 + fgam)*dx1/(b2_N+1.0);

      mean = mean + T1 + T2 + T3; 
    }
  */
  for (k=0; k<n2; k++) free(Yi[k]);
  free(Yi);

  return mean;
}
/* -------------------------------------------------------------------------- */

double calcVar(int *n, int index, int offset, double mean1, double mean2, 
               double *Y1, double *Y2, double **pdf)
{
  int    j, k, row, col;
  double var = 0.0;

  for (k=0; k<n[1]; k++)
    {
      row = k*n[0]*offset;
      for (j=0; j<n[0]; j++)
        {
          col = j*offset;
          var += (Y1[row+col+index]-mean1)*(Y2[row+col+index]-mean2)*pdf[k][j];
        } // end loop over columns
    } // end loop over rows

  return var;
}
/* -------------------------------------------------------------------------- */

void calcMarginalPDF(int dir, int *n, double **x, double *marginal, double **pdf)
{
  int i, j, k, row, col;

  // set the marginal to zero
  for (i=0; i<n[dir-1]; i++) marginal[i] = 0.0;

  for (k=0; k<n[1]; k++)
    {
      for (j=0; j<n[0]; j++)
        {
          if (dir == 1)
            marginal[j] += pdf[k][j];
          else if (dir == 2)
            marginal[k] += pdf[k][j];
        } // end loop over columns
    } // end loop over rows
}
/* -------------------------------------------------------------------------- */

double calcKullbackLeibler(int n1, int n2, double **P, double **Q)
{
  double KL = 0.0;
  int j, k;
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++)
        {
          if (P[k][j] > 0.0 && Q[k][j] > 0.0)
            KL += P[k][j]*log(P[k][j]/Q[k][j]);
        } // end loop over columns
    } // end loop over rows
  return KL;
}
/* -------------------------------------------------------------------------- */

double calcJensonShannon(int n1, int n2, double **P, double **Q)
{
  double M, JS = 0.0;
  int j, k;
  double tiny = 1.0e-6;
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++)
        {
          M = 0.5*(P[k][j] + Q[k][j]); // average distribution

          if (P[k][j] > 0.0 && M > 0.0) // D_KL(P||M)
            JS += 0.5*P[k][j]*log(P[k][j]/M);

          if (Q[k][j] > 0.0 && M > 0.0) // D_KL(Q||M)
            JS += 0.5*Q[k][j]*log(Q[k][j]/M);
        } // end loop over columns
    } // end loop over rows
  return JS;
}
/* -------------------------------------------------------------------------- */

double calcHellinger(int n1, int n2, double **P, double **Q)
{
  double tmp, H = 0.0;
  int j, k;
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++)
        {
          tmp = sqrt(P[k][j]) - sqrt(Q[k][j]);
          H  += tmp*tmp;
        } // end loop over columns
    } // end loop over rows
  return 0.5*H;
}
/* -------------------------------------------------------------------------- */

double calcChiSqDiv(int n1, int n2, double **P, double **Q)
{
  double chi = 0.0;
  int j, k;
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++)
        {
          if (P[k][j] > 0.0)
            chi  += Q[k][j]*Q[k][j]/P[k][j];
        } // end loop over columns
    } // end loop over rows
  return chi - 1;
}
/* -------------------------------------------------------------------------- */

double calcTVD(int n1, int n2, double **P, double **Q)
{
  double TVD = 0.0;
  int j, k;
  double tiny = 1.0e-6;
  for (k=0; k<n2; k++)
    {
      for (j=0; j<n1; j++)
        {
          TVD += fabs(P[k][j] - Q[k][j]);
        } // end loop over columns
    } // end loop over rows
  return 0.5*TVD;
}
/* -------------------------------------------------------------------------- */

int binPDF_1D(int nx, int ny, int nz, double ***var, int n1, double *x1,
              double *pdf)
{
  int i, j, k, jL;
  double plusPDF;

  plusPDF = 1.0/(double)(nx*ny*nz);

  // initialise to zero
  for (j=0; j<n1; j++) pdf[j] = 0.0;

  // loop over all data points
  for (k=0; k<nz; k++)
    {
      for (j=0; j<ny; j++)
        {
          for (i=0; i<nx; i++)
            {
              // get the index of the grid point to the left
              jL = findAdjacentGridPoint(n1, x1, var[k][j][i]);

              // assign to proper bin
              if (var[k][j][i] < 0.5*(x1[jL+1] - x1[jL]))
                  pdf[jL] += plusPDF;
              else
                  pdf[jL+1] += plusPDF;
            } // end loop over x
        }  // end loop over y
    } // end loop over z

  return 1;
}
/* -------------------------------------------------------------------------- */

int binPDF_2D(int nx, int ny, int nz, double ***var1, double ***var2,
              int *n, double **x, double **pdf)
{
  int i, j, k, jL, kB;
  double plusPDF, norm;

  // assing local grid properties
  int     n1 = n[0];
  int     n2 = n[1];
  double *x1 = x[0];
  double *x2 = x[1];

  plusPDF = 1.0;
  norm    = (double)(nx*ny*nz);

  // initialise to zero
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      pdf[k][j] = 0.0;

  // loop over all data points
  for (k=0; k<nz; k++)
    {
      for (j=0; j<ny; j++)
        {
          for (i=0; i<nx; i++)
            {
              // get the index of the grid point
              jL = findAdjacentGridPoint(n1, x1, var1[k][j][i]); // left
              kB = findAdjacentGridPoint(n2, x2, var2[k][j][i]); // below

              // assign to proper bin
              if (var1[k][j][i] < 0.5*(x1[jL+1] - x1[jL]))
                {
                  if (var2[k][j][i] < 0.5*(x2[kB+1] - x2[kB]))
                    pdf[kB][jL] += plusPDF;
                  else
                    pdf[kB+1][jL] += plusPDF;                    
                }
              else
                {
                  if (var2[k][j][i] < 0.5*(x2[kB+1] - x2[kB]))
                    pdf[kB][jL+1] += plusPDF;
                  else
                    pdf[kB+1][jL+1] += plusPDF;                    
                }

            } // end loop over x
        }  // end loop over y
    } // end loop over z

  // normalise
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      pdf[k][j] /= norm;

  return -1;
}
/* -------------------------------------------------------------------------- */

void unscalePDF(int str, int transform, int n1, int n2, double *x1, double *x2,
                double **pdf)
{
  int j, k, row, col;

  // grid info
  double dx1, dx2, dA; 

  for (k=0; k<n2; k++)
    {
      if (n2 == 1)
        dx2 = 1.0;
      else
        {
          if (k == 0)
            dx2 = 0.5*(x2[1]-x2[0]);
          else if (k == n2-1)
            dx2 = 0.5*(x2[k]-x2[k-1]);
          else
            dx2 = 0.5*(x2[k+1]-x2[k-1]);
        }

      for (j=0; j<n1; j++)
        {
          if (j == 0)
            dx1 = 0.5*(x1[1]-x1[0]);
          else if (j == n1-1)
            dx1 = 0.5*(x1[j]-x1[j-1]);
          else
            dx1 = 0.5*(x1[j+1]-x1[j-1]);

          // account for triangles if untransformed
          if (transform == 0)
            {
              if (j == n1-1-k)
                dA = 0.5*dx1*dx2;
              else if (j > n1-1-k)
                dA = 1.0;
              else
                dA = dx1*dx2;
            }
          else
            dA = dx1*dx2;

          // structure vs. unstructured array
          if (str)
            pdf[k][j] /= dA;
          else
            pdf[0][k*n1 + j] /= dA;

        } /* end loop over columns */
    } /* end loop over rows */
}
/* ------------------------------------------------------------------------- */

double gamma(double xx)
{
  if (xx > 0.0)
    return exp(lnGamma(xx));
  else if (xx < 0.0 && fmod(xx,1.0) != 0.0)
    return -3.1415926535897931/(xx*sin(3.1415926535897931*xx)*exp(lnGamma(-xx)));
  
  fprintf(stderr, "Gamma function not defined for x = %8.6f", xx);
  exit(2);
}
/* ------------------------------------------------------------------------- */

double gammaProd(int n, int d, double *num, double *den)
{
  int i, npDen = 0, npNum = 0;
  double PI = 3.1415926535897931;
  double fgam = 0.0, cof = 1.0, sumPoles;

  int *maskNum = (int*) calloc(n, sizeof(int));
  int *maskDen = (int*) calloc(d, sizeof(int));

  // find the number of poles in the numerator and denominator
  for (i=0; i<n; i++)
    {
      if (num[i] <= 0.0 && fmod(num[i],1.0) == 0.0)
        npNum++;
      else
        maskNum[i] = 1;
    }
  for (i=0; i<d; i++)
    {
      if (den[i] <= 0.0 && fmod(den[i],1.0) == 0.0)
        npDen++;
      else
        maskDen[i] = 1;
    }

  if (npNum < npDen) // 
    {
      free(maskNum);
      free(maskDen);
      return 0.0;
    }
  else if (npNum > npDen)
    {
      fprintf(stderr, "Gamma product not defined:\n\t");
      for (i=0; i<n; i++) fprintf(stderr, "GAMMA(%13.6e)", num[i]);
      fprintf(stderr, "\n\t");
      for (i=0; i<n*20; i++) fprintf(stderr, "-");
      fprintf(stderr, "\n\t");
      for (i=0; i<d; i++) fprintf(stderr, "GAMMA(%13.6e)", den[i]);
      fprintf(stderr, "\n");
      exit(2);
    } 
  /* Otherwise, either no poles or equal number in numerator and denominator,
     thus they cancel with:
          lim GAMMA(a)/GAMMA(b) = (-1)^(a+b)*GAMMA(1-b)/GAMMA(1-a)
     Note: this doesn't really help for non-complex numbers, as it will only
           be defined if the sum of all the poles is even, giving
           -1^(sum_poles) = 1
  */
  if (npNum > 0)
    {
      sumPoles = 0.0;
      for (i=0; i<npNum; i++) if (!maskNum[i]) sumPoles += num[i];
      for (i=0; i<npDen; i++) if (!maskDen[i]) sumPoles += den[i];
      if (fmod(sumPoles,2.0) != 0.0)
        {
          fprintf(stderr, "Poles for gamma product not even: complex number\n");
          exit(2);
        }
    }

  // add the contribution to the numerator
  for (i=0; i<n; i++)
    {
      if (maskNum[i]) // not a pole
        {
          if (num[i] > 0.0)
            fgam += lnGamma(num[i]);
          else
            {
              fgam -= lnGamma(-num[i]);
              cof  *= -PI/(num[i]*sin(num[i]*PI));
            }
        }
      else // pole
        fgam -= lnGamma(1.0 - num[i]);
    }

  // add the contribution to the denominator
  for (i=0; i<d; i++)
    {
      if (maskDen[i]) // not a pole
        {
          if (den[i] > 0.0)
            fgam -= lnGamma(den[i]);
          else if (den[i] < 0.0 && fmod(den[i],1.0) != 0.0)
            {
              fgam += lnGamma(-den[i]);
              cof  *= -den[i]*sin(den[i]*PI)/PI;
            }
        }
      else // pole
        fgam += lnGamma(1.0 - den[i]);
    }
  return cof*exp(fgam);
}
/* -------------------------------------------------------------------------- */

double lnGamma(double xx)
{
  int j;
  double x, tmp, y, ser;
  double cof[14] = {57.1562356658629235,     -59.5979603554754912,
                    14.1360979747417471,      -0.491913816097620199,
                     0.339946499848118887e-4,  0.465236289270485756e-4,
                    -0.983744753048795646e-4,  0.158088703224912494e-3,
                    -0.210264441724104883e-3,  0.217439618115212643e-3,
                    -0.164318106536763890e-3,  0.844182239838527433e-4,
                    -0.261908384015814087e-4,  0.368991826595316234e-5};
  y = x = xx;
  tmp = x + 5.24218750000000000;
  tmp = (x + 0.5)*log(tmp) - tmp;
  ser = 0.999999999999997092;
  for (j=0; j<14; j++) ser += cof[j]/++y;
  return tmp + log(2.5066282746310005*ser/x);
}
/* ------------------------------------------------------------------------- */

double normalDistribution(double x, double mean, double var)
{	
  return 0.398942280401/sqrt(fabs(var))*exp(-0.5*pow(x-mean,2)/var);
}
/* ------------------------------------------------------------------------- */
