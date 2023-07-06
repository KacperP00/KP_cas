#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "betaPDF.h"
#include "solveUtils.h"
#include "hypergeometric.h"
#ifdef MKL
#include "mkl_vml.h"
#endif

int calcBetaPDF(int n, double *x, double mean, double var, double *pdf,
                int *pdfBound)
{
  int    i, j, j0, j1, jMean, pdfFlag;
  double alpha, beta, tmp, fgam;    // beta parameters
  double sumPDF;
  double sumPDFI;		    // inverse for improved performance
  double truncate = 1.0e-6;         // stop computing pdf below this limit
  double dx;                        // grid spacing
  double dx0 = 0.5*(x[1]-x[0]);     // for boundary approximation
  double dxN = 0.5*(x[n-1]-x[n-2]);
  double dxI; 	 		    // inverse grid spacing for improved performance
     
  // zero out the pdf
  for (j=0; j<n; j++) pdf[j] = 0.0;

  // ensure that the mean and variance are properly bounded
  clipMoments1D(&mean, &var);

  // check for special cases of pdf, e.g. delta function
  pdfFlag = checkPDF_1D(n, x, mean, var, pdf, pdfBound);

  // if pdfFlag is greater than zero, special case was found, so exit
  if (pdfFlag > 0)
    return pdfFlag;      
        
  // compute the shape parameters for the beta distribution
  tmp   = mean*(1-mean)/var - 1;
  alpha = mean*tmp;
  beta  = (1-mean)*tmp;
  if (fabs(alpha - 1) < 1e-9) alpha = 1.0;
  if (fabs(beta - 1) < 1e-9) beta = 1.0;
       
  /* ----- compute a beta distribution based on alpha and beta ----- */
  sumPDF  = 0.0;
  if (alpha == 1 && beta == 1) // uniform distribution (not likely)
    {
      j0 = 0;
      j1 = n-1;
      dxI = 1.0 / (x[j1]-x[j0]);
      pdf[j0] = 0.5*(x[1]-x[0])*dxI;
      pdf[j1] = 0.5*(x[j1]-x[j1-1])*dxI;
      sumPDF = pdf[j0] + pdf[j1];
      for (j=j0+1; j < j1; j++)
        {
          dx    = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = dx*dxI;
          sumPDF += pdf[j];
        }
    }
  // u shaped distribution
  else if ( (alpha < 1.0 && beta < 1.0) )
    {
      // need to calculate over whole domain, ignore the ends
      j0 = 0;
      j1 = n-1;

      // compute the PDF over all the interior points
      for (j=1; j < j1; j++)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = betaDist(x[j], alpha, beta)*dx;
          sumPDF  += pdf[j];
        }
    }
  // decreasing distribution
  else if ( (alpha < 1.0 && beta >= 1.0) || (alpha == 1.0 && beta > 1.0) )
    {
      // start at the left boundary and march till less than truncate
      j0 = 0;
      for (j = j0+1; j<n-1; j++)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = betaDist(x[j], alpha, beta)*dx;
                                        
          if (pdf[j] < truncate)
            {
              pdf[j] = 0.0;
              break;
            }
                                        
          j1 = j;
          sumPDF  += pdf[j];
        }
    }
  // increasing distribution
  else if ( (alpha > 1.0 && beta <= 1.0) || (alpha == 1.0 && beta < 1.0) ) 
    {
      // start at the right boundary and march till less than truncate
      j1 = n-1;
      for (j = j1-1; j>0; j--)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = betaDist(x[j], alpha, beta)*dx;
                                        
          if (pdf[j] < truncate)
            {
              pdf[j] = 0.0;
              break;
            }
                                        
          j0 = j;
          sumPDF += pdf[j];
        }
    }
  // unimodal distribution
  else 
    {
      j0 = 0;
      j1 = n-1;

      // find the location of the mesh point before the mean
      jMean = -1;
      do { jMean++; } while (x[jMean] < mean);
      jMean--;

      // start from the mean and decrease x till the pdf is negligible
      for (j=jMean; j > 0; j--)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = betaDist(x[j], alpha, beta)*dx;
                                        
          if (pdf[j] < truncate)
            {
              pdf[j] = 0.0;
              break;
            }
                                        
          j0 = j;
          sumPDF += pdf[j];
        }
                        
      // start from the mean and increase x till the pdf is negligible
      for (j=jMean+1; j < n - 1; j++)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = betaDist(x[j], alpha, beta)*dx;
                                        
          if (pdf[j] < truncate)
            {
              pdf[j] = 0.0;
              break;
            }
                                        
          j1 = j;
          sumPDF += pdf[j];
        }
    } // end different cases
     
  /* compute the boundary points explicitly if necessary using approximate
     expansion */
  if (!(alpha == 1.0 && beta == 1.0))
    {
      fgam = lnGamma(alpha+beta) - lnGamma(alpha) - lnGamma(beta);
      if (j0 == 0 || j0 == 1)
        {
          j0 = 0;
          pdf[j0] = exp(alpha*log(dx0) + fgam)/alpha;
          sumPDF += pdf[j0];
        }

      if (j1 == n-1 || j1 == n-2)
        {
          j1 = n-1;
          pdf[j1] = exp(beta*log(dxN) + fgam)/beta;
          sumPDF += pdf[j1];
        }
    }

  sumPDFI = 1.0 / sumPDF;
  // normalize the PDF
  for (j=j0; j<=j1; j++)
    pdf[j] = pdf[j]*sumPDFI;

  // store the limits of the non-zero pdf
  pdfBound[0] = j0;
  pdfBound[1] = j1;

  return 0;
}
/* ------------------------------------------------------------------------- */
int calcBetaPDFall(int n, double *x,double mean, double var, double *pdf,
                int *pdfBound)
{
  int    i, j, j0, j1, jMean, pdfFlag;
  double alpha, beta, tmp, fgam;    // beta parameters
  double sumPDF;
  double sumPDFI;		    // inverse for improved performance
  double truncate = 1.0e-6;         // stop computing pdf below this limit
  double dx;                        // grid spacing
  double dx0 = 0.5*(x[1]-x[0]);     // for boundary approximation
  double dxN = 0.5*(x[n-1]-x[n-2]);
  double dxI; 	 		    // inverse grid spacing for improved performance
#ifdef MKL     
  vmlSetMode( VML_EP | VML_FTZDAZ_ON | VML_ERRMODE_IGNORE );
#endif /*MKL*/
  double *x1work = calloc(n, sizeof(double*));
  double *x2work = calloc(n, sizeof(double*));

  // zero out the pdf
  for (j=0; j<n; j++) pdf[j] = 0.0;

  // ensure that the mean and variance are properly bounded
  clipMoments1D(&mean, &var);

  // check for special cases of pdf, e.g. delta function
  pdfFlag = checkPDF_1D(n, x, mean, var, pdf, pdfBound);

  // if pdfFlag is greater than zero, special case was found, so exit
  if (pdfFlag > 0) {
    free(x1work); free(x2work);
    return pdfFlag;      
  }
  // compute the shape parameters for the beta distribution
  tmp   = mean*(1-mean)/var - 1;
  alpha = mean*tmp;
  beta  = (1-mean)*tmp;
  if (fabs(alpha - 1) < 1e-9) alpha = 1.0;
  if (fabs(beta - 1) < 1e-9) beta = 1.0;
       
  if (alpha < 1.0e-16) alpha = 1.0e-16;
  if (beta  < 1.0e-16) beta  = 1.0e-16;
  /* ----- compute a beta distribution based on alpha and beta ----- */
  sumPDF  = 0.0;
  if (alpha == 1 && beta == 1) // uniform distribution (not likely)
    {
      j0 = 0;
      j1 = n-1;
      dxI = 1.0 / (x[j1]-x[j0]);
      pdf[j0] = 0.5*(x[1]-x[0])*dxI;
      pdf[j1] = 0.5*(x[j1]-x[j1-1])*dxI;
      sumPDF = pdf[j0] + pdf[j1];
      for (j=j0+1; j < j1; j++)
        {
          dx    = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = dx*dxI;
          sumPDF += pdf[j];
        }
    }
  // u shaped distribution
  else if ( (alpha < 1.0 && beta < 1.0) )
    {
      // need to calculate over whole domain, ignore the ends
      j0 = 0;
      j1 = n-1;

      // compute the PDF over all the interior points
      betaDistAll(&x[1], &x1work[1], &x2work[1], alpha, beta, j1-1, &pdf[1]);
      for (j=1; j < j1; j++)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = pdf[j]*dx;
          sumPDF  += pdf[j];
        }
    }
  // decreasing distribution
  else if ( (alpha < 1.0 && beta >= 1.0) || (alpha == 1.0 && beta > 1.0) )
    {
      // start at the left boundary and march till less than truncate
      j0 = 0;

      betaDistAll(&x[j0+1], &x1work[j0+1], &x2work[j0+1], alpha, beta, n-1-(j0+1), &pdf[j0+1]);

      for (j = j0+1; j<n-1; j++)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = pdf[j]*dx;
                                        
          if (pdf[j] < truncate)
            {
              pdf[j] = 0.0;
              break;
            }
                                        
          j1 = j;
          sumPDF  += pdf[j];
        }
    }
  // increasing distribution
  else if ( (alpha > 1.0 && beta <= 1.0) || (alpha == 1.0 && beta < 1.0) ) 
    {
      // start at the right boundary and march till less than truncate
      j1 = n-1;
      betaDistAll(&x[1], &x1work[1] ,&x2work[1], alpha, beta, j1-1, &pdf[1]);
      for (j = j1-1; j>0; j--)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = pdf[j]*dx;
                                        
          if (pdf[j] < truncate)
            {
              pdf[j] = 0.0;
              break;
            }
                                        
          j0 = j;
          sumPDF += pdf[j];
        }
    }
  // unimodal distribution
  else 
    {
      j0 = 0;
      j1 = n-1;

      // find the location of the mesh point before the mean
      jMean = -1;
      do { jMean++; } while (x[jMean] < mean);
      jMean--;

      // start from the mean and decrease x till the pdf is negligible
      betaDistAll(&x[1], &x1work[1] ,&x2work[1], alpha, beta, jMean, &pdf[1]);
      for (j=jMean; j > 0; j--)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = pdf[j]*dx;
                                        
          if (pdf[j] < truncate)
            {
              pdf[j] = 0.0;
              break;
            }
                                        
          j0 = j;
          sumPDF += pdf[j];
        }
                        
      // start from the mean and increase x till the pdf is negligible
      betaDistAll(&x[jMean+1],&x1work[jMean+1] ,&x2work[jMean+1], alpha, beta, n-1 - (jMean+1), &pdf[jMean+1]);

      for (j=jMean+1; j < n - 1; j++)
        {
          dx = 0.5*(x[j+1]-x[j-1]);
          pdf[j] = pdf[j]*dx;
                                        
          if (pdf[j] < truncate)
            {
              pdf[j] = 0.0;
              break;
            }
                                        
          j1 = j;
          sumPDF += pdf[j];
        }
    } // end different cases
     
  /* compute the boundary points explicitly if necessary using approximate
     expansion */
  if (!(alpha == 1.0 && beta == 1.0))
    {
      fgam = lnGamma(alpha+beta) - lnGamma(alpha) - lnGamma(beta);
      if (j0 == 0 || j0 == 1)
        {
          j0 = 0;
          pdf[j0] = exp(alpha*log(dx0) + fgam)/alpha;
          sumPDF += pdf[j0];
        }

      if (j1 == n-1 || j1 == n-2)
        {
          j1 = n-1;
          pdf[j1] = exp(beta*log(dxN) + fgam)/beta;
          sumPDF += pdf[j1];
        }
    }

  sumPDFI = 1.0 / sumPDF;
  // normalize the PDF
  for (j=j0; j<=j1; j++)
    pdf[j] = pdf[j]*sumPDFI;

  //if (sumPDF>DBL_MIN) 
  //   {
  //	sumPDFI = 1.0 / sumPDF;
  //      // normalize the PDF
  //      for (j=j0; j<=j1; j++)
  //        pdf[j] = pdf[j]*sumPDFI;
  //   }
  //else
  //   {
  //      fprintf(stderr,"Error: Cannot scale pdf (sum is eq 0.0)\n");
  //      sumPDFI  = 1.0;
  //      // normalize the PDF
  //      for (j=j0; j<=j1; j++)
  //        pdf[j] = pdf[j]/sumPDF;
  //   }


  // store the limits of the non-zero pdf
  pdfBound[0] = j0;
  pdfBound[1] = j1;
  free(x1work); free(x2work);

  return 0;
}
/* ------------------------------------------------------------------------- */

int calcDirichletPDF(int *n,double **x, double mean1, double var1, double mean2,
                     double var2, double var, double **pdf, int **pdfBound,
                     int transform)
{
  int    i, j, k, j0, j1, k0, k1, pdfFlag;
  double dx1, dx2, dZ1, dZ2, sigma, covar;
  double beta1, beta2, beta3, tmp, fgam, B1, B2;
  double Z1, Z2;
  double mean, Q, S;
  double sumPDF, meanPDF;
  double truncate = 1.0e-6;
  double eps = 1.0e-4;

  double yc, xc;

  // set local grid parameters
  int    n1  = n[0];
  int    n2  = n[1];
  double *x1 = x[0];
  double *x2 = x[1];

  double dx1_0  = 0.5*(x1[1]-x1[0]);
  double dx1_N  = 0.5*(x1[n1-1]-x1[n1-2]);
  double dx2_0  = 0.5*(x2[1]-x2[0]);
  double dx2_N  = 0.5*(x2[n2-1]-x2[n2-2]);
  double eps0, epsN, b1_0, b1_N, b2_0, b2_N;

  // zero out the pdf
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      pdf[k][j] = 0.0;
           
  // ensure that the mean and variance are properly bounded
  clipMoments2D(&mean1, &var1, &mean2, &var2, &var, &covar);
  mean = mean1 + mean2; // total mean

  // treat special cases for the pdf
  pdfFlag = checkPDF_2D(n, x, mean1, var1, mean2, var2, var, pdf, pdfBound,
                        transform, 0);

  // exit if special case was found
  if (pdfFlag > 0)
    return pdfFlag;

  // compute mean and variance parameters
  S = mean1*mean1 + mean2*mean2 + (1-mean)*(1-mean);
  Q = var1 + var2 + var;

  // compute the shape parameters for the beta distribution
  tmp   = (1.0 - S)/Q - 1.0;
  beta1 = mean1*tmp;
  beta2 = mean2*tmp;
  beta3 = (1.0-mean)*tmp;
  if (fabs(beta1 - 1) < 1e-9) beta1 = 1.0;
  if (fabs(beta2 - 1) < 1e-9) beta2 = 1.0;
  if (fabs(beta3 - 1) < 1e-9) beta3 = 1.0;

  /* ----- compute the pdf in all interior points ----- */
  sumPDF = 0.0;
  for (k=1; k < n2-1; k++)
    {
      dx2 = 0.5*(x2[k+1]-x2[k-1]);
      for (j=1; j < n1-1; j++)
        {
          dx1 = 0.5*(x1[j+1]-x1[j-1]);
          if (transform == 0)
            {
              if (x1[j] + x2[k] < 1.0)
                pdf[k][j] = dirichletDist(x1[j], x2[k], beta1, beta2, beta3);
              else
                pdf[k][j] = 0.0;
            }
          else if (transform == 1)
            pdf[k][j] = dirichletDistDoran(x1[j], x2[k], beta1, beta2, beta3);
          else if (transform == 2)
            pdf[k][j] = dirichletDistHasse(x1[j], x2[k], beta1, beta2, beta3);

          // scale the pdf
          // NOTE: effect of grid on scaling taken into account in distribution
          pdf[k][j] = pdf[k][j]*dx1*dx2;

          sumPDF += pdf[k][j];
        } // end loop over columns
    } // end loop over rows

  /* ----- compute pdf on the boundaries ----- */
  fgam = lnGamma(beta1+beta2+beta3)-lnGamma(beta1)-lnGamma(beta2)-lnGamma(beta3);

  if (transform == 0 || transform == 1)
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

  // compute the corners
  pdf[0][0]       = exp(b1_0*log(dx1_0) + b2_0*log(dx2_0) + fgam)/(b1_0*b2_0);
  pdf[0][n1-1]    = exp(b1_N*log(dx1_N) + b2_0*log(dx2_0) + fgam)/(b1_N*b2_0);
  pdf[n2-1][0]    = exp(b1_0*log(dx1_0) + b2_N*log(dx2_N) + fgam)/(b1_0*b2_N);
  pdf[n2-1][n1-1] = exp(b1_N*log(dx1_N) + b2_N*log(dx2_N) + fgam)/(b1_N*b2_N);

  if (transform == 0)
    {
      // adjust corners according to transformation
      pdf[0][n1-1] = exp(b1_0*log(dx1_0) + b2_0*log(dx2_0/(1.0-0.5*dx1_N))
                         + fgam)/(b1_0*b2_0);
      pdf[n2-1][0]    = exp(b1_0*log(dx1_0) + b2_N*log(dx2_N/(1.0-0.5*dx1_0))
                            + fgam)/(b1_0*b2_N);
    }

  // add corners to sum of pdf
  sumPDF = sumPDF + pdf[0][0] + pdf[0][n1-1] + pdf[n2-1][0] + pdf[n2-1][n1-1];

  // left and right sides
  eps0 = b1_0*log(dx1_0);
  epsN = b1_N*log(dx1_N);
  for (k=1; k<n2-1; k++)
    {
      dx2 = 0.5*(x2[k+1]-x2[k-1]);
      B2 = (b2_0-1.0)*log(x2[k]) + (b2_N-1.0)*log(1.0-x2[k]);
      pdf[k][0]    = exp(eps0 + B2 + fgam)*dx2/b1_0;
      pdf[k][n1-1] = exp(epsN + B2 + fgam)*dx2/b1_N;

      sumPDF = sumPDF + pdf[k][0] + pdf[k][n1-1];
    }

  // top and bottom
  eps0 = b2_0*log(dx2_0);
  epsN = b2_N*log(dx2_N);
  for (j=1; j<n1-1; j++)
    {
      if (transform == 0)
        {
          eps0 = b2_0*log(dx2_0/(1.0-x1[j]));
          epsN = b2_N*log(dx2_N/(1.0-x1[j]));
        }

      dx1 = 0.5*(x1[j+1]-x1[j-1]);
      B1 = (b1_0-1.0)*log(x1[j])+(b1_N-1.0)*log(1.0-x1[j]);
      pdf[0][j]    = exp(eps0 + B1 + fgam)*dx1/b2_0;
      pdf[n2-1][j] = exp(epsN + B1 + fgam)*dx1/b2_N;

      sumPDF = sumPDF + pdf[0][j] + pdf[n2-1][j];
    }

   if (transform == 0)
      {
        /* remove the transform Z1+Z2 = 1
           NOTE: this only works for equal grid spacing in each direction */
        for (k=1; k<n2-1; k++)
          {
            //pdf[k][n1-1-k] = 0.5*pdf[n2-1][n1-1-k]/(1.0-x1[n1-1-k]);
            pdf[k][n1-1-k] = pdf[n2-1][n1-1-k];
            pdf[k][n1-1]   = pdf[n2-1][n1-1-k] = 0.0;
          }
        pdf[n2-1][n1-1] = 0.0;
      }

  /* normalise the PDF */
  for (k=0; k<n2; k++)
    {
      pdfBound[k][0] = 0;
      pdfBound[k][1] = n1-1;
      for(j=0; j<n1; j++)
        pdf[k][j] = pdf[k][j]/sumPDF;
    }

  return 0;
}
/* ------------------------------------------------------------------------- */

int calcJointBetaPDF(int *n,double **x, double mean1, double var1, double mean2,
                     double var2, double var, double **pdf, int **pdfBound,
                     int transform)
{
  int    i, j, k, j0, j1, k0, k1, pdfFlag, check = 0, flag;
  double dx1, dx2, dA, dZ1, dZ2;
  double beta[6], tmp, moms[6];
  double Z1, Z2;
  double mean, covar, P, Q, S;
  double sumPDF, meanPDF[3];
  double truncate = 1.0e-6;
  double eps = 1.0e-4;
  double x1c, x2c;
  double delta = 0.5;

  /* storage struct for data required by function and Jacobian evaluation */
  statStruct args;

  int    n1, n2;
  double *x1, *x2;
  double yc, xc;

  // assign parameters for computation
  args.nMom = 6;
  args.n1 = n1 = n[0];
  args.n2 = n2 = n[1];
  args.x1 = x1 = x[0];
  args.x2 = x2 = x[1];
  args.transform = transform;

  // zero out the pdf
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      pdf[k][j] = 0.0;
           
  // ensure that the mean and variance are properly bounded
  clipMoments2D(&mean1, &var1, &mean2, &var2, &var, &covar);
  mean = mean1 + mean2; // total mean

  // treat special cases for the pdf
  pdfFlag = checkPDF_2D(n, x, mean1, var1, mean2, var2, var, pdf, pdfBound,
                        transform, 0);

  // exit if special case was found
  if (pdfFlag > 0)
    return pdfFlag;
 
  // set the mean and variance vectors for passing
  moms[0] = 1.0;
  moms[1] = args.mean[0] = mean1;
  moms[2] = args.mean[1] = mean2;
  moms[3] = args.var[0]  = var1;
  moms[4] = args.var[1]  = var2;
  moms[5] = args.var[2]  = covar;

  // compute mean and variance parameters
  S = mean1*mean1 + mean2*mean2 + (1-mean)*(1-mean);
  Q = var1 + var2 + var;

  // compute the shape parameters for dirichlet distribution as initial guess
  tmp     = (1.0 - S)/Q - 1.0;
  beta[0] = 1.0;
  beta[1] = mean1*tmp;
  beta[2] = mean2*tmp;
  beta[3] = beta[4] = 1.0;
  beta[5] = (1.0-mean)*tmp;

  // use uniform distribution as initial guess
  beta[0] = beta[1] = beta[2] = beta[3] = beta[4] = beta[5] = 1.0;

#ifdef MDUC_PDF_KINSOL
  sdialStruct kin;  // allocate struct for holding KINSOL memory
  /* set up the solver */
  initKinJbeta(&kin, &args, beta);

  /* set the residual and Jacobian functions */
  args.fcn = &jointBetaRes;
  args.jac = &jointBetaJac;

  /* call the KINSOL solver */
  flag = KINSol(kin.mem, kin.Y, KIN_LINESEARCH, kin.Yscale, kin.Fscale);
  if (flag < 0)
    fprintf(stderr, "KINSOL Error: %i\n", flag);

  /* set the solution into the parameter array */
  for (i=0; i<6; i++)
    beta[i] = kin.Ydata[i];

  destroyKinJbeta(&kin);
#elif MDUC_PDF_BROYDEN
  /* use Broyden's method to get parameters */
  broyden(beta, 6, &check, &args, jointBetaRes, jointBetaJac);
#elif MDUC_PDF_LM
  int maxIts;
  double lmTol[4];
  double lmInfo[LM_INFO_SZ];

  maxIts   = 200;
  lmTol[0] = 1.0e-9;
  lmTol[1] = 1.0e-9;
  lmTol[2] = 1.0e-9;
  lmTol[3] = 1.0e-15;

  /* call the external Levenberg-Marquardt library */
  dlevmar_der(jointBetaFuncLM, jointBetaJacLM, beta, moms, 6, 6, maxIts, NULL,
              lmInfo, NULL, NULL, &args);

  if (lmInfo[6] == 3) check = 1;
#else

#ifdef MDUC_JBETA_ADAPT
  args.a  = beta;
  args.nxG    = 1;
  args.xG     = calloc(args.nxG, sizeof(double));
  args.nyG    = calloc(args.nxG, sizeof(int));
  args.nyG[0] = 1;
  args.yG     = calloc(args.nyG[0], sizeof(double*));

  /* use Newton solver to get parameters */
  newt(beta, 6, &check, &args, jointBetaFuncODE, jointBetaJacODE);
#else
  /* use Newton solver to get parameters */
  newt(beta, 6, &check, &args, jointBetaRes, jointBetaJac);
#endif /* MDUC_JBETA_ADAPT */
#endif /* MDUC_PDF_KINSOL */

  /* ----- compute the pdf based on the betas computed ----- */
  sumPDF = meanPDF[1] = meanPDF[2] = 0.0;
  for (k=0; k<n2; k++)
    {
      // get the current mesh point and spacing for ordinate
      if (k == 0)
        {
          dx2 = 0.5*(x2[1]-x2[0]);
          x2c = delta*dx2;
        }
      else if (k == n2-1)
        {
          dx2 = 0.5*(x2[k] - x2[k-1]);
          x2c = x2[k] - delta*dx2;
        }
      else
        {
          dx2 = 0.5*(x2[k+1] - x2[k-1]);
          x2c = x2[k];
        }

      // start loop over columns
      for (j=0; j<n1; j++)
        {
          // get the current mesh point and spacing for abscissa
          if (j == 0)
            {
              dx1 = 0.5*(x1[1]-x1[0]);
              x1c = delta*dx1;
            }
          else if (j == n1-1)
            {
              dx1 = 0.5*(x1[j] - x1[j-1]);
              x1c = x1[j] - delta*dx1;
            }
          else
            {
              dx1 = 0.5*(x1[j+1] - x1[j-1]);
              x1c = x1[j];
            }

          // get mixture fractions according to transform
          if (transform == 0)
            {
              if (j == n1-1-k)
                {
                  Z1 = x1c - dx1/6.0;
                  Z2 = x2c - dx2/6.0;
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
                  Z1 = x1c;
                  Z2 = x2c;
                  dA = dx1*dx2;
                }
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
 
          // get the value of the distribution
          P = jointBetaDist(Z1, Z2, beta);

          pdf[k][j] = P*dA;

          /* keep track of the computed means */
          meanPDF[1] += Z1*pdf[k][j];
          meanPDF[2] += Z2*pdf[k][j];
          sumPDF     += pdf[k][j];
        } // end loop over columns
    } // end loop over rows

  /* normalise the PDF */
  for (k=0; k<n2; k++)
    {
      for(j=0; j<n1; j++)
        {
          pdfBound[k][0] = 0;
          pdfBound[k][1] = n1-1;
          pdf[k][j] = pdf[k][j]/sumPDF;
        }
    }
  for (i=0; i<6; i++) fprintf(stderr, "b%i := %8.6f: ", i, beta[i]);
  fprintf(stderr, "\n");

#ifdef MDUC_JBETA_ADAPT
  

  FILE *pFile = fopen("jbetaGrid.dat", "w");

  int totPts = 0;
  for (i=0; i<args.nxG; i++) totPts += args.nyG[i];
  fprintf(pFile, "# total points: %i\n", totPts);
  int offset = 0;
  for (i=0; i<args.nxG; i++)
    {
      offset += args.nyG[i];
      if (i == 0)
        dx1 = 0.0;
      else
        dx1 = args.xG[i]-args.xG[i-1];

      for (j=0; j<args.nyG[i]; j++)
        {
          if (j == 0)
            dx2 = 0.0;
          else
            dx2 = args.yG[offset+j] - args.yG[offset+j-1];

          Z1 = fmin(0.999, fmax(0.001,args.xG[i]));
          Z2 = fmin(0.999, fmax(0.001,args.yG[offset+j]*(1.0-Z1)));
          dA = (1.0-Z1)*dx1*dx2;

          P = jointBetaDist(Z1, Z2, beta)*dA;

          fprintf(pFile, "%8.6f %8.6f %15.6e %8.6f\n",args.xG[i],
                  args.yG[offset+j],P,0.0);
        }
      fprintf(pFile, "\n");
    }
  fclose(pFile);
#endif /* MDUC_JBETA_ADAPT */

  if (check)
    return -10;

  return -1;
}
/* ------------------------------------------------------------------------- */

void jointBetaRes(int nF, double *a, double *f, void *args)
{
  // dereference variables passed from solver
  statStruct *beta = (statStruct*) args;

  /* compute the moments from the parameters */
  jointBetaFunc(nF, a, f, args);

  /* compute the residual */
  f[0] -= 1.0;
  f[1] -= beta->mean[0];
  f[2] -= beta->mean[1];
  f[3] -= beta->var[0];
  f[4] -= beta->var[1];
  f[5] -= beta->var[2];
}
/* ------------------------------------------------------------------------- */

void jointBetaFunc(int nF, double *a, double *f, void *args)
{
  int    i, j, k;
  int    n1, n2;
  double *x1, *x2, Z1, Z2;
  double *mean, *var;
  double dx1, dx2, dA;
  double tmp;
  double P;            // value of distribution
  double x1c, x2c; // current grid location
  double delta = 0.5;
  
  // dereference variables passed from solver
  statStruct *beta = (statStruct*) args;

  // reallocate to local variables for readability
  n1 = beta->n1;
  n2 = beta->n2;
  x1 = beta->x1;
  x2 = beta->x2;
  mean = beta->mean;
  var  = beta->var;

  // set all values to zero
  for (i=0; i<nF; i++) f[i] = 0.0;

  // sum over the interior points
  for (k=0; k<n2; k++)
    {
      // get the current mesh point and spacing for ordinate
      if (k == 0)
        {
          dx2 = 0.5*(x2[1]-x2[0]);
          x2c = delta*dx2;
        }
      else if (k == n2-1)
        {
          dx2 = 0.5*(x2[k] - x2[k-1]);
          x2c = x2[k] - delta*dx2;
        }
      else
        {
          dx2 = 0.5*(x2[k+1] - x2[k-1]);
          x2c = x2[k];
        }

      // start loop over columns
      for (j=0; j<n1; j++)
        {
          // get the current mesh point and spacing for abscissa
          if (j == 0)
            {
              dx1 = 0.5*(x1[1]-x1[0]);
              x1c = delta*dx1;
            }
          else if (j == n1-1)
            {
              dx1 = 0.5*(x1[j] - x1[j-1]);
              x1c = x1[j] - delta*dx1;
            }
          else
            {
              dx1 = 0.5*(x1[j+1] - x1[j-1]);
              x1c = x1[j];
            }

          // get mixture fractions according to transform
          if (beta->transform == 0)
            {
              /* if ( (k == 0 && j == n1-1) || (k == n2-1 && j == 0) ) */
              if (j == n1-1-k)
                {
                  Z1 = x1c - dx1/6.0;
                  Z2 = x2c - dx2/6.0;
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
                  Z1 = x1c;
                  Z2 = x2c;
                  dA = dx1*dx2;
                }
            }
          else if (beta->transform == 1)
            {
              Z1 = x1c;
              Z2 = x2c*(1-x1c);
              dA = dx1*dx2*(1-x1c);
            }
          else if (beta->transform == 2)
            {
              Z1 = x1c*(1.0 - x2c);
              Z2 = x1c*x2c;
              dA = dx1*dx2*x1c;
            }

          P = jointBetaDist(Z1, Z2, a);

          // normalisation condition
          f[0] += P*dA;

          // compute the means
          f[1] += P*Z1*dA; // mean1
          f[2] += P*Z2*dA; // mean2

          // compute the variances
          f[3] += P*(Z1-mean[0])*(Z1-mean[0])*dA; // variance 1
          f[4] += P*(Z2-mean[1])*(Z2-mean[1])*dA; // variance 2
          f[5] += P*(Z1-mean[0])*(Z2-mean[1])*dA; // co-variance
        } // end loop over columns
    } // end loop over rows
} 
/* -------------------------------------------------------------------------- */

void jointBetaJac(int nF, double *a, double **r, void *args)
{
  int    i, j, k, p, q;
  int    n1, n2;
  double *x1, *x2, Z1, Z2;
  double *mean, *var;
  double dx1, dx2, dA;
  double tmp, expTerm, d[6], f[6];
  double P;
  double x1c, x2c;
  double delta = 0.5;
  
  // dereference variables passed from solver
  statStruct *beta = (statStruct*) args;

  // reallocate to local variables for readability
  n1 = beta->n1;
  n2 = beta->n2;
  x1 = beta->x1;
  x2 = beta->x2;
  mean = beta->mean;
  var  = beta->var;

  // set all values to zero
  for (k=0; k<nF; k++)
    for (i=0; i<nF; i++)
      r[k][i] = 0.0;

  // sum over the interior points
  for (k=0; k<n2; k++)
    {
      // get the current mesh point and spacing for ordinate
      if (k == 0)
        {
          dx2 = 0.5*(x2[1]-x2[0]);
          x2c = delta*dx2;
        }
      else if (k == n2-1)
        {
          dx2 = 0.5*(x2[k] - x2[k-1]);
          x2c = x2[k] - delta*dx2;
        }
      else
        {
          dx2 = 0.5*(x2[k+1] - x2[k-1]);
          x2c = x2[k];
        }

      // start loop over columns
      for (j=0; j<n1; j++)
        {
          // get the current mesh point and spacing for abscissa
          if (j == 0)
            {
              dx1 = 0.5*(x1[1]-x1[0]);
              x1c = delta*dx1;
            }
          else if (j == n1-1)
            {
              dx1 = 0.5*(x1[j] - x1[j-1]);
              x1c = x1[j] - delta*dx1;
            }
          else
            {
              dx1 = 0.5*(x1[j+1] - x1[j-1]);
              x1c = x1[j];
            }

          // get mixture fractions according to transform
          if (beta->transform == 0)
            {
              if (j == n1-1-k)
                {
                  Z1 = x1c - dx1/6.0;
                  Z2 = x2c - dx2/6.0;
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
                  Z1 = x1c;
                  Z2 = x2c;
                  dA = dx1*dx2;
                }
            }
          else if (beta->transform == 1)
            {
              Z1 = x1c;
              Z2 = x2c*(1-x1c);
              dA = dx1*dx2*(1-x1c);
            }
          else if (beta->transform == 2)
            {
              Z1 = x1c*(1.0 - x2c);
              Z2 = x1c*x2c;
              dA = dx1*dx2*x1c;
            }

          // set components to zero
          P = 0.0;
          for (i=0; i<nF; i++) d[i] = f[i] = 0.0;
          
          // get the value of the probability
          P = jointBetaDist(Z1, Z2, a);

          // compute the extra terms for the derivatives
          d[0] = P*dA;
          d[1] = P*log(Z1)*dA;
          d[2] = P*log(Z2)*dA;
          d[3] = P*log(1.0-Z1)*dA;
          d[4] = P*log(1.0-Z2)*dA;
          d[5] = P*log(1.0-Z1-Z2)*dA;

          // compute the function terms
          f[0] = 1.0;
          f[1] = Z1;
          f[2] = Z2;
          f[3] = (Z1-mean[0])*(Z1-mean[0]);
          f[4] = (Z2-mean[1])*(Z2-mean[1]);
          f[5] = (Z1-mean[0])*(Z2-mean[1]);

          // fill the Jacobian
          for (p=0; p<nF; p++)
            for (q=0; q<nF; q++)
              r[p][q] += f[p]*d[q];
        } // end loop over columns
    } // end loop over rows
} 
/* -------------------------------------------------------------------------- */

/* wrapper function for internal Levenberg-Marquardt minimisation function */
void jointBetaLevMar(double *a, double *y, double **dyda,
                     int na, void *args)
{
  /* compute the values of the moments */
  jointBetaFunc(na, a, y, args);

  /* compute the Jacobian and set to dyda */
  jointBetaJac(na, a, dyda, args);
}
/* -------------------------------------------------------------------------- */

/* wrapper functions for external Levenberg-Marquardt library */
void jointBetaFuncLM(double *p, double *hx, int m, int n, void *args)
{
  jointBetaFunc(m, p, hx, args);
}
/* -------------------------------------------------------------------------- */

void jointBetaJacLM(double *p, double *j, int m, int n, void *args)
{
  int i,k;
  double **jac;
  jac = (double**) calloc(n, sizeof(double*));
  for (k=0; k<n; k++)
    jac[k] = (double*) calloc(m, sizeof(double));

  jointBetaJac(m, p, jac, args);

  for (k=0; k<n; k++)
    for (i=0; i<m; i++)
      j[k*m+i] = jac[k][i];

  for (k=0; k<n; k++) free(jac[k]); free(jac);
}
/* -------------------------------------------------------------------------- */

#ifdef MDUC_PDF_KINSOL
/* functions for KINSOL Newton solver */
void initKinJbeta(sdialStruct *kin, statStruct *stats, double *a)
{
  int i, flag;

  kin->nEq = stats->nMom;

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
  flag = KINMalloc(kin->mem, kinJbetaRes, kin->tmp);
  if (flag < 0)
    fprintf(stderr, "KINMalloc error: %i\n", flag);

  // attach a dense solver
  flag = KINDense(kin->mem, kin->nEq);
  if (flag < 0)
    fprintf(stderr, "KINDense error: %i\n", flag);

  // set pointers to data and function calls
  flag = KINSetFdata(kin->mem, stats);
  if (flag < 0)
    fprintf(stderr, "KINSetFdata error: %i\n", flag);
  flag = KINDenseSetJacFn(kin->mem, kinJbetaJac, stats); 
  if (flag < 0)
    fprintf(stderr, "KINDenseSetJacFn error: %i\n", flag);
#else

  flag = KINInit(kin->mem, kinJbetaRes, kin->tmp);
  flag = KINDense(kin->mem, kin->nEq);
  flag = KINSetUserData(kin->mem, stats);
  flag = KINDlsSetDenseJacFn(kin->mem, kinJbetaJac);

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
  kin->FscaleData[1] = fabs(1.0/(stats->mean[0]+1e-9));
  kin->FscaleData[2] = fabs(1.0/(stats->mean[1]+1e-9));
  kin->FscaleData[3] = fabs(1.0/(stats->var[0]+1e-9));
  kin->FscaleData[4] = fabs(1.0/(stats->var[1]+1e-9));
  kin->FscaleData[5] = fabs(1.0/(stats->var[2]+1e-9));

  // copy the initial guess into the solution vector
  for (i=0; i<kin->nEq; i++)
    kin->Ydata[i] = a[i];
}
/* -------------------------------------------------------------------------- */

void destroyKinJbeta(sdialStruct *kin)
{
  N_VDestroy_Serial(kin->Y);
  N_VDestroy_Serial(kin->tmp);
  N_VDestroy_Serial(kin->Yscale);
  N_VDestroy_Serial(kin->Fscale);
  KINFree(&kin->mem);
}
/* -------------------------------------------------------------------------- */

int kinJbetaRes(N_Vector u, N_Vector fval, void *f_data)
{
  // dereference storage variables
  statStruct *args = (statStruct*) f_data;
  realtype *a, *f;
  a = NV_DATA_S(u);
  f = NV_DATA_S(fval);

  jointBetaRes(args->nMom, a, f, args);
  return 0;
}
/* -------------------------------------------------------------------------- */

#if SUNDIALS_2_3
int kinJbetaJac(long int N, DenseMat J, N_Vector u, N_Vector fu,
                void *jac_data, N_Vector tmp1, N_Vector tmp2)
#else
int kinJbetaJac(int N, N_Vector u, N_Vector fu, DlsMat J, void *jac_data,
                N_Vector tmp1, N_Vector tmp2)
#endif /* SUNDIAL_2_3 */
{
  int i, j;
  // dereference storage variables
  statStruct *args = (statStruct*) jac_data;
  realtype *a, **Jtmp;

  // allocate memory for temporary Jacobian
  Jtmp = (double**) calloc(N, sizeof(double*));
  for (i=0; i<N; i++) Jtmp[i] = (double*) calloc(N, sizeof(double));

  /* pointers to data */
  a = NV_DATA_S(u);

  /* call the Jacobian function */
  jointBetaJac(args->nMom, a, Jtmp, args);

  // copy the computed Jacobian to J
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      DENSE_ELEM(J,i,j) = Jtmp[i][j];

  for (i=0; i<N; i++) free(Jtmp[i]); free(Jtmp);
  return 0;
}
/* -------------------------------------------------------------------------- */
#endif /* MDUC_PDF_KINSOL */

int calcAppellBetaPDF(int *n, double **x, double mean1,double var1,double mean2,
                      double var2, double var, double **pdf, int **pdfBound,
                      int transform)
{
  int    i, j, k, j0, j1, k0, k1, pdfFlag, check, appellF;
  double dx1, dx2, dA, dZ1, dZ2, C2, C3;
  double beta[6], AB[5], tmp, nGam[4],dGam[2];
  double Z1, Z2;
  double mean, covar, P, Q, S, rho;
  double sumPDF, meanPDF;
  double truncate = 1.0e-6;
  double eps = 1.0e-4;
  double x1c, x2c;
  double norm, delta = 0.25;

  statStruct args;
  void (*appellBetaFunc)(int, double *, double *, void *);

  int    n1, n2;
  double *x1, *x2;
  double yc, xc;

  // assign parameters for computation
  args.nMom = 6;
  args.n1 = n1 = n[0];
  args.n2 = n2 = n[1];
  args.x1 = x1 = x[0];
  args.x2 = x2 = x[1];
  args.transform = transform;

  // zero out the pdf
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      pdf[k][j] = 0.0;
           
  // ensure that the mean and variance are properly bounded
  clipMoments2D(&mean1, &var1, &mean2, &var2, &var, &covar);
  mean = mean1 + mean2; // total mean

  // treat special cases for the pdf
  pdfFlag = checkPDF_2D(n, x, mean1, var1, mean2, var2, var, pdf, pdfBound,
                        transform, 0);

  // exit if special case was found
  if (pdfFlag > 0)
    return pdfFlag;
 
  // set the mean and variance vectors for passing
  args.mean[0] = mean1;
  args.mean[1] = mean2;
  args.var[0]  = var1;
  args.var[1]  = var2;
  args.var[2]  = covar;
  /* args.tol     = tol; */
  rho = covar/sqrt(var1*var2);

  /* if (covar <= 0.0) */
    appellF = 3;
/*   else */
/*     appellF = 2; */

  if (appellF == 3)
    appellBetaFunc = &appellF3BetaFunc;
  else if (appellF == 2)
    appellBetaFunc = &appellF2BetaFunc;

  // compute mean and variance parameters
  S = mean1*mean1 + mean2*mean2 + (1-mean)*(1-mean);
  Q = var1 + var2 + var;

  // compute the shape parameters for dirichlet distribution as initial guess
  tmp     = (1.0 - S)/Q - 1.0;
  beta[0] = 1.0;
  AB[0] = beta[1] = mean1*tmp;
  AB[1] = beta[2] = mean2*tmp;
  AB[2] = AB[3] = beta[3] = beta[4] = 1.0;
  AB[4] = beta[5] = (1.0-mean)*tmp;

  // use uniform distribution as initial guess
  AB[0] = AB[1] = AB[2] = AB[3] = AB[4] = 0.99;

  // call Newton solver to get parameters
  newt(AB, 5, &check, &args, appellBetaFunc, NULL);
  for (i=1; i<6; i++) beta[i] = AB[i-1];

  if (appellF == 2)
    {
      nGam[0] = beta[1];
      nGam[1] = beta[2];
      nGam[2] = beta[3];
      nGam[3] = beta[4];
      dGam[0] = beta[1]+beta[3];
      dGam[1] = beta[2]+beta[4];
      C2 = 1/(gammaProd(4, 2, nGam, dGam)*
              appell_F2(1-beta[5],beta[1],beta[2],dGam[0],dGam[1],1,1,
                        5000,1.0e-9));
      fprintf(stderr, "C2 = %15.6e (%15.6e)\n", C2, log(C2));
      beta[0] = log(C2);
    }
  else if (appellF == 3)
    {
      nGam[0] = beta[1];
      nGam[1] = beta[2];
      nGam[2] = beta[5];
      tmp = beta[1]+beta[2]+beta[5];

      // check to make sure that the normalisation condition is satisfied
      C3 = 1/(gammaProd(3, 1, nGam, &tmp)*
              appell_F3(1-beta[3],1-beta[4],beta[1],beta[2],tmp,1,1,
                        1000,1.0e-9));
      fprintf(stderr, "C3 = %15.6e (%15.6e)\n", C3, log(C3));
      beta[0] = log(C3);
    }
  
  for (i=0; i<6; i++) fprintf(stderr, "b%i := %8.6f: ", i, beta[i]);
  fprintf(stderr, "\n");

  /* ----- compute the pdf based on the betas computed ----- */
  sumPDF = 0.0;
  for (k=0; k<n2; k++)
    {
      // get the current mesh point and spacing for ordinate
      if (k == 0)
        {
          dx2 = 0.5*(x2[1]-x2[0]);
          x2c = delta*dx2;
        }
      else if (k == n2-1)
        {
          dx2 = 0.5*(x2[k] - x2[k-1]);
          x2c = x2[k] - delta*dx2;
        }
      else
        {
          dx2 = 0.5*(x2[k+1] - x2[k-1]);
          x2c = x2[k];
        }

      // start loop over columns
      for (j=0; j<n1; j++)
        {
          // get the current mesh point and spacing for abscissa
          if (j == 0)
            {
              dx1 = 0.5*(x1[1]-x1[0]);
              x1c = delta*dx1;
            }
          else if (j == n1-1)
            {
              dx1 = 0.5*(x1[j] - x1[j-1]);
              x1c = x1[j] - delta*dx1;
            }
          else
            {
              dx1 = 0.5*(x1[j+1] - x1[j-1]);
              x1c = x1[j];
            }

          // get mixture fractions according to transform
          if (transform == 0)
            {
              if (j == n1-1-k)
                {
                  Z1 = x1c - dx1/6.0;
                  Z2 = x2c - dx2/6.0;
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
                  Z1 = x1c;
                  Z2 = x2c;
                  dA = dx1*dx2;
                }
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
 
          // get the value of the distribution
          P = jointBetaDist(Z1, Z2, beta);

          pdf[k][j] = P*dA;

          sumPDF += pdf[k][j];
        } // end loop over columns
    } // end loop over rows

  /* normalise the PDF */
  for (k=0; k<n2; k++)
    for(j=0; j<n1; j++)
      pdf[k][j] = pdf[k][j]/sumPDF;

  return -1;
}
/* ------------------------------------------------------------------------- */

void appellF2BetaFunc(int nF, double *a, double *f, void *args)
{
  int    i, j, k;
  int    nmax = 2000;
  double *mean, *var;
  double eps = 1.0e-9;
  
  // dereference variables passed from solver
  statStruct *beta = (statStruct*) args;

  // reallocate to local variables for readability
  mean = beta->mean;
  var  = beta->var;

  // set all values to zero
  for (i=0; i<nF; i++) f[i] = 0.0;

  // using Appell F2
  double alpha, b, bp, g, gp;
  double C2, F2;
  alpha = 1 - a[4];
  b     = a[0];
  bp    = a[1];
  g     = a[0] + a[2];
  gp    = a[1] + a[3];

  // standard F2 in denominator
  C2 = appell_F2(alpha, b, bp, g, gp, 1, 1, nmax, eps);

  // mean 1
  F2 = appell_F2(alpha, b+1, bp, g+1, gp, 1, 1, nmax, eps);
  f[0] = b/g*F2/C2;

  // mean 2
  F2 = appell_F2(alpha, b, bp+1, g, gp+1, 1, 1, nmax, eps);
  f[1] = bp/gp*F2/C2;

  // variance 1
  F2 = appell_F2(alpha, b+2, bp, g+2, gp, 1, 1, nmax, eps);
  f[2] = b*(b+1)/(g*(g+1))*F2/C2 - 2*f[0]*mean[0] + mean[0]*mean[0];

  // variance 2
  F2 = appell_F2(alpha, b, bp+2, g, gp+2, 1, 1, nmax, eps);
  f[3] = bp*(bp+1)/(gp*(gp+1))*F2/C2 - 2*f[1]*mean[1] + mean[1]*mean[1];

  // covariance
  F2 = appell_F2(alpha, b+1, bp+1, g+1, gp+1, 1, 1, nmax, eps);
  f[4] = b*bp/(g*gp)*F2/C2 - f[0]*mean[1] - f[1]*mean[0] + mean[0]*mean[1];

  f[0] -= mean[0];
  f[1] -= mean[1];
  f[2] -= var[0];
  f[3] -= var[1];
  f[4] -= var[2];
} 
/* -------------------------------------------------------------------------- */

void appellF3BetaFunc(int nF, double *a, double *f, void *args)
{
  int    i, j, k;
  int    nmax = 1000;
  double *mean, *var;
  double eps = 1.0e-9;
  
  // dereference variables passed from solver
  statStruct *beta = (statStruct*) args;

  // reallocate to local variables for readability
  mean = beta->mean;
  var  = beta->var;

  // set all values to zero
  for (i=0; i<nF; i++) f[i] = 0.0;

  // using Appell F3
  double alpha, alpha_p, b, bp, g;
  double F3, C3;
  alpha   = 1 - a[2];
  alpha_p = 1 - a[3];
  b       = a[0];
  bp      = a[1];
  g       = a[0] + a[1] + a[4];

  // standard F3 in denominator
  C3 = appell_F3(alpha, alpha_p, b, bp, g, 1, 1, nmax, eps);

  // mean 1
  F3 = appell_F3(alpha, alpha_p, b+1, bp, g+1, 1, 1, nmax, eps);
  f[0] = b/g*F3/C3;

  // mean 2
  F3 = appell_F3(alpha, alpha_p, b, bp+1, g+1, 1, 1, nmax, eps);
  f[1] = bp/g*F3/C3;

  // variance 1
  F3 = appell_F3(alpha, alpha_p, b+2, bp, g+2, 1, 1, nmax, eps);
  //f[2] = b*(b+1)/(g*(g+1))*F3/C3 - f[0]*f[0];
  f[2] = b*(b+1)/(g*(g+1))*F3/C3 - 2*f[0]*mean[0] + mean[0]*mean[0];

  // variance 2
  F3 = appell_F3(alpha, alpha_p, b, bp+2, g+2, 1, 1, nmax, eps);
  //f[3] = bp*(bp+1)/(g*(g+1))*F3/C3 - f[1]*f[1];
  f[3] = bp*(bp+1)/(g*(g+1))*F3/C3 - 2*f[1]*mean[1] + mean[1]*mean[1];

  // covariance
  F3 = appell_F3(alpha, alpha_p, b+1, bp+1, g+2, 1, 1, nmax, eps);
  //f[4] = b*bp/(g*(g+1))*F3/C3 - f[0]*f[1];
  f[4] = b*bp/(g*(g+1))*F3/C3 - f[0]*mean[1] - f[1]*mean[0] + mean[0]*mean[1];
  fprintf(stderr, ".");

  f[0] -= mean[0];
  f[1] -= mean[1];
  f[2] -= var[0];
  f[3] -= var[1];
  f[4] -= var[2];
} 
/* -------------------------------------------------------------------------- */

double betaDist(double x, double alpha, double beta)
{
  double tmp, fgam;
  
  fgam = lnGamma(alpha+beta) - lnGamma(alpha) - lnGamma(beta);
  tmp = (alpha-1.0)*log(x) + (beta-1.0)*log(1.0-x) + fgam;

  return exp(tmp);
}
void betaDistAll(double *x, double *x1work, double *x2work,double alpha, double beta,int n,  double *result)
{
  double tmp, fgam;
  int i;
 
#ifdef MKL
  for (i=0;i<n;i++) x2work[i] = 1.0 - x[i];
  vdLn(n,x,x1work);
  vdLn(n,x2work,x2work);
#else
  for (i=0;i<n;i++) {
	x1work[i] = log(x[i]);
   	x2work[i] = log(1.0-x[i]);
  }
#endif /*MKL*/  
  fgam = lnGamma(alpha+beta) - lnGamma(alpha) - lnGamma(beta);

  for (i=0;i<n;i++) {
    	result[i] = (alpha-1.0)*x1work[i] + (beta-1.0)*x2work[i] + fgam;
#ifndef MKL
	result[i] = exp(result[i]);
  }
#else
  }
  vdExp(n,result,result);
#endif /*MKL*/
  return;
}
/* ------------------------------------------------------------------------- */

double jointBetaDist(double x1, double x2, double *b)
{
  double tmp;

  tmp = b[0]
    + (b[1]-1.0)*log(x1)     + (b[2]-1.0)*log(x2) + (b[3]-1.0)*log(1.0-x1)
    + (b[4]-1.0)*log(1.0-x2) + (b[5]-1.0)*log(1.0-x1-x2);

  return exp(tmp);
}
/* ------------------------------------------------------------------------- */

double dirichletDist(double x1, double x2, double b1, double b2, double b3)
{
  double fgam, tmp;
  
  fgam = lnGamma(b1+b2+b3)-lnGamma(b1)-lnGamma(b2)-lnGamma(b3);
  tmp  = (b1-1.0)*log(x1) + (b2-1.0)*log(x2) + (b3-1.0)*log(1.0-x1-x2) + fgam;

  return exp(tmp);
}
/* ------------------------------------------------------------------------- */

double dirichletDistDoran(double x1,double x2,double b1,double b2,double b3)
{
  double tmp, fgam;
  
  fgam = lnGamma(b1+b2+b3)-lnGamma(b1)-lnGamma(b2)-lnGamma(b3);
  tmp =   (b1-1.0)*log(x1) + (b2+b3-1.0)*log(1.0-x1)
        + (b2-1.0)*log(x2) +    (b3-1.0)*log(1.0-x2) + fgam;

  return exp(tmp);
}
/* ------------------------------------------------------------------------- */

double dirichletDistHasse(double x1,double x2,double b1,double b2,double b3)
{
  double tmp, fgam;
  
  fgam = lnGamma(b1+b2+b3)-lnGamma(b1)-lnGamma(b2)-lnGamma(b3);
  tmp =   (b1+b2-1.0)*log(x1) + (b3-1.0)*log(1.0-x1)
        + (b2-1.0)*log(x2)    + (b1-1.0)*log(1.0-x2) + fgam;

  return exp(tmp);
}
/* ------------------------------------------------------------------------- */

#ifdef MDUC_JBETA_ADAPT
void jointBetaFuncODE(int nF, double *a, double *f, void *args)
{
  int    i, flag;
  int    itask,iret;
  long int nst;              // number of CVODE steps
  double dxPre, x = 0.0;

  // dereference variables passed from solver
  statStruct *beta = (statStruct*) args;

  double *xtmp = calloc(1000, sizeof(double));

/*   fprintf(stderr, "beta = {"); */
/*   for (i=0; i<nF; i++) fprintf(stderr, "%8.6f ", a[i]); */
/*   fprintf(stderr, "}\n"); */


  // CVODE tolerances
  double reltol = 1.0e-6;
  double abstol = 1.0e-9;

  // create CVODE vector for computed moments
  N_Vector nv_M = N_VNew_Serial(nF);
  realtype *M   = NV_DATA_S(nv_M);      // data access pointer
  for (i=0; i<nF; i++) M[i] = 0.0;

  // allocate and initialise CVODE object
  void *mem = CVodeCreate(CV_BDF, CV_NEWTON);
#ifdef SUNDIALS_2_3
  flag = CVodeMalloc(mem, jointBetaCVODE_dX, 0.0, nv_M, CV_SS, reltol, &abstol);
  flag = CVodeSetFdata(mem, args); // pass the data onward
  itask = CV_ONE_STEP_TSTOP;
  iret  = CV_TSTOP_RETURN;
#else
  flag  = CVodeInit(mem, jointBetaCVODE_dX, 0.0, nv_M);
  flag  = CVodeSStolerances(mem, reltol, abstol);
  flag  = CVodeSetUserData(mem, args);
  itask = CV_ONE_STEP;
  iret  = CV_TSTOP_RETURN;
#endif /* SUNDIALS_2_3 */
  CVDense(mem, nF);         // use dense matrix solver

  // set stop value for integration
  CVodeSetStopTime(mem, 1.0);
  /* CVodeSetInitStep(mem, 0.01); */

  // integrate the distribution in x-direction
  beta->iX = 1;
  do {
    flag = CVode(mem, 1.0, nv_M, &x, itask);
    if (flag < 0) exit(0);

    // get solver data about last time step
    CVodeGetNumSteps(mem, &nst);
    CVodeGetLastStep(mem, &dxPre);

    xtmp[nst] = x;
    beta->iX++;

    /* fprintf(stderr, "stp = %5li\tdx = %8.2e\tx = %7.4f\t y {", */
/*             nst, dxPre, x); */
/*     for (i=0; i<nF; i++) fprintf(stderr, "%15.6e", Ydata[i]); */
/*     fprintf(stderr, "}\n"); */
  } while (flag != iret);

  beta->nxG = nst+1;
  beta->xG  = (double*) realloc(beta->xG, beta->nxG*sizeof(double));
  for (i=0; i<beta->nxG; i++)
    beta->xG[i] = xtmp[i];

  fprintf(stderr, "moments = {");
  for (i=0; i<nF; i++) fprintf(stderr, "%8.6f ", M[i]);
  fprintf(stderr, "}\n");

  // for the newton solver
  f[0] = M[0] - 1.0;
  f[1] = M[1] - beta->mean[0];
  f[2] = M[2] - beta->mean[1];
  f[3] = M[3] - beta->var[0];
  f[4] = M[4] - beta->var[1];
  f[5] = M[5] - beta->var[2];

  // free the integrator
  N_VDestroy_Serial(nv_M);
  CVodeFree(&mem);
  free(xtmp);
} 
/* -------------------------------------------------------------------------- */

void jointBetaJacODE(int nF, double *a, double **r, void *args)
{
  int      i, j, flag;
  long int nst;              // number of CVODE steps
  double   dxPre, x = 0.0;

  // dereference variables passed from solver
  statStruct *beta = (statStruct*) args;

  // CVODE tolerances
  double reltol = 1.0e-6;
  double abstol = 1.0e-9;

  // create CVODE vector for computed moments
  N_Vector nv_J = N_VNew_Serial(nF*nF);
  realtype *J   = NV_DATA_S(nv_J);      // data access pointer
  for (i=0; i<nF*nF; i++) J[i] = 0.0;

  // allocate and initialise CVODE object
  void *mem = CVodeCreate(CV_BDF, CV_NEWTON);
#ifdef SUNDIALS_2_3
  CVodeMalloc(mem, jointBetaJacCVODE_dX, 0.0, nv_J, CV_SS, reltol, &abstol);
  CVodeSetFdata(mem, args); // pass the data onward
#else
  flag  = CVodeInit(mem, jointBetaJacCVODE_dX, 0.0, nv_J);
  flag  = CVodeSStolerances(mem, reltol, abstol);
  flag  = CVodeSetUserData(mem, args);
#endif /* SUNDIALS_2_3 */
  CVDense(mem, nF*nF);         // use dense matrix solver

  // set stop value for integration
  CVodeSetStopTime(mem, 1.0);
  /* CVodeSetInitStep(mem, 0.01); */

  // integrate the value of dM/dx from y = 0 --> 1
  flag = CVode(mem, 1.0, nv_J, &x, CV_NORMAL);

  // for the newton solver
  for (i=0; i<nF; i++)
    for (j=0; j<nF; j++)
      r[i][j] = J[i*nF + j];

  // free the integrator
  N_VDestroy_Serial(nv_J);
  CVodeFree(&mem);
} 
/* -------------------------------------------------------------------------- */

int jointBetaCVODE_dX(realtype t, N_Vector u, N_Vector udot, void *data)
{
  realtype *Y, *RHS;

  Y   = NV_DATA_S(u);
  RHS = NV_DATA_S(udot);

  realtype RHS_ODE[6], error[6], delta[6];

  jointBetaODE_RHS_X(t, Y, RHS, data);

  //jointBetaODE_RHS_Xgrid(t, Y, RHS, data);

  /* for (int i=0; i<6; i++) */
/*     { */
/*       delta[i] = RHS[i] - RHS_ODE[i]; */
/*       error[i] = fabs(delta[i])/RHS[i]*100; */
/*       if (error[i] > 1.0) */
/*         fprintf(stderr, "x = %8.6f eq[%i]: %15.6e - %15.6e (%6.4f)\n", t, i, */
/*                 RHS[i], delta[i], error[i]); */
/*     } */

  return 0;
}
/* -------------------------------------------------------------------------- */

int jointBetaCVODE_dY(realtype t, N_Vector u, N_Vector udot, void *data)
{
  realtype *Y, *RHS;

  Y   = NV_DATA_S(u);
  RHS = NV_DATA_S(udot);

  jointBetaODE_RHS_Y(t, Y, RHS, data);

  return 0;
}
/* -------------------------------------------------------------------------- */

int jointBetaJacCVODE_dX(realtype t, N_Vector u, N_Vector udot, void *data)
{
  realtype *Y, *RHS;

  Y   = NV_DATA_S(u);
  RHS = NV_DATA_S(udot);

  jointBetaJacODE_RHS_X(t, Y, RHS, data);

  return 0;
}
/* -------------------------------------------------------------------------- */

int jointBetaJacCVODE_dY(realtype t, N_Vector u, N_Vector udot, void *data)
{
  realtype *Y, *RHS;

  Y   = NV_DATA_S(u);
  RHS = NV_DATA_S(udot);

  jointBetaJacODE_RHS_Y(t, Y, RHS, data);

  return 0;
}
/* -------------------------------------------------------------------------- */

void jointBetaODE_RHS_Xgrid(double x, double *y, double *dydx, void *args)
{
  // just do standard integration
  int i, k, n2;
  double x2c, *x2, dx2, dA, Z1, Z2, P;
  statStruct *beta = (statStruct*) args;

  n2 = beta->n2;
  x2 = beta->x2;

  // set the current x value for the column solve
  beta->xi = fmin(0.999, fmax(x, 1.0e-3));

  for (i=0; i<beta->nMom; i++) dydx[i] = 0.0;

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
      
      Z1 = beta->xi;
      Z2 = x2c*(1.0 - beta->xi);
      dA = dx2*(1.0 - beta->xi);

      P = jointBetaDist(Z1, Z2, beta->a);

      // normalisation condition
      dydx[0] += P*dA;

      // compute the means
      dydx[1] += P*Z1*dA; // mean1
      dydx[2] += P*Z2*dA; // mean2

      // compute the variances
      dydx[3] += P*(Z1 - beta->mean[0])*(Z1 - beta->mean[0])*dA; // variance 1
      dydx[4] += P*(Z2 - beta->mean[1])*(Z2 - beta->mean[1])*dA; // variance 2
      dydx[5] += P*(Z1 - beta->mean[0])*(Z2 - beta->mean[1])*dA; // co-variance
    } // end loop over columns
}
/* -------------------------------------------------------------------------- */

void jointBetaODE_RHS_X(double x, double *y, double *dydx, void *args)
{
  int i, flag;
  long int nst;
  statStruct *beta = (statStruct*) args;
  double eta = 0.0, dEtaPre;

  double *ytmp = calloc(1000, sizeof(double));

  // set the current x value for the column solve
  beta->xi = fmin(0.999, fmax(x, 1.0e-3));

  // CVODE tolerances
  double reltol = 1.0e-6;
  double abstol = 1.0e-9;

  // set vector for value of dM/dx to be returned
  N_Vector nv_dMdx = N_VNew_Serial(beta->nMom);
  realtype *dMdx   = NV_DATA_S(nv_dMdx);
  for (i=0; i<beta->nMom; i++) dMdx[i] = 0.0;

  void *mem = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeMalloc(mem, jointBetaCVODE_dY, 0.0, nv_dMdx, CV_SS, reltol, &abstol);
  CVodeSetFdata(mem, args);
  CVDense(mem, beta->nMom);

  // set the value of y to integrate to
  CVodeSetStopTime(mem, 1.0);
/*   CVodeSetInitStep(mem, 0.01); */

  // integrate the value of dM/dx from y = 0..1
  /* flag = CVode(mem, 0.9999, cvodeY, &xcur, CV_NORMAL); */
  do {
    flag = CVode(mem, 1.0, nv_dMdx, &eta, CV_ONE_STEP_TSTOP);
    if (flag < 0) exit(0);

    // get solver data about last time step
    CVodeGetNumSteps(mem, &nst);
    CVodeGetLastStep(mem, &dEtaPre);

    ytmp[nst] = eta;

  } while (flag != CV_TSTOP_RETURN);

  int offset, nyGtot = 0;
  beta->nyG = realloc(beta->nyG, (beta->iX+1)*sizeof(int));
  beta->nyG[beta->iX] = nst+1;
  for (i=0; i<beta->iX; i++) nyGtot += beta->nyG[i];

  offset = nyGtot;
  nyGtot += beta->nyG[beta->iX];

  beta->yG = realloc(beta->yG, nyGtot*sizeof(double));
  for (i=0; i<beta->nyG[beta->iX]; i++) beta->yG[offset+i] = ytmp[i];

  for (i=0; i<beta->nMom; i++) dydx[i] = dMdx[i];

  // free the integrator memory
  N_VDestroy_Serial(nv_dMdx);
  CVodeFree(&mem);
  free(ytmp);
}
/* -------------------------------------------------------------------------- */

void jointBetaODE_RHS_Y(double x, double *y, double *dydx, void *args)
{
  statStruct *p = (statStruct*) args;
  double eta, Z1, Z2, P;
 
  // bound current y co-ordinate
  eta = fmin(0.999, fmax(x, 1.0e-3));

  // current values in untransformed space
  Z1 = p->xi;
  Z2 = eta*(1.0 - p->xi);
  // compute the distribution at this point
  P = (1.0 - p->xi)*jointBetaDist(Z1, Z2, p->a);

  // value of d/dy[dM/dx]
  dydx[0] = P;
  dydx[1] = P*Z1;
  dydx[2] = P*Z2;
  dydx[3] = P*(Z1 - p->mean[0])*(Z1 - p->mean[0]);
  dydx[4] = P*(Z2 - p->mean[1])*(Z2 - p->mean[1]);
  dydx[5] = P*(Z1 - p->mean[0])*(Z2 - p->mean[1]);
}
/* -------------------------------------------------------------------------- */

void jointBetaJacODE_RHS_X(double x, double *y, double *dydx, void *args)
{
  int i, flag;
  long int nst;
  statStruct *beta = (statStruct*) args;
  double eta = 0.0, dEtaPre;

  // set the current x value for the column solve
  beta->xi = fmin(0.999, fmax(x, 1.0e-3));

  // CVODE tolerances
  double reltol = 1.0e-6;
  double abstol = 1.0e-9;

  // set vector for value of dM/dx to be returned
  N_Vector nv_dJdx = N_VNew_Serial(beta->nMom*beta->nMom);
  realtype *dJdx   = NV_DATA_S(nv_dJdx);
  for (i=0; i<beta->nMom*beta->nMom; i++) dJdx[i] = 0.0;

  void *mem = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeMalloc(mem, jointBetaJacCVODE_dY, 0.0, nv_dJdx, CV_SS, reltol, &abstol);
  CVodeSetFdata(mem, args);
  CVDense(mem, beta->nMom*beta->nMom);

  // set the value of y to integrate to
  CVodeSetStopTime(mem, 1.0);
  /* CVodeSetInitStep(mem, 0.01); */

  // integrate the value of dM/dx from y = 0 --> 1
  flag = CVode(mem, 1.0, nv_dJdx, &eta, CV_NORMAL);

  for (i=0; i<beta->nMom*beta->nMom; i++) dydx[i] = dJdx[i];

  // free the integrator memory
  N_VDestroy_Serial(nv_dJdx);
  CVodeFree(&mem);
}
/* -------------------------------------------------------------------------- */

void jointBetaJacODE_RHS_Y(double x, double *y, double *dydx, void *args)
{
  int i, j;
  statStruct *p = (statStruct*) args;
  double eta, Z1, Z2, P, d[6], f[6];
 
  // bound current y co-ordinate
  eta = fmin(0.999, fmax(x, 1.0e-3));

  // current values in untransformed space
  Z1 = p->xi;
  Z2 = eta*(1.0 - p->xi);

  // compute the distribution at this point
  P = (1.0 - p->xi)*jointBetaDist(Z1, Z2, p->a);

  // compute distribution with extra terms from derivatives
  d[0] = P;
  d[1] = P*log(Z1);
  d[2] = P*log(Z2);
  d[3] = P*log(1.0 - Z1);
  d[4] = P*log(1.0 - Z2);
  d[5] = P*log(1.0 - Z1 - Z2);

  // compute function values
  f[0] = 1.0;
  f[1] = Z1;
  f[2] = Z2;
  f[3] = (Z1 - p->mean[0])*(Z1 - p->mean[0]);
  f[4] = (Z2 - p->mean[1])*(Z2 - p->mean[1]);
  f[5] = (Z1 - p->mean[0])*(Z2 - p->mean[1]);

  // fill the Jacobian --> d/dy[dJ/dx]
  for (i=0; i<6; i++)
    for (j=0; j<6; j++)
      dydx[6*i + j] = f[i]*d[j];
}
/* -------------------------------------------------------------------------- */
#endif /* MDUC_JBETA_ADAPT */
