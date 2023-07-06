#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "metric.h"
#include "solveUtils.h"
#include "utils.h"

void initMetric(flameStruct *flame)
{
  int i, j, k, jst, kst;
  double *x1 = flame->x1;
  double *x2 = flame->x2;

  // set the stencil width
  flame->stw = flame->order/2;

  // allocate and fill operators for x-coordinate
  flame->grad.x  = (double**) calloc(flame->n1, sizeof(double*));
  flame->grad.xx = (double**) calloc(flame->n1, sizeof(double*));
  for (j=0; j<flame->n1; j++)
    {
      flame->grad.x[j]  = (double*) calloc(flame->order+1, sizeof(double));
      flame->grad.xx[j] = (double*) calloc(flame->order+1, sizeof(double));

      jst = stencil_x(j, flame);
      fd_coeff(1, flame->order, &x1[jst], x1[j], flame->grad.x[j]);
      fd_coeff(2, flame->order, &x1[jst], x1[j], flame->grad.xx[j]);
    }

  /* exit here if only 1D */
  if (flame->nDim < 2)
    return;
  
  // allocate and fill operators for x-coordinate
  flame->grad.y  = (double**) calloc(flame->n2, sizeof(double*));
  flame->grad.yy = (double**) calloc(flame->n2, sizeof(double*));
  for (k=0; k<flame->n2; k++)
    {
      flame->grad.y[k]  = (double*) calloc(flame->order+1, sizeof(double));
      flame->grad.yy[k] = (double*) calloc(flame->order+1, sizeof(double));

      kst = stencil_y(k, flame);
      fd_coeff(1, flame->order, &x2[kst], x2[k], flame->grad.y[k]);
      fd_coeff(2, flame->order, &x2[kst], x2[k], flame->grad.yy[k]);
    } /* end loop over rows */
}
/* -------------------------------------------------------------------------- */

int stencil_x(int j, flameStruct *flame)
{
  return IMIN(flame->n1-1-flame->order, IMAX(0,j-flame->stw));
}
/* -------------------------------------------------------------------------- */

int stencil_y(int k, flameStruct *flame)
{
  return IMIN(flame->n2-1-flame->order, IMAX(0,k-flame->stw));
}
/* -------------------------------------------------------------------------- */

double grad_x(int st,int k,int j,int i,int nV,double *Y,flameStruct *flame)
{
  int    s, jp;
  double d = 0.0;

  for (s=0, jp=st; s<flame->order+1; s++, jp++)
    d += flame->grad.x[j][s]*Y[(k*flame->n1+jp)*nV + i];

  return d;
}
/* -------------------------------------------------------------------------- */

double grad_xx(int st,int k,int j,int i,int nV,double *Y,flameStruct *flame)
{
  int    s, jp;
  double d = 0.0;

  for (s=0, jp=st; s<flame->order+1; s++, jp++)
    d += flame->grad.xx[j][s]*Y[(k*flame->n1+jp)*nV + i];

  return d;
}
/* -------------------------------------------------------------------------- */

double grad_y(int st,int k,int j,int i,int nV,double *Y,flameStruct *flame)
{
  int    s, kp;
  double d = 0.0;

  for (s=0, kp=st; s<flame->order+1; s++, kp++)
    d += flame->grad.y[k][s]*Y[(kp*flame->n1+j)*nV + i];

  return d;
}
/* -------------------------------------------------------------------------- */

double grad_yy(int st,int k,int j,int i,int nV,double *Y,flameStruct *flame)
{
  int    s, kp;
  double d = 0.0;

  for (s=0, kp=st; s<flame->order+1; s++, kp++)
    d += flame->grad.yy[k][s]*Y[(kp*flame->n1+j)*nV + i];

  return d;
}
/* -------------------------------------------------------------------------- */

void fd_coeff(int deriv, int order, double *x1, double x, double *coeff)
{
  int i, j;
  int n = order+1;
  
  double **a = (double**) calloc(n, sizeof(double));
  double **b = (double**) calloc(n, sizeof(double));
  for (i=0; i<n; i++)
    {
      a[i] = (double*) calloc(n, sizeof(double));
      b[i] = (double*) calloc(n, sizeof(double));
    }

  /* fill the Taylor table matrix */
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      a[i][j] = pow(x1[i] - x, j)/factorial((IMAX(j,1)));

  /* get the inverse matrix */
  invertMatrix(n, a, b);

  for (i=0; i<n; i++)
    coeff[i] = b[deriv][i];

  for (i=0; i<n; i++)
    {
      free(a[i]);
      free(b[i]);
    }
  free(a);
  free(b);
}
/* -------------------------------------------------------------------------- */
