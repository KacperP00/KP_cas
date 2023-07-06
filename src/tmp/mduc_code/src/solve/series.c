#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "series.h"

double ShanksNext(double sum, int *n, int *ncv, double eps, double *e,
                  double *lastval, double *lasteps, int *cnvgd)
{
  int j;
  double diff, temp1, temp2, val;
  
  e[*n] = sum;
  temp2 = 0.0;
  for (j=*n; j>0; j--)
    {
      temp1 = temp2;
      temp2 = e[j-1];
      diff  = e[j] - temp2;
      if (fabs(diff) <= 10*DBL_MIN)
        e[j-1] = DBL_MAX;
      else
        e[j-1] = temp1 + 1.0/diff;
    }
  *n = *n + 1;
  val = (*n & 1) ? e[0] : e[1]; /* cases of n even or odd */
  if (fabs(val) > 0.01*DBL_MAX) val = *lastval;
  *lasteps = fabs(val - *lastval);
  if (*lasteps > eps) *ncv = 0;
  else *ncv = *ncv + 1;
  if (*ncv >= 3) *cnvgd = 1;
  return (*lastval = val);
}
/* -------------------------------------------------------------------------- */

double LevinNext(double sum, double omega, double beta, int *n, int *ncv,
                 double eps, double *num, double *den, double *lastval,
                 double *lasteps, int *cnvgd)
{
  int j;
  double fact,ratio,term,val;
  term = 1.0/(beta + *n);
  den[*n] = term/omega;
  num[*n] = sum*den[*n];
  if (n > 0)
    {
      ratio = (beta + *n - 1)*term;
      for (j=1; j<=*n; j++)
        {
          fact = (*n-j+beta)*term;
          num[*n-j] = num[*n-j+1] - fact*num[*n-j];
          den[*n-j] = den[*n-j+1] - fact*den[*n-j];
          term = term*ratio;
        }
    }
  *n = *n+1;
  val = fabs(den[0]) < 10.0*DBL_MIN ? *lastval : num[0]/den[0];
  *lasteps = fabs(val-*lastval);
  if (*lasteps <= eps) *ncv = *ncv + 1;
  if (*ncv >= 2) *cnvgd = 1;
  return (*lastval = val);
}
/* -------------------------------------------------------------------------- */
