#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hypergeometric.h"
#include "series.h"
#include "statsUtils.h"

double hyper(int p, int q, double *a, double *b, double x)
{
  int i, nmax = 1000;
  double h, eps = 1.0e-9;
 
  /* more special cases need to be implemented !!!*/
  /* handle special cases */
  if (p == 0)
    {
      if (q == 0) return exp(x);
      /* else if (q == 1) return hyper_0f1(b, x); */
      else
        h = hyperShanks(p, q, a, b, x, nmax, eps);
    }
  else if (p == 1)
    {
      if (q == 0)
        h = 1.0/pow(1.0-x, a[0]);
      else
        h = hyperShanks(p, q, a, b, x, nmax, eps);
    }
  else if (p == 2)
    {
      if (q == 1) // standard Gauss hypergeometric series
        {
          if (x > 0)
            h = hyper_2F1(a[0], a[1], b[0], x);
          else
            {
              /* use 1-z transformation to make positive 
                    Abramowitz and Stegun, 15.3.7
                 This is not quite right yet, as there are additional poles
                 that need to be accounted for */
              double aa = a[0], bb = a[1], c = b[0];
              double h1, h2;

              // first term
              double gn_1[2] = {c, c-aa-bb};
              double gd_1[2] = {c-aa, c-bb};

              h1  = hyper_2F1(aa, bb, aa+bb-c+1, 1-x);
              h1 *= gammaProd(2, 2, gn_1, gd_1);

              // second term
              double gn_2[2] = {c, aa+bb-c};
              double gd_2[2] = {aa, bb};

              h2 = hyper_2F1(c-aa, c-bb, c-aa-bb+1, 1-x);
              h1 *= pow(1-x,c-aa-bb)*gammaProd(2, 2, gn_2, gd_2);

              h = h1 + h2;
            }
        }
    }
  else if (p == q+1) // balanced series
    {
      if (p == 3)
        h = hyper_3F2(a[0], a[1], a[2], b[0], b[1], x);
      else /* general implementation */
        {
          if (x == 0.0)
            h = 1.0;
          else if (x >= 0 && x < 1.0)
            h = hyperShanks(p, q, a, b, x, nmax, eps);
          else
            {
              fprintf(stderr, "hyper: p=q+1 not yet implemented for x >= 1.0\n");
              exit(2);
            }
        }
    } /*  end p == q+1 */
  else if (p > q+1)
    {
      h = hyperShanks(p, q, a, b, x, nmax, eps);
    }
 
  return h;
}
/* -------------------------------------------------------------------------- */

double hyper_2F1(double a, double b, double c, double x)
{
  int    nmax = 1000;
  double h, s, eps = 1.0e-9;

  if (x == 0.0)
    h = 1.0;
  else if (x >= 0.0 && x < 1.0)
    {
      /* within unit circle, standard espilon algorithm with Shanks
         transformation works quite well*/
      double an[2] = {a, b};
      double bn[1] = {c};
      h = hyperShanks(2, 1, an, bn, x, nmax, eps);
    }
  else if (x > 1.0)
    {
      fprintf(stderr, "hyper_2F1: x > 1.0 doesn't work yet -- need complex\n");
      exit(2);
      /* use 1/z transformation
         Abramowitz and Stegun, 15.3.7 */
      double h1, h2;
      // first term
      double an[2] = {a, 1 - c + a};
      double bn[1] = {1 - b + a};

      // coefficients for gamma function
      double gNum[2] = {c, b - a};
      double gDen[2] = {b, c - a};
      h1  = hyperShanks(2, 1, an, bn, 1/x, nmax, eps);
      h1 *= pow(-x,-a)*gammaProd(2, 2, gNum, gDen);

      // second term
      an[0] = b;
      an[1] = 1 - c + b;
      bn[0] = 1 - a + b;

      // coefficients for gamma function
      gNum[0] = c;
      gNum[1] = a - b;
      gDen[0] = a;
      gDen[1] = c - b;
      h2  = hyperShanks(2, 1, an, bn, 1/x, nmax, eps);
      h2 *= pow(-x,-b)*gammaProd(2, 2, gNum, gDen);

      h = h1 + h2;
    }
  else /* x == 1.0 */
    {
      /* analytic continuation */
      s = c - a - b;
      double fNum[2] = {c, s};
      double fDen[2] = {c - a, c - b};
      
      if (s > 0)
        h = gammaProd(2, 2, fNum, fDen);
      else
        {
          fprintf(stderr, "hyper_2F1: s < 0 not implemented yet\n");
          exit(2);
        }
    }
  return h;
}
/* -------------------------------------------------------------------------- */

double hyper_3F2(double a, double b, double c, double e, double f, double x)
{
  int    nmax = 2000;
  double h, s, eps = 1.0e-9;

  double an[3] = {a, b, c};
  double bn[2] = {e, f};

  if (x == 0.0) /* simple case */
    h = 1.0;
  else if (x >= 0.0 && x < 1.0)
    {
      /* within unit circle, standard espilon algorithm with Shanks
         transformation works quite well*/
      h = hyperShanks(3, 2, an, bn, x, nmax, eps); 
    }
  else if (x > 1.0)
    {
      /* Need to implement transformation for x > 1.0, using x -> 1/x 
         Transformation exists
         http://functions.wolfram.com/HypergeometricFunctions/
         HypergeometricPFQ/06/01/05/02/0004/
      */
      fprintf(stderr, "hyper_3F2: not yet implemented for x > 1.0\n");
      exit(2);
    }
  else /* x == 1.0 */
    {
      // check whether if will converge
      s = e + f - a - b - c;

      if (s >= 0.0) /* logarithmic convergence on unit circle */
        {
          //h = hyperLevin(3, 2, an, bn, x, nmax, eps);
          //h = hyperShanks(3, 2, an, bn, x, nmax, eps);
          h = hypersum(3, 2, an, bn, x, nmax, eps);
        }
      else
        {
          /* not convergent without transformation */
          double fgam;

          // compute the coefficients for the unit argument transformation
          double ac[3] = {e - c, f - c, s};
          double bc[2] = {a + s, b + s};

          // compute the leading coefficient (rational of Gamma functions)
          double nGam[3] = {e, f, s};
          double dGam[3] = {a + s, b + s, c};
          fgam  = gammaProd(3, 3, nGam, dGam);

          h = hypersum(3, 2, ac, bc, x, nmax, eps);
          
          // compute the series using Levin acceleration
          //h  = hyperLevin(3, 2, ac, bc, x, nmax, eps);
          //h  = hyperShanks(3, 2, ac, bc, x, nmax, eps);
          h *= fgam;
        }
    }
  return h;
}
/* -------------------------------------------------------------------------- */

double appell_F1(double a, double b, double bp, double g, double x, double y,
                 int nmax, double eps)
{
  int    m, n, i;
  double term, sum_pre, sum_in, sum = 1.0;
  double a_out   = a;
  double b_out   = b;
  double g_out   = g;
  double a_in[2] = {a, bp};
  double b_in[1] = {g};

  sum  = sum_pre =  hyper(2, 1, a_in, b_in, y);
  term = 1.0;
  for (m=1; m<nmax; m++)
    {
      a_in[0] = a + m;
      a_in[1] = bp;
      b_in[0] = g + m;

      sum_in = hyper(2, 1, a_in, b_in, y);
      //term  *= (a+m-1)*(b+m-1)*x/(m*(g+m-1));
      term  *= a_out*b_out*x/(m*g_out);

      sum += term*sum_in;
      if (fabs(sum - sum_pre) < eps)
        {
          fprintf(stderr, "Appell converged to (%e) in %i steps\n", eps, m);
          return sum;
        }

      sum_pre = sum;
      a_out++;
      b_out++;
      g_out++;
    }
  fprintf(stderr, "Error: appell_F1 did not converge in %i iterations\n",nmax);
  return sum;
}
/* ------------------------------------------------------------------------- */

double appell_F2(double a, double b, double bp, double g, double gp,
                 double x, double y, int nmax, double eps)
{
  int    m, n, i;
  double term, fgam, sum_pre, sum_in, sum = 1.0;
  double a_in[2], b_in[1], a_out, b_out, g_out;

  if (x + y < 1) // convergent in the unit circle, |x| + |y| < 1
    {
      a_out = a;
      b_out = b;
      g_out = g;
      a_in[0] = a;
      a_in[1] = bp;
      b_in[0] = gp;

      // compute the first term in the series
      sum  = sum_pre =  hyper(2, 1, a_in, b_in, y);
      term = 1.0;
      for (m=1; m<nmax; m++)
        {
          // update the inner sum terms
          a_in[0] = a + m;
          a_in[1] = bp;
          b_in[0] = gp;

          // compute the inner sum
          sum_in = hyper(2, 1, a_in, b_in, y);

          // update the outer coefficient
          //term  *= (a+m-1)*(b+m-1)*x/(m*(g+m-1));
          term  *= a_out*b_out*x/(m*g_out);

          sum += term*sum_in; // compute the partial sum
          if (fabs(sum - sum_pre) < eps) // check for convergence
            {
              //fprintf(stderr, "Appell converged to (%e) in %i steps\n", eps, m);
              return sum;
            }

          // store previous sum and update outer parameters
          sum_pre = sum;
          a_out++;
          b_out++;
          g_out++;
        }
      fprintf(stderr, "Err: appell_F2 did not converge in %i iterations\n",nmax);
    }
/*   else if (x == 1.0) */
  else if (x > 1.0)
    {
      /* Use continuation from Hahne (Eq 12).  Only the first and third terms
         will be non-zero. */
      double h1, h2;

      // compute the coefficients for the first term
      double a_3F2_1[3] = {a + 1 - g + b - bp, gp - bp - a + b, gp - bp};
      double b_3F2_1[2] = {gp - bp + b, 1 - g + b + gp - bp};

      // leading values for gamma product
      double fNum_1[4] = {g, gp, g - b - gp + bp, gp - bp - a + b};
      double fDen_1[4] = {bp, g - b, gp - bp + b, g - a};

      // compute the series
      h1   = hyper(3, 2, a_3F2_1, b_3F2_1, x); // compute 3F2(1)
      fgam = gammaProd(4, 4, fNum_1, fDen_1);
      h1  *= fgam;//*pow(-1,a+bp-gp);
      //fprintf(stderr, "h1 = %.10f\n", h1);

      /* compute the coefficients for the third term (simultaneous interchange
         of b and bp, g and gp*/
      double a_3F2_2[3] = {a + 1 - gp + bp - b, g - b - a + bp, g - b};
      double b_3F2_2[4] = {g - b + bp, 1 - gp + bp + g - b};

      // leading values for gamma product
      double fNum_2[4] = {g, gp, gp - bp - g + b, g - b - a + bp};
      double fDen_2[4] = {b, gp - bp, g - b + bp, gp - a};

      // compute the series
      h2   = hyper(3, 2, a_3F2_2, b_3F2_2, x); // compute 3F2(1)
      fgam = gammaProd(4, 4, fNum_2, fDen_2);
      h2  *= fgam;//*pow(-1,a+b-g);
      //fprintf(stderr, "h2 = %.10f\n", h2);

      sum = h1 + h2;
    }
  else if (y == 1.0)
    {
      // compute the coefficients according to relation in Vidunas 2009
      double a_3F2[3] = {a, b, a - gp + 1};
      double b_3F2[2] = {g, a + bp - gp + 1};
/*       double a_3F2[3] = {a, b, gp - a - 1}; */
/*       double b_3F2[2] = {g, gp - a - bp - 1}; */

      // leading values for gamma product
      double fNum[2] = {gp, gp - a - bp};
      double fDen[2] = {gp - a, gp - bp};

      sum   = hyper(3, 2, a_3F2, b_3F2, x); // compute 3F2(1)
      fgam  = gammaProd(2, 2, fNum, fDen);
      sum  *= fgam;
    }
  return sum;
}
/* ------------------------------------------------------------------------- */

double appell_F3(double a, double ap, double b, double bp, double g,
                 double x, double y, int nmax, double eps)
{
  int    m, n, i;
  double term, sum_pre, sum_in, sum = 1.0;
  double a_in[2], b_in[1], a_out, b_out, g_out;

  if (x < 1.0 && y < 1.0) // absolutely convergent for |x| < 1, |y| < 1
    {
      a_out = a;
      b_out = b;
      g_out = g;
      a_in[0] = ap;
      a_in[1] = bp;
      b_in[0] = g;

      sum  = sum_pre =  hyper(2, 1, a_in, b_in, y);
      term = 1.0;
      for (m=1; m<nmax; m++)
        {
          a_in[0] = ap;
          a_in[1] = bp;
          b_in[0] = g + m;

          sum_in = hyper(2, 1, a_in, b_in, y);

          //term  *= (a+m-1)*(b+m-1)*x/(m*(g+m-1));
          term  *= a_out*b_out*x/(m*g_out);

          sum += term*sum_in;
          if (fabs(sum - sum_pre) < eps)
            return sum;

          sum_pre = sum;
          a_out++;
          b_out++;
          g_out++;
        }
      fprintf(stderr, "Err: appell_F3 did not converge in %i iterations\n",nmax);
    }
  else if (y == 1.0)
    {
      double fgam;

      // compute the coefficients according to relation in Vidunas 2009
      double a_3F2[3] = {a, b, g - ap - bp};
      double b_3F2[2] = {g - ap, g - bp};

      // leading values for gamma rational
      double fNum[2] = {g, g - ap - bp};
      double fDen[2] = {g - ap, g - bp};

      sum   = hyper(3, 2, a_3F2, b_3F2, x);
      fgam  = gammaProd(2, 2, fNum, fDen);
      sum  *= fgam;
    }
  return sum;
}
/* ------------------------------------------------------------------------- */

double appell_F4(double a, double b, double g, double gp, double x, double y,
                 int nmax, double eps)
{
  int    m, n, i;
  double term, sum_pre, sum_in, sum = 1.0;
  double a_in[2], b_in[1], a_out, b_out, g_out;

  a_out = a;
  b_out = b;
  g_out = g;
  a_in[0] = a;
  a_in[1] = b;
  b_in[0] = gp;

  sum  = sum_pre =  hyper(2, 1, a_in, b_in, y);
  term = 1.0;
  for (m=1; m<nmax; m++)
    {
      a_in[0] = a + m;
      a_in[1] = b + m;
      b_in[0] = gp;

      sum_in = hyper(2, 1, a_in, b_in, y);
      //term  *= (a+m-1)*(b+m-1)*x/(m*(g+m-1));
      term  *= a_out*b_out*x/(m*g_out);

      sum += term*sum_in;
      if (fabs(sum - sum_pre) < eps)
        return sum;

      sum_pre = sum;
      a_out++;
      b_out++;
      g_out++;
    }
  fprintf(stderr, "Error: appell_F4 did not converge in %i iterations\n",nmax);
  return sum;
}
/* ------------------------------------------------------------------------- */

double hypersum(int p, int q, double *a, double *b, double x, int nmax,
                double eps)
{
  int n, i;
  double num, den, term, sum_pre, sum = 1.0;
  double *atmp = (double*) calloc(p, sizeof(double));
  double *btmp = (double*) calloc(q, sizeof(double));
  // copy parameters so that they can be modified in the series loop
  for (i=0; i<p; i++) atmp[i] = a[i];
  for (i=0; i<q; i++) btmp[i] = b[i];

  term = sum_pre = 1.0;
  for (n=1; n<nmax; n++)
    {
      num = den = 1.0;
      for (i=0; i<p; i++) num *= atmp[i];
      for (i=0; i<q; i++) den *= btmp[i];

      term *= num*x/(n*den);
      sum  += term;
      //fprintf(stderr, "%8.6f\n", sum);
      if (fabs(sum - sum_pre) < eps)
        {
          free(atmp);
          free(btmp);
          return sum;
        }

      sum_pre = sum;
      for (i=0; i<p; i++) atmp[i] += 1.0;
      for (i=0; i<q; i++) btmp[i] += 1.0;
    }
  free(atmp);
  free(btmp);
  fprintf(stderr, "Error: hypersum did not converge in %i iterations\n",nmax);
  return sum;
}
/* ------------------------------------------------------------------------- */

double hyperShanks(int p, int q, double *a, double *b, double x, int nmax,
                   double eps)
{
  int    i, n, ncv, cvgd;
  double val, lastval, lasteps, term, sum, num, den;

  // temporary storage space
  double *atmp = (double*) calloc(p, sizeof(double));
  double *btmp = (double*) calloc(q, sizeof(double));
  double *e    = (double*) calloc(nmax, sizeof(double));

  // copy parameters so that they can be modified in the series loop
  for (i=0; i<p; i++) atmp[i] = a[i];
  for (i=0; i<q; i++) btmp[i] = b[i];

  ncv = cvgd = 0;
  n = 1;
  sum = term = e[0] = 1.0;
  while (!cvgd && n < nmax)
    {
      num = den = 1.0;
      for (i=0; i<p; i++) num *= atmp[i];
      for (i=0; i<q; i++) den *= btmp[i];
      if (n > 0) term *= num*x/(n*den);
      sum += term;
      val = ShanksNext(sum,&n,&ncv,eps,e,&lastval,&lasteps,&cvgd);

      for (i=0; i<p; i++) atmp[i] += 1.0;
      for (i=0; i<q; i++) btmp[i] += 1.0;
    }
  free(e);
  free(atmp);
  free(btmp);

  if (!cvgd)
    fprintf(stderr, "epsilon algorithm not converged in %i steps!\n", n);
  return val;
}
/* ------------------------------------------------------------------------- */

double hyperLevin(int p, int q, double *a, double *b, double x, int nmax,
                  double eps)
{
  int    i, n, ncv, cvgd;
  double val, omega, lastval, lasteps, term, sum, num, den;
  double beta = 1.0;

  // temporary storage space
  double *atmp  = (double*) calloc(p, sizeof(double));
  double *btmp  = (double*) calloc(q, sizeof(double));
  double *numer = (double*) calloc(nmax, sizeof(double));
  double *denom = (double*) calloc(nmax, sizeof(double));

  // copy parameters so that they can be modified in the series loop
  for (i=0; i<p; i++) atmp[i] = a[i];
  for (i=0; i<q; i++) btmp[i] = b[i];

  ncv = cvgd = 0;
  n = 1;
  sum = term = numer[0] = denom[0] = 1.0;
  while (!cvgd && n < nmax)
    {
      num = den = 1.0;
      for (i=0; i<p; i++) num *= atmp[i];
      for (i=0; i<q; i++) den *= btmp[i];
      if (n > 0) term *= num*x/(n*den);
      sum += term;
      omega = (beta + n)*term;     // u transformation
      val = LevinNext(sum,omega,beta,&n,&ncv,eps,numer,denom,&lastval,&lasteps,
                      &cvgd);

      for (i=0; i<p; i++) atmp[i] += 1.0;
      for (i=0; i<q; i++) btmp[i] += 1.0;
    }
  free(numer);
  free(denom);
  free(atmp);
  free(btmp);

  if (!cvgd)
    fprintf(stderr, "Levin algorithm not converged in %i steps!\n", n);
  return val;
}
/* ------------------------------------------------------------------------- */
