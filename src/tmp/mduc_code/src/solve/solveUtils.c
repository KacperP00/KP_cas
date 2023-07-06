#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "solveUtils.h"

/* DEFINE GLOBAL VARIABLES */
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SWAP(a,b) {temp=(a); (a)=(b); (b)=temp;}

static int    MAXITS = 500;
static double STPMAX = 100.0;
static double TINY   = 1.0e-20;

static double EPS    = 1.0e-9;
static double TOLF   = 1.0e-4;
static double TOLX   = 1.0e-7;
static double TOLMIN = 1.0e-9;
static double TOLSVD = 1.0e-5;
static double ALF    = 1.0e-4;

/* Function returns the central difference of 1st derivative 
 (2nd order accurate) */
double cGrad(double y_im1, double y_i, double y_ip1, double dx1, double dx2)
{
  double den, dx1_2, dx2_2;
        
  dx1_2 = dx1*dx1;
  dx2_2 = dx2*dx2;
        
  return (y_ip1*dx1_2 + y_i*(dx2_2-dx1_2) - y_im1*(dx2_2))/(dx1*dx2*(dx1+dx2));
}
/* -------------------------------------------------------------------------- */


/* Function returns the forward difference of 1st derivative 
 (1st order accurate) */
double fGrad(double y_i, double y_ip1, double dx)
{
  return (y_ip1 - y_i)/dx;
}
/* -------------------------------------------------------------------------- */


/* Function returns the backward difference of 1st derivative 
 (1st order accurate) */
double bGrad(double y_i, double y_im1, double dx)
{
  return (y_i - y_im1)/dx;
}
/* -------------------------------------------------------------------------- */


/* Function returns the central difference of 2nd derivative 
 (2nd order accurate) */
double cGrad2(double y_im1, double y_i, double y_ip1, double dx1, double dx2)
{
  return 2.0/(dx1+dx2)*((y_ip1-y_i)/dx2 - (y_i-y_im1)/dx1);
}
/* -------------------------------------------------------------------------- */


/* Function returns the forward difference of 2nd derivative 
 (1st order accurate) */
double fGrad2(double y_i, double y_ip1, double y_ip2, double dx1, double dx2)
{
  return 2.0/(dx1+dx2)*((y_ip2-y_ip1)/dx2 - (y_ip1-y_i)/dx1);
}
/* -------------------------------------------------------------------------- */


/* Function returns the backward difference of 2nd derivative 
 (1st order accurate) */
double bGrad2(double y_i, double y_im1, double y_im2, double dx1, double dx2)
{
  return 2.0/(dx1+dx2)*((y_i-y_im1)/dx1 - (y_im1-y_im2)/dx2);
}
/* -------------------------------------------------------------------------- */

void solveTriDiag(int n, double *a, double *b, double *c, double *u, double *r)
{
  int i;
  double m;
  
  // eliminate a in forward substiution
  for (i=1; i < n; i++)
    {
      m = a[i]/b[i-1];
      
      b[i] -= m*c[i-1];
      r[i] -= m*r[i-1];
    }
  
  // compute the value of u_n
  u[n-1] = r[n-1]/b[n-1];
  
  // update solution with backward substitution
  for (i=n-2; i > 0; i--)
    u[i] = (r[i] - c[i]*u[i+1])/b[i];
}
/* -------------------------------------------------------------------------- */

void triDiagNR(int n, double *a, double *b, double *c, double *u, double *r)
{
  int i;
  double bet, *gam;
  
  gam = (double*) calloc(n, sizeof(double));
  
  if (b[0] == 0.0)
    {
      fprintf(stderr, "ERROR: trivial row in tridiagonal solve\n");
      exit(2);
    }
  
  u[0] = r[0]/(bet=b[0]);
  
  // forward substitution
  for (i=1; i<n; i++)
    {
      gam[i] = c[i-1]/bet;
      bet = b[i] - a[i]*gam[i];
      if (bet == 0.0)
        {
          fprintf(stderr, "ERROR: zero pivot in tridiagonal solve\n");
          exit(2);
        }
      u[i] = (r[i] - a[i]*u[i-1])/bet;
    }
  
  for (i=n-2; i>0; i--)
    u[i] -= gam[i+1]*u[i+1];
  
  free(gam);
}
/* -------------------------------------------------------------------------- */

double brent(double (*funcEval)(double, void *), double x1, double x2,
             double tol, void *args)
{
  int j, JMAX=100;
  double EPS=3.0e-8;
  double a=x1, b=x2, c=x2, d,e,min1,min2;
  double fa=(*funcEval)(a,args),fb=(*funcEval)(b,args),fc,p,q,r,s,tol1,xm;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    {
      fprintf(stderr, "root must be bracketed in brent\n");
      exit(0);
    }
      
  fc=fb;
  for (j=1; j < JMAX; j++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    } /* end if */
    if (fabs(fc) < fabs(fb)) {
      a=b; b=c; c=a;
      fa=fb; fb=fc; fc=fa;
    } /* end if */
    tol1=2.0*EPS*fabs(b)+0.5*tol; /* check convergence */
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
        p=2.0*xm*s;
        q=1.0-s;
      } else {
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      } /* end if */
      if (p > 0.0) q=-q; // check whether in bounds
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=d;
        d=p/q;
      } else {
        d=xm;
        e=d;
      } /* end if */
    } else {
      d=xm;
      e=d;
    } /* end if */
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb=(*funcEval)(b,args);
  } /* end for j */
  fprintf(stderr, "Maximum number of iterations in brent exceeded\n");
  exit(0);
  return 0.0;
}
/* -------------------------------------------------------------------------- */

void newt(double *x, int n, int *check, void *args, 
          void (*vecFunc)(int, double *, double *, void *),
          void (*jacFunc)(int, double *, double **, void *))
{
  int i, its, j, *indx;
  double d, den, f, fold, stpmax, sum, temp, test;
  double *g, *p, *xold, *fvec, **fjac;
        
  /* allocate memory for vectors */
  indx = (int*)     calloc(n, sizeof(int));
  g    = (double*)  calloc(n, sizeof(double));
  p    = (double*)  calloc(n, sizeof(double));
  xold = (double*)  calloc(n, sizeof(double));
  fvec = (double*)  calloc(n, sizeof(double));
  fjac = (double**) calloc(n, sizeof(double*));
  for (i=0; i<n; i++)
    fjac[i] = (double*) calloc(n, sizeof(double));


  f = dotF(n, x, fvec, args, vecFunc);
  test = 0.0;
  for (i=0;i<n;i++)
    if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
  if (test < 0.01*TOLF)
    {
      *check = 0;
      freeNewt(n, indx, g, p, xold, fvec, fjac);
      return;
    }
  for (sum=0.0,i=0;i<n;i++) sum+= x[i]*x[i];
  stpmax=STPMAX*fmax(sqrt(sum),(double)n);
        
  for (its=0; its<MAXITS; its++) /* begin iteration loop */
    {
      if (jacFunc == NULL) // use finite-difference Jacobian approximation
        fdjac(n,x,fvec,fjac,args,vecFunc);  // (re)initialise Jacobian in r
      else // compute jacobian from function
        jacFunc(n,x,fjac,args);

      for (i=0;i<n;i++) /* compute grad(f) for the line search */
        {
          for (sum=0.0,j=0;j<n;j++) sum += fjac[j][i]*fvec[j];
          g[i] = sum;
        }
      for (i=0;i<n;i++) xold[i]=x[i]; /* store x */
      fold=f;
      for (i=0;i<n;i++) p[i] = -fvec[i]; /* RHS for linear equations */
      ludcmp(fjac,n,indx,&d);
      lubksb(fjac,n,indx,p);
      lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fvec,args,vecFunc,dotF);
                        
      test = 0.0;
      for (i=0;i<n;i++)
        if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
      if (test < TOLF)
        {
          *check=0;
          freeNewt(n, indx, g, p, xold, fvec, fjac);
          return;
        }
      if (*check)
        {
          test=0.0;
          den=fmax(f,0.5*n);
          for (i=0;i<n;i++)
            {
              temp=fabs(g[i])*fmax(fabs(x[i]),1.0)/den;
              if (temp > test) test=temp;
            }
          *check=(test < TOLMIN ? 1 : 0);
          freeNewt(n, indx, g, p, xold, fvec, fjac);
          return;
        }
      test=0.0;
      for (i=0;i<n;i++)
        {
          temp=(fabs(x[i]-xold[i]))/fmax(fabs(x[i]),1.0);
          if (temp > test) test=temp;
        }
      if (test < TOLX)
        {
          freeNewt(n, indx, g, p, xold, fvec, fjac);
          return;
        }
    } /* end iteration loop */
  *check = 1;
  freeNewt(n, indx, g, p, xold, fvec, fjac);
  fprintf(stderr,"ERROR: MAXITS exceeded in newt\n");
}
/* -------------------------------------------------------------------------- */

void broyden(double *x, int n, int *check, void *args, 
             void (*vecFunc)(int, double *, double *, void *),
             void (*jacFunc)(int, double *, double **, void *))
{
  int i, j, k;
  int its, restart, sing, skip;
  double f, fold, maxStep, sum, den, temp, test, *c, *d, *fvcold;
  double *g, *p, **qt, **r, *s, *t, *w, *xold, *fvec;

  // allocate memory
  c      = (double*)  calloc(n, sizeof(double));
  d      = (double*)  calloc(n, sizeof(double));
  g      = (double*)  calloc(n, sizeof(double));
  p      = (double*)  calloc(n, sizeof(double));
  s      = (double*)  calloc(n, sizeof(double));
  t      = (double*)  calloc(n, sizeof(double));
  w      = (double*)  calloc(n, sizeof(double));
  xold   = (double*)  calloc(n, sizeof(double));
  fvec   = (double*)  calloc(n, sizeof(double));
  fvcold = (double*)  calloc(n, sizeof(double));
  qt     = (double**) calloc(n, sizeof(double*));
  r      = (double**) calloc(n, sizeof(double*));

  for (i=0; i<n; i++)
    {
      qt[i] = (double*) calloc(n, sizeof(double));
      r[i]  = (double*) calloc(n, sizeof(double));
    }

  f = dotF(n, x, fvec, args, vecFunc);
  test = 0.0;
  for (i=0; i<n; i++)
    if (fabs(fvec[i]) > test) test=fabs(fvec[i]);

  /* test whether the intial guess is a root using more stringent test
     If so, free memory and exit */
  if (test < 0.01*TOLF)
    {
      *check = 0;
      freeBroyden(n,c,d,g,p,s,t,w,xold,fvec,fvcold,qt,r);
      return;
    }

  // calculate maximum step for line searches
  for (sum=0.0, i=0; i<n; i++) sum += x[i]*x[i];
  maxStep = STPMAX*fmax(sqrt(sum), (double)(n));
  restart = 1;

  for (its=1; its < MAXITS; its++)
    {
      if (restart)
        {
          if (jacFunc == NULL) // use finite-difference Jacobian approximation
            fdjac(n,x,fvec,r,args,vecFunc);  // (re)initialise Jacobian in r
          else // compute jacobian from function
            jacFunc(n,x,r,args);
          qrdcmp(r,n,c,d,&sing);      // QR decomposition of Jacobian
          if (sing)
            {
              fprintf(stderr, "Error: singular Jacobian in broyden\n");
              exit(2);
            }

          // compute Q^T explicitly
          for (i=0; i<n; i++)
            {
              for (j=0; j<n; j++)
                qt[i][j] = 0.0;
              qt[i][i] = 1.0;    // set diagonal
            } // end loop over x

          for (k=0; k<n; k++)
            {
              if (c[k])
                {
                  for (j=0; j<n; j++)
                    {
                      sum = 0.0;
                      for (i=k; i<n; i++)
                        sum += r[i][k]*qt[i][j];

                      sum /= c[k];
                      for (i=k; i<n; i++)
                        qt[i][j] -= sum*r[i][k];
                    }
                }
            } // end loop over ?

          for (i=0; i<n; i++)
            {
              r[i][i] = d[i];
              for (j=0; j<i; j++)
                r[i][j] = 0.0;
            } // end loop over ?
        } // end restart
      else // carry out Broyden update
        {
          for (i=0; i<n; i++) s[i] = x[i] - xold[i]; // s = delta*x
          for (i=0; i<n; i++)                        // t = R*s
            {
              for (sum=0.0, j=i; j<n; j++) sum += r[i][j]*s[j];
              t[i] = sum;
            }
          skip = 1;

          // w = delta*F - B*s
          for (i=0; i<n; i++)
            {
              for (sum=0.0, j=0; j<n; j++) sum += qt[j][i]*t[j];
              w[i] = fvec[i] - fvcold[i] - sum;

              // Don't update with noisy components of w
              if (fabs(w[i]) >= EPS*(fabs(fvec[i]) + fabs(fvcold[i])))
                skip = 0;
              else
                w[i] = 0.0;
            }
          if (!skip)
            {
              for (i=0; i<n; i++)
                {
                  for (sum=0.0, j=0; j<n; j++) sum += qt[i][j]*w[j];
                  t[i] = sum;
                }
              for (den=0.0, i=0; i<n; i++) den += s[i]*s[i];
              for (i=0; i<n; i++) s[i] /= den; // store s/(s*s) in s
              qrupdate(r,qt,n,t,s);   // update R and Q^T
              for (i=0; i<n; i++)
                {
                  if (r[i][i] == 0.0)
                    {
                      fprintf(stderr, "Error: r singular in broyden\n");
                      exit(2);
                    }
                  d[i] = r[i][i]; // diagonal of R stored in d
                }
            } // end if skip
        } // end restart

      // r.h.s for linear equations -Q^T*F
      for (i=0; i<n; i++)
        {
          for (sum=0.0, j=0; j<n; j++) sum += qt[i][j]*fvec[j];
          p[i] = -sum;
        }

      // compute grad(f) = transpose(Q*R)*F for line search
      for (i=n-1; i>=0; i--)
        {
          for (sum=0.0, j=0; j<n; j++) sum -= r[j][i]*p[j];
          g[i] = sum;
        }

      // store x and F
      for (i=0; i<n; i++)
        {
          xold[i]   = x[i];
          fvcold[i] = fvec[i];
        }

      fold = f;
      rsolve(r,n,d,p); // solve linear equations

      /* lnsrch returns new x and f, also calculates fvec at the new x when
         it calls dotF */
      lnsrch(n,xold,fold,g,p,x,&f,maxStep,check,fvec,args,vecFunc,dotF);

      test = 0.0;
      for (i=0; i<n; i++)
        if (fabs(fvec[i]) > test)
          test = fabs(fvec[i]);

      if (test < TOLF)
        {
          *check = 0;
          freeBroyden(n,c,d,g,p,s,t,w,xold,fvec,fvcold,qt,r);
          return;
        }

      if (*check)
        {
          if (restart)
            {
              freeBroyden(n,c,d,g,p,s,t,w,xold,fvec,fvcold,qt,r);
              return;
            }
          else
            {
              test = 0.0;
              den = fmax(f, 0.5*n);
              for (i=0; i<n; i++)
                {
                  temp = fabs(g[i])*fmax(fabs(x[i]), 1.0)/den;
                  if (temp > test)
                    test = temp;
                }
              if (test < TOLMIN)
                {
                  freeBroyden(n,c,d,g,p,s,t,w,xold,fvec,fvcold,qt,r);
                  return;
                }
              else
                restart = 1;
            }
        }
      else
        {
          restart = 0;
          test = 0.0;
          for (i=0; i<n; i++)
            {
              temp = (fabs(x[i]-xold[i]))/fmax(fabs(x[i]), 1.0);
              if (temp > test)
                test = temp;
            }
          if (test < TOLX)
            {
              freeBroyden(n,c,d,g,p,s,t,w,xold,fvec,fvcold,qt,r);
              return;
            }
        }
    } // end iterations
  *check = 1;
  freeBroyden(n,c,d,g,p,s,t,w,xold,fvec,fvcold,qt,r);
  fprintf(stderr, "Warning: max iterations exceeded in broyden\n");
}
/* -------------------------------------------------------------------------- */

void lnsrch(int n, double *xold, double fold, double *g, double *p,
            double *x, double *f, double stpmax, int *check, double *fvec,
            void *args, void (*vecFunc)(int, double *, double *, void *), 
            double (*func)(int, double *, double *, void *,
                           void(*)(int, double *, double *, void *)))
{
  int i;
  double a, alam, alam2, alamin, b, disc, f2, rhs1, rhs2;
  double slope, sum, temp, test, tmplam;
        
  *check = 0;
  for (sum = 0.0, i = 0; i < n; i++)
    sum += p[i]*p[i];
  sum = sqrt(sum);
  // scale if step used is too large
  if (sum > stpmax)
    for (i = 0; i < n; i++)
      p[i] *= stpmax/sum;
        
  for (slope = 0.0, i = 0; i < n; i++)
    slope += g[i]*p[i];
                
  if (slope >= 0.0) {
    fprintf(stderr, "\nRoundoff problem in lineSearch\n");
    exit(2);
  }
        
  test = 0.0;
  for (i = 0; i < n; i++) {
    temp = fabs(p[i])/fmax(fabs(xold[i]),1.0);
    if (temp > test)
      test = temp;
  }
        
  alamin = TOLX/test;
  /*TODO: alam2 & f2 previously not initialized */
  /*The following values has been assumed by MD */
  /*but needs to be double checked!*/

  alam   = alam2 = 1.0; 
  f2  	 = *f;
  for (;;) {
    for (i = 0; i < n; i++)
      x[i] = xold[i] + alam*p[i];
    *f = (*func)(n, x, fvec, args, vecFunc);
    if (alam < alamin) {
      for (i = 0; i < n; i++)
        x[i] = xold[i];
      *check = 1;
      return;
    } else if (*f <= fold+ALF*alam*slope)
      return;
    else {
      if (alam == 1.0)
        tmplam = -slope/(2.0*(*f-fold-slope));
      else {
        rhs1 = *f - fold - alam*slope;
        rhs2 = f2 - fold - alam2*slope;
        a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
        b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
        if (a == 0.0)
          tmplam = -slope/(2.0*b);
        else {
          disc = b*b-3.0*a*slope;
          if (disc < 0.0)
            tmplam = 0.5*alam;
          else if (b <= 0.0)
            tmplam = (-b + sqrt(disc))/(3.0*a);
          else
            tmplam = -slope/(b+sqrt(disc));
        }
        if (tmplam > 0.5*alam)
          tmplam = -0.5*alam;
      }
    }
    alam2 = alam;
    f2 = *f;
    alam = fmax(tmplam,0.1*alam);
  } // end for
} // end lineSearch
/* -------------------------------------------------------------------------- */

void fdjac(int n, double *x, double *fvec, double **df, void *args,
           void (*vecFunc)(int, double *, double *, void *))
{
  int i, j;
  double h, temp, *f;

  f = (double*) calloc(n, sizeof(double));
  for (j=0; j<n; j++)
    {
      temp = x[j];
      h = EPS*fabs(temp);
      if (h == 0.0) h = EPS;
      x[j] = temp + h;
      h = x[j] - temp;
      (*vecFunc)(n, x, f, args);
      x[j] = temp;
      for (i=0; i<n; i++)
        df[i][j] = (f[i] - fvec[i])/h;
    }
  free(f);
}
/* -------------------------------------------------------------------------- */

double dotF(int n, double *x, double *fvec, void *args,
            void(*vecFunc)(int, double *, double *, void *))
{
  int i;
  double sum;

  (*vecFunc)(n,x,fvec,args);
  for (sum=0.0, i=0; i<n; i++) sum += fvec[i]*fvec[i];
  return 0.5*sum;
}
/* -------------------------------------------------------------------------- */

void freeNewt(int n, int *a, double *b, double *c, double *d, double *e,
              double **aa)
{
  int ii;
  free(a); free(b); free(c); free(d); free(e);
  for (ii=0; ii<n; ii++) free(aa[ii]); free(aa);
}
/* -------------------------------------------------------------------------- */

void freeBroyden(int n, double *a, double *b, double *c, double *d, double *e,
                 double *f, double *g, double *h, double *i, double *j,
                 double **aa, double **bb)
{
  int ii;
  free(a); free(b); free(c); free(d); free(e); free(f); free(g); free(h);
  free(i); free(j);
  for (ii=0; ii<n; ii++) { free(aa[ii]); free(bb[ii]); }
  free(aa); free(bb);
}
/* -------------------------------------------------------------------------- */

void mrqmin(double *x, double *y, double *sig, int ndat, double *a, int *ia,
            int ma, double **covar, double **alpha, double *chisq, void *args,
            void (*funcs)(double *, double *, double **, int, void *))
{
  int j,k,l,iter,done=0;
  static int NDONE = 4;
  static int mfit;
  static double alamda, ochisq,*atry,*beta,*da,**oneda,**temp;

  atry = (double*) calloc(ma, sizeof(double));
  beta = (double*) calloc(ma, sizeof(double));
  da   = (double*) calloc(ma, sizeof(double));
  
  mfit = 0;
  for (j=0; j<ma; j++) if (ia[j]) mfit++;
  oneda = (double**) calloc(mfit, sizeof(double*));
  temp  = (double**) calloc(mfit, sizeof(double*));
  for (j=0; j<mfit; j++)
    {
      oneda[j] = (double*) calloc(1, sizeof(double));
      temp[j]  = (double*) calloc(mfit, sizeof(double));
    }
  alamda = 0.001;
  mrqcof(x,y,sig,ndat,a,ia,ma,alpha,beta,chisq,args,funcs);
  for (j=0; j<ma; j++) atry[j] = a[j];
  ochisq = (*chisq);

  for (iter=0; iter<MAXITS; iter++)
    {
      if (done == NDONE) alamda = 0.0; /* last pass, use zero alamda */
      for (j=0; j<mfit; j++)
        {
          for (k=0; k<mfit; k++) covar[j][k] = alpha[j][k];
          covar[j][j] = alpha[j][j]*(1.0+alamda);
          for (k=0; k<mfit; k++) temp[j][k]  = covar[j][k];
          oneda[j][0] = beta[j];
        }
      gaussj(temp,mfit,oneda,1);
      for (j=0; j<mfit; j++)
        {
          for (k=0; k<mfit; k++) covar[j][k]=temp[j][k];
          da[j] = oneda[j][0];
        }
      if (done == NDONE)
        {
          covsrt(covar, ma, ia, mfit);
          covsrt(alpha, ma, ia, mfit);
          free(atry); free(beta); free(da);
          for (j=0; j<mfit; j++) { free(oneda[j]); free(temp[j]); }
          free(oneda); free(temp);
          return;
        }
      for (j=0, l=0; l<ma; l++)
        if (ia[l]) atry[l] = a[l]+da[j++];
      mrqcof(x,y,sig,ndat,atry,ia,ma,covar,da,chisq,args,funcs);
      if (fabs(*chisq-ochisq) < fmax(1e-9, 1e-9*(*chisq))) done++;
      if (*chisq < ochisq)
        {
          alamda *= 0.1;
          ochisq = *chisq;
          for (j=0; j<mfit; j++)
            {
              for (k=0; k<mfit; k++) alpha[j][k] = covar[j][k];
              beta[j] = da[j];
            }
          for (l=0; l<ma; l++) a[l] = atry[l];
        }
      else
        {
          alamda *= 10.0;
          *chisq = ochisq;
        }
    } /* end iterations */

  free(atry); free(beta); free(da);
  for (j=0; j<mfit; j++) { free(oneda[j]); free(temp[j]); }
  free(oneda); free(temp);
}
/* -------------------------------------------------------------------------- */

void mrqcof(double *x, double *y, double *sig, int ndat, double *a, int *ia, 
            int ma, double **alpha, double *beta, double *chisq, void *args,
            void (*funcs)(double *, double *, double **, int, void *))
{
  int i,j,k,l,m,mfit=0;
  double *ymod,wt,sig2i,dy,**dyda;
  ymod = (double*) calloc(ma, sizeof(double));
  dyda = (double**) calloc(ma, sizeof(double));
  for (j=0; j<ma; j++)
    dyda[j] = (double*) calloc(ma, sizeof(double));
  for (j=0; j<ma; j++) if (ia[j]) mfit++;
  for (j=0; j<mfit; j++)
    {
      for (k=0; k<=j; k++) alpha[j][k] = 0.0;
      beta[j] = 0.0;
    }
  *chisq = 0.0;
  (*funcs)(a,ymod,dyda,ma,args);
  for (i=0; i<ndat; i++)
    {
      sig2i = 1.0/(sig[i]*sig[i]);
      dy = y[i] - ymod[i];
      for (j=0, l=0; l<ma; l++)
        {
          if (ia[l])
            {
              wt = dyda[i][l]*sig2i;
              for (k=0, m=0; m<l+1; m++)
                if (ia[m]) alpha[j][k++] += wt*dyda[i][m];
              beta[j++] += dy*wt;
            }
        }
      *chisq += dy*dy*sig2i;
    }
  for (j=1; j<mfit; j++)
    for (k=0; k<j; k++) alpha[k][j] = alpha[j][k];
  free(ymod);
  for(j=0; j<ma; j++) { free(dyda[j]); }
  free(dyda);
}
/* -------------------------------------------------------------------------- */

void covsrt(double **covar, int ma, int *ia, int mfit)
{
  int i,j,k;
  double temp;
  for (i=mfit; i<ma; i++)
    for (j=0; j<i+1; j++) covar[i][j] = covar[j][i] = 0.0;
  k=mfit-1;
  for (j=ma-1; j>=0; j--)
    {
      if (ia[j])
        {
          for (i=0; i<ma; i++) SWAP(covar[i][k], covar[i][j]);
          for (i=0; i<ma; i++) SWAP(covar[k][i], covar[j][i]);
          k--;
        }
    }
}
/* -------------------------------------------------------------------------- */

/* functions for QR decomposition */
void qrdcmp(double **a, int n, double *c, double *d, int *sing)
{
  int i, j, k;
  double scale, sigma, sum, tau;
  
  *sing = 0;
  for (k=0; k<n-1; k++)
    {
      scale = 0.0;
      for (i=k; i<n; i++) scale = fmax(scale, fabs(a[i][k]));
      if (scale == 0.0)
        {
          *sing = 1;
          c[k] = d[k] = 0.0;
        }
      else
        {
          for (i=k; i<n; i++) a[i][k] /= scale;
          for (sum=0.0, i=k; i<n; i++) sum += a[i][k]*a[i][k];
          sigma = SIGN(sqrt(sum), a[k][k]);
          a[k][k] += sigma;
          c[k] = sigma*a[k][k];
          d[k] = -scale*sigma;
          for (j=k+1; j<n; j++)
            {
              for (sum=0.0, i=k; i<n; i++) sum += a[i][k]*a[i][j];
              tau = sum/c[k];
              for (i=k; i<n; i++) a[i][j] -= tau*a[i][k];
            }
        }
    }
  d[n-1] = a[n-1][n-1];
  if (d[n-1] == 0.0) *sing = 1;
}
/* -------------------------------------------------------------------------- */

void qrsolve(double **a, int n, double *c, double *d, double *b)
{
  int i, j;
  double sum, tau;

  // form (Q^trans, b)
  for (j=0; j<n; j++)
    {
      for (sum=0.0, i=j; i<n; i++) sum += a[i][j]*b[i];
      tau = sum/c[j];
      for (i=j; i<n; i++) b[i] -= tau*a[i][j];
    }
  rsolve(a,n,d,b); // solve (R, x) = (Q^trans, b)
}
/* -------------------------------------------------------------------------- */

void rsolve(double **a, int n, double *d, double *b)
{
  int i, j;
  double sum;

  b[n-1] /= d[n-1];
  for (i=n-2; i>=0; i--)
    {
      for (sum=0.0, j=i+1; j<n; j++) sum += a[i][j]*b[j];
      b[i]=(b[i]-sum)/d[i];
    }
}
/* -------------------------------------------------------------------------- */

void qrupdate(double **r, double **qt, int n, double *u, double *v)
{
  int i, j, k;

  // find largest k such that u[k] != 0.0
  for (k=n-1; k>=0; k--)
    {
      if (u[k]) break;
    }

  if (k<0) k = 0;

  // transform R + (u X v) to upper Hessenberg
  for (i=k-1; i>=0; i--)
    {
      rotate(r,qt,n,i,u[i], -u[i+1]);
      if (u[i] == 0.0) u[i] = fabs(u[i+1]);
      else if (fabs(u[i]) > fabs(u[i+1]))
        u[i] = fabs(u[i])*sqrt(1.0 + u[i+1]*u[i+1]/(u[i]*u[i]));
      else u[i] = fabs(u[i+1])*sqrt(1.0 + u[i]*u[i]/(u[i+1]*u[i+1]));
    }
  for (j=0; j<n; j++) r[0][j] += u[0]*v[j];

  // transform upper Hessenberg matrix to upper triangular
  for (i=0; i<k; i++)
    rotate(r,qt,n,i,r[i][i], -r[i+1][i]);
}
/* -------------------------------------------------------------------------- */

void rotate(double **r, double **qt, int n, int i, double a, double b)
{
  int j;
  double c, fact, s, w, y;

  // avoid unnecessary overflow or underflow
  if (a == 0.0)
    {
      c = 0.0;
      s = (b >= 0.0 ? 1.0 : -1.0);
    }
  else if (fabs(a) > fabs(b))
    {
      fact = b/a;
      c = SIGN(1.0/sqrt(1.0+(fact*fact)), a);
      s = fact*c;
    }
  else
    {
      fact = a/b;
      s = SIGN(1.0/sqrt(1.0+(fact*fact)), b);
      c = fact*s;
    }

  // pre-multiply r by Jacobi rotation
  for (j=i; j<n; j++)
    {
      y = r[i][j];
      w = r[i+1][j];
      r[i][j]   = c*y - s*w;
      r[i+1][j] = s*y + c*w;
    }

  // pre-multiply r by Jacobi rotation
  for (j=i; j<n; j++)
    {
      y = qt[i][j];
      w = qt[i+1][j];
      qt[i][j]   = c*y - s*w;
      qt[i+1][j] = s*y + c*w;
    }
}
/* -------------------------------------------------------------------------- */

void gaussj(double **a, int n, double **b, int m)
{
  int *indxc, *indxr, *ipiv;
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv,temp;

  indxc = (int*) calloc(n, sizeof(int));
  indxr = (int*) calloc(n, sizeof(int));
  ipiv  = (int*) calloc(n, sizeof(int));
  for (j=0; j<n; j++) ipiv[j] = 0;
  for (i=0; i<n; i++)
    {
      big = 0.0;
      for (j=0; j<n; j++)
        if (ipiv[j] != 1)
          for (k=0; k<n; k++)
            {
              if (ipiv[k] == 0)
                {
                  if (fabs(a[j][k]) >= big)
                    {
                      big = fabs(a[j][k]);
                      irow = j;
                      icol = k;
                    }
                }
            }
      ++(ipiv[icol]);
      if (irow != icol)
        {
          for (l=0; l<n; l++) SWAP(a[irow][l],a[icol][l]);
          for (l=0; l<m; l++) SWAP(b[irow][l],b[icol][l]);
        }
      indxr[i] = irow;
      indxc[i] = icol;
      if (a[icol][icol] == 0.0)
        {
          fprintf(stderr, "gaussj: singular matrix\n");
          exit(2);
        }
      pivinv = 1.0/a[icol][icol];
      a[icol][icol] = 1.0;
      for (l=0; l<n; l++) a[icol][l] *= pivinv;
      for (l=0; l<m; l++) b[icol][l] *= pivinv;
      for (ll=0; ll<n; ll++)
        if (ll != icol)
          {
            dum = a[ll][icol];
            a[ll][icol] = 0.0;
            for (l=0; l<n; l++) a[ll][l] -= a[icol][l]*dum;
            for (l=0; l<m; l++) b[ll][l] -= b[icol][l]*dum;
          }
    }
  for (l=n-1; l>=0; l--)
    {
      if (indxr[l] != indxc[l])
        for (k=0; k<n; k++)
          SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }
  free(ipiv); free(indxr); free(indxc);
}
/* -------------------------------------------------------------------------- */

void invertMatrix(int n, double **a, double **a_inv)
{
  int i,j;
  double d;
  int   *indx = (int*)    calloc(n, sizeof(double));
  double *col = (double*) calloc(n, sizeof(double));
  
  // make sure inverse is zero
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      a_inv[i][j] = 0.0;

  ludcmp(a,n,indx,&d);
  for (j=0; j<n; j++)
    {
      for (i=0; i<n; i++) col[i] = 0.0;
      col[j] = 1.0;
      lubksb(a,n,indx,col);
      for (i=0; i<n; i++) a_inv[i][j] = col[i];
    }
  free(indx);
  free(col);
}
/* -------------------------------------------------------------------------- */

void ludcmp(double **a, int n, int *indx, double *d)
{
  int i, imax, j, k;
  double big, dum, sum, temp;
  double *vv; // stores the implicit scaling of each row
        
  vv = (double*) calloc(n, sizeof(double));
  *d = 1.0; // no row interchanges yet
  // loop over rows to get implicit scaling information
  for (i = 0; i < n; i++)
    {
      big = 0.0;
      for (j = 0; j < n; j++)
        if ((temp =fabs(a[i][j])) > big)
          big = temp;
      if (big == 0.0)
        {
          fprintf(stderr,"Singular matrix in routine LU Decomposition\n");
          exit(1);
        }
      vv[i] = 1.0/big; // save scaling
    } // end implicit scaling loop
        
  // loop over columns for Crouts method
  for (j = 0; j < n; j++)
    {
      for (i = 0; i < j; i++)
        {
          sum = a[i][j];
          for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
          a[i][j] = sum;
        } // end i loop (1)
      big = 0.0;
      for (i = j; i < n; i++)
        {
          sum = a[i][j];
          for (k = 0; k < j; k++)
            sum -= a[i][k]*a[k][j];
          a[i][j] = sum;
          // check to see if figure of merit for the pivot is best so far
          if ( (dum=vv[i]*fabs(sum)) >= big)
            {
              big = dum;
              imax = i;
            } // end if
        } // end i loop (j)
      // check to see if rows need to be interchanged
      if (j != imax)
        { // if yes, do so
          for (k = 0; k < n; k++)
            {
              dum = a[imax][k];
              a[imax][k] = a[j][k];
              a[j][k] = dum;
            } // end k loop
          *d = -(*d); // change parity of d
          vv[imax] = vv[j]; // interchange scale factor
        } // end if
      indx[j] = imax;
      // if the pivot elemet is zero the matrix is singular
      if (a[j][j] == 0.0) a[j][j] = TINY;
      if (j != n-1) // divide by the pivot element
        {
          dum = 1.0/(a[j][j]);
          for (i = j+1; i < n; i++) a[i][j] *= dum;
        } // end if
    } // end j loop: go back for next column in the reduction
  free(vv);
}
/* -------------------------------------------------------------------------- */

void lubksb(double **a, int n, int *indx, double *b)
{
  int i, ii=-1, ip, j;
  double sum;
        
  for (i = 0; i < n; i++)
    {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii>=0)
        for (j = ii; j <= i-1; j++)
          sum -= a[i][j]*b[j];
      else if (sum)
        ii = i;
      b[i] = sum;
    } // end forward loop
  for (i = n-1; i >= 0; i--)
    {
      sum = b[i];
      for (j = i+1; j < n; j++)
        sum -= a[i][j]*b[j];
      b[i] = sum/a[i][i];
    } // end backward subsitution loop
}
/* -------------------------------------------------------------------------- */

void svdcmp(double **a, int m, int n, double *w, double **v)
{
  int flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

  rv1 = (double*) calloc(n, sizeof(double));
  g = scale = anorm = 0.0;
  for (i=0; i<n; i++) /* Householder reduction to bidiagonal form */
    {
      l = i+1;
      rv1[i] = scale*g;
      g = s = scale = 0.0;
      if (i < m) /* !!!!! < or <= ????? */
        {
          for (k=i; k<m; k++) scale += fabs(a[k][i]);
          if (scale)
            {
              for (k=i; k<m; k++)
                {
                  a[k][i] /= scale;
                  s += a[k][i]*a[k][i];
                }
              f = a[i][i];
              g = -SIGN(sqrt(s),f);
              h = f*g - s;
              a[i][i] = f - g;
              for (j=l; j<n; j++)
                {
                  for (s=0.0, k=i; k<m; k++) s += a[k][i]*a[k][j];
                  f = s/h;
                  for (k=i; k<m; k++) a[k][j] += f*a[k][i];
                }
              for (k=i; k<m; k++) a[k][i] *= scale;
            }
        }
      w[i] = scale*g;
      g = s = scale = 0.0;
      if (i < m && i != n-1) /* i !=n or i != n-1 ???? */
        {
          for (k=l; k<n; k++) scale += fabs(a[i][k]);
          if (scale)
            {
              for (k=l; k<n; k++)
                {
                  a[i][k] /= scale;
                  s += a[i][k]*a[i][k];
                }
              f = a[i][l];
              g = -SIGN(sqrt(s),f);
              h = f*g - s;
              a[i][l] = f - g;
              for (k=l; k<n; k++) rv1[k] = a[i][k]/h;
              for (j=l; j<m; j++)
                {
                  for (s=0.0, k=l; k<n; k++) s += a[j][k]*a[i][k];
                  for (k=l; k<n; k++) a[j][k] += s*rv1[k];
                }
              for (k=l; k<n; k++) a[i][k] *= scale;
            }
        }
      anorm = fmax(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
  for (i=n-1; i>=0; i--) /* Accumulation of right-hand transformations */
    {
      if (i < n-1)
        {
          if (g)
            {
              for (j=l; j<n; j++) /* double division to avoid underflow */
                v[j][i] = (a[i][j]/a[i][l])/g;
              for (j=l; j<n; j++)
                {
                  for (s=0.0, k=l; k<n; k++) s += a[i][k]*v[k][j];
                  for (k=l; k<n; k++) v[k][j] += s*v[k][i];
                }
            }
          for (j=l; j<n; j++) v[i][j] = v[j][i] = 0.0;
        }
      v[i][i] = 1.0;
      g = rv1[i];
      l = i;
    }
  for (i=IMIN(m,n)-1; i>=0; i--) /* Accumulation of left-hand transformations */
    {
      l = i+1;
      g = w[i];
      for (j=l; j<n; j++) a[i][j] = 0.0;
      if (g)
        {
          g = 1.0/g;
          for (j=l; j<n; j++)
            {
              for (s=0.0, k=l; k<m; k++) s += a[k][i]*a[k][j];
              f = (s/a[i][i])*g;
              for (k=i; k<m; k++) a[k][j] += f*a[k][i];
            }
          for (j=i; j<m; j++) a[j][i] *= g;
        }
      else
        for (j=i; j<m; j++) a[j][i] = 0.0;
      ++a[i][i];
    }
  for (k=n-1; k>=0; k--) /* diagonalization of the bidiagonal form */
    {
      for (its=1; its<=30; its++)
        {
          flag = 1;
          for (l=k; l>=0; l--) /* test for splitting */
            {                  /*  note that rv1[0] is always zero */
              nm = l-1; 
              if ((double)(fabs(rv1[l])+anorm) == anorm)
                {
                  flag = 0;
                  break;
                }
              if ((double)(fabs(w[nm])+anorm) == anorm) break;
            }
          if (flag)
            {
              c = 0.0;   /* cancellation of rv1[l], if l > 1 */
              s = 1.0;
              for (i=l; i<=k; i++)
                {
                  f = s*rv1[i];
                  rv1[i] = c*rv1[i];
                  if ((double)(fabs(f)+anorm) == anorm) break;
                  g = w[i];
                  h = pythag(f,g);
                  w[i] = h;
                  h = 1.0/h;
                  c =  g*h;
                  s = -f*h;
                  for (j=0; j<m; j++)
                    {
                      y = a[j][nm];
                      z = a[j][i];
                      a[j][nm] = y*c + z*s;
                      a[j][i]  = z*c - y*s;
                    }
                }
            }
          z = w[k];
          if (l == k) /* convergence */
            {
              if (z < 0.0) /* singular value is non-negative */
                {
                  w[k] = -z;
                  for (j=0; j<n; j++) v[j][k] = -v[j][k];
                }
              break;
            }
          if (its == 30)
            {
              fprintf(stderr, "no convergence in 30 svdcmp iterations\n");
              exit(2);
            }
          x = w[l];
          nm = k-1;
          y = w[nm];
          g = rv1[nm];
          h = rv1[k];
          f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
          g = pythag(f,1.0);
          f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
          c = s = 1.0; /* next QR transformation */
          for (j=l; j<=nm; j++) /* !!!! < or <=  ???? */
            {
              i = j+1;
              g = rv1[i];
              y = w[i];
              h = s*g;
              g = c*g;
              z = pythag(f,h);
              rv1[j] = z;
              c = f/z;
              s = h/z;
              f = x*c + g*s;
              h = y*s;
              y *= c;
              for (jj=0; jj<n; jj++)
                {
                  x = v[jj][j];
                  z = v[jj][i];
                  v[jj][j] = x*c + z*s;
                  v[jj][i] = z*c - x*s;
                }
              z = pythag(f,h);
              w[j] = z; /* Rotation can be arbitrary if z = 0 */
              if (z)
                {
                  z = 1.0/z;
                  c = f*z;
                  s = h*z;
                }
              f = c*g + s*y;
              x = c*y - s*g;
              for (jj=0; jj<m; jj++)
                {
                  y = a[jj][j];
                  z = a[jj][i];
                  a[jj][j] = y*c + z*s;
                  a[jj][i] = z*c - y*s;
                }
            }
          rv1[l] = 0.0;
          rv1[k] = f;
          w[k]   = x;
        }
    }
  free(rv1);
}
/* -------------------------------------------------------------------------- */

void svbksb(double **u, double *w, double **v, int m, int n, double *b,
            double *x)
{
  int jj,j,i;
  double s, *tmp;

  tmp = (double*) calloc(n, sizeof(double));
  for (j=0; j<n; j++) /* calculate U^T*B */
    {
      s = 0.0;
      if (w[j]) /* non-zero result only if w_j is non-zero */
        {
          for (i=0; i<m; i++) s += u[i][j]*b[i];
          s /= w[j]; /* this is the divide by w_j */
        }
      tmp[j] = s;
    }
  for (j=0; j<n; j++) /* matrix multiply by V to get answer */
    {
      s = 0.0;
      for (jj=0; jj<n; jj++) s += v[j][jj]*tmp[jj];
      x[j] = s;
    }
  free(tmp);
}
/* -------------------------------------------------------------------------- */

double pythag(double a, double b)
{
  double absa, absb;
  absa = fabs(a);
  absb = fabs(b);
  if (absa > absb) return absa*sqrt(1.0+absb*absb/(absa*absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+absa*absa/(absb*absb)));
}
/* -------------------------------------------------------------------------- */

void lfit(double *x, double *y, double *sig, int ndat, double *a, int *ia,
          int ma, double **covar, double *chisq,
          void (*funcs)(double, double *, int))
{
  int i,j,k,l,m,mfit=0;
  double ym,wt,sum,sig2i,**beta,*afunc;

  beta = (double**) calloc(ma, sizeof(double*));
  for (i=0; i<ma; i++) beta[i] = (double*) calloc(1, sizeof(double));
  afunc = (double*) calloc(ma, sizeof(double));
  for (j=0; j<ma; j++)
    if (ia[j]) mfit++;
  if (mfit == 0)
    {
      fprintf(stderr, "lfit: no parameters to be fitted\n");
      exit(2);
    }
  for (j=0; j<mfit; j++) /* initialise the symmetric matrix */
    {
      for (k=0; k<mfit; k++) covar[j][k] = 0.0;
      beta[j][0] = 0.0;
    }
  for (i=0; i<ndat; i++) /* loop over data to accumulate coefficients of */
    {                    /*   the normal matrix */
      (*funcs)(x[i],afunc,ma);
      ym = y[i];
      if (mfit < ma) /* subtract off dependencies on known pieces of the */
        {            /*   fitting function */
          for (j=0; j<ma; j++)
            if (!ia[j]) ym -= a[j]*afunc[j];
        }
      sig2i = 1.0/(sig[i]*sig[i]);
      for (j=-1, l=0; l<ma; l++)
        {
          if (ia[l])
            {
              wt = afunc[l]*sig2i;
              for (j++, k=-1, m=0; m<=l; m++)
                if (ia[m]) covar[j][++k] += wt*afunc[m];
              beta[j][0] += ym*wt;
            }
        }
    }
  for (j=1; j<mfit; j++) /* fill in above diagonal from symmetry */
    for (k=0; k<j; k++)
      covar[k][j] = covar[j][k];
  gaussj(covar,mfit,beta,1); /* matrix solution */
  for (j=-1, l=0; l<ma; l++)
    if (ia[l]) a[l] = beta[++j][0]; /* partition solution to appropriate coeff*/
  *chisq = 0.0;
  for (i=0; i<ndat; i++)
    {
      (*funcs)(x[i],afunc,ma);
      for (sum=0.0, j=0; j<ma; j++) sum += a[j]*afunc[j];
      *chisq += (y[i]-sum)*(y[i]-sum)/(sig[i]*sig[i]);
    }
/*   covsrt(covar,ma,ia,mfit); */
  free(afunc);
  for (i=0; i<ma; i++) free(beta[i]);
  free(beta);
}
/* -------------------------------------------------------------------------- */

void svdfit(double *x, double *y, double *sig, int ndata, double *a, int ma,
            double **u, double **v, double *w, double *chisq,
            void (*funcs)(double, double *, int))
{
  int i,j;
  double wmax,tmp,thresh,sum,*b,*afunc;

  b     = (double*) calloc(ndata, sizeof(double));
  afunc = (double*) calloc(ma,    sizeof(double));

  for (i=0; i<ndata; i++) /* accumulate coefficients of the fitting matrix */
    {
      (*funcs)(x[i],afunc,ma);
      tmp = 1.0/sig[i];
      for (j=0; j<ma; j++) u[i][j] = afunc[j]*tmp;
      b[i] = y[i]*tmp;
    }
  svdcmp(u,ndata,ma,w,v); /* singular value decomposition */
  wmax = 0.0;
  for (j=0; j<ma; j++)
    if (w[j] > wmax) wmax = w[j];
  thresh = TOLSVD*wmax;
  for (j=0; j<ma; j++)
    if (w[j] < thresh) w[j] = 0.0;
  svbksb(u,w,v,ndata,ma,b,a);
  *chisq = 0.0;
  for (i=0; i<ndata; i++)
    {
      (*funcs)(x[i],afunc,ma);
      for (sum=0.0, j=0; j<ma; j++) sum += a[j]*afunc[j];
      *chisq += (tmp=(y[i]-sum)/sig[i], tmp*tmp);
    }
  free(b);
  free(afunc);
}
/* -------------------------------------------------------------------------- */

void fitPoly(double x, double *p, int np)
{
  int j;

  p[0] = 1.0;
  for (j=1; j<np; j++) p[j]=p[j-1]*x;
}
/* -------------------------------------------------------------------------- */

void Linbcg(int n, double *b, double *x, int itol, double tol, int itmax,
            int *iter, double *err, void *args, 
            void (*atimes)(int, double* , double *, int, void *), 
            void (*asolve)(int, double* , double *, int, void *))
{
  double ak,akden,bk,bkden=1.0,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
  double eps_bcg = 1.0e-14;
  int j;
  double *p, *pp, *r, *rr, *z, *zz;
  p  = (double*) calloc(n, sizeof(double));
  pp = (double*) calloc(n, sizeof(double));
  r  = (double*) calloc(n, sizeof(double));
  rr = (double*) calloc(n, sizeof(double));
  z  = (double*) calloc(n, sizeof(double));
  zz = (double*) calloc(n, sizeof(double));

  *iter = 0;
  (*atimes)(n,x,r,0,args); // calculate the initial residual
  for (j=0; j<n; j++)
    {
      r[j]  = b[j]-r[j];
      rr[j] = r[j];
    }
  // (*atimes)(r,rr,0); /* uncomment this line to get the "minimum residual" */
  if (itol == 1)
    {
      bnrm = snrm(n, b, itol); // !!!! function
      (*asolve)(n,r,z,0,args);
    }
  else if (itol == 2)
    {
      (*asolve)(n,b,z,0,args);
      bnrm = snrm(n, z, itol);
      (*asolve)(n,r,z,0,args);
    }
  else if (itol == 3 || itol == 4)
    {
      (*asolve)(n,b,z,0,args);
      bnrm = snrm(n, z, itol);
      (*asolve)(n,r,z,0,args);
      znrm = snrm(n, z, itol);
    }
  else
    {
      fprintf(stderr, "Error: illegal itol in linbcg\n");
      exit(2);
    }

  while (*iter < itmax) /* main loop */
    {
      *iter = *iter + 1;
      (*asolve)(n,rr,zz,1,args);
      for (bknum=0.0,j=0; j<n; j++) bknum += z[j]*rr[j];
      /* calculate coefficient bk in the direction vectors p and pp */
      if (*iter == 1)
        {
          for (j=0; j<n; j++)
            {
              p[j]  = z[j];
              pp[j] = zz[j];
            }
        }
      else
        {
          bk = bknum/bkden;
          for (j=0; j<n; j++)
            {
              p[j]  = bk*p[j]+z[j];
              pp[j] = bk*pp[j]+zz[j];
            }
        }
      bkden = bknum;
      (*atimes)(n,p,z,0,args);
      for (akden=0.0,j=0; j<n; j++) akden += z[j]*pp[j];
      ak = bknum/akden;
      (*atimes)(n,pp,zz,1,args);
      for (j=0; j<n; j++)
        {
          x[j]  += ak*p[j];
          r[j]  -= ak*z[j];
          rr[j] -= ak*zz[j];
        }
      (*asolve)(n,r,z,0,args);
      if (itol == 1)
        *err = snrm(n, r, itol)/bnrm;
      else if (itol == 2)
        *err = snrm(n, z, itol)/bnrm;
      else if (itol == 3 || itol ==4)
        {
          zm1nrm = znrm;
          znrm   = snrm(n, z, itol);
          if (fabs(zm1nrm-znrm) > eps_bcg*znrm)
            {
              dxnrm = fabs(ak)*snrm(n, p, itol);
              *err   = znrm/fabs(zm1nrm-znrm)*dxnrm;
            }
          else
            {
              *err = znrm/bnrm;
              continue;
            }
          xnrm = snrm(n, x, itol);
          if (*err <= 0.5*xnrm) *err /= xnrm;
          else
            {
              *err = znrm/bnrm;
              continue;
            }
        }
      if (*err <= tol) break;
    } // end iteration

  free(p);
  free(pp);
  free(r);
  free(rr);
  free(z);
  free(zz);
}
/* -------------------------------------------------------------------------- */

double snrm(int n, double *sx, int itol)
{
  int i, isamax;
  double ans;
  if (itol <= 3)
    {
      ans = 0.0;
      for (i=0; i<n; i++) ans += sx[i]*sx[i];
      return sqrt(ans);
    }
  else
    {
      isamax = 0;
      for (i=0; i<n; i++)
        {
          if (fabs(sx[i]) > fabs(sx[isamax])) isamax = i;
        }
      return fabs(sx[isamax]);
    }
}
/* -------------------------------------------------------------------------- */

void LaplaceInterp_asolve(int n, double *b, double *x, int itrnsp, void *args)
{
  int i;
  for (i=0; i<n; i++) x[i] = b[i];
}
/* -------------------------------------------------------------------------- */

void LaplaceInterp_atimes_1D(int n, double *u, double *r, int itrnsp,void *args)
{
  int    i,j;
  double del;
  
  LapIntStruct *arg = (LapIntStruct*) args;
  int nx    = arg->n1;
  double *x = arg->x1;

  for (j=0; j<n; j++) r[j] = 0.0;
  for (j=0; j<nx; j++)
    {
      if (arg->mask[j]) // measured point
        {
          r[j] += u[j];
        }
      else if (j == 0) // left boundary
        {
          if (itrnsp)
            {
              r[j]   += u[j];
              del = -1.0/(x[j+1]-x[j])*u[j];
              r[j+1] +=   del*(x[j+2] - x[j]);
              r[j+2] += - del*(x[j+2] - x[j+1]);
            }
          else
            {
              r[j] = u[j] - 1.0/(x[j+1]-x[j])*( (x[j+2]-x[j])*u[j+1] -
                                                (x[j+2]-x[j+1])*u[j+2] );
            }
        }
      else if (j == nx-1) // right boundary
        {
          if (itrnsp)
            {
              r[j]   += u[j];
              del = -1.0/(x[j]-x[j-1])*u[j];
              r[j-1] +=   del*(x[j]   - x[j-2]);
              r[j-2] += - del*(x[j-1] - x[j-2]);
            }
          else
            {
              r[j] = u[j] - 1.0/(x[j]-x[j-1])*( (x[j]-x[j-2])*u[j-1] -
                                                (x[j-1]-x[j-2])*u[j-2] );
            }
        }
      else // interior point
        {
          if (itrnsp)
            {
              r[j]   += u[j];
              del = -1.0/(x[j+1]-x[j-1])*u[j];
              r[j+1] += del*(x[j]   - x[j-1]);
              r[j-1] += del*(x[j+1] - x[j]);
            }
          else
            {
              r[j] = u[j] - 1.0/(x[j+1]-x[j-1])*( (x[j]-x[j-1])*u[j+1] + 
                                                  (x[j+1]-x[j])*u[j-1] );
            }
        }
    }
}
/* -------------------------------------------------------------------------- */

void LaplaceInterp_atimes_2D(int n, double *u, double *r, int itrnsp,void *args)
{
  int    i, ic, ip, im, j, jc, jp, jm, k, kc, kp, km, f;
  double del, den, dx, dxm, dxp, dy, dym, dyp;
  
  LapIntStruct *arg = (LapIntStruct*) args;
  int nx    = arg->n1;
  int ny    = arg->n2;
  double *x = arg->x1;
  double *y = arg->x2;

  for (j=0; j<n; j++) r[j] = 0.0;
  for (k=0; k<ny; k++)
    {
      i = k*nx;
      if (k == 0)
        {
          km  = k;
          kc  = k+1;
          kp  = k+2;
        }
      else if (k == ny-1)
        {
          km  = k-2;
          kc  = k-1;
          kp  = k;
        }
      else
        {
          km  = k-1;
          kc  = k;
          kp  = k+1;
        }
      dym = y[kc] - y[km];
      dy  = y[kp] - y[km];
      dyp = y[kp] - y[kc];
      im  = km*nx;
      ic  = kc*nx;
      ip  = kp*nx;
      

      for (j=0; j<nx; j++)
        {
          if (j == 0)
            {
              jm = j;
              jc = j+1;
              jp = j+2;
            }
          else if (j == nx-1)
            {
              jm = j-2;
              jc = j-1;
              jp = j;
            }
          else
            {
              jm = j-1;
              jc = j;
              jp = j+1;
            }
          dxm = x[jc] - x[jm];
          dx  = x[jp] - x[jm];
          dxp = x[jp] - x[jc];

          if (arg->mask[i+j]) // measured point
            {
              r[i+j] += u[i+j];
            }
          else
            {
              if (itrnsp)
                {
                  r[im+jc] += -1/(dx*dxm)*u[i+j];
                  r[ic+jc] += 1/dx*(1/dxp+1/dxm)*u[i+j];
                  r[ip+jc] += -1/(dx*dxp)*u[i+j];
                  r[ic+jm] += -1/(dy*dym)*u[i+j];
                  r[ic+jc] += 1/dy*(1/dyp+1/dym)*u[i+j];
                  r[ic+jp] += -1/(dy*dyp)*u[i+j];
                }
              else
                r[i+j] =   1/dx*( (1/dxp+1/dxm)*u[ic+jc] - u[ic+jp]/dxp -
                                  u[ic+jm]/dxm)
                         + 1/dy*( (1/dyp+1/dym)*u[ic+jc] - u[ip+jc]/dyp -
                                  u[im+jc]/dym);
            }
        } // end loop over columns
    } // end loop over rows
}
/* -------------------------------------------------------------------------- */
