#include <stdio.h>
#include <limits.h>
#include "polyInt.h"

#define bigReal 1.E38  // FLT_MAX, DBL_MAX if double above

typedef struct{gPoint min; gPoint max;} box;
typedef long long hp;
typedef struct{long x; long y;} ipoint;
typedef struct{long mn; long mx;} rng;
typedef struct{ipoint ip; rng rx; rng ry; short in;} vertex;

void bd(real *X, real y){ *X = (*X < y) ? *X : y; }
void bu(real *X, real y){ *X = (*X > y) ? *X : y; }
/* -------------------------------------------------------------------------- */

void range(gPoint *x, int c, box *B)
{
  while(c--) {
    bd(&B->min.x, x[c].x);
    bu(&B->max.x, x[c].x);
    bd(&B->min.y, x[c].y);
    bu(&B->max.y, x[c].y);
  }
}
/* -------------------------------------------------------------------------- */

void fit(gPoint *x, int cx, vertex *ix, box B, real sclx, real scly,
         real mid, int fudge)
{
  int c;
  c = cx;
  /* loop over all vertices */
  while(c--) {
    ix[c].ip.x = (long)((x[c].x - B.min.x)*sclx - mid)&~7|fudge|c&1;
    ix[c].ip.y = (long)((x[c].y - B.min.y)*scly - mid)&~7|fudge;
  }
  ix[0].ip.y += cx&1;
  ix[cx].ip.x = ix[0].ip.x;
  ix[cx].ip.y = ix[0].ip.y;
  
  /* find the min/max co-ordinates of the polygon */
  c = cx;
  while(c--) {
    ix[c].rx = (ix[c].ip.x < ix[c+1].ip.x) ?
      (rng){ix[c].ip.x,ix[c+1].ip.x} : (rng){ix[c+1].ip.x,ix[c].ip.x};
    ix[c].ry = (ix[c].ip.y < ix[c+1].ip.y) ?
      (rng){ix[c].ip.y,ix[c+1].ip.y} : (rng){ix[c+1].ip.y,ix[c].ip.y};
    ix[c].in=0;
  }
  ix[cx].rx.mn = ix[0].rx.mn;
  ix[cx].rx.mx = ix[0].rx.mx;
  ix[cx].ry.mn = ix[0].ry.mn;
  ix[cx].ry.mx = ix[0].ry.mx;
  ix[cx].in    = ix[0].in;
}
/* -------------------------------------------------------------------------- */

hp area(ipoint a, ipoint p, ipoint q) {
  return (hp)p.x*q.y - (hp)p.y*q.x + (hp)a.x*(p.y - q.y) + (hp)a.y*(q.x - p.x);
}
/* -------------------------------------------------------------------------- */

void cntrib(ipoint f, ipoint t, short w, hp *s)
{s[0] += (hp)w*(t.x-f.x)*(t.y+f.y)/2;}

/* determine whether there is overlap */
int ovl(rng p, rng q){return p.mn < q.mx && q.mn < p.mx;}
/* -------------------------------------------------------------------------- */

void cross(vertex *a, vertex *b, vertex *c, vertex *d,
           double a1, double a2, double a3, double a4, hp *s)
{
  real r1=a1/((real)a1+a2), r2 = a3/((real)a3+a4);
  cntrib((ipoint){a->ip.x + r1*(b->ip.x-a->ip.x),
        a->ip.y + r1*(b->ip.y-a->ip.y)}, b->ip, 1, s);
  cntrib(d->ip, (ipoint){c->ip.x + r2*(d->ip.x - c->ip.x), 
        c->ip.y + r2*(d->ip.y - c->ip.y)}, 1, s);
  ++a->in;
  --c->in;
}
/* -------------------------------------------------------------------------- */

void inness(vertex *P, int cP, vertex * Q, int cQ, hp *s)
{
  int j, c=cQ;
  int ss = 0;
  ipoint p = P[0].ip;
  int sgn;
  while(c--) {
    if(Q[c].rx.mn < p.x && p.x < Q[c].rx.mx)
      {
        sgn = 0 < area(p, Q[c].ip, Q[c+1].ip);
        ss += sgn != Q[c].ip.x < Q[c+1].ip.x ? 0 : (sgn?-1:1);
      }
  }

  for (j=0; j<cP; ++j)
    {
      if (ss)
        cntrib(P[j].ip, P[j+1].ip, ss, s);
      ss += P[j].in;
    }
}
/* -------------------------------------------------------------------------- */

real polyIntersect(int na, gPoint *a, int nb, gPoint *b)
{
  int j, k;
  vertex ipa[na+1], ipb[nb+1];
  double ascale;
  box B = {{bigReal, bigReal}, {-bigReal, -bigReal}};

  /* exit if either a or b represents a line */
  if(na < 3 || nb < 3) return 0;

  /* get the bounding box around the polygons */
  range(a, na, &B);
  range(b, nb, &B);

  const real gamut = 500000000., mid = gamut/2.;
  real rngx = B.max.x - B.min.x, sclx = gamut/rngx;
  real rngy = B.max.y - B.min.y, scly = gamut/rngy;
    
  fit(a, na, ipa, B, sclx, scly, mid, 0);
  fit(b, nb, ipb, B, sclx, scly, mid, 2);
  ascale = sclx*scly;

  int o;
  hp a1, a2, a3, a4;
  hp s = 0;
  for (j=0; j<na; ++j)
    for(k=0; k<nb; ++k)
      if (ovl(ipa[j].rx, ipb[k].rx) && ovl(ipa[j].ry, ipb[k].ry))
        {
          a1 = -area(ipa[j].ip,   ipb[k].ip, ipb[k+1].ip);
          a2 =  area(ipa[j+1].ip, ipb[k].ip, ipb[k+1].ip);
          o = a1<0;
          if (o == a2<0)
            {
              a3 =  area(ipb[k].ip,   ipa[j].ip, ipa[j+1].ip);
              a4 = -area(ipb[k+1].ip, ipa[j].ip, ipa[j+1].ip);
              if (a3<0 == a4<0)
                {
                  if (o)
                    cross(&ipa[j],&ipa[j+1],&ipb[k],&ipb[k+1],a1,a2,a3,a4,&s);
                  else
                    cross(&ipb[k],&ipb[k+1],&ipa[j],&ipa[j+1],a3,a4,a1,a2,&s);
                }
            }
        }
  inness(ipa, na, ipb, nb, &s);
  inness(ipb, nb, ipa, na, &s);
  return s/ascale;
}
/* -------------------------------------------------------------------------- */
