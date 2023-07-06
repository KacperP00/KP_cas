#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "engineUtils.h"

double cad_to_time(double CAD, engineStruct *engine)
{
  return (CAD - engine->CAD_0)/(6*engine->RPM);
}
/* -------------------------------------------------------------------------- */

double time_to_cad(double t, engineStruct *engine)
{
  return t*6*engine->RPM + engine->CAD_0;
}
/* -------------------------------------------------------------------------- */

double calcCylLength(double CAD, engineStruct *engine)
{
  double PI    = 3.141592653589793;
  double theta = CAD*PI/180.0;      /* convert crank angle to radians */
  double rod   = engine->rod;
  double r_cam = 0.5*engine->stroke;
  double tmp1  = r_cam*sin(theta);
  double tmp2  = pow(rod*rod - tmp1*tmp1,0.5);

  return rod + r_cam - r_cam*cos(theta) - tmp2;
}
/* -------------------------------------------------------------------------- */

double calcPistonVelocity(double CAD, engineStruct *engine)
{
  double PI    = 3.141592653589793;
  double theta = CAD*PI/180.0;          /* convert crank angle to radians */
  double omega = 2*PI*engine->RPM/60.0; /* angular velocity */
  double rod   = engine->rod;
  double r_cam = 0.5*engine->stroke;
  double tmp1  = r_cam*sin(theta);
  double tmp2  = pow(rod*rod - tmp1*tmp1,0.5);

  return -omega*(tmp1 + r_cam*tmp1*cos(theta)/tmp2);
}
/* -------------------------------------------------------------------------- */


double calcCylVol(double CAD, engineStruct *engine)
{
  return engine->V_TDC + engine->Acyl*calcCylLength(CAD, engine);
}
/* -------------------------------------------------------------------------- */

double calcValveArea(double lift, valveStruct *v)
{
  if (lift <= 0.0) 
    return 0.0;

  double PI     = 3.141592653589793;
  double D_m    = v->D_head - v->W;
  double A_open = v->D_port*v->D_port - v->D_stem*v->D_stem;

  if (lift < v->W/(sin(v->theta)*cos(v->theta)) )
    return PI*lift*cos(v->theta)*(v->D_head-2*v->W +0.5*lift*sin(2*v->theta));
  else if (lift < pow(0.25*A_open/D_m - v->W*v->W,0.5)+v->W*tan(v->theta) )
    return PI*D_m*pow( pow(lift - v->W*tan(v->theta), 2) + v->W*v->W, 0.5);
  else
    return 0.25*PI*A_open;
}
/* -------------------------------------------------------------------------- */

double calcValveSteadyMassFlow(double gamma, double rhoi, double Pu, double Pd,
                               double Af)
{
  double term_sqrt, Cd;
  double mdot;
  double dir = 1.0;
  double pd_pu = Pd/Pu;

  double fm, fg1, fg2;

  fg1 = 2.0/(gamma+1);
  fg2 = gamma/(gamma-1.0);

  if (pd_pu <= pow(fg1,fg2)) { /* choked */
    fm = pow(gamma,0.5)*pow(fg1, (gamma+1.0)/(2.0*(gamma-1.0)));
  }
  else { /* unchoked */
    term_sqrt = 2*fg2*(1.0-pow(pd_pu,1.0/fg2));
    term_sqrt = fmax(term_sqrt, 0.0);
    fm = pow(pd_pu,1.0/gamma)*pow(term_sqrt,0.5);
  }

  /* compute discharge coefficient based on simple model */
  Cd = 0.85 - 0.25*pow(pd_pu,2);

  return Cd*Af*pow(Pu*rhoi,0.5)*fm;
}
/* -------------------------------------------------------------------------- */

double calcEngine_dPdt(double t1, double t2, double P1, double gamma,
                       engineStruct *eng)
{
  double V1, V2;
  V1 = calcCylVol(time_to_cad(t1,eng),eng);
  V2 = calcCylVol(time_to_cad(t2,eng),eng);

  if (t2 > t1)
    return P1*(pow(V1/V2,gamma)-1)/(t2-t1);
  else
    return 0.0;
}
/* -------------------------------------------------------------------------- */

double calcWoschniQdot(double t, double Pc, double Tc, double gamma, 
                       engineStruct *engine)
{
  double PI = 3.141592653589793;
  double C1; 
  double V, L, A, Vref;   /* geometric properties of the engine*/
  double omega;           /* engine speed in [rad/s] */
  double S_p;             /* mean piston velocity */
  double u_t;             /* turbulent velocity */
  double P_mot;           /* pressure assuming adiabatic comp. (motoring) */
  double h_e;             /* effective heat transfer coefficient */
  woschniStruct *w;
  if (engine->htModel == 1)
    w = &engine->hohen;
  else
    w = &engine->woschni;

  /* calculate the current engine geometry */
  omega = 2*PI/60.0*engine->RPM;
  /* L     = calcCylLength(time_to_cad(t,engine), engine); */
  V     = calcCylVol(time_to_cad(t,engine), engine);
  L     = V/engine->Acyl;
  A     = PI*engine->bore*L + 2*engine->Acyl;
  Vref  = engine->V_BDC;

  /* compute the turbulent velocity */
  S_p   = 2*engine->stroke*omega/(2*PI);
  w->C1 = 2.28 + 0.308*PI*engine->bore*omega/S_p;
  P_mot = w->Pref*pow(Vref/V, gamma);
  u_t   = w->C1*S_p + w->C2/6*engine->V_SWP*w->Tref*(Pc - P_mot)/(w->Pref*Vref);
  if (engine->htModel == 1) {
    L = V;
    u_t += w->C3;
  }
  u_t   = fmax(u_t, 0.0);

  /* compute the effective heat transfer coefficient */
  h_e = w->alpha*pow(L,w->a)*pow(Pc/w->Pscale,w->b)*pow(Tc,w->c)*pow(u_t,w->d);

  return -h_e*A*(Tc - engine->T_wall);
}
/* -------------------------------------------------------------------------- */

void calcEGR(double phi, double psi, double EGR, double M_fuel, int nC, int nH,
             int nO, double *nb, double *nu)
{
  double y, eps, zeta;

  /* default molar N/O ratio */
  if (psi == -1)
    psi = 3.773;

  /* compute the necessary constants */
  y   = nH/nC;
  eps = 4.0/(4.0+y);

  /* account for any oxygenated fuel */
  zeta = 2.0/(2.0-eps*nO*(1.0-phi));
  phi *= zeta;
  psi *= (1.0-0.5*eps*nO)*zeta;

  if (phi > 1.0) { /* fuel rich */
    double aa,bb,cc,tmp,c0,c1;
    double c, K;
    
    /* compute equilibrium constant Eq. (4.5) in Heywood */
    /* TODO: assumed constant burned temperature */
    K = calcEGR_K(1740);

    /* solve quadratic Eq. (4.6) in Heywood for c */
    c = 0.0;
    aa = K-1;
    bb = -1.0*( K*(2.0*(phi-1.0)+eps*phi) + 2.0*(1.0-eps*phi) );
    cc = 2.0*K*eps*phi*(phi-1.0);
    tmp = bb*bb - 4*aa*cc;
    if (tmp < 0.0) {
      fprintf(stderr, "Error: rich mixture computation failed");
      exit(2);
    }
    c0 = (-bb - pow(tmp,0.5))/(2*aa);
    c1 = (-bb + pow(tmp,0.5))/(2*aa);
    if (c1 > 0.99)
      c = c0;
    else if (c1 < 0.99 && c0 > 0.99)
      c = c1;
    else {
      fprintf(stderr, "Error: no constant in right range");
      exit(2);
    }

    /* compute the composition of the burned mixture (in moles) */
    nb[0] = 0.0;                                   /* Fuel */
    nb[1] = 0.0;                                   /* O2 */
    nb[2] = psi;                                   /* N2 */
    nb[3] = eps*phi - c;                           /* CO2 */
    nb[4] = 2*(1.0-eps*phi) + c;                   /* H2O */
    nb[5] = c;                                     /* CO */
    nb[6] = 2*(phi-1.0) - c;                       /* H2 */
    
    /* compute the composition of the unburned mixture (in moles) */
    nu[0] = 4.0*(1.0-EGR)*(1+2.0*eps)*phi/M_fuel;  /* Fuel */
    nu[1] = 1 - EGR;                               /* O2 */
    nu[2] = psi;                                   /* N2 */
    nu[3] = nb[3]*EGR;                             /* CO2 */
    nu[4] = nb[4]*EGR;                             /* H2O */
    nu[5] = nb[5]*EGR;                             /* CO */
    nu[6] = nb[6]*EGR;                             /* H2 */
  }
  else { /* fuel lean/stoichiometric */
    /* compute the composition of the burned mixture (in moles) */
    nb[0] = 0.0;                                   /* Fuel */
    nb[1] = 1 - phi;                               /* O2 */
    nb[2] = psi;                                   /* N2 */
    nb[3] = eps*phi;                               /* CO2 */
    nb[4] = 2*(1.0-eps)*phi;                       /* H2O */
    nb[5] = 0.0;                                   /* CO */
    nb[6] = 0.0;                                   /* H2 */
    
    /* compute the composition of the unburned mixture (in moles) */
    nu[0] = 4.0*(1.0-EGR)*(1+2.0*eps)*phi/M_fuel;  /* Fuel */
    nu[1] = 1 - EGR*phi;                           /* O2 */
    nu[2] = psi;                                   /* N2 */
    nu[3] = nb[3]*EGR;                             /* CO2 */
    nu[4] = nb[4]*EGR;                             /* H2O */
    nu[5] = 0.0;                                   /* CO */
    nu[6] = 0.0;                                   /* H2 */
  }
}
/* -------------------------------------------------------------------------- */

double calcEGR_K(double Tb)
{
  double tmp;
  tmp = 2.743 - 1.0/Tb*(1.761e3 + 1.0/Tb*(1.611e6 - 0.2803e9/Tb));
  return exp(tmp);
}
/* -------------------------------------------------------------------------- */

double calcFlameSpeedHeywood(char *f,double phi,double xb,double P,double T)
{
  int form = 0;
  double SL0,phiM,Bm,Bphi,alpha,beta;
  double Tref = 298.0;
  double Pref = 101325;

  /* select constants for relevant fuel */
  if (!strcmp(f,"methanol")) {
    phiM = 1.11;
    Bm   =   36.9; /* cm/s */
    Bphi = -140.5; /* cm/s */
  }
  else if (!strcmp(f,"propane")) {
    phiM = 1.08;
    Bm   =   34.2; /* cm/s */
    Bphi = -138.7; /* cm/s */
  }
  else if (!strcmp(f,"iso-octane")) {
    phiM = 1.13;
    Bm   =   26.3; /* cm/s */
    Bphi =  -84.7; /* cm/s */
  }
  else if (!strcmp(f,"gasoline")) {
    phiM = 1.21;
    Bm   =   30.5; /* cm/s */
    Bphi =  -54.9; /* cm/s */
    form = 1;
  }
  else {
    fprintf(stderr, "Warning: Heywood flame speed not available for %s\n",f);
    fprintf(stderr, "         available fuels: methanol, propane, iso-octane, "
            "gasoline\n");
  }

  /* compute the P/T exponent factors */
  if (form == 1) {
    alpha =  2.4   - 0.271*pow(phi,3.51);
    beta  = -0.357 + 0.140*pow(phi,2.77);
  }
  else {
    alpha =  2.18 - 0.80*(phi - 1.0);
    beta  = -0.16 + 0.22*(phi - 1.0);
  }

  /* compute reference laminar flame speed for this equivalence ratio */
  SL0  = Bm + Bphi*(phi-phiM)*(phi-phiM);

  /* scale for diluent */
  SL0 *= (1.0 - 2.06*pow(xb,0.77));

  /* return flame speed scaled by pressure and temperature*/
  return SL0*pow(T/Tref,alpha)*pow(P/Pref,beta);
}
/* -------------------------------------------------------------------------- */
