#include <stdlib.h>
#include <math.h>
#include "mducConst.h"
#include "transport.h"

double calcViscosity(double Mw, double eps_k, double sigma, double T)
{
  double k_b = 1.3806488e-23;
  double PI  = 3.141592653589793;
  double cst = 6.5511372043853852e-13; /* 5/16*sqrt(pi*k_b)/pi */

  double num, den;
  double T_red;     /* reduced temperature [-] */
  double omega;     /* collision integral  [-] */

  /* compute the reduced temperature */
  T_red = T/eps_k;

  /* compute the collision integral */
  omega = calcOmegaMu(T_red);
  
  return cst*pow(Mw*T,0.5)/(sigma*sigma*omega);
}
/* -------------------------------------------------------------------------- */

double calcViscCoeff(double Mw, double sigma)
{
  /* FIXME: check the exponent in the power function */
  return 2.6693e-6*pow(Mw,0.5)/(sigma*sigma);
}
/* -------------------------------------------------------------------------- */

double calcMixViscosity(int nSpecies, double T, double *Y, chemStruct *chem)
{
  int i,j;
  double Tr, omega, *mu05, *mui05,mixMu, mixDen, tmp, phi_ij;
  double sqrt_T = pow(T,0.5);

  /* ----- compute the species viscosities ----- */
  for (i=0; i<nSpecies; i++) {
    /* compute the reduced temperature */
    Tr = T/chem->eps_k[i];

    /* compute the collision integral */
    omega = calcOmegaMu(Tr);
    
    /* add to the overall viscosity */
    chem->mu05[i]  = sqrt(chem->coeffMu[i]*sqrt_T/omega);
    chem->mu05i[i] = 1.0 / chem->mu05[i];
  }

  /* ----- compute the mixture viscosity ----- */
  /* TODO: use chem->transport flag to ensure no invalid numbers */
  mixMu = 0.0;
  for (i=0; i<nSpecies; i++) {
    mixDen = 0.0;
    for (j=0; j<nSpecies; j++) {
      tmp     = 1.0 + chem->mu05[i]*chem->mu05i[j]*chem->WjOverWi[i][j];
      phi_ij  = chem->phiMixMu_a[i][j]*tmp*tmp;
      mixDen += Y[j]*chem->Wi[j]*phi_ij;
    }
    mixMu += chem->mu05[i]*chem->mu05[i]*Y[i]/(chem->W[i]*mixDen);
  }

  return mixMu;
}
/* -------------------------------------------------------------------------- */

double calcThermalCond()
{
  double lambda = 0.0;

  return lambda;
}
/* -------------------------------------------------------------------------- */

void fillBinDiffMatrix(double T, double P, double **Dij, chemStruct *chem)
{
  int i,j,n,p;
  double cst = 7.674784138632968e-36; /* 3/16*sqrt(2*pi*kb^3)/pi */
  double Tr, T3, delta;
  double m_ij, sig_ij, eps_ij, muSq_ij;
  double xi, alpha_n, mu_p, omega, tmp;

  T3 = pow(T,3.0);

  for (i=0; i<chem->nSpecies; i++) {
    for (j=0; j<chem->nSpecies; j++) {
      /* compute properties of species pairs */
      m_ij    = chem->W[i]*chem->W[j]/(chem->W[i]+chem->W[j]);
      eps_ij  = pow(chem->eps_k[i]*chem->eps_k[j],0.5);
      sig_ij  = 0.5*(chem->sigma[i]+chem->sigma[j]);
      muSq_ij = chem->muDipole[i]*chem->muDipole[j];

      /* correct if polarity of molecules is different */
      if ((chem->muDipole[i] > 0.0 && chem->muDipole[j] < 1.0e-6) ||
          (chem->muDipole[j] > 0.0 && chem->muDipole[i] < 1.0e-6) ) {
        xi = 1.0;

        if (chem->muDipole[i] > 0.0) { /* i is polar, j is non-polar */
          p = i;
          n = j;
        }
        else { /* j is polar, i is non-polar */
          p = j;
          n = i;
        }

        /* compute reduced polarizability and dipole moments */
        alpha_n = chem->alphaPolar[n]/pow(chem->sigma[n],3);
        mu_p    = chem->muDipole[p]*pow(chem->eps_k[p]*c_kb*
                                        pow(chem->sigma[p],3),-0.5);
        
        /* correction factor */
        xi += 0.25*alpha_n*mu_p*pow(chem->eps_k[p]/chem->eps_k[n], 0.5);

        /* apply correction */
        eps_ij *= xi*xi;
        sig_ij *= pow(xi,-1.0/6.0);
        muSq_ij = 0.0;
      }

      /* calculate the reduced quantities */
      Tr    = T/eps_ij;
      delta = 0.5*muSq_ij;

      /* compute the collision integral */
      omega = calcOmegaD(Tr);

      /* binary diffusion coefficent */
      Dij[i][j] = Dij[j][i] = cst*pow(T3/m_ij,0.5)/(P*sig_ij*sig_ij*omega);
    }
  }
}
/* -------------------------------------------------------------------------- */

double calcOmegaMu(double Tr)
{
  /* set the coefficients */
  static double m1 = 3.3530622607;
  static double m2 = 2.53272006;
  static double m3 = 2.9024238575;
  static double m4 = 0.11186138893;
  static double m5 = 0.8662326188;      /* = -0.1337673812 + 1.0 */
  static double m6 = 1.3913958626;
  static double m7 = 3.158490576;
  static double m8 = 0.18973411754;
  static double m9 = 0.00018682962894;
  
  double num, den;
  num = m1 + Tr*(m2 + Tr*(m3 + Tr*m4));
  den = m5 + Tr*(m6 + Tr*(m7 + Tr*(m8 + Tr*m9)));
  return num/den;
}
/* -------------------------------------------------------------------------- */

double calcOmegaD(double Tr)
{               
  /* set the coefficients */    
  static double m1 = 6.8728271691;
  static double m2 = 9.4122316321;
  static double m3 = 7.7442359037;
  static double m4 = 0.23424661229;
  static double m5 = 1.45337701568;     /* = 1.0 + 0.45337701568 */
  static double m6 = 5.2269794238;
  static double m7 = 9.7108519575;
  static double m8 = 0.46539437353;
  static double m9 = 0.00041908394781;
  
  double num, den;
  num = m1 + Tr*(m2 + Tr*(m3 + Tr*m4));
  den = m5 + Tr*(m6 + Tr*(m7 + Tr*(m8 + Tr*m9)));
  return num/den;
}        
/* -------------------------------------------------------------------------- */
