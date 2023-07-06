#ifndef MDUC_SOOT_HMOM_H
#define MDUC_SOOT_HMOM_H

/* ---------------------------------------------------------- */
/* Header file for the HMOM soot model for detailed chemistry */
/*    Author: Michael E. Mueller                              */
/*    Date: 08/30/12                                          */
/*    Version: 1.1                                            */
/* ---------------------------------------------------------- */

#include <math.h>

/* --------------------------------------------------------- */
/* Number of Soot Moments and Moment Orders                  */
/*   Needs to be either 4 (3+1) for first-order moments only */
/*     or 7 (6+1) for second-order moments                   */
/* --------------------------------------------------------- */
static       int    nSootMoments;
const static double MomOrderV[6] = {0.0, 1.0, 0.0, 2.0, 1.0, 0.0};
const static double MomOrderS[6] = {0.0, 0.0, 1.0, 0.0, 1.0, 2.0};
double sootMom[7];

double oxidfact[4];
/* -------------- */
/* Soot Constants */
/* -------------- */
#define SootDensity 1800.0    /* kg/m^3 */
//#define SootMolarMass 12.0e-3 /* kg/mol */ /* FIXME: 24.0e-3 in FlameMaster */
//#define SootMolarMass 24.0e-3 /* kg/mol */ 
//#define Avogadro 6.022e23     /* 1/mol */
#define SootMolarMass 24.0 /* kg/kmol */ 
#define Avogadro 6.022e26     /* 1/kmol */
#define Pi (4.0*atan(1.0))
//#define Rgas 8.314            /* J/mol/K */
#define Rgas 8.314e3            /* J/kmol/K */
/* Soot Fractal Dimension */
#define SootDf 1.8
/* Surface Density of Hydrogenated Sites */
#define SootChi 1.7e19        /* 1/m^2 */
/* Coefficients for Fit to Small Surface Area Change */
#define SootFitC (2.0/3.0)
#define SootFitE -0.2043
/* Free Molecular Size Exponents */
#define SootAv (1.0-(2.0/SootDf))
#define SootAs ((3.0/SootDf)-1.0)

/* Surface chemistry rate parameters    */
/*   Units are in J, m^3, mol, K, and s */
/*   Reference: G. Blanquart and H. Pitsch, A joint volume-surface-
     hydrogen multi-variate model for soot formation, in: H. Bockhorn,
     A. D'Anna, A. Sarofim, H. Wang (Eds.), Combustion Generated Fine
     Carbonaceous Particles, Karlsruhe University Press, 2009, 439-466 */
#define SootA1f 6.72e-2
#define Sootn1f 3.33
#define SootE1f 6.09e6
#define SootA1b 6.44e-4
#define Sootn1b 3.79
#define SootE1b 27.96e6
#define SootA2f 1.00e+5
#define Sootn2f 1.80
#define SootE2f 68.42e6
#define SootA2b 8.68e+1
#define Sootn2b 2.36
#define SootE2b 25.46e6
#define SootA3f 1.13e16 //1.13e13 
#define Sootn3f -0.06
#define SootE3f 476.05e6
#define SootA3b 4.17e10
#define Sootn3b 0.15
#define SootE3b 0.0
#define SootA4  2.52e6
#define Sootn4  1.10
#define SootE4  17.13e6
#define SootA5  2.20e9
#define Sootn5  0.0
#define SootE5  31.38e6
#define SootGamma6 0.13
enum HMOM_RXN
{
  HMOM_IRXN_PAH,          /* A4 + A4        -> DIMER */
  HMOM_IRXN_SG_OH_to_H2O, /* Soot-H + OH   <-> Soot-* + H2O */
  HMOM_IRXN_SG_H_to_H2,   /* Soot-H + H    <-> Soot-* + H2 */
  HMOM_IRXN_SG_H_abs,     /* Soot-H        <-> Soot-* + H */
  HMOM_IRXN_SG_C2H2_add,  /* Soot-* + C2H2  -> Soot-H */
  HMOM_IRXN_OX_O2,        /* Soot-* + O2    -> Soot-* + 2CO */
  HMOM_IRXN_OX_OH,        /* Soot-H + OH    -> Soot-H + CO */
  HMOM_NRXN
};

/* Globale variables used for nucleation and condensation with more than one PAH */

  //List of Species for Nucleation & Condensation
  int nPAH;
  double *PAHNbrC2;
  double *PAHNbrH;
  double *PAHStick;

  double DimerVol;
  double NuclVol;
  double NuclSurf;
  double LABSTOL;

/* function to set soot order */
void sootHMOMsetOrder(int order);
void sootHMOMsetTolerance(double tol);

/* ------------------------------ */
/* Function Prototypes: Interface */
/* ------------------------------ */
void sootHMOM(double Temp, double Viscosity, double Density, double MolarMass,
              double *c_species, double *moments, double *src, double *rxnrates,
              double *rhodot, double *srcP, int *srcActive,double ABSTOL);
/* The main interface function takes the temperature, a few local properties,a
   few local species concentrations, and the local moments and returns the 
   moment source terms and reaction rates for PAH dimerization, surface growth,
   and oxidation. */
/* Note: There is no initialization function that needs to be called
   before the main interface. */
/* Units (SI) for interface function:
     T: K
     Dynamic Viscosity: kg / m-s
     Gas Density: kg / m^3
     Mixture MolarMass: kg / mol
     C_i: mol / m^3
     M_xy: m^(x+2/3y) / m^3 (e.g. 1/m^3, m^3/m^3 for M00, M10)
     src_xy: m^(x+2/3y) / m^3 (divided into production and consumption)
     Reaction Rates: mol / m^3-s
     Continuity Source Term: kg / m^3-s
*/
/* Internally, for soot, the units are converted to counting in mols of
   C atoms.  However, this is completely shielded from the user by
   converting units on entry/exit with the main interface. */

void sootHMOM_GasPhase_Coeff(double Temp, double Viscosity, double Density, double MolarMass,
              double *c_species, double *ksg, double *kox, double *ko2, double *ProdRate);

/* ----------------------------------------- */
/* Function Prototypes: Moment Interpolation */
/* ----------------------------------------- */
double FracMom(double VolOrd, double SurfOrd, double *moments);
double FracMomDensity(double VolOrd, double SurfOrd, double *moments, double density);
double FracMomLarge(double VolOrd, double SurfOrd, double *moments);

/* ------------------------------------------------------- */
/* Function Prototypes: Dimerization and Surface Chemistry */
/* ------------------------------------------------------- */
//double Dimerization(double Temp, double Density, double *moments);
double Dimerization(double Temp, double *C_PAH, double *moments,int *srcActive);
void SurfaceChemistryRates(double Temp, double C_OH, double C_H2O, double C_H, 
                           double C_H2, double C_C2H2, double C_O2,
                           double *ksg, double *kox, double *ko2);
void ComputeProdRate(double *C_PAH, double Temp, double *ProdRate);

/* -------------------------------------- */
/* Function Prototypes: Soot Source Terms */
/* -------------------------------------- */
void Nucleation(double Temp, double DimerConc, double *moments, double *src);
void Condensation(double Temp, double DimerConc, double *moments, double *src);
void Coagulation(double Temp,double Viscosity,double Density,double MolarMass,
                 double *moments, double *src);
void SurfaceGrowth(double ksg, double *moments, double *src);
void Oxidation(double kox, double *moments, double *src);
void Fragmentation(double ko2, double *moments, double *src);

/* --------------------------------------------------------------------- */
/* Function Prototypes: Collision Coeffs for Nucleation and Condensation */
/* --------------------------------------------------------------------- */
double GetBetaNucl(double Temp);
double GetBetaCond(double Temp, double *moments);
double GetBetaDimer(double Temp, double i);
double GetC(double temp);
/* ----------------------------------------------- */
/* Function Prototypes: Free Molecular Coagulation */
/* ----------------------------------------------- */
double FMCoagSL(int i, double *moments);
double PsiSL(double x, double y, double a, double b, double *moments);
double FMCoagLL(int i, double *moments);
double PsiLL(double x, double y, double a, double b, double *moments);

/* ------------------------------------------ */
/* Function Prototypes: Continuum Coagulation */
/* ------------------------------------------ */
double CNCoagSL(int i, double *moments, double Temp, double Viscosity,
                double Density, double MolarMass);
double CNCoagLL(int i, double *moments, double Temp, double Viscosity,
                double Density, double MolarMass);

/* --------------------------------------------- */
/* Function Prototypes: Gas-Phase Reaction Rates */
/* --------------------------------------------- */
void GasPhaseRxnRates(double Temp, double *C_PAH, double C_OH, double C_H2O, 
                      double C_H, double C_H2, double C_C2H2, double C_O2, 
                      double *moments, double *rxnrates, int *srcActive);

/* ------------------------------------------- */
/* Function Prototypes: Continuity Source Term */
/* ------------------------------------------- */
double ContinuitySource(double *rxnrates);

#endif /* MDUC_SOOT_HMOM_H */
