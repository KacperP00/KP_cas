/* ----------------------------------------------------------- */
/* Source code for the HMOM soot model for detailed chemisitry */
/*    Author: Michael E. Mueller                               */
/*    Date: 08/30/12                                           */
/*    Version: 1.1                                             */
/* ----------------------------------------------------------- */

#include "sootHMOM.h"
#include <stdio.h>

void sootHMOMsetOrder(int order)
{
  nSootMoments = order;
}

void sootHMOMsetTolerance(double tol)
{
  LABSTOL = tol;
  //LABSTOL = 1.0e-15;
}
/* -------------------------------------------------------------------------- */

/* --------------------------------------------- */
/* Main Interface Function                       */
/*   Inputs: Array of soot moments in SI units   */
/*   Outputs: Soot source terms in SI units      */
/* --------------------------------------------- */
void sootHMOM(double Temp, double Viscosity, double Density, double MolarMass,
              double *c_species, double *moments, double *src, double *rxnrates,
              double *rhodot, double *srcP, int *srcActive,double ABSTOL)
{
  int i,j,k;
  double DimerConc;
  double ksg, kox, ko2, tmp;
  double C_PAH[nPAH];
  double V0 = 2.0 * NuclVol;
  //double dV = NuclVol;
  double dV = 1.0;
  double S0 = pow(V0,2.0/3.0);
  double dS = dV * S0 / V0;
  double N0,oxidmin,Mmin;
  double x;

  for (i=0;i<nSootMoments;i++)
      sootMom[i] = moments[i];

// clipping
  Mmin = LABSTOL*Density;
  if (moments[nSootMoments-1] > 0.0) x = fmax(0.0,(moments[nSootMoments-1] - Mmin))/moments[nSootMoments-1];
  else x = 0.0;
  oxidfact[nSootMoments-1]=-2.0 * pow(x,3.0) + 3.0*pow(x,2.0);    
  moments[nSootMoments-1] = fmax(moments[nSootMoments-1],Mmin);

  Mmin = moments[nSootMoments-1]+LABSTOL*Density;
  if (moments[0] > 0.0) x = fmax(0.0,(moments[0] - Mmin))/moments[0];
  else x = 0.0;
  oxidfact[0]=-2.0 * pow(x,3.0) + 3.0*pow(x,2.0);    
  moments[0] = fmax(moments[0],Mmin);

  Mmin = (moments[0]-moments[nSootMoments-1])*(V0+dV) + moments[nSootMoments-1]*V0;
  if (moments[1] > 0.0) x = fmax(0.0,(moments[1] - Mmin))/moments[1];
  else x = 0.0;
  oxidfact[1]=-2.0 * pow(x,3.0) + 3.0*pow(x,2.0);    
  moments[1] = fmax(moments[1],Mmin);

  Mmin = (moments[0]-moments[nSootMoments-1])*(S0+dS) + moments[nSootMoments-1]*S0;
  if (moments[2] > 0.0) x = fmax(0.0,(moments[2] - Mmin))/moments[2];
  else x = 0.0;
  oxidfact[2]=-2.0 * pow(x,3.0) + 3.0*pow(x,2.0);    
  moments[2] = fmax(moments[2],Mmin);

  if (nSootMoments > 4) 
     fprintf(stderr, "Error: Moment clipping only available for 1st order moments\n");

  for (i=0; i<nPAH; i++)
     C_PAH[i] = c_species[i];

  double C_OH   = c_species[nPAH];
  double C_H2O  = c_species[nPAH+1];
  double C_H    = c_species[nPAH+2];
  double C_H2   = c_species[nPAH+3];
  double C_C2H2 = c_species[nPAH+4];
  double C_O2   = c_species[nPAH+5];

  /* Zero the source terms */
  for (i=0; i<nSootMoments; i++)   src[i]  = 0.0;
  for (i=0; i<6*nSootMoments; i++) srcP[i] = 0.0;

  /* Compute the dimerization rate */
  DimerConc = Dimerization(Temp, C_PAH, moments, srcActive);

  /* Compute the surface chemistry rate coefficients */
  SurfaceChemistryRates(Temp,C_OH,C_H2O,C_H,C_H2,C_C2H2,C_O2,&ksg,&kox,&ko2);

  /* Compute the soot source terms */
  k = 0;
  if (srcActive[0])
    Nucleation(Temp, DimerConc, moments, &srcP[0*nSootMoments]);
  if (srcActive[1])
    Condensation(Temp, DimerConc, moments, &srcP[1*nSootMoments]);
  if (srcActive[2])
    Coagulation(Temp,Viscosity,Density,MolarMass,moments,&srcP[2*nSootMoments]);
  if (srcActive[3])
    SurfaceGrowth(ksg, moments, &srcP[3*nSootMoments]);
  if (srcActive[4])
    Oxidation(kox, moments,&srcP[4*nSootMoments]);
  if (srcActive[5])
    Fragmentation(ko2, moments,&srcP[5*nSootMoments]);

  /* accumulate the total source terms */
  for (j=0; j<6; j++)
    for (i=0; i<nSootMoments; i++)
      src[i] += srcP[j*nSootMoments+i];

  /* Compute the gas-phase reaction rates for PAH dimerization, surface growth,
     and oxidation */
  GasPhaseRxnRates(Temp,C_PAH,C_OH,C_H2O,C_H,C_H2,C_C2H2,C_O2,moments,rxnrates,srcActive);

  /* Compute the density source term from the gas-phase reaction rates */
  *rhodot = ContinuitySource(rxnrates);
#ifndef RHODOT
  for (i=0; i<nSootMoments; i++) {
      src[i] -= *rhodot * moments[i] / Density * oxidfact[i];
  }
#endif /* RHODOT*/
#ifdef DEBUG
  for (i=0; i<nSootMoments; i++) {
      if (isinf(src[i]) || isnan(src[i])) {
	fprintf(stderr,"HMOM source term NaN or inf: %i %g\n",i, src[i]);
	fprintf(stderr,"moments: %i %g\n",i, moments[i]);
	fprintf(stderr,"oxidfact: %i %g\n",i, oxidfact[i]);
        for (j=0; j<6; j++)
	    fprintf(stderr,"srcP: %i %i %g\n",i,j, srcP[j*nSootMoments+i]);
      }
  }
#endif /* DEBUG */
 
}
/* -------------------------------------------------------------------------- */

/* ---------------------- */
/* Nucleation Source Term */
/* ---------------------- */
void Nucleation(double Temp, double DimerConc, double *moments, double *src)
{
  int i;

  double BetaNucl = GetBetaNucl(Temp);
  double V0 = 2.0*NuclVol;
  double S0 = pow(V0,2.0/3.0);
  /* Total Moments */
  for (i=0; i<nSootMoments-1; i++) {
    src[i] += 0.5*BetaNucl*pow(V0,MomOrderV[i])*
      pow(S0,MomOrderS[i])*DimerConc*DimerConc;
  }

  /* Weight of the Delta Function */
  src[nSootMoments-1] += 0.5*BetaNucl*DimerConc*DimerConc;
}
/* -------------------------------------------------------------------------- */

/* ------------------------------------------------- */
/* Interface Function                                */
/*   Inputs: Gas phase variables in SI units         */
/*   Outputs: Thermodynamic source term coefficients */
/* ------------------------------------------------- */
void sootHMOM_GasPhase_Coeff(double Temp, double Viscosity, double Density, double MolarMass,
              double *c_species, double *ksg, double *kox, double *ko2, double *ProdRate)
{
  int i,j,k;
  double C_PAH[nPAH];
  for (i=0; i<nPAH; i++)
     C_PAH[i] = c_species[i];
  double C_OH   = c_species[nPAH];
  double C_H2O  = c_species[nPAH+1];
  double C_H    = c_species[nPAH+2];
  double C_H2   = c_species[nPAH+3];
  double C_C2H2 = c_species[nPAH+4];
  double C_O2   = c_species[nPAH+5];
  double ProdRate_tmp[3];

  /* convert the molar mass: g/mol -> kg/mol */
  /* MolarMass /= 1000.0;*/

  /* Compute the surface chemistry rate coefficients */
  SurfaceChemistryRates(Temp,C_OH,C_H2O,C_H,C_H2,C_C2H2,C_O2,ksg,kox,ko2);

  /* Compute Dimer production rate */
  ComputeProdRate(C_PAH, Temp, ProdRate_tmp);
  if (ProdRate_tmp[0] != 0.0)
    ProdRate_tmp[0] = ProdRate_tmp[1] / NuclVol;

  /* TODO: Check if this is correct!*/
  *ProdRate = ProdRate_tmp[0];
}
/* -------------------------------------------------------------------------- */
/* ------------------------ */
/* Condensation Source Term */
/* ------------------------ */
void Condensation(double Temp, double DimerConc, double *moments, double *src)
{
  int i;
  double VolTerm, SurfTerm;

  /* Constant for Free Molecular Collisions */
  /*   No van der Waals enhancement         */
  double C = GetC(Temp) / 2.2;

  DimerVol = NuclVol;
  double V0 = 2.0*NuclVol;
  /* Total Moments */
  for (i=0; i<nSootMoments-1; i++) {
    VolTerm = FracMom(MomOrderV[i]+2.0*SootAv-1.0,
                      MomOrderS[i]+2.0*SootAs,moments)*pow(DimerVol, 3.0/6.0)
      +   2.0*FracMom(MomOrderV[i]+SootAv-1.0,   
                      MomOrderS[i]+SootAs,    moments)*pow(DimerVol, 5.0/6.0)
      +       FracMom(MomOrderV[i]-1.0,           
                      MomOrderS[i],           moments)*pow(DimerVol, 7.0/6.0)
      +   0.5*FracMom(MomOrderV[i]+2.0*SootAv-2.0,
                      MomOrderS[i]+2.0*SootAs,moments)*pow(DimerVol, 9.0/6.0)
      +       FracMom(MomOrderV[i]+SootAv-2.0,    
                      MomOrderS[i]+SootAs,    moments)*pow(DimerVol,11.0/6.0)
      +   0.5*FracMom(MomOrderV[i]-2.0,           
                      MomOrderS[i],           moments)*pow(DimerVol,13.0/6.0);
    
    SurfTerm = FracMom(MomOrderV[i]-2.0*SootFitE+2.0*SootAv-1.0,
                       MomOrderS[i]+3.0*SootFitE+2.0*SootAs,
                       moments)*pow(DimerVol, 3.0/6.0)
      +    2.0*FracMom(MomOrderV[i]-2.0*SootFitE+SootAv-1.0,
                       MomOrderS[i]+3.0*SootFitE+SootAs,
                       moments)*pow(DimerVol, 5.0/6.0)
      +        FracMom(MomOrderV[i]-2.0*SootFitE-1.0,
                       MomOrderS[i]+3.0*SootFitE,moments)*pow(DimerVol, 7.0/6.0)
      +    0.5*FracMom(MomOrderV[i]-2.0*SootFitE+2.0*SootAv-2.0,
                       MomOrderS[i]+3.0*SootFitE+2.0*SootAs,
                       moments)*pow(DimerVol, 9.0/6.0)
      +        FracMom(MomOrderV[i]-2.0*SootFitE+SootAv-2.0,
                       MomOrderS[i]+3.0*SootFitE+SootAs,
                       moments)*pow(DimerVol,11.0/6.0)
      +    0.5*FracMom(MomOrderV[i]-2.0*SootFitE-2.0,
                       MomOrderS[i]+3.0*SootFitE,
                       moments)*pow(DimerVol,13.0/6.0);

    src[i] += C*(MomOrderV[i]*VolTerm+SootFitC*MomOrderS[i]*SurfTerm)*DimerConc;
  }

  /* Weight of the Delta Function */
  src[nSootMoments-1] +=  - C * sqrt((1.0/V0) + (1.0/DimerVol))* 
    pow((pow(V0, 1.0/3.0) + pow(DimerVol, 1.0/3.0)), 2.0)*
    DimerConc*moments[nSootMoments-1]*oxidfact[nSootMoments-1];
}
/* -------------------------------------------------------------------------- */

/* ----------------------- */
/* Coagulation Source Term */
/* ----------------------- */
void Coagulation(double Temp,double Viscosity,double Density,double MolarMass,
                 double *moments, double *src)
{
  int i;
  double ss,ss_fm,ss_cn,sl,sl_fm,sl_cn,ll,ll_fm,ll_cn;

  /* Free molecular collision coefficient (w/ van der Waals enhancement) */
  double C_fm = GetC(Temp);

  /* Continuum collision coefficient */
  double C_cn = 8.0*Rgas*Temp/(3.0*Viscosity);

  /* Mean free path for finite Knudsen number correction in continuum regime */
  double Lambda = 3.0*Viscosity/Density * sqrt(Pi*MolarMass/(8.0*Rgas*Temp))/
    pow(6.0*SootMolarMass/(Pi*SootDensity*Avogadro),1.0/3.0);

  double V0 = 2.0*NuclVol;  

  /* Total Moments */
  for (i=0; i<nSootMoments-1; i++) {
    /* Collisions between two first mode particles */
    /*     Collision Model: Pure coalescence       */
    /* - Free molecular regime */
    ss_fm = C_fm*pow(2.0, 2.5)*
      (pow(2.0, MomOrderV[i]+(2.0/3.0)*MomOrderS[i]-1.0) - 1.0)*
      pow(V0, MomOrderV[i]+(2.0/3.0)*MomOrderS[i]+(1.0/6.0))*
      moments[nSootMoments-1]*moments[nSootMoments-1];
    /* - Continuum regime */
    ss_cn = 4.0*C_cn*(1 + 1.257*Lambda*pow(V0,-1.0/3.0))*
      (pow(2.0, MomOrderV[i]+(2.0/3.0)*MomOrderS[i]-1.0) - 1.0)*
      pow(V0, MomOrderV[i]+(2.0/3.0)*MomOrderS[i])*
      moments[nSootMoments-1]*moments[nSootMoments-1];
    /* - Harmonic mean for transitional regime */
    if (ss_fm == 0.0 || ss_cn == 0.0)
      ss = 0.0;
    else
      ss = (ss_fm*ss_cn)/(ss_fm + ss_cn);

    /* Collision between a particle in each mode */
    /*     Collision model: "Splashing"          */
    /* - Free molecular regime */
    sl_fm = C_fm*FMCoagSL(i,moments);
    /* - Continuum regime */
    sl_cn = C_cn*CNCoagSL(i,moments,Temp,Viscosity,Density,MolarMass);
    /* - Harmonic mean for transitional regime */
    if (sl_fm == 0.0 || sl_cn == 0.0)
      sl = 0.0;
    else
      sl = (sl_fm*sl_cn)/(sl_fm + sl_cn);

    /* Collision between two second mode particles */
    /*     Collision Model: Pure aggregation       */
    /* - Free molecular regime */
    ll_fm = C_fm*FMCoagLL(i,moments);
    /* - Continuum regime */
    ll_cn = C_cn*CNCoagLL(i,moments,Temp,Viscosity,Density,MolarMass);
    /* - Harmonic mean for transitional regime */
    if (ll_fm == 0.0 || ll_cn == 0.0)
      ll = 0.0;
    else
      ll = (ll_fm*ll_cn)/(ll_fm + ll_cn);

    src[i] += (ss+sl+ll)*oxidfact[i];
  }

  /* Weight of the Delta Function */
  /* - Free molecular regime */
  ss_fm = -C_fm*pow(2.0,2.5)*pow(1.0/V0,0.5)*pow(V0,2.0/3.0)*
    moments[nSootMoments-1]*moments[nSootMoments-1];
  /* - Continuum regime */
  ss_cn = -4.0*C_cn*(1 + 1.257*Lambda*pow(V0, -1.0/3.0))*
    moments[nSootMoments-1]*moments[nSootMoments-1];
  /* - Harmonic mean for transitional regime */
  if (ss_fm == 0.0 || ss_cn == 0.0)
    ss = 0.0;
  else
    ss = (ss_fm*ss_cn)/(ss_fm + ss_cn);

  /* - Free molecular regime */
  sl_fm = C_fm*FMCoagSL(nSootMoments-1,moments);
  /* - Continuum regime */
  sl_cn = C_cn*CNCoagSL(nSootMoments-1,moments,Temp,Viscosity,Density,
                        MolarMass);
  /* - Harmonic mean for transitional regime */
  sl = (sl_fm * sl_cn)/(sl_fm + sl_cn);

  src[nSootMoments-1] += (ss + sl)*oxidfact[nSootMoments-1];
}
/* -------------------------------------------------------------------------- */

/* -------------------------- */
/* Surface Growth Source Term */
/* -------------------------- */
void SurfaceGrowth(double ksg, double *moments, double *src)
{
  int i;
  /* Convert surface density from SI to mol of C */
  double Chi = SootChi * pow(36.0*Pi,1.0/3.0) * 
    pow(SootMolarMass/(Avogadro*SootDensity),2.0/3.0);
  /* Increase in number of carbon atoms for surface reaction */
  //double dV = NuclVol;
  double dV = 1.0;
  double V0 = 2.0*NuclVol;
  double S0 = pow(V0,2.0/3.0); 
  double sgfact=1.0;
  /* Total Moments */
  for (i=0; i<nSootMoments-1; i++) {
    src[i] += ksg*Chi*dV*
      (  MomOrderV[i]*FracMom(MomOrderV[i]-1.0, MomOrderS[i]+1.0, moments)
       + MomOrderS[i]*SootFitC*FracMom(MomOrderV[i]-1.0-2.0*SootFitE,
                                       MomOrderS[i]+1.0+3.0*SootFitE, moments));
  }

  /* Weight of the Delta Function */
  src[nSootMoments-1] -= ksg*Chi*S0*moments[nSootMoments-1]*oxidfact[nSootMoments-1];
}
/* -------------------------------------------------------------------------- */

/* --------------------- */
/* Oxidation Source Term */
/* --------------------- */
void Oxidation(double kox, double *moments, double *src)
{
  int i;
  /* Convert surface density from SI to mol of C */
  double Chi = SootChi * pow(36.0*Pi,1.0/3.0) * 
    pow(SootMolarMass/(Avogadro*SootDensity),2.0/3.0);
  /* Decrease in number of carbon atoms for oxidation */
  //double dV = NuclVol; /* FlameMaster: 2.0*dimer_nbrC2 */
  double dV = 1.0; /* FlameMaster: 2.0*dimer_nbrC2 */
  double V0 = 2.0*NuclVol;
  double S0 = pow(V0,2.0/3.0); 
  /* Oxidation of the smaller particles */
  double small;
  /* Oxidation of the larger particles */
  double large;
  /* Total Moments */
  for (i=0; i<nSootMoments-1; i++) {
      /* Oxidation of the small particles */
      small = - kox*Chi*pow(V0,MomOrderV[i]-1.0+(2.0/3.0)*
                (MomOrderS[i]+1.0))*moments[nSootMoments-1];
      /* Oxidation of larger particles */
      large = - kox*Chi*(MomOrderV[i]+(2.0/3.0)*MomOrderS[i])*
        FracMomLarge(MomOrderV[i]-1.0,MomOrderS[i]+1.0,moments);

      src[i] += (small + large)*oxidfact[i];
  }

   /* Weight of the Delta Function */
   small = - kox*Chi*pow(V0,-1.0/3.0)*moments[nSootMoments-1];
   double inter = V0/(FracMomLarge(1.0,0.0,moments)/
                           FracMomLarge(0.0,0.0,moments));
   large = inter*kox*Chi*dV*FracMomLarge(-1.0,1.0,moments);

   src[nSootMoments-1] += small *oxidfact[nSootMoments-1]+large*oxidfact[0];
}
/* -------------------------------------------------------------------------- */

/* ------------------------- */
/* Fragmentation Source Term */
/* ------------------------- */
void Fragmentation(double ko2, double *moments, double *src)
{
  int i;
  /* Convert surface density from SI to mol of C */
  double Chi = SootChi*pow(36.0*Pi,1.0/3.0)* 
    pow(SootMolarMass/(Avogadro*SootDensity),2.0/3.0);
  /* Decrease in number of carbon atoms for oxidation */
  //double dV = NuclVol;
  double dV = 1.0;
  double V0 = 2.0*NuclVol;
  double fragfact;
  double TOL=LABSTOL;

  fragfact =(moments[1]-TOL*V0)/moments[1]; 

  /* Total Moments */
  for (i = 0; i < nSootMoments-1; i++) {
    src[i] += (pow(2.0,1.0-MomOrderV[i]-MomOrderS[i])-1.0)*
      2.0*ko2*Chi*dV*FracMomLarge(MomOrderV[i]-1.0,MomOrderS[i]+1.0,moments)*fragfact;
  }

  /* Weight of the Delta Function */
#ifdef FM_VS_MDUC
  double inter = 0.01; // FIXME: from FlameMaster
  //double inter = V0/(FracMomLarge(1.0,0.0,moments)/
  //                           FracMomLarge(0.0,0.0,moments));
#else
  double inter = V0/(FracMomLarge(1.0,0.0,moments)/
                             FracMomLarge(0.0,0.0,moments));
#endif
  src[nSootMoments-1] += inter*2.0*ko2*Chi*dV*FracMomLarge(-1.0,1.0,moments)*fragfact;
}
/* -------------------------------------------------------------------------- */

/* ------------------------------------------------------ */
/* Dimerization: Calculation of dimer concentration       */
/*   Note: Chemical source term for A4 accounted for with */
/*     reaction mechanism                                 */
/* ------------------------------------------------------ */
double Dimerization(double Temp, double *C_PAH, double *moments,int *srcActive)
{
  /* Dimer production rate from pyrene */
  double ProdRate[3];

  ComputeProdRate(C_PAH, Temp, ProdRate);

  if (ProdRate[0] != 0.0)
    ProdRate[0] = ProdRate[1] / NuclVol;

  double BetaNucl=0.0;
  double BetaCond=0.0;
  /* Nucleation and condenstion collision coefficients */
  if (srcActive[0]) BetaNucl = GetBetaNucl(Temp);
  if (srcActive[1]) BetaCond = GetBetaCond(Temp, moments);

  /* QSS for pyrene dimers */
  /*   ProdRate = BetaNucl*[Dimer]*[Dimer] + BetaCond*[Dimer] */
  double Delta = BetaCond*BetaCond + 4.0*BetaNucl*ProdRate[0];
  return (sqrt(Delta)-BetaCond)/(2.0*BetaNucl);
}
/* -------------------------------------------------------------------------- */

/*--------------------------------------------------------*/
/* Computes dimer production rates			  */
/*--------------------------------------------------------*/
void ComputeProdRate(double *C_PAH, double Temp, double *ProdRate)
{
  int i;

  ProdRate[0] = 0.0;
  ProdRate[1] = 0.0;
  ProdRate[2] = 0.0;

  //Frenklach path
  for (i = 0; i < nPAH; i++)
  {
    //double aromconc = density * MAX(Y[AromInd[i]],1.0e-60) / molarMass[AromInd[i]]; //Mueller
    double wDimer = 0.5 * GetBetaDimer(Temp,PAHNbrC2[i]) * C_PAH[i] * C_PAH[i];

    if (wDimer < 0.0)
      wDimer = 0.0;

     //Sticky coeffs
     wDimer *= PAHStick[i];

     ProdRate[0] += wDimer; //Number density
     ProdRate[1] += 2.0 * PAHNbrC2[i] * wDimer; //nbrC2 density
     ProdRate[2] += 2.0 * PAHNbrH[i] * wDimer; //nbrH density
  }
  
  if (ProdRate[0] == 0.0)
    ProdRate[0] = 1.0e-60;
  
}
/*----------------------------------------------------------------------------*/

/* -------------------------------------------------- */
/* Rate coefficients for surface growth and oxidation */
/* -------------------------------------------------- */
void SurfaceChemistryRates(double Temp, double C_OH, double C_H2O,
                           double C_H, double C_H2, double C_C2H2,
                           double C_O2, double *ksg, double *kox, 
                           double *ko2)
{
  /* Rate coefficients */
  double k1f = SootA1f*pow(Temp,Sootn1f)*exp(-SootE1f/(Rgas*Temp));
  double k1b = SootA1b*pow(Temp,Sootn1b)*exp(-SootE1b/(Rgas*Temp));
  double k2f = SootA2f*pow(Temp,Sootn2f)*exp(-SootE2f/(Rgas*Temp));
  double k2b = SootA2b*pow(Temp,Sootn2b)*exp(-SootE2b/(Rgas*Temp));
  double k3f = SootA3f*pow(Temp,Sootn3f)*exp(-SootE3f/(Rgas*Temp));
  double k3b = SootA3b*pow(Temp,Sootn3b)*exp(-SootE3b/(Rgas*Temp));
  double k4  = SootA4 *pow(Temp,Sootn4 )*exp(-SootE4 /(Rgas*Temp));

  /* Fraction of hydrogenated sites that are radical sites */
  /*   QSS for surface radical sites on soot               */
  double fSootStar = (k1f*C_OH + k2f*C_H + k3f) /
    fmax(k1b*C_H2O + k2b*C_H2 + k3b*C_H + k4*C_C2H2, 1.0e-30);
  fSootStar = fSootStar/(1.0 + fSootStar);

  /* Surface chemistry rates [1/s]                             */
  /*   Used in expression for moment source terms              */
  /*   Multiply by concentration of hydrogenated surface sites */
  *ksg = SootA4*pow(Temp,Sootn4)*exp(-SootE4/(Rgas*Temp))*C_C2H2*fSootStar;
  *kox = SootA5*pow(Temp,Sootn5)*exp(-SootE5/(Rgas*Temp))*C_O2*fSootStar
    + (0.5*Avogadro/SootChi)*8.94*SootGamma6*sqrt(Temp)*C_OH;
  *ko2 = SootA5*pow(Temp,Sootn5)*exp(-SootE5/(Rgas*Temp))*C_O2*fSootStar;
}
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* Gas phase source terms from dimerization, surface growth, and oxidation    */
/*   Note: Reaction rates must be multiplied by species stoichiometric        */
/*     coefficients to get species production rates                           */
/* -------------------------------------------------------------------------- */
void GasPhaseRxnRates(double Temp, double *C_PAH, double C_OH, double C_H2O, 
                      double C_H, double C_H2, double C_C2H2, double C_O2, 
                      double *moments, double *rxnrates, int *srcActive)
{
  int i;

  if (srcActive[0]) {/* Nucleation */
#ifndef ONEWAY
     for (i=0; i<nPAH; i++) {
      double wDimer = 0.5*GetBetaDimer(Temp, PAHNbrC2[i])*C_PAH[i]*C_PAH[i];
      
      rxnrates[HMOM_IRXN_PAH+i] = PAHStick[i]*wDimer;
     } 
  }
  else {
#else
     for (i=0; i<nPAH; i++) 
      rxnrates[HMOM_IRXN_PAH+i] = 0.0;
#endif /* ONEWAY */
  }

  /* Convert surface density from SI to mol of C */
  double Chi = SootChi*pow(36.0*Pi,1.0/3.0)*
	pow(SootMolarMass/(Avogadro*SootDensity),2.0/3.0);

  /* Rate coefficients */
  double k1f = SootA1f*pow(Temp,Sootn1f)*exp(-SootE1f/(Rgas*Temp));
  double k1b = SootA1b*pow(Temp,Sootn1b)*exp(-SootE1b/(Rgas*Temp));
  double k2f = SootA2f*pow(Temp,Sootn2f)*exp(-SootE2f/(Rgas*Temp));
  double k2b = SootA2b*pow(Temp,Sootn2b)*exp(-SootE2b/(Rgas*Temp));
  double k3f = SootA3f*pow(Temp,Sootn3f)*exp(-SootE3f/(Rgas*Temp));
  double k3b = SootA3b*pow(Temp,Sootn3b)*exp(-SootE3b/(Rgas*Temp));
  double k4  = SootA4 *pow(Temp,Sootn4 )*exp(-SootE4 /(Rgas*Temp));
  double k5  = SootA5 *pow(Temp,Sootn5 )*exp(-SootE5 /(Rgas*Temp));
  double k6  = 8.94*SootGamma6*sqrt(Temp)*Avogadro;

  /* Concentration of hydrogentated and radical soot surface sites */
  /*   QSS for surface radical sites on soot                       */
  double fSootStar = (k1f*C_OH + k2f*C_H + k3f)/
    fmax(k1b*C_H2O + k2b*C_H2 + k3b*C_H + k4*C_C2H2, 1.0e-30);
  fSootStar = fSootStar/(1.0 + fSootStar);
  //double C_Soot     = Chi*moments[2];
  double C_Soot     = Chi * FracMom(0.0,1.0,moments);
  double C_SootStar = fSootStar*C_Soot;
  C_Soot *= (1.0-fSootStar);

  /* Individual reaction rates */
  double w1f  = k1f*C_Soot*C_OH;
  double w1b  = k1b*C_SootStar*C_H2O;
  double w2f  = k2f*C_Soot*C_H;
  double w2b  = k2b*C_SootStar*C_H2;
  double w3f  = k3f*C_Soot;
  double w3b  = k3b*C_SootStar*C_H;
  double w4   = k4*C_SootStar*C_C2H2;
  double w5   = k5*C_SootStar*C_O2;
  double Surf = pow(36.0*Pi,1.0/3.0)*
    pow(SootMolarMass/(Avogadro*SootDensity),2.0/3.0)*FracMom(0.0,1.0,moments);
  double w6   = k6*C_OH*Surf;
  double fact = 1.0; 

  fact = fmin(fact,oxidfact[0]);
  fact = fmin(fact,oxidfact[1]);
  fact = fmin(fact,oxidfact[2]);
  fact = fmin(fact,oxidfact[3]);
  /* Total reaction rates (forward - backward) [mol / m^3-s] */
  if (srcActive[4]) { /* SurfaceGrowth */
#ifndef ONEWAYSR
     rxnrates[HMOM_IRXN_SG_OH_to_H2O+nPAH-1] = fact*w1f - w1b;
     rxnrates[HMOM_IRXN_SG_H_to_H2+nPAH-1]   = fact*w2f - w2b;
     rxnrates[HMOM_IRXN_SG_H_abs+nPAH-1]     = fact*w3f - w3b;
     rxnrates[HMOM_IRXN_SG_C2H2_add+nPAH-1]  = w4;
  }
  else {
#else
     rxnrates[HMOM_IRXN_SG_OH_to_H2O+nPAH-1] = 0.0;
     rxnrates[HMOM_IRXN_SG_H_to_H2+nPAH-1]   = 0.0;
     rxnrates[HMOM_IRXN_SG_H_abs+nPAH-1]     = 0.0;
     rxnrates[HMOM_IRXN_SG_C2H2_add+nPAH-1]  = 0.0;
#endif
  }
  if (srcActive[3]) { /* Oxidation */ 
#ifndef ONEWAYSR
     rxnrates[HMOM_IRXN_OX_O2+nPAH-1]        = fact*w5;
     rxnrates[HMOM_IRXN_OX_OH+nPAH-1]        = fact*w6;
  }
  else {
#else
     rxnrates[HMOM_IRXN_OX_O2+nPAH-1]        = 0.0;
     rxnrates[HMOM_IRXN_OX_OH+nPAH-1]        = 0.0;
#endif
  }

}
/* -------------------------------------------------------------------------- */

/* ------------------------------- */
/* Continuity Equation Source Term */
/* ------------------------------- */
double ContinuitySource(double *rxnrates)
{
  int i;
  double rhodot = 0.0;

  /* Molar masses of relevant species [kg / mol] */
  //double W_A4   = 0.20224;
  double W_C    = 0.0120107E3;
  double W_OH   = 0.017008E3;
  double W_H2O  = 0.018016E3;
  double W_H    = 0.001008E3;
  double W_H2   = 0.002016E3;
  double W_C2H2 = 0.026036E3;
  double W_O2   = 0.032E3;
  double W_CO   = 0.02801E3;

  /* PAH */
  for (i=0; i<nPAH; i++) {
    double W_PAH = 2.0*W_C*PAHNbrC2[i]+W_H*PAHNbrH[i];
    rhodot -= 2.0*W_PAH*rxnrates[HMOM_IRXN_PAH+i];
  }
  /* OH */
  rhodot += W_OH*(-rxnrates[HMOM_IRXN_SG_OH_to_H2O+nPAH-1]-rxnrates[HMOM_IRXN_OX_OH+nPAH-1]);
  /* H2O */
  rhodot += W_H2O*(rxnrates[HMOM_IRXN_SG_OH_to_H2O+nPAH-1]);
  /* H */
  rhodot += W_H*(-rxnrates[HMOM_IRXN_SG_H_to_H2+nPAH-1]+rxnrates[HMOM_IRXN_SG_H_abs+nPAH-1]);
  /* H2 */
  rhodot += W_H2*(rxnrates[HMOM_IRXN_SG_H_to_H2+nPAH-1]);
  /* C2H2 */
  rhodot += W_C2H2*(-rxnrates[HMOM_IRXN_SG_C2H2_add+nPAH-1]);
  /* O2 */
  rhodot += W_O2*(-rxnrates[HMOM_IRXN_OX_O2+nPAH-1]);
  /* CO */
  rhodot += W_CO*(2.0*rxnrates[HMOM_IRXN_OX_O2+nPAH-1]+rxnrates[HMOM_IRXN_OX_OH+nPAH-1]);
  
  return rhodot;
}
/* -------------------------------------------------------------------------- */

/* ------------------------------------------------------ */
/* Collision coefficients for nucleation and condensation */
/* ------------------------------------------------------ */
double GetBetaNucl(double Temp)
{
  /* Collision frequency (coefficient) between two dimers in the */
  /* free molecular regime with van der Waals enhancement        */
  /*   Units: m^3 / mol-s                                        */
  //return 2.2*4.0*sqrt(2.0)*pow(NuclVol,1.0/6.0)*pow(6.0/Pi,2.0/3.0)* 
  //  sqrt(Pi*Rgas*Temp/(2.0*Avogadro*SootDensity))*
  //  pow(SootMolarMass/(Avogadro*SootDensity),1.0/6.0)*Avogadro;
  double C = GetC(Temp);

  return C * (4.0 * sqrt(2.0) * pow(NuclVol,1.0/6.0));
}
/* -------------------------------------------------------------------------- */

double GetBetaCond(double Temp, double *moments)
{
  /* Collision frequency between a dimer and a particle in the */
  /*    free molecular regime                                  */
  /*   No van der Waals enhancement                            */
  /*   Units: 1 / s                                            */
  //double Cfm = pow(6.0/Pi,2.0/3.0)*
  //  sqrt(Pi*Rgas*Temp/(2.0*Avogadro*SootDensity))*
  //  pow(SootMolarMass/(Avogadro*SootDensity),1.0/6.0)*Avogadro;
  double Cfm = GetC(Temp) / 2.2;

  double S0 = FracMom(2.0*SootAv    ,2.0*SootAs,moments)*pow(NuclVol,-3.0/6.0)
    +     2.0*FracMom(    SootAv    ,    SootAs,moments)*pow(NuclVol,-1.0/6.0)
    +         FracMom(           0.0,       0.0,moments)*pow(NuclVol, 1.0/6.0)
    +     0.5*FracMom(2.0*SootAv-1.0,2.0*SootAs,moments)*pow(NuclVol, 3.0/6.0)
    +         FracMom(    SootAv-1.0,    SootAs,moments)*pow(NuclVol, 5.0/6.0)
    +     0.5*FracMom(          -1.0,       0.0,moments)*pow(NuclVol, 7.0/6.0);
  return Cfm*S0;
}
/* -------------------------------------------------------------------------- */

double GetBetaDimer(double Temp, double i)
{

    double C = GetC(Temp);

    return C * (4.0 * sqrt(2.0) * pow(i,1.0/6.0));
}
double GetC(double temp)
{
  //Units of C are [m^3 / (kmole-s)]
  
  double V1 = SootMolarMass / (Avogadro * SootDensity);
  double A = pow(6.0/Pi, 1.0/3.0);
  
  double C = 2.2 * A * A * sqrt((Pi * Rgas * temp) / (2.0 * Avogadro * SootDensity)) * pow(V1, 1.0/6.0);
  
  return C * Avogadro;
}
/* -------------------- */
/* Moment interpolation */
/* -------------------- */
double FracMom(double VolOrd, double SurfOrd, double *moments)
{
  double M00,M10,M01,M20,M11,M02,N0;
  double fracmom,firstpeak,bothpeaks;
  double fact1,fact2,fact3,fact4,fact5,fact6;
  double V0 = 2.0*NuclVol;
  double S0 = pow(V0,2.0/3.0);
  double dV = 1.0;
  double dS = dV * S0 / V0;
  double NL;
  /* First order interpolation: M00, M10, M01 */
  if (nSootMoments == 4) {
    /* Compute the moments of the second mode */
    N0  = fmax(moments[3]        ,LABSTOL);
    M00 = fmax(moments[0] - N0   ,LABSTOL);
    M10 = fmax(moments[1] - V0*N0,M00*(V0+dV));
    M01 = fmax(moments[2] - S0*N0,M00*(S0+dS));

    /* ensure that large particles are acutally larger than N0 */
    if (M10 < M01 * V0 / S0) M10 = M01 * V0 / S0;

    firstpeak = N0*pow(V0, VolOrd)*pow(S0, SurfOrd);
    /* Contribution from both modes */
    bothpeaks = firstpeak + pow(M00, 1.0-VolOrd-SurfOrd)*pow(M10, VolOrd)*pow(M01, SurfOrd);

    /* If the moments are effectively zero, then use only one mode */
#ifndef BOTHPEAKS
    if (M00 < LABSTOL || M10 < V0*LABSTOL || M01 < S0*LABSTOL)
      fracmom = firstpeak;
    else
#endif
      fracmom = bothpeaks;
    return fracmom;
  }
//  /* First order interpolation: M00, M10, M01 */
//  if (nSootMoments == 4) {
//    /* Compute the moments of the second mode */
//    N0  = fmax(moments[nSootMoments-1], 0.0);
//
//    firstpeak = N0*pow(V0, VolOrd)*pow(S0, SurfOrd);
//
//    NL  = fmax(moments[0] - N0   ,0.0);
//    if (NL>0.0) {
//       M10 = fmax((moments[1] - V0*N0)/NL,0.0);
//       M01 = fmax((moments[2] - S0*N0)/NL,0.0);
//    }
//    else 
//       return firstpeak;
//
//    if (M01 <= 0.0) 
//        return firstpeak;
//    if (M10 <= 0.0) 
//        return firstpeak;
//    /* Contribution from both modes */
//    bothpeaks = firstpeak + NL*pow(M10, VolOrd)*pow(M01, SurfOrd);
//
//    /* If the moments are effectively zero, then use only one mode */
//    return bothpeaks;
//  }
  
  /* Second order interpolation: M00, M10, M01, M20, M11, M02 */
  if (nSootMoments == 7) {
    /* Compute the moments of the second mode */
    M00 = moments[0] - moments[6];
    M10 = moments[1] - V0*moments[6];
    M01 = moments[2] - S0*moments[6];
    M20 = moments[3] - V0*V0*moments[6];
    M11 = moments[4] - V0*S0*moments[6];
    M02 = moments[5] - S0*S0*moments[6];

    /* Contribution from only one mode */
    fact1 = moments[0];
    fact2 = pow(pow(moments[1],  2.0)*pow(moments[0], -1.5)*
                pow(moments[3], -0.5), VolOrd);
    fact3 = pow(pow(moments[2],  2.0)*pow(moments[0], -1.5)*
                pow(moments[5], -0.5), SurfOrd);
    fact4 = pow(pow(moments[3],  0.5)*pow(moments[0],  0.5)*
                pow(moments[1], -1.0), VolOrd*VolOrd);
    fact5 = pow((moments[4]*moments[0])/(moments[1]*moments[2]),VolOrd*SurfOrd);
    fact6 = pow(pow(moments[5],  0.5)*pow(moments[0], 0.5)*
                pow(moments[2], -1.0), SurfOrd*SurfOrd);
    firstpeak = fact1*fact2*fact3*fact4*fact5*fact6;

    /* contribution from both modes */
    fact1 = M00;
    fact2 = pow(pow(M10, 2.0)*pow(M00, -1.5)*pow(M20, -0.5), VolOrd);
    fact3 = pow(pow(M01, 2.0)*pow(M00, -1.5)*pow(M02, -0.5), SurfOrd);
    fact4 = pow(pow(M20, 0.5)*pow(M00,  0.5)*pow(M10, -1.0), VolOrd*VolOrd);
    fact5 = pow((M11*M00)/(M10*M01), VolOrd*SurfOrd);
    fact6 = pow(pow(M02, 0.5)*pow(M00,  0.5)*pow(M01, -1.0), SurfOrd*SurfOrd);
    bothpeaks =   moments[6]*pow(V0,VolOrd)*pow(S0,SurfOrd)
                + fact1*fact2*fact3*fact4*fact5*fact6;

    /* If the moments are effectively zero, then use only one mode */
    if (M00 < 1.0e-30 || M10 < 1.0e-30 || M01 < 1.0e-30 || M20 < 1.0e-30 ||
        M11 < 1.0e-30 || M02 < 1.0e-30)
      fracmom = firstpeak;
    else
      fracmom = bothpeaks;

    return fracmom;
  }

  return 0.0;
}
double FracMomDensity(double VolOrd, double SurfOrd, double *moments, double density)
{
  double M00,M10,M01,M20,M11,M02,N0;
  double fracmom,firstpeak,bothpeaks;
  double fact1,fact2,fact3,fact4,fact5,fact6;
  double V0 = 2.0*NuclVol;
  //double dV = NuclVol;
  double dV = 1.0;
  double S0 = pow(V0,2.0/3.0);
  double dS = dV * S0 / V0;
  int    i;
  double Mmin,x,np,np0,dm,dm0,NL;
  double Df = 1.8;

  for (i=0;i<nSootMoments;i++)
    sootMom[i] = moments[i]*density;

  /* First order interpolation: M00, M10, M01 */
  if (nSootMoments == 4) {
    /* Compute the moments of the second mode */
    N0  = fmax(sootMom[nSootMoments-1], 0.0);

    firstpeak = N0*pow(V0, VolOrd)*pow(S0, SurfOrd);

    NL  = fmax(sootMom[0] - N0   ,0.0);
    if (NL>0.0) {
       M10 = fmax((sootMom[1] - V0*N0)/NL,0.0);
       M01 = fmax((sootMom[2] - S0*N0)/NL,0.0);
    }
    else 
       return firstpeak;

    if (M10 <= V0) 
	M10 = V0;
    if (M01 <= S0) 
	M01 = S0;

    /* Contribution from both modes */
    bothpeaks = firstpeak + NL*pow(M10, VolOrd)*pow(M01, SurfOrd);

    /* If the moments are effectively zero, then use only one mode */
    return bothpeaks;
  }
  
  return 0.0;
}
/* -------------------------------------------------------------------------- */
double FracMomLarge(double VolOrd, double SurfOrd, double *moments)
{
  double V0 = 2.0*NuclVol;
  double N0 =  fmax(moments[nSootMoments-1],LABSTOL);
  /* Subtract out the contribution from the first mode */
  double fracmom = FracMom(VolOrd, SurfOrd, moments)
    - N0*pow(V0,VolOrd)*pow(pow(V0,2.0/3.0),SurfOrd);

  /* If the moment is negative, then return some small (consistent) value */
  double fracmom_min = pow(V0,VolOrd)*pow(pow(V0,2.0/3.0),SurfOrd)*LABSTOL;
  if (fracmom <= fracmom_min)
    return fracmom_min;

  return fracmom;
}
/* -------------------------------------------------------------------------- */

/* ----------------------------------------- */
/* Free molecular coagulation source terms   */
/*     Small-Large: "Splashing"              */
/*     Large-Large: Pure Aggregation         */
/* - Generalized grid function follows terms */
/* ----------------------------------------- */
double FMCoagSL(int i, double *moments)
{
  if (i == nSootMoments-1) /* N0 */
    return -PsiSL(0.0, 0.0, 0.0, 0.0, moments);

  switch (i)
  {
  case 0: /* M00 */
    return -PsiSL(0.0, 0.0, 0.0, 0.0, moments);
  case 1: /* M10 */
    return 0.0;
  case 2: /* M01 */
    return SootFitC*PsiSL(-2.0*SootFitE-1.0,3.0*SootFitE+1.0,1.0,0.0,moments)
      - PsiSL(0.0, 0.0, 0.0, 1.0, moments);
  case 3: /* M20 */
    return 2.0*PsiSL(1.0, 0.0, 1.0, 0.0, moments);
  case 4: /* M11 */
    return SootFitC*PsiSL(-2.0*SootFitE, 3.0*SootFitE+1.0, 1.0, 0.0, moments)
      + PsiSL(0.0, 1.0, 1.0, 0.0, moments)
      + SootFitC*PsiSL(-2.0*SootFitE-1.0, 3.0*SootFitE+1.0, 2.0, 0.0, moments)
      - PsiSL(0.0, 0.0, 1.0, 1.0, moments);
  case 5: /* M02 */
    return  2.0*SootFitC*PsiSL(-2.0*SootFitE-1.0, 3.0*SootFitE+2.0, 1.0, 0.0,
                               moments)
      + SootFitC*SootFitC*PsiSL(-4.0*SootFitE-2.0,6.0*SootFitE+2.0,2.0,0.0,
                                moments)
      - PsiSL(0.0, 0.0, 0.0, 2.0, moments);
  }

  return 0.0;
}
/* -------------------------------------------------------------------------- */

double FMCoagLL(int i, double *moments)
{
  switch (i) {
  case 0: /* M00 */
    return -0.5*PsiLL(0.0, 0.0, 0.0, 0.0, moments);
  case 1: /* M10 */
    return 0.0;
  case 2: /* M01 */
    return 0.0;
  case 3: /* M20 */
    return PsiLL(1.0, 0.0, 1.0, 0.0, moments);
  case 4: /* M11 */
    return PsiLL(1.0, 0.0, 0.0, 1.0, moments);
  case 5: /* M02 */
    return PsiLL(0.0, 1.0, 0.0, 1.0, moments);
  }

  return 0.0;
}
/* -------------------------------------------------------------------------- */

double PsiSL(double x, double y, double a, double b, double *moments)
{
  /* Interpolation for the reduced mass term (square root of sum) in the */
  /*     collision kernel for collision between a particle in each mode  */
  /* - Only two grid functions used for all moments                      */
  /* - Limited sensitivity to increasing the number of grid functions    */

  double N_0  = moments[nSootMoments-1]; /* weight of the delta function */
  double m    = a+2.0/3.0*b;
  double f[3] = { x, SootAv+x, 2.0*SootAv+x };
  double g[3] = { y, SootAs+y, 2.0*SootAs+y };
  double psi1, psi2;
  double V0 = 2.0*NuclVol;

  psi1 =  pow(V0,  1.0/6.0+m)*N_0*FracMomLarge(f[0]-0.5, g[0], moments)
    + 2.0*pow(V0, -1.0/6.0+m)*N_0*FracMomLarge(f[1]-0.5, g[1], moments)
    +     pow(V0, -3.0/6.0+m)*N_0*FracMomLarge(f[2]-0.5, g[2], moments);

  psi2 =  pow(V0,  7.0/6.0+m)*N_0*FracMomLarge(f[0]-0.5, g[0], moments)
    + 2.0*pow(V0,  5.0/6.0+m)*N_0*FracMomLarge(f[1]-0.5, g[1], moments)
    +     pow(V0,  3.0/6.0+m)*N_0*FracMomLarge(f[2]-0.5, g[2], moments)
    +     pow(V0,  1.0/6.0+m)*N_0*FracMomLarge(f[0]+0.5, g[0], moments)
    + 2.0*pow(V0, -1.0/6.0+m)*N_0*FracMomLarge(f[1]+0.5, g[1], moments)
    +     pow(V0, -3.0/6.0+m)*N_0*FracMomLarge(f[2]+0.5, g[2], moments);

  return sqrt(psi1*psi2);
}
/* -------------------------------------------------------------------------- */

double PsiLL(double x, double y, double a, double b, double *moments)
{
  /* Interpolation for the reduced mass (square root of sum) in the collision */
  /*     kernel for collision between two particles in the second mode        */
  /* - Only two grid functions used for all moments                           */
  /* - Limited sensitivity to increasing the number of grid functions         */

  double f[3] = { x, SootAv+x, 2.0*SootAv+x };
  double g[3] = { y, SootAs+y, 2.0*SootAs+y };
  double p[3] = { a, SootAv+a, 2.0*SootAv+a };
  double q[3] = { b, SootAs+b, 2.0*SootAs+b };
  double psi1, psi2;

  psi1 =  FracMomLarge(f[2]-.5,g[2],moments)*FracMomLarge(p[0]-.5,q[0],moments)
    + 2.0*FracMomLarge(f[1]-.5,g[1],moments)*FracMomLarge(p[1]-.5,q[1],moments)
    +     FracMomLarge(f[0]-.5,g[0],moments)*FracMomLarge(p[2]-.5,q[2],moments);

  psi2 =  FracMomLarge(f[2]+.5,g[2],moments)*FracMomLarge(p[0]-.5,q[0],moments)
    + 2.0*FracMomLarge(f[1]+.5,g[1],moments)*FracMomLarge(p[1]-.5,q[1],moments)
    +     FracMomLarge(f[0]+.5,g[0],moments)*FracMomLarge(p[2]-.5,q[2],moments)
    +     FracMomLarge(f[2]-.5,g[2],moments)*FracMomLarge(p[0]+.5,q[0],moments)
    + 2.0*FracMomLarge(f[1]-.5,g[1],moments)*FracMomLarge(p[1]+.5,q[1],moments)
    +     FracMomLarge(f[0]-.5,g[0],moments)*FracMomLarge(p[2]+.5,q[2],moments);

  return sqrt(psi1*psi2);
}
/* -------------------------------------------------------------------------- */

/* ---------------------------------- */
/* Continuum coagulation source terms */
/*     Small-Large: "Splashing"       */
/*     Large-Large: Pure Aggregation  */
/* ---------------------------------- */
double CNCoagSL(int i, double *moments, double Temp, double Viscosity,
                double Density, double MolarMass)
{
  double *m = moments;
  /* weight of the delta function */
  double N_0 = moments[nSootMoments-1];

  /* Mean free path for finite Knudsen number correction in continuum regime */
  double Lambda = 3.0*Viscosity/Density*sqrt(Pi*MolarMass/(8.0*Rgas*Temp))/
    pow(6.0*SootMolarMass/(Pi*SootDensity*Avogadro),1.0/3.0);

  double V0 = 2.0*NuclVol;
  if (i == nSootMoments-1) /* N_0 */
    return -N_0*
      (2.0*FracMomLarge( 0.0       , 0.0       ,m)
       + pow(V0, 1.0/3.0)*FracMomLarge(    -SootAv,    -SootAs,m)
       + pow(V0,-1.0/3.0)*FracMomLarge(     SootAv,     SootAs,m)
       + 1.257*Lambda*
       (  pow(V0,-1.0/3.0)*FracMomLarge(0.0        ,0.0        ,m)
        +                       FracMomLarge(    -SootAv,    -SootAs,m)
        + pow(V0,-2.0/3.0)*FracMomLarge(     SootAv,     SootAs,m)
        + pow(V0, 1.0/3.0)*FracMomLarge(-2.0*SootAv,-2.0*SootAs,m) )
       );
  
  switch (i) {
  case 0: /* M00 */
    return -N_0*
      (2.0*FracMomLarge(0.0,0.0,m)
       + pow(V0, 1.0/3.0)*FracMomLarge(-SootAv,-SootAs,m)
       + pow(V0,-1.0/3.0)*FracMomLarge( SootAv, SootAs,m) +
       1.257*Lambda*
       (  pow(V0,-1.0/3.0)*FracMomLarge( 0.0       , 0.0       ,m)
        +                       FracMomLarge(    -SootAv,    -SootAs,m)
	      + pow(V0,-2.0/3.0)*FracMomLarge(     SootAv,     SootAs,m)
        + pow(V0, 1.0/3.0)*FracMomLarge(-2.0*SootAv,-2.0*SootAs,m) )
       );
  case 1: /* M10 */
    return 0.0;
  case 2: /* M01 */
       return SootFitC*N_0*
       (2.0*V0*FracMomLarge(-2.0*SootFitE-1.0,3.0*SootFitE+1.0,m)
        + pow(V0, 4.0/3.0)*FracMomLarge(-SootAv-2.0*SootFitE-1.0,
                                             -SootAs+3.0*SootFitE+1.0,m)
        + pow(V0, 2.0/3.0)*FracMomLarge( SootAv-2.0*SootFitE-1.0,
                                              SootAs+3.0*SootFitE+1.0,m)
        + 1.257*Lambda*
        (  pow(V0, 2.0/3.0)*FracMomLarge(-2.0*SootFitE-1.0,
                                              3.0*SootFitE+1.0,m)
         +               V0*FracMomLarge(-SootAv-2.0*SootFitE-1.0,
                                              -SootAs+3.0*SootFitE+1.0,m)
         + pow(V0, 1.0/3.0)*FracMomLarge(SootAv-2.0*SootFitE-1.0,
                                              SootAs+3.0*SootFitE+1.0,m)
         + pow(V0, 4.0/3.0)*FracMomLarge(-2.0*SootAv-2.0*SootFitE-1.0,
                                              -2.0*SootAs+3.0*SootFitE+1.0, m) )
        )
       - N_0*
       (2.0*pow(V0, 2.0/3.0)*FracMomLarge(0.0    ,0.0,    m)
        +   pow(V0, 3.0/3.0)*FracMomLarge(-SootAv,-SootAs,m)
        +   pow(V0, 1.0/3.0)*FracMomLarge( SootAv, SootAs,m)
        + 1.257*Lambda*
        (  pow(V0, 1.0/3.0)*FracMomLarge(0.0        ,        0.0,m)
         + pow(V0, 2.0/3.0)*FracMomLarge(    -SootAv,    -SootAs,m)
         + pow(V0, 0.0/3.0)*FracMomLarge(     SootAv,     SootAs,m)
         + pow(V0, 3.0/3.0)*FracMomLarge(-2.0*SootAv,-2.0*SootAs,m) )
        );
  case 3: /* M20 */
    return 2.0*N_0*
       (           2.0*V0*FracMomLarge(1.0,0.0,m)
        + pow(V0,4.0/3.0)*FracMomLarge(-SootAv+1.0,-SootAs,m)
        + pow(V0,2.0/3.0)*FracMomLarge( SootAv+1.0, SootAs,m)
        + 1.257*Lambda*
        (  pow(V0,2.0/3.0)*FracMomLarge(1.0, 0.0,m)
         +              V0*FracMomLarge(    -SootAv+1.0,-SootAs,m)
         + pow(V0,1.0/3.0)*FracMomLarge(     SootAv+1.0,     SootAs,m)
         + pow(V0,4.0/3.0)*FracMomLarge(-2.0*SootAv+1.0,-2.0*SootAs,m) )
        );
  case 4: /* M11 */
    return SootFitC*N_0*
       (2.0*V0*FracMomLarge(-2.0*SootFitE,3.0*SootFitE+1.0,m)
        + pow(V0, 4.0/3.0)*FracMomLarge(-SootAv-2.0*SootFitE,
                                             -SootAs+3.0*SootFitE+1.0,m)
        + pow(V0, 2.0/3.0)*FracMomLarge( SootAv-2.0*SootFitE,
                                              SootAs+3.0*SootFitE+1.0,m)
        + 1.257*Lambda*
        (  pow(V0, 2.0/3.0)*FracMomLarge(-2.0*SootFitE,3.0*SootFitE+1.0,m)
         +               V0*FracMomLarge(-SootAv-2.0*SootFitE,
                                              -SootAs+3.0*SootFitE+1.0,m)
         + pow(V0, 1.0/3.0)*FracMomLarge( SootAv-2.0*SootFitE,
                                               SootAs+3.0*SootFitE+1.0,m)
         + pow(V0, 4.0/3.0)*FracMomLarge(-2.0*SootAv-2.0*SootFitE,
                                              -2.0*SootAs+3.0*SootFitE+1.0,m) )
        )
        + N_0*
        (2.0*V0*FracMomLarge(0.0,1.0,m)
         + pow(V0, 4.0/3.0)*FracMomLarge(-SootAv,-SootAs+1.0,m)
         + pow(V0, 2.0/3.0)*FracMomLarge( SootAv, SootAs+1.0,m)
         + 1.257*Lambda*
         (  pow(V0, 2.0/3.0)*FracMomLarge(0.0,1.0,m)
          +               V0*FracMomLarge(-SootAv,-SootAs+1.0,m)
          + pow(V0, 1.0/3.0)*FracMomLarge( SootAv, SootAs+1.0,m)
          + pow(V0, 4.0/3.0)*FracMomLarge(-2.0*SootAv,-2.0*SootAs+1.0,m) )
         )
        + SootFitC*N_0*
        (2.0*pow(V0, 6.0/3.0)*FracMomLarge(-2.0*SootFitE-1.0,
                                                 3.0*SootFitE+1.0,m)
         +   pow(V0, 7.0/3.0)*FracMomLarge(-SootAv-2.0*SootFitE-1.0,
                                                -SootAs+3.0*SootFitE+1.0,m)
         +   pow(V0, 5.0/3.0)*FracMomLarge( SootAv-2.0*SootFitE-1.0,
                                                 SootAs+3.0*SootFitE+1.0,m)
         + 1.257*Lambda*
         (  pow(V0, 5.0/3.0)*FracMomLarge(-2.0*SootFitE-1.0,
                                                3.0*SootFitE+1.0,m)
          + pow(V0, 6.0/3.0)*FracMomLarge(-SootAv-2.0*SootFitE-1.0,
                                               -SootAs+3.0*SootFitE+1.0,m)
          + pow(V0, 4.0/3.0)*FracMomLarge( SootAv-2.0*SootFitE-1.0,
                                                SootAs+3.0*SootFitE+1.0,m)
          + pow(V0, 7.0/3.0)*FracMomLarge(-2.0*SootAv-2.0*SootFitE-1.0,
                                               -2.0*SootAs+3.0*SootFitE+1.0,m) )
         )
        - N_0*
        (2.0*pow(V0, 5.0/3.0)*FracMomLarge(0.0,0.0,m)
         + pow(V0, 6.0/3.0)*FracMomLarge(-SootAv,-SootAs,m)
         + pow(V0, 4.0/3.0)*FracMomLarge( SootAv, SootAs,m)
         + 1.257*Lambda*
         (  pow(V0, 4.0/3.0)*FracMomLarge(0.0,0.0,m)
          + pow(V0, 5.0/3.0)*FracMomLarge(    -SootAv,    -SootAs,m)
          + pow(V0, 3.0/3.0)*FracMomLarge(     SootAv,     SootAs,m)
          + pow(V0, 6.0/3.0)*FracMomLarge(-2.0*SootAv,-2.0*SootAs,m) )
         );
  case 5: /* M02 */
    return 2.0*N_0*
        (2.0*SootFitC*
         (2.0*V0*FracMomLarge(-2.0*SootFitE-1.0,3.0*SootFitE+2.0,m)
          + pow(V0, 4.0/3.0)*FracMomLarge(-SootAv-2.0*SootFitE-1.0,
                                               -SootAs+3.0*SootFitE+2.0,m)
          + pow(V0, 2.0/3.0)*FracMomLarge( SootAv-2.0*SootFitE-1.0,
                                                SootAs+3.0*SootFitE+2.0,m)
          + 1.257*Lambda*
          (  pow(V0, 2.0/3.0)*FracMomLarge(-2.0*SootFitE-1.0,
                                                 3.0*SootFitE+2.0,m)
           +               V0*FracMomLarge(    -SootAv-2.0*SootFitE-1.0,
                                                    -SootAs+3.0*SootFitE+2.0,m)
           + pow(V0, 1.0/3.0)*FracMomLarge(     SootAv-2.0*SootFitE-1.0,
                                                     SootAs+3.0*SootFitE+2.0,m)
           + pow(V0, 4.0/3.0)*FracMomLarge(-2.0*SootAv-2.0*SootFitE-1.0,
                                                -2.0*SootAs+3.0*SootFitE+2.0,m))
          )
         + SootFitC*SootFitC*N_0*
         (2.0*pow(V0, 6.0/3.0)*FracMomLarge(-4.0*SootFitE-2.0,
                                                  6.0*SootFitE+2.0,m)
          + pow(V0, 7.0/3.0)*FracMomLarge(-SootAv-4.0*SootFitE-2.0,
                                               -SootAs+6.0*SootFitE+2.0,m)
          + pow(V0, 5.0/3.0)*FracMomLarge(SootAv-4.0*SootFitE-2.0,
                                               SootAs+6.0*SootFitE+2.0,m)
          + 1.257*Lambda*
          (  pow(V0, 5.0/3.0)*FracMomLarge(-4.0*SootFitE-2.0,
                                         6.0*SootFitE+2.0,m)
           + pow(V0, 6.0/3.0)*FracMomLarge(-SootAv-4.0*SootFitE-2.0,
                                         -SootAs+6.0*SootFitE+2.0,m)
           + pow(V0, 4.0/3.0)*FracMomLarge( SootAv-4.0*SootFitE-2.0,
                                          SootAs+6.0*SootFitE+2.0,m)
           + pow(V0, 7.0/3.0)*FracMomLarge(-2.0*SootAv-4.0*SootFitE-2.0,
                                                -2.0*SootAs+6.0*SootFitE+2.0,m))
          )
         - N_0*
         (2.00*pow(V0, 4.0/3.0)*FracMomLarge(0.0,0.0,m)
          +    pow(V0, 5.0/3.0)*FracMomLarge(-SootAv,-SootAs,m)
          +    pow(V0, 3.0/3.0)*FracMomLarge( SootAv, SootAs,m)
          + 1.257*Lambda*
          (  pow(V0, 3.0/3.0)*FracMomLarge(0.0,0.0,m)
           + pow(V0, 4.0/3.0)*FracMomLarge(-SootAv,-SootAs,m)
           + pow(V0, 2.0/3.0)*FracMomLarge( SootAv, SootAs,m)
           + pow(V0, 5.0/3.0)*FracMomLarge(-2.0*SootAv,-2.0*SootAs,m) )
          )
         );
  }

  return 0.0;
}
/* -------------------------------------------------------------------------- */

double CNCoagLL(int i, double *moments, double Temp, double Viscosity,
                double Density, double MolarMass)
{
  double *m = moments;
  /* weight of the delta function */
  double N_0 = moments[nSootMoments-1];

  /* Mean free path for finite Knudsen number correction in continuum regime */
  double Lambda = 3.0*Viscosity/Density * sqrt(Pi*MolarMass/(8.0*Rgas*Temp))/
    pow(6.0*SootMolarMass/(Pi*SootDensity*Avogadro),1.0/3.0);

  switch (i) {
  case 0: /* M00 */
    return -0.5*
      (2.0*FracMomLarge(    0.0,    0.0,m)*FracMomLarge(    0.0,    0.0,m)
       +   FracMomLarge( SootAv, SootAs,m)*FracMomLarge(-SootAv,-SootAs,m)
       +   FracMomLarge(-SootAv,-SootAs,m)*FracMomLarge( SootAv, SootAs,m)
       + 1.257*Lambda*
       (  FracMomLarge(-SootAv,-SootAs,m)*FracMomLarge(0.0,0.0,m)
        + FracMomLarge(0.0,0.0,m)*FracMomLarge(-SootAv,-SootAs,m)
        + FracMomLarge(-2.0*SootAv,-2.0*SootAs,m)*FracMomLarge(SootAv,SootAs,m)
        + FracMomLarge(SootAv,SootAs,m)*FracMomLarge(-2.0*SootAv,-2.0*SootAs,m))
       );
  case 1: /* M10 */
    return 0.0;
  case 2: /* M01 */
    return 0.0;
  case 3: /* M20 */
    return 2.0*FracMomLarge(1.0,0.0,m)*FracMomLarge(1.0,0.0,m)
      + FracMomLarge( SootAv+1.0, SootAs,m)*FracMomLarge(-SootAv+1.0,-SootAs,m)
      + FracMomLarge(-SootAv+1.0,-SootAs,m)*FracMomLarge( SootAv+1.0, SootAs,m)
      + 1.257*Lambda*
      (  FracMomLarge(-SootAv+1.0,-SootAs,m)*FracMomLarge(1.0,0.0,m)
       + FracMomLarge(1.0,0.0,m)*FracMomLarge(-SootAv+1.0,-SootAs,m)
       + (FracMomLarge(-2.0*SootAv+1.0,-2.0*SootAs,m)*
          FracMomLarge(SootAv+1.0,SootAs,m))
       + (FracMomLarge(SootAv+1.0,SootAs,m)*
          FracMomLarge(-2.0*SootAv+1.0,-2.0*SootAs,m))
         );
      (  FracMomLarge(-SootAv+1.0,-SootAs,m)*FracMomLarge(1.0,0.0,m)
       + FracMomLarge(1.0,0.0,m)*FracMomLarge(-SootAv+1.0,-SootAs,m)
       + (FracMomLarge(-2.0*SootAv+1.0,-2.0*SootAs,m)*
          FracMomLarge(SootAv+1.0,SootAs,m))
       + (FracMomLarge(SootAv+1.0,SootAs,m)*
          FracMomLarge(-2.0*SootAv+1.0,-2.0*SootAs,m))
         );
  case 4: /* M11 */
    return 2.0*FracMomLarge(1.0,0.0,m)*FracMomLarge(0.0,1.0,m)
      + FracMomLarge( SootAv+1.0, SootAs,m)*FracMomLarge(-SootAv,-SootAs+1.0,m)
      + FracMomLarge(-SootAv+1.0,-SootAs,m)*FracMomLarge( SootAv, SootAs+1.0,m)
      + 1.257*Lambda*
      (  FracMomLarge(-SootAv+1.0,-SootAs,m)*FracMomLarge(0.0,1.0,m)
       + FracMomLarge(1.0,0.0,m)*FracMomLarge(-SootAv,-SootAs+1.0,m)
       + (FracMomLarge(-2.0*SootAv+1.0,-2.0*SootAs,m)*
          FracMomLarge(SootAv,SootAs+1.0,m))
       + (FracMomLarge(SootAv+1.0,SootAs,m)*
          FracMomLarge(-2.0*SootAv,-2.0*SootAs+1.0,m))
         );
  case 5: /* M02 */
    return 2.0*FracMomLarge(0.0,1.0,m)*FracMomLarge(0.0,1.0,m)
      + FracMomLarge( SootAv, SootAs+1.0,m)*FracMomLarge(-SootAv,-SootAs+1.0,m)
      + FracMomLarge(-SootAv,-SootAs+1.0,m)*FracMomLarge( SootAv, SootAs+1.0,m)
      + 1.257*Lambda*
      (  FracMomLarge(-SootAv,-SootAs+1.0,m)*FracMomLarge(0.0,1.0,m)
       + FracMomLarge(0.0,1.0,m)*FracMomLarge(-SootAv,-SootAs+1.0,m)
       + (FracMomLarge(-2.0*SootAv,-2.0*SootAs+1.0,m)*
          FracMomLarge(SootAv,SootAs+1.0,m))
       + (FracMomLarge(SootAv,SootAs+1.0,m)*
          FracMomLarge(-2.0*SootAv,-2.0*SootAs+1.0,m))
         );
  }

  return 0.0;
}
/* -------------------------------------------------------------------------- */
