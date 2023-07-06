#include "udfWiebe.h"

extern udfStruct udf;
static wiebeStruct *wiebe = NULL;

void UDFinitWiebe()
{ 
  Message0("===== Initializing Wiebe combustion ===== \n\n");

  /* allocate the memory */
  wiebe = (wiebeStruct*) calloc(1, sizeof(wiebeStruct));
  
#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      UDFinitWiebeFile(wiebe);
    }

  /* ----- fuel properties -----*/
  udf.LHV  = 42.01e6;
  if (RP_is_Defined("mduc/lower-heating-value"))
    udf.LHV = RP_Get_Real("mduc/lower-heating-value");
  udf.m_fuel = 9.6e-6/2.0;
  if (RP_is_Defined("mduc/mass-fuel"))
    udf.m_fuel = RP_Get_Real("mduc/mass-fuel");

  Message0("Fuel Properties:\n");
  Message0("  injected mass: %15.6e [kg]\n", udf.m_fuel);
  Message0("            LHV: %15.6e [J/kg]\n", udf.LHV);
  Message0("\n");

  /* ----- coefficients for the Wiebe function ----- */
  if (RP_is_Defined("wiebe/coeff-a"))
    wiebe->a = RP_Get_Real("wiebe/coeff-a");
  if (wiebe->a < 0.0) {
    Message0("Warning: using default wiebe->a\n");
    wiebe->a = 3.9544;
  }

  if (RP_is_Defined("wiebe/coeff-m"))
    wiebe->m = RP_Get_Real("wiebe/coeff-m");
  if (wiebe->m < 0.0) {
    Message0("Warning: using default wiebe->m\n");
    wiebe->m = 1.8537;
  }
  if (RP_is_Defined("wiebe/soc")) /* start of combustion */
    wiebe->SOC = RP_Get_Real("wiebe/soc");

  /* compute the burn duration */
  if (RP_is_Defined("wiebe/soc-alpha"))
    wiebe->alpha = RP_Get_Real("wiebe/soc-alpha");
  if (wiebe->alpha < 0.0) {
    Message0("Warning: using default wiebe->alpha\n");
    wiebe->alpha = 0.76;
  }
  if (RP_is_Defined("wiebe/soc-beta"))
    wiebe->beta = RP_Get_Real("wiebe/soc-beta");
  if (wiebe->beta < 0.0) {
    Message0("Warning: using default wiebe->beta\n");
    wiebe->beta = 20;
  }
  wiebe->dt = wiebe->alpha*wiebe->SOC + wiebe->beta;

  wiebe->sFuel = wiebe->sO2 = wiebe->sCO2 = wiebe->sH2O = wiebe->sHTOT = 0.0;
  wiebe->bFuel = wiebe->bO2 = wiebe->bCO2 = wiebe->bH2O = wiebe->bHTOT = 0.0;

  Message0("Wiebe parameters:\n");
  Message0("                    a: %8.6f [-]\n", wiebe->a);
  Message0("                    m: %8.6f [-]\n", wiebe->m);
  Message0("                alpha: %8.6f [-]\n", wiebe->alpha);
  Message0("                 beta: %8.6f [-]\n", wiebe->beta);
  Message0("  start of combustion: %8.6f [CAD]\n", wiebe->SOC);
  Message0("  ccmbustion duration: %8.6f [CAD]\n", wiebe->dt);
  Message0("\n");


  if (RP_is_Defined("wiebe/rxn-c"))
    wiebe->nC = RP_Get_Real("wiebe/rxn-c");
  if (wiebe->nC < 0.0) {
    Message0("Warning: using default wiebe->rxn-c\n");
    wiebe->nC = 6.9;
  }
  if (RP_is_Defined("wiebe/rxn-h"))
    wiebe->nH = RP_Get_Real("wiebe/rxn-h");
  if (wiebe->nH < 0.0) {
    Message0("Warning: using default wiebe->rxn-h\n");
    wiebe->nH = 14.0;
  }

  Message0("Reaction Properties:\n");
  Message0("  C: %8.3f [-]\n", wiebe->nC);
  Message0("  H: %8.3f [-]\n", wiebe->nH);
  Message0("\n");

  /* udf.Mw_fuel = wiebe->nC*12 + wiebe->nH*1.0; */

  /* compute the coefficients for the species equations [kg/mol] */
  wiebe->aFuel = udf.Mw_fuel/1000; 
  wiebe->aO2   = udf.Mw_O2/1000*(wiebe->nC+0.25*wiebe->nH);
  wiebe->aH2O  = udf.Mw_H2O/1000*0.5*wiebe->nH;
  wiebe->aCO2  = udf.Mw_CO2/1000*wiebe->nC;

  Message0("===== finished Wiebe initialization ===== \n");
}
/* -------------------------------------------------------------------------- */

void UDFinitWiebeFile(wiebeStruct* w)
{
  /* generate the file name */
  int i, nT, cur_n;
  char fileName[512], *id, fileLine[2056], dummy[2056], *prev;
  long nPrev;
  FILE *pFile;

  /* open state file and write header */
  sprintf(fileName, "%s/monitor/wiebe.%s", udf.caseDir,
          RP_Get_String("mduc/case-id"));

  /* check to see if the file exists */
  pFile = NULL;
  pFile = fopen(fileName, "r");
  if (pFile) {
    Message0("  wiebe monitor %s exists\n", fileName);
    
    /* read in the file and search for the time step */
    cur_n = -1;
    while (!feof(pFile)) {
      nPrev = ftell(pFile);
      fgets(fileLine, 2056, pFile);
      
      if (fileLine[0] != '#' && fileLine[0] != '\n')
        sscanf(fileLine, "%i %*s", &cur_n);
      
      if (cur_n >= N_TIME)
        break;
    } /* end loop over file */
    
      /* go to beginning of file */
    fseek (pFile, 0, SEEK_SET);
    
    /* read in everything up until the current time step */
    prev = (char*) calloc(nPrev, sizeof(char));
    fread(prev,1,nPrev,pFile);
    
    /* reopen the file and write the previous information */
    freopen(fileName, "w", pFile);
    fwrite(prev,1,nPrev,pFile);
    
    free(prev);
  }
  else { /* create the file new and write the header */
    pFile = fopen(fileName, "w");
    
    nT = 3;
    fprintf(pFile, "# nt<1> time<2> ");
    if (in_cylinder_p) {
      fprintf(pFile, "CAD<3>  CAD_rel<4> ");
      nT += 2;
    }

    fprintf(pFile, "mass[kg]<%i> ",             nT); nT++;
    fprintf(pFile, "m_fuel[kg]<%i> ",           nT); nT++;
    fprintf(pFile, "dNb/dt[mol/s]<%i> ",        nT); nT++;
    fprintf(pFile, "srcFuelAvg[kg/m^3-s]<%i> ", nT); nT++;
    fprintf(pFile, "srcO2Avg[kg/m^3-s]<%i> ",   nT); nT++;
    fprintf(pFile, "srcH2OAvg[kg/m^3-s]<%i> ",  nT); nT++;
    fprintf(pFile, "srcCO2Avg[kg/m^3-s]<%i> ",  nT); nT++;
    fprintf(pFile, "srcEnthAvg[W/m^3-s]<%i> ",  nT); nT++;
    fprintf(pFile, "fracAvg[]<%i> ",            nT); nT++;
    fprintf(pFile, "srcFuelInt[kg/s]<%i> ",     nT); nT++;
    fprintf(pFile, "srcO2Int[kg/s]<%i> ",       nT); nT++;
    fprintf(pFile, "srcH2OInt[kg/s]<%i> ",      nT); nT++;
    fprintf(pFile, "srcCO2Int[kg/s]<%i> ",      nT); nT++;
    fprintf(pFile, "srcEnthInt[W/s]<%i> ",      nT); nT++;
    fprintf(pFile, "fracInt[]<%i> ",            nT); nT++;
    fprintf(pFile, "FuelInt[kg]<%i> ",          nT); nT++;
    fprintf(pFile, "O2Int[kg]<%i> ",            nT); nT++;
    fprintf(pFile, "H2OInt[kg]<%i> " ,          nT); nT++;
    fprintf(pFile, "CO2Int[kg]<%i> ",           nT); nT++;
    fprintf(pFile, "EnthInt[W]<%i> ",           nT); nT++;
    fprintf(pFile, "Fuel[kg]<%i> ",             nT); nT++;
    fprintf(pFile, "O2[kg]<%i> ",               nT); nT++;
    fprintf(pFile, "H2O[kg]<%i> " ,             nT); nT++;
    fprintf(pFile, "CO2[kg]<%i> ",              nT); nT++;
    fprintf(pFile, "Enth[W]<%i> ",              nT); nT++;
    fprintf(pFile, "\n");
  }
  fflush(pFile);

  w->oFile = pFile;
}
/* -------------------------------------------------------------------------- */

void WiebeCalcSrc(Domain *d)
{
  int    iT;
  Thread *t;
  cell_t  c;

  /* calculate Wiebe function source */
  double CAD     = RP_Get_Real("dynamesh/in-cyn/crank-angle");
  double CAD_rel = fmod(CAD,720);

  /* ensure that relative crank angle is defined centred around combustion at
     0CAD, i.e. -360 < CAD_rel < 360 */
  if (CAD_rel > 360.0)
    CAD_rel -= 720.0;

  /* keep track of total masses in the cylinder */
  wiebe->mf_n = 0.0;
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      wiebe->mf_n += C_YI(c,t,udf.iFuel)*C_R(c,t)*C_VOLUME(c,t);
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

#if RP_NODE /* parallel */
  wiebe->mf_n = PRF_GRSUM1(wiebe->mf_n);
#endif /* RP_NODE */
  
  /* find out if combustion has started */
  wiebe->src = 0.0;
  if ((CAD_rel >= wiebe->SOC)) {
    wiebe->src=UDFcalcWiebe_dt(wiebe->a,wiebe->m,wiebe->SOC,wiebe->dt,CAD_rel);

    /* convert from d/dCAD to d/dt */
    wiebe->src *= 6*in_cyn_ca_rpm;

    /* account for initial moles of fuel in the system */
    wiebe->src *= udf.m_fuel/udf.Mw_fuel*1000;
  }
  else {
   
  }
}
/* -------------------------------------------------------------------------- */

void writeWiebe(Domain *d)
{
  int    iT;
  Thread *t;
  cell_t  c;
  double frac, aO2, aCO2, aH2O, mass, m_fuel, totMass, totVol;
  double avgFuel = 0.0, avgO2 = 0.0, avgH2O = 0.0, avgCO2 = 0.0, avgHTOT = 0.0;
  double iFuel = 0.0, iO2 = 0.0, iH2O = 0.0, iCO2 = 0.0, iHTOT = 0.0;
  double sFuel, sO2, sH2O, sCO2, sHTOT;
  double avgFrac = 0.0, iFrac;

  /* calculate Wiebe function source */
  double CAD     = RP_Get_Real("dynamesh/in-cyn/crank-angle");
  double CAD_rel = fmod(CAD,720);

  /* ensure that relative crank angle is defined centred around combustion at
     0CAD, i.e. -360 < CAD_rel < 360 */
  if (CAD_rel > 360.0)
    CAD_rel -= 720.0;

   /* store the source terms to output */
  totMass = m_fuel = totVol = 0.0;
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;

    begin_c_loop_int(c,t) {
      mass = C_R(c,t)*C_VOLUME(c,t);

      /* compute the coefficients */
      frac = C_YI(c,t,udf.iFuel)*mass/udf.m_fuel;
      frac = C_YI(c,t,udf.iFuel)*mass/wiebe->mf_n;

      /* compute the source terms */
      sFuel = -frac*wiebe->src*wiebe->aFuel/C_VOLUME(c,t); 
      sO2   = -frac*wiebe->src*wiebe->aO2/C_VOLUME(c,t);
      sH2O  =  frac*wiebe->src*wiebe->aH2O/C_VOLUME(c,t);
      sCO2  =  frac*wiebe->src*wiebe->aCO2/C_VOLUME(c,t);
      sHTOT =  frac*wiebe->aFuel*udf.LHV*wiebe->src/C_VOLUME(c,t);

      /* keep track of the integral */
      iFuel += C_VOLUME(c,t)*sFuel;
      iO2   += C_VOLUME(c,t)*sO2;
      iH2O  += C_VOLUME(c,t)*sH2O;
      iCO2  += C_VOLUME(c,t)*sCO2;
      iHTOT += C_VOLUME(c,t)*sHTOT;
      iFrac += frac;

      /* keep track of the mass average */
      avgFuel += mass*sFuel;
      avgO2   += mass*sO2;
      avgH2O  += mass*sH2O;
      avgCO2  += mass*sCO2;
      avgHTOT += mass*sHTOT;
      avgFrac += mass*frac;

      totMass += mass;
      m_fuel  += mass*C_YI(c,t,udf.iFuel);
      totVol  += C_VOLUME(c,t);
    } end_c_loop_int(c,t);
  } /* end loop over threads */

#if RP_NODE /* parallel */
  avgFuel = PRF_GRSUM1(avgFuel);
  avgO2   = PRF_GRSUM1(avgO2);
  avgH2O  = PRF_GRSUM1(avgH2O);
  avgCO2  = PRF_GRSUM1(avgCO2);
  avgHTOT = PRF_GRSUM1(avgHTOT);
  avgFrac = PRF_GRSUM1(avgFrac);
  iFuel   = PRF_GRSUM1(iFuel);
  iO2     = PRF_GRSUM1(iO2);
  iH2O    = PRF_GRSUM1(iH2O);
  iCO2    = PRF_GRSUM1(iCO2);
  iHTOT   = PRF_GRSUM1(iHTOT);
  iFrac   = PRF_GRSUM1(iFrac);
  totMass = PRF_GRSUM1(totMass);
  m_fuel  = PRF_GRSUM1(m_fuel);
  totVol  = PRF_GRSUM1(totVol);
#endif /* RP_NODE */
  avgFuel /= totMass;
  avgO2   /= totMass;
  avgH2O  /= totMass;
  avgCO2  /= totMass;
  avgHTOT /= totMass;
  avgFrac /= totMass;

  /* keep track of the integral values */
  wiebe->bFuel += CURRENT_TIMESTEP*iFuel;
  wiebe->bO2   += CURRENT_TIMESTEP*iO2;
  wiebe->bH2O  += CURRENT_TIMESTEP*iH2O;
  wiebe->bCO2  += CURRENT_TIMESTEP*iCO2;
  wiebe->bHTOT += CURRENT_TIMESTEP*iHTOT;

  wiebe->sFuel += CURRENT_TIMESTEP*wiebe->src*wiebe->aFuel;
  wiebe->sO2   += CURRENT_TIMESTEP*wiebe->src*wiebe->aO2;
  wiebe->sH2O  += CURRENT_TIMESTEP*wiebe->src*wiebe->aH2O;
  wiebe->sCO2  += CURRENT_TIMESTEP*wiebe->src*wiebe->aCO2;
  wiebe->sHTOT += CURRENT_TIMESTEP*wiebe->src*wiebe->aFuel*udf.LHV;

  /* write the info to the output file */
#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      fprintf(wiebe->oFile, "%6i %15.6e %15.6e %15.6e %15.6e %15.6e  %15.6e ",
              N_TIME, CURRENT_TIME, CAD, CAD_rel, totMass, m_fuel, wiebe->src);
      fprintf(wiebe->oFile, "%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e ",
              avgFuel, avgO2, avgH2O, avgCO2, avgHTOT, avgFrac);
      fprintf(wiebe->oFile, "%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e ",
              iFuel, iO2, iH2O, iCO2, iHTOT, iFrac);
      fprintf(wiebe->oFile, "%15.6e %15.6e %15.6e %15.6e %15.6e ",
              wiebe->bFuel, wiebe->bO2, wiebe->bH2O, wiebe->bCO2, wiebe->bHTOT);
      fprintf(wiebe->oFile, "%15.6e %15.6e %15.6e %15.6e %15.6e ",
              wiebe->sFuel, wiebe->sO2, wiebe->sH2O, wiebe->sCO2, wiebe->sHTOT);
      fprintf(wiebe->oFile, "\n");
      fflush(wiebe->oFile);
    }
}
/* -------------------------------------------------------------------------- */

double UDFcalcWiebe(double a, double m, double t0, double dt, double t)
{
  if (t > t0)
    return 1.0-exp(-a*pow((t-t0)/dt, m+1));
  else
    return 0.0;
}
/* ------------------------------------------------------------------------- */

double UDFcalcWiebe_dt(double a, double m, double t0, double dt, double t)
{
  if (t > t0)
    return (1.0-UDFcalcWiebe(a,m,t0,dt,t))*a*(m+1)/dt*pow((t-t0)/dt,m);
  else 
    return 0.0;
}
/* ------------------------------------------------------------------------- */

/* Define the source terms for the major species to be used */
DEFINE_SOURCE(Wiebe_Y_fuel, c, t, dS, eqn)
{
  double f = 0.0;
  if (wiebe->mf_n > 1.0e-7)
    f = C_YI(c,t,udf.iFuel)*C_R(c,t)*C_VOLUME(c,t)/wiebe->mf_n;
  return -f*wiebe->src*wiebe->aFuel/C_VOLUME(c,t);
}
/* ------------------------------------------------------------------------- */

DEFINE_SOURCE(Wiebe_Y_O2, c, t, dS, eqn)
{
  double f = 0.0;
  if (wiebe->mf_n > 1.0e-7)
    f = C_YI(c,t,udf.iFuel)*C_R(c,t)*C_VOLUME(c,t)/wiebe->mf_n;
  return -f*wiebe->src*wiebe->aO2/C_VOLUME(c,t);
}
/* ------------------------------------------------------------------------- */

DEFINE_SOURCE(Wiebe_Y_CO2, c, t, dS, eqn)
{
  double f = 0.0;
  if (wiebe->mf_n > 1.0e-7)
    f = C_YI(c,t,udf.iFuel)*C_R(c,t)*C_VOLUME(c,t)/wiebe->mf_n;
  return f*wiebe->src*wiebe->aCO2/C_VOLUME(c,t);
}
/* ------------------------------------------------------------------------- */

DEFINE_SOURCE(Wiebe_Y_H2O, c, t, dS, eqn)
{
  double f = 0.0;
  if (wiebe->mf_n > 1.0e-7)
    f = C_YI(c,t,udf.iFuel)*C_R(c,t)*C_VOLUME(c,t)/wiebe->mf_n;
  return f*wiebe->src*wiebe->aH2O/C_VOLUME(c,t);
}
/* ------------------------------------------------------------------------- */

DEFINE_SOURCE(Wiebe_Enthalpy, c, t, dS, eqn)
{
  double f = 0.0;
  if (wiebe->mf_n > 1.0e-7)
    f = C_YI(c,t,udf.iFuel)*C_R(c,t)*C_VOLUME(c,t)/wiebe->mf_n;
/*   return f*udf.LHV*udf.m_fuel/C_VOLUME(c,t)*wiebe->src*wiebe->aFuel; */
  return f*wiebe->aFuel*udf.LHV*wiebe->src/C_VOLUME(c,t);
}
/* ------------------------------------------------------------------------- */

