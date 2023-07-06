#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "engine.h"
#include "engineUtils.h"
#include "accessMDUC.h"
#include "chemMech.h"
#include "fileIO.h"
#include "utils.h"
#include "parser.h"
#include "wiebe.h"

void defaultEngine(engineStruct *engine)
{
  /* general for error checking */
  strcpy(engine->manifoldFile, "");
  strcpy(engine->monFile,  "engine_monitor");
  strcpy(engine->type, "HCCI");
  engine->CAD_0    =  0.0;         // starting crank angle
  engine->CAD_end  = -1.0; 
  engine->Dt_CAD   = -1.0;
  engine->vol_prof = NULL;
  engine->bore     = -1.0;
  engine->stroke   = -1.0;
  engine->CR       = -1.0;
  engine->Acyl = engine->V_TDC = engine->V_BDC = engine->V_SWP = -1.0;
  engine->mass = engine->mfuel = engine->mair  = engine->megr  = -1.0;
  
  /* general for performance parameters */
  engine->CA50 = -1000;

  /* intake valve */
  strcpy(engine->intake.P_file, "");
  strcpy(engine->intake.liftFile, "");
  strcpy(engine->intake.mdotFile, "");
  engine->intake.open   = engine->intake.close  = -1.0;
  engine->intake.CA_L   = engine->intake.CA_P   = engine->intake.CA_mdot  =NULL;
  engine->intake.L_prof = engine->intake.P_prof = engine->intake.mdotProf =NULL;
  engine->intake.iY = NULL;
  engine->intake.Y  = NULL;

  /* exhaust valve */
  strcpy(engine->exhaust.P_file, "");
  strcpy(engine->exhaust.liftFile, "");
  engine->exhaust.open   = engine->exhaust.close  = -1.0;
  engine->exhaust.CA_L   = engine->exhaust.CA_P   = NULL;
  engine->exhaust.L_prof = engine->exhaust.P_prof = NULL;
  engine->exhaust.iY = NULL;
  engine->exhaust.Y  = NULL;

  /* injection */
  engine->inj = NULL;
/*   engine->inj->iY = NULL; */
/*   engine->inj->t  = engine->inj->mdot = engine->inj->Y = NULL; */

  /* spark */
  strcpy(engine->sparkModel, "wiebe");
  engine->spark  = engine->ignited = engine->fullBurn = 0;
  engine->Dflame = engine->Aflame  = 0.0;
  engine->sparkCA = 0.0;
  engine->sparkD  = 0.001;
  engine->m_prem  = 0.0;
  engine->dVbdt   = 0.0;
  engine->nPremSp = -1;
  engine->premSp  = NULL;
  engine->iPremSp = NULL;
  engine->Yburnt  = engine->Yprem = NULL;
  engine->qdotPrem = engine->qPrem = engine->qPrem_m1 = 0.0;
  engine->nBurnRate = -1;
  engine->burnRateScale = 1.0;
  strcpy(engine->burnRateFile, "");

  engine->htModel = 0;
  engine->woschni.Pscale = 1.0;
  engine->woschni.Pref   = 101325;
  engine->woschni.Tref   = 400.0;
  engine->woschni.C1     = 2.28;
  engine->woschni.C2     = 0.0032;
  engine->woschni.C3     = 0.0;
  engine->hohen.Pscale   = 1.0;
  engine->hohen.Pref     = 101325;
  engine->hohen.Tref     = 400.0;
  engine->hohen.C1       = 0.0;
  engine->hohen.C2       = 0.0;
  engine->hohen.C3       = 1.4;
}
/* -------------------------------------------------------------------------- */

void writeEngineManifold(double t, engineStruct *engine)
{
  fprintf(engine->pMan, "%6.2f %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e "
          "%15.6e\n", time_to_cad(t,engine),engine->intake.P, engine->intake.L,
          engine->intake.A_f, engine->intake.mdot, engine->intake.Pr,
          engine->intake.gamma, engine->intake.rho);

  fflush(engine->pMan);
}
/* -------------------------------------------------------------------------- */

void writeEngineMonitorHeader(ioStruct *io, engineStruct *engine)
{
  int i = 1;
  char monFile[256];
  if (strlen(io->dataDir) > 0)
    sprintf(monFile, "%s/%s", io->dataDir, engine->monFile);
  else
    strcpy(monFile, engine->monFile);
  engine->monitor = fopen(monFile, "w");
  if (!engine->monitor) {
    fprintf(stderr, "Error: cannot open engine monitor file %s\n", monFile);
    exit(2);
  }

  fprintf(engine->monitor, "# time[s]<%i> ", i); i++;
  fprintf(engine->monitor, "CAD[deg]<%i> ", i); i++;
  fprintf(engine->monitor, "m_cyl[kg]<%i> ", i); i++;
  fprintf(engine->monitor, "m_air[kg]<%i> ", i); i++;
  fprintf(engine->monitor, "m_fuel[kg]<%i> ", i); i++;
  fprintf(engine->monitor, "m_egr[kg]<%i> ", i); i++;
  fprintf(engine->monitor, "m_burned[kg]<%i> ", i); i++;
  fprintf(engine->monitor, "wiebe_src[1/s]<%i> ", i); i++;
  fprintf(engine->monitor, "rho_b[kg/m^3]<%i> ", i); i++;
  fprintf(engine->monitor, "T_b[K]<%i> ", i); i++;
  fprintf(engine->monitor, "vol_b[m^3]<%i> ", i); i++;
  fprintf(engine->monitor, "cp_b[J/(kg-K)]<%i> ", i); i++;
  fprintf(engine->monitor, "h_b[J/kg]<%i> ", i); i++;
  fprintf(engine->monitor, "h_u[J/kg]<%i> ", i); i++;
  fprintf(engine->monitor, "Yfuel[-]<%i> ", i); i++;

  fprintf(engine->monitor, "\n");
  fflush(engine->monitor);
}
/* -------------------------------------------------------------------------- */

void writeEngineMonitor(double t, ioStruct *io, chemStruct *chem, 
                        engineStruct *engine, flameStruct *flame)
{
  int i;
  double m_cyl = 0.0;
  for (i=0; i<4; i++) m_cyl += flame->Y[flame->sM+i];
  fprintf(engine->monitor, "%15.6e %15.6e ", t, time_to_cad(t,engine));
  fprintf(engine->monitor, "%15.6e ", m_cyl);
  fprintf(engine->monitor, "%15.6e ", flame->Y[flame->sMa]);
  fprintf(engine->monitor, "%15.6e ", flame->Y[flame->sMf]);
  fprintf(engine->monitor, "%15.6e ", flame->Y[flame->sMe]);
  fprintf(engine->monitor, "%15.6e ", flame->Y[flame->sMb]);
  fprintf(engine->monitor, "%15.6e ", engine->wiebeSI.src);
  fprintf(engine->monitor, "%15.6e ", engine->rhoPrem);
  fprintf(engine->monitor, "%15.6e ", engine->Tprem);
  fprintf(engine->monitor, "%15.6e ", engine->Vprem);
  fprintf(engine->monitor, "%15.6e ", engine->cpPrem);
  fprintf(engine->monitor, "%15.6e ", engine->hB);
  fprintf(engine->monitor, "%15.6e ", engine->hU);
  fprintf(engine->monitor, "%15.6e ", engine->Y_fuel);

  fprintf(engine->monitor, "\n");
  fflush(engine->monitor);
}
/* -------------------------------------------------------------------------- */

void initEngine(ioStruct *io, chemStruct *chem, engineStruct *engine)
{
  int i;
  double PI    = 3.141592653589793;
  double cpmix, wmix, R, ysum, tmp;
  FILE *pF;

  /* compute the geometry */
  engine->Acyl  = 0.25*PI*engine->bore*engine->bore;
  engine->V_SWP = engine->Acyl*engine->stroke;

  if (engine->V_TDC > 0.0) {
    engine->V_BDC = engine->V_TDC + engine->V_SWP;
    engine->CR    = engine->V_BDC/engine->V_TDC;
  }
  else {
    engine->V_TDC = engine->V_SWP/(engine->CR-1);
    engine->V_BDC = engine->CR*engine->V_TDC;
  }

  /* calculate the piston location at TDC */
  engine->zTDC = 4*engine->V_TDC/(PI*engine->bore*engine->bore);

  /* only initialize geometry if not actually engine mode */
  if (strcmp(io->timeMode, "engine") != 0)
    return;

  writeEngineMonitorHeader(io, engine);

  /* compute the volume profile and it's derivative in CAD */
  double dCAD, dt;
  engine->nCAD     = 2880;
  dCAD = 720.0/(engine->nCAD-1); /* assume 1 cycle of 720 degrees */
  dt   = dCAD/(6*engine->RPM);
  engine->CAD_vec  = (double*) calloc(engine->nCAD, sizeof(double));
  engine->vol_prof = (double*) calloc(engine->nCAD, sizeof(double));
  engine->dVdt     = (double*) calloc(engine->nCAD, sizeof(double));
  for (i=0; i<engine->nCAD; i++) { /* volume */
    engine->CAD_vec[i]  = i*dCAD;
    engine->vol_prof[i] = calcCylVol(engine->CAD_vec[i], engine);
  }
  /* volume rate */
  engine->dVdt[0] = engine->dVdt[engine->nCAD-1] = 0.0;
  for (i=1; i<engine->nCAD-1; i++)
    engine->dVdt[i]=(engine->vol_prof[i+1]-engine->vol_prof[i-1])/(2*dt);

  /* write profile to file */
  pF = fopen("engine_profile.dat", "w");
  fprintf(pF, "# CAD[deg]<1> volume[m^3]<2> dVdt[m^3/s]<3> ");
  fprintf(pF, "piston_vel[m/s]<4>\n");
  for (i=0; i<engine->nCAD; i++)
    fprintf(pF, "%6.1f %15.6e %15.6e %15.6e\n",engine->CAD_vec[i],
            engine->vol_prof[i], engine->dVdt[i], 
            calcPistonVelocity(engine->CAD_vec[i],engine));
  fclose(pF);

  /* determine whether spark should be computed */
  if (!strcmp(engine->type, "SI") || !strcmp(engine->type, "SACI") || 
      !strcmp(engine->type, "si") || !strcmp(engine->type, "saci") ) {

    /* assign the timing */
    engine->sparkCA = 720 - engine->sparkAdv;

    /* set up specific engine parameters */
    if (!strcmp(engine->sparkModel, "wiebe")) {
      engine->spark = 1;
      engine->wiebe = &engine->wiebeSI;
      if (chem->args)
        chem->args[0] = engine->wiebe;
      else {
        fprintf(stderr, "Warning: arguments not allocated in for chemistry\n");
        chem->args = (void**) calloc(1, sizeof(void*));
        chem->args[0] = engine->wiebe;
/*         exit(2); */
      }
      engine->wiebeSI.t0 = 0.2516*(720-engine->sparkAdv) - 194.2;
      engine->wiebeSI.dt =  engine->wiebeSI.alpha*engine->wiebeSI.t0
                        + engine->wiebeSI.beta
                        - (-1.8802*0.3*100 + 51.74);
      /* debug */
      engine->wiebeSI.t0 = -engine->sparkAdv;
      engine->wiebeSI.a = engine->wiebeSI.a;
      engine->wiebeSI.m = engine->wiebeSI.m;
      engine->wiebeSI.dt = engine->wiebeSI.alpha*engine->wiebeSI.t0
        + engine->wiebeSI.beta;
      fprintf(stderr, "spark:    %g\n", engine->sparkCA);
      fprintf(stderr, "SOC:      %g\n", engine->wiebeSI.t0);
      fprintf(stderr, "duration: %g\n", engine->wiebeSI.dt);
    }
    else if (!strcmp(engine->sparkModel, "flame-speed"))
      engine->spark = 2;

    /* check to see if species have been initialized */
    if (engine->nPremSp > 0) {
      engine->iPremSp = (int*) calloc(engine->nPremSp, sizeof(int));
      if (engine->Yburnt) {
        

      }
      else {
        /* TODO: compute burned species from equivalence ratio */
      }
    }

    engine->Yprem = (double*) calloc(chem->nSpecies, sizeof(double));
    fprintf(stderr, "Premixed zone species:\n");
    engine->MwPrem = 0.0;
    for (i=0; i<engine->nPremSp; i++) {
      engine->iPremSp[i] = getSpeciesIndex(chem->nSpecies, chem->species,
                                           engine->premSp[i]);
      engine->Yprem[engine->iPremSp[i]] = engine->Yburnt[i];
      engine->MwPrem += engine->Yburnt[i]*chem->W[engine->iPremSp[i]];
      fprintf(stderr, "%18s (%3i): %8.6f\n",engine->premSp[i],
              engine->iPremSp[i], engine->Yburnt[i]);
    } /* end loop over burnt species */
    engine->Rprem = 8314.34/engine->MwPrem;
    fprintf(stderr, "premixed gas constant: %8.4f (%8.4f)", engine->Rprem,
            engine->MwPrem);
  } /* end SI setup (TODO: consider separate function for this if stat.) */
  else if (!strcmp(engine->type, "hcci") || !strcmp(engine->type, "HCCI")) {
    if (!strcmp(chem->mechType, "wiebe")) {
      fprintf(stderr, "Engine: using wiebe function w/ paraemeters\n");
      engine->wiebe = &engine->wiebeHCCI;
      chem->args[0] = engine->wiebe; /* assign into the chemistry struct */
      engine->wiebe->dt = engine->wiebe->alpha*engine->wiebe->t0 + engine->wiebe->beta;
      engine->wiebe->tConv = 6*engine->RPM;
      fprintf(stderr, "        a: %6.2f\n", engine->wiebe->a);
      fprintf(stderr, "        m: %6.2f\n", engine->wiebe->m);
      fprintf(stderr, "      SOC: %6.2f\n", engine->wiebe->t0);
      fprintf(stderr, "    alpha: %6.2f\n", engine->wiebe->alpha);
      fprintf(stderr, "     beta: %6.2f\n", engine->wiebe->beta);
      fprintf(stderr, " duration: %6.2f\n", engine->wiebe->dt);
    }
  }

  engine->pMan = fopen("manifold_info.dat", "w");
  fprintf(engine->pMan, "# CAD<1> intake_pressure[Pa]<2> "
          "exhaust_pressure[Pa]<3> \n");
  fflush(engine->pMan);

  /* set up the injections */
  injStruct *inj;
  if (engine->nInj > 0) {
    engine->inj = (injStruct*) calloc(engine->nInj, sizeof(injStruct));
    for (i=0; i<engine->nInj; i++) {
      inj = &engine->inj[i];
      strcpy(inj[i].injFile, engine->injProf[i]); 
      inj->nMdot = readTabFileVar(engine->injProf[i],1,&inj->t,2,&inj->mdot);
      
      inj->CA0 = engine->inj_t[i]; /* must be in cad */
      inj->dCA = inj->t[inj->nMdot-1]/(6*engine->RPM);
      inj->CA1 = inj->CA0 + inj->dCA;
    } /* end loop over injections */
  }

  /* ----- setup intake ----- */
  valveStruct *v;
  v = &engine->intake;
  /* finalize geometry */
  v->W = 0.5*(v->D_head - v->D_seat);
  v->theta *= PI/180.0;
  /* read in the profiles */
  if (strlen(v->P_file) > 0) {
    v->nP = readTabFileVar(v->P_file,   1, &v->CA_P, 2, &v->P_prof);
  }
  for (i=0; i<v->nP; i++) /* convert to absolute pressure from gauge */
    v->P_prof[i] += 101325;
  if (strlen(v->liftFile) > 0) {
    v->nL = readTabFileVar(v->liftFile, 1, &v->CA_L, 2, &v->L_prof);
    /* find the first instance where the valve starts to open */
    v->open = v->close = -1.0;
    v->L_max = -1e20;
    for (i=0; i<v->nL-1; i++) {
      if (v->L_prof[i] > v->L_max)
        v->L_max = v->L_prof[i];
      if (v->open < 0.0 && v->L_prof[i+1] > 1.0e-5) {
        v->open = v->CA_L[i];
      }
      else if (v->open > 0.0 && v->close < 0.0 && v->L_prof[i+1] < 1.0e-5) {
        v->close = v->CA_L[i+1];
        break;
      }
    }
    engine->IVO = v->open;
    engine->IVC = v->close;
  }
  else {
    v->open  = engine->IVO;
    v->close = engine->IVC;
  }
  if (strlen(v->mdotFile) > 0) {
    v->nMdot = readTabFileVar(v->mdotFile, 1, &v->CA_mdot, 5, &v->mdotProf);
  }

  /* set-up the intake species to be the oxidizer boundary species */
  v->nY = io->nOx;
  v->T  = io->oxT;
  chemComputeThermoData(chem->h, chem->cp, v->T, chem);
  v->iY    = (int*)    calloc(v->nY, sizeof(int));
  v->Y     = (double*) calloc(v->nY, sizeof(double));
  v->maskY = (int*)    calloc(chem->nSpecies, sizeof(int));
  ysum = tmp = cpmix = 0.0;
  for (i=0; i<chem->nSpecies; i++) v->maskY[i] = 0;
  for (i=0; i<v->nY; i++) {
    v->iY[i] = io->oxIndex[i];
    v->Y[i]  = io->oxY[i];
    v->maskY[v->iY[i]] = 1;
    cpmix += v->Y[i]*chem->cp[v->iY[i]];
    tmp  += fmin(fmax(v->Y[i],0.0),1.0)/chem->W[v->iY[i]];
    ysum += fmin(fmax(v->Y[i],0.0),1.0);
  }
  wmix     = ysum/tmp;
  v->R     = 8314.34/wmix;
  v->gamma = cpmix/(cpmix - v->R);
  fprintf(stderr, "inlet gamma: %g\n", v->gamma);

  fprintf(stderr, "air species:\n");
  for (i=0; i<v->nY; i++)
    fprintf(stderr, "%i %s\n", v->iY[i], chem->species[v->iY[i]]); 

  FILE *pFt = fopen("test_pressure.dat", "w");
  fprintf(stderr, "# CAD p_intake\n");
  for (i=0; i<v->nP; i++) {
    tmp = i*1.0;
    v->P = lookupProfile(tmp, v->nP, v->CA_P, v->P_prof);
    fprintf(pFt, "%6.1f %15.6e %6.1f %15.6e\n", tmp, v->P, v->CA_P[i],
            v->P_prof[i]);
  }
  fclose(pFt);

  /* ----- setup exhaust ----- */
  fprintf(stderr, "intake P_file: %s\n", v->P_file);
  v = &engine->exhaust;
  fprintf(stderr, "exhaust P_file: %s\n", v->P_file);
  /* finalize geometry */
  v->W = 0.5*(v->D_head - v->D_seat);
  v->theta *= PI/180.0;
  /* read in the profiles */
  if (strlen(v->P_file) > 0) {
    v->nP = readTabFileVar(v->P_file,   1, &v->CA_P, 3, &v->P_prof);
  }
  for (i=0; i<v->nP; i++) /* convert to absolute pressure from gauge */
    v->P_prof[i] += 101325;
  if (strlen(v->liftFile) > 0) {
    v->nL = readTabFileVar(v->liftFile, 1, &v->CA_L, 2, &v->L_prof);
    /* find the first instance where the valve starts to open */
    v->open = v->close = -1.0;
    v->L_max = -1e20;
    for (i=0; i<v->nL-1; i++) {
      if (v->L_prof[i] > v->L_max)
        v->L_max = v->L_prof[i];
      if (v->open < 0.0 && v->L_prof[i+1] > 1.0e-5) {
        v->open = v->CA_L[i];
      }
      else if (v->open > 0.0 && v->close < 0.0 && v->L_prof[i+1] < 1.0e-5) {
        v->close = v->CA_L[i+1];
        break;
      }
    }
    engine->EVO = v->open;
    engine->EVC = v->close;
  }
  else {
    v->open  = engine->EVO;
    v->close = engine->EVC;
  }

  /* set-up the exhaust species to be the all species, initialized from the
     boundary definition */
  v->nY = chem->nSpecies;
  v->T  = io->fuelT1;
  chemComputeThermoData(chem->h, chem->cp, v->T, chem);
  v->iY    = (int*)    calloc(v->nY, sizeof(int));
  v->Y     = (double*) calloc(v->nY, sizeof(double));
  v->maskY = (int*)    calloc(chem->nSpecies, sizeof(int));
  for (i=0; i<chem->nSpecies; i++) v->maskY[i] = 0;
  ysum = tmp = cpmix = 0.0;
  for (i=0; i<io->nFuel1; i++) {
    v->Y[io->fuelIndex1[i]] = io->fuelY1[i];
  }
  for (i=0; i<v->nY; i++) {
    v->iY[i] = i;
    v->maskY[v->iY[i]] = 1;
    cpmix += v->Y[i]*chem->cp[i];
    tmp   += fmin(fmax(v->Y[i],0.0),1.0)/chem->W[i];
    ysum  += fmin(fmax(v->Y[i],0.0),1.0);
  }
  wmix = ysum/tmp;
  v->R = 8314.34/wmix;
  v->gamma = cpmix/(cpmix - v->R);

  if (strlen(engine->burnRateFile) > 0) {
    engine->nBurnRate =
      readTabFileVar(engine->burnRateFile, 1, &engine->burnRateCA, 16,
                     &engine->burnFrac);
    engine->burnRate = (double*) calloc(engine->nBurnRate, sizeof(double));

    for (i=0; i<engine->nBurnRate; i++) {
      engine->burnRateCA[i] += 720;
      if (engine->burnFrac[i] > 0.1)
        engine->burnFrac[i] *= 1e-2;
      else
        engine->burnFrac[i] = 0.0;
    }

    /* burn rate */
    engine->burnRate[0] = engine->burnRate[engine->nBurnRate-1] = 0.0;
    for (i=1; i<engine->nBurnRate-1; i++)
      engine->burnRate[i]=(engine->burnFrac[i+1]-engine->burnFrac[i-1])/
        (engine->burnRateCA[i+1]-engine->burnRateCA[i-1]);

    pFt = fopen("burn_rate.dat", "w");
    fprintf(pFt, "# CA[deg]<1> burn_frac[perc]<2> burn_rate[]<3>\n");
    for (i=0; i<engine->nBurnRate; i++)
      fprintf(pFt, "%8.2f %15.6e %15.6e\n", engine->burnRateCA[i],
              engine->burnFrac[i], engine->burnRate[i]);
    fclose(pFt);
  }
  
  fprintf(stderr, " ----- Valve Timing Summary ----- \n");
  fprintf(stderr, "   IVO: %6.1f  (%6.1f) CAD\n", engine->intake.open,
          engine->intake.open-720); 
  fprintf(stderr, "   IVC: %6.1f  (%6.1f) CAD\n", engine->intake.close,
          engine->intake.close-720); 
  fprintf(stderr, "   EVO: %6.1f  (%6.1f) CAD\n", engine->exhaust.open, 
          engine->exhaust.open+720); 
  fprintf(stderr, "   EVC: %6.1f  (%6.1f) CAD\n", engine->exhaust.close,
          engine->exhaust.close+720); 
  fprintf(stderr, " --------------------------------\n");

  fprintf(stderr, " Valve information \n");
  fprintf(stderr, " %10s | %10s | %10s | \n", "", "intake", "exhaust");
  fprintf(stderr, " %10s | %10.2f | %10.2f | mm\n", "D",
          engine->intake.D_seat*1e3, engine->exhaust.D_seat*1e3);
  fprintf(stderr, " %10s | %10.2f | %10.2f | mm\n", "D_s",
          engine->intake.D_stem*1e3, engine->exhaust.D_stem*1e3);
  fprintf(stderr, " %10s | %10.2f | %10.2f | mm\n", "D_v",
          engine->intake.D_head*1e3, engine->exhaust.D_head*1e3);
  fprintf(stderr, " %10s | %10.2f | %10.2f | mm\n", "D_p",
          engine->intake.D_port*1e3, engine->exhaust.D_port*1e3);
  fprintf(stderr, " %10s | %10.2f | %10.2f | mm\n", "w",
          engine->intake.W*1e3, engine->exhaust.W*1e3);
  fprintf(stderr, " %10s | %10.1f | %10.1f | deg\n", "beta",
          engine->intake.theta*180/PI, engine->exhaust.theta*180/PI);
}
/* -------------------------------------------------------------------------- */

double getEngineMass(char *type, void *mem)
{
  // dereference internal flamelet solver structs
  mducStruct   *mduc   = (mducStruct*) mem;
/*   ioStruct     *io     = mduc->io; */
/*   driveStruct  *driver = mduc->driver; */
/*   chemStruct   *chem   = mduc->chem; */
/*   timeStruct   *time   = mduc->time; */
/*   flowStruct   *flow   = mduc->flow; */
  flameStruct  *flame  = mduc->flame;
  engineStruct *engine = mduc->engine;
  
  if (!strcmp(type, "air"))
    return flame->Y[flame->sMa];
  else if (!strcmp(type, "fuel"))
    return flame->Y[flame->sMf];
  else if (!strcmp(type, "egr"))
    return flame->Y[flame->sMe];
  else if (!strcmp(type, "all")) {
/*     fprintf(stderr, "masses: %15.6e %15.6e %15.6e\n", flame->Y[flame->sM], */
/*             flame->Y[flame->sM+1], flame->Y[flame->sM+2]); */
    return flame->Y[flame->sMa] + flame->Y[flame->sMf] + flame->Y[flame->sMe];
  }
  return -1.0;
}
/* -------------------------------------------------------------------------- */

double calcEnginePressure(double t, engineStruct *engine, flowStruct *flow,
                          flameStruct *flame)
{
  int i,j,k,r,row,col;
  double m, V, rho, P;
  double Tmean, Rmean;
  m = 0.0;
  for (i=0; i<3; i++) 
    m += flame->Y[flame->sM+i];

  V = calcCylVol(time_to_cad(t,engine),engine);
  if (engine->spark && engine->ignited) {
    if (1.0-engine->Vprem/V < 0.01 || engine->fullBurn) {
      rho = flame->Y[flame->sMb]/V;
      engine->fullBurn = 1; /* force burned in future */
      return rho*engine->Rprem*flame->Y[flame->sTb];
    }
    V -= engine->Vprem;
  }

  if (V < 0.0) {
    fprintf(stderr, "Error: unburned volume is negative\n");
    exit(2);
  }

  /* average density based on mass and volume */
  rho = m/(V+1e-12);

  /* compute mean temperature and gas constant */
  Tmean = Rmean = 0.0;
  for (k=0; k<flame->n2; k++) {
    r   = k*flame->n1;
    row = k*flame->n1*flame->nVars;
    for (j=0; j<flame->n1; j++) {
      col = j*flame->nVars;
      Tmean += flow->pdf[r+j]*flame->Y[row+col+flame->sT];
      Rmean += flow->pdf[r+j]*flame->R[r+j];
    }
  }

  return rho*Rmean*Tmean;
}
/* -------------------------------------------------------------------------- */

void initEngineMass(void *mem)
{
  mducStruct   *mduc   = (mducStruct*) mem;
  driveStruct  *driver = mduc->driver;
  ioStruct     *io     = mduc->io;
  chemStruct   *chem   = mduc->chem;
  timeStruct   *time   = mduc->time;
  flowStruct   *flow   = mduc->flow;
  flameStruct  *flame  = mduc->flame;
  engineStruct *engine = mduc->engine;

  int i,j,k,r,row,col,pdfBound[2];

  /* make sure pdf is initialized */
  if (flame->nDim == 0)
    flow->pdf[0] = 1.0;
  else if (flame->nDim == 1)
    MDUCcalcPDF_1D(io->pdfMode,flame->n1,flame->x1,flow->phiMean[1],
                   flow->phiVar[1],flow->pdf,pdfBound);

  /* ensure the local gas properties are up-to-date */
  double Rmean, Tmean;
  Rmean = Tmean = 0.0;
  for (k=0; k<flame->n2; k++) {
    r   = k*flame->n1;
    row = k*flame->n1*flame->nVars;
    for (j=0; j<flame->n1; j++) {
      col = j*flame->nVars;
      Tmean += flow->pdf[r+j]*flame->Y[row+col+flame->sT];
      flame->Mw[r+j] = calcMixtureW(chem->nSpecies,&flame->Y[row+col],chem->W);
      flame->R[r+j]  = flame->RGAS/flame->Mw[r+j];
      Rmean += flow->pdf[r+j]*flame->R[r+j];
    }
  }
  
  /* compute total mass assuming that pressure and CAD_0 have been set */
  double V0 = calcCylVol(time_to_cad(time->t, engine), engine);
  double m0 = flow->P*V0/(Rmean*Tmean);

  /* assign the masses assuming mass fractions have been initialized */
  engine->mair0 = engine->mair  = (1.0-flow->phiMean[1]-flow->phiMean[2])*m0;
  engine->mf0   = engine->mfuel = flow->phiMean[1]*m0;
  engine->megr0 = engine->megr  = flow->phiMean[2]*m0;
  engine->m0    = engine->mass  = engine->mair + engine->mfuel + engine->megr;
  flame->Y[flame->sMa]    = engine->mair;
  flame->Y[flame->sMf]    = engine->mfuel;
  flame->Y[flame->sMe]    = engine->megr;
  flame->Y_m1[flame->sMa] = flame->Y[flame->sMa];
  flame->Y_m1[flame->sMf] = flame->Y[flame->sMf];
  flame->Y_m1[flame->sMe] = flame->Y[flame->sMe];

  fprintf(stderr, "initial masses:\n");
  fprintf(stderr, " air:   %g\n", flame->Y[flame->sMa]);
  fprintf(stderr, " fuel:  %g\n", flame->Y[flame->sMf]);
  fprintf(stderr, " egr:   %g\n", flame->Y[flame->sMe]);
  fprintf(stderr, " total: %g\n", engine->mass);

  if (!strcmp(chem->mechType, "wiebe")) {
    /* initialize the number of fuel moles in the system */
    engine->wiebe->mol0 = engine->wiebe->molF = 0.0;
    for (i=0; i<chem->nSpecies; i++)
      engine->wiebe->mol0 += flame->Y[i]*engine->mass*chem->W[i];
    for (i=0; i<chem->nFuel; i++) 
      engine->wiebe->molF +=
        flame->Y[chem->sFuel[i]]*engine->mass*chem->W[chem->sFuel[i]];
    
/*     engine->wiebe->mol0 = engine->mfuel*chem->fuelMw; */
    fprintf(stderr, " total moles for Wiebe: %g\n", engine->wiebe->mol0);
    fprintf(stderr, " fuel moles for Wiebe: %g\n", engine->wiebe->molF);
  }

  if (engine->spark) {
    flame->Y[flame->sMb] = flame->Y_m1[flame->sMb] = 0.0;
    flame->Y[flame->sTb] = flame->Y_m1[flame->sTb] = flame->Y[flame->sT];
    engine->Tprem = flame->Y[flame->sTb];
    flame->Y[flame->sDF] = flame->Y_m1[flame->sDF] = 0.0;
    engine->rhoPrem = flow->P/(engine->Rprem*flame->Y[flame->sTb]);
    fprintf(stderr, "burnt zone initialization:\n");
    fprintf(stderr, " mass:        %g\n", flame->Y[flame->sMb]);
    fprintf(stderr, " pressure:    %g\n", flow->P);
    fprintf(stderr, " temperature: %g\n", flame->Y[flame->sTb]);
    fprintf(stderr, " density:     %g\n", engine->rhoPrem);
    fprintf(stderr, " flame:       %g\n", flame->Y[flame->sDF]);
  }
}
/* -------------------------------------------------------------------------- */

void calcEngineRHS_0D(double t, timeStruct *time, chemStruct *chem,
                      flowStruct *flow, flameStruct *flame,engineStruct *engine,
                      double *Y, double *RHS)
{
  engineStruct *e = engine;
  valveStruct  *v;
  injStruct    *inj;
  int    i,j,k;
  double CAD, CAD_rel;
  double t_inj, A_f;
  double m_fuel, m_air, m_egr, m_tot, m_prem;
  double mdot_fuel, mdot_air, mdot_egr, mdot;

  /* zero out the mass flow rates */
  for (i=0; i<3; i++) RHS[flame->sM+i] = 0.0;

  /* get the current crank-angle (absolute and relative) */
  CAD     = time_to_cad(t, engine);
  CAD_rel = fmod(CAD, 720);
  /* fprintf(stderr, "CAD_rel: %g\n", CAD_rel); */
  /* ensure that relative crank angle is defined centred around combution */
/*   if (CAD_rel > 360.0) */
/*     CAD_rel -= 720.0; */

  /* retrieve the masses from the solution array */
  m_air  = Y[flame->sMa];
  m_fuel = Y[flame->sMf];
  m_egr  = Y[flame->sMe];
  m_tot  = m_air + m_fuel + m_egr;

  mdot_fuel = mdot_air = mdot_egr = 0.0;

  /* mass flow rate from injection(s) */
  for (i=0; i<e->nInj; i++) {
    inj = &e->inj[i];
    if (CAD_rel > inj->CA0 && CAD_rel < inj->CA1) {
      /* get the current injection position
         NOTE: assumes injection profile is in seconds starting at 0 */
      t_inj = (CAD - inj->CA0)/(6*engine->RPM);
      mdot_fuel += lookupProfile(t_inj, inj->nMdot, inj->t, inj->mdot);
    }
  } /* end loop over injections */

  /* intake mass flow rate */
  if (CAD_rel >= e->IVO && CAD_rel <= e->IVC) {
    /* fprintf(stderr, "computing mass flow rate intake\n"); */
    v = &e->intake; /* reference the intake valve */
    /* lookup the current pressure */
    v->P = lookupProfile(CAD_rel, v->nP, v->CA_P, v->P_prof);

    /* lookup the current lift and compute valve area */
    v->L = lookupProfile(CAD_rel, v->nL, v->CA_L, v->L_prof);
    /* fprintf(stderr, "%g %g %g\n", CAD_rel, v->P, v->L); */
    A_f  = calcValveArea(v->L, v); 

    /* compute the current inlet density */
    v->rho = v->P/(v->R*v->T);
    if (v->CA_mdot != NULL) {
      mdot_air = lookupProfile(CAD_rel, v->nMdot, v->CA_mdot, v->mdotProf);
    }
    else if (v->P/flow->P >= 1.0) { /* inflow */
      mdot_air = calcValveSteadyMassFlow(v->gamma, v->rho, v->P, flow->P, A_f);
    }
    else { /* outflow */
      mdot_air = calcValveSteadyMassFlow(flame->gamma[0], flame->rho[0],
                                         flow->P, v->P, A_f);
      mdot_air *= -1.0;
    }
    if (mdot_air < 0.0 && m_air <= 0.0) {
      mdot_egr = mdot_air;
      mdot_air = 0.0;
    }
    v->A_f = A_f;
    v->mdot = mdot_air;
    v->Pr = flow->P/v->P;
  }

  /* exhaust mass flow rate */
  if (CAD_rel > e->EVO && CAD_rel < e->EVO) {
    v = &e->exhaust; /* reference the intake valve */
    /* lookup the current pressure */
    v->P = lookupProfile(CAD_rel, v->nP, v->CA_P, v->P_prof);

    /* lookup the current lift and compute valve area */
    v->L = lookupProfile(CAD_rel, v->nL, v->CA_L, v->L_prof);
    A_f  = calcValveArea(v->L, v); 

    /* compute the current cylinder density */
    v->rho   = flame->rho[0];
    v->gamma = flame->gamma[0];
    mdot_egr = calcValveSteadyMassFlow(v->gamma, v->rho, flow->P, v->P, A_f);
  }

  /* total mass flow */
  mdot             = mdot_fuel+mdot_air+mdot_egr;
  /* fprintf(stderr, "mdot in: %15.6e\n", mdot_air); */
  RHS[flame->sMa] = mdot_air;
  RHS[flame->sMf] = mdot_fuel;
  RHS[flame->sMe] = mdot_egr;
  /* fprintf(stderr, "%g %g %g\n", mdot_air, mdot_fuel, mdot_egr); */

  /* add appropriately to the */
  double invM  = 1.0/m_tot;
  double tmp   = invM*mdot;
  double srcY  = 0.0;
  double srcH  = 0.0;
  double cpmix = 0.0;

  if (mdot != 0.0) {
/*     fprintf(stderr, "fuel: %g\n", mdot_fuel); */
/*     inj = &e->inj[0]; */

/*     /\* compute the enthalpy at the fuel temperature *\/ */
/*     /\* TODO: make this the surface temperature, account for droplet heatup *\/ */
/*     chemComputeThermoData(chem->h, chem->cp, inj->T, chem); */

/*     /\* compute the fuel species source term, add up enthalpy source *\/ */
/*     srcH = cpmix = 0.0; */
/*     for (i=0; i<inj->nY; i++) { */
/*       srcY   = inj->Y[i]*invM*(mdot_fuel+m_fuel*tmp); */
/*       cpmix += inj->Y[i]*chem->cp[inj->iY[i]]; */
/*       srcH  += srcY*chem->h[inj->iY[i]]; */
/*       /\* RHS[inj->iY[i]] += srcY; *\/ */
/*     } */
/*     RHS[flame->sT] += srcH/cpmix; */

    v = &e->intake;
    /* compute the enthalpy at the air temperature */
    chemComputeThermoData(chem->h, chem->cp, v->T, chem);

    /* compute the air species source term, add up enthalpy source */
    double sumY = 0.0;
    double dMdt = invM*(mdot_air+m_air*tmp);
    sumY = srcH = cpmix = 0.0;
/*     for (i=0; i<v->nY; i++) { */
/*       sumY  += Y[v->iY[i]]; */
/*       srcY   = v->Y[i]*dMdt; */
/*       cpmix += v->Y[i]*chem->cp[v->iY[i]]; */
/*       srcH  += srcY*chem->h[v->iY[i]]; */
/*       RHS[v->iY[i]] += srcY; */
/*     } */
/*     for (i=0; i<chem->nSpecies; i++) { */
/*       if (v->maskY[i] < 1) { */
/*         RHS[i] -= Y[i]*(1.0-sumY)*dMdt; */
/*       } */
/*     } */
/*     RHS[flame->sT] += srcH/cpmix; */

    if (mdot_air < 0.0) mdot_air = 0.0;

    /* use the dY/dt formulation */
    srcH = cpmix = 0.0;
    for (i=0; i<v->nY; i++) {
      srcY  += v->Y[i]*invM*mdot_air;
      cpmix += v->Y[i]*chem->cp[v->iY[i]];
      srcH  += srcY*chem->h[v->iY[i]];
      RHS[v->iY[i]] += srcY;
    }
    RHS[flame->sT] += srcH/cpmix;
    /* fprintf(stderr, "air: %g %g %g\n", mdot_air*1e3, srcY, srcH); */

/*     fprintf(stderr, "egr: %g\n", mdot_egr); */
    v = &e->exhaust;
    /* compute the enthalpy at the cylinder temperature */
    chemComputeThermoData(chem->h, chem->cp, flame->Y[flame->sT], chem);

 /*    srcH = cpmix = 0.0; */
/*     for (i=0; i<v->nY; i++) { */
/*       srcY  += v->Y[i]*invM*mdot_egr; */
/*       cpmix += v->Y[i]*chem->cp[v->iY[i]]; */
/*       srcH  += srcY*chem->h[v->iY[i]]; */
/*       RHS[v->iY[i]] += srcY; */
/*     } */
/*     RHS[flame->sT] += srcH/cpmix; */
/*     fprintf(stderr, "air: %g %g %g\n", mdot_egr*1e3, srcY, srcH); */

    /* account for all Y changes */
    for (i=0; i<chem->nSpecies; i++)
      RHS[i] -= invM*Y[i]*mdot;

 /*    /\* compute the egr species source term, add up enthalpy source *\/ */
/*     dMdt = invM*(mdot_egr+m_egr*tmp); */
/*     srcH = cpmix = sumY = 0.0; */
/*     for (i=0; i<v->nY; i++) { */
/*       sumY  += Y[v->iY[i]]; */
/*       srcY   = v->Y[i]*dMdt; */
/*       cpmix += v->Y[i]*chem->cp[v->iY[i]]; */
/*       srcH  += srcY*chem->h[v->iY[i]]; */
/*       RHS[v->iY[i]] += srcY; */
/*     } */
/*     for (i=0; i<chem->nSpecies; i++) { */
/*       if (v->maskY[i] < 1) { */
/*         RHS[i] -= Y[i]*(1.0-sumY)*dMdt; */
/*       } */
/*     } */
/*     RHS[flame->sT] += srcH/cpmix; */
/*     fprintf(stderr, "air: %g %g %g\n", mdot_egr*1e3, srcY, srcH/cpmix); */
  }

  /* ----- update the pressure ----- */
  /* double vol = calcCylVol(CAD,engine); */
/*   /\* fprintf(stderr, "m_tot: %g, vol: %g\n", m_tot, vol); *\/ */
/*   flow->P = m_tot/vol*flame->R[0]*Y[flame->sT]; */
  
/*   if (isnan(flow->P)) { */
/*     fprintf(stderr, "volume:  %g\n", vol); */
/*     fprintf(stderr, "R:       %g\n", flame->R[0]); */
/*     fprintf(stderr, "T:       %g\n", Y[flame->sT]); */
/*     fprintf(stderr, "m_air:   %g\n", Y[flame->sM]); */
/*     fprintf(stderr, "m_fuel:  %g\n", Y[flame->sM+1]); */
/*     fprintf(stderr, "m_egr:   %g\n", Y[flame->sM+2]); */
/*     fprintf(stderr, "m_tot:   %g\n", m_tot); */
/*     fprintf(stderr, "flow->P: %g\n", flow->P); */
/*     exit(2); */
/*   } */
}
/* -------------------------------------------------------------------------- */

void calcEngineSpark(double t, driveStruct *driver, timeStruct *time,
                     chemStruct *chem, flowStruct *flow, flameStruct *flame,
                     engineStruct *engine, double *Y, double *RHS)
{
  /* this function needs to return a dm/dt term for the burnt region and a
     volume rate of change for the unburnt zone for the dPdt */
  double PI = 3.141592653589793;
  int i,j,k,r,iCAD;
  double CAD     = time_to_cad(t, engine);
  double CAD_rel = fmod(CAD, 720);
  double CAD_wiebe, src;
  double P, dPdt, V;
  double dmbdt = 0.0, dVbdt = 0.0, dTbdt = 0.0;

  /* zero out derivative terms */
  RHS[flame->sMb] = RHS[flame->sTb] = RHS[flame->sDF] = 0.0;
  V = calcCylVol(CAD, engine);
  engine->qdotPrem = 0.0;

  if (driver->PRES_MODE == 3) {
    P    = Y[flame->sP];
    dPdt = RHS[flame->sP];
  }
  else {
    P    = flow->P;
    dPdt = flow->Dt_P;
  }

  /* initialize everything needed in subsequent steps to zero */
  engine->dVbdt = 0.0;

  double m_unb = 0.0, m_b = 0.0, m_cyl = 0.0, Y_fuel = 0.0;
  for (i=0; i<3; i++) m_unb += Y[flame->sM+i];
  m_b = Y[flame->sMb];
  m_cyl = m_unb+m_b;
  
  /* check to see if there is any fuel left */
  Y_fuel = 0.0;
  for (k=0; k<flame->n2; k++) {
    r = k*flame->n1;
    for (j=0; j<flame->n1; j++) {
      for (i=0; i<chem->nFuel; i++) {
        Y_fuel += flow->pdf[r+j]*Y[(r + j)*flame->nVars + chem->sFuel[i]];
      }
    }
  }
  engine->Y_fuel = Y_fuel;

  /* if here ignited */
  double Dflame = Y[flame->sDF];

  /* ----- get the flame properties ----- */
  double sL,sT;         /* laminar and turbulent flame speeds */
  double equivRatio,xb; /* equivalence ratio and diluent */

  engine->Tprem = Y[flame->sTb];
  if (!engine->ignited) {
    /* if it hasn't actually burned yet, use unburned properties */
    engine->rhoPrem = engine->cpPrem = 0.0;
    for (k=0; k<flame->n2; k++) {
      r   = k*flame->n1;
      for (j=0; j<flame->n1; j++) {
        engine->rhoPrem += flow->pdf[r+j]*flame->rho[r+j];
        engine->cpPrem  += flow->pdf[r+j]*flame->cpmix[r+j];
      }
    }
  }
  else { /* burned mass */
    /* compute the thermo properties of the burned zone */
    chemComputeThermoData(chem->h, chem->cp, Y[flame->sTb], chem);
    engine->cpPrem = engine->hB = 0.0;
    for (i=0; i<chem->nSpecies; i++) {
      engine->cpPrem += engine->Yprem[i]*chem->cp[i];
      engine->hB     += engine->Yprem[i]*chem->h[i];
    }

    /* compute the thermo properties of the unburned zone */
    chemComputeThermoData(chem->h, chem->cp, Y[flame->sT], chem);
    engine->hU = 0.0;
    for (i=0; i<chem->nSpecies; i++)
      engine->hU += Y[i]*chem->h[i];

/*     fprintf(stderr, "rhoPrem: %15.6e (%15.6e %15.6e %15.6e)\n", */
/*             engine->rhoPrem, P, engine->Rprem, Y[flame->sTb]); */
    engine->Vprem  = V;
    engine->Vprem /= 1+( (m_unb*engine->MwPrem*flame->Y[flame->sT])/
                         (Y[flame->sMb]*flame->Mw[0]*Y[flame->sTb]) );
    if (engine->Vprem > 0.0)
      engine->rhoPrem = Y[flame->sMb]/engine->Vprem;
  }

  RHS[flame->sTb] = 1.0/(engine->cpPrem*engine->rhoPrem)*dPdt;

  if (m_unb <= 0.0) {
    engine->Vprem = V;
    return;
  }

  /* check whether any fuel left in the burned zone */
  if (Y_fuel < 1.0e-6) return;

  /* ignite if necessary */
  if (!engine->ignited) {
    if (CAD >= engine->sparkCA) { /* ignite */
      Y[flame->sDF] = engine->sparkD;
      engine->ignited = 1;
      fprintf(stderr, "igniting spark to %6.1f mm\n", engine->sparkD*1e3);
    }
    else
      return;
  }

  /* ----- calculate the mass and volume rates of change ----- */
  engine->Dflame = Dflame;
  engine->Aflame = PI*Dflame*Dflame;
  
  double psi = engine->mf0/(engine->mair0+engine->megr0);

  if (engine->spark == 1) { /* wiebe */
    CAD_wiebe = CAD_rel;
    if (CAD_wiebe > 360)
      CAD_wiebe -= 720;
    if (engine->nBurnRate > 0) {
      iCAD = gridLocate(engine->nBurnRate,engine->burnRateCA,CAD);
      src  = linInterp(CAD,engine->burnRateCA[iCAD],engine->burnRate[iCAD],
                       engine->burnRateCA[iCAD+1], engine->burnRate[iCAD+1]);
      engine->wiebeSI.src = fmax(0.0, src);
    }
    else {
      engine->wiebeSI.src = chemCalcWiebe_dt(CAD_wiebe,engine->wiebeSI); /* 1/CAD */
    }
    engine->wiebeSI.src *= 6*engine->RPM;      /* 1/s */
    dmbdt = engine->mf0*engine->wiebeSI.src;   /* kg/s */

    /* account for total mass flow rate based on fuel/air ratio */
    dmbdt *= (1+1/psi);

    dmbdt *= engine->burnRateScale;

    /* TODO: update this based on equation above */
    dVbdt  = dmbdt/engine->rhoPrem; /* m^3/s */

    /* compute the heat release rate for the premixed zone */
    engine->qdotPrem = -dmbdt*(engine->hU - engine->hB);
    
    /* TODO: add heat loss for premixed zone */
    dTbdt = dmbdt/((Y[flame->sMb]+1e-8)*engine->cpPrem)*(engine->hU-engine->hB);
  }
  else if (engine->spark == 2) { /* flame speed */
    equivRatio = Y[flame->sMf]/Y[flame->sMa]/14.7;
    xb = Y[flame->sMe]/Y[flame->sM];
    sL  = calcFlameSpeedHeywood("gasoline",equivRatio,xb,P,Y[flame->sT]);
    sT = 2.0*sL;
    dVbdt = engine->Aflame*sT;
    dmbdt = flame->rho[0]*dVbdt;
    /* dTbdt = 0.0; /\* TODO: fixed burned temperature equation for flame speed *\/ */
    /* compute the heat release rate for the premixed zone */
    engine->qdotPrem = -dmbdt*(engine->hU - engine->hB);

    dTbdt = dmbdt/((Y[flame->sMb]+1e-8)*engine->cpPrem)*(engine->hU-engine->hB);
    /* engine->Vprem = engine->Aflame*Dflame/6.0; */
  }

  /* adjust the diameter based on the flame speed */
  RHS[flame->sDF] = 2.0*sT; /* multiply 2 for diameter */
  /* if (Y_fuel < 1.0e-6) return; */
  RHS[flame->sTb] += dTbdt;

  /* provide factor to adjust unburned mass proportionally */
  RHS[flame->sMa] -= Y[flame->sMa]/m_unb*dmbdt;
  RHS[flame->sMf] -= Y[flame->sMf]/m_unb*dmbdt;
  RHS[flame->sMe] -= Y[flame->sMe]/m_unb*dmbdt;

  /* rate of increase of burned mass */
  RHS[flame->sMb] = dmbdt;

/*   writeEngineMonitor(t, NULL, chem, engine, flame); */
/*   fprintf(stderr, "%8.1f (%8.1f): %15.6e %15.6e %15.6e %15.6e %15.6e " */
/*           "%15.6e\n", */
/*           CAD, CAD_rel, engine->wiebeSI.src, dmbdt, dVbdt, Dflame, */
/*           engine->Vprem, Y[flame->sTb]); */

  /* set the volume rate of change of the burnt region to be added to the
     overall dVdt term in RHS pressure calculation */
  engine->dVbdt = dVbdt;
}
/* -------------------------------------------------------------------------- */

double engineCalcFlameArea(double t, engineStruct *engine)
{
  double PI = 3.141592653589793;
  double V,z;
  double CAD     = time_to_cad(t, engine);
  double c,c_sq,l_s,alpha,beta,a;
  double Df = engine->Dflame;
  double rf = 0.5*Df;
  double Af;

  /* calculate the engine volume */
  V = calcCylVol(CAD, engine);
  
  /* calculate the piston position */
  z = 4*V/(PI*engine->bore*engine->bore);

  /* calculate the intersection of the flame diameter with the cylinder */
  c_sq = engine->sparkX*engine->sparkX + engine->sparkY*engine->sparkY;
  c    = pow(c_sq,0.5);

  if (c + rf < 0.5*engine->bore) {
    /* half sphere*/
    Af = 2*PI*rf*rf;

    if (rf > z) { /* flame intersects piston */
      Af -= PI*Df*(Df-z);
    }
  }
  else { /* flame base intersects cylinder wall */
    /* distance to circle intersection line */
    l_s  = (c_sq - rf*rf + 0.25*engine->bore*engine->bore)/(2*c);
    alpha = 2*acos((l_s-c)/rf);

    /* check if it intersects the piston */
    if (rf > z) {
      /* calculate flame/piston intersection radius */
      a = pow(rf*rf - z*z, 0.5);
      if (c + a < 0.5*engine->bore) {
        Af  = 2*rf*rf*(2*PI-alpha);
        Af -= PI*Df*(rf-z);
      }
      else {
        beta = 0.5*PI + acos(a/rf);
        Af = rf*rf*(2*PI-alpha)*(1-cos(beta));
      }
    }
    else {
      Af = 2*rf*rf*(2*PI-alpha);
    }
  }

  return Af;
}
/* -------------------------------------------------------------------------- */

void destroyEngine(engineStruct *engine)
{
  int i;
  if (engine->vol_prof) free(engine->vol_prof);
  if (engine->dVdt) free(engine->dVdt);
  if (engine->CAD_vec) free(engine->CAD_vec);
  if (engine->iPremSp) free(engine->iPremSp);
  if (engine->premSp) {
    for (i=0; i<engine->nPremSp; i++) free(engine->premSp[i]);
    free(engine->premSp);
  }
  if (engine->Yprem) free(engine->Yprem);
  if (engine->Yburnt) free(engine->Yburnt);
  if (engine->burnRate) free(engine->burnRate);
  if (engine->burnFrac) free(engine->burnFrac);
  if (engine->burnRateCA) free(engine->burnRateCA);
  if (engine->inj_t) free(engine->inj_t);
  if (engine->injDur) free(engine->injDur);
  if (engine->injProf) {
    for (i=0; i<engine->nInj; i++) free(engine->injProf[i]);
    free(engine->injProf);
  }
  if (engine->inj) destroyInjection(engine->inj);
}
/* -------------------------------------------------------------------------- */

void destroyValve(valveStruct *valve)
{
  if (valve->CA_P) free(valve->CA_P);
  if (valve->CA_L) free(valve->CA_L);
  if (valve->L_prof) free(valve->L_prof);
  if (valve->P_prof) free(valve->P_prof);
  if (valve->iY) free(valve->iY);
  if (valve->Y) free(valve->Y);
  if (valve->maskY) free(valve->maskY);
  if (valve->CA_mdot) free(valve->CA_mdot);
  if (valve->mdotProf) free(valve->mdotProf);
}
/* -------------------------------------------------------------------------- */

void destroyInjection(injStruct *inj)
{
  if (inj->t) free(inj->t);
  if (inj->mdot) free(inj->mdot);
  if (inj->iY) free(inj->iY);
  if (inj->Y) free(inj->Y);
 }
/* -------------------------------------------------------------------------- */
