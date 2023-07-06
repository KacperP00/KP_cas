#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "mducMain.h"
#include "postProcess.h"
#include "accessMDUC.h"
#include "chemMech.h"
#include "scalarDiss.h"
#include "parser.h"
#include "errFunc.h"
#include "utils.h"
#include "display.h"
#include "engine.h"
#include "engineUtils.h"
#include "accessChem.h"

void postProcess1D(char *postFile, ioStruct *io, driveStruct *driver,
                   void *memMDUC)
{
  int i, j, k, count;
  
  /* variables for tracking programme run time */
#ifdef ACCESS_MPI
  double startRun, endRun;
#else
  time_t startRun, endRun;
#endif /* ACCESS_MPI */
  double elapsedTime;
  
  // time variables for CFD
  int    nt, ntSteps, *ntIn;
  double *timeIn, curTime, CAD, dt, timeStart, timeEnd;
  
  /* ----- allocate flamelet variables for use in CFD ----- */
  double **Y, *Temp;  // storage for mass fraction and temperature
  int     nPhi;       // number of grid point in flamelet solution
  double *phi;        // flamelet grid
  int    nSpecies;    // number of species

  /* ----- local flow variables ----- */
  int    *pdfBound, pdfFlag = -1, pdfDefined = 0;
  double *chi, *pdf;
  double meanRho, meanRho_0;   // density
  double Tmean;                // mean temperature in domain
  double Tmax;                 // track maximum temperature in domain
  double mass, mass_0;         // system mass
  double vol, vol_0, vol_pre;  // system volume
  double volChem, volPrem;
  double dVdt;
  double gamma;                // polytropic coefficient
  double Rmean;                // mean gas constant

  /* ----- declare flow variables from post-processing ----- */
  double P, *PIn, dPdt, P_pre; // pressure
  double **TbIn;               // time varying boundary temperature
  double **Zmean, **Zvar;      // mixture fraction mean and variance
  double **chiRef;             // reference scalar dissipation rate
  double **ZRef;               // reference mixture fraction
  double strain;
  double **xmin, **xmax;       // min and max of flamelet points
  double tmp;
  /* heat release rates */
  double qdot, Q_tot;          // total heat release rate
  double qdotChem, Q_chem;     // heat release rate from kinetics
  double qdotPrem, Q_prem;     // heat release rate from premixed burn (option)
  double qdotWall, Q_wall;     // wall heat transfer rate
  double Q_h, Q_int;           // copmuted from volume/pressure change

  /* ----- parallel variables ----- */
  int myid = 0;
#ifdef ACCESS_MPI
  mpiStruct   *mpi;
  getMPIStruct(&mpi, memMDUC);
  MPI_Comm_rank(mpi->world, &myid); // get processor ID
  if (myid == 0)
    fprintf(stderr, "Post-processing file: %s\n", postFile);

#endif /* ACCESS_MPI */

  /* get the engine struct just in case */
  engineStruct *engine = MDUCgetEngineStruct(memMDUC);

  /* write the output to a file rather than stderr */
  FILE *oFile;
  char outFileName[256];
  if (strlen(io->dataDir) > 0)
    sprintf(outFileName, "%s/mduc.out", io->dataDir);
  else
    sprintf(outFileName, "mduc.out");
  oFile = fopen(outFileName, "w");
  if (oFile == NULL) fprintf(stderr, "Warning: could not open monitor file\n");

  // set all flow arrays to NULL pointers
  timeIn = PIn = chi = NULL;
  TbIn = chiRef = ZRef = Zmean = Zvar = xmin = xmax = NULL;

  // get the flow variable vectors (over time)
  if (driver->TIME_MODE == 0) {
    /* post-process from variables specified in input file */
    ntSteps = inputPostProcess(&ntIn, &timeIn, &PIn, &chiRef, &ZRef,
                               &Zmean, &Zvar, &xmin, &xmax, memMDUC);
  }
  else if (driver->TIME_MODE == 1) {
    /* read flow variables from individual files output from DNS */
    ntSteps = readFileList(io->postFile, io->fileMin, io->fileMax, &ntIn,
                           &timeIn, &PIn, &chiRef, &ZRef, &Zmean, &Zvar,
                           &xmin, &xmax, memMDUC);
  }
  else if (driver->TIME_MODE == 2) {
    /* read flow variables from scalar file */
    if (!strcmp(postFile, "none")) strcpy(postFile, io->postFile);
    ntSteps = readScalarFile(postFile, &ntIn, &timeIn, &PIn, &TbIn, &chiRef,
                             &ZRef, &Zmean, &Zvar, &xmin, &xmax, memMDUC);
  }
  else if (driver->TIME_MODE == 3) { 
    /* post-process from engine variables specified in input file */
    timeStruct   *time   = MDUCgetTimeStruct(memMDUC);
    
    if (engine->Dt_CAD > 0.0)
      time->dt = cad_to_time(engine->CAD_0+engine->Dt_CAD, engine);
    if (engine->CAD_end > 0.0)
      time->t_end = cad_to_time(engine->CAD_end, engine);
    
    if (!strcmp(postFile, "none"))
      ntSteps = inputPostProcess(&ntIn, &timeIn, &PIn, &chiRef, &ZRef,
                                 &Zmean, &Zvar, &xmin, &xmax, memMDUC);
    else
      ntSteps = readScalarFile(postFile, &ntIn, &timeIn, &PIn, &TbIn, &chiRef,
                               &ZRef, &Zmean, &Zvar, &xmin, &xmax, memMDUC);
  }
  
  /* get initial variables */
  timeStart = timeIn[0];
  dt = timeIn[1] - timeIn[0];
  P = P_pre = PIn[0];
  
  /* ----- initialise the combustion solver ----- */
  nSpecies = MDUCgetNspecies(memMDUC);
  nPhi     = MDUCgetGrid(&phi, memMDUC);

  /* allocate flamelet solution variables */
  Y    = (double**) calloc(nPhi, sizeof(double*));
  Temp = (double*)  calloc(nPhi, sizeof(double));
  for (j=0; j<nPhi; j++)
    Y[j] = (double*) calloc(nSpecies, sizeof(double));

  // allocate local flow variables  
  chi      = (double*) calloc(nPhi, sizeof(double));
  pdf      = (double*) calloc(nPhi, sizeof(double));
  pdfBound = (int*)    calloc(2, sizeof(int));

  // get Zmax for scalar dissipation rate calculation
  double chiZ1max;
  chiZ1max = MDUCgetChiZ1max(memMDUC);

  /* ------------------------------------------ */
  /* ----- initialise the flow parameters ----- */
  /* ------------------------------------------ */
  MDUCinitFlow(timeStart, P, xmin[0][1], xmax[0][1], 0.0, 0.0, memMDUC);

  // read in initial chi/pdf from the flow file if needed
  if (driver->CHI_MODE < 0 || driver->PDF_MODE == -1) {
    MDUCreadFlowData(&chi, NULL, NULL, &pdf, NULL, ntIn[0], memMDUC);
    chi[0] = chi[nPhi-1] = 0;
    
    /* set pdf bounds to whole grid to be safe (only changed if pdf computed) */
    pdfBound[0] = 0;
    pdfBound[1] = nPhi-1;
    pdfDefined  = 1;
    
    /* interpolate to ensure continuous distribution */
    if (driver->CHI_INTERP == 1)
      MDUCinterpChiDNS(&chi, NULL, NULL, memMDUC);
    else if (driver->CHI_INTERP == 2)
      MDUCfitChi1D(nPhi, xmin[0][1], xmax[0][1], phi, chi);
    else if (driver->CHI_INTERP == 3)
      MDUClaplaceChi1D(nPhi, xmin[0][1], xmax[0][1], phi, chi);
  }

  /* ----- calculate the pdf if not read in ----- */
  if (driver->PDF_MODE != -1 && Zmean[0][1] >= 0.0 && 
      MDUCgetNdim(memMDUC) > 0) {
    pdfFlag = MDUCcalcPDF_1D(io->pdfMode,nPhi,phi,Zmean[0][1],Zvar[0][1],
                             pdf,pdfBound);
    pdfDefined = 1;
  }
  else if (MDUCgetNdim(memMDUC) == 0) {
    pdf[0] = 1.0;
    pdfBound[0] = 0;
    pdfBound[1] = 1;
    pdfDefined = 1;
  }

  if (pdfDefined)
    MDUCsetPDF(nPhi,1,&pdf,memMDUC);

  /* ----- scalar dissipation rate ----- */
  if (driver->CHI_MODE >= 0) { /* not read from file */
    if (driver->CHI_MODE < 2) {
      /* initialise which scalar dissipation form to use */
      if (driver->CHI_MODE == 0)      calcChiForm = &calcLogForm;
      else if (driver->CHI_MODE == 1) calcChiForm = &calcErfcForm;
          
      /* compute dissipation according to specified form */
      for (i=0; i<nPhi; i++)
        chi[i] = chiRef[0][1]*calcChiForm(phi[i], ZRef[0][1], chiZ1max);
      
      /* strain = MDUCgetStrainRate(1,memMDUC); */
      /* calcChi1DcounterflowDiffFlame(nPhi,phi,0,strain,chi); */
    }
    else {
      /* specify initial chi here */
      strain = MDUCgetStrainRate(1,memMDUC);
      calcChi1DcounterflowDiffFlame(nPhi,phi,0,strain,chi);
      
      /* for (i=0; i<nPhi; i++) chi[i] = 0.0; */
    }
  }
  /* set scalar dissipation rate into pre timestep */
  MDUCsetScalarDiss(nPhi, 1, &chi, NULL, NULL, memMDUC);

  /* compute mean density if necessary */
  if (driver->ENGINE) {
    /* initEngineMass(memMDUC); */
    mass = mass_0 = getEngineMass("all", memMDUC);
    meanRho_0 = MDUCcalcMeanDensity(P,NULL,&pdf,&pdfBound,pdfFlag,memMDUC);
    Tmean = MDUCcalcMeanVar("temp",Zmean[0][1],Zvar[0][1],-1.0,-1.0,-1.0,
                              &pdf,&pdfBound,pdfFlag, memMDUC);
    vol = vol_0 = vol_pre = calcCylVol(time_to_cad(timeStart, engine),engine);
    gamma = MDUCcalcMeanGamma(&pdf, &pdfBound, pdfFlag, &Rmean, memMDUC);
    
  }
  else if (driver->PRES_MODE > 0) {
    if (pdfDefined) {
      /* get mean properties */
      meanRho_0=MDUCcalcMeanDensity(P,NULL,&pdf,&pdfBound,pdfFlag,memMDUC);
      Tmean = MDUCcalcMeanVar("temp",Zmean[0][1],Zvar[0][1],-1.0,-1.0,-1.0,
                              &pdf,&pdfBound,pdfFlag, memMDUC);
    }
    
    /* define initial volume */
    if (driver->PRES_MODE == 1) { /* constant volume case */
      vol = vol_0 = vol_pre = 1.0;
    }
    else { /* compression case */
      vol_0 = vol_pre = calcCylVol(time_to_cad(timeStart, engine),engine);
      gamma = MDUCcalcMeanGamma(&pdf, &pdfBound, pdfFlag, &Rmean, memMDUC);
    }

    /* calculate initial mass to be conserved in closed volume case */
    mass = mass_0 = meanRho_0*vol_0;
  }

  /* write error file header */
  if (driver->PRES_MODE == 1) {
    writeErrorFileHeader(io->errorFile, 1);
    writeErrorFile(io->errorFile,timeStart,P,PIn[0],xmin[0][1],xmax[0][1],
                   0.0,0.0,memMDUC);
  }
  MDUCgetSolution(&Y, &Temp, memMDUC);

  monStruct *mon = monitorInit("test_monitor");
  nt      = 0;
  curTime = timeStart;
  if (myid == 0) {
    for (j=0, Tmax=1.0e-20; j<nPhi; j++) if (Temp[j] > Tmax) Tmax = Temp[j];
    monitorRegister("step","",&nt,0,mon);
    monitorRegister("time","s",&curTime,-1,mon);
    if (driver->PRES_MODE == 2 || driver->ENGINE)
      monitorRegister("CA","deg",&CAD,-1,mon);
    monitorRegister("P","bar",&P,-1,mon);
    if (driver->PRES_MODE == 2 || driver->ENGINE) {
      monitorRegister("T","K",&Tmean,-1,mon);
      monitorRegister("mass","kg",&mass,-1,mon);
      monitorRegister("gamma","-",&gamma,-1,mon);
      monitorRegister("volume","m^3",&vol,-1,mon);
      monitorRegister("qdot","J/s",&qdot,-1,mon);
      monitorRegister("Q","J",&Q_tot,-1,mon);
      monitorRegister("qdotChem","J/s",&qdotChem,-1,mon);
      monitorRegister("Q_chem","J",&Q_chem,-1,mon);
      monitorRegister("qdotPrem","J/s",&qdotPrem,-1,mon);
      monitorRegister("Q_prem","J",&Q_prem,-1,mon);
      monitorRegister("qdotPres","J/s",&Q_h,-1,mon);
      monitorRegister("Q_pres","J",&Q_int,-1,mon);
      monitorRegister("qdotWall","J/s",&qdotWall,-1,mon);
      monitorRegister("Q_wall","J",&Q_wall,-1,mon);
    }
    else {
      monitorRegister("Tmax","K",&Tmax,-1,mon);
    }
    monitorWriteHead(mon);
    monitorWriteLine(mon);


/*  fprintf(oFile, "%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e" */
/*               "%15.6e %15.6e %15.6e %15.6e %15.6e", P, Tmean, mass, gamma, vol, */
/*               qdot, Q_tot, qdotChem, Q_chem, qdotPrem, Q_prem, Q_h, Q_int); */
    /* print output header */
    int iC = 2;
    fprintf(oFile, "%6s<1> %15s<2> ", "step", "time[s]"); iC++;
    if (driver->PRES_MODE == 2 || driver->ENGINE) {
      fprintf(oFile, "%15s<%i> ", "CAD[deg]",      iC); iC++;
      fprintf(oFile, "%15s<%i> ", "P[bar]",        iC); iC++;
      fprintf(oFile, "%15s<%i> ", "T[K]",          iC); iC++;
      fprintf(oFile, "%15s<%i> ", "mass[kg]",      iC); iC++;
      fprintf(oFile, "%15s<%i> ", "gamma[-]",      iC); iC++;
      fprintf(oFile, "%15s<%i> ", "volume[m^3]",   iC); iC++;
      fprintf(oFile, "%15s<%i> ", "qdot[J/s]",     iC); iC++;
      fprintf(oFile, "%15s<%i> ", "Q[J]",          iC); iC++;
      fprintf(oFile, "%15s<%i> ", "qdotChem[J/s]", iC); iC++;
      fprintf(oFile, "%15s<%i> ", "Q_chem[J]",     iC); iC++;
      fprintf(oFile, "%15s<%i> ", "qdotPrem[J/s]", iC); iC++;
      fprintf(oFile, "%15s<%i> ", "Q_prem[J]",     iC); iC++;
      fprintf(oFile, "%15s<%i> ", "qdotPres[J/s]", iC); iC++;
      fprintf(oFile, "%15s<%i> ", "Q_pres[J]",     iC); iC++;
      fprintf(oFile, "%15s<%i> ", "qdotWall[J/s]", iC); iC++;
      fprintf(oFile, "%15s<%i> ", "Q_wall[J]",     iC); iC++;
    }
    else {
      fprintf(oFile, "%15s<%i> ", "P[bar]",       iC); iC++;
      fprintf(oFile, "%15s<%i> ", "T_max[K]",     iC); iC++;
      fprintf(oFile, "%15s<%i> ", "Z_min[-]",     iC); iC++;
      fprintf(oFile, "%15s<%i> ", "Z_max[-]",     iC); iC++;
    }
    fprintf(oFile, "\n");
    fprintf(oFile, "%6i %15.6e ", 0, timeStart);
    if (driver->PRES_MODE == 2 || driver->ENGINE) {
      fprintf(oFile, "%15.6e ", time_to_cad(timeStart, engine));
      fprintf(oFile, "%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e"
              "%15.6e %15.6e %15.6e %15.6e %15.6e ", P, Tmean, mass_0, gamma,
              vol_0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
    else {
      fprintf(oFile,"%15.6e %15.6e %15.6e %15.6e ",P,Tmax,
              xmin[0][1],xmax[0][1]);
    }
    fprintf(oFile, "\n"); fflush(oFile);
  }
    
#ifdef ACCESS_MPI
  startRun = MPI_Wtime();
#else
  time(&startRun); // start recording run time
#endif /* ACCESS_MPI */
  /* -------------------------------------- */
  /* ----- begin time integration CFD ----- */
  /* -------------------------------------- */
  Q_int = Q_chem = Q_prem = Q_tot = Q_wall = 0.0;
  MDUCwriteFlamelet("initial_solution.rif", memMDUC);
  for (nt=1; nt<ntSteps; nt++) {
    // dereference data from input file
    curTime = timeIn[nt];
    dt      = timeIn[nt] - timeIn[nt-1];
    
    // read data from flow file
    if (driver->CHI_MODE < 0 || driver->PDF_MODE == -1) {
      MDUCreadFlowData(&chi, NULL, NULL, &pdf, NULL, ntIn[nt], memMDUC);
      chi[0] = chi[nPhi-1] = 0;
      pdfDefined = 1;
      
      /* set pdf bounds to whole grid to be safe (only changed if pdf 
         computed) */
      pdfBound[0] = 0;
      pdfBound[1] = nPhi-1;
      
      /* interpolate to ensure continuous distribution */
      if (driver->CHI_INTERP == 1)
        MDUCinterpChiDNS(&chi, NULL, NULL, memMDUC);
      else if (driver->CHI_INTERP == 2)
        MDUCfitChi1D(nPhi, xmin[nt][1], xmax[nt][1], phi, chi);
      else if (driver->CHI_INTERP == 3)
        MDUClaplaceChi1D(nPhi,xmin[nt][1],xmax[nt][1],phi,chi);
    }

    /* ----- calculate the pdf if not read in ----- */
    if (driver->PDF_MODE != -1 && Zmean[nt][1] >= 0.0
        && MDUCgetNdim(memMDUC) > 0) {
      pdfFlag = MDUCcalcPDF_1D(io->pdfMode,nPhi,phi,Zmean[nt][1],Zvar[nt][1],
                               pdf,pdfBound);
      pdfDefined = 1;
      MDUCsetPDF(nPhi,1,&pdf,memMDUC);
    }

    /* ----- calculate scalar dissipation rate ----- */
    if (driver->CHI_MODE >= 0 && driver->CHI_MODE < 2) {
      for (i=0; i<nPhi; i++)
        chi[i] = chiRef[nt][1]*calcChiForm(phi[i], ZRef[nt][1], chiZ1max);
      
      /* calcChi1DcounterflowDiffFlame(nPhi,phi,0,strain,chi); */
    }

    if (nt > 1 && TbIn != NULL) {
      if (myid == 0)
        fprintf(stderr, "Setting temperature bounds: %8.1f %8.1f\n",
                TbIn[nt][0], TbIn[nt][1]);
      MDUCsetTempBounds(TbIn[nt][0], TbIn[nt][1], -1.0, memMDUC);
      MDUCresetODE(memMDUC);
    }

    /* ----- main combustion solve ----- */
    MDUCsolve1D(dt, P, xmin[nt][1], xmax[nt][1], chi, Y, Temp, memMDUC);

    /* update pressure if reading from file */
    if (driver->PRES_MODE == -1) {
      // update pressure and rate
      P    = PIn[nt];
      dPdt = (P - P_pre)/dt;
      
      /* update the temperature base on the current dP/dt and
         temperature from last iteration */
      MDUCupdateTempDPDT(P, dPdt, &Temp, memMDUC);
      
      /* copy the final temperature back into solution array and
         reinitialise SUNDIALS */
      MDUCsetSolution(curTime, P, &Y, &Temp, memMDUC);
    }

    /* ----- calculate Favre average quantities ----- */
    if (pdfDefined || driver->PDF_MODE == -1) {
      // recompute pressure if constant volume

      if (driver->ENGINE) {
        P       = MDUCgetPressure(0, memMDUC);
        meanRho = MDUCcalcMeanDensity(-1,&Temp,&pdf,&pdfBound,pdfFlag,memMDUC);
        vol     = calcCylVol(time_to_cad(curTime,engine),engine);
        Tmean   = MDUCcalcMeanVar("temp",Zmean[nt][1],Zvar[nt][1],-1.0,-1.0,
                                  -1.0,&pdf,&pdfBound,pdfFlag,memMDUC);
        gamma   = MDUCcalcMeanGamma(&pdf,&pdfBound,pdfFlag,&Rmean,memMDUC);
        mass    = getEngineMass("all", memMDUC);

        /* compute the heat release rates (all should be J/s) */
        volChem = vol;
        if (engine->spark)
          volChem -= engine->Vprem;
        qdotChem  = MDUCgetHRrate(volChem, memMDUC);
        qdotPrem  = engine->qdotPrem;
        qdotWall  = MDUCgetHeatLoss(memMDUC);
        qdot      = qdotChem + engine->qdotPrem - qdotWall;

        /* compute the heat release rate based on pressure/volume */
        dPdt = (P - P_pre)/dt;
        dVdt = (vol - vol_pre)/(dt);
        Q_h  = gamma/(gamma-1.0)*P*dVdt + 1.0/(gamma-1.0)*vol*dPdt;
        
        /* integrate the heat release rates */
        if (nt > 2) {
          Q_int  += Q_h*dt;
          Q_chem += qdotChem*dt;
/*           Q_prem += qdotPrem*dt; */
          Q_prem  = engine->qPrem;
          Q_wall  = qdotWall*dt;
          Q_tot  += qdot*dt;
        }

        /* MDUCsetSolution(curTime, P, &Y, &Temp, memMDUC); */
        writeEngineManifold(timeIn[nt], engine);
      }
      else if (driver->PRES_MODE > 0) {
        /* iterate to converge the pressure/temperature that gives the
           correct density */
        count   = 0;
        meanRho = MDUCcalcMeanDensity(P,&Temp,&pdf,&pdfBound,pdfFlag,memMDUC);

        /* update volume if compression */
        if (driver->PRES_MODE == 2)
          vol = calcCylVol(time_to_cad(curTime,engine),engine);
        
        /* compute current mass */
        mass = meanRho*vol;
        
        while (fabs(1.0 - mass/mass_0) > 1.0e-12) {
          // update pressure and rate
          P    -= P*(mass/mass_0 - 1.0);
          dPdt  = (P - P_pre)/dt;
          
          /* update the temperature base on the current dP/dt and
             temperature from last iteration */
          MDUCupdateTempDPDT(P, dPdt, &Temp, memMDUC);
          
          /* compute the new mean density from the flamelet at current
             pressure and temperature */
          meanRho = MDUCcalcMeanDensity(P,&Temp,&pdf,&pdfBound,pdfFlag,memMDUC);
          /* update mass */
          mass    = meanRho*vol;
          
          if (++count > 100)
            mducErr("constant volume pressure not converged","postProcess1D");
        }

        /* copy the final temperature back into solution array and
           reinitialise SUNDIALS */
        MDUCsetSolution(curTime, P, &Y, &Temp, memMDUC);
        
        /* compute the mean properties */
        Tmean = MDUCcalcMeanVar("temp",Zmean[nt][1],Zvar[nt][1],-1.0,-1.0,
                                -1.0,&pdf,&pdfBound,pdfFlag,memMDUC);
        meanRho = MDUCcalcMeanDensity(P,&Temp,&pdf,&pdfBound,pdfFlag,memMDUC);
        gamma = MDUCcalcMeanGamma(&pdf,&pdfBound,pdfFlag,&Rmean,memMDUC);
        mass  = meanRho*vol;
        
        dVdt = (vol - vol_pre)/(dt);
        Q_h = gamma/(gamma-1.0)*P*dVdt + 1.0/(gamma-1.0)*vol*dPdt;
        
        if (nt > 2)
          Q_int += Q_h*dt;
      }
      
      // output average scalar quantities
      MDUCwriteScalars(curTime,Zmean[nt],Zvar[nt],&pdf,&pdfBound,pdfFlag,
                       memMDUC);
    }

#ifdef ACCESS_MPI
      /* send correct pressure to all processes */
    MPI_Bcast(&P, 1, MPI_DOUBLE_PRECISION, 0, mpi->world);
    MDUCsetPressure(-1, P, memMDUC);
#endif /* ACCESS_MPI */

//    // Dump the flow file
//    if (io->writeFlow && MDUCgetNdim(memMDUC) > 0)
//      MDUCwriteFlowData(&chi, NULL, &pdf, NULL, memMDUC);
    
    // write error file
    if (driver->PRES_MODE == 1)
      writeErrorFile(io->errorFile,timeIn[nt],P,PIn[nt],xmin[nt][1],
                     xmax[nt][1],0.0,0.0,memMDUC);
    
    // write out the current step
    if (myid == 0) {
      for (j=0, Tmax=1.0e-20; j<nPhi; j++) if (Temp[j] > Tmax) Tmax = Temp[j];
      monitorWriteLine(mon);
      fprintf(oFile, "%6i %15.6e ", nt, curTime);
      if (driver->PRES_MODE == 2 || driver->ENGINE) {
        fprintf(oFile, "%15.6e ", time_to_cad(curTime, engine));
        fprintf(oFile, "%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e"
                "%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e",P,Tmean,mass,
                gamma,vol,qdot,Q_tot,qdotChem,Q_chem,qdotPrem,Q_prem,Q_h,Q_int,
                qdotWall, Q_wall);
      }
      else {
        fprintf(oFile, "%15.6e %15.6e %15.6e %15.6e ",  P, Tmax,
                xmin[nt][1], xmax[nt][1]);
      }
      fprintf(oFile, "\n"); fflush(oFile);
    }

    // store the old pressure
    P_pre = P;
    if (driver->PRES_MODE > 0) 
    vol_pre = vol;
  }
  /* ----- end time integration ----- */
#ifdef ACCESS_MPI
  endRun = MPI_Wtime();
  elapsedTime = endRun - startRun;
#else
  time(&endRun);
  elapsedTime = difftime(endRun, startRun);
#endif /* ACCESS_MPI */
  fclose(oFile);
  if (myid == 0)
    fprintf(stderr, "Total time elapsed: %5.3f\n", elapsedTime);
  
  // free memory allocated by reading file
  if (ntIn)   free(ntIn);
  if (timeIn) free(timeIn);
  if (PIn)    free(PIn);
  if (chiRef) { for (i=0; i<ntSteps; i++) free(chiRef[i]); free(chiRef); }
  if (ZRef)   { for (i=0; i<ntSteps; i++) free(ZRef[i]);   free(ZRef); }
  if (Zmean)  { for (i=0; i<ntSteps; i++) free(Zmean[i]);  free(Zmean); }
  if (Zvar)   { for (i=0; i<ntSteps; i++) free(Zvar[i]);   free(Zvar); } 
  if (xmin)   { for (i=0; i<ntSteps; i++) free(xmin[i]);   free(xmin); }
  if (xmax)   { for (i=0; i<ntSteps; i++) free(xmax[i]);   free(xmax); }

  // free locally allocated memory
  free(chi); free(pdf); free(pdfBound);
  for (j=0; j<nPhi; j++) free(Y[j]);
  free(Y); free(Temp); free(phi);

  destroyMonitor(mon); free(mon);

  // free all flamelet memory
  MDUCdestroy(memMDUC);

#ifdef ACCESS_MPI
  MPI_Finalize();
#endif /* ACCESS_MPI */
}
/* -------------------------------------------------------------------------- */

void postProcess2D(char *postFile, ioStruct *io, driveStruct *driver, 
                   void *memMDUC)
{
  int i, j, k, count, trans;
  
  // variables for tracking programme run time
#ifdef ACCESS_MPI
  double startRun, endRun;
#else
  time_t startRun, endRun;
#endif /* ACCESS_MPI */
  double elapsedTime;
  
  // time variables for CFD
  int    nt, ntSteps, *ntIn;
  double *timeIn, curTime, dt, timeStart, timeEnd;
  
  /* ----- allocate flamelet variables for use in CFD ----- */
  double ***Y, **Temp;       // storage for mass fraction and temperature
  int    *nPhi;              // number of flamelet grid points
  double **phi;              // flamelet grid
  int    nSpecies;           // number of species

  /* ----- local flow variables ----- */
  int    **pdfBound, pdfFlag = -1, pdfDefined = 0;
  double **chi0, **chi1, **chi2, **chi12, **pdf;
  double meanRho, meanRho_0;     // mean density
  double mass, mass_0;           // system mass
  double vol, vol_0;             // system volume
  double Tmean;                  // mean temperature in domain
  double Tmax;                   // maximum temperature in domain  
  double gamma;                  // polytropic coefficient
  double Rmean;
  
  /* ----- declare flow variables from post-processing ----- */
  double P, *PIn, dPdt, P_pre;   // pressure
  double **TbIn;                 // time varying boundary temperature
  double **Zmean, **Zvar;        // mixture fraction mean and variance
  double **chiRef;               // reference scalar dissipation rate
  double **ZRef;                 // reference mixture fraction
  double **xmin, **xmax;         // min and max of flamelet points
  double *Sdot;                  // conditional source term
  int     nZ;
  int dir, dirZ, dir2;
  
  /* ----- parallel variables ----- */
  int myid = 0;
#ifdef ACCESS_MPI
  mpiStruct   *mpi;
  getMPIStruct(&mpi, memMDUC);
  MPI_Comm_rank(mpi->world, &myid); // get processor ID
  if (myid == 0)
    fprintf(stderr, "Post-processing file %s:\n", postFile);
#endif // ACCESS_MPI
  
  /* get the engine struct just in case */
  engineStruct *engine = MDUCgetEngineStruct(memMDUC);

  // set all flow arrays to NULL pointers
  timeIn = PIn = Sdot = NULL;
  TbIn = chiRef = ZRef = Zmean = Zvar = xmin = xmax = NULL;
  chi1 = chi12 = chi2 = pdf = NULL;

  // read in the flow variables from post-process file
  if (driver->TIME_MODE == 0) {
    // post-process from variables specified in input file
    ntSteps = inputPostProcess(&ntIn, &timeIn, &PIn, &chiRef, &ZRef,
                               &Zmean, &Zvar, &xmin, &xmax, memMDUC);
  }
  else if (driver->TIME_MODE == 1) {
    // read flow variables from individual files output from DNS
    ntSteps = readFileList(io->postFile, io->fileMin, io->fileMax, &ntIn,
                           &timeIn, &PIn, &chiRef, &ZRef, &Zmean, &Zvar,
                           &xmin, &xmax, memMDUC);
  }
  else if (driver->TIME_MODE == 2) { // read flow variables from scalar file
    if (!strcmp(postFile, "none")) strcpy(postFile, io->postFile);
    ntSteps = readScalarFile(postFile, &ntIn, &timeIn, &PIn, &TbIn, &chiRef,
                             &ZRef, &Zmean, &Zvar, &xmin, &xmax, memMDUC);
  }
  else if (driver->TIME_MODE == 3) {
    // post-process from engine variables specified in input file
    timeStruct   *time   = MDUCgetTimeStruct(memMDUC);
    
    if (engine->Dt_CAD > 0.0)
      time->dt = cad_to_time(engine->CAD_0+engine->Dt_CAD, engine);
    if (engine->CAD_end > 0.0)
      time->t_end = cad_to_time(engine->CAD_end, engine);
    
    ntSteps = inputPostProcess(&ntIn, &timeIn, &PIn, &chiRef, &ZRef,
                               &Zmean, &Zvar, &xmin, &xmax, memMDUC);
  }
  
  // get initial variables
  timeStart = timeIn[0];
  dt = timeIn[1] - timeIn[0];
  P  = P_pre = PIn[0];
  
  /* ----- initialise the combustion solver ----- */
  nSpecies = MDUCgetNspecies(memMDUC);
  nPhi     = MDUCgetGrid2D(&phi, memMDUC);
  
  // allocate local memory
  Y        = (double***) calloc(nPhi[1], sizeof(double**));
  Temp     = (double**)  calloc(nPhi[1], sizeof(double*));
  chi0     = (double**)  calloc(nPhi[1], sizeof(double*));
  chi1     = (double**)  calloc(nPhi[1], sizeof(double*));
  chi12    = (double**)  calloc(nPhi[1], sizeof(double*));
  chi2     = (double**)  calloc(nPhi[1], sizeof(double*));
  pdf      = (double**)  calloc(nPhi[1], sizeof(double*));
  pdfBound = (int**)     calloc(nPhi[1], sizeof(int*));
  for (k=0; k < nPhi[1]; k++) {
    Y[k]        = (double**) calloc(nPhi[0], sizeof(double*));
    Temp[k]     = (double*)  calloc(nPhi[0], sizeof(double));
    chi0[k]     = (double*)  calloc(nPhi[0], sizeof(double));
    chi1[k]     = (double*)  calloc(nPhi[0], sizeof(double));
    chi12[k]    = (double*)  calloc(nPhi[0], sizeof(double));
    chi2[k]     = (double*)  calloc(nPhi[0], sizeof(double));
    pdf[k]      = (double*)  calloc(nPhi[0], sizeof(double));
    pdfBound[k] = (int*)     calloc(2, sizeof(int));
    for (j=0; j < nPhi[0]; j++)
      Y[k][j] = (double*) calloc(nSpecies, sizeof(double));
  }

  // Find out if H direction. If yes, find grid points in Z
  nZ = 0;
  if (driver->H_DIM == 1) {
    nZ = nPhi[1];
    dirZ = 1;
  }
  if (driver->H_DIM == 2) {
    nZ = nPhi[0];
    dirZ = 0;
  }

  // allocate Sdot
  if (driver->H_DIM != 0)
    Sdot = (double*) calloc(nZ, sizeof(double));
  
  /* ------------------------------------------ */
  /* ----- initialise the flow parameters ----- */
  /* ------------------------------------------ */
  MDUCinitFlow(timeStart,P,xmin[0][1],xmax[0][1],xmin[0][2],xmax[0][2],memMDUC);

  // read in initial data from the flow file
  if (driver->CHI_MODE == -1 || driver->PDF_MODE == -1) {
    MDUCreadFlowData(chi1,chi12,chi2,pdf,Sdot,ntIn[0],memMDUC);
    
    // force chi1 to be zero at boundaries Z1=0 and Z1=1
    for (i=0; i<nPhi[1]; i++)
      chi1[i][0] = chi1[i][nPhi[0]-1] = 0;
    
    // force chi2 to be zero at boundaries Z2=0 and Z2=1
    for (i=0; i<nPhi[0]; i++)
      chi2[0][i] = chi2[nPhi[1]-1][i] = 0;
    
    // interpolate scalar dissipation rate
    if (driver->CHI_INTERP == 1)
      MDUCinterpChiDNS(chi1, chi12, chi2, memMDUC);
    else if (driver->CHI_INTERP == 3)
      MDUCinterpChiLap(chi1, chi12, chi2, memMDUC);
    
    // transform variables if 2D mixture fraction
    if (driver->H_DIM == 0 && driver->TRANS_Z2 > 0) {
      MDUCtransformVar(1,driver->TRANS_Z2,nPhi[0],nPhi[1],phi[0],phi[1],chi1);
      MDUCtransformVar(1,driver->TRANS_Z2,nPhi[0],nPhi[1],phi[0],phi[1],chi2);
      MDUCtransformVar(1,driver->TRANS_Z2,nPhi[0],nPhi[1],phi[0],phi[1],chi12);
      MDUCtransformVar(1,driver->TRANS_Z2,nPhi[0],nPhi[1],phi[0],phi[1],pdf);
    }
    
    /* set pdf bounds to whole grid to be safe (only changed if pdf computed) */
    pdfDefined = 1;
    for (j = 0; j<nPhi[1]; j++) {
      pdfBound[j][0] = 0;
      pdfBound[j][1] = nPhi[0]-1;
    }
  }
  
  /* ----- calculate the pdf if not read in ----- */
  if (driver->PDF_MODE != -1 && Zmean) {
    pdfFlag = MDUCcalcPDF_2D(io->pdfMode,nPhi,phi,Zmean[0][1],
                             Zvar[0][1],Zmean[0][2],Zvar[0][2],
                             Zvar[0][0],pdf,pdfBound,driver->TRANS_Z2);
    pdfDefined = 1;
  } 

  /* ----- scalar dissipation rate ----- */
  if (driver->CHI_MODE != -1) // compute from reference value
    calcChi2D(io->chiMode, nPhi, phi, ZRef[0], chiRef[0], chi1, chi12, chi2,
              driver->TRANS_Z2, driver->CHI_SCALE);
  
  /* set scalar dissipation rate into pre timestep */
  MDUCsetScalarDiss(nPhi[0], nPhi[1], chi1, chi12, chi2, memMDUC);

  /* compute mean density if necessary */
  if (driver->PRES_MODE > 0) {
    if (pdfDefined) {
      /* get mean properties */
      meanRho_0 = MDUCcalcMeanDensity(P,NULL,pdf,pdfBound,pdfFlag,memMDUC);
      Tmean     = MDUCcalcMeanVar("temp",Zmean[0][1],Zvar[0][1],Zmean[0][2],
                                  Zvar[0][2],Zvar[0][0],pdf,pdfBound,pdfFlag,
                                  memMDUC);
    }

    /* define initial volume */
    if (driver->PRES_MODE == 1) { /* constant volume case */
      vol = vol_0 = 1.0;
    }
    else { /* compression case */
      vol_0 = calcCylVol(time_to_cad(timeStart, engine),engine);
      gamma = MDUCcalcMeanGamma(pdf, pdfBound, pdfFlag, &Rmean, memMDUC);
    }
    
    /* calculate initial mass to be conserved in closed volume case */
    mass = mass_0 = meanRho_0*vol_0;
  }

  /* write error file header */
  if (driver->PRES_MODE == 1 && myid == 0) {
    writeErrorFileHeader(io->errorFile, 2);
    writeErrorFile(io->errorFile,timeStart,P,PIn[0],xmin[0][1],xmax[0][1],
                   xmin[0][2],xmax[0][2],memMDUC);
  }

  MDUCgetSolution(Y, Temp, memMDUC);
  if (myid == 0) {
    /* print output header */
    fprintf(stderr, "%6s %15s %15s %15s %15s %15s %15s %15s\n","step","time[s]",
            "P[bar]","Tmax[K]","phi1_min","phi1_max","phi2_min","phi2_max");
    for (k=0, Tmax=1.0e-20; k<nPhi[1]; k++)
      for (j=0; j<nPhi[0]; j++) if (Temp[k][j] > Tmax) Tmax = Temp[k][j];
    fprintf(stderr, "%6i %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n", 0,
            timeStart, P, Tmax, xmin[0][1], xmax[0][1], xmin[0][2], xmax[0][2]);
  }
#ifdef ACCESS_MPI
  startRun = MPI_Wtime();
#else
  time(&startRun); // start recording run time
#endif /* ACCESS_MPI */
  /* -------------------------------------- */
  /* ----- begin time integration CFD ----- */
  /* -------------------------------------- */

  for (nt=1; nt<ntSteps; nt++) {
    // dereference data from input file
    curTime  = timeIn[nt];
    dt       = timeIn[nt] - timeIn[nt-1];
    
    // read data from flow file
    if (driver->CHI_MODE == -1 || driver->PDF_MODE == -1) {
      MDUCreadFlowData(chi1, chi12, chi2, pdf, Sdot, ntIn[nt], memMDUC);
          
      // force chi1 to be zero at boundaries Z1=0 and Z1=1
      for (i=0; i<nPhi[1]; i++)
        chi1[i][0] = chi1[i][nPhi[0]-1] = 0;
      
      // force chi2 to be zero at boundaries Z2=0 and Z2=1
      for (i=0; i<nPhi[0]; i++)
        chi2[0][i] = chi2[nPhi[1]-1][i] = 0;
      
      // interpolate scalar dissipation rate
      if (driver->CHI_INTERP == 1)
        MDUCinterpChiDNS(chi1, chi12, chi2, memMDUC);
      else if (driver->CHI_INTERP == 3)
        MDUCinterpChiLap(chi1, chi12, chi2, memMDUC);
      
      // transform variables if 2D mixture fraction
      if (driver->H_DIM == 0 && driver->TRANS_Z2 > 0) {
        MDUCtransformVar(1,driver->TRANS_Z2,nPhi[0],nPhi[1],phi[0],phi[1],chi1);
        MDUCtransformVar(1,driver->TRANS_Z2,nPhi[0],nPhi[1],phi[0],phi[1],chi2);
        MDUCtransformVar(1,driver->TRANS_Z2,nPhi[0],nPhi[1],phi[0],phi[1],chi12);
        MDUCtransformVar(1,driver->TRANS_Z2,nPhi[0],nPhi[1],phi[0],phi[1],pdf);
      }
    }
    
    /* ----- calculate the pdf if not read in ----- */
    if (driver->PDF_MODE != -1 && Zmean) {
      pdfFlag = MDUCcalcPDF_2D(io->pdfMode,nPhi,phi,Zmean[nt][1],Zvar[nt][1],
                               Zmean[nt][2],Zvar[nt][2],Zvar[nt][0],pdf,
                               pdfBound,driver->TRANS_Z2);
      pdfDefined = 1;
    } 
    
    /* ----- calculate scalar dissipation rate ----- */
    if (driver->CHI_MODE != -1) // compute from reference value
      calcChi2D(io->chiMode,nPhi,phi,ZRef[nt],chiRef[nt],chi1,chi12,chi2,
                driver->TRANS_Z2,driver->CHI_SCALE);

    if (nt > 1 && TbIn != NULL) {
      if (myid == 0)
        fprintf(stderr, "Setting temperature bounds: %8.1f %8.1f %8.1f\n",
                TbIn[nt][0], TbIn[nt][1], TbIn[nt][2]);
      MDUCsetTempBounds(TbIn[nt][0], TbIn[nt][1], TbIn[nt][2], memMDUC);
      MDUCresetODE(memMDUC);
    }
    
    /* ----- main combustion solve ----- */
    MDUCflamelet2D(dt,P,xmin[nt][1],xmax[nt][1],xmin[nt][2],xmax[nt][2],
                   chi1,chi12,chi2,0,Y,Temp,memMDUC);
    
    /* update pressure if reading from file */
    if (myid == 0 && driver->PRES_MODE == -1) {
      // update pressure and rate
      P    = PIn[nt];
      dPdt = (P - P_pre)/dt;
      
      /* update the temperature base on the current dP/dt and
         temperature from last iteration */
      MDUCupdateTempDPDT(P, dPdt, Temp, memMDUC);
      
      /* copy the final temperature back into solution array and
         reinitialise SUNDIALS */
      MDUCsetSolution(curTime, P, Y, Temp, memMDUC);
    }
    
    /* ----- calculate Favre average quantities ----- */
    if (myid == 0 && pdfDefined) { 
      /* recompute pressure if constant volume */
      if (driver->PRES_MODE > 0) {
        /* iterate to converge the pressure/temperature that gives the
           correct density */
        count   = 0;
        meanRho = MDUCcalcMeanDensity(P, Temp, pdf, pdfBound, pdfFlag, memMDUC);
        
        /* update volume if compression */
        if (driver->PRES_MODE == 2)
          vol = calcCylVol(time_to_cad(curTime,engine),engine);
        
        /* compute current mass */
        mass = meanRho*vol;
        
        while (fabs(1.0 - mass/mass_0) > 1.0e-12) {
          /* update pressure and rate */
          P    -= P*(mass/mass_0 - 1.0);
          dPdt  = (P - P_pre)/dt;
          
          /* update the temperature base on the current dP/dt and
             temperature from last iteration */
          MDUCupdateTempDPDT(P, dPdt, Temp, memMDUC);
          
          /* compute the new mean density from the flamelet at current
             pressure and temperature */
          meanRho = MDUCcalcMeanDensity(P,Temp,pdf,pdfBound,pdfFlag,memMDUC);

          /* udpate mass */
          mass = meanRho*vol;
          
          if (++count > 100)
            mducErr("constant volume pressure not converged","postProcess2D");
        }
        
        /* copy the final temperature back into solution array and
           reinitialise SUNDIALS */
        MDUCsetSolution(curTime, P, Y, Temp, memMDUC);
      }
      
      /* output scalars */
      MDUCwriteScalars(curTime,Zmean[nt],Zvar[nt],pdf,pdfBound,pdfFlag,memMDUC);
    }
      
#ifdef ACCESS_MPI
      /* send correct pressure to all processes */
      MPI_Bcast(&P, 1, MPI_DOUBLE_PRECISION, 0, mpi->world);
      MDUCsetPressure(-1, P, memMDUC);
#endif /* ACCESS_MPI */
      
//      /* dump the flow file */
//      if (io->writeFlow && myid == 0)
//        MDUCwriteFlowData(chi1, chi2, pdf, Sdot, memMDUC);
  
      /* write error file */
      if (driver->PRES_MODE == 1 && myid == 0)
        writeErrorFile(io->errorFile,timeIn[nt],P,PIn[nt],xmin[nt][1],
                       xmax[nt][1],xmin[nt][2],xmax[nt][2],memMDUC);

      for (k=0, Tmax=1.0e-20; k<nPhi[1]; k++)
        for (j=0; j<nPhi[0]; j++) if (Temp[k][j] > Tmax) Tmax = Temp[k][j];

      if (myid == 0) {
        fprintf(stderr,"%6i %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n",
                nt,curTime,P,Tmax,xmin[nt][1],xmax[nt][1],
                xmin[nt][2],xmax[nt][2]);
      }

      /* store pressure */
      P_pre = P;
    } /* end time integration */
#ifdef ACCESS_MPI
  endRun = MPI_Wtime();
  elapsedTime = endRun - startRun;
#else
  time(&endRun);
  elapsedTime = difftime(endRun, startRun);
#endif /* ACCESS_MPI */
  if (myid == 0)
    fprintf(stderr, "Total time elapsed: %5.3f\n", elapsedTime);
  
  // free memory allocated by reading file
  if (ntIn)   free(ntIn);
  if (timeIn) free(timeIn);
  if (PIn)    free(PIn);
  for (nt=0; nt<ntSteps; nt++) {
    if (TbIn)   free(TbIn[nt]);
    if (chiRef) free(chiRef[nt]);
    if (ZRef)   free(ZRef[nt]);
    if (Zmean)  free(Zmean[nt]);
    if (Zvar)   free(Zvar[nt]);
    if (xmin)   free(xmin[nt]);
    if (xmax)   free(xmax[nt]);
  }
  if (TbIn)   free(TbIn);
  if (chiRef) free(chiRef);
  if (ZRef)   free(ZRef);
  if (Zmean)  free(Zmean);
  if (Zvar)   free(Zvar);
  if (xmin)   free(xmin);
  if (xmax)   free(xmax);
  
  /* free locally allocated memory */
  for (k=0; k<nPhi[1]; k++) {
    for (j=0; j<nPhi[0]; j++) free(Y[k][j]);
    free(Y[k]); free(Temp[k]);
    free(pdf[k]); free(pdfBound[k]);
    free(chi0[k]); free(chi1[k]); free(chi2[k]); free(chi12[k]);
  }
  free(Y); free(Temp);
  free(pdf); free(pdfBound);
  free(chi0); free(chi1); free(chi2); free(chi12);
  free(nPhi); free(phi[0]); free(phi[1]); free(phi);
  
  /* free all flamelet memory */
  if (myid == 0) fprintf(stderr, "freeing flamelet memory...\n");
  MDUCdestroy(memMDUC);
  
#ifdef ACCESS_MPI
  /* call MPI_Init before handling any other options */
  MPI_Finalize();
#endif /* ACCESS_MPI */
}
/* -------------------------------------------------------------------------- */

void postProcessMulti1D(char *postFile, ioStruct *io, driveStruct *driver, 
                        void *mem, void *mem1, void *mem2)
{
  int i, j, k, count, trans;
  
  // variables for tracking programme run time
#ifdef ACCESS_MPI
  double startRun, endRun;
#else
  time_t startRun, endRun;
#endif /* ACCESS_MPI */
  double elapsedTime;
  
  // time variables for CFD
  int    nt, ntSteps, *ntIn;
  double *timeIn, curTime, dt, timeStart, timeEnd;
  
  /* ----- allocate flamelet variables for use in CFD ----- */
  double ***Y, **Temp;       // storage for mass fraction and temperature
  double ***Y_pre, **T_pre;
  double **H_tot, **HR, Tmax;
  double **Y_1, **Y_2;
  double *T_1, *T_2;
  int    *nPhi;              // number of flamelet grid points
  double **phi;              // flamelet grid
  int    nSpecies;           // number of species
  int    nSolve;

  /* ----- local flow variables ----- */
  int    **pdfBound, pdfFlag = -1, *pdfB_1, *pdfB_2;
  double **chi0, **chi1, **chi2, **chi12, **pdf;
  double *chiM_1, *chiM_2, *pdfM_1, *pdfM_2;
  double Z1, Z2;
  double meanRho, meanRho_0;     // mean density
  
  /* ----- declare flow variables from post-processing ----- */
  double P, *PIn, dPdt, P_pre;   // pressure
  double **TbIn;                 // variable boundary temperature
  double **Zmean, **Zvar;        // mixture fraction mean and variance
  double **chiRef;               // reference scalar dissipation rate
  double **ZRef;                 // reference mixture fraction
  double **xmin, **xmax;         // min and max of flamelet points
  double *Sdot;                  // conditional source term
  double *h_i, *cp_i, *h_pre, *cp_pre, *Mw, wmix, wmix_pre, hr_local;
  double rho_L, rho_L_pre;
  int    nZ;
  int dir, dirZ, dir2;
  double RGAS = 8314.34;
  
  /* ----- parallel variables ----- */
  int myid = 0;
#ifdef ACCESS_MPI
  mpiStruct   *mpi;
  getMPIStruct(&mpi, mem);
  MPI_Comm_rank(mpi->world, &myid); /* get processor ID */
  if (myid == 0)
    fprintf(stderr, "Post-processing file %s:\n", postFile);
#endif /* ACCESS_MPI */
  
  /* set all flow arrays to NULL pointers */
  timeIn = PIn = Sdot = NULL;
  chiRef = ZRef = Zmean = Zvar = xmin = xmax = NULL;
  chi1 = chi12 = chi2 = pdf = NULL;

  /* read in the flow variables from post-process file */
  if (driver->TIME_MODE == 0) {
    /* post-process from variables specified in input file */
    ntSteps = inputPostProcess(&ntIn, &timeIn, &PIn, &chiRef, &ZRef,
                               &Zmean, &Zvar, &xmin, &xmax, mem);
  }
  else if (!strcmp(postFile, "none")) {
    /* read flow variables from individual files output from DNS */
    ntSteps = readFileList(io->postFile, io->fileMin, io->fileMax, &ntIn,
                           &timeIn, &PIn, &chiRef, &ZRef, &Zmean, &Zvar,
                           &xmin, &xmax, mem);
  }
  else { /* read flow variables from scalar file */
    if (!strcmp(postFile, "none")) strcpy(postFile, io->postFile);
    ntSteps = readScalarFile(postFile, &ntIn, &timeIn, &PIn, &TbIn, &chiRef,
                             &ZRef, &Zmean, &Zvar, &xmin, &xmax, mem);
  }

  /* get initial variables */
  timeStart = timeIn[0];
  dt = timeIn[1] - timeIn[0];
  P  = P_pre = PIn[0];
  
  /* ----- initialise the combustion solver ----- */
  nSpecies = MDUCgetNspecies(mem);
  nSolve   = MDUCgetNspecs(mem);
  /* for now, just make sure that */
  nPhi = (int*) calloc(2, sizeof(int));

/*   nPhi[0] = MDUCgetGrid(phi[0], mem1); */
/*   nPhi[1] = MDUCgetGrid(phi[1], mem2); */
  nPhi = MDUCgetGrid2D(&phi, mem);
  int n1 = nPhi[0], n2 = nPhi[1];
  int jL, kB;

  h_i    = (double*) calloc(nSpecies, sizeof(double));
  cp_i   = (double*) calloc(nSpecies, sizeof(double));
  h_pre  = (double*) calloc(nSpecies, sizeof(double));
  cp_pre = (double*) calloc(nSpecies, sizeof(double));
  Mw     = (double*) calloc(nSpecies, sizeof(double));
  MDUCgetMolWeightAll(Mw, mem);
  
  // allocate local memory
  Y        = (double***) calloc(nPhi[1], sizeof(double**));
  Y_pre    = (double***) calloc(nPhi[1], sizeof(double**));
  Temp     = (double**)  calloc(nPhi[1], sizeof(double*));
  T_pre    = (double**)  calloc(nPhi[1], sizeof(double*));
  H_tot    = (double**)  calloc(nPhi[1], sizeof(double*));
  HR       = (double**)  calloc(nPhi[1], sizeof(double*));
  chi0     = (double**)  calloc(nPhi[1], sizeof(double*));
  chi1     = (double**)  calloc(nPhi[1], sizeof(double*));
  chi12    = (double**)  calloc(nPhi[1], sizeof(double*));
  chi2     = (double**)  calloc(nPhi[1], sizeof(double*));
  pdf      = (double**)  calloc(nPhi[1], sizeof(double*));
  pdfBound = (int**)     calloc(nPhi[1], sizeof(int*));
  pdfB_1   = (int*)      calloc(2, sizeof(int));
  pdfB_2   = (int*)      calloc(2, sizeof(int));
  for (k=0; k < nPhi[1]; k++) {
    Y[k]        = (double**) calloc(nPhi[0], sizeof(double*));
    Y_pre[k]    = (double**) calloc(nPhi[0], sizeof(double*));
    Temp[k]     = (double*)  calloc(nPhi[0], sizeof(double));
    T_pre[k]    = (double*)  calloc(nPhi[0], sizeof(double));
    H_tot[k]    = (double*)  calloc(nPhi[0], sizeof(double));
    HR[k]       = (double*)  calloc(nPhi[0], sizeof(double));
    chi0[k]     = (double*)  calloc(nPhi[0], sizeof(double));
    chi1[k]     = (double*)  calloc(nPhi[0], sizeof(double));
    chi12[k]    = (double*)  calloc(nPhi[0], sizeof(double));
    chi2[k]     = (double*)  calloc(nPhi[0], sizeof(double));
    pdf[k]      = (double*)  calloc(nPhi[0], sizeof(double));
    pdfBound[k] = (int*)     calloc(2, sizeof(int));
    for (j=0; j < nPhi[0]; j++) {
      Y[k][j]     = (double*) calloc(nSpecies, sizeof(double));
      Y_pre[k][j] = (double*) calloc(nSpecies, sizeof(double));
    }
  }

  double **Y_1_m1, *T_1_m1, *HR_1, *HR_1_est;

  pdfM_1 = (double*) calloc(nPhi[0], sizeof(double));
  chiM_1 = (double*) calloc(nPhi[0], sizeof(double));
  pdfM_2 = (double*) calloc(nPhi[1], sizeof(double));
  chiM_2 = (double*) calloc(nPhi[1], sizeof(double));
  Y_1    = (double**) calloc(nPhi[0], sizeof(double*));
  Y_1_m1 = (double**) calloc(nPhi[0], sizeof(double*));
  T_1    = (double*) calloc(nPhi[0], sizeof(double));
  T_1_m1 = (double*) calloc(nPhi[0], sizeof(double));
  HR_1     = (double*) calloc(nPhi[0], sizeof(double));
  HR_1_est = (double*) calloc(nPhi[0], sizeof(double));
  Y_2    = (double**) calloc(nPhi[1], sizeof(double*));
  T_2    = (double*) calloc(nPhi[1], sizeof(double));
  double HRavg1, HRavg1_est;

  for (j=0; j<nPhi[0]; j++) {
    Y_1[j] = (double*) calloc(nSpecies, sizeof(double));
    Y_1_m1[j] = (double*) calloc(nSpecies, sizeof(double));
  }
  for (k=0; k<nPhi[1]; k++) Y_2[k] = (double*) calloc(nSpecies, sizeof(double));

  double h_n, h_nm1, rho_n, rho_nm1;

  // retrieve the 1D solutions and make sure that the 2D is up-to-date
  MDUCgetSolution(&Y_1, &T_1, mem1);
  for (j=0; j<n1; j++) {
    MDUCcomputeThermoData(h_i, cp_i, T_1[j], mem);
    for (i=0, H_tot[0][j]=0.0; i<nSpecies; i++) {
      Y[0][j][i]   = Y_1[j][i];
      Temp[0][j]   = T_1[j];
      H_tot[0][j] += Y_1[j][i]*h_i[i];
    }
  }
  MDUCgetSolution(&Y_2, &T_2, mem2);
  for (k=0; k<n2; k++) {
    MDUCcomputeThermoData(h_i, cp_i, T_2[k], mem2);
    for (i=0, H_tot[k][0]=0.0; i<nSpecies; i++) {
      Y[k][0][i]   = Y_2[k][i];
      Temp[k][0]   = T_2[k];
      H_tot[k][0] += Y_2[k][i]*h_i[i];
    }
  }
  
  double Y_Z1, Y_Z2, H_Z1, H_Z2, T_Z1, T_Z2, Z0;
  double *x1 = phi[0], *x2 = phi[1];
  /* reconstruct a 2D field */
  for (k=1; k<n2; k++) {
    for (j=1; j<n1; j++) {
      getZfromGrid(driver->TRANS_Z2, phi[0][j], phi[1][k], &Z1, &Z2);
      Z0 = Z1+Z2;
      jL = findAdjacentGridPoint(n1,phi[0],Z0);
      kB = findAdjacentGridPoint(n2,phi[1],Z0);
      
      // linearly interpolate species and enthalpy
      for (i=0; i<nSpecies; i++) {
        Y_Z1 = linInterp(Z0,x1[jL],Y_1[jL][i],x1[jL+1],Y_1[jL+1][i]);
        Y_Z2 = linInterp(Z0,x2[kB],Y_2[kB][i],x2[kB+1],Y_2[kB+1][i]);
        Y[k][j][i] = Y_Z1 + Z2/Z0*(Y_Z2 - Y_Z1);
      }
      
      H_Z1 = linInterp(Z0,x1[jL],H_tot[0][jL],x1[jL+1],H_tot[0][jL+1]);
      H_Z2 = linInterp(Z0,x2[kB],H_tot[kB][0],x2[kB+1],H_tot[kB+1][0]);
      H_tot[k][j] = H_Z1 + Z2/Z0*(H_Z2 - H_Z1);
      
      
      /*  T_Z1 = linInterp(Z0,x1[jL],T_1[jL],x1[jL+1],T_1[jL+1]); */
      /*           T_Z2 = linInterp(Z0,x2[kB],T_2[kB],x2[kB+1],T_2[kB+1]); */
      /*           Temp[k][j] = T_Z1 + Z2/Z0*(T_Z2 - T_Z1); */
      // compute temperature
      Temp[k][j] = MDUCcalcTfromH(nSolve,H_tot[k][j],Y[k][j],1000,mem);
    }
  }
  MDUCsetSolution(timeStart, P, Y, Temp, mem);

  for (k=0; k<n2; k++) {
    for (j=0; j<n1; j++) {
      for (i=0; i<nSpecies; i++)
        Y_pre[k][j][i] = Y[k][j][i];
      T_pre[k][j] = Temp[k][j];
    }
  }

  // Find out if H direction. If yes, find grid points in Z
  nZ = 0;
  if (driver->H_DIM == 1) {
    nZ = nPhi[1];
    dirZ = 1;
  }
  if (driver->H_DIM == 2) {
    nZ = nPhi[0];
    dirZ = 0;
  }

  // allocate Sdot
  if (driver->H_DIM != 0)
    Sdot = (double*) calloc(nZ, sizeof(double));
  
  /* ------------------------------------------ */
  /* ----- initialise the flow parameters ----- */
  /* ------------------------------------------ */
  MDUCinitFlow(timeStart,P,xmin[0][1],xmax[0][1],0.0,0.0,mem1);
  MDUCinitFlow(timeStart,P,xmin[0][2],xmax[0][2],0.0,0.0,mem2);
  
  // read in initial data from the flow file
  if (driver->CHI_MODE == -1 || driver->PDF_MODE == -1) {
    MDUCreadFlowData(chi1,chi12,chi2,pdf,Sdot,ntIn[0],mem);
    
    // force chi1 to be zero at boundaries Z1=0 and Z1=1
    for (i=0; i<nPhi[1]; i++)
      chi1[i][0] = chi1[i][nPhi[0]-1] = 0;
    
    // force chi2 to be zero at boundaries Z2=0 and Z2=1
    for (i=0; i<nPhi[0]; i++)
      chi2[0][i] = chi2[nPhi[1]-1][i] = 0;
    
    for (k=0; k<n2; k++) { pdfM_2[k] = 0.0; chiM_2[k] = 0.0; }
    for (j=0; j<n1; j++) { pdfM_1[j] = 0.0; chiM_1[j] = 0.0; }
    for (k=0; k<n2; k++) {
      for (j=0; j<n1; j++) {
        pdfM_1[j] += pdf[k][j];
        pdfM_2[k] += pdf[k][j];
      } /* end loop over columns */
    } /* end loop over rows */ 

    for (k=0; k<n2; k++) {
      for (j=0; j<n1; j++) {
        if (pdfM_1[j] > 0.0)
          chiM_1[j] += chi1[k][j]*pdf[k][j]/pdfM_1[j];
        if (pdfM_2[k] > 0.0)
          chiM_2[k] += chi2[k][j]*pdf[k][j]/pdfM_2[k];
      } /* end loop over columns */
    } /* end loop over rows */ 
    
      // transform variables if 2D mixture fraction
      /* if (driver->H_DIM == 0 && driver->TRANS_Z2 > 0) */
/*         { */
/*           MDUCtransformVar(1, driver->TRANS_Z2, nPhi[0], nPhi[1], */
/*                            phi[0], phi[1], chi1); */
/*           MDUCtransformVar(1, driver->TRANS_Z2, nPhi[0], nPhi[1], */
/*                            phi[0], phi[1], chi2); */
/*           MDUCtransformVar(1, driver->TRANS_Z2, nPhi[0], nPhi[1], */
/*                            phi[0], phi[1], chi12); */
/*           MDUCtransformVar(1, driver->TRANS_Z2, nPhi[0], nPhi[1], */
/*                            phi[0], phi[1], pdf); */
/*         } */

      // interpolate scalar dissipation rate      
      if (driver->CHI_INTERP == 1) {
        MDUCinterpChiDNS(&chiM_1, NULL, NULL, mem1);
        MDUCinterpChiDNS(&chiM_2, NULL, NULL, mem2);
      }
      else if (driver->CHI_INTERP == 2) {
        MDUCfitChi1D(nPhi[0], xmin[0][1], xmax[0][1], phi[0], chiM_1);
        MDUCfitChi1D(nPhi[1], xmin[0][2], xmax[0][2], phi[1], chiM_2);
      }
      else if (driver->CHI_INTERP == 3) {
        MDUClaplaceChi1D(nPhi[0],xmin[0][1],xmax[0][1],phi[0],chiM_1);
        MDUClaplaceChi1D(nPhi[1],xmin[0][2],xmax[0][2],phi[1],chiM_2);
      }
      
      /* set pdf bounds to whole grid to be safe (only changed if pdf
	 computed) */
      for (j = 0; j<nPhi[1]; j++) {
        pdfBound[j][0] = 0;
        pdfBound[j][1] = nPhi[0]-1;
      }
      pdfB_1[0] = pdfB_2[1] = 0;
      pdfB_1[1] = n1-1;
      pdfB_2[1] = n2-1;
  }
  
  /* ----- scalar dissipation rate ----- */
  if (driver->CHI_MODE != -1) // compute from reference value
    calcChi2D(io->chiMode, nPhi, phi, ZRef[0], chiRef[0], chi1, chi12, chi2,
              driver->TRANS_Z2, driver->CHI_SCALE);

  // set scalar dissipation rate into pre timestep
  MDUCsetScalarDiss(nPhi[0], 1, &chiM_1, NULL, NULL, mem1);
  MDUCsetScalarDiss(nPhi[1], 1, &chiM_2, NULL, NULL, mem2);
  
  // compute mean density if necessary
  if (driver->PRES_MODE == 1) {
    if (driver->PDF_MODE != -1 && Zmean) // compute PDF from moments
      pdfFlag = MDUCcalcPDF_2D(io->pdfMode,nPhi,phi,Zmean[0][1],Zvar[0][1],
                               Zmean[0][2],Zvar[0][2],Zvar[0][0],
                               pdf,pdfBound,driver->TRANS_Z2);
    
    // get mean density to conserve for constant volume cases
    meanRho_0 = MDUCcalcMeanDensity(P,NULL,pdf,pdfBound,pdfFlag,mem);
  }

  // write error file header
  if (driver->PRES_MODE == 1 && myid == 0) {
    writeErrorFileHeader(io->errorFile, 2);
    writeErrorFile(io->errorFile,timeStart,P,PIn[0],xmin[0][1],xmax[0][1],
                   xmin[0][2],xmax[0][2],mem1);
  }
  
  FILE *hrFile = fopen("HR_compare.dat", "w");
  fprintf(hrFile, "# time[s] HR_cvode[J/s] HR_dt[J/s]\n");
  fflush(hrFile);

#ifdef ACCESS_MPI
  startRun = MPI_Wtime();
#else
  time(&startRun); // start recording run time
#endif /* ACCESS_MPI */
  /* -------------------------------------- */
  /* ----- begin time integration CFD ----- */
  /* -------------------------------------- */
  fprintf(stderr, "%6s %15s %15s %15s %15s\n", "iter.", "time[s]", "dt[s]",
          "P[Pa]", "Tmax[K]");
  for (k=0, Tmax=-1.0e20; k<n2; k++)
    for (j=0; j<n1; j++) if (Temp[k][j] > Tmax) Tmax = Temp[k][j];
  fprintf(stderr, "%6i %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %6i\n", 0, timeStart, dt,
          P, Tmax, xmin[0][1], xmax[0][1], xmin[0][2], xmax[0][2],driver->TRANS_Z2);
  for (nt=1; nt<ntSteps; nt++) {
    // dereference data from input file
    curTime  = timeIn[nt];
    dt       = timeIn[nt] - timeIn[nt-1];
    
    // read data from flow file
    if (driver->CHI_MODE == -1 || driver->PDF_MODE == -1) {
      MDUCreadFlowData(chi1, chi12, chi2, pdf, Sdot, ntIn[nt], mem);
      
      // force chi1 to be zero at boundaries Z1=0 and Z1=1
      for (i=0; i<nPhi[1]; i++)
        chi1[i][0] = chi1[i][nPhi[0]-1] = 0;
      
      // force chi2 to be zero at boundaries Z2=0 and Z2=1
      for (i=0; i<nPhi[0]; i++)
        chi2[0][i] = chi2[nPhi[1]-1][i] = 0;
      
      for (k=0; k<n2; k++) { pdfM_2[k] = 0.0; chiM_2[k] = 0.0; }
      for (j=0; j<n1; j++) { pdfM_1[j] = 0.0; chiM_1[j] = 0.0; }
      for (k=0; k<n2; k++) {
        for (j=0; j<n1; j++) {
          pdfM_1[j] += pdf[k][j];
          pdfM_2[k] += pdf[k][j];
        } /* end loop over columns */
      } /* end loop over rows */ 
      
      for (k=0; k<n2; k++) {
        for (j=0; j<n1; j++) {
          if (pdfM_1[j] > 0.0)
            chiM_1[j] += chi1[k][j]*pdf[k][j]/pdfM_1[j];
          if (pdfM_2[k] > 0.0)
            chiM_2[k] += chi2[k][j]*pdf[k][j]/pdfM_2[k];
        } /* end loop over columns */
      } /* end loop over rows */ 
      
      // interpolate scalar dissipation rate      
      if (driver->CHI_INTERP == 1) {
        MDUCinterpChiDNS(&chiM_1, NULL, NULL, mem1);
        MDUCinterpChiDNS(&chiM_2, NULL, NULL, mem2);
      }
      else if (driver->CHI_INTERP == 2) {
        MDUCfitChi1D(nPhi[0], xmin[nt][1], xmax[nt][1], phi[0], chiM_1);
        MDUCfitChi1D(nPhi[1], xmin[nt][2], xmax[nt][2], phi[1], chiM_2);
      }
      else if (driver->CHI_INTERP == 3) {
        MDUClaplaceChi1D(nPhi[0],xmin[nt][1],xmax[nt][1],phi[0],chiM_1);
        MDUClaplaceChi1D(nPhi[1],xmin[nt][2],xmax[nt][2],phi[1],chiM_2);
      }
    }
    
    /* ----- calculate scalar dissipation rate ----- */
    if (driver->CHI_MODE != -1) // compute from reference value
      calcChi2D(io->chiMode,nPhi,phi,ZRef[nt],chiRef[nt],chi1,chi12,chi2,
                driver->TRANS_Z2,driver->CHI_SCALE);
    
    /* ----- main combustion solve ----- */
    MDUCsolve1D(dt,P,xmin[nt][1],xmax[nt][1],chiM_1,Y_1,T_1,mem1);
    MDUCsolve1D(dt,P,xmin[nt][2],xmax[nt][2],chiM_2,Y_2,T_2,mem2);
    
    for (j=0; j<n1; j++) {
      MDUCcomputeThermoData(h_i, cp_i, T_1[j], mem1);
      Temp[0][j] = T_1[j];
      for (i=0, H_tot[0][j]=0.0; i<nSpecies; i++) {
        Y[0][j][i]   = Y_1[j][i];
        H_tot[0][j] += Y_1[j][i]*h_i[i];
      }
    }
    for (k=0; k<n2; k++) {
      MDUCcomputeThermoData(h_i, cp_i, T_2[k], mem2);
      Temp[k][0] = T_2[k];
      for (i=0, H_tot[k][0]=0.0; i<nSpecies; i++) {
        Y[k][0][i]   = Y_2[k][i];
        H_tot[k][0] += Y_2[k][i]*h_i[i];
      }
    }
    
      /* reconstruct a 2D field */
/*       for (k=1; k<n2; k++) */
/*         { */
/*           for (j=1; j<n1; j++) */
/*             { */
/*               getZfromGrid(driver->TRANS_Z2, phi[0][j], phi[1][k], &Z1, &Z2); */

/*               if (Z1 + Z2 > 1.0) continue; */

/*               // linearly interpolate species and enthalpy */
/*               for (i=0; i<nSpecies; i++) */
/*                 Y[k][j][i] = Y_1[j][i] + Z2/(Z1+Z2)*(Y_2[k][i] - Y_1[j][i]); */
/*               H_tot[k][j] = H_tot[0][j] + Z2/(Z1+Z2)*(H_tot[k][0]-H_tot[0][j]); */

/*               // compute temperature */
/*               Temp[k][j] = calcTfromH(nSpecies,nSolve,H_tot[k][j],Y[k][j],1000); */
/*             } */
/*         } */
    for (k=1; k<n2; k++) {
      for (j=1; j<n1; j++) {
        getZfromGrid(driver->TRANS_Z2, phi[0][j], phi[1][k], &Z1, &Z2);
        Z0 = Z1+Z2;
        jL = findAdjacentGridPoint(n1,phi[0],Z0);
        kB = findAdjacentGridPoint(n2,phi[1],Z0);
        
        // linearly interpolate species and enthalpy
        for (i=0; i<nSpecies; i++) {
          Y_Z1 = linInterp(Z0,x1[jL],Y_1[jL][i],x1[jL+1],Y_1[jL+1][i]);
          Y_Z2 = linInterp(Z0,x2[kB],Y_2[kB][i],x2[kB+1],Y_2[kB+1][i]);
          Y[k][j][i] = Y_Z1 + Z2/Z0*(Y_Z2 - Y_Z1);
        }
        
        H_Z1 = linInterp(Z0,x1[jL],H_tot[0][jL],x1[jL+1],H_tot[0][jL+1]);
        H_Z2 = linInterp(Z0,x2[kB],H_tot[kB][0],x2[kB+1],H_tot[kB+1][0]);
        H_tot[k][j] = H_Z1 + Z2/Z0*(H_Z2 - H_Z1);
        
        /* T_Z1 = linInterp(Z0,x1[jL],T_1[jL],x1[jL+1],T_1[jL+1]); */
          /*               T_Z2 = linInterp(Z0,x2[kB],T_2[kB],x2[kB+1],T_2[kB+1]); */
          /*               Temp[k][j] = T_Z1 + Z2/Z0*(T_Z2 - T_Z1); */
          // compute temperature
        Temp[k][j] = MDUCcalcTfromH(nSolve,H_tot[k][j],Y[k][j],1000,mem);
      }
    }
    MDUCsetSolution(curTime, P, Y, Temp, mem);
      
    /* update pressure if reading from file */
    if (myid == 0 && driver->PRES_MODE == -1) {
      // update pressure and rate
      P    = PIn[nt];
      dPdt = (P - P_pre)/dt;
      
      /* update the temperature base on the current dP/dt and
         temperature from last iteration */
      MDUCupdateTempDPDT(P, dPdt, Temp, mem);
      
      /* copy the final temperature back into solution array and
         reinitialise SUNDIALS */
      MDUCsetSolution(curTime, P, Y, Temp, mem);
    }
    
    /* ----- calculate Favre average quantities ----- */
    if (myid == 0 && (Zmean || driver->PDF_MODE == -1)) { 
      if (driver->PDF_MODE != -1)
        pdfFlag = MDUCcalcPDF_2D(io->pdfMode,nPhi,phi,Zmean[nt][1],
                                 Zvar[nt][1],Zmean[nt][2],Zvar[nt][2],
                                 Zvar[nt][0],pdf,pdfBound,driver->TRANS_Z2);
      
      // recompute pressure if constant volume
      if (driver->PRES_MODE == 1) {
	      /* iterate to converge the pressure/temperature that gives the
           correct density */
	      count   = 0;
	      meanRho = MDUCcalcMeanDensity(P, Temp, pdf, pdfBound, pdfFlag,
                                      mem);
	      while (fabs(1.0 - meanRho/meanRho_0) > 1.0e-12) {
          // update pressure and rate
          P    -= P*(meanRho/meanRho_0 - 1.0);
          dPdt  = (P - P_pre)/dt;
          
          /* update the temperature base on the current dP/dt and
             temperature from last iteration */
          MDUCupdateTempDPDT(P, dPdt, Temp, mem);
          
          /* compute the new mean density from the flamelet at current
             pressure and temperature */
          meanRho = MDUCcalcMeanDensity(P, Temp, pdf, pdfBound, pdfFlag,
                                        mem);
          
          if (++count > 100)
            mducErr("constant volume pressure not converged",
                    "postProcessMulti1D");
        }
        
	      /* copy the final temperature back into solution array and
           reinitialise SUNDIALS */
	      MDUCsetSolution(curTime, P, Y, Temp, mem);
	    }

      for (j=0; j<n1; j++) T_1[j] = Temp[0][j];
      MDUCsetSolution(curTime, P, &Y_1, &T_1, mem1);
      for (k=0; k<n2; k++) T_2[k] = Temp[k][0];
      MDUCsetSolution(curTime, P, &Y_2, &T_2, mem2);
      
      // Output scalars
      MDUCupdateTime(curTime, dt, nt, mem);
      MDUCwriteFlamelet(NULL, mem);
      MDUCwriteScalars(curTime,Zmean[nt],Zvar[nt],pdf,pdfBound,pdfFlag,mem);
      MDUCwriteScalars(curTime,Zmean[nt],Zvar[nt],&pdfM_1,&pdfB_1,1,mem1);
      MDUCwriteScalars(curTime,Zmean[nt],Zvar[nt],&pdfM_2,&pdfB_2,1,mem2);
    }
    
    /* compute the heat release */
    for (k=0; k<n2; k++) {
      for (j=0; j<n1; j++) {
        MDUCcomputeThermoData(h_i,   cp_i,   Temp[k][j], mem1);
        MDUCcomputeThermoData(h_pre, cp_pre, T_pre[k][j], mem1);
        wmix     = calcMixtureW(nSpecies,Y[k][j],Mw);
        wmix_pre = calcMixtureW(nSpecies,Y_pre[k][j],Mw);
        rho_n    = wmix*P/(RGAS*Temp[k][j]);
        rho_nm1  = wmix_pre*P_pre/(RGAS*T_pre[k][j]);
        
        HR[k][j] = 0.0;
        for (i=0, h_n=0.0, h_nm1=0.0; i<nSpecies; i++) {
          h_n   += Y[k][j][i]*h_i[i];
          h_nm1 += Y_pre[k][j][i]*h_pre[i];
          HR[k][j] += (Y[k][j][i] - Y_pre[k][j][i])/dt*h_pre[i];
        }
        
        /* HR[k][j] = ( 0.5*(rho_n+rho_nm1)*(h_n-h_nm1) + */
        /*                            0.5*(h_n+h_nm1)*(rho_n-rho_nm1) )/dt; */
        
        /*               HR[k][j] = hr_local/dt; */
        HR[k][j] *= -rho_nm1;
      }
    }
    MDUCupdateHRrate(HR, H_tot, mem);
    /* compute the heat release */
    for (j=0; j<n1; j++) {
      MDUCcomputeThermoData(h_i,   cp_i,   T_1[j], mem);
      MDUCcomputeThermoData(h_pre, cp_pre, T_1_m1[j], mem);
      wmix     = calcMixtureW(nSpecies,Y_1[j],Mw);
      wmix_pre = calcMixtureW(nSpecies,Y_1_m1[j],Mw);
      rho_n    = wmix*P/(RGAS*T_1[j]);
      rho_nm1  = wmix_pre*P_pre/(RGAS*T_1_m1[j]);
      
      HR_1_est[j] = 0.0;
      for (i=0, h_n=0.0, h_nm1=0.0; i<nSpecies; i++) {
        h_n   += Y_1[j][i]*h_i[i];
        h_nm1 += Y_1_m1[j][i]*h_pre[i];
        HR_1_est[j] += (Y_1[j][i] - Y_1_m1[j][i])/dt*h_i[i];
      }
      
/*           HR_1_est[j] = ( 0.5*(rho_n+rho_nm1)*(h_n-h_nm1) + */
/*                       0.5*(h_n+h_nm1)*(rho_n-rho_nm1) )/dt; */
/*           HR_1_est[j] = rho_nm1*h_nm1*(h_n-h_nm1)/dt; */
      HR_1_est[j] *= -rho_nm1;
    }
    
#ifdef ACCESS_MPI
      // send correct pressure to all processes
    MPI_Bcast(&P, 1, MPI_DOUBLE_PRECISION, 0, mpi->world);
    MDUCsetPressure(-1, P, mem);
#endif /* ACCESS_MPI */

//    // dump the flow file
//    if (io->writeFlow && myid == 0)
//      MDUCwriteFlowData(chi1, chi2, pdf, Sdot, mem);
    
    // write error file
    if (driver->PRES_MODE == 1 && myid == 0)
      writeErrorFile(io->errorFile,timeIn[nt],P,PIn[nt],xmin[nt][1],
                     xmax[nt][1],xmin[nt][2],xmax[nt][2],mem);
    
    for (k=0, Tmax=-1.0e20; k<n2; k++)
      for (j=0; j<n1; j++) if (Temp[k][j] > Tmax) Tmax = Temp[k][j];
    fprintf(stderr, "%6i %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n", nt, curTime, dt,
            P, Tmax, xmin[nt][1], xmax[nt][1], xmin[nt][2], xmax[nt][2]);

    MDUCgetVar("HRrate", &HR_1, mem1);
    // compute the average heat release
    for (j=0, HRavg1=0.0, HRavg1_est=0.0; j<n1; j++) {
      HRavg1     += HR_1[j]*pdfM_1[j];
      HRavg1_est += HR_1_est[j]*pdfM_1[j];
    }
    fprintf(hrFile, "%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n", curTime, HRavg1, HRavg1_est, rho_n, rho_nm1, h_n, h_nm1);
    fflush(hrFile);

    // store pressure
    P_pre = P;
    for (k=0; k<n2; k++) {
      for (j=0; j<n1; j++) {
        for (i=0; i<nSpecies; i++)
          Y_pre[k][j][i] = Y[k][j][i];
        T_pre[k][j] = Temp[k][j];
        
        for (i=0; i<nSpecies; i++)
          Y_1_m1[j][i] = Y_1[j][i];
        T_1_m1[j] = T_1[j];
      }
    }
    
  } // end time integration
#ifdef ACCESS_MPI
  endRun = MPI_Wtime();
  elapsedTime = endRun - startRun;
#else
  time(&endRun);
  elapsedTime = difftime(endRun, startRun);
#endif /* ACCESS_MPI */
  if (myid == 0)
    fprintf(stderr, "Total time elapsed: %5.3f\n", elapsedTime);

  fclose(hrFile);
  
  // free memory allocated by reading file
  if (ntIn)   free(ntIn);
  if (timeIn) free(timeIn);
  if (PIn)    free(PIn);
  for (nt=0; nt<ntSteps; nt++) {
    if (chiRef[nt]) free(chiRef[nt]);
    if (ZRef[nt])   free(ZRef[nt]);
    if (Zmean[nt])  free(Zmean[nt]);
    if (Zvar[nt])   free(Zvar[nt]);
    if (xmin[nt])   free(xmin[nt]);
    if (xmax[nt])   free(xmax[nt]);
  }
  if (chiRef) free(chiRef);
  if (ZRef)   free(ZRef);
  if (Zmean)  free(Zmean);
  if (Zvar)   free(Zvar); 
  if (xmin)   free(xmin);
  if (xmax)   free(xmax);
  
  // free locally allocated memory
  for (k=0; k<nPhi[1]; k++) {
    for (j=0; j<nPhi[0]; j++) { free(Y[k][j]); free(Y_pre[k][j]); }
    free(Y[k]); free(Temp[k]);
    free(Y_pre[k]); free(T_pre[k]);
    free(H_tot[k]); free(HR[k]);
    free(pdf[k]); free(pdfBound[k]);
    free(chi0[k]); free(chi1[k]); free(chi2[k]); free(chi12[k]);
  }
  free(Y); free(Temp);
  free(Y_pre); free(T_pre); free(H_tot); free(HR);
  free(pdf); free(pdfBound);
  free(chi0); free(chi1); free(chi2); free(chi12);
  free(nPhi); free(phi[0]); free(phi[1]); free(phi);
  free(h_i); free(cp_i); free(h_pre); free(cp_pre); free(Mw);
  free(pdfB_1); free(pdfB_2);
  free(pdfM_1); free(pdfM_2); free(chiM_1); free(chiM_2);
  free(T_1); free(T_2); free(T_1_m1); free(HR_1); free(HR_1_est);
  for (j=0; j<n1; j++) { free(Y_1[j]); free(Y_1_m1[j]); }
  free(Y_1);
  for (k=0; k<n2; k++) free(Y_2[k]); free(Y_2); 
  
  // free all flamelet memory
  if (myid == 0) fprintf(stderr, "freeing flamelet memory...\n");
  MDUCdestroy(mem);
  MDUCdestroy(mem1);
  MDUCdestroy(mem2);
  
#ifdef ACCESS_MPI
  // call MPI_Init before handling any other options
  MPI_Finalize();
#endif // ACCESS_MPI
}
/* -------------------------------------------------------------------------- */

int inputPostProcess(int **nt, double **tIn, double **P, double ***chiR,
                     double ***ZR, double ***Zmean, double ***Zvar,
                     double ***xmin, double ***xmax, void *mem)
{
  int i, n;
  timeStruct *time = MDUCgetTimeStruct(mem);

  if (time->nIter > 0) {
    if (time->t_end > 0.0)
      time->dt = (time->t_end - time->t0)/time->nIter;
    else
      time->t_end = time->t0 + time->nIter*time->dt;
  }
  else
    time->nIter = (time->t_end - time->t0)/time->dt;

  time->nIter++;

  /* allocate memory */
  nt[0]    = (int*)     calloc(time->nIter, sizeof(int));
  tIn[0]   = (double*)  calloc(time->nIter, sizeof(double));
  P[0]     = (double*)  calloc(time->nIter, sizeof(double));
  chiR[0]  = (double**) calloc(time->nIter, sizeof(double*));
  ZR[0]    = (double**) calloc(time->nIter, sizeof(double*));
  Zmean[0] = (double**) calloc(time->nIter, sizeof(double*));
  Zvar[0]  = (double**) calloc(time->nIter, sizeof(double*));
  xmin[0]  = (double**) calloc(time->nIter, sizeof(double*));
  xmax[0]  = (double**) calloc(time->nIter, sizeof(double*));

  for (i=0; i<time->nIter; i++) {
    chiR[0][i]  = (double*) calloc(3, sizeof(double));
    ZR[0][i]    = (double*) calloc(3, sizeof(double));
    Zmean[0][i] = (double*) calloc(3, sizeof(double));
    Zvar[0][i]  = (double*) calloc(3, sizeof(double));
    xmin[0][i]  = (double*) calloc(3, sizeof(double));
    xmax[0][i]  = (double*) calloc(3, sizeof(double));
  }

  /* set values into first array value */
  nt[0][0]  = 0;
  tIn[0][0] = time->t0;
  P[0][0]   = MDUCgetPressure(0, mem);
  MDUCgetGridBound(&xmin[0][0][1],&xmax[0][0][1],&xmin[0][0][2],&xmax[0][0][2],
                   mem);
  for (i=0; i<3; i++) {
    chiR[0][0][i] = MDUCgetChiRef(i, mem);
    ZR[0][0][i]   = MDUCgetZref(i, mem);
    MDUCgetZmoments(i, &Zmean[0][0][i], &Zvar[0][0][i], mem);
  }

  /* fill the arrays */
  for (n=1; n<time->nIter; n++) {
    nt[0][n]  = n;
    tIn[0][n] = tIn[0][n-1] + time->dt;
    P[0][n]   = P[0][0];
    for (i=0; i<3; i++) {
      chiR[0][n][i] = chiR[0][0][i];
      ZR[0][n][i]   = ZR[0][0][i];
      xmin[0][n][i] = xmin[0][0][i];
      xmax[0][n][i] = xmax[0][0][i];
      MDUCgetZmoments(i, &Zmean[0][n][i], &Zvar[0][n][i], mem);
    }
  }
  return time->nIter;
}
/* -------------------------------------------------------------------------- */

int readFluentData(char *dataFile, int nYFav, char **species, double **P,
                   double **RHO, double ***YFav, double **HFav, double **K,
                   double **D, double **Z1Fav, double **Z1Var, double **Z2Fav,
                   double **Z2Var, double **Z12Var, double **V)
{
  printf("Reading file: %s\n", dataFile);
  int i, j, k, m, nCells, nCols;
  float *f_tmp;
  char cTest;

  // integers for storing the column index of CFD variables
  int iP, iRho, *iY, iH, iK, iD, iZ1, iZ1Var, iZ2, iZ2Var, iZ12Var, iV;

  // strings to store the tag for each line and comments (if any)
  char  tag[50], rawLine[1000];
  char **lowSpecies;
  
  // open the file for reading
  FILE *pFile = fopen(dataFile, "r");
  if (!pFile)
    {
      fprintf(stderr, "Error: %s file not present!\n", dataFile);
      exit(2);
    }

  // read the last line to determine how many cells there are to read
  while (!feof(pFile))
    fgets(rawLine, 1000, pFile);
  sscanf(rawLine, "%d", &nCells); // read in number of cells from last line

  // go back to the beginning of the file
  rewind(pFile);

  // allocate all memory
  iY = (int*) calloc(nYFav, sizeof(int));
  lowSpecies = (char**) calloc(nYFav, sizeof(char*));

  if (YFav[0] == NULL)
    YFav[0]  = (double**) calloc(nYFav, sizeof(double*));

  P[0]      = (double*) realloc(P[0], nCells*sizeof(double));
  RHO[0]    = (double*) realloc(RHO[0], nCells*sizeof(double));
  HFav[0]   = (double*) realloc(HFav[0], nCells*sizeof(double));
  K[0]      = (double*) realloc(K[0], nCells*sizeof(double));
  D[0]      = (double*) realloc(D[0], nCells*sizeof(double));
  Z1Fav[0]  = (double*) realloc(Z1Fav[0], nCells*sizeof(double));
  Z1Var[0]  = (double*) realloc(Z1Var[0], nCells*sizeof(double));
  Z2Fav[0]  = (double*) realloc(Z2Fav[0], nCells*sizeof(double));
  Z2Var[0]  = (double*) realloc(Z2Var[0], nCells*sizeof(double));
  Z12Var[0] = (double*) realloc(Z2Var[0], nCells*sizeof(double));
  V[0]      = (double*) realloc(V[0], nCells*sizeof(double));

  for (i=0; i<nYFav; i++)
    {
      YFav[0][i] = (double*) realloc(YFav[0][i], nCells*sizeof(double));
      lowSpecies[i] = (char*) calloc(25, sizeof(char));
    }

  // copy all the species into lower case for Fluent format
  for (i=0; i<nYFav; i++)
    {
      strcpy(lowSpecies[i], species[i]);
      for (j=0; j<strlen(species[i]); j++)
        if (isalpha(species[i][j]))
          lowSpecies[i][j] = tolower(species[i][j]);
    }

  // read in the header and get the number of columns
  // also determine the index of each parameter needed
  nCols = 0;
  do {
    fscanf(pFile, "%s", tag);
    if (!strcmp(tag, "absolute-pressure"))
      iP = nCols;
    else if (!strcmp(tag, "density"))
      iRho = nCols;
    else if (!strcmp(tag, "turb-kinetic-energy"))
      iK = nCols;
    else if (!strcmp(tag, "turb-diss-rate"))
      iD = nCols;
    else if (!strcmp(tag, "uds-0-scalar"))
      iZ1 = nCols;
    else if (!strcmp(tag, "uds-1-scalar"))
      iZ1Var = nCols;
    else if (!strcmp(tag, "uds-2-scalar" ))
      iZ2 = nCols;
    else if (!strcmp(tag, "uds-3-scalar"))
      iZ2Var = nCols;
    else if (!strcmp(tag, "uds-4-scalar"))
      iZ12Var = nCols;
    else if (!strcmp(tag, "uds-5-scalar"))
      iH = nCols;
    else if (!strcmp(tag, "cell-volume"))
      iV = nCols;
    else
      {
        for (i=0; i<nYFav; i++)
          if (!strcmp(tag, lowSpecies[i]))
            iY[i] = nCols;
      }

    nCols++;

    // read in all whitespace until next header
    do {
      cTest = fgetc(pFile);
    } while (cTest == ' ');
    ungetc(cTest, pFile);
  }
  while (cTest != '\n');

  // create an array to store the contents of a line
  f_tmp = (float*) calloc(nCells, sizeof(float));

  double Pavg = 0;
  double mass = 0;
  double mFuel = 0;
  double volume = 0;
  for (i=0; i<nCells; i++)
    {
      for (j=0; j<nCols; j++)
        fscanf(pFile, "%f", &f_tmp[j]);

      P[0][i]      = (double)(f_tmp[iP]);
      RHO[0][i]    = (double)(f_tmp[iRho]);
      HFav[0][i]   = (double)(f_tmp[iH]);
      K[0][i]      = (double)(f_tmp[iK]);
      D[0][i]      = (double)(f_tmp[iD]);
      Z1Fav[0][i]  = (double)(f_tmp[iZ1]);
      Z1Var[0][i]  = (double)(f_tmp[iZ1Var]);
      Z2Fav[0][i]  = (double)(f_tmp[iZ2]);
      Z2Var[0][i]  = (double)(f_tmp[iZ2Var]);
      Z12Var[0][i] = (double)(f_tmp[iZ12Var]);
      V[0][i]      = (double)(f_tmp[iV]);

      for (k=0; k<nYFav; k++)
        YFav[0][k][i] = (double)(f_tmp[iY[k]]);
     
      mass += RHO[0][i]*V[0][i];
      Pavg += P[0][i]*V[0][i];
      mFuel += YFav[0][0][i]*RHO[0][i]*V[0][i];
      volume += V[0][i];
    } // end loop over file

  printf("Total mass: %15.6e\n", mass);
  printf("Fuel mass: %15.6e\n", mFuel);
  printf("The average pressure is: %15.6e\n", Pavg/volume);

  fclose(pFile);

  free(iY);
  for (i=0; i<nYFav; i++) free(lowSpecies[i]);
  free(lowSpecies);

  printf("%i cells successfully read.\n", nCells);
  return nCells;
}
/* -------------------------------------------------------------------------- */

int readFileList(char *fileForm, int iStart, int iStop, int **nt, double **tIn,
                 double **P, double ***chiR, double ***ZR, double ***Zmean,
                 double ***Zvar, double ***xmin, double ***xmax, void *mem)
{
  int i, n, len, isochor = 0, rot;

  char  file[200], buffer[200], prefix[200], suffix[20], *head;
  FILE  *pFile;

  double tmp;
  mducStruct  *mduc   = (mducStruct*) mem;
  ioStruct    *io     = mduc->io;
  driveStruct *driver = mduc->driver;

  int myid = 0;
#ifdef ACCESS_MPI
  mpiStruct   *mpi;
  getMPIStruct(&mpi, mem);
  MPI_Comm_rank(mpi->world, &myid); // get processor ID
#endif // ACCESS_MPI

  /* ----- check to see what files are there ----- */
  // store prefix/suffix
  head = strtok(fileForm, ".");
  strcpy(prefix, head);
  head = strtok(NULL, "\n");
  if (head == NULL)
    strcpy(suffix, "");
  else
    strcpy(suffix, head);
  
  len = 0;
  if (myid == 0) fprintf(stderr, "checking for files: ");
  nt[0] = (int*) calloc(1, sizeof(int));
  // loop over files to see if they exist
  rot = 0;
  for (i=iStart; i<=iStop; i++)
    {
      sprintf(buffer, "%0*i", io->fileNumPrec, i);
      strcpy(file, prefix);
      strcat(file, buffer);
      if (head != NULL)
        {
          strcat(file, ".");
          strcat(file, suffix);
        }

      pFile = fopen(file, "r");
      if (pFile)
        {
          len++;
          if (myid == 0) fprintf(stderr, ".");
          nt[0]     = (int*) realloc(nt[0], len*sizeof(int));
          nt[0][len-1] = i;
          fclose(pFile);
        }
      else
        {
          if (myid == 0) fprintf(stderr, "-");
        }
    } // end file loop
  if (myid == 0) fprintf(stderr, "\n%i files found\n\n", len);
  
  // allocate memory (some if it is dummy)
  tIn[0] = (double*) calloc(len, sizeof(double));
  P[0]   = (double*) calloc(len, sizeof(double));

  // check the first file to see if necessary information is available
  if (driver->PRES_MODE == 1)
    {
      // set pressure mode temporarily so that P is set in solver
      isochor = 1;
      driver->PRES_MODE = -1;
    }
  driver->readBody = 0; // set so that only the header info will be read in
  tIn[0][0] = MDUCreadFlowData(NULL, NULL, NULL, NULL, NULL, nt[0][0], mem);

  if (tIn[0][0] < 0.0)
    mducErr("time information not available in file","readFileList");

  P[0][0] = MDUCgetPressure(0, mem);
  if (P[0][0] == 0.0)
    mducErr("pressure information not available in file","readFileList");
  
  chiR[0]  = (double**) calloc(len, sizeof(double*));
  ZR[0]    = (double**) calloc(len, sizeof(double*));
  Zmean[0] = (double**) calloc(len, sizeof(double*));
  Zvar[0]  = (double**) calloc(len, sizeof(double*));
  xmin[0]  = (double**) calloc(len, sizeof(double*));
  xmax[0]  = (double**) calloc(len, sizeof(double*));
  for (i=0; i<len; i++)
    {
      chiR[0][i]  = (double*) calloc(3, sizeof(double));
      ZR[0][i]    = (double*) calloc(3, sizeof(double));
      Zmean[0][i] = (double*) calloc(3, sizeof(double));
      Zvar[0][i]  = (double*) calloc(4, sizeof(double));
      xmin[0][i]  = (double*) calloc(3, sizeof(double));
      xmax[0][i]  = (double*) calloc(3, sizeof(double));
    }

  MDUCgetGridBound(&xmin[0][0][1],&xmax[0][0][1],&xmin[0][0][2],&xmax[0][0][2],
                   mem);

  // read in the rest of the data
  for (n=1; n<len; n++)
    {
      tIn[0][n] = MDUCreadFlowData(NULL,NULL,NULL,NULL,NULL,nt[0][n], mem);
      P[0][n]   = MDUCgetPressure(0, mem);
      MDUCgetGridBound(&xmin[0][n][1], &xmax[0][n][1],
                       &xmin[0][n][2], &xmax[0][n][2], mem);
    }

  /* if isochoric, reset the pressure mode appropriately so pressure is no
     longer stored when read in */
  if (isochor)
    driver->PRES_MODE = 0;

  // return to reading entire file
  driver->readBody = 1;

  return len;
}
/* -------------------------------------------------------------------------- */

int readScalarFile(char *file, int **nt, double **tIn, double **P, double ***T,
                   double ***chiR, double ***ZR, double ***Zmean,
                   double ***Zvar, double ***xmin, double ***xmax, void *mem)
{
  int     i, j, n, len, CA = 0, nHead, nVars, index, var_alloc;
  int     iTime,iP,iT[3],iZR[3],iChiR[3],iMean[3],iVar[4],iMin[3],iMax[3];
  double  tmp;
  
  /* strings to store the tag for each line and comments (if any) */
  char  tag[50], comment[1000], header[1000], *token;
  FILE  *pFile = fopen(file, "r");
  
  if (!pFile) {
    fprintf(stderr, "Error: post-processing file %s not found\n", file);
    exit(2);
  }
  
  /* check to see if this is a legacy CA file */
  fscanf(pFile, "%s", tag);
  if (!strcmp(tag, "RPM")) CA = 1;
  rewind(pFile);

  if (!CA) { /* read in number of header lines */
    nHead = -1;
    do {
      nHead++;
      strcpy(comment, "\n");
      fgets(comment, 1000, pFile);
    } while(comment[0] == '#');
    if (nHead == 0) {
      fprintf(stderr, "Error: post-process file needs header\n");
      exit(2);
    }
    rewind(pFile);
  }
  else
    nHead = 3;

  /* get the number of data lines */
  len = 0;
  while (!feof(pFile)) {
    strcpy(comment, "\n");
    fgets(comment, 1000, pFile);
    if (strcmp(comment, "\n"))
      len++;
  }
  len -= nHead;
  rewind(pFile);

  /* get the header line */
  for (i=0; i<nHead; i++) fgets(header,1000,pFile);

  /* parse the header to get the indices */
  if (CA) index = 0;
  else    index = -1;

  iTime = iP = iVar[3] = -1;
  for (i=0; i<3; i++)
    iT[i]=iChiR[i]=iZR[i]=iMean[i]=iVar[i]=iMin[i]=iMax[i] = -1;
  token = strtok(header, " \t\n");
  while (token != NULL) {
    /* convert string to lower case */
    strcpy(tag, token);
    for (i=0; i<strlen(tag); i++) tag[i] = tolower(tag[i]);
    
    if      (!strncmp(tag, "time", 4))      iTime    = index;
    else if (!strncmp(tag, "pressure", 5))  iP       = index;
    else if (!strcmp(tag, "t1"))            iT[0]    = index;
    else if (!strcmp(tag, "t2"))            iT[1]    = index;
    else if (!strcmp(tag, "t3"))            iT[2]    = index;
    else if (!strncmp(tag, "chi0ref", 4))   iChiR[0] = index;
    else if (!strncmp(tag, "chi1ref", 5) || !strncmp(tag,"chiref", 4) ||
             !strncmp(tag, "sci", 3) )      
      iChiR[1] = index;
    else if (!strncmp(tag, "chi2ref", 5))   iChiR[2] = index;
    else if (!strncmp(tag, "z0ref", 3))     iZR[0] = index;
    else if (!strncmp(tag, "z1ref", 3) || !strncmp(tag, "zref", 2))
      iZR[1] = index;
    else if (!strncmp(tag, "z2ref", 3) || !strncmp(tag, "href", 4) ) 
      iZR[2] = index;
    else if (!strncmp(tag, "z0mean", 4))
      iMean[0] = index;
    else if (!strncmp(tag, "z1mean", 4) || !strcmp(tag, "zmean"))
      iMean[1] = index;
    else if (!strncmp(tag, "z2mean", 4) || !strncmp(tag, "hmean", 5))
      iMean[2] = index;
    else if (!strncmp(tag, "z0var", 3))
      iVar[0] = index;
    else if (!strncmp(tag, "z1var", 3) || !strcmp(tag,"zvar"))
      iVar[1] = index;
    else if (!strncmp(tag, "z2var", 3) || !strncmp(tag, "hvar", 4))
      iVar[2] = index;
    else if (!strncmp(tag, "z12var", 4))
      iVar[3] = index;
    else if (!strncmp(tag, "z0min", 4))
      iMin[0] = index;
    else if (!strncmp(tag, "z1min", 4) || !strcmp(tag, "zmin"))
      iMin[1] = index;
    else if (!strncmp(tag, "z2min", 4) || !strncmp(tag, "hmin", 4))
      iMin[2] = index;
    else if (!strncmp(tag, "z0max", 4))
      iMax[0] = index;
    else if (!strncmp(tag, "z1max", 4) || !strcmp(tag, "zmax"))
      iMax[1] = index;
    else if (!strncmp(tag, "z2max", 4) || !strncmp(tag, "hmax", 4))
      iMax[2] = index;
    
    token = strtok(NULL, " \t\n");
    index++;
  }
  nVars = index;

  if (iTime < 0) {
    fprintf(stderr, "Error: time not found in file\n");
    exit(2);
  }
  else {
    fprintf(stderr, "=============== Scalar File Summary ===============\n");
    fprintf(stderr, "  Header lines:     %3i\n", nHead);
    fprintf(stderr, "  No. of columns:   %3i\n", nVars);
    fprintf(stderr, "  No. of timesteps: %3i\n", len);
    fprintf(stderr, "  Variables:\n");
    fprintf(stderr, "  %15s %5i\n","time",iTime);
  }
  tIn[0] = (double*) calloc(len, sizeof(double));
  nt[0]  = (int*)    calloc(len, sizeof(int));
  
  /* allocate additional memory as necessary */
  if (iP > 0 || MDUCgetPressure(0, mem) > 0.0) {
    P[0]     = (double*) calloc(len, sizeof(double));
    fprintf(stderr, "  %15s %5i\n", "pressure", iP);
  }
  if (iT[0] > 0 || iT[1] > 0 || iT[2] > 0) {
    T[0]  = (double**) calloc(len, sizeof(double*));
    if (iT[0] > 0) fprintf(stderr, "  %15s %5i\n", "T1", iT[0]);
    if (iT[1] > 0) fprintf(stderr, "  %15s %5i\n", "T2", iT[1]);
    if (iT[2] > 0) fprintf(stderr, "  %15s %5i\n", "T3", iT[2]);
  }
  if (iChiR[0] > 0 || iChiR[1] > 0 || iChiR[2] > 0) {
    chiR[0]  = (double**) calloc(len, sizeof(double*));
    if (iChiR[0] > 0) fprintf(stderr, "  %15s %5i\n", "chi0ref", iChiR[0]);
    if (iChiR[1] > 0) fprintf(stderr, "  %15s %5i\n", "chi1ref", iChiR[1]);
    if (iChiR[2] > 0) fprintf(stderr, "  %15s %5i\n", "chi2ref", iChiR[2]);
  }
  if (iZR[0] > 0 || iZR[1] > 0 || iZR[2] > 0) {
    ZR[0]    = (double**) calloc(len, sizeof(double*));
    if (iZR[0] > 0) fprintf(stderr, "  %15s %5i\n", "z0ref", iZR[0]);
    if (iZR[1] > 0) fprintf(stderr, "  %15s %5i\n", "z1ref", iZR[1]);
    if (iZR[2] > 0) fprintf(stderr, "  %15s %5i\n", "z2ref", iZR[2]);
  }
  if (iMean[0] > 0 || iMean[1] > 0 || iMean[2] > 0) {
    Zmean[0] = (double**) calloc(len, sizeof(double*));
    if (iMean[0] > 0) fprintf(stderr, "  %15s %5i\n", "z0mean", iMean[0]);
    if (iMean[1] > 0) fprintf(stderr, "  %15s %5i\n", "z1mean", iMean[1]);
    if (iMean[2] > 0) fprintf(stderr, "  %15s %5i\n", "z2mean", iMean[2]);
  }
  if (iVar[0] > 0 || iVar[1] > 0 || iVar[2] > 0 || iVar[3] > 0) {
    Zvar[0]  = (double**) calloc(len, sizeof(double*));
    if (iVar[0] > 0) fprintf(stderr, "  %15s %5i\n", "z0var", iVar[0]);
    if (iVar[1] > 0) fprintf(stderr, "  %15s %5i\n", "z1var", iVar[1]);
    if (iVar[2] > 0) fprintf(stderr, "  %15s %5i\n", "z2var", iVar[2]);
    if (iVar[3] > 0) fprintf(stderr, "  %15s %5i\n", "z2var", iVar[3]);
  }

  /* always need min/max so allocate and set here */
  xmin[0] = (double**) calloc(len, sizeof(double*));
  xmax[0] = (double**) calloc(len, sizeof(double*));

  for (i=0; i<len; i++) {
    if (T[0])     T[0][i]     = (double*) calloc(3, sizeof(double));
    if (chiR[0])  chiR[0][i]  = (double*) calloc(3, sizeof(double));
    if (ZR[0])    ZR[0][i]    = (double*) calloc(3, sizeof(double));
    if (Zmean[0]) Zmean[0][i] = (double*) calloc(3, sizeof(double));
    if (Zvar[0])  Zvar[0][i]  = (double*) calloc(4, sizeof(double));
    xmin[0][i] = (double*)  calloc(3, sizeof(double));
    xmax[0][i] = (double*)  calloc(3, sizeof(double));
  }

  /* set Zmean to be negative if not present */
  if (iMean[0] > 0 || iMean[1] > 0 || iMean[2] > 0) {
    for (n=0; n<len; n++)
      for (j=0; j<3; j++)
        Zmean[0][n][j] = -1.0;
  }

  /* loop over data lines in file */
  for (n=0; n<len; n++) {
    nt[0][n] = n;
    /* loop over variables in line */
    for (i=0; i<nVars; i++) {
      /* get current variable */
      tmp = readReal(pFile);
      
      /* assign appropriately */
      if      (i == iTime) tIn[0][n] = tmp;
      else if (i == iP)    P[0][n]   = tmp;
      else {
        for (j=0; j<3; j++) {
          if      (i == iT[j])    T[0][n][j]     = tmp;
          else if (i == iChiR[j]) chiR[0][n][j]  = tmp;
          else if (i == iZR[j])   ZR[0][n][j]    = tmp;
          else if (i == iMean[j]) Zmean[0][n][j] = tmp;
          else if (i == iVar[j])  Zvar[0][n][j]  = tmp;
          else if (i == iMin[j])  xmin[0][n][j]  = tmp;
          else if (i == iMax[j])  xmax[0][n][j]  = tmp;
        }
      }
      
      /* if extent data not present, set to default (only possible for Z) */
      for (j=0; j<3; j++) {
        if (iMin[j] < 0) xmin[0][n][j] = 0.0;
        if (iMax[j] < 0) xmax[0][n][j] = 1.0;
      }
    } /* end loop over variables */
  } /* end loop over data lines */

  /* if pressure is not set, try to read from flow file */
  if (iP < 0 && MDUCgetPressure(0, mem) > 0.0) {
    tmp = MDUCgetPressure(0,mem);
    for (j=0; j<len; j++) P[0][j] = tmp;
  }

  /* check whether some variables should be read as constants from input */
  for (i=0; i<3; i++) {
    if (MDUCgetChiRef(i, mem) >= 0.0) {
      if (iChiR[i] < 0 && !chiR[0]) {
        chiR[0] = (double**) calloc(len, sizeof(double*));
        for (j=0; j<len; j++)
          chiR[0][j] = (double*) calloc(3, sizeof(double));
      }
      tmp = MDUCgetChiRef(i, mem);
      for (j=0; j<len; j++) chiR[0][j][i] = tmp;
    }
    if (MDUCgetZref(i, mem) >= 0.0) {
      if (iZR[i] < 0 && !ZR[0]) {
        ZR[0] = (double**) calloc(len, sizeof(double*));
        for (j=0; j<len; j++)
          ZR[0][j] = (double*) calloc(3, sizeof(double));
      }
      tmp = MDUCgetZref(i, mem);
      for (j=0; j<len; j++) ZR[0][j][i] = tmp;
    }
  }

  fclose(pFile);
  fprintf(stderr, "===================================================\n");
  return len;
}
/* -------------------------------------------------------------------------- */

void writeErrorFileHeader(char *fileName, int dim)
{
  FILE *pFile = fopen(fileName, "w");
  // time
  fprintf(pFile, "# %15s ", "time");
  // pressure
  fprintf(pFile, "%15s %15s %15s ", "P_MDUC[Pa]","P_DNS[Pa]", "P_err[%]");
  // grid bounds
  fprintf(pFile, "%15s %15s %15s %15s %15s %15s ","phi1min_MDUC",
          "phi1min_DNS", "phi1min_err[%]","phi1max_MDUC", "phi1max_DNS",
          "phi1max_err[%]");
  if (dim == 2)
    fprintf(pFile, "%15s %15s %15s %15s %15s %15s ","phi2min_MDUC",
          "phi2min_DNS", "phi2min_err[%]","phi2max_MDUC", "phi2max_DNS",
          "phi2max_err[%]");

  fprintf(pFile, "\n");
 
  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

void writeErrorFile(char *fileName, double t, double P, double P_DNS,
                    double dns1min, double dns1max, double dns2min,
                    double dns2max, void *mem)
{
  double phi1min, phi1max, phi2min, phi2max, err;
  mducStruct  *mduc   = (mducStruct*) mem;
  driveStruct *driver = mduc->driver;
  chemStruct  *chem   = mduc->chem;
  flameStruct *flame  = mduc->flame;
  flowStruct  *flow   = mduc->flow;

  /* Set data row  */
  FILE *pFile = fopen(fileName, "a");

  if (!pFile)
    {
      fprintf(stderr, "Error: could not open error file %s\n", fileName);
      exit(2);
    }

  // by default, set to values in solver (should be correct if mixture fraction)
  phi1min = flow->phi1min;
  phi1max = flow->phi1max;
  phi2min = flow->phi2min;
  phi2max = flow->phi2max;

  // if enthalpy direction, re-compute from the solution and set appropriately
  if (driver->H_DIM > 0)
    {
      calcTotalEnthalpy(flame->n1,flame->n2,flame->nVars,flame->sT,flame->Y,
                        flame->H,chem);

      if (driver->H_DIM == 1)
        {
          phi1min = getMinValue(flame->n1,flame->n2,0,1,flame->H);
          phi1max = getMaxValue(flame->n1,flame->n2,0,1,flame->H);
        }
      else if (driver->H_DIM == 2)
        {
          phi2min = getMinValue(flame->n1,flame->n2,0,1,flame->H);
          phi2max = getMaxValue(flame->n1,flame->n2,0,1,flame->H);
        }
    }

  // write the data line
  fprintf(pFile, "%15.6e ", t); // time
  // pressure
  fprintf(pFile, "%15.6e %15.6e %15.6e ", P, P_DNS, fabs(P-P_DNS)/P_DNS*100);
  // grid bounds
  err = dns1min > 0 ? fabs(phi1min-dns1min)/dns1min*100 : 0.0;
  fprintf(pFile, "%15.6e %15.6e %15.6e ", phi1min, dns1min, err);
  fprintf(pFile, "%15.6e %15.6e %15.6e ", phi1max, dns1max,
          fabs(phi1max-dns1max)/dns1max*100);
  if (flame->nDim == 2)
    {
      err = dns2min > 0 ? fabs(phi2min-dns2min)/dns2min*100 : 0.0;
      fprintf(pFile, "%15.6e %15.6e %15.6e ", phi2min, dns2min, err);
      fprintf(pFile, "%15.6e %15.6e %15.6e ", phi2max, dns2max,
              fabs(phi2max-dns2max)/dns2max*100);
    }
  
  fprintf(pFile, "\n");
  fclose(pFile);
}
/* -------------------------------------------------------------------------- */
