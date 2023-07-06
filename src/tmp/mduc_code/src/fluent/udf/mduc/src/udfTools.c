#include "udfTools.h"
#include "udfPremixed.h"
#if RampantReleaseMajor < 15
#include "shash.h"
#endif /* RampantReleaseMajor < 15 */
#include "accessMDUC.h"

extern udfStruct udf;

#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : fabs(b))

void initFluentOptions(udfStruct *udf)
{
  int i,j;
  Message0("initializing Fluent options... \n");

  /* ----- map the scalars/memory according to input ----- */
  UDFsetScalarIndex();
  UDFsetMemoryIndex();

  /* ----- set up the RIF/statistics co-ordinates ----- */
  /* set defaults */
  udf->iPhi[1]  = udf->iPhi[2]  = -1;
  udf->iPhiV[1] = udf->iPhiV[2] = udf->iPhiV[0] = -1;
  if (udf->iZ1  >= 0) udf->iPhi[1]  = udf->iZ1;  /* Z1 */
  if (udf->iZ1V >= 0) udf->iPhiV[1] = udf->iZ1V; /* Z1-variance */
  if (udf->iZ2  >= 0) udf->iPhi[2]  = udf->iZ2;  /* Z2 */
  if (udf->iZ2V >= 0) udf->iPhiV[2] = udf->iZ2V; /* Z2-variance */
  if (udf->iZ0V >= 0) udf->iPhiV[0] = udf->iZ0V; /* Z0-variance */
  udf->rifEnthalpy = 0; /* default to mixture flamelet */
  udf->rifDim = 1;      /* default to 1D flamelet */

  /* assign according to input */
  if (RP_is_Defined("mduc/rif/dimensions"))
    udf->rifDim = RP_Get_Integer("mduc/rif/dimensions");
  if (RP_is_Defined("mduc/rif/x-coord"))
    checkRIF_coord(RP_Get_String("mduc/rif/x-coord"),
                   &udf->iPhi[1],&udf->iPhiV[1]);
  if (RP_is_Defined("mduc/rif/y-coord"))
    checkRIF_coord(RP_Get_String("mduc/rif/y-coord"),
                   &udf->iPhi[2],&udf->iPhiV[2]);

  Message0("  RIF co-ordinates set to:\n");
  if (udf->iPhi[1] >= 0)
    Message0("    %i (%s)\n", udf->iPhi[1], udf->udsNames[udf->iPhi[1]]);
  if (udf->iPhi[2] >= 0)
    Message0("    %i (%s)\n", udf->iPhi[2], udf->udsNames[udf->iPhi[2]]);

  for (i=0; i<3; i++) udf->phiBuf[i] = 0.0;
  if (RP_is_Defined("mduc/enthalpy-buffer")) {
    if (RP_Get_Real("mduc/enthalpy-buffer") > 0.0) {
      Message0("  Setting enthalpy buffer to +/- %15.6e\n", 
               RP_Get_Real("mduc/enthalpy-buffer"));
      if (udf->iPhi[1] == udf->iH)
        udf->phiBuf[udf->iPhi[1]] = RP_Get_Real("mduc/enthalpy-buffer");
      else if (udf->iPhi[2] == udf->iH)
        udf->phiBuf[udf->iPhi[2]] = RP_Get_Real("mduc/enthalpy-buffer");
      else
        Message0("  Warning: no enthalpy co-ordinate specified\n");
    }
  }

  /* ----- determine which models to employ ----- */
  char *combModel;
  if (RP_is_Defined("mduc/combustion-model"))
    combModel = RP_Get_String("mduc/combustion-model");
  else {
    combModel = (char*) calloc(32, sizeof(char));
    strcpy(combModel, "none");
  }

  udf->rif = udf->premixed = udf->rifPremixed = udf->wiebe = 0;
  if (!strcmp(combModel, "rif"))
    udf->rif = 1;
  else if (!strcmp(combModel, "premixed"))
    udf->premixed = 1;
  else if (!strcmp(combModel, "saci") || !strcmp(combModel, "rif-premixed")) {
    udf->rif = 1;
    udf->rifPremixed = 1;
  }
  else if (!strcmp(combModel, "wiebe"))
    udf->wiebe = 1;
  /* CHANGE MARCO */
  udf->avgMix = 0;
  if (RP_is_Defined("mduc/avg-mixture")) {
    if (udf->rifPremixed && RP_Get_Boolean("mduc/avg-mixture"))
      udf->avgMix = 1;
  }
  else {
    Message0("!!! mduc/avgMixture not defined !!!\n");
  }
  /* check the single droplet model */
  udf->sdm_model = 1;
  if (RP_is_Defined("mduc/sdm-model"))
    udf->sdm_model = RP_Get_Boolean("mduc/sdm-model");

  /* ----- further MDUC information ----- */
  udf->mduc = NULL;
  udf->transform = 1;
  if (RP_is_Defined("mduc/transformation"))
    udf->transform = RP_Get_Integer("mduc/transformation");

  udf->caseDir    = RP_Get_String("mduc/directory");
  udf->filePrefix = RP_Get_String("mduc/file-prefix");

  /* explicitly specify when to initialise injection */
  udf->tPilot = udf->tMain = -1;
  if (RP_is_Defined("mduc/init-pilot-time"))
    udf->tPilot = RP_Get_Real("mduc/init-pilot-time");
  if (RP_is_Defined("mduc/init-main-time"))
    udf->tMain = RP_Get_Real("mduc/init-main-time");

  /* second injection info */
  udf->nZ2 = -1;
  udf->nFuel2 = -1;
  udf->Tfuel2 = -1;
  udf->Yfuel2 = udf->fuel2 = NULL;
  if (RP_is_Defined("mduc/inj2/n-grid"))
    udf->nZ2 = RP_Get_Integer("mduc/inj2/n-grid");
  if (RP_is_Defined("mduc/inj2/n-fuel"))
    udf->nFuel2 = RP_Get_Integer("mduc/inj2/n-fuel");
  if (RP_is_Defined("mduc/inj2/temp"))
    udf->Tfuel2 = RP_Get_Real("mduc/inj2/temp");
  if (RP_is_Defined("mduc/inj2/fuel"))
    readFluentString(udf->nFuel2, RP_Get_String("mduc/inj2/fuel"), &udf->fuel2);
  if (RP_is_Defined("mduc/inj2/mass-fraction"))
    readFluentReal(udf->nFuel2, RP_Get_String("mduc/inj2/mass-fraction"),
                   &udf->Yfuel2);

  /* EGR info */
  udf->nEGR = -1;
  if (RP_is_Defined("mduc/egr-species"))
    udf->nEGR = parseFluentString(RP_Get_String("mduc/egr-species"),
                                  &udf->EGRspecies);

  j = -1;
  char **tmpY = NULL;
  if (udf->nEGR > 0) {
    if (RP_is_Defined("mduc/egr-mass-fractions"))
      j = parseFluentString(RP_Get_String("mduc/egr-mass-fractions"), &tmpY);

    if (j != udf->nEGR)
      Error("EGR species defined with no mass fractions!");
    
    udf->Y_EGR = (double*) calloc(udf->nEGR, sizeof(double));
    Message0("Initializing %i EGR species:\n", udf->nEGR);
    for (i=0; i<udf->nEGR; i++) {
      udf->Y_EGR[i] = atof(tmpY[i]);
      Message0(" %s: %8.4f\n", udf->EGRspecies[i], udf->Y_EGR[i]);
      free(tmpY[i]);
    }
    free(tmpY);
  }

  /* ----- information about threads in case ----- */
  udf->nThreadSolve = -1;
  if (RP_is_Defined("mduc/solve-zones")) {
    udf->nThreadSolve = parseFluentString(RP_Get_String("mduc/solve-zones"),
                                          &udf->solve_thread_names);
    for (i=0; i<udf->nThreadSolve; i++)
      UDFstripStringBraces(udf->solve_thread_names[i]);
  }
  udf->nThreadIntake = -1;
  if (RP_is_Defined("mduc/intake-zones")) {
    udf->nThreadIntake = parseFluentString(RP_Get_String("mduc/intake-zones"),
                                           &udf->intake_thread_names);
    for (i=0; i<udf->nThreadIntake; i++)
      UDFstripStringBraces(udf->intake_thread_names[i]);
  }
  udf->nThreadExhaust = -1;
  if (RP_is_Defined("mduc/exhaust-zones")) {
    udf->nThreadExhaust = parseFluentString(RP_Get_String("mduc/exhaust-zones"),
                                           &udf->exhaust_thread_names);
    for (i=0; i<udf->nThreadExhaust; i++)
      UDFstripStringBraces(udf->exhaust_thread_names[i]);
  }

  /* ----- check monitor ----- */
  char **tmp = NULL;
  udf->monitor = 0;
  if (RP_is_Defined("mduc/monitor/zone-names"))
    j = parseFluentString(RP_Get_String("mduc/monitor/zone-names"), &tmp);
  if (j > 0)
    udf->monitor = 1;

  /* ----- flow variables ----- */
  /* initial pressure rate */
  udf->dPdt = 0.0; /* default */
  if (RP_is_Defined("mduc/initial-pressure-rate"))
    udf->dPdt = RP_Get_Real("mduc/initial-pressure-rate");

  /* ----- spray variables ----- */
  /* liquid drop reference temperature */
  udf->T_drop_ref = 298.15; /* default */
  if (RP_is_Defined("mduc/drop-ref-temp"))
    udf->T_drop_ref = RP_Get_Real("mduc/drop-ref-temp");

  /* get the injection information */
  udf->nNoz = udf->nInj = 0;
  if (RP_is_Defined("spray/nozzles"))
    udf->nNoz = RP_Get_Integer("spray/nozzles");
  if (RP_is_Defined("spray/injection-ids"))
    udf->nInj = parseFluentString(RP_Get_String("spray/injection-ids"),
                                  &udf->inject_id);

  /* remove braces if necessary */
  for (i=0; i<udf->nInj; i++)
    UDFstripStringBraces(udf->inject_id[i]);
  udf->stream = (Injection***) calloc(udf->nInj, sizeof(Injection**));
  for (i=0; i<udf->nInj; i++)
    getInjection(udf->inject_id[i], udf->nNoz, &udf->stream[i]);

  /* ----- scalar dissipation rate settings ----- */
  udf->chiForm1D = udf->chiForm2D = udf->chiMode1D = udf->chiMode2D = NULL;
  udf->chiNorm1D = udf->chiNorm2D = NULL;
  if (RP_is_Defined("mduc/chi-1d-form"))
    udf->chiForm1D = RP_Get_String("mduc/chi-1d-form");
  if (RP_is_Defined("mduc/chi-2d-form"))
    udf->chiForm2D = RP_Get_String("mduc/chi-2d-form");
  if (RP_is_Defined("mduc/chi-1d-mode"))
    udf->chiMode1D = RP_Get_String("mduc/chi-1d-mode");
  if (RP_is_Defined("mduc/chi-2d-mode"))
    udf->chiMode2D = RP_Get_String("mduc/chi-2d-mode");
  if (RP_is_Defined("mduc/chi-1d-scaling"))
    udf->chiNorm1D = RP_Get_String("mduc/chi-1d-scaling");
  if (RP_is_Defined("mduc/chi-2d-scaling"))
    udf->chiNorm2D = RP_Get_String("mduc/chi-2d-scaling");

  /* ----- pdf settings ----- */
  udf->pdf1D = udf->pdf2D = NULL;
  udf->varScale = udf->chiScale = 1.0;
  if (RP_is_Defined("mduc/pdf-1d"))
    udf->pdf1D = RP_Get_String("mduc/pdf-1d");
  if (RP_is_Defined("mduc/pdf-2d"))
    udf->pdf2D = RP_Get_String("mduc/pdf-2d");
  if (RP_is_Defined("mduc/variance-scale"))
      udf->varScale = RP_Get_Real("mduc/variance-scale");
  if (RP_is_Defined("mduc/variance-scale"))
    udf->varScale = RP_Get_Real("mduc/variance-scale");
  if (RP_is_Defined("mduc/dissipation-scale"))
    udf->chiScale = RP_Get_Real("mduc/dissipation-scale");

  /* set defaults if nothing read in */
  if (udf->pdf1D == NULL) {
    udf->pdf1D = (char*) calloc(50, sizeof(char));
    strcpy(udf->pdf1D, "beta");
  }
  if (udf->pdf2D == NULL) {
    udf->pdf2D = (char*) calloc(50, sizeof(char));
    strcpy(udf->pdf2D, "dirichlet");
  }

  /* if enthalpy is a co-ordinate, ensure the settings are correct */
  if (udf->rifEnthalpy == 1) {
    udf->transform = 0;
    strcpy(udf->pdf2D, "independent");
  }

  Message0("done\n");
}
/* ------------------------------------------------------------------------- */

int RP_is_Defined(char *var_name)
{
  int rpDefined = 0;

#if !RP_NODE
  if(NNULLP(seem_hash_lookup(var_name)))
    rpDefined = 1;
#endif
  host_to_node_int_1(rpDefined);

  return rpDefined;
}
/* ------------------------------------------------------------------------- */

void readFluentString(int n, char *input, char ***species)
{
  int i;
  char *name;

  if (n < 1) return;

  species[0] = (char**) calloc(n, sizeof(char*));
  for (i=0; i<n; i++) {
    species[0][i] = (char*) calloc(50, sizeof(char));
    /* get current species name from input list using tokens */
    if (i == 0)
      name = strtok(input, " ");
    else
      name = strtok(NULL, " ");
    strcpy(species[0][i], name);
  }
}
/* ------------------------------------------------------------------------- */

int parseFluentString(char *input, char ***var)
{
  int i, n;
  char sWork[1024], *p;

  if (input == NULL)
    return 0;

  /* copy the input string */
  strcpy(sWork, input);

  var[0] = (char**) calloc(1, sizeof(char*));
  p = strtok(sWork, " ,:");
  i = 0;
  while(p != NULL) {
    var[0]    = realloc(var[0], (i+1)*sizeof(char*));
    var[0][i] = calloc(128, sizeof(char));
    strcpy(var[0][i], p);
    i++;
    p = strtok(NULL, " ,:");
  }
  return i;
}
/* ------------------------------------------------------------------------- */

int parseFluentReal(char *input, double **var)
{
  int i;
  char sWork[1024], tmp[32], *p;

  if (input == NULL)
    return 0;

  /* copy the input string */
  strcpy(sWork, input);

  var[0] = (double*) calloc(1, sizeof(double));
  p = strtok(sWork, " ,:");
  i = 0;
  while(p != NULL) {
    var[0] = realloc(var[0], (i+1)*sizeof(double));
    strcpy(tmp, p);
    var[0][i] = atof(tmp);
    i++;
    p = strtok(NULL, " ,:");
  }
  return i;
}
/* ------------------------------------------------------------------------- */

void UDFstripStringBraces(char *name)
{
  int len;
  char cStart, cEnd, buffer[512];

  len    = strlen(name);
  cStart = name[0];
  cEnd   = name[len-1];

  /* removing leading or trailing brackets*/
  if (cStart == '(') {
    strcpy(buffer, &name[1]);
    strcpy(name, buffer);
    len  = strlen(name);
    cEnd = name[len-1];
  }
  if (cEnd == ')') {
    strncpy(buffer, name, len-1);
    buffer[len-1] = '\0';
    strcpy(name, buffer);
  }
}
/* ------------------------------------------------------------------------- */

void readFluentReal(int n, char *input, double **var)
{
  int i;
  char *name, tmp_name[50];
  float tmp;

  if (n < 1) return;

  var[0] = (double*) calloc(n, sizeof(double));
  for (i=0; i<n; i++) {
    /* get current species name from input list using tokens */
    if (i == 0)
      name = strtok(input, " ");
    else
      name = strtok(NULL, " ");
    strcpy(tmp_name, name);
    sscanf(tmp_name, "%f", &tmp);
    var[0][i] = (double)(tmp);
  }
}
/* ------------------------------------------------------------------------- */

int UDFgetThreads(int n, char **names, Thread **threads)
{
  Domain *d = Get_Domain(1);
  Thread *t;
  int i, iT, nFound;

  /* set all threads to null */
  for (i=0; i<n; i++) threads[i] = NULL;

  /* loop over all threads */
  i = 0;
  nFound = 0;
  thread_loop(t,d,iT) {
    for (i=0; i<n; i++) {
      if (!strcmp(names[i], THREAD_NAME(t))) {
        threads[i] = t;
        nFound++;
      }
    }
  } /* end loop over all threads */
  return nFound;
}
/* -------------------------------------------------------------------------- */

/* ##### start functions that only work for serial/compute node ##### */
#if !RP_HOST
void initThreadSolve(Domain *d)
{
  Message0(" -----  Initializing threads ----- \n");
  Thread *t;
  int i,j,iT;
  
  /* ----- determine which threads to solve ----- */
  if (udf.nThreadSolve < 0) { /* default to all fluid threads */
    Message0("    Defaulting to solve all fluid threads\n");
    /* count the number of fluid cells */
    udf.nThreadSolve = 0;
    thread_loop(t,d,iT) {
      if (FLUID_CELL_THREAD_P(t))
        udf.nThreadSolve++;
    } /* end fluid thread */
    Message0("    Fluid threads found: %i\n", udf.nThreadSolve);
    
    /* allocate the memory for names */
    udf.solve_thread_names = (char**) calloc(udf.nThreadSolve, sizeof(char*));
    
    /* store the thread names */
    i=0; 
    thread_loop(t,d,iT) {
      if (FLUID_CELL_THREAD_P(t)) {
        udf.solve_thread_names[i] = (char*) calloc(128, sizeof(char));
        strcpy(udf.solve_thread_names[i], THREAD_NAME(t));
        i++;
      }
    } /* end fluid thread */
  } /* end default threads */

  /* allocate additional memory */
  udf.solve_thread_ids  = (int*)     calloc(udf.nThreadSolve,  sizeof(int));
  udf.solve_threads     = (Thread**) calloc(udf.nThreadSolve,  sizeof(Thread*));
  udf.intake_thread_ids = (int*)     calloc(udf.nThreadIntake, sizeof(int));
  udf.intake_threads    = (Thread**) calloc(udf.nThreadIntake, sizeof(Thread*));
  udf.exhaust_thread_ids = (int*)    calloc(udf.nThreadExhaust, sizeof(int));
  udf.exhaust_threads = (Thread**)  calloc(udf.nThreadExhaust, sizeof(Thread*));

  for (i=0; i<udf.nThreadSolve; i++) {
    udf.solve_threads[i]    = NULL;
    udf.solve_thread_ids[i] = -1;
  }
  for (i=0; i<udf.nThreadIntake; i++)
    udf.intake_threads[i] = NULL;
  for (i=0; i<udf.nThreadExhaust; i++)
    udf.exhaust_threads[i] = NULL;

  int nFound = 0, nSolve = 0, nIntake = 0, nExhaust = 0;
  thread_loop(t,d,iT) {
    for (i=0; i<udf.nThreadSolve; i++) {
      if (!strcmp(udf.solve_thread_names[i], THREAD_NAME(t))) {
        udf.solve_thread_ids[i] = THREAD_ID(t);
        udf.solve_threads[i] = t;
        nFound++;
        nSolve++;
      }
    } /* end loop over threads to solve */
    for (i=0; i<udf.nThreadIntake; i++) {
      if (!strcmp(udf.intake_thread_names[i], THREAD_NAME(t))) {
        udf.intake_thread_ids[i] = THREAD_ID(t);
        udf.intake_threads[i] = t;
        nFound++;
        nIntake++;
      }
    } /* end loop over intake threads */
    for (i=0; i<udf.nThreadExhaust; i++) {
      if (!strcmp(udf.exhaust_thread_names[i], THREAD_NAME(t))) {
        udf.exhaust_thread_ids[i] = THREAD_ID(t);
        udf.exhaust_threads[i] = t;
        nFound++;
        nExhaust++;
      }
    } /* end loop over exhaust threads */
  } /* end thread loop */
  
  Message0("\n     Total threads found: %i\n", nFound);
  Message0("       Reacting threads (%i/%i) --> ( ", nSolve, udf.nThreadSolve);
  for (i=0; i<udf.nThreadSolve; i++)
    Message0("%s ", udf.solve_thread_names[i]);
  Message0(")\n");
  if (udf.nThreadIntake > 0) {
    Message0("       Intake threads (%i) --> ( ", udf.nThreadIntake);
    for (i=0; i<udf.nThreadIntake; i++)
      Message0("%s ", udf.intake_thread_names[i]);
    Message0(")\n");
  }
  if (udf.nThreadExhaust > 0) {
    Message0("       Exhaust threads (%i) --> ( ", udf.nThreadExhaust);
    for (i=0; i<udf.nThreadExhaust; i++)
      Message0("%s ", udf.exhaust_thread_names[i]);
    Message0(")\n");
  }
  Message0("\n");
  Message0(" -----  Threads successfully initialized ----- \n\n");
}
/* ------------------------------------------------------------------------- */

void checkThreadSolve(Domain *d)
{
  Thread *t;
  int i, iT, nFound;

  Message0("  Updating reacting zones -> ");
  for (i=0; i<udf.nThreadSolve; i++)
    udf.solve_threads[i] = NULL;

  nFound = 0;
  thread_loop(t,d,iT) {
    for (i=0; i<udf.nThreadSolve; i++) {
      if (!strcmp(udf.solve_thread_names[i], THREAD_NAME(t))) {
        udf.solve_thread_ids[i] = THREAD_ID(t);
        udf.solve_threads[i] = t;
        nFound++;
      }
    } /* end loop over threads to solve */
  } /* end thread loop */
  Message0("%i/%i found\n", nFound, udf.nThreadSolve);
}
/* ------------------------------------------------------------------------- */

void UDFupdateSolverThreadNames(Domain *d)
{
  Thread *t;
  int i, iT;
  int n_names, n_ids;
  char **t_names = NULL, **t_ids = NULL;

  Message0("***** Thread info in rpvar *****\n");
  /* get the solver names  */
  if (RP_Get_String("mduc/thread-names")) {
    Message0("   mduc/thread-names: %s\n", RP_Get_String("mduc/thread-names") );
    n_names = parseFluentString(RP_Get_String("mduc/thread-names"), &t_names);
  }
  if (t_names == NULL || n_names == 0)
    Message0("Warning: mduc/thread-names not set\n");

  /* get the solver ids */
  if (RP_Get_String("mduc/thread-ids")) {
    Message0("   mduc/thread-ids: %s\n", RP_Get_String("mduc/thread-ids") );
    n_ids = parseFluentString(RP_Get_String("mduc/thread-ids"), &t_ids);
  }
  if (t_names == NULL || n_names == 0)
    Message0("Warning: mduc/thread-ids not set\n");

  /* check to ensure that names/ids are consistent */
  if (n_ids != n_names) {
    Message0("Warning: number of cylinder names not equal to cylinder ids\n");
    Message0("  ids: %i --> names: %i\n", n_ids, n_names);
  }

  Message0("   matching:\n");
  for (i=0; i<n_ids; i++) {
    UDFstripStringBraces(t_names[i]);
    UDFstripStringBraces(t_ids[i]);
    Message0("   %5i --> %s\n", atoi(t_ids[i]), t_names[i]);
  }
  Message0("********************************\n\n");

  /* rename based on current names */
  thread_loop(t,d,iT) {
    for (i=0; i<n_ids; i++) {
      if (THREAD_ID(t) == atoi(t_ids[i])) {
        Message0("   Renaming thread %i: %s --> %s\n", THREAD_ID(t),
                 THREAD_NAME(t), t_names[i]);
        strcpy(THREAD_NAME(t), t_names[i]);
      }
    } /* end loop over threads */
  } /* end thread loop */

  Message0("\n   Updated thread names:\n");
  thread_loop(t,d,iT) {
    Message0("   Thread %i: %s\n", THREAD_ID(t), THREAD_NAME(t));
  } /* end thread loop */
}
/* ------------------------------------------------------------------------- */

double calcAvgPressure(Domain *d)
{
  Thread *t;
  cell_t c;
        
  double PAvg = 0.0;
  double totVolume = 0.0;

  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      PAvg      += C_VOLUME(c,t)*ABS_P(C_P(c,t), op_pres);
      totVolume += C_VOLUME(c,t);
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

#if RP_NODE /* parallel */
  totVolume = PRF_GRSUM1(totVolume);
  PAvg      = PRF_GRSUM1(PAvg);
#endif /* RP_NODE */

  return PAvg/totVolume;
}
/* ------------------------------------------------------------------------- */

double calcAvgTemp(Domain *d)
{
  int iT;
  Thread *t;
  cell_t c;
        
  double TAvg = 0.0;
  double mass, totMass = 0.0;
  
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;

    begin_c_loop_int(c,t) {
      mass     = C_R(c,t)*C_VOLUME(c,t);
      TAvg    += mass*C_T(c,t);
      totMass += mass;
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

#if RP_NODE /* parallel */
  totMass = PRF_GRSUM1(totMass);
  TAvg    = PRF_GRSUM1(TAvg);
#endif /* RP_NODE */

  return TAvg/totMass;
}
/* ------------------------------------------------------------------------- */

void UDFupdateEnthalpy(Domain *d, void *mem)
{
  int iT;
  Thread *t;
  cell_t c;

  int    i;
  int    nSp = MDUCgetNspecies(mem);
  double *h  = (double*) calloc(nSp, sizeof(double));
  double *cp = (double*) calloc(nSp, sizeof(double));

  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      /* compute the total enthalpy for this cell */
      MDUCcomputeThermoData(h, cp, C_T(c,t), mem);
      C_HTOT(c,t) = 0.0;
      for (i=0; i<udf.nY; i++) {
        if (udf.indexY[i] >= 0)
          C_HTOT(c,t) += C_YI(c,t,i)*h[udf.indexY[i]];
      }
      if (udf.iHV >= 0)
        C_HVAR(c,t) = 0.0;
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

  free(h); free(cp);
}
/* ------------------------------------------------------------------------- */

void clipVariance(Domain *d)
{
  Thread *t;
  cell_t c;

  double var0Max, var1Max, var2Max, var0Min, egrVarMax;

  /* loop over all cells */
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      /* ensure that the means are properly bounded */
      if (udf.iZ1 >= 0) {
        C_Z1MEAN(c,t) = MIN_D(1.0, MAX_D(0.0, C_Z1MEAN(c,t)));
        var1Max = C_Z1MEAN(c,t)*(1.0 - C_Z1MEAN(c,t));
        C_Z1VAR(c,t) = MAX_D(0.0, MIN_D(C_Z1VAR(c,t), var1Max));
      }
      if (udf.iZ2 >= 0) {
        C_Z2MEAN(c,t) = MIN_D(1.0, MAX_D(0.0, C_Z2MEAN(c,t)));
        var2Max = C_Z2MEAN(c,t)*(1.0 - C_Z2MEAN(c,t));
        C_Z2VAR(c,t) = MAX_D(0.0, MIN_D(C_Z2VAR(c,t), var2Max));
      }
      if (udf.iEGR >= 0) {
        C_EGR(c,t) = MIN_D(1.0, MAX_D(0.0, C_EGR(c,t)));
        egrVarMax = C_EGR(c,t)*(1.0 - C_EGR(c,t));
        C_EGRVAR(c,t) = MAX_D(0.0, MIN_D(C_EGRVAR(c,t), egrVarMax));
      }
      if (udf.iProg >= 0) {
        C_PROG(c,t) = MAX_D(0.0, C_PROG(c,t));
      }
      
      if (udf.iZ1 >= 0 && udf.iZ2 >= 0 && udf.iZ0V >= 0) {
        var0Min = C_Z1VAR(c,t)+C_Z2VAR(c,t)-2*sqrt(C_Z1VAR(c,t)*C_Z2VAR(c,t));
        var0Max = C_Z1VAR(c,t)+C_Z2VAR(c,t)+2*sqrt(C_Z1VAR(c,t)*C_Z2VAR(c,t));
        C_Z0VAR(c,t) = MAX_D(var0Min, MIN_D(C_Z0VAR(c,t),  var0Max));
        C_Z0VAR(c,t) = MAX_D(0.0, C_Z0VAR(c,t)); /* ensure > 0 */
      }
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
}
/* ------------------------------------------------------------------------- */

void clipEnthalpyVariance(Domain *d)
{
  Thread *t;
  cell_t c;

  double HvarMax;
  /* loop over all cells */
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      /* make sure the varianc is properly bounded */
      HvarMax     = C_HTOT(c,t)*(udf.Hmin + udf.Hmax - C_HTOT(c,t))
        - udf.Hmin*udf.Hmax;
      C_HVAR(c,t) = MAX_D(0.0, MIN_D(C_HVAR(c,t), HvarMax));
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
}
/* ------------------------------------------------------------------------- */

void getPhiBounds(Domain *d)
{
  Thread *t;
  cell_t c;

  int i, iT;
  double phiMin[3]  = { 1.0e20,  1.0e20,  1.0e20};
  double phiMax[3]  = {-1.0e20, -1.0e20, -1.0e20};
  double phiMinV[3] = { 1.0e20,  1.0e20,  1.0e20};
  double phiMaxV[3] = {-1.0e20, -1.0e20, -1.0e20};
  double Tmin =  1.0e20;
  double Tmax = -1.0e20;
  double rWork[3];
 /* CHANGE MARCO */
  double rho0;
  double rho1;
  rho0 = MDUCgetRhobounds(-1,udf.mduc);
  rho1 = MDUCgetRhobounds(1,udf.mduc);
  /* END CHANGE MARCO */ 
  /* loop over all cells */
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;
    
    begin_c_loop_int(c,t) {
      /* CHANGE MARCO - Do not consider cells that burned premixed */
      if (udf.rifPremixed) {
        /* mean values */
        if (C_PROG(c,t)<0.5) {
          if (udf.iPhi[1] >= 0) {
            if (C_PHI1(c,t) < phiMin[1]) phiMin[1] = C_PHI1(c,t);
            if (C_PHI1(c,t) > phiMax[1]) phiMax[1] = C_PHI1(c,t);
          }
          
          if (udf.iPhi[2] >= 0) {
            if (C_PHI2(c,t) < phiMin[2]) phiMin[2] = C_PHI2(c,t);
            if (C_PHI2(c,t) > phiMax[2]) phiMax[2] = C_PHI2(c,t);
          }
          
          if (udf.iPhi[1] >= 0 && udf.iPhi[2] >= 0) {
            if (C_PHI1(c,t)+C_PHI2(c,t) < phiMin[0])
              phiMin[0] = C_PHI1(c,t)+C_PHI2(c,t);
            if (C_PHI1(c,t)+C_PHI2(c,t) > phiMax[0])
              phiMax[0] = C_PHI1(c,t)+C_PHI2(c,t);
          }
          
          /* variances */
          for (i=0; i<3; i++) {
            if (udf.iPhiV[i] >= 0) {
              if (C_UDSI(c,t,udsOffset+udf.iPhiV[i]) < phiMinV[i])
                phiMinV[i] = C_UDSI(c,t,udsOffset+udf.iPhiV[i]);
              if (C_UDSI(c,t,udsOffset+udf.iPhiV[i]) > phiMaxV[i])
                phiMaxV[i] = C_UDSI(c,t,udsOffset+udf.iPhiV[i]);
            }
          }
          /* keep track of temperature as well */
          if (C_T(c,t) < Tmin) Tmin = C_T(c,t);
          if (C_T(c,t) > Tmax) Tmax = C_T(c,t);
	}
      }
      else {
        /* mean values */
        if (udf.iPhi[1] >= 0) {
          if (C_PHI1(c,t) < phiMin[1]) phiMin[1] = C_PHI1(c,t);
          if (C_PHI1(c,t) > phiMax[1]) phiMax[1] = C_PHI1(c,t);
        }
        
        if (udf.iPhi[2] >= 0) {
          if (C_PHI2(c,t) < phiMin[2]) phiMin[2] = C_PHI2(c,t);
          if (C_PHI2(c,t) > phiMax[2]) phiMax[2] = C_PHI2(c,t);
        }
        
        if (udf.iPhi[1] >= 0 && udf.iPhi[2] >= 0) {
          if (C_PHI1(c,t)+C_PHI2(c,t) < phiMin[0])
            phiMin[0] = C_PHI1(c,t)+C_PHI2(c,t);
          if (C_PHI1(c,t)+C_PHI2(c,t) > phiMax[0])
            phiMax[0] = C_PHI1(c,t)+C_PHI2(c,t);
        }
        
        /* variances */
        for (i=0; i<3; i++) {
          if (udf.iPhiV[i] >= 0) {
            if (C_UDSI(c,t,udsOffset+udf.iPhiV[i]) < phiMinV[i])
              phiMinV[i] = C_UDSI(c,t,udsOffset+udf.iPhiV[i]);
            if (C_UDSI(c,t,udsOffset+udf.iPhiV[i]) > phiMaxV[i])
              phiMaxV[i] = C_UDSI(c,t,udsOffset+udf.iPhiV[i]);
          }
        }
        /* keep track of temperature as well */
      	if (C_T(c,t) < Tmin) Tmin = C_T(c,t);
      	if (C_T(c,t) > Tmax) Tmax = C_T(c,t);
      }

    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

  /* get the min/max across all processors */
#if RP_NODE /* parallel */
  PRF_GRLOW(phiMin,   3, rWork);
  PRF_GRHIGH(phiMax,  3, rWork);
  PRF_GRLOW(phiMinV,  3, rWork);
  PRF_GRHIGH(phiMaxV, 3, rWork);
  Tmin = PRF_GRLOW1(Tmin);
  Tmax = PRF_GRHIGH1(Tmax);
#endif /* RP_NODE */
  
  /* apply buffer if specified */
  for (i=0; i<3; i++) {
    phiMin[i] -= udf.phiBuf[i];
    phiMax[i] += udf.phiBuf[i];
  }
    
  /* assign to global memory */
  for (i=0; i<3; i++) {
    udf.phiMin[i]    = phiMin[i];
    udf.phiMax[i]    = phiMax[i];
    udf.phiMinV[i]   = phiMinV[i];
    udf.phiMaxV[i]   = phiMaxV[i];
    udf.deltaPhi[i]  = phiMax[i] - phiMin[i];
    udf.deltaPhi2[i] = udf.deltaPhi[i]*udf.deltaPhi[i];
  }
  udf.Tmin = Tmin;
  udf.Tmax = Tmax;
  udf.delT = Tmax - Tmin;
  /* CHANGE MARCO */
  /* udf.Hmin & udf.Hmax are used to clip enthalpy variance
     --> only set by getEnthalpyBounds yet, which is only called in udfRIFinit */
  if (udf.rifPremixed) {
    udf.Hmin = udf.phiMin[1];
    udf.Hmax = udf.phiMax[1];
  }
  /* END CHANGE MARCO */
}
/* ------------------------------------------------------------------------- */

void getEnthalpyBounds(Domain *d)
{
  Thread *t;
  cell_t c;

  int    iT;
  double Hmin =  1e20;
  double Hmax = -1e20;
/* CHANGE MARCO */
  double rho0;
  double rho1;
  rho0 = MDUCgetRhobounds(-1,udf.mduc);
  rho1 = MDUCgetRhobounds(1,udf.mduc);
/* END CHANGE MARCO */

  for (iT=0; iT<udf.nThreadSolve; iT++) { 
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;

    begin_c_loop_int(c,t) {
      /* CHANGED MARCO */
      /* DO NOT CONSIDER CELLS, WHICH ARE ALREADY BURN PREMIXED IF RIF*/
      if (udf.rifPremixed) {
      /*	if (C_HTOT(c,t) < Hmin && !C_REGIME(c,t)) Hmin = C_HTOT(c,t);*/
      	if ((C_HTOT(c,t) < Hmin && C_PROG(c,t) < 0.5)) Hmin = C_HTOT(c,t);
      	if ((C_HTOT(c,t) > Hmax && C_PROG(c,t) < 0.5)) Hmax = C_HTOT(c,t);
      /*	if (C_HTOT(c,t) < Hmin && C_PROG(c,t)<0.5) Hmin = C_HTOT(c,t);
      	if (C_HTOT(c,t) > Hmax && C_PROG(c,t)<0.5) Hmax = C_HTOT(c,t);*/
      }
      else {
        if (C_HTOT(c,t) < Hmin) Hmin = C_HTOT(c,t);
        if (C_HTOT(c,t) > Hmax) Hmax = C_HTOT(c,t);
      }
      /* END CHANGE MARCO */
    } end_c_loop_int(c,t); /* end loop over all cells */
  } /* end loop over cell threads */

  /* get the min/max across all processors */
#if RP_NODE /* parallel */
  udf.Hmin = PRF_GRLOW1(Hmin);
  udf.Hmax = PRF_GRHIGH1(Hmax);
#endif /* RP_NODE */

  udf.delH = udf.Hmax - udf.Hmin;
}
/* ------------------------------------------------------------------------- */

void adjustMass(double mSet, Domain *d)
{
  Message0(" Adjusting mass in cylinder to %15.6e\n", mSet);
  Thread *t;
  cell_t c;
        
  int    iT, nCells = 0;
  double delM, totMass = 0.0, totVol = 0.0;

  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL)
      continue;
    
    begin_c_loop_int(c,t) {
      totVol  += C_VOLUME(c,t);
      totMass += C_R(c,t)*C_VOLUME(c,t);
      nCells++;
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

#if RP_NODE
  /* Get the total average quantities from all compute nodes */
  totVol  = PRF_GRSUM1(totVol);
  totMass = PRF_GRSUM1(totMass);
  nCells  = PRF_GISUM1(nCells);
#endif /* RP_NODE */

  /* compute the required mass adjustment */
  delM = mSet - totMass;
  Message0(" Mass adjustment: %15.6e (%6.4f \%)\n", delM,
           ABS_D(delM)/totMass*100.0);
  
  /* update the density by removing mass over entire volume */
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL)
      continue;
    
    begin_c_loop_int(c,t) {
      C_R(c,t) += delM/totVol;
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
}
/* ------------------------------------------------------------------------- */

void adjustTemp(double Tset, Domain *d)
{
  Message0("\n Adjusting average temperature in cylinder to %6.1f\n", Tset);
  Thread *t;
  cell_t c;
        
  int    iT, nCells = 0;
  double delT, mass, totMass = 0.0, TFav = 0.0;

  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL)
      continue;
    
    begin_c_loop_int(c,t) {
      mass     = C_R(c,t)*C_VOLUME(c,t);
      TFav    += mass*C_T(c,t);
      totMass += mass;
      nCells++;
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */

#if RP_NODE
  /* Get the total average quantities from all compute nodes */
  TFav    = PRF_GRSUM1(TFav);
  totMass = PRF_GRSUM1(totMass);
  nCells  = PRF_GISUM1(nCells);
#endif /* RP_NODE */
  TFav   /= totMass;
  
  /* compute the required mass adjustment */
  delT = Tset - TFav;
  Message0(" Temp adjustment: %8.6f (%6.4f \%)\n",delT,ABS_D(delT)/TFav*100.0);
  
  /* update the density by removing mass over entire volume */
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL)
      continue;
    
    begin_c_loop_int(c,t) {
      C_T(c,t) += delT;/**C_R(c,t)*C_VOLUME(c,t)/totMass;*/
    } end_c_loop_int(c,t);
  } /* end loop over cell threads */
}
/* ------------------------------------------------------------------------- */

void UDFwritePDF_2D(Domain *d, int local)
{
  Thread *t;
  cell_t c;

  int    i, j, k, iT, iX, iY, iZ;
  double TINY = 1.0e-12;
  double mass, totMass, tmp;
  double x, y, z, dx, dy, dz;

  /* allocate the binning grid */
  int n1 = 101, n2 = 101, nT = 101;
  double *x1    = (double*)  calloc(n1, sizeof(double));
  double *x2    = (double*)  calloc(n2, sizeof(double));
  double *xT    = (double*)  calloc(nT, sizeof(double));
  double **pdf  = (double**) calloc(n2, sizeof(double*));
  double **pTZ1 = (double*)  calloc(nT, sizeof(double));
  double **pTZ2 = (double*)  calloc(nT, sizeof(double));
  double *pdfZ1 = (double*)  calloc(n1, sizeof(double));
  double *pdfZ2 = (double*)  calloc(n2, sizeof(double));
  double *pdfT  = (double*)  calloc(nT, sizeof(double));
  double **T    = (double**) calloc(n2, sizeof(double*));
  double **Tmin = (double**) calloc(n2, sizeof(double*));
  double **Tmax = (double**) calloc(n2, sizeof(double*));
  double *T1    = (double*)  calloc(n1, sizeof(double));
  double *T2    = (double*)  calloc(n2, sizeof(double));
  double *T1min = (double*)  calloc(n1, sizeof(double));
  double *T1max = (double*)  calloc(n1, sizeof(double));
  double *T2min = (double*)  calloc(n2, sizeof(double));
  double *T2max = (double*)  calloc(n2, sizeof(double));
  double **H    = (double**) calloc(n2, sizeof(double*));
  double **Hmin = (double**) calloc(n2, sizeof(double*));
  double **Hmax = (double**) calloc(n2, sizeof(double*));
  double *H1    = (double*)  calloc(n1, sizeof(double));
  double *H2    = (double*)  calloc(n2, sizeof(double));
  double *H1min = (double*)  calloc(n1, sizeof(double));
  double *H1max = (double*)  calloc(n1, sizeof(double));
  double *H2min = (double*)  calloc(n2, sizeof(double));
  double *H2max = (double*)  calloc(n2, sizeof(double));
  double *rwork = (double*)  calloc(n1, sizeof(double));
  for (k=0; k<n2; k++) {
    pdf[k]  = (double*) calloc(n1, sizeof(double));
    pTZ1[k] = (double*) calloc(n1, sizeof(double));
    pTZ2[k] = (double*) calloc(n2, sizeof(double));
    T[k]    = (double*) calloc(n1, sizeof(double));
    Tmin[k] = (double*) calloc(n1, sizeof(double));
    Tmax[k] = (double*) calloc(n1, sizeof(double));
    H[k]    = (double*) calloc(n1, sizeof(double));
    Hmin[k] = (double*) calloc(n1, sizeof(double));
    Hmax[k] = (double*) calloc(n1, sizeof(double));
  }

  /* fill the grid */
  double dx1 = 1.0/(n1-1);
  double dx2 = 1.0/(n2-1);
  double dT  = 1.0/(nT-1);
  x1[0] = x2[0] = xT[0] = 0.0;
  for (j=1; j<n1; j++) x1[j] = x1[j-1] + dx1;
  for (j=1; j<n2; j++) x2[j] = x2[j-1] + dx2;
  for (j=1; j<n2; j++) xT[j] = xT[j-1] + dT;

  /* file variables */
  char  pdfFileName[1024], suffix[64];
  FILE *pFile;

  /* make sure everything is set to zero*/
  for (k=0; k<n2; k++) {
    for (j=0; j<n1; j++) {
      pdf[k][j] = T[k][j] = H[k][j] = 0.0;
      Tmin[k][j] = Hmin[k][j] =  1.0e20;
      Tmax[k][j] = Hmax[k][j] = -1.0e20;
    }
  }

  /* ----- compute the means ----- */
  totMass = 0.0;
  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;
    
    begin_c_loop_int(c,t) {
      mass = C_R(c,t)*C_VOLUME(c,t);     /* cell mass */
      /* only consider a PDF over the cells with non-zero mixture */
      if ( (C_PHI1(c,t) > 0.0) || (C_PHI2(c,t) > 0.0) ) {
        x = (C_PHI1(c,t)  - udf.phiMin[1])/udf.deltaPhi[1];
        y = (C_PHI2(c,t)  - udf.phiMin[2])/udf.deltaPhi[2];
        z = (C_T(c,t)     - udf.Tmin)/udf.delT;
        
        /* search for the bin that it fits in */
        j = UDFgridLocate(n1, x1, x);
        k = UDFgridLocate(n2, x2, y);
        i = UDFgridLocate(nT, xT, z);
        
        if (x < 0.5*(x1[j]+x1[j+1])) {
          iX = j;
          if (y < 0.5*(x2[k]+x2[k+1]))
            iY = k;
          else
            iY = k+1;
        }
        else {
          iX = j+1;
          if (y < 0.5*(x2[k]+x2[k+1]))
            iY = k;
          else
            iY = k+1;
        }

        if (z < 0.5*(xT[i]+xT[i+1]))
          iZ = i;
        else
          iZ = i+1;

        /* assign the variables */
        pdf[iY][iX]  += mass;
        pTZ1[iZ][iX] += mass;
        pTZ2[iZ][iY] += mass;

        /* temperature */
        T[iY][iX] += C_T(c,t)*mass;
        if (C_T(c,t) < Tmin[iY][iX]) Tmin[iY][iX] = C_T(c,t);
        if (C_T(c,t) > Tmax[iY][iX]) Tmax[iY][iX] = C_T(c,t);
        if (udf.iH >= 0) {
          H[iY][iX] += C_HTOT(c,t)*mass;
          if (C_HTOT(c,t) < Hmin[iY][iX]) Hmin[iY][iX] = C_HTOT(c,t);
          if (C_HTOT(c,t) > Hmax[iY][iX]) Hmax[iY][iX] = C_HTOT(c,t);
        }

      } /* end pdf computation */
      totMass += mass;
    } end_c_loop_int(c,t);
  }  /* end loop over cell threads */

#if RP_NODE
  /* if parallel, get the contribution to from each node */
  for (k=0; k<n2; k++) {
    PRF_GRSUM(pdf[k], n1, rwork);
    PRF_GRSUM(T[k], n1, rwork);
    PRF_GRSUM(H[k], n1, rwork);
    PRF_GRLOW(Tmin[k], n1, rwork);
    PRF_GRHIGH(Tmax[k], n1, rwork);
  }
  for (k=0; k<nT; k++) {
    PRF_GRSUM(pTZ1[k], n1, rwork);
    PRF_GRSUM(pTZ2[k], n2, rwork);
  } 
  totMass = PRF_GRSUM1(totMass);
#endif /* RP_NODE */

  /* normalize the pdf itself */
  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      pdf[k][j] /= totMass;

  for (i=0; i<nT; i++) {
    for (j=0; j<n1; j++) pTZ1[i][j] /= totMass;
    for (k=0; k<n2; k++) pTZ2[i][k] /= totMass;
  }
    
  /* compute the marginal pdfs */  
  for (j=0; j<n1; j++) {
    pdfZ1[j] = T1[j] = H1[j] = 0.0;
    T1min[j] = H1min[j] =  1.0e20;
    T1max[j] = H1max[j] = -1.0e20;
  }
  for (k=0; k<n2; k++) {
    pdfZ2[k] = T2[k] = H2[k] = 0.0;
    T2min[k] = H2min[k] =  1.0e20;
    T2max[k] = H2max[k] = -1.0e20;
  }
  for (k=0; k<n2; k++) {
    for (j=0; j<n1; j++) {
      pdfZ1[j] += pdf[k][j];
      pdfZ2[k] += pdf[k][j];
      T1[j]    += T[k][j];
      T2[k]    += T[k][j];
      H1[j]    += H[k][j];
      H2[k]    += H[k][j];

      if (Tmin[k][j] < T1min[j]) T1min[j] = Tmin[k][j];
      if (Tmin[k][j] < T2min[k]) T2min[k] = Tmin[k][j];
      if (Tmax[k][j] > T1max[j]) T1max[j] = Tmax[k][j];
      if (Tmax[k][j] > T2max[k]) T2max[k] = Tmax[k][j];

      if (Hmin[k][j] < H1min[j]) H1min[j] = Hmin[k][j];
      if (Hmin[k][j] < H2min[k]) H2min[k] = Hmin[k][j];
      if (Hmax[k][j] > H1max[j]) H1max[j] = Hmax[k][j];
      if (Hmax[k][j] > H2max[k]) H2max[k] = Hmax[k][j];
    } /* end loop over columns */
  } /* end loop over rows */

  for (i=0; i<nT; i++) {
    for (j=0; j<n1; j++)
      pdfT[i] += pTZ1[i][j];
  }

  /* normalize the conditional variables */
  for (k=0; k<n2; k++) {
    for (j=0; j<n1; j++) {
      if (pdf[k][j] > 0.0) {
        T[k][j] /= pdf[k][j]*totMass;
        H[k][j] /= pdf[k][j]*totMass;
      }
      else {
        Tmin[k][j] = Tmax[k][j] = 0.0;
        Hmin[k][j] = Hmax[k][j] = 0.0;
      }
    } /* end loop over columns */
  } /* end loop over rows */

  for (j=0; j<n1; j++) {
    if (pdfZ1[j] > 0.0) {
      T1[j] /= pdfZ1[j]*totMass;
      H1[j] /= pdfZ1[j]*totMass;
    }
    else {
      T1min[j] = T1max[j] = H1min[j] = H1max[j] = 0.0;
    }
  }
  for (k=0; k<n2; k++) {
    if (pdfZ2[k] > 0.0) {
      T2[k] /= pdfZ2[k]*totMass;
      H2[k] /= pdfZ2[k]*totMass;
    }
    else {
      T2min[k] = T2max[k] = H2min[k] = H2max[k] = 0.0;
    }
  }

  /* ----- compute the variances ----- */
  double **Tv = (double**) calloc(n2, sizeof(double*));
  double *T1v = (double*)  calloc(n1, sizeof(double));
  double *T2v = (double*)  calloc(n2, sizeof(double));
  double **Hv = (double**) calloc(n2, sizeof(double*));
  double *H1v = (double*)  calloc(n1, sizeof(double));
  double *H2v = (double*)  calloc(n2, sizeof(double));
  for (k=0; k<n2; k++) {
    Tv[k] = (double*) calloc(n1, sizeof(double));
    Hv[k] = (double*) calloc(n1, sizeof(double));
  }

  for (k=0; k<n2; k++)
    for (j=0; j<n1; j++)
      Tv[k][j] = Hv[k][j] = 0.0;

  for (iT=0; iT<udf.nThreadSolve; iT++) {
    t = udf.solve_threads[iT]; /* current thread */
    if (t == NULL) continue;
    
    begin_c_loop_int(c,t) {
      mass = C_R(c,t)*C_VOLUME(c,t);     /* cell mass */
      /* only consider a PDF over the cells with non-zero mixture */
      if ( (C_PHI1(c,t) > 0.0) || (C_PHI2(c,t) > 0.0) ) {
        x = (C_PHI1(c,t)  - udf.phiMin[1])/udf.deltaPhi[1];
        y = (C_PHI2(c,t)  - udf.phiMin[2])/udf.deltaPhi[2];
        
        /* search for the bin that it fits in */
        j = UDFgridLocate(n1, x1, x);
        k = UDFgridLocate(n2, x2, y);
        
        if (x < 0.5*(x1[j]+x1[j+1])) {
          iX = j;
          if (y < 0.5*(x2[k]+x2[k+1]))
            iY = k;
          else
            iY = k+1;
        }
        else {
          iX = j+1;
          if (y < 0.5*(x2[k]+x2[k+1]))
            iY = k;
          else
            iY = k+1;
        }

        /* compute the variance */
        Tv[iY][iX] += SQR(C_T(c,t)-T[iY][iX])*mass;
        if (udf.iH >= 0)
          Hv[iY][iX] += SQR(C_HTOT(c,t)-H[iY][iX])*mass;

      } /* end pdf computation */
    } end_c_loop_int(c,t);
  }  /* end loop over cell threads */
#if RP_NODE
  /* if parallel, get the contribution to from each node */
  for (k=0; k<n2; k++) {
    PRF_GRSUM(Tv[k], n1, rwork);
    PRF_GRSUM(Hv[k], n1, rwork);
  }
#endif /* RP_NODE */

  /* compute the marginal pdfs */  
  for (j=0; j<n1; j++) T1v[j] = H1v[j] = 0.0;
  for (k=0; k<n2; k++) T2v[k] = H2v[k] = 0.0;
  for (k=0; k<n2; k++) {
    for (j=0; j<n1; j++) {
      T1v[j] += Tv[k][j];
      T2v[k] += Tv[k][j];
      H1v[j] += Hv[k][j];
      H2v[k] += Hv[k][j];
    } /* end loop over columns */
  } /* end loop over rows */

  /* normalize the conditional variance */
  for (k=0; k<n2; k++) {
    for (j=0; j<n1; j++) {
      if (pdf[k][j] > 0.0) {
        Tv[k][j] /= pdf[k][j]*totMass;
        Hv[k][j] /= pdf[k][j]*totMass;
      }
    } /* end loop over columns */
  } /* end loop over rows */

  /* marginal 1 */
  for (j=0; j<n1; j++) {
    if (pdfZ1[j] > 0.0) {
      T1v[j] /= pdfZ1[j]*totMass;
      H1v[j] /= pdfZ1[j]*totMass;
    }
  }
  /* marginal 2 */
  for (k=0; k<n2; k++) {
    if (pdfZ2[k] > 0.0) {
      T2v[k] /= pdfZ2[k]*totMass;
      H2v[k] /= pdfZ2[k]*totMass;
    }
  }

  /* ----- write the results to a file ----- */
#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      /* generate the file suffix */
      if (in_cylinder_p) {
        sprintf(suffix,"%08.3fCAD",RP_Get_Real("dynamesh/in-cyn/crank-angle"));
      }
      else
        sprintf(suffix, "%05i", N_TIME);

      /* ----- write the joint pdf ----- */
      /* setup the file name */
      if (local) {
        sprintf(pdfFileName, "pdf_%s.dat", suffix);
        Message0(" writing pdf to file: %s\n", pdfFileName);
      }
      else
        sprintf(pdfFileName, "%s/pdf/pdf-%s.dat", udf.caseDir, suffix);
      pFile = fopen(pdfFileName, "w");
      
      /* write the header */
      fprintf(pFile,"# phi1 %i %15.6e %15.6e\n",n1,udf.phiMin[1],udf.phiMax[1]);
      fprintf(pFile,"# phi2 %i %15.6e %15.6e\n",n2,udf.phiMin[2],udf.phiMax[2]);
      fprintf(pFile,"# Temp %i %15.6e %15.6e\n",nT,udf.Tmin,udf.Tmax);
      fprintf(pFile,"# phi1<1> phi2<2> pdf_1_2<3> pdf_T_1<4> pdf_T_2<5> ");
      fprintf(pFile,"T[K]<6> Tvar[K]<7> Tmin[K]<8> Tmax[K]<9> ");
      if (udf.iH >= 0)
        fprintf(pFile, "H<10> Hvar<11> Hmin<12> Hmax<13> ");
      fprintf(pFile, "\n");

      /* write the 2D data */
      for (k=0; k<n2; k++) {
        for (j=0; j<n1; j++) {
          fprintf(pFile, "%15.6e %15.6e %15.6e", x1[j], x2[k], pdf[k][j]);
          fprintf(pFile, "%15.6e %15.6e", pTZ1[k][j], pTZ2[k][j]);
          fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e", T[k][j], Tv[k][j], 
                  Tmin[k][j], Tmax[k][j]);
          if (udf.iH >= 0)
            fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e", H[k][j], Hv[k][j],
                    Hmin[k][j], Hmax[k][j]);
          fprintf(pFile, "\n");
        } /* end loop over columns */
        fprintf(pFile, "\n");
      } /* end loop over rows */
      fclose(pFile);

      /* ----- write the marginal pdfs ----- */
      if (local) {
        sprintf(pdfFileName, "pdf-marg-%s.dat", suffix);
        Message0(" writing marginal to file: %s\n", pdfFileName);
      }
      else
        sprintf(pdfFileName, "%s/pdf/pdf-marg-%s.dat", udf.caseDir, suffix);
      pFile = fopen(pdfFileName, "w");

      /* marginal for phi1 */
      fprintf(pFile, "# phi1/phi2<1> pdf<2> ");
      fprintf(pFile, "T[K]<3> Tvar[K^2]<4> Tmin[K]<5> Tmax[K]<6> ");
      if (udf.iH >= 0) 
        fprintf(pFile, "H<7> Hvar<8> Hmin<9> Hmax<10> ");
      fprintf(pFile, "\n");
      for (j=0; j<n1; j++) {
        fprintf(pFile, "%15.6e %15.6e ", x1[j], pdfZ1[j]);
        fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e ", T1[j], T1v[j],
                T1min[j], T1max[j]);
        if (udf.iH >= 0)
          fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e ", H1[j], H1v[j],
                  H1min[j], H1max[j]);
        fprintf(pFile, "\n");
      }

      /* marginal for phi2 */
      fprintf(pFile, "\n\n");
      for (k=0; k<n2; k++) {
        fprintf(pFile, "%15.6e %15.6e ", x2[k], pdfZ2[k]);
        fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e ", T2[k], T2v[k],
                T2min[k], T2max[k]);
        if (udf.iH >= 0)
          fprintf(pFile, "%15.6e %15.6e %15.6e %15.6e ", H2[k], H2v[k],
                  H2min[k], H2max[k]);
        fprintf(pFile, "\n");
      }

      fprintf(pFile, "\n\n");
      for (i=0; i<nT; i++) {
        fprintf(pFile, "%15.6e %15.6e ", xT[i], pdfT[i]);
        fprintf(pFile, "\n");
      }
      fclose(pFile);
    }

  free(x1); free(x2);
  for (k=0; k < n2; k++) {
    free(pdf[k]);
    free(T[k]); free(Tmin[k]); free(Tmax[k]); free(Tv[k]);
    free(H[k]); free(Hmin[k]); free(Hmax[k]); 
  }
  for (i=0; i<nT; i++) {
    free(pTZ1[i]); free(pTZ2[i]);
  }
  free(pdf); free(pdfZ1); free(pdfZ2); free(pdfT);
  free(pTZ1); free(pTZ2); free(rwork);
  free(T); free(T1); free(T2);
  free(Tmin); free(Tmax); free(T1min); free(T1max); free(T2min); free(T2max);
  free(Tv); free(T1v); free(T2v);
  free(H); free(H1); free(H2);
  free(Hmin); free(Hmax); free(H1min); free(H1max); free(H2min); free(H2max);
  free(Hv); free(H1v); free(H2v);
  Message0(" Exiting write pdf\n");
}
/* ------------------------------------------------------------------------- */
#endif /* !RP_HOST */
/* ##### end functions that only work for serial/compute node ##### */

int calcFactorial(int x)
{
  int i, fact;

  fact = x;
  for (i=x-1; i > 1; --i)
    fact *= i;

  return fact;
}
/* ------------------------------------------------------------------------- */

void write1DArray(char *fileName,char *header,int n1,double *x1,double *var)
{
  int j;
  FILE *pFile = fopen(fileName, "w");

  fprintf(pFile, "%s\n", header);

  for (j=0; j<n1; j++) {
    fprintf(pFile, "%8.6f %15.6e\n", x1[j], var[j]);
  }
  
  fclose(pFile);
}
/* ------------------------------------------------------------------------- */

void write2DArray(char *fileName, char *header, int n1, int n2, double *x1, 
                  double *x2, double **var)
{
  int j, k;
  FILE *pFile = fopen(fileName, "w");

  fprintf(pFile, "%s\n", header);

  for (k=0; k<n2; k++) {
    for (j=0; j<n1; j++) {
      fprintf(pFile, "%8.6f %8.6f %15.6e\n", x1[j], x2[k], var[k][j]);
    }
    fprintf(pFile, "\n");
  }
  
  fclose(pFile);
}
/* ------------------------------------------------------------------------- */

int rootBrac(double (*func)(double, void *), double *x1, double *x2, void *args)
{
  int j, NTRY=50;
  double f1, f2, FACTOR=1.6;

  if (*x1 == *x2) {
    fprintf(stderr, "bad initial range in brac\n");
    exit(0);
  }

  f1 = (*func)(*x1, args);
  f2 = (*func)(*x2, args);

  for (j=1; j<=NTRY; j++) {
    if (f1*f2 < 0.0) return 1;
    if (fabs(f1) < fabs(f2))
      f1 = (*func)(*x1 += FACTOR*(*x1-*x2), args);
    else
      f2 = (*func)(*x2 += FACTOR*(*x2-*x1), args);
  }

  return 0;
}
/* -------------------------------------------------------------------------- */

int rootBrak(double (*func)(double, void *), double x1, double x2, int n, 
             double *xb1, double *xb2, int *nb, void *args)
{
  int nbb, i;
  double x, fp, fc, dx;

  nbb = 0;
  dx = (x2-x1)/n;
  fp = (*func)(x=x1, args);
  for (i=1; i<=n; i++) {
    fc = (*func)(x += dx, args);
    if (fc*fp <= 0.0) {
      xb1[++nbb] = x-dx;
      xb2[nbb] = x;
      if (*nb == nbb) return;
    }
    fp = fc;
  }
  *nb = nbb;
}
/* -------------------------------------------------------------------------- */

double rootBisect(double (*func)(double, void *), double x1, double x2,
                  double tol, void *args)
{
  int j, JMAX=100;
  double dx, f, fmid, xmid, rtb;

  f    = (*func)(x1, args);
  fmid = (*func)(x2, args);
  if (f*fmid >= 0.0) {
    fprintf(stderr, "root must be bracketed in bisection\n");
    exit(0);
  }

  /* orient the search so that f>0 lies at x+dx */
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);

  for (j=1; j<=JMAX; j++) { /* bisection loop */
    fmid = (*func)(xmid=rtb+(dx *= 0.5), args);
    if (fmid <= 0.0) rtb = xmid;
    if (fabs(dx) < tol || fmid == 0.0) return rtb;
  }
  fprintf(stderr, "Maximum number of iterations in bisect exceeded\n");
  exit(0);
  return 0.0;
}
/* -------------------------------------------------------------------------- */

double rootBrent(double (*funcEval)(double, void *), double x1, double x2,
                 double tol, void *args)
{
  int j, JMAX=100;
  double EPS=3.0e-8;
  double a=x1, b=x2, c=x2, d,e,min1,min2;
  double fa=(*funcEval)(a,args),fb=(*funcEval)(b,args),fc,p,q,r,s,tol1,xm;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
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
      if (p > 0.0) q=-q; /* check whether in bounds */
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

void indexIntArray(int n, int *arr, int *indx)
{
  int M = 7;
  int NSTACK = 64;
  int i, indxt, ir=n-1, itemp, j, k, l=0;
  int jstack=-1, *istack;
  int a;

  istack = (int*) calloc(NSTACK, sizeof(int));

  for (j=0; j<n; j++) indx[j] = j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1; j<=ir; j++) {
        indxt = indx[j];
        a = arr[indxt];
        for (i=j-1; i>=l; i--) {
          if (arr[indx[i]] <= a) break;
          indx[i+1] = indx[i];
        }
        indx[i+1] = indxt;
      }
      if (jstack < 0) break;
      ir = istack[jstack--];
      l  = istack[jstack--];
    }
    else {
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1]);
      /* swapInt(&indx[k],&indx[l+1]); */
      if (arr[indx[l]] > arr[indx[ir]]) {
        SWAP(indx[l],indx[ir]);
        /* swapInt(&indx[l],&indx[ir]); */
      }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
        SWAP(indx[l+1],indx[ir]);
        /* swapInt(&indx[l+1],&indx[ir]); */
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
        SWAP(indx[l],indx[l+1]);
        /* swapInt(&indx[l],&indx[l+1]); */
      }
      i = l+1;
      j = ir;
      indxt = indx[l+1];
      a = arr[indxt];
      for (;;) {
        do i++; while (arr[indx[i]] < a);
        do j--; while (arr[indx[j]] > a);
        if (j < i) break;
        SWAP(indx[i],indx[j]);
        /* swapInt(&indx[i],&indx[j]); */
      }
      indx[l+1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (jstack >= NSTACK) {
        fprintf(stderr, "NSTACK too small in index.\n");
        exit(2);
      }
      if (ir-i+1 >= j-l) {
        istack[jstack]   = ir;
        istack[jstack-1] = i;
        ir = j-1;
      }
      else {
        istack[jstack] = j-1;
        istack[jstack-1] = l;
        l = i;
      }
    }
  }
  free(istack);
}
/* -------------------------------------------------------------------------- */

int UDFgridLocate(int n1, double *x1, double x)
{
  int jL = 0, jR = n1-1, jM;
  int ascnd = (x1[n1-1] >= x1[0]);  /* whether grid is ascending or not */

  /* use bisection to find nearest grid point */
  while (jR - jL > 1) {
    /* compute a midpoint */
    jM = (jL+jR)/2;
    if (x >= x1[jM] == ascnd)
      jL = jM;
    else
      jR = jM;
  }
  return MAX_I(0, MIN_I(n1-2, jL));
}
/* -------------------------------------------------------------------------- */

double MAX_D(double x, double y)
{
  if (x > y) return x;
  else return y;
}
/* ------------------------------------------------------------------------- */

double MIN_D(double x, double y)
{
  if (x < y) return x;
  else return y;
}
/* ------------------------------------------------------------------------- */

double ABS_D(double x)
{
  if (x < 0.0) x = -x;
  return x;
}
/* ------------------------------------------------------------------------- */

int MAX_I(int x, int y)
{
  if (x > y) return x;
  else return y;
}
/* ------------------------------------------------------------------------- */

int MIN_I(int x, int y)
{
  if (x < y) return x;
  else return y;
}
/* ------------------------------------------------------------------------- */

double SGN_D(double x)
{
  if (x > 0.0) return 1.0;
  else return -1.0;
}
/* ------------------------------------------------------------------------- */
