#include "exportUDF.h"
#include "ensightIO.h"
#include "vtkIO.h"

expStruct ex;

DEFINE_ON_DEMAND(init)
{
  int i, j, line;

  /* initialize the thread info */
  ex.nCells = initThreads(&ex.cell_th, &ex.part_id, &ex.part_name,
                          &ex.part_type, &ex.part_nv, &ex.part_verts);

  /* initialize the variables to export */
  UDFinitExport();

  /* display summary*/
  Message0("\n ================== Custom Export Summary =================\n\n");
  if (ex.binary)
    Message0("  Format:             binary\n");
  else
    Message0("  Format:             ASCII\n");
  if (ex.byteOrder == 0)
    Message0("  Native Byte Ord:    Little Endian\n");
  else if (ex.byteOrder == 1)
    Message0("  Native Byte Ord:    Big Endian\n");
  Message0("  Export directory:   %s\n", ex.dir);
  Message0("  File prefix:        %s\n", ex.prefix);
  Message0("  Frequency:          %i\n", ex.freq);

  Message0("  Fluid Zones:        %i\n%22s", ex.nCells, " ");
  for (i=0; i<ex.nCells; i++) {
    Message0("%s ( ", ex.part_name[i]);
    for (j=0; j<ex.part_nv[i]; j++)
      Message0("%s ", ex.part_type[i][j]);
    Message0(")\n%22s", " ");
  }
  Message0("\n");

  Message0("  Specified Zones:    %i\n%22s", ex.nZone, " ");
  for (i=0; i<ex.nZone; i++) {
    Message0("%s ( ", ex.zones[i]);
    for (j=0; j<ex.nZoneThreads[i]; j++)
      Message0("%s ", ex.zThreadNames[i][j]);
    Message0(")\n%22s", " ");
  }
  Message0("\n");

  Message0("  Specified Walls:    %i\n%22s", ex.nSurf, " ");
  for (i=0; i<ex.nSurf; i++) {
    Message0("%s ( ", ex.surfs[i]);
    for (j=0; j<ex.nSurfThreads[i]; j++)
      Message0("%s ", ex.sThreadNames[i][j]);
    Message0(")\n%22s", " ");
  }
  Message0("\n");

  Message0("  Specified Surfaces:    %i\n%22s", ex.nCxSurf, " ");
  for (i=0; i<ex.nCxSurf; i++) {
    Message0("%s ( ", ex.cxSurfs[i]);
    for (j=0; j<ex.nCxSurfThreads[i]; j++)
      Message0("%s ", ex.cxThreadNames[i][j]);
    Message0(")\n%22s", " ");
  }
  Message0("\n");


  Message0("  Standard variables: ");
  for (i=0, line=1; i<ex.nSc; i++, line++) {
    Message0("%s ", ex.sc_names[i]);
    if (line%5 == 0)
      Message0("\n%22s", " ");
  }
  Message0("\n");
  if (ex.dump_species) {
    Message0("  Species:            ");
    for (i=0; i<ex.nY; i++)
      Message0("%s ", ex.Y_names[i]);
    Message0("\n");
  }
  if (ex.dump_particle) {
    Message0("  Particle data:      \n");
    Message0("     streams:         %i -> ", ex.nInj);
    for (i=0; i<ex.nInj; i++)
      Message0("%s ", ex.inject_id[i]);
    Message0("\n");
    Message0("     nozzles:         %i\n", ex.nNoz);
  }
  Message0("  SDM data:           %i\n", ex.dump_sdm);
#if PARALLEL
#if UDF_MPI_WRITE
  Message0("  Parallel write:     MPI (needs testing)\n");
#else /* ROOT WRITE */
  Message0("  Parallel write:     root process\n");
#endif /* MDUC_MPI_WRITE */
#endif /* PARALLEL */
  Message0("\n ==========================================================\n\n");
}
/* -------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(writeData)
{
#if !RP_HOST
  Message0("\n =================== Custom Data Export ===================\n\n");

  /* check to see if number of fluid zones has changed */
  checkCurrentZones(0);

  Message0("  Exporting data to: %s\n", ex.dir);
  if (ex.format == 0)
    UDFexportEnsight(ex.dir, ex.prefix);
  else if (ex.format == 1)
    VTKexport(ex.dir, ex.prefix);

  Message0("\n ==========================================================\n\n");
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(checkZones)
{
#if !RP_HOST
  Message0("\n ================== Checking Export Zones =================\n\n");

  /* check to see if number of fluid zones has changed */
  checkCurrentZones(1);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

DEFINE_EXECUTE_AT_END(writeFreq)
{
#if !RP_HOST
  /* check if data should be dumped at this time step */
  if (N_TIME%ex.freq != 0)
    return;

  /* check to see if number of fluid zones has changed */
  checkCurrentZones(0);

  Message0("\n ================= Custom Data Export =================\n\n");
  Message0("  Exporting data to: %s\n", ex.dir);
  
  if (ex.format == 0)
    UDFexportEnsight(ex.dir, ex.prefix);
  else if (ex.format == 1)
    VTKexport(ex.dir, ex.prefix);

  Message0("\n =========================================================\n\n");
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void UDFinitExport(void)
{
  /* set standard options */
  ex.single_zone   = 0; /* dump entire volume as single part */
  ex.binary        = 1; /* write in binary format */
  ex.init          = 0; /* initialize the case */
  ex.freq          = 1; /* export frequency */
  ex.dump_velocity = 1;
  ex.dump_particle = 0; /* dump particle history */
  ex.dump_sdm      = 0; /* don't dump single droplet model variables */
  ex.dump_species  = 0; /* don't dump any species */
  ex.single_dir    = 1; /* write everything to single data directory */
  /* standard files */
  strcpy(ex.dir, ".");
  strcpy(ex.prefix, "fluent_export");

  int i, j, k, nYi, nSDM = 18, nBase = 7;
  char buf[1000];

  /* retrieve the ensight directory */
  if (EX_RP_is_Defined("export/directory"))
    strcpy(ex.dir, RP_Get_String("export/directory"));
  if (EX_RP_is_Defined("export/file-prefix"))
    strcpy(ex.prefix, RP_Get_String("export/file-prefix"));
  if (EX_RP_is_Defined("export/init"))
    ex.init = RP_Get_Boolean("export/init");
  if (EX_RP_is_Defined("export/freq"))
    ex.freq = RP_Get_Integer("export/freq");
  if (EX_RP_is_Defined("export/binary"))
    ex.binary = RP_Get_Boolean("export/binary");
  ex.nZone = 0;
  if (EX_RP_is_Defined("export/zones"))
    ex.nZone = parseFluentString_EX(RP_Get_String("export/zones"), &ex.zones);
  ex.nSurf = 0;
  if (EX_RP_is_Defined("export/wall-surfaces"))
    ex.nSurf = parseFluentString_EX(RP_Get_String("export/wall-surfaces"),
                                    &ex.surfs);
  ex.nCxSurf = 0;
  if (EX_RP_is_Defined("export/surfaces"))
    ex.nCxSurf = parseFluentString_EX(RP_Get_String("export/surfaces"),
                                      &ex.cxSurfs);

  ex.byteOrder = EX_checkEndian();

  int fluentZones = 0;
  if (ex.nZone == 0) { /* default to all fluid threads */
    Message0("  Exporting all Fluent zones\n");
    fluentZones = 1;
    char **zNames;
    ex.nZone = getFluidThreadNames(&zNames);
    Message0(" %i zones found\n", ex.nZone);
    ex.zones = (char**) calloc(ex.nZone, sizeof(char*));
    for (i=0; i<ex.nZone; i++) {
      ex.zones[i] = (char*)  calloc(128, sizeof(char));
      strcpy(ex.zones[i], zNames[i]);
      free(zNames[i]);
    }
    free(zNames);
  }

  int nZt;
  char **tNames;
  if (EX_RP_is_Defined("export/zone-threads")) {
    nZt = parseFluentString_EX(RP_Get_String("export/zone-threads"), &tNames);

    if (nZt == 0) { /* default to all fluid threads */
      nZt = getFluidThreadNames(&tNames);
    }

    /* allocate memory for zones */
    ex.nZoneThreads = (int*)    calloc(ex.nZone, sizeof(int));
    ex.zThreadNames = (char***) calloc(ex.nZone, sizeof(char**));
    
    if (fluentZones) { /* each zone has one thread */
      for (i=0; i<ex.nZone; i++) {
        ex.nZoneThreads[i]    = 1;
        ex.zThreadNames[i]    = (char**) calloc(1,   sizeof(char*));
        ex.zThreadNames[i][0] = (char*)  calloc(128, sizeof(char));
        strcpy(ex.zThreadNames[i][0], ex.zones[i]);
      }
    }
    else {
      EXparseThreadNames(nZt, tNames, ex.nZoneThreads, ex.zThreadNames);
    }
  }

  int fluentWalls = 0;
  if (ex.nSurf == 0) { /* default to all wall threads */
   Message0("  Exporting all Fluent zones\n");
    fluentWalls = 1;
    char **wNames;
    ex.nSurf = getWallThreadNames(&wNames);;
    ex.surfs = (char**) calloc(ex.nSurf, sizeof(char*));
    for (i=0; i<ex.nSurf; i++) {
      ex.surfs[i] = (char*)  calloc(128, sizeof(char));
      strcpy(ex.surfs[i], wNames[i]);
      free(wNames[i]);
    }
    free(wNames);
  }

  int nSt;
  char **sNames;
  char *defaultEngThreads;
  if (ex.nSurf > 0) {
    if (!strcmp(ex.surfs[0], "default-engine-symmetry")) {
      Message0("Using default symmetry engine surfaces\n");
      /* define the default surfaces for symmetry engine */
      ex.nSurf = 7;
      ex.surfs = (char**) realloc(ex.surfs, ex.nSurf*sizeof(char*));
      for (i=1; i<ex.nSurf; i++)
        ex.surfs[i] = (char*) calloc(234, sizeof(char));
      defaultEngThreads = (char*) calloc(2048, sizeof(char));

      strcpy(ex.surfs[0], "cylinder");
      strcpy(defaultEngThreads,
             "(cyl-head wall-fluid-ch cyl-tri cyl-layer cyl-piston) ");
      strcpy(ex.surfs[1], "piston");
      strcat(defaultEngThreads, "(piston) ");
      strcpy(ex.surfs[2], "intake-manifold");
      strcat(defaultEngThreads,
             "(invalve-1-port invalve1-seat ice-inlet-invalve-1-port) ");
      strcpy(ex.surfs[3], "exhaust-manifold");
      strcat(defaultEngThreads,
             "(exvalve-1-port exvalve1-seat ice-outlet-exvalve-1-port id-38) ");
      strcpy(ex.surfs[4], "intake-valve");
      strcat(defaultEngThreads,
             "(invalve1-ch invalve1-ib invalve1-ob invalve1-stem) ");
      strcpy(ex.surfs[5], "exhaust-valve");
      strcat(defaultEngThreads,
             "(exvalve1-ch exvalve1-ib exvalve1-ob exvalve1-stem) ");
      strcpy(ex.surfs[6], "symmetry");
      strcat(defaultEngThreads, "(symmetry) ");

      nSt = parseFluentString_EX(defaultEngThreads, &sNames);
    }
    else if (EX_RP_is_Defined("export/wall-surface-threads")) {
      nSt = parseFluentString_EX(RP_Get_String("export/wall-surface-threads"),
                                 &sNames);
    }
    
    if (nSt == 0) { /* default to all wall threads */
      /* TODO: I think this can only work with one zone */
     /*  nSt = getWallThreadNames(&sNames); */
    }
    
    /* allocate memory for zones */
    ex.nSurfThreads = (int*)    calloc(ex.nSurf, sizeof(int));
    ex.sThreadNames = (char***) calloc(ex.nSurf, sizeof(char**));

    if (fluentWalls) { /* each wall has one thread */
      for (i=0; i<ex.nSurf; i++) {
        ex.nSurfThreads[i]    = 1;
        ex.sThreadNames[i]    = (char**) calloc(1,   sizeof(char*));
        ex.sThreadNames[i][0] = (char*)  calloc(128, sizeof(char));
        strcpy(ex.sThreadNames[i][0], ex.surfs[i]);
      }
    }
    else {
      EXparseThreadNames(nSt, sNames, ex.nSurfThreads, ex.sThreadNames);
    }
    /* clear memory allocated in parseFluentString */
    /* for (i=0; i<nSt; i++) free(sNames); */
/*     free(sNames); */
  }

  /* interior defined surfaces through CX */
  if (ex.nCxSurf > 0) {
    nSt = 0;
    if (EX_RP_is_Defined("export/surface-threads")) {
      nSt = parseFluentString_EX(RP_Get_String("export/surface-threads"),
                                 &sNames);
    }

    /* allocate memory for zones */
    ex.nCxSurfThreads = (int*)    calloc(ex.nCxSurf, sizeof(int));
    ex.cxThreadNames = (char***) calloc(ex.nCxSurf, sizeof(char**));
    
    if (nSt > 0) {
      EXparseThreadNames(nSt, sNames, ex.nCxSurfThreads, ex.cxThreadNames);
    }
    else { /* default one thread per surface */
      for (i=0; i<ex.nCxSurf; i++) {
        ex.nCxSurfThreads[i]    = 1;
        ex.cxThreadNames[i]    = (char**) calloc(1,   sizeof(char*));
        ex.cxThreadNames[i][0] = (char*)  calloc(128, sizeof(char));
        strcpy(ex.cxThreadNames[i][0], ex.cxSurfs[i]);
      }
    }
    /* clear memory allocated in parseFluentString */
    /* for (i=0; i<nSt; i++) free(sNames); */
/*     free(sNames); */

    /* need to get the reference id values for the surfaces from RP, as
       I can't figure out how to get them from the code */
    if (EX_RP_is_Defined("export/surface-ref-names")) {
      ex.nCxRef =
        parseFluentString_EX(RP_Get_String("export/surface-ref-names"),
                             &ex.cxRefNames);
    }
    else {
      Message0("Warning: couldn't find reference surface names\n");
    }
    if (EX_RP_is_Defined("export/surface-ref-ids")) {
      nSt = parseFluentString_EX(RP_Get_String("export/surface-ref-ids"),
                                 &sNames);
      if (nSt != ex.nCxRef)
        Message0("Warning: no. of reference surf ids and names don't match\n");
      ex.cxRefIds = (int*) calloc(ex.nCxRef, sizeof(int));
      for (i=0; i<ex.nCxRef; i++) {
        ex.cxRefIds[i] = atoi(sNames[i]);
        /* free(sNames[i]); */
      }
      /* free(sNames); */
    }
    else {
      Message0("Warning: couldn't find reference surface ids\n");
    }
    

  }

  ex.format = 0;
  if (!strcmp(strrchr(ex.dir, '/')+1, "paraview") ||
      !strcmp(strrchr(ex.dir, '/')+1, "vtk") ) {
    ex.inLine = 0;
    ex.format = 1;
  }

  if (ex.binary) {
    writeDataStr  = &writeBinStr;
    writeDataInt  = &writeBinInt;
    writeDataReal = &writeBinReal;
    strcpy(ex.vtkFormat, "raw");
  }
  else {
    writeDataStr  = &writeAsciiStr;
    writeDataInt  = &writeAsciiInt;
    writeDataReal = &writeAsciiReal;
    strcpy(ex.vtkFormat, "ascii");
  }

  /* assign the data write method */
  writeBinaryData = &fwrite;
/*   writeBinaryData = &BINARY_FILE_WRITE; */
  
  /* standard variables */
  ex.nSc = nBase;

  /* check for species export */
  ex.nY = 0;
  if (EX_RP_is_Defined("export/species"))
    ex.nY = parseFluentString_EX(RP_Get_String("export/species"), &ex.Y_names);
  if (ex.nY > 0) {
    ex.dump_species = 1;
    ex.Y_indx = (int*) calloc(ex.nY, sizeof(int));
    
    /* get the species indices */
#if !RP_HOST
    Domain *d;
    Material *mat, *sp;
    d   = Get_Domain(1);        /* get the domain from Fluent */
    mat = mixture_material(d);  /* get domain material */
    for (k=0; k<ex.nY; k++) {
      mixture_species_loop(mat,sp,i) {
        strcpy(buf, MIXTURE_SPECIE_NAME(mat,i));
        /* convert species to upper case */
        for (j=0; j<25; j++) buf[j] = toupper(buf[j]);
        if (!strcmp(ex.Y_names[k], buf))
          ex.Y_indx[k] = i;
      }
    }
#endif /* !RP_HOST */

    /* output the gas constant and specific heat */
    ex.nSc += 2;
  }

  /* find out whether to dump SDM variables */
  if (EX_RP_is_Defined("export/dump-sdm"))
    ex.dump_sdm = RP_Get_Boolean("export/dump-sdm");

  /* check whether particles are dumped */
  if (EX_RP_is_Defined("export/dump-particle"))
    ex.dump_particle = RP_Get_Boolean("export/dump-particle");

  /* update the total number of scalars */
  ex.nSc += N_UDS;
  ex.nSc += N_UDM;
  if (ex.dump_sdm)     ex.nSc += nSDM;
  if (ex.dump_species) ex.nSc += ex.nY;

  /* allocate memory */
  ex.sc_indx  = (int*)   calloc(ex.nSc, sizeof(int));
  ex.sc_type  = (int*)   calloc(ex.nSc, sizeof(int));
  ex.sc_names = (char**) calloc(ex.nSc, sizeof(char*));
  ex.sc_abbrv = (char**) calloc(ex.nSc, sizeof(char*));
  for (i=0; i<ex.nSc; i++) {
    ex.sc_names[i] = (char*) calloc(80, sizeof(char));
    ex.sc_abbrv[i] = (char*) calloc(25, sizeof(char));
  }

  /* ----- set the standard variables ----- */
  k = 0;
  ex.sc_indx[0] = SV_P;                     ex.sc_type[0] = 0;
  strcpy(ex.sc_names[0], "Pressure");       strcpy(ex.sc_abbrv[0], "P");
  ex.sc_indx[1] = SV_T;                     ex.sc_type[1] = 1;
  strcpy(ex.sc_names[1], "Temperature");    strcpy(ex.sc_abbrv[1], "T");
  ex.sc_indx[2] = SV_DENSITY;               ex.sc_type[2] = 1;
  strcpy(ex.sc_names[2], "Density");        strcpy(ex.sc_abbrv[2], "RHO");
  ex.sc_indx[3] = SV_K;                     ex.sc_type[3] = 1;
  strcpy(ex.sc_names[3], "T_Kin_Energy");   strcpy(ex.sc_abbrv[3], "TKE");
  ex.sc_indx[4] = SV_D;                     ex.sc_type[4] = 1;
  strcpy(ex.sc_names[4], "T_Dissipation");  strcpy(ex.sc_abbrv[4], "EPS");
  ex.sc_indx[5] = SV_MU_T;                  ex.sc_type[5] = 1;
  strcpy(ex.sc_names[5], "T_Viscosity");    strcpy(ex.sc_abbrv[5], "MU_T");
  ex.sc_indx[6] = SV_VOLUME;                ex.sc_type[6] = 1;
  strcpy(ex.sc_names[6], "Cell_Volume");    strcpy(ex.sc_abbrv[6], "VOL");
  k = nBase;

  /* output the gas constant if species present */
  if (ex.nY > 0) {
    ex.sc_indx[k] = SV_CP;                    ex.sc_type[k] = 1;
    strcpy(ex.sc_names[k], "Specific_Heat");  strcpy(ex.sc_abbrv[k], "CP");
    k++;
    ex.sc_indx[k] = SV_RGAS;                ex.sc_type[k] = 1;
    strcpy(ex.sc_names[k], "Gas_Constant"); strcpy(ex.sc_abbrv[k],"RGAS");
    k++;
  }

  /* ----- initialize user-defined-memory export ----- */
  for (i=0; i<N_UDS; i++) {
    ex.sc_indx[k+i] = SV_UDS_I(i);
    ex.sc_type[k+i] = 1;
    if (strncmp(user_scalar_vars[i].name, "User", 4)) {
      strcpy(ex.sc_names[k+i], user_scalar_vars[i].name);
      strcpy(ex.sc_abbrv[k+i], user_scalar_vars[i].name);
    }
    else {
      sprintf(buf, "UDS_%i", i);
      strcpy(ex.sc_names[k+i], buf);
      strcpy(ex.sc_abbrv[k+i], buf);
    }
  }
  k += N_UDS;

  /* ----- initialize user-defined-memory export ----- */
  for (i=0; i<N_UDM; i++) {
    ex.sc_indx[k+i] = user_memory_vars[i].index;
    ex.sc_type[k+i] = 2;
    if (strncmp(user_memory_vars[i].name, "User", 4)) {
      strcpy(ex.sc_names[k+i], user_memory_vars[i].name);
      strcpy(ex.sc_abbrv[k+i], user_memory_vars[i].name);
    }
    else {
      sprintf(buf, "UDM_%i", i);
      strcpy(ex.sc_names[k+i], buf);
      strcpy(ex.sc_abbrv[k+i], buf);
    }
  }
  k += N_UDM;

  /* ----- initialise species export ----- */
  if (ex.dump_species) {
    for (i=0; i<ex.nY; i++) {
      ex.sc_indx[k+i] = ex.Y_indx[i];
      ex.sc_type[k+i] = 3;
      strcpy(ex.sc_names[k+i], "Y_");
      strcat(ex.sc_names[k+i], ex.Y_names[i]);
      strcpy(ex.sc_abbrv[k+i], ex.Y_names[i]);
    }
  }
  k += ex.nY;

  /* retrieve the information about the particle streams */
  ex.streams = NULL;
  if (ex.dump_particle) {
    if (EX_RP_is_Defined("spray/nozzles"))
      ex.nNoz = RP_Get_Integer("spray/nozzles");
    if (EX_RP_is_Defined("spray/injection-ids"))
      ex.nInj = parseFluentString_EX(RP_Get_String("spray/injection-ids"),
                                     &ex.inject_id);

    ex.streams     = (Injection***) calloc(ex.nInj, sizeof(Injection**));
    ex.all_streams = (Injection**)  calloc(ex.nInj*ex.nNoz, sizeof(Injection*));
    
    for (i=0; i<ex.nInj; i++) {
      EXstripStringBraces(ex.inject_id[i]);
      getInjection(ex.inject_id[i], ex.nNoz, &ex.streams[i]);
      for (j=0; j<ex.nNoz; j++)
        ex.all_streams[i*ex.nNoz+j] = ex.streams[i][j];
    }
  } /* end particle dump */
}
/* -------------------------------------------------------------------------- */

void checkCurrentZones(int force)
{
  /* check to see that the number of zones hasn't changed */
  int i, j, currentCells;
  currentCells = countFluidThreads();
  
  /* reinitialize if cell count has changed */
  if ((currentCells != ex.nCells) || force == 1) {
    Message0("  Reinitializing fluid zones:\n");
    destroyExportStruct(&ex);
    ex.nCells = initThreads(&ex.cell_th, &ex.part_id, &ex.part_name,
                            &ex.part_type, &ex.part_nv, &ex.part_verts);
    
    Message0("  New Fluid Zones:        %i\n%22s", ex.nCells, " ");
    for (i=0; i<ex.nCells; i++) {
      Message0("%s ( ", ex.part_name[i]);
      for (j=0; j<ex.part_nv[i]; j++)
        Message0("%s ", ex.part_type[i][j]);
      Message0(")\n%22s", " ");
    }
    Message0("\n");
  }
}
/* -------------------------------------------------------------------------- */

int countFluidThreads()
{
  Domain *d = Get_Domain(1);
  Thread *t;
  int iT, nT = 0;
  thread_loop(t,d,iT) {
    if (FLUID_CELL_THREAD_P(t)) {
      nT++;
    }
  } /* end loop over all threads */
  return nT;
}
/* -------------------------------------------------------------------------- */

void destroyExportStruct(expStruct *e)
{
  int i, j;
  for (i=0; i<e->nCells; i++) {
    for (j=0; j<e->part_nv[i]; j++) {
      if (e->part_type[i][j]) free(e->part_type[i][j]);
    }
    if (e->part_type[i])  free(e->part_type[i]);
    if (e->part_verts[i]) free(e->part_verts[i]);
    if (e->part_name[i])  free(e->part_name[i]);
  }
  if (e->part_type)  free(e->part_type);
  if (e->part_verts) free(e->part_verts);
  if (e->part_name)  free(e->part_name);
  if (e->part_nv)    free(e->part_nv);
  if (e->part_id)    free(e->part_id);
  if (e->cell_th)    free(e->cell_th);
}
/* -------------------------------------------------------------------------- */

int initThreads(Thread ***tF, int **ids, char ***name, char ****type, int **nv,
                int ***verts)
{
  Domain *d = Get_Domain(1);
  Thread *t;
  int i, iP, iT, nT;
  char tmp[200];

  tF[0]    = (Thread**) calloc(1, sizeof(Thread*));
  ids[0]   = (int*)     calloc(1, sizeof(int));
  nv[0]    = (int*)     calloc(1, sizeof(int));
  verts[0] = (int**)    calloc(1, sizeof(int*));
  name[0]  = (char**)   calloc(1, sizeof(char*));
  type[0]  = (char***)  calloc(1, sizeof(char**));

  /* loop over all threads */
  nT = 0;  
  i = 0;
  thread_loop(t,d,iT) {
    if (FLUID_CELL_THREAD_P(t)) {
      nT++;
      iP = nT - 1;
      
      /* store the thread pointer */
      tF[0]     = (Thread**) realloc(tF[0], nT*sizeof(Thread*));
      tF[0][iP] = t;
      
      /* store the thread id */
      ids[0]     = (int*) realloc(ids[0], nT*sizeof(int));
      ids[0][iP] = THREAD_ID(t);
      
      /* store the thread name */
      name[0]     = (char**) realloc(name[0], nT*sizeof(char**));
      name[0][iP] = (char*)  calloc(256, sizeof(char));
      strcpy(name[0][iP], THREAD_NAME(t));
      if (!strcmp(THREAD_NAME(t), "")) {
        sprintf(tmp, "fluid-%i", THREAD_ID(t));
        strcpy(name[0][iP], tmp);
      }
      
      /* keep track of the number of vertices */
      nv[0]    = (int*)    realloc(nv[0],    nT*sizeof(int));
      verts[0] = (int**)   realloc(verts[0], nT*sizeof(int*));
      type[0]  = (char***) realloc(type[0],  nT*sizeof(char**));
      if (THREAD_ELEMENT_TYPE(t) != MIXED_CELL) {
        nv[0][iP] = 1;
        verts[0][iP]   = (int*)   calloc(1,  sizeof(int));
        type[0][iP]    = (char**) calloc(1,  sizeof(char*));
        type[0][iP][0] = (char*)  calloc(30, sizeof(char));
      }
      
      /* figure out the cell type */
      switch (THREAD_ELEMENT_TYPE(t)) {
      case MIXED_CELL:
        nv[0][iP] = getNumberCellTypes(t, &verts[0][iP], &type[0][iP]);
        break;
        /* 2D cell types */
      case TRI_CELL:
        strcpy(type[0][iP][0], "tria3");
        break;
      case QUAD_CELL:
        strcpy(type[0][iP][0], "quad4");
        break;
        /* 3D cell types */
      case TET_CELL:
        strcpy(type[0][iP][0], "tetra4");
        verts[0][iP][0] = 4;
        break;
      case PYRAMID_CELL:
        strcpy(type[0][iP][0], "pyramid5");
        verts[0][iP][0] = 5;
        break;
      case WEDGE_CELL:
        strcpy(type[0][iP][0], "penta6");
        verts[0][iP][0] = 6;
        break;
      case HEX_CELL:
        strcpy(type[0][iP][0], "hexa8");
        verts[0][iP][0] = 8;
        break;
      case POLY_CELL:
        strcpy(type[0][iP][0], "nfaced");
        ex.poly_mesh = 1;
        break;
      } /* end switch */
    }
  } /* end loop over all threads */

  return nT;
}
/* -------------------------------------------------------------------------- */

int getFluidThreadNames(char ***names)
{
  Domain *d = Get_Domain(1);
  Thread *t;
  int nT, iT, iP;
  char tmp[256];

  names[0] = (char**) calloc(1, sizeof(char*));

  nT = 0;  
  thread_loop(t,d,iT) {
    if (FLUID_CELL_THREAD_P(t)) {
      nT++;
      iP = nT - 1;

      /* store the thread name */
      names[0]     = (char**) realloc(names[0], nT*sizeof(char**));
      names[0][iP] = (char*)  calloc(256, sizeof(char));
      strcpy(names[0][iP], THREAD_NAME(t));
      if (!strcmp(THREAD_NAME(t), "")) {
        sprintf(tmp, "fluid-%i", THREAD_ID(t));
        strcpy(names[0][iP], tmp);
      }
    }
  } /* end loop over all threads */

  return nT;
}
/* -------------------------------------------------------------------------- */

int getWallThreadNames(char ***names)
{
  Domain *d = Get_Domain(1);
  Thread *t;
  int nT, iT, iP;
  char tmp[256];

  names[0] = (char**) calloc(1, sizeof(char*));

  nT = 0;  
  thread_loop(t,d,iT) {
    if (THREAD_TYPE(t) ==  THREAD_F_WALL) {
      nT++;
      iP = nT - 1;

      /* store the thread name */
      names[0]     = (char**) realloc(names[0], nT*sizeof(char**));
      names[0][iP] = (char*)  calloc(256, sizeof(char));
      strcpy(names[0][iP], THREAD_NAME(t));
      if (!strcmp(THREAD_NAME(t), "")) {
        sprintf(tmp, "wall-%i", THREAD_ID(t));
        strcpy(names[0][iP], tmp);
      }
    }
  } /* end loop over all threads */

  return nT;
}
/* -------------------------------------------------------------------------- */

int getNumberCellTypes(Thread *t, int **verts, char ***types)
{
  int totTypes, uniqueTypes;
#if !RP_HOST
  cell_t c;
  int i, j, *nTypes, *vertsL, *all_verts, *indx, *tOff, *iWork, newType;
  int *cTypes, *all_types, *typeN;

  vertsL = (int*) calloc(1, sizeof(int));
  cTypes = (int*) calloc(1, sizeof(int));
  nTypes = (int*) calloc(compute_node_count, sizeof(int));
  tOff   = (int*) calloc(compute_node_count, sizeof(int));
  iWork  = (int*) calloc(compute_node_count, sizeof(int));
  
#if !PARALLEL
  int myid = 0;
#endif /* !PARALLEL */

  /* loop over all cells to count the number of types on this process */
  nTypes[myid] = 0;
  begin_c_loop_int(c,t) {
    newType = 1;
    for (i=0; i<nTypes[myid]; i++)
      if (C_TYPE(c,t) == cTypes[i])
        newType = 0;

    if (C_TYPE(c,t) == 0)
      Message("Warning: Cell type for %i is 0: verts: %i\n", c, C_NNODES(c,t));

    if (C_TYPE(c,t) == 7)
      ex.poly_mesh = 1;

    if (newType) {
      nTypes[myid]++;
      vertsL = (int*) realloc(vertsL, nTypes[myid]*sizeof(int));
      cTypes = (int*) realloc(cTypes, nTypes[myid]*sizeof(int));
      vertsL[nTypes[myid]-1] = C_NNODES(c,t);
      cTypes[nTypes[myid]-1] = C_TYPE(c,t);
    }
  } end_c_loop_int(c,t);

#if PARALLEL
  /* find out if any processor had a polyhedron */
  ex.poly_mesh = PRF_GIHIGH1(ex.poly_mesh);

  /* accumulate number of types onto all nodes */
  PRF_GISUM(nTypes, compute_node_count, iWork);

  /* get the number of types found over all processors */
  for (i=0, totTypes=0; i<compute_node_count; i++)
    totTypes += nTypes[i];

  /* find the offset for each node when transferring to the root */
  for (i=1, tOff[0]=0; i<compute_node_count; i++)
    for (j=0, tOff[i]=0; j<i; j++)
      tOff[i] += nTypes[j];

  /* allocate the memory on all processors */
  all_verts = (int*) calloc(totTypes, sizeof(int));
  all_types = (int*) calloc(totTypes, sizeof(int));
  indx      = (int*) calloc(totTypes, sizeof(int));

  if (myid == node_zero) {
    /* fill types found on root into global array */
    for (i=0; i<nTypes[myid]; i++) {
      all_verts[i] = vertsL[i];
      all_types[i] = cTypes[i];
    }
    
    /* fill types found on all other processors */
    for (i=1; i<compute_node_count; i++) {
      PRF_CRECV_INT(i, &all_verts[tOff[i]], nTypes[i], myid);
      PRF_CRECV_INT(i, &all_types[tOff[i]], nTypes[i], myid);
    }
    
    /* sort types found */
    indexIntArray(totTypes, all_types, indx);
    
    /* count the number of unique types found across all nodes */
    uniqueTypes = 1;
    for (i=1; i<totTypes; i++) {
      if (all_types[indx[i]] != all_types[indx[i-1]]) {
        uniqueTypes++;
      }
    }
    
    /* send the number of types to all processors */
    for (i=1; i<compute_node_count; i++) 
      PRF_CSEND_INT(i, &uniqueTypes, 1, myid);
  }
  else {
    PRF_CSEND_INT(node_zero, vertsL, nTypes[myid], node_zero);
    PRF_CSEND_INT(node_zero, cTypes, nTypes[myid], node_zero);
    PRF_CRECV_INT(node_zero, &uniqueTypes, 1, node_zero);
  }
  
  /* allocate the memory on all processes */
  verts[0] = (int*)   calloc(uniqueTypes, sizeof(int));
  typeN    = (int*)   calloc(uniqueTypes, sizeof(int));
  types[0] = (char**) calloc(uniqueTypes, sizeof(char*));

  if (myid == node_zero) { 
    /* assign only unique types to the final list of vertices */
    verts[0][0] = all_verts[indx[0]];
    typeN[0]    = all_types[indx[0]];
    for (i=1, j=1; i<totTypes; i++) {
      if (all_types[indx[i]] != all_types[indx[i-1]]) {
        verts[0][j] = all_verts[indx[i]];
        typeN[j]    = all_types[indx[i]];
        j++;
      }
    }
    
    /* send out the updated list to all nodes */
    for (i=1; i<compute_node_count; i++) {
      PRF_CSEND_INT(i, verts[0], uniqueTypes, myid);
      PRF_CSEND_INT(i, typeN,    uniqueTypes, myid);
    }
  }
  else {
    PRF_CRECV_INT(node_zero, verts[0], uniqueTypes, node_zero);
    PRF_CRECV_INT(node_zero, typeN,    uniqueTypes, node_zero);
  }
  
  free(all_verts); free(all_types); free(indx);
#else /* SERIAL */
  uniqueTypes = nTypes[myid];
  verts[0] = (int*)   calloc(uniqueTypes, sizeof(int));
  types[0] = (char**) calloc(uniqueTypes, sizeof(char*));
  for (i=0; i<uniqueTypes; i++)
    verts[0][i] = vertsL[i];
#endif /* PARALLEL */

  /* assign the type names based on number of vertices found */
  for (i=0; i<uniqueTypes; i++) {
    types[0][i] = (char*) calloc(30, sizeof(char));
    switch (typeN[i]) {
    case TRI_CELL:
      strcpy(types[0][i], "tria3");
      break;
    case QUAD_CELL:
      strcpy(types[0][i], "quad4");
      break;
    case TET_CELL:
      strcpy(types[0][i], "tetra4");
      break;
    case PYRAMID_CELL:
      strcpy(types[0][i], "pyramid5");
      break;
    case WEDGE_CELL:
      strcpy(types[0][i], "penta6");
      break;
    case HEX_CELL:
      strcpy(types[0][i], "hexa8");
      break;
    case POLY_CELL:
      strcpy(types[0][i], "polyhedron");
      ex.poly_mesh = 1;
      break;
    }
  }
  
  free(vertsL); free(typeN); free(cTypes); free(nTypes); free(tOff);
  free(iWork);
#endif /* !RP_HOST */

  return uniqueTypes;
}
/* -------------------------------------------------------------------------- */

int initSpecThreadNames(int **ids, char ***name, char ***type)
{
  /* This function requires that a string variable with the zones has been 
     defined in the case file as: thread-names 
     See journal/zones.scm for the scheme function */
  int i, i_pre, ni, nt, nn, nThreads = 0;
  char *t_ids, *t_names, *t_types, s[1000];
  char *c_ids, *c_names, *c_types;  /* for cells */
  char *w_ids, *w_names, *w_types;  /* for walls */
  t_ids   = (char*) calloc(2000, sizeof(char));
  t_names = (char*) calloc(2000, sizeof(char));
  t_types = (char*) calloc(2000, sizeof(char));
  c_ids   = RP_Get_String("cell-thread-ids");
  c_names = RP_Get_String("cell-thread-names");
  c_types = RP_Get_String("cell-thread-types");
  w_ids   = RP_Get_String("wall-thread-ids");
  w_names = RP_Get_String("wall-thread-names");
  w_types = RP_Get_String("wall-thread-types");
  strcpy(t_ids, c_ids); strcat(t_ids, w_ids);
  strcpy(t_names, c_names); strcat(t_names, w_names);
  strcpy(t_types, c_types); strcat(t_types, w_types);

  ids[0]  = (int*)   calloc(1, sizeof(int));
  name[0] = (char**) calloc(1, sizeof(char*));
  type[0] = (char**) calloc(1, sizeof(char*));

  while(t_ids[0]) {
    /* read in the first zone id and increment */
    sscanf(t_ids,"%s%n",s,&ni); t_ids += ni;
    i = atoi(s);
    if (i == i_pre)
      break;
    if (i > nThreads) {
      nThreads++;
      name[0]    = realloc(name[0], (i+1)*sizeof(char*));
      type[0]    = realloc(type[0], (i+1)*sizeof(char*));
      name[0][i] = calloc(100, sizeof(char));
      type[0][i] = calloc(30, sizeof(char));
    }
    else {
      name[0][i] = calloc(100, sizeof(char));
      type[0][i] = calloc(30, sizeof(char));
    }
    ids[0] = realloc(ids[0],nThreads*sizeof(int));
    ids[0][nThreads-1] = i;
    sscanf(t_names,"%s%n",name[0][i],&nn); t_names += nn;
    sscanf(t_types,"%s%n",type[0][i],&nt); t_types += nt;
    Message("%i %s %s\n", i, name[0][i], type[0][i]);
    i_pre = i;
  }
  /*  free(t_ids); free(t_names); free(t_types); */
  return nThreads;
}
/* -------------------------------------------------------------------------- */

#if PARALLEL && UDF_MPI_WRITE
void writeBinStr(MPI_File *pFile, int len, char *str)
{
#if !RP_HOST
  char *buf = (char*) calloc(len, sizeof(char));
  sprintf(buf, "%s", str);
  MPI_Status status;
  MPI_File_write(pFile, buf, len, MPI_CHAR, &status);
  free(buf);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */
void writeBinInt(MPI_File *pFile, int len, int *var)
{
#if !RP_HOST
  MPI_Status status;
  MPI_File_write(pFile, var, len, MPI_INT, &status);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */
void writeBinReal(MPI_File *pFile, int len, float *var)
{
#if !RP_HOST
  MPI_Status status;
  MPI_File_write(pFile, var, len, MPI_FLOAT, &status);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */
#else /* SERIAL || PARALLEL ROOT WRITE */
void writeBinStr(FILE *pFile, int len, char *str)
{
  char *buf = (char*) calloc(len, sizeof(char));
  sprintf(buf, "%s", str);
/*   writeBinaryData(buf, sizeof(char), len, pFile); */
  fwrite(buf, sizeof(char), len, pFile);
  free(buf);
}
void writeBinInt(FILE *pFile, int len, int *var, int perLine)
{
  fwrite(var, sizeof(int), len, pFile);
/*   writeBinaryData(var, sizeof(int), len, pFile); */
}
void writeBinReal(FILE *pFile, int len, float *var, int perLine)
{
  fwrite(var, sizeof(float), len, pFile);
/*   writeBinaryData(var, sizeof(float), len, pFile); */
}

void writeAsciiStr(FILE *pFile, int len, char *str)
{
  char *buf = (char*) calloc(len, sizeof(char));
  sprintf(buf, "%s", str);
  fprintf(pFile, "%s\n", buf);
  free(buf);
}
void writeAsciiInt(FILE *pFile, int len, int *var, int perLine)
{
  int i, line = 1;
  for (i=0; i<len; i++)
    {
      fprintf(pFile, "%10d", var[i]);
      if (line%perLine == 0) fprintf(pFile, "\n");
      line++;
    }
}
void writeAsciiReal(FILE *pFile, int len, float *var, int perLine)
{
  int i, line = 1;
  for (i=0; i<len; i++)
    {
      fprintf(pFile, "%12.5e", var[i]);
      if (line%perLine == 0) fprintf(pFile, "\n");
      line++;
    }
}
/* -------------------------------------------------------------------------- */
#endif /* PARALLEL && UDF_MPI_WRITE */

int BINARY_FILE_WRITE(void *buffer, int size, int count, FILE *pFile)
{
  int i,n;
  unsigned char *tmp;
  char *pos;

  tmp = malloc(size);
  pos = (char*) buffer;
   
  for (n=0; n<count; n++) {
    memcpy(tmp,pos,size);
    for (i=0; i<size; i++) {
      fputc(tmp[i], pFile);
    }
    pos += size;
  }
   
  free(tmp);
  return 0;
} 
/* ------------------------------------------------------------------------- */

int parseFluentString_EX(char *input, char ***var)
{
  int i, n, cur, len;
  char *name, *sWork, s[1024];

  if (input == NULL)
    return 0;

  /* copy the input string */
  sWork = (char*) calloc(1024, sizeof(char));
  strcpy(sWork, input);

  var[0] = (char**) calloc(1, sizeof(char*));
  i = cur = 0;
  len = strlen(sWork);
  while(cur < len) {
    /* read in the first zone id and increment */
    sscanf(&sWork[cur], "%s", s);
    cur += (strlen(s)+1);
    
    /* allocate additional memory */
    var[0]    = realloc(var[0], (i+1)*sizeof(char*));
    var[0][i] = calloc(1024, sizeof(char));
    
    /* copy variable */
    strcpy(var[0][i], s);
    i++;
  }
  free(sWork);
  return i;
}
/* ------------------------------------------------------------------------- */

void getInjection(char *name, int nNozzle, Injection ***I)
{
  Message0("Getting injection: %s\n", name);
  int i;
  char ident[80];
  I[0] = (Injection**) calloc(nNozzle, sizeof(Injection*));

  for (i=0; i<nNozzle; i++) {
    sprintf(ident, "%s-%i", name, i+1);
    I[0][i] = Pick_Injection(ident);
  }
}
/* ------------------------------------------------------------------------- */

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

int gridLocateIntIndx(int n1, int *x1, int *indx, int x)
{
  int jL = 0, jR = n1-1, jM;
 /* whether grid is ascending or not */
  int ascnd = (x1[indx[n1-1]] >= x1[indx[0]]); 

  /* use bisection to find nearest grid point */
  while (jR - jL > 1) {
    /* compute a midpoint */
    jM = (jL+jR)/2;
    if (x >= x1[indx[jM]] == ascnd)
      jL = jM;
    else
      jR = jM;
  }
  
  return IMAX(0, IMIN(n1-2, jL));
}
/* -------------------------------------------------------------------------- */

void EXstripString(char *str)
{
  int pos, len;
  char loc[100];

  /* remove any leading whitespace */
  pos = 0;
  while (loc[pos] == ' ' && pos < strlen(str)) ++pos;
  strcpy(loc, &str[pos]);

  /* remove any trailing whitespace */
  pos = strlen(loc) - 2;
  if (pos < 0)
    loc[0] = '\0'; /* empty string */
  else {
    while (loc[pos] == ' ' && pos > 0) --pos;
    loc[pos+1] = '\0';
  }
  
  /* assign original string to modified local string */
  strcpy(str,loc);
}
/* -------------------------------------------------------------------------- */

void EXstripStringBraces(char *name)
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

int EXparseThreadNames(int nRaw, char **raw, int *nT, char ***names)
{
  int i, iZ, nZ, iT, len, tot;
  char cStart, cEnd;

  iZ = nZ = 0;
  for (i=0; i<nRaw; i++) {
    nT[iZ] = nT[iZ] + 1;
    iT     = nT[iZ] - 1;

    names[iZ]     = (char**) realloc(names[iZ], nT[iZ]*sizeof(char*));
    names[iZ][iT] = (char*)   calloc(512, sizeof(char));

    /* get the current string characteristics */
    len    = strlen(raw[i]);
    cStart = raw[i][0];
    cEnd   = raw[i][len-1];

    if (cStart == '(') {
      if (cEnd == ')') { /* only one thread in zone */
        strncpy(names[iZ][iT], &raw[i][1], len-2);
        names[iZ][iT][len-2] = '\0';
        iZ++;
        nZ++;
      }
      else {
        strcpy(names[iZ][iT], &raw[i][1]);
      }
    }
    else if (cEnd == ')') {
      strncpy(names[iZ][iT], raw[i], len-1);
      iZ++;
      nZ++;
    }
    else {
      strcpy(names[iZ][iT], raw[i]);
    }
  } /* end loop over threads */

  return nZ;
}
/* -------------------------------------------------------------------------- */

int EX_getCxId(char *name)
{
  int i;
  for (i=0; i<ex.nCxRef; i++) {
    if (!strcmp(name, ex.cxRefNames[i]))
      return ex.cxRefIds[i];
  }
  Message0("Warning: couldn't find surface %s\n",name);
  return -1;
}
/* -------------------------------------------------------------------------- */

int EX_RP_is_Defined(char *var_name)
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

void EX_facet_area_3D( Surface *s, int i, double area[] )
{
  double v[MAX_FACE_NODES][3];
  int k, m;
  area[0]=area[1]=area[2]=0.0;
  
  for (k=0; k < s->facets[i]; k++) {
    v[k][0]=Get_Surface_Point_Coord(&(s->points[s->facets[i+k+1]]),X_DIM); 
    v[k][1]=Get_Surface_Point_Coord(&(s->points[s->facets[i+k+1]]),Y_DIM); 
    v[k][2]=Get_Surface_Point_Coord(&(s->points[s->facets[i+k+1]]),Z_DIM);
  }

  for (k=0; k < s->facets[i]; k++) {
    m = (k+1)%(s->facets[i]);
    area[0]+= 0.5*((v[m][1] + v[k][1]) * (v[m][2] - v[k][2]));
    area[1]+= 0.5*((v[m][2] + v[k][2]) * (v[m][0] - v[k][0]));
    area[2]+= 0.5*((v[m][0] + v[k][0]) * (v[m][1] - v[k][1]));
  }
}
/* ------------------------------------------------------------------------- */

int EX_checkEndian()
{
  int n = 1;
  if(*(char *)&n == 1) return 0; /* little endian */
  return 1; /* big endian */
}

int IMIN(int x, int y)
{
  if (x < y) return x;
  else return y;
}
/* ------------------------------------------------------------------------- */

int IMAX(int x, int y)
{
  if (x > y) return x;
  else return y;
}
/* ------------------------------------------------------------------------- */

double ABS_D(double x)
{
  if (x < 0.0) x = -x;
  return x;
}
/* ------------------------------------------------------------------------- */
