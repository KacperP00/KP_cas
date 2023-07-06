#include "vtkIO.h"

extern expStruct ex;

void VTKexport()
{
  int i;
#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL && RP_NODE */
    VTKwriteTimeInfo(ex.dir, ex.prefix);

  /* write the multiblock data set stuff */
  VTKwriteMultiBlock(ex.nZone, ex.zones);

  /* write the geometry for this time-step */
  for (i=0; i<ex.nZone; i++)
    VTKwriteUnstructData(ex.zones[i],ex.nZoneThreads[i],ex.zThreadNames[i]);

  if (ex.dump_particle) {
    for (i=0; i<ex.nInj; i++)
      VTKwriteParticleData(ex.inject_id[i], ex.nNoz, ex.streams[i]);
  }

  /* write any walls defined */
  ex.nAllSurfThreads = 0;
  ex.sThreadIDs = (int*) calloc(1, sizeof(int));
  for (i=0; i<ex.nSurf; i++)
    VTKwriteWallSurf(ex.surfs[i], ex.nSurfThreads[i], ex.sThreadNames[i]);
  free(ex.sThreadIDs);

  /* write any surfaces defined */
  ex.nAllCxThreads = 0;
  ex.cxThreadIDs = (int*) calloc(1, sizeof(int));
  for (i=0; i<ex.nCxSurf; i++)
    VTKwriteCxSurf(ex.cxSurfs[i],ex.nCxSurfThreads[i],ex.cxThreadNames[i]);
  free(ex.cxThreadIDs);
}
/* -------------------------------------------------------------------------- */

void VTKwriteTimeInfo(char *dir, char *prefix)
{
#if !RP_HOST

  int  i, j, nt, iStart, increment, iL, line, perLine = 5, iCstart, iCend;
  int lineL = 1024;
  long double r;
  double *tVals, cur_time;
  float tLine;
  char fileName[1024], buffer[1024], tag[64], *ident, **dataLines, *pch;
  char *cQ0, *cQ1, *str_time;
  int nXML, rXML, maxXML;
  char *xmlLine, *xmlType, **xmlRaw, **xmlTags, **xmlVals;
  sprintf(fileName, "%s/%s.pvd", dir, prefix);
  Message0("  %-46s  %s\n", "Writing vtk case file:", strrchr(fileName, '/')+1);

  FILE *pFile = NULL;
  pFile = fopen(fileName, "r");

  char *buff = (char*) calloc(lineL, sizeof(char));

  if (!pFile)
    Message0("Error: could not open file %s\n", fileName);
   
  if (in_cylinder_p)
    cur_time = RP_Get_Real("dynamesh/in-cyn/crank-angle");
  else
    cur_time = CURRENT_TIME;

  if (pFile && ex.init == 1) { /* file exists, read in contents to get time */
    /* read in until the time section is reached */
    iCstart = 0;
    while (!feof(pFile)) { 
      fgets(buff, 1024, pFile);  /* read in line */
      iCstart++;
      pch = strstr(buff, "<Collection>");
      if (pch) break;
    }
    
    dataLines = (char**) calloc(1, sizeof(char*));

    maxXML = 128; /* note actual number of tags is maxXML/2 */
    xmlRaw  = (char**) calloc(maxXML, sizeof(char*));
    xmlTags = (char**) calloc(maxXML, sizeof(char*));
    xmlVals = (char**) calloc(maxXML, sizeof(char*));
    for (i=0; i<maxXML; i++) {
      xmlRaw[i]  = (char*) calloc(512, sizeof(char));
      xmlTags[i] = (char*) calloc(512, sizeof(char));
      xmlVals[i] = (char*) calloc(512, sizeof(char));
    }

    /* figure out how many time steps exist in the file */
    nt = iL = 0;
    while (!feof(pFile)) { 
      fgets(buff, lineL, pFile);  /* read in line */
      pch = strstr(buff, "</Collection>");
      if (pch) break;
      nt++;
      dataLines     = (char**) realloc(dataLines, nt*sizeof(char*));
      dataLines[iL] = (char*)  calloc(lineL, sizeof(char));
      strcpy(dataLines[iL], buff);

      /* copy the current data line for XML processing */
      xmlLine = (char*) calloc(lineL, sizeof(char));
      strcpy(xmlLine, buff);

      /* split the string into tokens */
      rXML = 0;
      pch = strtok(buff, " =");
      while (pch != NULL) {
        rXML++;
        strcpy(xmlRaw[rXML-1], pch);
        pch = strtok(NULL, " =");
      }
      
      /* make sure spaces are there so that the first and last are
         not arguments */
      nXML = (rXML-2)/2;
      for (i=0, j=1; i<nXML; i++, j+=2) {
        strcpy(xmlTags[i], xmlRaw[j]);
        strncpy(xmlVals[i], &xmlRaw[j+1][1], strlen(xmlRaw[j+1])-2);
      }

      /* extract the current time */
      tLine = -1;
      for (j=0; j<nXML; j++) {
        if (!strcmp(xmlTags[j], "timestep"))
          tLine = atof(xmlVals[j]);
      }

      free(xmlLine);
      if (tLine >= cur_time)
        break;
      iL++;
    } /* end while over datalines */
    iCend = iCstart + nt;
    ex.fileNo = nt-1;

    if (tLine < cur_time) {
      /* write the new dataline */
      nt++;
      dataLines       = (char**) realloc(dataLines, nt*sizeof(char*));
      dataLines[nt-1] = (char*)  calloc(lineL, sizeof(char));
      ex.fileNo = nt-1;

      sprintf(buffer, "data/%s_%05i.vtm", prefix, ex.fileNo);
      VTKstrDataSet(dataLines[nt-1], 4, cur_time, "", "0", buffer);
    }
    
    /* free restart memory */
    for (i=0; i<maxXML; i++) {
      free(xmlRaw[i]); free(xmlTags[i]); free(xmlVals[i]);
    }
    free(xmlRaw); free(xmlTags); free(xmlVals);
  }
  else {
    nt        = 1;
    iL        = nt-1;
    dataLines     = (char**) calloc(1,     sizeof(char*));
    dataLines[iL] = (char*)  calloc(lineL, sizeof(char));
    iStart    = 0;
    ex.fileNo = iStart;
    ex.init = 1; /* now restart */

    sprintf(buffer, "data/%s_%05i.vtm", prefix, ex.fileNo);
    VTKstrDataSet(dataLines[iL], 4, cur_time, "", "0", buffer);
  }

  /* open the file for writing */
  pFile = fopen(fileName, "w");

  /* file info */
  fprintf(pFile, "<VTKFile type=\"Collection\" version=\"0.1\" ");
  if (ex.byteOrder == 0)
    fprintf(pFile, "byte_order=\"LittleEndian\"");
  else if (ex.byteOrder == 1)
    fprintf(pFile, "byte_order=\"BigEndian\"");
  fprintf(pFile, ">\n");

  /* Data collection */
  fprintf(pFile, "  <Collection>\n");
  for (i=0; i<nt; i++) {
    fprintf(pFile, "%s", dataLines[i]);
  }
  fprintf(pFile, "  </Collection>\n");

  fprintf(pFile, "</VTKFile>\n");
  fclose(pFile);

  for (i=0; i<nt; i++) free(dataLines[i]);
  free(dataLines);
  free(buff);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

/* write the multiblock data set stuff */
void VTKwriteMultiBlock(int nZone, char **zoneNames)
{
#if !RP_HOST
  int i, j, tab, nBlock;

  /* generate the filename */
  char fileName[1024];
  sprintf(fileName,"%s/data/%s_%05i.vtm",ex.dir, ex.prefix, ex.fileNo);
  Message0("  %-46s  %s\n", "Writing vtk multiblock file:",
           strrchr(fileName, '/')+1);

  FILE *mbF = NULL;
#if !PARALLEL 
  int myid      = 0;
  int node_zero = 0;
#endif /* !PARALLEL  */

  if (myid == node_zero) {
    mbF = fopen(fileName, "w");
    if (!mbF)
      Message0("Error: could not open file %s\n", fileName);
    
    /* write the header information */
    /* fprintf(mbF, "<?xml version=\"1.0\"?>\n"); */
    fprintf(mbF, "<VTKFile type=\"vtkMultiBlockDataSet\" ");
    fprintf(mbF, "version=\"1.0\" ");
    if (ex.byteOrder == 0)
      fprintf(mbF, "byte_order=\"LittleEndian\"");
    else if (ex.byteOrder == 1)
      fprintf(mbF, "byte_order=\"BigEndian\"");
    fprintf(mbF, ">\n");
    tab = 2;
    fprintf(mbF, "%*s<vtkMultiBlockDataSet>\n", tab, "");
    tab += 2;

    nBlock = 0;
    if (ex.dump_particle || ex.nSurf > 0) { /* put fluid region in own block */
      fprintf(mbF, "%*s<Block index=\"%i\" name=\"fluid\">\n", tab, "", nBlock);
      tab += 2;
    }

    /* write the individual unstructured data zones */
    for (i=0; i<nZone; i++) {
      fprintf(mbF, "%*s<DataSet index=\"%i\" ", tab, "", i);
      fprintf(mbF, "name=\"%s\" ", zoneNames[i]);
      fprintf(mbF, "file=\"%s-%s_%05i.vtu\" />\n", ex.prefix, zoneNames[i],
              ex.fileNo);
    }

    tab -= 2;
    if (ex.dump_particle || ex.nSurf > 0) {
      fprintf(mbF, "%*s</Block>\n", tab, "");
    }
   
    /* dump the spray */
    if (ex.dump_particle) {
      nBlock++;
      fprintf(mbF, "%*s<Block index=\"%i\" name=\"spray\">\n", tab, "", nBlock);
      tab += 2;
      for (i=0; i<ex.nInj; i++) {
        fprintf(mbF, "%*s<Block index=\"%i\" name=\"%s\">\n", tab, "", i,
                ex.inject_id[i]);
        tab += 2;

        for (j=0; j<ex.nNoz; j++) {
          fprintf(mbF, "%*s<DataSet index=\"%i\" ", tab, "", j);
          fprintf(mbF, "name=\"nozzle_%02i\" ", j+1);
          fprintf(mbF, "file=\"%s-spray-%s-%02i_%05i.vtp\" />\n", ex.prefix, 
                  ex.inject_id[i], j+1, ex.fileNo);
        }
        tab -= 2;
        fprintf(mbF, "%*s</Block>\n", tab, "");
      }

      tab -= 2;
      fprintf(mbF, "%*s</Block>\n", tab, "");
    } /* end particles */

    /* dump the walls */
    if (ex.nSurf > 0) {
      nBlock++;
      fprintf(mbF, "%*s<Block index=\"%i\" name=\"walls\">\n", tab, "", nBlock);
      tab += 2;

      for (i=0; i<ex.nSurf; i++) {
        fprintf(mbF, "%*s<DataSet index=\"%i\" ", tab, "", i);
        fprintf(mbF, "name=\"%s\" ", ex.surfs[i]);
        fprintf(mbF, "file=\"%s-%s_%05i.vtp\" />\n", ex.prefix, ex.surfs[i],
                ex.fileNo);
      }

      tab -= 2;
      fprintf(mbF, "%*s</Block>\n", tab, "");
    } /* end walls */

    /* dump the walls */
    if (ex.nCxSurf > 0) {
      nBlock++;
      fprintf(mbF, "%*s<Block index=\"%i\" name=\"surfaces\">\n",tab,"",nBlock);
      tab += 2;

      for (i=0; i<ex.nCxSurf; i++) {
        fprintf(mbF, "%*s<DataSet index=\"%i\" ", tab, "", i);
        fprintf(mbF, "name=\"%s\" ", ex.cxSurfs[i]);
        fprintf(mbF, "file=\"%s-%s_%05i.vtp\" />\n", ex.prefix, ex.cxSurfs[i],
                ex.fileNo);
      }

      tab -= 2;
      fprintf(mbF, "%*s</Block>\n", tab, "");
    } /* end walls */

    tab -= 2;
    fprintf(mbF, "%*s</vtkMultiBlockDataSet>\n", tab, "");
    fprintf(mbF, "</VTKFile>\n");

    fclose(mbF);
  }
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

int VTKsetZones(int nZ, char **zNames, Thread ***zT)
{
#if !RP_HOST
  Domain *d = Get_Domain(1);
  Thread *t, *tf;
  cell_t  c;
  face_t  f;
  Node   *v;
  int i, j, n, nf, iT, iP, nT, nC, nN, nCV;

  zT[0] = (Thread**) calloc(1, sizeof(Thread*));

  /* get the threads that correspond to the current zones */
  nT = 0;  
  i = 0;
  nC = nN = 0;
  thread_loop(t,d,iT) {
    if (FLUID_CELL_THREAD_P(t)) {
      /* check whether thread name exists */
      if (!strcmp(THREAD_NAME(t), "")) {
        Message0("Warning: thread names not set\n");
        /* assume zNames contains thread ids instead of names */
        /* if (!strcmp(zNames, atoi(THREAD_ID(t)))) { */
          
/*         } */
      }
      else {
        /* check whether thread belongs to zone */
        for (j=0; j<nZ; j++) {
          if (!strcmp(zNames[j], THREAD_NAME(t))) {
            nT++;
            iP = nT - 1;

            /* store the thread pointer */
            zT[0]     = (Thread**) realloc(zT[0], nT*sizeof(Thread*));
            zT[0][iP] = t;
          }
        } /* end loop over zone names */
      }
        
    } /* end fluid thread */
  } /* end loop over all threads */

  return nT;
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

int VTKsetSurfs(int nZ, char **zNames, Thread ***zT)
{
#if !RP_HOST
  Domain *d = Get_Domain(1);
  Thread *t, *tf;
  cell_t  c;
  face_t  f;
  Node   *v;
  int i, j, k, m, n, nf, iT, iP, nT, nC, nN, nCV, tID, iID, match;

  /* allocate memory for thread pointers */
  zT[0]   = (Thread**) calloc(1, sizeof(Thread*));

  /* get the threads that correspond to the current zones */
  nT = 0;  
  i = 0;
  nC = nN = 0;
  thread_loop(t,d,iT) {
    /* check whether thread belongs to zone */
    match = 0;
    for (j=0; j<nZ; j++) {
      if (!strcmp(zNames[j], "symmetry")) {
        if (THREAD_TYPE(t) == THREAD_F_SYMMETRIC)
          match = 1;
        
        /* check to see if already included elsewhere */
        for (k=0; k<ex.nAllSurfThreads; k++) {
          if (THREAD_ID(t) == ex.sThreadIDs[k])
            match = 0;
        }
      }
      else if (!strncmp(zNames[j], "id-", 3)) {
        iID = strcspn(zNames[j], "-") + 1;
        tID = atoi(&zNames[j][iID]);
        if (tID == THREAD_ID(t)) {
          match = 1;
        }
      }
      else if (!strcmp(zNames[j], THREAD_NAME(t))) 
        match = 1;
    } /* end loop over zone names */
    
      /* store thread if matched */
    if (match) {
      nT++;
      iP = nT - 1;
      
      /* store the thread pointer */
      zT[0]     = (Thread**) realloc(zT[0], nT*sizeof(Thread*));
      zT[0][iP] = t;
      
    }
  } /* end loop over all threads */

  return nT;
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

int VTKsetCx(int nZ, char **zNames, int **zT)
{
#if !RP_HOST
  int i,j,iP,match;

  /* allocate memory for thread pointers */
  zT[0]    = (int*) calloc(1, sizeof(int));
  zT[0][0] = -1;

  /* get the threads that correspond to the current zones */
  int nT = 0;  
  i = 0;
  for (i=0; i<ex.nCxRef; i++) {
    match = 0;
    for (j=0; j<nZ; j++) {
      if (!strcmp(zNames[j], ex.cxRefNames[i])) {
        match = 1;
      }
    } /* end loop over surface names */

    if (match) {
      nT++;
      iP = nT - 1;
      
      /* store the thread pointer */
      zT[0]     = (int*) realloc(zT[0], nT*sizeof(int));
      zT[0][iP] = ex.cxRefIds[i];
    }
  } /* end loop over all surfaces */

  return nT;
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void VTKwriteUnstructData(char *zone, int nZ, char **zNames)
{
#if !RP_HOST
  FILE *dF;
  int i, j;
  long long int offset;

#if !PARALLEL 
  myid          = 0;
  int node_zero = 0;
#endif /* !PARALLEL  */

  int n, nf, nC, nN, nCV, nFCV;
  int iC, iN, iCV, iFace;
  int totCells, totNodes, totCV, totFCV, nZt, iT, sizeFaceArray;
  int bytes;
  Thread *t, *tf, **Zt;
  cell_t  c;
  face_t  f;
  Node   *v;

  /* get the info about the zones */
  nZt = VTKsetZones(nZ, zNames, &Zt);

  /* generate the filename */
  char fileName[1024];
  sprintf(fileName,"%s/data/%s-%s_%05i.vtu",ex.dir,ex.prefix,zone,ex.fileNo);
  Message0("  %-36s (%3i/%3i)  %s\n", "Writing vtk unstructured data file:",nZt,
           nZ, strrchr(fileName, '/')+1);

  /* get the information about the cells/nodes on each processor and define
     the appropriate offsets for parallel reduction */
  int *nCells = (int*) calloc(compute_node_count, sizeof(int));
  int *nNodes = (int*) calloc(compute_node_count, sizeof(int));
  int *nCVs   = (int*) calloc(compute_node_count, sizeof(int));
  int *cOff   = (int*) calloc(compute_node_count, sizeof(int));
  int *cxOff  = (int*) calloc(compute_node_count, sizeof(int));
  int *nOff   = (int*) calloc(compute_node_count, sizeof(int));
  int *nxOff  = (int*) calloc(compute_node_count, sizeof(int));
  int *cvOff  = (int*) calloc(compute_node_count, sizeof(int));
  int *iWork  = (int*) calloc(compute_node_count, sizeof(int));

  /* support for polyhedra */
  int  nPtFtmp[128];
  int *nFaceCV = (int*) calloc(compute_node_count, sizeof(int));
  int *fcvOff  = (int*) calloc(compute_node_count, sizeof(int));

  /* mark nodes */
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_c_loop_int(c,t) {
      c_node_loop(c,t,n) {
        v = C_NODE(c,t,n);
        NODE_ITMP2(v) = -1;
      }
    } end_c_loop_int(c,t);
  } /* end loop over zone threads */

  /* get the total number of nodes and cells */
  nCells[myid] = nNodes[myid] = nCVs[myid] = nFaceCV[myid] = 0;
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_c_loop_int(c,t) { /* loop over all interior cells */
      i=0;
      c_face_loop(c,t,n) {  /* loop over this cell's faces */
        f  = C_FACE(c,t,n);        /* get the current face */
        tf = C_FACE_THREAD(c,t,n); /* get the associated face thread */
        nPtFtmp[i] = F_NNODES(f,tf);
        i++;
        if (PRINCIPAL_FACE_P(f,tf)) {
          f_node_loop(f,tf,nf) { /* loop over nodes of principal face */
            v = F_NODE(f,tf,nf);
            if (NODE_ITMP2(v) == -1) {
              NODE_ITMP2(v) = 1;
              nNodes[myid]++;
            }
          } /* end loop over face nodes */
        } /* end principal face */
      } /* end loop over cell faces */
      nCells[myid]++;
      nCVs[myid] += C_NNODES(c,t);
      nFaceCV[myid] += 1 + C_NFACES(c,t);
      for (i=0; i<C_NFACES(c,t); i++)
        nFaceCV[myid] += nPtFtmp[i];
    } end_c_loop_int(c,t);
  } /* end loop over zone threads */

#if PARALLEL
  PRF_GISUM(nCells, compute_node_count, iWork);
  for (i=0, totCells=0; i<compute_node_count; i++) totCells += nCells[i];

  PRF_GISUM(nNodes, compute_node_count, iWork);
  for (i=0, totNodes=0; i<compute_node_count; i++) totNodes += nNodes[i];

  PRF_GISUM(nCVs, compute_node_count, iWork);
  for (i=0, totCV=0; i<compute_node_count; i++)    totCV += nCVs[i];

  PRF_GISUM(nFaceCV, compute_node_count, iWork);
  for (i=0, totFCV=0; i<compute_node_count; i++)   totFCV += nFaceCV[i];

  nC = nCells[myid];
  /* find the offset for each compute node when transferring to the root */
  for (i=1, cOff[0]=0, cxOff[0]=0; i<compute_node_count; i++) {
    for (j=0, cOff[i]=0, cxOff[i]=0; j<i; j++) {
      cOff[i]  += nCells[j];
      cxOff[i] += 3*nCells[j];
    }
  }

  nN = nNodes[myid];
  /* find the offset for each compute node when transferring to the root */
  for (i=1, nOff[0]=0, nxOff[0]=0; i<compute_node_count; i++) {
    for (j=0, nOff[i]=0, nxOff[i]=0; j<i; j++) {
      nOff[i]  += nNodes[j];
      nxOff[i] += 3*nNodes[j];
    }
  }

  nCV = nCVs[myid];
  /* find the offset for each compute node when transferring to the root */
  for (i=1, cvOff[0]=0; i<compute_node_count; i++)
    for (j=0, cvOff[i]=0; j<i; j++)
      cvOff[i] += nCVs[j];

  nFCV = nFaceCV[myid];
  /* find the offset for each compute node when transferring to the root */
  for (i=1, fcvOff[0]=0; i<compute_node_count; i++)
    for (j=0, fcvOff[i]=0; j<i; j++)
      fcvOff[i] += nFaceCV[j];

  if (myid == node_zero) {
    nC   = totCells;
    nN   = totNodes;
    nCV  = totCV;
    nFCV = totFCV;
  }
#else /* SERIAL */
  totCells = nCells[myid];
  totNodes = nNodes[myid];
  totCV    = nCVs[myid];
  totFCV   = nFaceCV[myid];
  nC       = totCells;
  nN       = totNodes;
  nCV      = totCV;
  nFCV     = totFCV;
#endif /* PARALLEL */
  if (myid == node_zero) {
    dF = fopen(fileName, "w");
    if (!dF)
      Message0("Error: could not open file %s\n", fileName);
    
    /* write the header information */
    fprintf(dF, "<?xml version=\"1.0\"?>\n");
    fprintf(dF, "<VTKFile type=\"UnstructuredGrid\" ");
    fprintf(dF, "version=\"0.1\" ");
    if (ex.byteOrder == 0)
      fprintf(dF, "byte_order=\"LittleEndian\"");
    else if (ex.byteOrder == 1)
      fprintf(dF, "byte_order=\"BigEndian\"");
    fprintf(dF, ">\n");
    
    /* -----  mesh info ----- */
    fprintf(dF, "%2s<UnstructuredGrid>\n", "");
    fprintf(dF, "%4s<Piece NumberOfPoints=\"%i\" ", "", totNodes);
    fprintf(dF, "NumberOfCells=\"%i\">\n", totCells);
    
    offset = 0;
    /* ----- co-ordinate info ----- */
    fprintf(dF, "%6s<Points>\n", "");
    VTKinfoDataArray(dF,8,"Float32","coordinates",3,"appended",offset);
    offset += sizeof(int) + 3*totNodes*sizeof(float);
    fprintf(dF, "%6s</Points>\n", "");
    
    /* ----- cell info ----- */
    fprintf(dF, "%6s<Cells>\n", "");
    VTKinfoDataArray(dF,8,"Int32","connectivity",1,"appended",offset);
    offset += sizeof(int) + totCV*sizeof(int);
    VTKinfoDataArray(dF,8,"Int32","offsets",1,"appended",offset);
    offset += sizeof(int) + totCells*sizeof(int);
    VTKinfoDataArray(dF,8,"Int8","types",1,"appended",offset);
    offset += sizeof(int) + totCells*sizeof(char);
    if (ex.poly_mesh == 1) {
      VTKinfoDataArray(dF,8,"UInt32","faces",1,"appended",offset);
      offset += sizeof(int) + totFCV*sizeof(int);
      VTKinfoDataArray(dF,8,"UInt32","faceoffsets",1,"appended",offset);
      offset += sizeof(int) + totCells*sizeof(int);
    }
    fprintf(dF, "%6s</Cells>\n", "");
    
    /* ----- point data data ----- */
    fprintf(dF, "%6s<PointData>\n", "");
    /* put data output here */
    /* VTKinfoDataArray(dF,8,"Float32","dummy-points",1,"appended",offset); */
    /* offset += sizeof(int) + totNodes*sizeof(float); */
    fprintf(dF, "%6s</PointData>\n", "");
    
    /* ----- cell data ----- */
    fprintf(dF, "%6s<CellData>\n", "");
      /* velocity */
      VTKinfoDataArray(dF,8,"Float32","velocity",3,"appended",offset);
      offset += sizeof(int) + 3*totCells*sizeof(float);
      
      /* scalars */
      for (j=0; j<ex.nSc; j++) {
        VTKinfoDataArray(dF,8,"Float32",ex.sc_names[j],1,"appended",offset);
        offset += sizeof(int) + totCells*sizeof(float);
      }
      fprintf(dF, "%6s</CellData>\n", "");
      
      fprintf(dF, "%4s</Piece>\n", "");
      fprintf(dF, "%2s</UnstructuredGrid>\n", "");
      fprintf(dF, "%2s<AppendedData encoding=\"raw\">\n_", "");
  } /* end header write */

  /* ========== start the appended data ========== */

  int *nodeID   = (int*)   calloc(totNodes, sizeof(int));
  int *nodeIndx = (int*)   calloc(totNodes, sizeof(int));
  float *nodes  = (float*) calloc(3*nN,     sizeof(float));
  int *faceCV;
  if (ex.poly_mesh == 1)
    faceCV = (int*) calloc(nFCV, sizeof(int));
  
  /* fill all the node data */
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_c_loop_int(c,t) {
      c_node_loop(c,t,n) {
        v = C_NODE(c,t,n);
        NODE_ITMP2(v) = -1;
      }
    } end_c_loop_int(c,t);
  } /* end loop over threads */

  /* store the nodes */
  iN = 0;
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_c_loop_int(c,t) {
      c_face_loop(c,t,n) { /* loop over this cell's faces */
        f  = C_FACE(c,t,n);        /* get the current face */
        tf = C_FACE_THREAD(c,t,n); /* get the associated face thread */
        if (PRINCIPAL_FACE_P(f,tf)) {
          f_node_loop(f,tf,nf) {  /* loop over nodes of principal face */
            v = F_NODE(f,tf,nf);
            if (NODE_ITMP2(v) == -1) {
              /* store IDs and co-ordinates */
#if PARALLEL
              nodeID[iN]    = NODE_ID(v);
#else /* SERIAL */
              nodeID[iN]    = iN;
              NODE_ID(v)    = iN;
#endif /* PARALLEL */
              nodes[3*iN+0] = NODE_X(v);
              nodes[3*iN+1] = NODE_Y(v);
              nodes[3*iN+2] = NODE_Z(v);
              NODE_ITMP2(v) = 1;
              iN++;
            }
          } /* end loop over nodes */
        } /* end principal face */
      } /* end loop over cell faces */
    } end_c_loop_int(c,t);
  } /* end loop over threads */

#if PARALLEL /* reduce everything to root process */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++) {
      PRF_CRECV_INT(i,   &nodeID[nOff[i]], nNodes[i],   myid);
      PRF_CRECV_FLOAT(i, &nodes[nxOff[i]], 3*nNodes[i], myid);
    }
  }
  else {
    PRF_CSEND_INT(node_zero,  nodeID, nNodes[myid],   node_zero);
    PRF_CSEND_FLOAT(node_zero, nodes, 3*nNodes[myid], node_zero);
  }
#endif /* PARALLEL */

  /* write the data from the root */
  if (myid == node_zero) {
    /* index array in sorted order of node ident if necessary */
    indexIntArray(totNodes, nodeID, nodeIndx);

    /* write the data size */
    bytes = sizeof(int) + 3*totNodes*sizeof(float);
    writeDataInt(dF, 1, &bytes, 1);

    /* write the node coordinates */
    writeDataReal(dF, 3*totNodes, nodes, 3);
  }

#if PARALLEL 
  if (myid == node_zero) {
    /* send the sorted list and index out to the compute nodes */
    for (i=1; i<compute_node_count; i++) {
      PRF_CSEND_INT(i, nodeID,   totNodes, myid);
      PRF_CSEND_INT(i, nodeIndx, totNodes, myid);
    }
  }
  else {
    /* receive entire sorted list from root */
    PRF_CRECV_INT(node_zero, nodeID,   totNodes, node_zero);
    PRF_CRECV_INT(node_zero, nodeIndx, totNodes, node_zero);
  }
#endif /* PARALLEL */

  /* set all the tmp1 values to zero */
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_c_loop_int(c,t) {
      c_node_loop(c,t,n) {
        v = C_NODE(c,t,n);
        NODE_ITMP2(v) = -1;
      }
    } end_c_loop_int(c,t);
  } /* end loop over threads */

  /* store the order node numbering in temporary value used later in 
     unstructured connectivity write (VTKwriteUnstructGeom) */
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_c_loop_int(c,t) {
      c_node_loop(c,t,n) {
        v = C_NODE(c,t,n);
        if (NODE_ITMP2(v) == -1) {
          j = gridLocateIntIndx(totNodes,nodeID,nodeIndx,NODE_ID(v));
          /* account for the grid locator always returning the value
             to the left */
          if (NODE_ID(v) == nodeID[nodeIndx[totNodes-1]]) j++;

          /* store the correct node id for use in connectivity */
          NODE_ITMP2(v) = nodeIndx[j];
        }
      }
    } end_c_loop_int(c,t);
  } /* end loop over threads */

  char *types;
  int *offs;
  int *cv    = (int*) calloc(nCV, sizeof(int));  /* connectivity */
  int *nVert = (int*) calloc(nC,  sizeof(int));  /* number of vertices */
  int *cType = (int*) calloc(nC,  sizeof(int));  /* cell type */
  int *fOffs;
  if (ex.poly_mesh == 1)
    fOffs = (int*) calloc(nC,  sizeof(int));     /* face memory offsets */

  /* fill the connectivity array */
  iC = iCV = iFace = 0;
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_c_loop_int(c,t) {

      /* get the cell node numbers */
      iN = 0;
      c_node_loop(c,t,n) {
        v = C_NODE(c,t,n);
        cv[iCV+iN] = NODE_ITMP2(v);
        iN++;
      } /* end loop over nodes */
      
      /* store the number of vertices for this cell */
      nVert[iC] = C_NNODES(c,t);

      /* store the cell type */
      cType[iC] = C_TYPE(c,t);

      if (ex.poly_mesh == 1) {
        /* store the face information if polyhedron */
        faceCV[iFace] = C_NFACES(c,t);
        iFace++;
        c_face_loop(c,t,n) { /* loop over this cell's faces */
          f  = C_FACE(c,t,n);        /* get the current face */
          tf = C_FACE_THREAD(c,t,n); /* get the associated face thread */
          faceCV[iFace] = F_NNODES(f,tf);
          iFace++;
          f_node_loop(f,tf,nf) {
            v = F_NODE(f,tf,nf);
            faceCV[iFace] = NODE_ITMP2(v);
            iFace++;
          } /* end loop over face nodes*/
        } /* end loop over cell faces */
        fOffs[iC] = iFace;
      }
      
      iC++;
      iCV += C_NNODES(c,t);
    } end_c_loop_int(c,t);
  } /* end loop over threads */

  /* add the offset of the previous processors to the offset */
  if (ex.poly_mesh == 1) {
    for (i=0; i<nCells[myid]; i++)
      fOffs[i] += fcvOff[myid];
  }
 
#if PARALLEL /* reduce everything to root process */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++) {
      PRF_CRECV_INT(i, &cv[cvOff[i]],      nCVs[i],    myid);
      PRF_CRECV_INT(i, &nVert[cOff[i]],    nCells[i],  myid);
      PRF_CRECV_INT(i, &cType[cOff[i]],    nCells[i],  myid);
      if (ex.poly_mesh == 1) {
        PRF_CRECV_INT(i, &faceCV[fcvOff[i]], nFaceCV[i], myid);
        PRF_CRECV_INT(i, &fOffs[cOff[i]],    nCells[i],  myid);
      }
    }
  }
  else {
    PRF_CSEND_INT(node_zero, cv,     nCVs[myid],    node_zero);
    PRF_CSEND_INT(node_zero, nVert,  nCells[myid],  node_zero);
    PRF_CSEND_INT(node_zero, cType,  nCells[myid],  node_zero);
    if (ex.poly_mesh == 1) {
      PRF_CSEND_INT(node_zero, faceCV, nFaceCV[myid], node_zero);
      PRF_CSEND_INT(node_zero, fOffs,  nCells[myid],  node_zero);
    }
  }
#endif /* PARALLEL */

  /* write <Cells> */
  if (myid == node_zero) {
    offs  = (int*)  calloc(totCells, sizeof(int));
    types = (char*) calloc(totCells, sizeof(char));

    /* compute the offsets and get cell types */
    offs[0]  = nVert[0];
    types[0] = VTKcellType(cType[0]);
    for (j=1; j<totCells; j++) {
      offs[j]  = offs[j-1] + nVert[j];
      types[j] = VTKcellType(cType[j]);
    }
    
    /* cell connectivity */
    bytes = sizeof(int) + totCV*sizeof(int);
    writeDataInt(dF, 1, &bytes, 1);
    writeDataInt(dF, totCV, cv, totCV); 
    
    /* cell offsets */
    bytes = sizeof(int) + totCells*sizeof(int);
    writeDataInt(dF, 1, &bytes, 1);
    writeDataInt(dF, totCells, offs, totCells);
    
    /* cell types */
    bytes = sizeof(int) + totCells*sizeof(char);
    writeDataInt(dF, 1, &bytes, 1);
    fwrite(types, sizeof(char), totCells, dF);

    if (ex.poly_mesh == 1) {
      /* faces */
      bytes = sizeof(int) + totFCV*sizeof(int);
      writeDataInt(dF, 1, &bytes, 1);
      writeDataInt(dF, totFCV, faceCV, totFCV);
      
      /* face offsets */
      bytes = sizeof(int) + totCells*sizeof(int);
      writeDataInt(dF, 1, &bytes, 1);
      writeDataInt(dF, totCells, fOffs, totCells);
    }

    free(offs); free(types);
  } /* end <Cells> write */
  free(nodes); free(nodeID); free(nodeIndx);
  free(cv); free(nVert); free(cType);
  if (ex.poly_mesh == 1) {
    free(faceCV); free(fOffs);
  }

  /* ----- write <PointData> ----- */
  /* VTKwriteDummyData(totNodes, dF); */

  /* ----- write <CellData> ----- */
  /* vectors (velocity for now) */
  float *cVec = (float*) calloc(3*nC, sizeof(float));
  iC = 0;
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_c_loop_int(c,t) {
      cVec[iC*3+0] = (float)(C_U(c,t));
      cVec[iC*3+1] = (float)(C_V(c,t));
      cVec[iC*3+2] = (float)(C_W(c,t));
      iC++;
    } end_c_loop_int(c,t);
  } /* end loop over threads */

#if PARALLEL
  /* transfer all the cell information to the root */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++)
      PRF_CRECV_FLOAT(i, &cVec[cxOff[i]], 3*nCells[i], myid);
  }
  else {
    PRF_CSEND_FLOAT(node_zero, cVec, 3*nCells[myid], node_zero);
  }
#endif /* PARALLEL */
  
  if (myid == node_zero) {
    bytes = sizeof(int) + 3*totCells*sizeof(float);
    writeDataInt(dF, 1, &bytes, 1);
    writeDataReal(dF, 3*totCells, cVec, 1);
  }
  free(cVec);

  /* scalars */
  /* array for scalar storage */
  float *cSc = (float*) calloc(nC, sizeof(float));
  for (j=0; j<ex.nSc; j++) {
    iC = 0;
    for (iT=0; iT<nZt; iT++) {
      t = Zt[iT];
      if (t == NULL) break;
      begin_c_loop_int(c,t) {
        if (ex.sc_type[j] == 0)
          cSc[iC] = (float)(ABS_P(C_P(c,t),op_pres));        /* abs. pressure */
        else if (ex.sc_type[j] == 1)
          cSc[iC] = (float)(C_STORAGE_R(c,t,ex.sc_indx[j])); /* standard */
        else if (ex.sc_type[j] == 2)
          cSc[iC] = (float)(C_UDMI(c,t,ex.sc_indx[j]));      /* udm */
        else if (ex.sc_type[j] == 3)
          cSc[iC] = (float)(C_YI(c,t,ex.sc_indx[j]));        /* species */
        iC++;
      } end_c_loop_int(c,t);
    } /* end loop over threads */
    
#if PARALLEL
    /* transfer all the cell information to the root */
    if (myid == node_zero) {
      for (i=1; i<compute_node_count; i++)
        PRF_CRECV_FLOAT(i, &cSc[cOff[i]], nCells[i], myid);
    }
    else {
      PRF_CSEND_FLOAT(node_zero, cSc, nCells[myid], node_zero);
    }
#endif /* PARALLEL */

    if (myid == node_zero) {
      bytes = sizeof(int) + totCells*sizeof(float);
      writeDataInt(dF, 1, &bytes, 1);
      writeDataReal(dF, totCells, cSc, 1);
    }
  } /* end loop over cell scalars */
  free(cSc);

  /* ========== end the appended data ========== */
  if (myid == node_zero) {
    fprintf(dF, "\n%2s</AppendedData>\n", "");
    fprintf(dF, "</VTKFile>\n");
    fclose(dF);
  }

  free(nCells); free(nNodes); free(nCVs); free(nFaceCV); free(iWork);
  free(cOff); free(cxOff); free(nOff); free(nxOff); free(cvOff); free(fcvOff);
  free(Zt);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void VTKwriteParticleData(char *inj_id, int nNoz, Injection **streams)
{
#if !RP_HOST
  Injection *I;
  Particle  *p;
  char nozID[128];
  int i, j, k, bytes;
  long long int offset;
 
  /* generate the filename */
  char fileName[1000];
  Message0("  %-25s %-20s\n", "Writing vtk spray data:", inj_id);

  FILE *dF;
#if !PARALLEL 
  myid          = 0;
  int node_zero = 0;
#endif /* !PARALLEL  */

  /* parallel counting variables */
  int totPart, nP, iP;
  int *nPart = (int*) calloc(compute_node_count, sizeof(int));
  int *pOff  = (int*) calloc(compute_node_count, sizeof(int));
  int *pxOff = (int*) calloc(compute_node_count, sizeof(int));
  int *iWork = (int*) calloc(compute_node_count, sizeof(int));

  /* storage for output variables */
  int   *vertIds;
  float *pNodes, *pTemp, *pDiam, *pMass, *pVel;

  /* write the nozzle data files */
  for (k=0; k<nNoz; k++) {

    /* make sure everything is reset for this nozzle */
    for (i=0; i<compute_node_count; i++)
      nPart[i] = pOff[i] = pxOff[i] = iWork[i] = 0.0;

    /* get the injection for this nozzle */
    sprintf(nozID, "%s-%i", inj_id, k+1);
    I = Pick_Injection(nozID);
  
    /* create the file name */
    sprintf(fileName, "%s/data/%s-spray-%s-%02i_%05i.vtp", ex.dir, ex.prefix,
            inj_id, k+1, ex.fileNo);
    Message0("  %-35s  nozzle %02i: %s\n", "", k+1, strrchr(fileName, '/')+1);

    /* determine the number of particles for this nozzle */
    nPart[myid] = 0;
    loop(p, I->p) {
      nPart[myid]++; /* number of particles on this processor */
    }

#if PARALLEL
    /* get total number of particles across all processors */
    PRF_GISUM(nPart, compute_node_count, iWork);
    for (i=0, totPart=0; i<compute_node_count; i++)
      totPart += nPart[i];
    
    /* find the offset for each compute node when transferring to the root */
    for (i=1, pOff[0]=0, pxOff[0]=0; i<compute_node_count; i++) {
      for (j=0, pOff[i]=0, pxOff[i]=0; j<i; j++) {
        pOff[i]  +=   nPart[j];
        pxOff[i] += 3*nPart[j];
      }
    }

    /* local storage needed */
    nP = nPart[myid];
    if (myid == node_zero)
      nP = totPart;
#else /* SERIAL */
    totPart = nPart[myid];
    nP      = totPart;
#endif /* PARALLEL */

    /* write the head information */
    if (myid == node_zero) {
      /* open the file for writing */
      dF = fopen(fileName, "w");
      if (!dF)
        Message0("Error: could not open file %s\n", fileName);
      
      /* write the header information */
      fprintf(dF, "<?xml version=\"1.0\"?>\n");
      fprintf(dF, "<VTKFile type=\"PolyData\" ");
      fprintf(dF, "version=\"0.1\" ");
      if (ex.byteOrder == 0)
        fprintf(dF, "byte_order=\"LittleEndian\"");
      else if (ex.byteOrder == 1)
        fprintf(dF, "byte_order=\"BigEndian\"");
      fprintf(dF, ">\n");
      
      /* -----  mesh info ----- */
      fprintf(dF, "%2s<PolyData>\n", "");
      fprintf(dF, "%4s<Piece NumberOfPoints=\"%i\" ", "", totPart);
      fprintf(dF, "NumberOfVerts=\"%i\" ", totPart);
      fprintf(dF, "NumberOfLines=\"%i\" ", 0);
      fprintf(dF, "NumberOfStrips=\"%i\" ", 0);
      fprintf(dF, "NumberOfPolys=\"%i\">\n", 0);
      
      offset = 0;
      /* ----- co-ordinate info ----- */
      fprintf(dF, "%6s<Points>\n", "");
      VTKinfoDataArray(dF,8,"Float32","coordinates",3,"appended",offset);
      offset += sizeof(int) + 3*totPart*sizeof(float);
      fprintf(dF, "%6s</Points>\n", "");
      
      /* ----- write the info for vertices ----- */
      fprintf(dF, "%6s<Verts>\n", "");
      VTKinfoDataArray(dF,8,"Float32","vertices",1,"appended",offset);
      offset += sizeof(int) + totPart*sizeof(int);
      fprintf(dF, "%6s</Verts>\n", "");
      
      /* ----- point data data ----- */
      fprintf(dF, "%6s<PointData>\n", "");
      VTKinfoDataArray(dF,8,"Float32","velocity",3,"appended",offset);
      offset += sizeof(int) + 3*totPart*sizeof(float);
      VTKinfoDataArray(dF,8,"Float32","temperature",1,"appended",offset);
      offset += sizeof(int) + totPart*sizeof(float);
      VTKinfoDataArray(dF,8,"Float32","diameter",1,"appended",offset);
      offset += sizeof(int) + totPart*sizeof(float);
      VTKinfoDataArray(dF,8,"Float32","mass",1,"appended",offset);
      offset += sizeof(int) + totPart*sizeof(float);
      fprintf(dF, "%6s</PointData>\n", "");
      
      fprintf(dF, "%4s</Piece>\n", "");
      fprintf(dF, "%2s</PolyData>\n", "");
      fprintf(dF, "%2s<AppendedData encoding=\"raw\">\n_", "");
    }

    /* storage for co-ordinates */
    pNodes = (float*) calloc(3*nP, sizeof(float));
    pVel   = (float*) calloc(3*nP, sizeof(float));
    pTemp  = (float*) calloc(nP, sizeof(float));
    pDiam  = (float*) calloc(nP, sizeof(float));
    pMass  = (float*) calloc(nP, sizeof(float));

    /* write the co-ordinates */
    iP = 0;
    loop(p, I->p) {
      pNodes[3*iP]   = (float)(P_POS(p)[0]);
      pNodes[3*iP+1] = (float)(P_POS(p)[1]);
      pNodes[3*iP+2] = (float)(P_POS(p)[2]);
      pVel[3*iP]     = (float)(P_VEL(p)[0]);
      pVel[3*iP+1]   = (float)(P_VEL(p)[1]);
      pVel[3*iP+2]   = (float)(P_VEL(p)[2]);
      pTemp[iP]      = (float)(P_T(p));
      pDiam[iP]      = (float)(P_DIAM(p));
      pMass[iP]      = (float)(P_MASS(p));
      iP++;
    }

#if PARALLEL
    /* transfer all the cell information to the root */
    if (myid == node_zero) {
      for (i=1; i<compute_node_count; i++) {
        PRF_CRECV_FLOAT(i, &pNodes[pxOff[i]], 3*nPart[i], myid);
        PRF_CRECV_FLOAT(i, &pVel[pxOff[i]],   3*nPart[i], myid);
        PRF_CRECV_FLOAT(i, &pTemp[pOff[i]],     nPart[i], myid);
        PRF_CRECV_FLOAT(i, &pDiam[pOff[i]],     nPart[i], myid);
        PRF_CRECV_FLOAT(i, &pMass[pOff[i]],     nPart[i], myid);
      }
    }
    else {
      PRF_CSEND_FLOAT(node_zero, pNodes, 3*nPart[myid], node_zero);
      PRF_CSEND_FLOAT(node_zero, pVel,   3*nPart[myid], node_zero);
      PRF_CSEND_FLOAT(node_zero, pTemp,    nPart[myid], node_zero);
      PRF_CSEND_FLOAT(node_zero, pDiam,    nPart[myid], node_zero);
      PRF_CSEND_FLOAT(node_zero, pMass,    nPart[myid], node_zero);
    }
#endif /* PARALLEL */

    if (myid == node_zero) {
      /* write the nodes */
      bytes = sizeof(int) + 3*totPart*sizeof(float);
      writeDataInt(dF,  1, &bytes, 1);
      writeDataReal(dF, 3*totPart, pNodes, 3);
      
      /* also write to the vertices so that I can see them */
      vertIds = (int*) calloc(totPart, sizeof(int));
      bytes = sizeof(int) + totPart*sizeof(int);
      for (j=0; j<totPart; j++) vertIds[j] = j;
      writeDataInt(dF,  1, &bytes, 1);
      writeDataReal(dF, totPart, vertIds, 3);
      free(vertIds);

      /* write the velocity */
      bytes = sizeof(int) + 3*totPart*sizeof(float);
      writeDataInt(dF,  1, &bytes, 1);
      writeDataReal(dF, 3*totPart, pVel, 3);
      
      /* write the scalars */
      bytes = sizeof(int) + totPart*sizeof(float);
      writeDataInt(dF,  1, &bytes, 1);
      writeDataReal(dF, totPart, pTemp, 3);
      writeDataInt(dF,  1, &bytes, 1);
      writeDataReal(dF, totPart, pDiam, 3);
      writeDataInt(dF,  1, &bytes, 1);
      writeDataReal(dF, totPart, pMass, 3);

      fprintf(dF, "\n%2s</AppendedData>\n", "");
      fprintf(dF, "</VTKFile>\n");
      fclose(dF);
    }

    PRF_GSYNC();
    free(pNodes); free(pVel); free(pTemp); free(pDiam); free(pMass);
  } /* end loop over nozzles */

  PRF_GSYNC();
  free(nPart); free(pOff); free(pxOff); free(iWork);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void VTKwriteWallSurf(char *zone, int nZ, char **zNames)
{
#if !RP_HOST
  FILE *dF;
  int i, j, k;
  long long int offset;


#if !PARALLEL 
  myid          = 0;
  int node_zero = 0;
#endif /* !PARALLEL  */

  int n, nf, nF, nN, nCV;
  int iF, iN, iCV;
  int totFaces, totNodes, totCV, nZt, iT;
  int bytes;
  Thread *t, *tf, **Zt;
  face_t  f;
  Node   *v;

  /* get the info about the zones */
  nZt = VTKsetSurfs(nZ, zNames, &Zt);

  /* store the thread ids of the ones just found */
  k = ex.nAllSurfThreads;
  ex.nAllSurfThreads += nZt;
  ex.sThreadIDs = (int*) realloc(ex.sThreadIDs, ex.nAllSurfThreads*sizeof(int));
  for (i=k, j=0; i<ex.nAllSurfThreads; i++, j++) {
    ex.sThreadIDs[i] = THREAD_ID(Zt[j]);
  }

  /* generate the filename */
  char fileName[1024];
  sprintf(fileName,"%s/data/%s-%s_%05i.vtp",ex.dir,ex.prefix,zone,ex.fileNo);
  Message0("  %-36s (%3i/%3i)  %s\n", "Writing vtk poly surface file:", nZt,nZ, 
           strrchr(fileName, '/')+1);

  /* get the information about the cells/nodes on each processor and define
     the appropriate offsets for parallel reduction */
  int *nFaces = (int*) calloc(compute_node_count, sizeof(int));
  int *nNodes = (int*) calloc(compute_node_count, sizeof(int));
  int *nCVs   = (int*) calloc(compute_node_count, sizeof(int));
  int *fOff   = (int*) calloc(compute_node_count, sizeof(int));
  int *nOff   = (int*) calloc(compute_node_count, sizeof(int));
  int *nxOff  = (int*) calloc(compute_node_count, sizeof(int));
  int *cvOff  = (int*) calloc(compute_node_count, sizeof(int));
  int *iWork  = (int*) calloc(compute_node_count, sizeof(int));

  /* mark nodes */
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_f_loop_int(f,t) {
      f_node_loop(f,t,n) {
        v = F_NODE(f,t,n);
        NODE_ITMP2(v) = -1;
      }
    } end_f_loop_int(f,t);
  } /* end loop over zone threads */

  /* get the total number of nodes and cells */
  nFaces[myid] = nNodes[myid] = nCVs[myid] = 0;
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_f_loop_int(f,t) { /* loop over all interior cells */
      if (PRINCIPAL_FACE_P(f,t)) {
        f_node_loop(f,t,nf) { /* loop over nodes of principal face */
          v = F_NODE(f,t,nf);
          if (NODE_ITMP2(v) == -1) {
            NODE_ITMP2(v) = 1;
            nNodes[myid]++;
          }
        } /* end loop over face nodes */
        nFaces[myid]++;
        nCVs[myid] += F_NNODES(f,t);
      } /* end principal face */
    } end_f_loop_int(f,t);
  } /* end loop over zone threads */
  
#if PARALLEL
  PRF_GISUM(nFaces, compute_node_count, iWork);
  for (i=0, totFaces=0; i<compute_node_count; i++) totFaces += nFaces[i];

  PRF_GISUM(nNodes, compute_node_count, iWork);
  for (i=0, totNodes=0; i<compute_node_count; i++) totNodes += nNodes[i];

  PRF_GISUM(nCVs, compute_node_count, iWork);
  for (i=0, totCV=0; i<compute_node_count; i++)    totCV    += nCVs[i];

  nF = nFaces[myid];
  /* find the offset for each compute node when transferring to the root */
  for (i=1, fOff[0]=0; i<compute_node_count; i++)
    for (j=0, fOff[i]=0; j<i; j++)
      fOff[i] += nFaces[j];

  nN = nNodes[myid];
  /* find the offset for each compute node when transferring to the root */
  for (i=1, nOff[0]=0, nxOff[0]=0; i<compute_node_count; i++) {
    for (j=0, nOff[i]=0, nxOff[i]=0; j<i; j++) {
      nOff[i]  += nNodes[j];
      nxOff[i] += 3*nNodes[j];
    }
  }

  nCV = nCVs[myid];
  /* find the offset for each compute node when transferring to the root */
  for (i=1, cvOff[0]=0; i<compute_node_count; i++)
    for (j=0, cvOff[i]=0; j<i; j++)
      cvOff[i] += nCVs[j];

  if (myid == node_zero) {
    nF  = totFaces;
    nN  = totNodes;
    nCV = totCV;
  }
#else /* SERIAL */
  totFaces = nFaces[myid];
  totNodes = nNodes[myid];
  totCV    = nCVs[myid];
  nF       = totFaces;
  nN       = totNodes;
  nCV      = totCV;
#endif /* PARALLEL */

  if (myid == node_zero) {
    dF = fopen(fileName, "w");
    if (!dF)
      Message0("Error: could not open file %s\n", fileName);
    
    /* write the header information */
    fprintf(dF, "<?xml version=\"1.0\"?>\n");
    fprintf(dF, "<VTKFile type=\"PolyData\" ");
    fprintf(dF, "version=\"0.1\" ");
    if (ex.byteOrder == 0)
      fprintf(dF, "byte_order=\"LittleEndian\"");
    else if (ex.byteOrder == 1)
      fprintf(dF, "byte_order=\"BigEndian\"");
    fprintf(dF, ">\n");
    
    /* -----  mesh info ----- */
    fprintf(dF, "%2s<PolyData>\n", "");
    fprintf(dF, "%4s<Piece NumberOfPoints=\"%i\" ", "", totNodes);
    fprintf(dF, "NumberOfVerts=\"%i\" ", 0);
    fprintf(dF, "NumberOfLines=\"%i\" ", 0);
    fprintf(dF, "NumberOfStrips=\"%i\" ", 0);
    fprintf(dF, "NumberOfPolys=\"%i\">\n", totFaces);
    
    offset = 0;
    /* ----- co-ordinate info ----- */
    fprintf(dF, "%6s<Points>\n", "");
    VTKinfoDataArray(dF,8,"Float32","coordinates",3,"appended",offset);
    offset += sizeof(int) + 3*totNodes*sizeof(float);
    fprintf(dF, "%6s</Points>\n", "");
    
    /* ----- write poly connectivity info ----- */
    fprintf(dF, "%6s<Polys>\n", "");
    VTKinfoDataArray(dF,8,"Int32","connectivity",1,"appended",offset);
    offset += sizeof(int) + totCV*sizeof(int);
    VTKinfoDataArray(dF,8,"Int32","offsets",1,"appended",offset);
    offset += sizeof(int) + totFaces*sizeof(int);
    fprintf(dF, "%6s</Polys>\n", "");
    
    /* ----- point data data ----- */
    fprintf(dF, "%6s<PointData>\n", "");
    fprintf(dF, "%6s</PointData>\n", "");

    /* ----- face data here ----- */
    fprintf(dF, "%6s<CellData>\n", "");
    VTKinfoDataArray(dF,8,"Float32","heat_flux",1,"appended",offset);
    offset += sizeof(int) + totFaces*sizeof(float);
    if (dpm_par.wall_film_active) {
      VTKinfoDataArray(dF,8,"Float32","wall_film_height",1,"appended",offset);
      offset += sizeof(int) + totFaces*sizeof(float);
      VTKinfoDataArray(dF,8,"Float32","wall_film_mass",1,"appended",offset);
      offset += sizeof(int) + totFaces*sizeof(float);
      VTKinfoDataArray(dF,8,"Float32","wall_film_temp",1,"appended",offset);
      offset += sizeof(int) + totFaces*sizeof(float);
    }
    fprintf(dF, "%6s</CellData>\n", "");
    
    /* ----- cell data ----- */
         
    fprintf(dF, "%4s</Piece>\n", "");
    fprintf(dF, "%2s</PolyData>\n", "");
    fprintf(dF, "%2s<AppendedData encoding=\"raw\">\n_", "");
  } /* end header write */

  /* ========== start the appended data ========== */

  int *nodeID   = (int*)   calloc(totNodes, sizeof(int));
  int *nodeIndx = (int*)   calloc(totNodes, sizeof(int));
  float *nodes  = (float*) calloc(3*nN,     sizeof(float));
  
  /* mark nodes */
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_f_loop_int(f,t) {
      f_node_loop(f,t,n) {
        v = F_NODE(f,t,n);
        NODE_ITMP2(v) = -1;
      }
    } end_f_loop_int(f,t);
  } /* end loop over zone threads */

  /* store the nodes */
  iN = 0;
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_f_loop_int(f,t) {
      if (PRINCIPAL_FACE_P(f,t)) {
        f_node_loop(f,t,nf) {  /* loop over nodes of principal face */
          v = F_NODE(f,t,nf);
          if (NODE_ITMP2(v) == -1) {
            /* store IDs and co-ordinates */
#if PARALLEL
            nodeID[iN]    = NODE_ID(v);
#else /* SERIAL */
            nodeID[iN]    = iN;
            NODE_ID(v)    = iN;
#endif /* PARALLEL */
            nodes[3*iN+0] = NODE_X(v);
            nodes[3*iN+1] = NODE_Y(v);
            nodes[3*iN+2] = NODE_Z(v);
            NODE_ITMP2(v) = 1;
            iN++;
          }
        } /* end loop over nodes */
      } /* end principal face */
    } end_f_loop_int(f,t);
  } /* end loop over threads */

#if PARALLEL /* reduce everything to root process */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++) {
      PRF_CRECV_INT(i,   &nodeID[nOff[i]], nNodes[i],   myid);
      PRF_CRECV_FLOAT(i, &nodes[nxOff[i]], 3*nNodes[i], myid);
    }
  }
  else {
    PRF_CSEND_INT(node_zero,  nodeID, nNodes[myid],   node_zero);
    PRF_CSEND_FLOAT(node_zero, nodes, 3*nNodes[myid], node_zero);
  }
#endif /* PARALLEL */

  /* write the data from the root */
  if (myid == node_zero) {
    /* index array in sorted order of node ident if necessary */
    indexIntArray(totNodes, nodeID, nodeIndx);

    /* write the data size */
    bytes = sizeof(int) + 3*totNodes*sizeof(float);
    writeDataInt(dF, 1, &bytes, 1);

    /* write the node coordinates */
    writeDataReal(dF, 3*totNodes, nodes, 3);
  }

#if PARALLEL 
  if (myid == node_zero) {
    /* send the sorted list and index out to the compute nodes */
    for (i=1; i<compute_node_count; i++) {
      PRF_CSEND_INT(i, nodeID,   totNodes, myid);
      PRF_CSEND_INT(i, nodeIndx, totNodes, myid);
    }
  }
  else {
    /* receive entire sorted list from root */
    PRF_CRECV_INT(node_zero, nodeID,   totNodes, node_zero);
    PRF_CRECV_INT(node_zero, nodeIndx, totNodes, node_zero);
  }
#endif /* PARALLEL */

  /* mark nodes */
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_f_loop_int(f,t) {
      f_node_loop(f,t,n) {
        v = F_NODE(f,t,n);
        NODE_ITMP2(v) = -1;
      }
    } end_f_loop_int(f,t);
  } /* end loop over zone threads */

  /* store the order node numbering in temporary value used later in 
     unstructured connectivity write (VTKwriteUnstructGeom) */
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_f_loop_int(f,t) {
      if (PRINCIPAL_FACE_P(f,t)) {
        f_node_loop(f,t,n) {
          v = F_NODE(f,t,n);
          if (NODE_ITMP2(v) == -1) {
            j = gridLocateIntIndx(totNodes,nodeID,nodeIndx,NODE_ID(v));
            /* account for the grid locator always returning the value
               to the left */
            if (NODE_ID(v) == nodeID[nodeIndx[totNodes-1]]) j++;
            
            /* store the correct node id for use in connectivity */
            NODE_ITMP2(v) = nodeIndx[j];
          }
        } /* end loop over face nodes */
      } /* end principal face */
    } end_f_loop_int(f,t);
  } /* end loop over threads */

  int   *offs;
  int   *cv    = (int*)   calloc(nCV, sizeof(int));
  int   *nVert = (int*)   calloc(nF,  sizeof(int));
  float *hFlux = (float*) calloc(nF, sizeof(float));
  float *wfHeight, *wfMass, *wfTemp;

  if (dpm_par.wall_film_active) {
    wfHeight = (float*) calloc(nF, sizeof(float));
    wfMass   = (float*) calloc(nF, sizeof(float));
    wfTemp   = (float*) calloc(nF, sizeof(float));
  }

  /* fill the connectivity array */
  iF = iCV = 0;
  for (iT=0; iT<nZt; iT++) {
    t = Zt[iT];
    if (t == NULL) break;
    begin_f_loop_int(f,t) {
      if (PRINCIPAL_FACE_P(f,t)) {
        
        /* get the cell node numbers */
        iN = 0;
        f_node_loop(f,t,n) {
          v = F_NODE(f,t,n);
          cv[iCV+iN] = NODE_ITMP2(v);
          iN++;
        } /* end loop over nodes */
        
        /* store the number of vertices for this cell */
        nVert[iF] = F_NNODES(f,t);

        /* store face scalars */
        if (THREAD_TYPE(t) != THREAD_F_SYMMETRIC)
          hFlux[iF] = BOUNDARY_HEAT_FLUX(f,t);

        if (dpm_par.wall_film_active && WALL_FILM_THREAD_P(t)) {
          wfHeight[iF] = F_WALL_FILM_HEIGHT(f,t);
          wfMass[iF]   = F_WALL_FILM_MASS(f,t);
          wfTemp[iF]   = F_WALL_FILM_TEMP(f,t);
        }
        
        iF++;
        iCV += F_NNODES(f,t);
      } /* end principal face */
    } end_f_loop_int(f,t);
  } /* end loop over threads */

#if PARALLEL /* reduce everything to root process */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++) {
      PRF_CRECV_INT(i, &cv[cvOff[i]],   nCVs[i],   myid);
      PRF_CRECV_INT(i, &nVert[fOff[i]], nFaces[i], myid);
      PRF_CRECV_FLOAT(i, &hFlux[fOff[i]], nFaces[i], myid);
      if (dpm_par.wall_film_active) {
        PRF_CRECV_FLOAT(i, &wfHeight[fOff[i]], nFaces[i], myid);
        PRF_CRECV_FLOAT(i, &wfMass[fOff[i]], nFaces[i], myid);
        PRF_CRECV_FLOAT(i, &wfTemp[fOff[i]], nFaces[i], myid);
      }
    }
  }
  else {
    PRF_CSEND_INT(node_zero, cv,    nCVs[myid],   node_zero);
    PRF_CSEND_INT(node_zero, nVert, nFaces[myid], node_zero);
    PRF_CSEND_FLOAT(node_zero, hFlux, nFaces[myid], node_zero);
    if (dpm_par.wall_film_active) {
      PRF_CSEND_FLOAT(node_zero, wfHeight, nFaces[myid], node_zero);
      PRF_CSEND_FLOAT(node_zero, wfMass, nFaces[myid], node_zero);
      PRF_CSEND_FLOAT(node_zero, wfTemp, nFaces[myid], node_zero);
    }
  }
#endif /* PARALLEL */

  /* write <Cells> */
  if (myid == node_zero) {
    offs  = (int*)  calloc(totFaces, sizeof(int));

    /* compute the offsets and get cell types */
    offs[0]  = nVert[0];
    for (j=1; j<totFaces; j++) {
      offs[j]  = offs[j-1] + nVert[j];
    }
    
    /* cell connectivity */
    bytes = sizeof(int) + totCV*sizeof(int);
    writeDataInt(dF, 1, &bytes, 1);
    writeDataInt(dF, totCV, cv, totCV); 
    
    /* cell offsets */
    bytes = sizeof(int) + totFaces*sizeof(int);
    writeDataInt(dF, 1, &bytes, 1);
    writeDataInt(dF, totFaces, offs, totFaces);

    free(offs);

    /* cell data */
    bytes = sizeof(int) + totFaces*sizeof(float);
    writeDataInt(dF, 1, &bytes, 1);
    writeDataInt(dF, totFaces, hFlux, totFaces);

    if (dpm_par.wall_film_active) {
      bytes = sizeof(int) + totFaces*sizeof(float);
      writeDataInt(dF, 1, &bytes, 1);
      writeDataInt(dF, totFaces, wfHeight, totFaces);

      bytes = sizeof(int) + totFaces*sizeof(float);
      writeDataInt(dF, 1, &bytes, 1);
      writeDataInt(dF, totFaces, wfMass, totFaces);

      bytes = sizeof(int) + totFaces*sizeof(float);
      writeDataInt(dF, 1, &bytes, 1);
      writeDataInt(dF, totFaces, wfTemp, totFaces);
    }

  } /* end <Cells> write */
  free(nodes); free(nodeID); free(nodeIndx);
  free(cv); free(nVert);
  free(hFlux);
  if (dpm_par.wall_film_active) {
    free(wfHeight); free(wfMass); free(wfTemp);
  }

  /* ========== end the appended data ========== */
  if (myid == node_zero) {
    fprintf(dF, "\n%2s</AppendedData>\n", "");
    fprintf(dF, "</VTKFile>\n");
    fclose(dF);
  }

  free(nFaces); free(nNodes); free(nCVs); free(iWork);
  free(fOff); free(nOff); free(nxOff); free(cvOff);
  free(Zt);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void VTKwriteCxSurf(char *zone, int nZ, char **zNames)
{
#if !RP_HOST
  FILE *dF;
  int i, j, k;
  long long int offset;

#if !PARALLEL 
  myid          = 0;
  int node_zero = 0;
#endif /* !PARALLEL  */

  int n, nf, nF, nN, nCV;
  int iF, iN, iCV;
  int totFaces, totNodes, totCV, nZt, iT;
  int bytes;
  int sid, *Zt;
  Thread *t, *tf;
  Surface *s;
  cell_t  c;
  face_t  f;
  Node   *v;

  /* get the info about the zones */
  nZt = VTKsetCx(nZ, zNames, &Zt);

  /* store the thread ids of the ones just found */
  k = ex.nAllCxThreads;
  ex.nAllCxThreads += nZt;
  ex.cxThreadIDs = (int*) realloc(ex.cxThreadIDs,ex.nAllCxThreads*sizeof(int));
  for (i=k, j=0; i<ex.nAllCxThreads; i++, j++) {
    ex.cxThreadIDs[i] = Zt[i];
  }

  /* generate the filename */
  char fileName[1024];
  sprintf(fileName,"%s/data/%s-%s_%05i.vtp",ex.dir,ex.prefix,zone,ex.fileNo);
  Message0("  %-36s (%3i/%3i)  %s\n", "Writing vtk poly surface file:", nZt,nZ, 
           strrchr(fileName, '/')+1);

  /* get the information about the cells/nodes on each processor and define
     the appropriate offsets for parallel reduction */
  int *nFaces = (int*) calloc(compute_node_count, sizeof(int));
  int *nNodes = (int*) calloc(compute_node_count, sizeof(int));
  int *nCVs   = (int*) calloc(compute_node_count, sizeof(int));
  int *fOff   = (int*) calloc(compute_node_count, sizeof(int));
  int *fxOff  = (int*) calloc(compute_node_count, sizeof(int));
  int *nOff   = (int*) calloc(compute_node_count, sizeof(int));
  int *nxOff  = (int*) calloc(compute_node_count, sizeof(int));
  int *cvOff  = (int*) calloc(compute_node_count, sizeof(int));
  int *iWork  = (int*) calloc(compute_node_count, sizeof(int));

  /* get the number of faces, nodes, and size of connectivity */
  double vv[3];
  nFaces[myid] = nNodes[myid] = nCVs[myid] = 0;
  for (iT=0; iT<nZt; iT++) {
    sid = Zt[iT];
    if (sid < 0) break;
    s = SurfaceList + sid;
    nFaces[myid] += s->nf;
    nNodes[myid] += s->np;
    /* this list of connectivity information somehow contains an additional
       integer on top of the number of nodes.  I believe it is the id of the 
       correct node itself */
    /* nCVs[myid]   += s->nfl; */

    /* loop over the facets of this surface */
    for (k=0, i=0; k<s->nf; i+=s->facets[i]+1,k++) {
      /* make sure it is at least a tri (i.e. not a line) */
      if (s->facets[i] >= 3) {
        nCVs[myid] += s->facets[i];
      }
    }
  } /* end loop over zone threads */

#if PARALLEL
  PRF_GISUM(nFaces, compute_node_count, iWork);
  for (i=0, totFaces=0; i<compute_node_count; i++) totFaces += nFaces[i];

  PRF_GISUM(nNodes, compute_node_count, iWork);
  for (i=0, totNodes=0; i<compute_node_count; i++) totNodes += nNodes[i];

  PRF_GISUM(nCVs, compute_node_count, iWork);
  for (i=0, totCV=0; i<compute_node_count; i++)    totCV    += nCVs[i];

  nF = nFaces[myid];
  /* find the offset for each compute node when transferring to the root */
  for (i=1, fOff[0]=0, fxOff[0]=0; i<compute_node_count; i++) {
    for (j=0, fOff[i]=0, fxOff[i]=0; j<i; j++) {
      fOff[i]  += nFaces[j];
      fxOff[i] += 3*nFaces[j];
    }
  }

  nN = nNodes[myid];
  /* find the offset for each compute node when transferring to the root */
  for (i=1, nOff[0]=0, nxOff[0]=0; i<compute_node_count; i++) {
    for (j=0, nOff[i]=0, nxOff[i]=0; j<i; j++) {
      nOff[i]  += nNodes[j];
      nxOff[i] += 3*nNodes[j];
    }
  }

  nCV = nCVs[myid];
  /* find the offset for each compute node when transferring to the root */
  for (i=1, cvOff[0]=0; i<compute_node_count; i++)
    for (j=0, cvOff[i]=0; j<i; j++)
      cvOff[i] += nCVs[j];

  if (myid == node_zero) {
    nF  = totFaces;
    nN  = totNodes;
    nCV = totCV;
  }
#else /* SERIAL */
  totFaces = nFaces[myid];
  totNodes = nNodes[myid];
  totCV    = nCVs[myid];
  nF       = totFaces;
  nN       = totNodes;
  nCV      = totCV;
#endif /* PARALLEL */

  if (myid == node_zero) {
    dF = fopen(fileName, "w");
    if (!dF)
      Message0("Error: could not open file %s\n", fileName);
    
    /* write the header information */
    fprintf(dF, "<?xml version=\"1.0\"?>\n");
    fprintf(dF, "<VTKFile type=\"PolyData\" ");
    fprintf(dF, "version=\"0.1\" ");
    if (ex.byteOrder == 0)
      fprintf(dF, "byte_order=\"LittleEndian\"");
    else if (ex.byteOrder == 1)
      fprintf(dF, "byte_order=\"BigEndian\"");
    fprintf(dF, ">\n");
    
    /* -----  mesh info ----- */
    fprintf(dF, "%2s<PolyData>\n", "");
    fprintf(dF, "%4s<Piece NumberOfPoints=\"%i\" ", "", totNodes);
    fprintf(dF, "NumberOfVerts=\"%i\" ", 0);
    fprintf(dF, "NumberOfLines=\"%i\" ", 0);
    fprintf(dF, "NumberOfStrips=\"%i\" ", 0);
    fprintf(dF, "NumberOfPolys=\"%i\">\n", totFaces);
    
    offset = 0;
    /* ----- co-ordinate info ----- */
    fprintf(dF, "%6s<Points>\n", "");
    VTKinfoDataArray(dF,8,"Float32","coordinates",3,"appended",offset);
    offset += sizeof(int) + 3*totNodes*sizeof(float);
    fprintf(dF, "%6s</Points>\n", "");
    
    /* ----- write poly connectivity info ----- */
    fprintf(dF, "%6s<Polys>\n", "");
    VTKinfoDataArray(dF,8,"Int32","connectivity",1,"appended",offset);
    offset += sizeof(int) + totCV*sizeof(int);
    VTKinfoDataArray(dF,8,"Int32","offsets",1,"appended",offset);
    offset += sizeof(int) + totFaces*sizeof(int);
    fprintf(dF, "%6s</Polys>\n", "");
    
    /* ----- point data data ----- */
    fprintf(dF, "%6s<PointData>\n", "");
    fprintf(dF, "%6s</PointData>\n", "");

    /* ----- face data here ----- */
    fprintf(dF, "%6s<CellData>\n", "");
    VTKinfoDataArray(dF,8,"Float32","Velocity",3,"appended",offset);
    offset += sizeof(int) + 3*totFaces*sizeof(float);
    for (j=0; j<ex.nSc; j++) {
      VTKinfoDataArray(dF,8,"Float32",ex.sc_names[j],1,"appended",offset);
      offset += sizeof(int) + totFaces*sizeof(float);
    }
    fprintf(dF, "%6s</CellData>\n", "");
    
    fprintf(dF, "%4s</Piece>\n", "");
    fprintf(dF, "%2s</PolyData>\n", "");
    fprintf(dF, "%2s<AppendedData encoding=\"raw\">\n_", "");
  } /* end header write */

  /* ========== start the appended data ========== */

  int *nodeID   = (int*)   calloc(totNodes, sizeof(int));
  int *nodeIndx = (int*)   calloc(totNodes, sizeof(int));
  float *nodes  = (float*) calloc(3*nN,     sizeof(float));
  
  iN = 0;
  for (iT=0; iT<nZt; iT++) {
    sid = Zt[iT];
    if (sid < 0) break;
    s = SurfaceList + sid;

    for (k=0; k<s->np; k++) {
      nodeID[iN]    = nOff[myid]+iN;
      nodes[3*iN+0] = Get_Surface_Point_Coord(&(s->points[k]), X_DIM);
      nodes[3*iN+1] = Get_Surface_Point_Coord(&(s->points[k]), Y_DIM);
      nodes[3*iN+2] = Get_Surface_Point_Coord(&(s->points[k]), Z_DIM);
      iN++;
    }
  } /* end loop over zone threads */

#if PARALLEL /* reduce everything to root process */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++) {
      PRF_CRECV_INT(i,   &nodeID[nOff[i]], nNodes[i],   myid);
      PRF_CRECV_FLOAT(i, &nodes[nxOff[i]], 3*nNodes[i], myid);
    }
  }
  else {
    PRF_CSEND_INT(node_zero,  nodeID, nNodes[myid],   node_zero);
    PRF_CSEND_FLOAT(node_zero, nodes, 3*nNodes[myid], node_zero);
  }
#endif /* PARALLEL */

  /* write the data from the root */
  if (myid == node_zero) {
    /* index array in sorted order of node ident if necessary */
    indexIntArray(totNodes, nodeID, nodeIndx);

    /* TODO: should check for duplicate nodes */
    /* for (i=0; i<totNodes-1; i++) { */
    /*   for (j=i+1; j<totNodes; j++) { */
    /*     if ( (ABS_D(nodes[3*i+0] - nodes[3*j+0]) < 1.0e-9) && */
    /*          (ABS_D(nodes[3*i+1] - nodes[3*j+1]) < 1.0e-9) && */
    /*          (ABS_D(nodes[3*i+2] - nodes[3*j+2]) < 1.0e-9) ) { */
    /*       Message("%i: Nodes the same at %i/%i\n",myid,i,j); */
    /*     } */
    /*   } */
    /* } */

    /* write the data size */
    bytes = sizeof(int) + 3*totNodes*sizeof(float);
    writeDataInt(dF, 1, &bytes, 1);

    /* write the node coordinates */
    writeDataReal(dF, 3*totNodes, nodes, 3);
  }

#if PARALLEL 
  if (myid == node_zero) {
    /* send the sorted list and index out to the compute nodes */
    for (i=1; i<compute_node_count; i++) {
      PRF_CSEND_INT(i, nodeID,   totNodes, myid);
      PRF_CSEND_INT(i, nodeIndx, totNodes, myid);
    }
  }
  else {
    /* receive entire sorted list from root */
    PRF_CRECV_INT(node_zero, nodeID,   totNodes, node_zero);
    PRF_CRECV_INT(node_zero, nodeIndx, totNodes, node_zero);
  }
#endif /* PARALLEL */

  int   *offs;
  int   *cv    = (int*)   calloc(nCV,  sizeof(int));
  int   *nVert = (int*)   calloc(nF,   sizeof(int));
  float *fVel  = (float*) calloc(3*nF, sizeof(float));

  /* fill the connectivity array */
  iF = iCV = iN = 0;
  for (iT=0; iT<nZt; iT++) {
    sid = Zt[iT];
    if (sid < 0) break;
    s = SurfaceList + sid;

    /* loop over the facets of this surface */
    for (k=0, i=0; k<s->nf; i+=s->facets[i]+1,k++) {
      /* make sure it is at least a tri (i.e. not a line) */
      if (s->facets[i] >= 3) {
        /* loop over nodes of face */
        for (j=0; j<s->facets[i]; j++) {
          cv[iCV+j] = nOff[myid] + iN + s->facets[i+j+1];
        } /* end loop over nodes */

        /* store the number of vertices */
        nVert[iF] = s->facets[i];

        /* store cell scalars here */
        c = RP_CELL(s->cells+k);
        t = RP_THREAD(s->cells+k);
        fVel[3*iF+0] = (float)(C_U(c,t));
        fVel[3*iF+1] = (float)(C_V(c,t));
        fVel[3*iF+2] = (float)(C_W(c,t));

        iCV += s->facets[i];
        iF++;
      }
    } /* end loop over faces*/
    iN += s->np;
  } /* end loop over zone threads */

#if PARALLEL /* reduce everything to root process */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++) {
      PRF_CRECV_INT(i, &cv[cvOff[i]],   nCVs[i],   myid);
      PRF_CRECV_INT(i, &nVert[fOff[i]], nFaces[i], myid);
      PRF_CRECV_FLOAT(i, &fVel[fxOff[i]], 3*nFaces[i], myid);
    }
  }
  else {
    PRF_CSEND_INT(node_zero, cv,    nCVs[myid],   node_zero);
    PRF_CSEND_INT(node_zero, nVert, nFaces[myid], node_zero);
    PRF_CSEND_FLOAT(node_zero, fVel, 3*nFaces[myid], node_zero);
  }
#endif /* PARALLEL */
    
  /* write <Cells> */
  if (myid == node_zero) {
    offs  = (int*)  calloc(totFaces, sizeof(int));

    /* compute the offsets and get cell types */
    offs[0]  = nVert[0];
    for (j=1; j<totFaces; j++) {
      offs[j]  = offs[j-1] + nVert[j];
    }
    
    /* cell connectivity */
    bytes = sizeof(int) + totCV*sizeof(int);
    writeDataInt(dF, 1, &bytes, 1);
    writeDataInt(dF, totCV, cv, totCV); 
    
    /* cell offsets */
    bytes = sizeof(int) + totFaces*sizeof(int);
    writeDataInt(dF, 1, &bytes, 1);
    writeDataInt(dF, totFaces, offs, totFaces);

    free(offs);

    /* cell data (velocity) */
    bytes = sizeof(int) + 3*totFaces*sizeof(float);
    writeDataInt(dF, 1, &bytes, 1);
    writeDataInt(dF, 3*totFaces, fVel, 1);
  } /* end <Cells> write */
  free(fVel);

  /* write all the scalars
     NOTE: this should probably not be in a loop like this, but ... */
  float *fSc  = (float*) calloc(nF,  sizeof(float));
  for (j=0; j<ex.nSc; j++) {
    iF = 0;
    for (iT=0; iT<nZt; iT++) {
      sid = Zt[iT];
      if (sid < 0) break;
      s = SurfaceList + sid;
      
      /* loop over the facets of this surface */
      for (k=0, i=0; k<s->nf; i+=s->facets[i]+1,k++) {
        /* make sure it is at least a tri (i.e. not a line) */
        if (s->facets[i] >= 3) {
          /* store cell scalars here */
          c = RP_CELL(s->cells+k);
          t = RP_THREAD(s->cells+k);
          
          if (ex.sc_type[j] == 0)
            fSc[iF] = (float)(ABS_P(C_P(c,t),op_pres));        /* pressure */
          else if (ex.sc_type[j] == 1)
            fSc[iF] = (float)(C_STORAGE_R(c,t,ex.sc_indx[j])); /* standard */
          else if (ex.sc_type[j] == 2)
            fSc[iF] = (float)(C_UDMI(c,t,ex.sc_indx[j]));      /* udm */
          else if (ex.sc_type[j] == 3)
            fSc[iF] = (float)(C_YI(c,t,ex.sc_indx[j]));        /* species */
          
          iF++;
        }
      } /* end loop over faces*/
    } /* end loop over zone threads */

#if PARALLEL /* reduce everything to root process */
    if (myid == node_zero) {
      for (i=1; i<compute_node_count; i++)
        PRF_CRECV_FLOAT(i, &fSc[fOff[i]], nFaces[i], myid);
    }
    else {
      PRF_CSEND_FLOAT(node_zero, fSc, nFaces[myid], node_zero);
    }
#endif /* PARALLEL */

    if (myid == node_zero) {
      bytes = sizeof(int) + totFaces*sizeof(float);
      writeDataInt(dF, 1, &bytes, 1);
      writeDataInt(dF, totFaces, fSc, totFaces);
    }
  }

  free(nodes); free(nodeID); free(nodeIndx);
  free(cv); free(nVert); free(fSc);

  /* ========== end the appended data ========== */
  if (myid == node_zero) {
    fprintf(dF, "\n%2s</AppendedData>\n", "");
    fprintf(dF, "</VTKFile>\n");
    fclose(dF);
  }

  free(nFaces); free(nNodes); free(nCVs); free(iWork);
  free(fOff); free(fxOff); free(nOff); free(nxOff); free(cvOff);
  free(Zt);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void VTKwriteDummyData(int n, FILE *pFile)
{
#if !RP_HOST
  Message0("writing dummy data of length: %i\n", n);
  int bytes;
  float *dummy = (float*) calloc(n, sizeof(float));
  bytes = sizeof(int) + n*sizeof(float);

  writeDataInt(pFile, 1, &bytes, 1);
  writeDataReal(pFile, n, dummy, 1);

  free(dummy);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

int VTKcellType(int verts)
{
  switch (verts) {   /* 3D cell types */
  case TET_CELL:     /* VTK_TETRA */
    return 10; 
  case PYRAMID_CELL: /* VTK_PYRAMID */
    return 14;
  case WEDGE_CELL:   /* VTK_WEDGE */
    return 13;
  case HEX_CELL:     /* VTK_HEXAHEDRON */
    return 12;
  case POLY_CELL:    /* VTK_POLYHEDRON */
    return 42;
  } /* end switch */
}
/* -------------------------------------------------------------------------- */

void VTKinfoDataSet(FILE *p,int tab,double t,char *group,char *part,char *file)
{
  char *tmp = (char*) calloc(1024, sizeof(char));
  VTKstrDataSet(tmp,tab,t,group,part,file);
  fprintf(p, "%s", tmp);
  free(tmp);
}
/* -------------------------------------------------------------------------- */

void VTKstrDataSet(char *p,int tab,double t,char *group,char *part,char *file)
{
  char *tmp = (char*) calloc(128, sizeof(char));
  char *out = (char*) calloc(1024, sizeof(char));
  sprintf(out, "%*s<DataSet ", tab, ""); 
  sprintf(tmp, "timestep=\"%.6e\" ", t);   strcat(out, tmp);
  sprintf(tmp, "part=\"%s\" ", part);      strcat(out, tmp);
  sprintf(tmp, "group=\"%s\" ", group);    strcat(out, tmp);
  sprintf(tmp, "file=\"%s\" ", file);      strcat(out, tmp);
  sprintf(tmp, "/>\n");                    strcat(out, tmp);
  strcpy(p, out);
  free(tmp); free(out);
}
/* -------------------------------------------------------------------------- */

void VTKinfoDataArray(FILE *p,int tab,char *type,char *name,int n,char *format,
                      long long int offset)
{
  char *tmp = (char*) calloc(1024, sizeof(char));
  VTKstrDataArray(tmp,tab,type,name,n,format,offset);
  fprintf(p, "%s", tmp);
  free(tmp);
}
/* -------------------------------------------------------------------------- */

void VTKstrDataArray(char *p,int tab,char *type,char *name,int n,char *format, 
                     long long int offset)
{
  char *tmp = (char*) calloc(128, sizeof(char));
  char *out = (char*) calloc(1024, sizeof(char));
  sprintf(out, "%*s<DataArray ", tab, "");
  sprintf(tmp, "type=\"%s\" ", type);              strcat(out, tmp);
  sprintf(tmp, "Name=\"%s\" ", name);              strcat(out, tmp);
  if (n > 1) {
    sprintf(tmp, "NumberOfComponents=\"%i\" ", n); strcat(out, tmp);
  }
  sprintf(tmp, "format=\"%s\" ", format);          strcat(out, tmp);
  if (!strcmp(format, "appended")) {
    sprintf(tmp, "offset=\"%i\" ", offset);        strcat(out, tmp);
    sprintf(tmp, "/>\n");                          strcat(out, tmp);
  }
  else {
    sprintf(tmp, ">\n");                           strcat(out, tmp);
  }
  strcpy(p, out);
  free(tmp); free(out);
}
/* -------------------------------------------------------------------------- */

int XMLparseLine(char *xmlIn, char **type, char ***tags, char ***vals)
{
  int i, j, nRaw, nArgs;
  char *pch;

  /* Message0("parsing xml line: %s\n", xmlIn); */

  type[0] = (char*)  calloc(128, sizeof(char));

  char **raw = (char**) calloc(1, sizeof(char*));
  char seps[2] = " =";

  nRaw = 0;
  pch = strtok(xmlIn, seps);
  while (pch != NULL) {
    nRaw++;
    raw         = (char**) realloc(raw, nRaw*sizeof(char*));
    raw[nRaw-1] = (char*)   calloc(128, sizeof(char));
    strcpy(raw[nRaw-1], pch);
    pch = strtok(NULL, seps);
  }

  /* copy is in the first one, take away the < */
  strcpy(type[0], &raw[0][1]);
  
  /* make sure spaces are there so that the first and last are not arguments */
  nArgs = (nRaw-2)/2;
  tags[0] = (char**) calloc(nArgs,   sizeof(char*));
  vals[0] = (char**) calloc(nArgs,   sizeof(char*));
  for (i=0, j=1; i<nArgs; i++, j+=2) {
    tags[0][i] = (char*) calloc(256, sizeof(char));
    vals[0][i] = (char*) calloc(256, sizeof(char));
    /* Message0("i: %i, j: %i\n", i, j); */
    strcpy(tags[0][i], raw[j]);
    strncpy(vals[0][i], &raw[j+1][1], strlen(raw[j+1])-2);
  }
  
  for (i=0; i<nRaw; i++) free(raw[i]);
  free(raw);

  return nArgs;
}
/* -------------------------------------------------------------------------- */
