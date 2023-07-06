#include "ensightIO.h"

extern expStruct ex;

void UDFexportEnsight(char *dir, char *prefix)
{
  int i;
#if PARALLEL && RP_NODE
  if (myid == node_zero)
#endif /* PARALLEL && RP_NODE */
    UDFwriteEnsightCase(dir, prefix);

  /* write the geometry for this time-step */
  UDFwriteEnsightGeom(dir, prefix);

  /* write the velocity */
  if (ex.dump_velocity)
    UDFwriteEnsightVelocity(dir, prefix);

  /* write all the scalars */
  for (i=0; i<ex.nSc; i++)
    UDFwriteEnsightScalar(dir, prefix, ex.sc_names[i], ex.sc_abbrv[i],
                          ex.sc_indx[i], ex.sc_type[i]);

  /* write the particle data */
  if (ex.dump_particle)
    UDFwriteEnsightParticle(dir, prefix, ex.nInj*ex.nNoz, ex.all_streams);
}
/* -------------------------------------------------------------------------- */

void UDFwriteEnsightCase(char *dir, char *prefix)
{
#if !RP_HOST

  int  i, nt, iStart, increment, line, perLine = 5;
  long double r;
  double *tVals, cur_time;
  char fileName[1000], buffer[1000], tag[50], *ident;
  sprintf(fileName, "%s/%s.case", dir, prefix);
  Message0("  %-46s  %s\n", "Writing case file:", strrchr(fileName, '/')+1);

  FILE *pFile = fopen(fileName, "r");

  if (in_cylinder_p)
    cur_time = RP_Get_Real("dynamesh/in-cyn/crank-angle");
  else
    cur_time = CURRENT_TIME;

  if (pFile && ex.init == 1) { /* file exists, read in contents to get time */
    /* read in until the time section is reached */
    while (!feof(pFile)) { 
      fscanf(pFile, "%s", tag);  /* read in string */
      if (!strcmp(tag, "TIME"))
        break;
    }
    
    fgets(buffer, 80, pFile); /* time set line */
    fgets(buffer, 80, pFile); /* time set */
    fgets(buffer, 80, pFile); /* number of steps */
    ident = strtok(buffer, ":");
    nt = atoi(strtok(NULL,""));
    fgets(buffer, 80, pFile); /* file start number */
    ident = strtok(buffer, ":");
    iStart = atoi(strtok(NULL,""));
    fgets(buffer, 80, pFile); /* file increment */
    ident = strtok(buffer, ":");
    increment = atoi(strtok(NULL,""));
    fgets(buffer, 80, pFile); /* time values: (assumes on single line) */
    
    tVals = (double*) calloc(nt+1, sizeof(double));
    for (i=0; i<nt; i++) {
      fscanf(pFile, "%Lf\n", &r);
      tVals[i] = (double)(r);
      if (tVals[i] > cur_time) {
        nt = i-1;
        break;
      }
    }
    nt++;
    tVals[nt-1] = cur_time;
    ex.fileNo = nt;
    
    fclose(pFile);
  }
  else {
    nt        = 1;
    iStart    = 1;
    ex.fileNo = iStart;
    tVals    = (double*) calloc(nt, sizeof(double));
    tVals[0] = cur_time;
    ex.init = 1; /* now restart */
  }
  
  /* open the file for writing */
  pFile = fopen(fileName, "w");

  fprintf(pFile, "FORMAT\n");
  fprintf(pFile, "type:\tensight gold\n");

  /* write out the geometry */
  fprintf(pFile, "GEOMETRY\n");
  if (ex.single_dir)
    fprintf(pFile, "model: 1 data/%s-geo.*****\n", prefix);
  else
    fprintf(pFile, "model: 1 geometry/%s-geo.*****\n", prefix);
  if (ex.dump_particle)
    fprintf(pFile, "measured: 1 particle/%s-mpg.*****\n", prefix);

  /* write out all the variables */
  fprintf(pFile, "VARIABLE\n");
  if (ex.dump_velocity) {
    if (ex.single_dir)
      fprintf(pFile, "vector per element: 1 Velocity data/%s-vel.*****\n",
              prefix);
    else
      fprintf(pFile, "vector per element: 1 Velocity velocity/%s-vel.*****\n",
              prefix);
  }

  /* scalars */
  for (i=0; i<ex.nSc; i++) {
    if (ex.single_dir)
      fprintf(pFile, "scalar per element: 1 %s data/%s-%s.*****\n",
              ex.sc_names[i], prefix, ex.sc_abbrv[i]);
    else
      fprintf(pFile, "scalar per element: 1 %s %s/%s-%s.*****\n",
              ex.sc_names[i], ex.sc_abbrv[i], prefix, ex.sc_abbrv[i]);
  }

  /* particle variables */ 
  if (ex.dump_particle) {
    fprintf(pFile, "scalar per measured node: 1 Drop_velocity particle/%s-vel.*****\n", prefix);
    fprintf(pFile, "scalar per measured node: 1 Drop_temperature particle/%s-temp.*****\n", prefix);
    fprintf(pFile, "scalar per measured node: 1 Drop_diameter particle/%s-diam.*****\n", prefix);
  }

  /* write the time */
  fprintf(pFile, "TIME\n");
  fprintf(pFile, "time set: 1\n");
  fprintf(pFile, "number of steps: %i\n", nt);
  fprintf(pFile, "filename start number: %i\n", iStart);
  fprintf(pFile, "filename increment: 1\n"); /* using individual files */
  fprintf(pFile, "time values:\n");
  for (i=0, line=1; i<nt; i++) {
    fprintf(pFile, "%15.6e", tVals[i]);
    if (line%perLine == 0)
      fprintf(pFile, "\n");
    line++;
  }
  fclose(pFile);
  free(tVals);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

#if PARALLEL
void UDFwriteEnsightGeom_SP(char *dir, char *prefix, char *ident)
{
#if !RP_HOST
  Domain *d;
  Thread *t;
  cell_t c;
  face_t f;
  Node   *node;
  double x[3];

  d = Get_Domain(1); /* get the domain from FLUENT */

  int i, j, n, iC, iN, *nCells, *nNodes, *iWork, iTmp = 1;
  int nL, nC, *nOff,  *nodeMask, *nodeID, *nodeIndx, nodeDupl;
  float *nodeX, *nodeY, *nodeZ;
  int flag;
  int totNodes, totCells;
  char buffer[80];
  float fTmp, *fTmpA;

  char fileName[1000];
  if (ex.single_dir)
    sprintf(fileName, "%s/data/%s-geo.%05i", dir, prefix, ex.fileNo);
  else
    sprintf(fileName, "%s/geometry/%s-geo.%05i", dir, prefix, ex.fileNo);
  Message0("  %-46s  %s\n", "Writing geometry file:", strrchr(fileName, '/')+1);

#if UDF_MPI_WRITE
  MPI_Status status;
  MPI_File *pFile;
  flag = MPI_File_open(MPI_COMM_SELF, fileName, MPI_MODE_WRONLY |
                       MPI_MODE_CREATE,MPI_INFO_NULL,&pFile);
  MPI_Offset *cOff=(MPI_Offset*) calloc(compute_node_count, sizeof(MPI_Offset));
#else /* ROOT WRITE */
  FILE *pFile = fopen(fileName, "wb");
  int *cOff = (int*) calloc(compute_node_count, sizeof(int));
#endif /* UDF_MPI_WRITE */

  nNodes = (int*) calloc(compute_node_count, sizeof(int));
  nCells = (int*) calloc(compute_node_count, sizeof(int));
  nOff   = (int*) calloc(compute_node_count, sizeof(int));
  iWork  = (int*) calloc(compute_node_count, sizeof(int));

  /* get the total number of nodes and cells internal to the process */
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      c_node_loop(c,t,n) {
        node = C_NODE(c,t,n);
        NODE_ITMP1(node) = 0;
      }
    } end_c_loop_int(c,t);
  }
  nNodes[myid] = nCells[myid] = 0;
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      c_node_loop(c,t,n) {
        node = C_NODE(c,t,n);
        if (!NODE_ITMP1(node)) {
          NODE_ITMP1(node) = 1;
          nNodes[myid]++;
        }
      }
      nCells[myid]++;
    } end_c_loop_int(c,t);
  }

  /* accumlate to all processes */
  PRF_GISUM(nNodes, compute_node_count, iWork);
  for (i=0, totNodes=0; i<compute_node_count; i++)
    totNodes += nNodes[i];

  PRF_GISUM(nCells, compute_node_count, iWork);
  for (i=0, totCells=0; i<compute_node_count; i++)
    totCells += nCells[i];

  /* allocate enough memory on the root for all nodes, otherwise only need
     for local */
  if (myid == node_zero)
    nL = totNodes;
  else
    nL = nNodes[myid];

  /* allocate the necessary storage for the nodes */
  nodeID   = (int*)   calloc(nL, sizeof(int));
  nodeIndx = (int*)   calloc(nL, sizeof(int));
  nodeMask = (int*)   calloc(nL, sizeof(int));
  nodeX    = (float*) calloc(nL, sizeof(float));
  nodeY    = (float*) calloc(nL, sizeof(float));
  nodeZ    = (float*) calloc(nL, sizeof(float));

  /* find the offset for each node when transferring to the root */
  for (i=1, nOff[0]=0; i<compute_node_count; i++)
    for (j=0, nOff[i]=0; j<i; j++)
      nOff[i] += nNodes[j];

  /* fill all the node data */
  iN = 0;
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      c_node_loop(c,t,n) {
        node = C_NODE(c,t,n);
        if (NODE_ITMP1(node)) {
          nodeID[iN] = NODE_ID(node);
          nodeX[iN]  = NODE_X(node);
          nodeY[iN]  = NODE_Y(node);
          nodeZ[iN]  = NODE_Z(node);
          NODE_ITMP1(node) = 0;
          iN++;
        }
      }
    } end_c_loop_int(c,t);
  }

  /* transfer all the node data to the root */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++) {
      PRF_CRECV_INT(i, &nodeID[nOff[i]], nNodes[i], myid);
      PRF_CRECV_FLOAT(i, &nodeX[nOff[i]], nNodes[i], myid);
      PRF_CRECV_FLOAT(i, &nodeY[nOff[i]], nNodes[i], myid);
      PRF_CRECV_FLOAT(i, &nodeZ[nOff[i]], nNodes[i], myid);
    }
  }
  else {
    PRF_CSEND_INT(node_zero, nodeID, nNodes[myid], node_zero);
    PRF_CSEND_FLOAT(node_zero, nodeX, nNodes[myid], node_zero);
    PRF_CSEND_FLOAT(node_zero, nodeY, nNodes[myid], node_zero);
    PRF_CSEND_FLOAT(node_zero, nodeZ, nNodes[myid], node_zero);
  }
  
  int increasing = 1;
  if (myid == node_zero) {
    /* check to see if the nodes are already in increasing order */
    for (i=0, increasing = 1; i<totNodes; i++)
      if (nodeID[nodeIndx[i]] < nodeID[nodeIndx[i-1]])
        increasing = 0;
    
    /* index array in sorted order of node ident if necessary */
    indexIntArray(totNodes, nodeID, nodeIndx);
    
    /* find and mask duplicate nodes */
    for (i=1, nodeDupl=0; i<totNodes; i++) {
      nodeMask[nodeIndx[i]] = 0;
      if (nodeID[nodeIndx[i]] == nodeID[nodeIndx[i-1]]) {
        nodeMask[nodeIndx[i]] = 1;
        nodeDupl++;
      }
    }
    
    /* update the total nodes on all processors */
    totNodes = totNodes - nodeDupl;
    for (i=1; i<compute_node_count; i++)
      PRF_CSEND_INT(i, &totNodes, 1, myid);
    
    /* write the header information */
    writeDataStr(pFile, 80, "C Binary");
    writeDataStr(pFile, 80, "Model Geometry File exported from Fluent");
    writeDataStr(pFile, 80, "Ensight-Gold format");
    writeDataStr(pFile, 80, "node id assign");
    writeDataStr(pFile, 80, "element id assign");
    writeDataStr(pFile, 80, "part");
    iTmp = 1;
    writeDataInt(pFile, 1, &iTmp, 1);
    sprintf(buffer, "%s_volume", ident);
    writeDataStr(pFile, 80, buffer);
    writeDataStr(pFile, 80, "coordinates");
    writeDataInt(pFile, 1, &totNodes, 1);
    
    /* write the node coordinates */
    fTmpA = (float*) calloc(totNodes, sizeof(float));
    for (i=0, j=0; i<totNodes+nodeDupl; i++) {
      if (nodeMask[nodeIndx[i]] == 0) {
        fTmpA[j] = nodeX[nodeIndx[i]];
        j++;
      }
    }
    writeDataReal(pFile, totNodes, fTmpA, 1);
    for (i=0, j=0; i<totNodes+nodeDupl; i++) {
      if (nodeMask[nodeIndx[i]] == 0) {
        fTmpA[j] = nodeY[nodeIndx[i]];
        j++;
      }
    }
    writeDataReal(pFile, totNodes, fTmpA, 1);
    for (i=0, j=0; i<totNodes+nodeDupl; i++) {
      if (nodeMask[nodeIndx[i]] == 0) {
        fTmpA[j] = nodeZ[nodeIndx[i]];
        j++;
      }
    }
    writeDataReal(pFile, totNodes, fTmpA, 1);
    free(fTmpA);
    
    /* write the element information */
    writeDataStr(pFile, 80, "hexa8");
    writeDataInt(pFile, 1, &totCells, 1);
  }
  else {
    PRF_CRECV_INT(node_zero, &totNodes, 1, node_zero);
  }
  free(nodeID); free(nodeMask); free(nodeIndx);
  free(nodeX); free(nodeY); free(nodeZ);

  nC = 8*nCells[myid];
  /* find the offset for each node when transferring to the root */
  for (i=1, cOff[0]=0; i<compute_node_count; i++)
    for (j=0, cOff[i]=0; j<i; j++)
      cOff[i] += 8*nCells[j];

#if UDF_MPI_WRITE
  /* change offset to bytes */
  cOff[myid] *= sizeof(int);

  /* advance the front material */
  cOff[myid] += 8*80*sizeof(char) + 2*sizeof(int);     /* header */
  cOff[myid] += 3*totNodes*sizeof(float);              /* coordinates */
  cOff[myid] += 80*sizeof(char) + 1*sizeof(int);       /* element type */

  MPI_File_set_view(pFile,cOff[myid],MPI_INT,MPI_INT,"native",MPI_INFO_NULL);
#else /* ROOT WRITE */
  if (myid == node_zero)
    nC = 8*totCells;
#endif /* UDF_MPI_WRITE */
  /* allocate an array for the connectivities on this processor */
  int *conn = (int*) calloc(nC, sizeof(int));

  /* write the connectivity */
  iC = 0;
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      iN = 0;
      c_node_loop(c,t,n) {
        node = C_NODE(c,t,n);
        conn[8*iC+iN] = NODE_ID(node);
        iN++;
      }
      iC++;
    } end_c_loop_int(c,t);
  }

#if UDF_MPI_WRITE
  MPI_File_write(pFile, conn, 8*nCells[myid], MPI_INT, &status);
  MPI_File_close(&pFile);
#else /* ROOT WRITE */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++)
      PRF_CRECV_INT(i, &conn[cOff[i]], 8*nCells[i], myid);
    
    writeDataInt(pFile, 8*totCells, conn, 8);
    fclose(pFile);
  }
  else {
    PRF_CSEND_INT(node_zero, conn, 8*nCells[myid], node_zero);
  }
#endif /* UDF_MPI_WRITE */
  free(conn);
  free(nNodes); free(nCells); free(nOff); free(cOff); free(iWork);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */
#else /* SERIAL */
void UDFwriteEnsightGeom_SP(char *dir, char *prefix, char *ident)
{
#if !RP_HOST
  Domain *d;
  Thread *t;
  cell_t c;
  face_t f;
  Node   *node;

  d = Get_Domain(1); /* get the domain from FLUENT */

  int i, n, iC, iN, nCells, nNodes, iTmp = 1, *iTmpA;
  float fTmp, *fTmpA;

  /* generate the filename */
  char fileName[1000], buffer[80];
  if (ex.single_dir)
    sprintf(fileName, "%s/data/%s-geo.%05i", dir, prefix, ex.fileNo);
  else
    sprintf(fileName, "%s/geometry/%s-geo.%05i", dir, prefix, ex.fileNo);
  Message0("  %-46s  %s\n", "Writing geometry file:", strrchr(fileName, '/')+1);

  FILE *pFile = fopen(fileName, "wb");

  /* write the header information */
  writeDataStr(pFile, 80, "C Binary");
  writeDataStr(pFile, 80, "Model Geometry File exported from Fluent");
  writeDataStr(pFile, 80, "Ensight-Gold format");
  writeDataStr(pFile, 80, "node id assign");
  writeDataStr(pFile, 80, "element id assign");
  writeDataStr(pFile, 80, "part");
  iTmp = 1;
  writeDataInt(pFile, 1, &iTmp, 1);
  sprintf(buffer, "%s_volume", ident);
  writeDataStr(pFile, 80, buffer);
  writeDataStr(pFile, 80, "coordinates");
  nNodes = d->node_count;
  writeDataInt(pFile, 1, &nNodes, 1);

  /* output the node co-ordinates */
  fTmpA = (float*) calloc(nNodes, sizeof(float));
  /* x co-ordinates */
  iN = 0;
#if RampantReleaseMajor >= 15
  begin_node_loop(node,d) {
#else
  begin_node_loop(node) {
#endif /* RampantReleaseMajor >= 15 */
    fTmpA[iN] = NODE_X(node); iN++;
    NODE_ITMP1(node) = iN;
  } end_node_loop(node)
  writeDataReal(pFile, nNodes, fTmpA, 1);
  /* y co-ordinates */
  iN = 0;
#if RampantReleaseMajor >= 15
  begin_node_loop(node,d) { 
#else
  begin_node_loop(node) { 
#endif /* RampantReleaseMajor >= 15 */
    fTmpA[iN] = NODE_Y(node); iN++;
  } end_node_loop(node)
  writeDataReal(pFile, nNodes, fTmpA, 1);
  /* Z co-ordinates */
  iN = 0;
#if RampantReleaseMajor >= 15
  begin_node_loop(node,d) {
#else
  begin_node_loop(node) { 
#endif /* RampantReleaseMajor >= 15 */
    fTmpA[iN] = NODE_Z(node); iN++;
  } end_node_loop(node)
  writeDataReal(pFile, nNodes, fTmpA, 1);
  free(fTmpA);

  /* write the cell connectivity */
  nCells=local_count_c(d);
  writeDataStr(pFile, 80, "hexa8");
  writeDataInt(pFile, 1, &nCells, 1);

  iTmpA = (int*) calloc(8*nCells, sizeof(int));
  iC = iN = 0;
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      iN = 0;
      c_node_loop(c,t,n) {
        node = C_NODE(c,t,n);
        iTmpA[8*iC+iN] = NODE_ITMP1(node);
        iN++;
      }
      iC++;
    } end_c_loop_int(c,t);
  }
  writeDataInt(pFile, 8*nCells, iTmpA, 8);
  free(iTmpA);
  fclose(pFile);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */
#endif /* PARALLEL */

void UDFwriteEnsightVelocity_SP(char *dir, char *prefix)
{
#if !RP_HOST
  Domain *d;
  Thread *t;
  cell_t c;
  int i, j, id, nC, nCells, *nCellsL, totCells, *iWork;
  float *cU, *cV, *cW;

  d = Get_Domain(1); /* get the domain from FLUENT */

  char fileName[1000];
  if (ex.single_dir)
    sprintf(fileName, "%s/data/%s-vel.%05i", dir, prefix, ex.fileNo);
  else
    sprintf(fileName, "%s/velocity/%s-vel.%05i", dir, prefix, ex.fileNo);
  Message0("  %-46s  %s\n", "Writing velocity file:", strrchr(fileName, '/')+1);

#if PARALLEL
#if UDF_MPI_WRITE
  int flag;
  MPI_Status status;
  MPI_File *pFile;
  flag = MPI_File_open(MPI_COMM_SELF, fileName, MPI_MODE_WRONLY |
                       MPI_MODE_CREATE,MPI_INFO_NULL,&pFile);
  MPI_Offset *cOff=(MPI_Offset*) calloc(compute_node_count, sizeof(MPI_Offset));
#else /* ROOT WRITE */
  FILE *pFile;
  if (myid == node_zero)
    pFile = fopen(fileName, "wb");
  int *cOff = (int*) calloc(compute_node_count, sizeof(int));
#endif /* UDF_MPI_WRITE */
#else /* SERIAL */
  FILE *pFile = fopen(fileName, "wb");
#endif /* PARALLEL */

  if (pFile == NULL)
    Message0("file not opened\n");

#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      /* write the header */
      writeDataStr(pFile, 80, "Velocity");
      writeDataStr(pFile, 80, "part");      /* beginning of part */
      id = 1;
      writeDataInt(pFile, 1, &id, 1);
      writeDataStr(pFile, 80, "hexa8");      /* element type */
    }

  /* count the number of cells */
  nCells = 0;
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      nCells++;
    } end_c_loop_int(c,t);
  }

  nC = nCells;
#if PARALLEL
  nCellsL = (int*) calloc(compute_node_count, sizeof(int));
  iWork   = (int*) calloc(compute_node_count, sizeof(int));
  nCellsL[myid] = nCells;
  totCells = PRF_GRSUM1(nCells);
  PRF_GISUM(nCellsL, compute_node_count, iWork);

  /* compute the offset */
  for (i=1, cOff[0]=0; i<compute_node_count; i++)
    for (j=0, cOff[i]=0; j<i; j++)
      cOff[i] += nCellsL[j];

#if UDF_MPI_WRITE
  cOff[myid] *= sizeof(float);

  cOff[myid] += 3*80*sizeof(char) + 1*sizeof(int);
#else
  if (myid == node_zero)
    nC = totCells;
#endif /* UDF_MPI_WRITE */
#endif /* PARALLEL */

  /* store the velocities in an array */
  cU = (float*) calloc(nC, sizeof(float));
  cV = (float*) calloc(nC, sizeof(float));
  cW = (float*) calloc(nC, sizeof(float));

  i = 0;
  thread_loop_c(t,d) {
    begin_c_loop_int(c,t) {
      cU[i] = (float)(C_U(c,t));
      cV[i] = (float)(C_V(c,t));
      cW[i] = (float)(C_W(c,t));
      i++;
    } end_c_loop_int(c,t);
  }
  
#if PARALLEL
#if UDF_MPI_WRITE
  /* write the velocities */
  MPI_File_set_view(pFile,cOff[myid],MPI_INT,MPI_INT,"native",MPI_INFO_NULL);
  MPI_File_write(pFile, cU, nCells, MPI_INT, &status);

  cOff[myid] += totCells*sizeof(float);
  MPI_File_set_view(pFile,cOff[myid],MPI_INT,MPI_INT,"native",MPI_INFO_NULL);
  MPI_File_write(pFile, cV, nCells, MPI_INT, &status);

  cOff[myid] += totCells*sizeof(float);
  MPI_File_set_view(pFile,cOff[myid],MPI_INT,MPI_INT,"native",MPI_INFO_NULL);
  MPI_File_write(pFile, cW, nCells, MPI_INT, &status);
  
  MPI_File_close(&pFile);
#else /* ROOT WRITE */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++) {
      PRF_CRECV_FLOAT(i, &cU[cOff[i]], nCellsL[i], myid);
      PRF_CRECV_FLOAT(i, &cV[cOff[i]], nCellsL[i], myid);
      PRF_CRECV_FLOAT(i, &cW[cOff[i]], nCellsL[i], myid);
    }
    
    writeDataReal(pFile, totCells, cU, 1);
    writeDataReal(pFile, totCells, cV, 1);
    writeDataReal(pFile, totCells, cW, 1);
    fclose(pFile);
  }
  else {
    PRF_CSEND_FLOAT(node_zero, cU, nCellsL[myid], node_zero);
    PRF_CSEND_FLOAT(node_zero, cV, nCellsL[myid], node_zero);
    PRF_CSEND_FLOAT(node_zero, cW, nCellsL[myid], node_zero);
  }
#endif /* UDF_MPI_WRITE */
  free(nCellsL); free(iWork); free(cOff);
#else /* SERIAL */
  fwrite(cU, sizeof(float), nCells, pFile);
  fwrite(cV, sizeof(float), nCells, pFile);
  fwrite(cW, sizeof(float), nCells, pFile);
  fclose(pFile);
#endif /* PARALLEL */
  free(cU); free(cV); free(cW);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void UDFwriteEnsightVelocity(char *dir, char *prefix)
{
#if !RP_HOST
  int i, j, id;

  char fileName[1000];
  if (ex.single_dir)
    sprintf(fileName, "%s/data/%s-vel.%05i", dir, prefix, ex.fileNo);
  else
    sprintf(fileName, "%s/velocity/%s-vel.%05i", dir, prefix, ex.fileNo);
  Message0("  %-46s  %s\n", "Writing velocity file:", strrchr(fileName, '/')+1);

#if PARALLEL
#if UDF_MPI_WRITE
  int flag;
  MPI_Status status;
  MPI_File *pFile;
  flag = MPI_File_open(MPI_COMM_SELF, fileName, MPI_MODE_WRONLY |
                       MPI_MODE_CREATE,MPI_INFO_NULL,&pFile);
#else /* ROOT WRITE */
  FILE *pFile;
  if (myid == node_zero)
    pFile = fopen(fileName, "wb");
#endif /* UDF_MPI_WRITE */
#else /* SERIAL */
  FILE *pFile = fopen(fileName, "wb");
#endif /* PARALLEL */

  if (pFile == NULL)
    Message0("file not opened\n");

#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    { /* write the header */
      writeDataStr(pFile, 80, "Velocity");
    }

  /* write the cell geometry */
  for (i=0; i<ex.nCells; i++) {
#if PARALLEL
    if (myid == node_zero)
#endif /* PARALLEL */
      { /* write the part header */
        writeDataStr(pFile, 80, "part");      /* beginning of part */
        id = i+1;
        writeDataInt(pFile, 1, &id, 1);
      }
    
    /* write the cell connectivity for the current thread */
    switch (THREAD_ELEMENT_TYPE(ex.cell_th[i])) {
    case MIXED_CELL:
      for (j=0; j<ex.part_nv[i]; j++)
        writeEnsightVector(ex.cell_th[i], pFile, ex.part_verts[i][j],
                           ex.part_type[i][j]);
      break;
    default:
      writeEnsightVector(ex.cell_th[i], pFile, ex.part_verts[i][0],
                         ex.part_type[i][0]);
      break;
    } /* end switch */
  }

#if PARALLEL
#if UDF_MPI_WRITE
  MPI_File_close(&pFile);
#else /* ROOT WRITE */
  if (myid == node_zero) fclose(pFile);
#endif /* UDF_MPI_WRITE */
#else /* SERIAL */
  fclose(pFile);
#endif
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void writeEnsightVector(Thread *t, FILE *pFile, int verts, char *cell)
{
#if !RP_HOST
  cell_t c;
  int i, j, nC, nCells, *nCellsL, totCells, *iWork;
  float *cU, *cV, *cW;

#if PARALLEL
#if UDF_MPI_WRITE
  int flag;
  MPI_Status status;
  MPI_Offset *cOff=(MPI_Offset*) calloc(compute_node_count, sizeof(MPI_Offset));
#else /* ROOT WRITE */
  int *cOff = (int*) calloc(compute_node_count, sizeof(int));
#endif /* UDF_MPI_WRITE */
  if (myid == node_zero)
#endif /* PARALLEL */
    { /* write the element type */
      writeDataStr(pFile, 80, cell);
    }

  /* count the number of cells */
  nCells = 0;
  begin_c_loop_int(c,t) {
    if (C_NNODES(c,t) == verts)
      nCells++;
  } end_c_loop_int(c,t);

  nC = nCells;
#if PARALLEL
  nCellsL = (int*) calloc(compute_node_count, sizeof(int));
  iWork   = (int*) calloc(compute_node_count, sizeof(int));
  nCellsL[myid] = nCells;
  totCells = PRF_GRSUM1(nCells);
  PRF_GISUM(nCellsL, compute_node_count, iWork);

  /* compute the offset */
  for (i=1, cOff[0]=0; i<compute_node_count; i++)
    for (j=0, cOff[i]=0; j<i; j++)
      cOff[i] += nCellsL[j];

#if UDF_MPI_WRITE
  cOff[myid] *= sizeof(float); /* convert to bytes */
  cOff[myid] += 3*80*sizeof(char) + 1*sizeof(int); /* add header */
#else
  if (myid == node_zero)
    nC = totCells;
#endif /* UDF_MPI_WRITE */
#endif /* PARALLEL */

  /* store the velocities in an array */
  cU = (float*) calloc(nC, sizeof(float));
  cV = (float*) calloc(nC, sizeof(float));
  cW = (float*) calloc(nC, sizeof(float));

  i = 0;
  begin_c_loop_int(c,t) {
    if (C_NNODES(c,t) == verts) {
      cU[i] = (float)(C_U(c,t));
      cV[i] = (float)(C_V(c,t));
      cW[i] = (float)(C_W(c,t));
      i++;
    }
  } end_c_loop_int(c,t);

#if PARALLEL
#if UDF_MPI_WRITE
  /* write the velocities */
  MPI_File_set_view(pFile,cOff[myid],MPI_INT,MPI_INT,"native",MPI_INFO_NULL);
  MPI_File_write(pFile, cU, nCells, MPI_INT, &status);

  cOff[myid] += totCells*sizeof(float);
  MPI_File_set_view(pFile,cOff[myid],MPI_INT,MPI_INT,"native",MPI_INFO_NULL);
  MPI_File_write(pFile, cV, nCells, MPI_INT, &status);

  cOff[myid] += totCells*sizeof(float);
  MPI_File_set_view(pFile,cOff[myid],MPI_INT,MPI_INT,"native",MPI_INFO_NULL);
  MPI_File_write(pFile, cW, nCells, MPI_INT, &status);
#else /* ROOT WRITE */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++) {
      PRF_CRECV_FLOAT(i, &cU[cOff[i]], nCellsL[i], myid);
      PRF_CRECV_FLOAT(i, &cV[cOff[i]], nCellsL[i], myid);
      PRF_CRECV_FLOAT(i, &cW[cOff[i]], nCellsL[i], myid);
    }

    writeDataReal(pFile, totCells, cU, 1);
    writeDataReal(pFile, totCells, cV, 1);
    writeDataReal(pFile, totCells, cW, 1);
  }
  else {
    PRF_CSEND_FLOAT(node_zero, cU, nCellsL[myid], node_zero);
    PRF_CSEND_FLOAT(node_zero, cV, nCellsL[myid], node_zero);
    PRF_CSEND_FLOAT(node_zero, cW, nCellsL[myid], node_zero);
  }
#endif /* UDF_MPI_WRITE */
  free(nCellsL); free(iWork); free(cOff);
#else /* SERIAL */
  writeDataReal(pFile, nCells, cU, 1);
  writeDataReal(pFile, nCells, cV, 1);
  writeDataReal(pFile, nCells, cW, 1);
#endif /* PARALLEL */
  free(cU); free(cV); free(cW);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void UDFwriteEnsightScalar_SP(char *dir, char *prefix, char *name, char *abbrv,
                              int ind, int type)
{
#if !RP_HOST
  Domain *d;
  Thread *t;
  cell_t c;
  int i, j, id, nC, nCells, totCells, *nCellsL, *iWork;
  float *cSc;

  d = Get_Domain(1); /* get the domain from FLUENT */

  char fileName[1000];
  if (ex.single_dir)
    sprintf(fileName, "%s/data/%s-%s.%05i", dir, prefix, abbrv, ex.fileNo);
  else
    sprintf(fileName, "%s/%s/%s-%s.%05i", dir, abbrv, prefix, abbrv, ex.fileNo);
  Message0("  Writing scalar file for: %20s%3s%s\n",
           name, " ", strrchr(fileName,'/')+1);

#if PARALLEL 
#if UDF_MPI_WRITE
  int flag;
  MPI_Status status;
  MPI_File *pFile;
  flag = MPI_File_open(MPI_COMM_SELF, fileName, MPI_MODE_WRONLY |
                       MPI_MODE_CREATE,MPI_INFO_NULL,&pFile);
  MPI_Offset *cOff=(MPI_Offset*) calloc(compute_node_count, sizeof(MPI_Offset));
#else
  FILE *pFile;
  if (myid == node_zero)
    pFile = fopen(fileName, "wb");
  int *cOff = (int*) calloc(compute_node_count, sizeof(int));
#endif /* UDF_MPI_WRITE */
#else /* SERIAL OR PARALLEL SINGLE WRITE */
  FILE *pFile = fopen(fileName, "wb");
#endif /* PARALLEL && UDF_MPI_WRITE */

  if (pFile == NULL)
    Message0("file not opened\n");

  id = 1;
#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      /* write the header */
      writeDataStr(pFile, 80, name);
      writeDataStr(pFile, 80, "part");      /* beginning of part */
      writeDataInt(pFile, 1, &id, 1);
      writeDataStr(pFile, 80, "hexa8");      /* element type */
    }

  /* count the number of cells */
  nCells = 0;
  thread_loop_c(t,d)
    {
      begin_c_loop_int(c,t)
        {
          nCells++;
        }
      end_c_loop_int(c,t)
    }

  nC = nCells;
#if PARALLEL
  nCellsL = (int*) calloc(compute_node_count, sizeof(int));
  iWork   = (int*) calloc(compute_node_count, sizeof(int));
  totCells = PRF_GRSUM1(nCells); /* get the total number of cells */
  /* accumulate all information about cells on every node */
  nCellsL[myid] = nCells;
  PRF_GISUM(nCellsL, compute_node_count, iWork);

  /* compute the offset */
  for (i=1, cOff[0]=0; i<compute_node_count; i++)
    for (j=0, cOff[i]=0; j<i; j++)
      cOff[i] += nCellsL[j];

#if UDF_MPI_WRITE
  cOff[myid] *= sizeof(float);     /* adjust offset to reflect byte size */
  cOff[myid] += 3*80*sizeof(char) + 1*sizeof(int); /* add the header material */

  /* set the file view */
  MPI_File_set_view(pFile,cOff[myid],MPI_INT,MPI_INT,"native",MPI_INFO_NULL);
#else /* PARALLEL SINGLE WRITE */
  if (myid == node_zero) /* make sure root has enough memory for all cells */
    nC = totCells;
#endif /* UDF_MPI_WRITE */
#endif /* PARALLEL */

  /* store the velocities in an array */
  cSc = (float*) calloc(nC, sizeof(float));
  i = 0;
  thread_loop_c(t,d)
    {
      begin_c_loop_int(c,t)
        {
          if (type == 0)
            cSc[i] = (float)(ABS_P(C_P(c,t),op_pres)); /* absolute pressure */
          else if (type == 1)
            cSc[i] = (float)(C_STORAGE_R(c,t,ind));   /* standard */
          else if (type == 2)
            cSc[i] = (float)(C_UDMI(c,t,ind)); /* udm */
          i++;
        }
      end_c_loop_int(c,t)
    }

#if PARALLEL
#if UDF_MPI_WRITE
  /* write the velocities */
  MPI_File_write(pFile, cSc, nCells, MPI_INT, &status);
  MPI_File_close(&pFile);
#else /* PARALLEL SINGLE WRITE */
  /* transfer all the cell information to the root */
  if (myid == node_zero)
    {
      for (i=1; i<compute_node_count; i++)
        PRF_CRECV_FLOAT(i, &cSc[cOff[i]], nCellsL[i], myid);

      writeDataReal(pFile, totCells, cSc, 1);
      fclose(pFile);
    }
  else
    {
      PRF_CSEND_FLOAT(node_zero, cSc, nCellsL[myid], node_zero);
    }
#endif /* UDF_MPI_WRITE */
  free(nCellsL); free(iWork); free(cOff);
#else /* SERIAL */
  writeDataReal(pFile, nCells, cSc, 1);
  fclose(pFile);
#endif /* PARALLEL */
  free(cSc);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void UDFwriteEnsightGeom(char *dir, char *prefix)
{
#if !RP_HOST
  int i, j;

  /* generate the filename */
  char fileName[1000];
  if (ex.single_dir)
    sprintf(fileName, "%s/data/%s-geo.%05i", dir, prefix, ex.fileNo);
  else
    sprintf(fileName, "%s/geometry/%s-geo.%05i", dir, prefix, ex.fileNo);
  Message0("  %-46s  %s\n", "Writing geometry file:", strrchr(fileName, '/')+1);

#if PARALLEL 
#if UDF_MPI_WRITE
  int flag;
  MPI_File *geoFile;
  flag = MPI_File_open(MPI_COMM_SELF, prefix, MPI_MODE_RDWR | MPI_MODE_CREATE,
                       MPI_INFO_NULL,&geoFile);

  MPI_File_set_view(geoFile,(MPI_Offset)0,MPI_INT,MPI_INT,"native",MPI_INFO_NULL);
#else
  FILE *geoFile;
  if (myid == node_zero)
    geoFile = fopen(fileName, "wb");
#endif /* UDF_MPI_WRITE */
#else /* SERIAL */
  FILE *geoFile = fopen(fileName, "wb");
#endif /* PARALLEL  */

#if RP_NODE
  if (myid == node_zero)
#endif /* RP_NODE */
    {
      if (!geoFile)
        Message0("Error: could not open file %s\n", fileName);

      /* write the header information */
      if (ex.binary)
        writeDataStr(geoFile, 80, "C Binary");
      writeDataStr(geoFile, 80, "Model Geometry File exported from Fluent");
      writeDataStr(geoFile, 80, "Ensight-Gold format");
      writeDataStr(geoFile, 80, "node id assign");
      writeDataStr(geoFile, 80, "element id assign");
    }

  /* write the cell geometry */
  for (i=0; i<ex.nCells; i++) {
    /* write the node data for the current thread */
    UDFwriteEnsightNodes(ex.cell_th[i], geoFile, i+1);

    /* write the cell connectivity for the current thread */
    switch (THREAD_ELEMENT_TYPE(ex.cell_th[i])) {
    case MIXED_CELL:
      /* TODO: generalize for poly cells as well */
      for (j=0; j<ex.part_nv[i]; j++)
        UDFwriteEnsightUnstructGeom(ex.cell_th[i], geoFile,
                                    ex.part_verts[i][j],ex.part_type[i][j]);
      break;
      /* 2D cell types */
    case TRI_CELL:
      break;
    case QUAD_CELL:
      break;
      /* 3D cell types */
    case TET_CELL:
      UDFwriteEnsightUnstructGeom(ex.cell_th[i], geoFile, 4, "tetra4");
      break;
    case PYRAMID_CELL:
      UDFwriteEnsightUnstructGeom(ex.cell_th[i], geoFile, 5, "pyramid5");
      break;
    case WEDGE_CELL:
      UDFwriteEnsightUnstructGeom(ex.cell_th[i], geoFile, 6, "penta6");
      break;
    case HEX_CELL:
      UDFwriteEnsightUnstructGeom(ex.cell_th[i], geoFile, 8, "hexa8");
      break;
    case POLY_CELL:
      break;
    } /* end switch */
  }

#if PARALLEL 
#if UDF_MPI_WRITE
  MPI_File_close(geoFile);
#else /* PARALLEL SINGLE WRITE */
  if (myid == node_zero) fclose(geoFile);
#endif /* UDF_MPI_WRITE */
#else /* SERIAL */
  fclose(geoFile);
#endif /* PARALLEL */
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

#if PARALLEL
void UDFwriteEnsightNodes(Thread *t, FILE *pFile, int id)
{
#if !RP_HOST
  Thread *tf;
  cell_t c;
  face_t f;
  Node *v;

  /* counters */
  int i, j, n, nf, iN, nL, totNodes, maxNodes; 
  int *nNodes, *nOff, *iWork;

  /* storage for the nodal data */
  int   *nodeID, *nodeIndx, *iTmpA;
  float *nodeX, *nodeY, *nodeZ;

  float *fTmpA;
  FILE *idFile;
  char buffer[200];

  nNodes = (int*) calloc(compute_node_count, sizeof(int));
  nOff   = (int*) calloc(compute_node_count, sizeof(int));
  iWork  = (int*) calloc(compute_node_count, sizeof(int));

  /* get the total number of nodes and cells internal to the process */
  begin_c_loop_int(c,t) {
    c_node_loop(c,t,n) {
      v = C_NODE(c,t,n);
      NODE_ITMP1(v) = 0;
    }
  } end_c_loop_int(c,t);
    
  /* get the total number of nodes and cells internal to the process */
  for (i=0; i<compute_node_count; i++) nNodes[i] = 0;
  begin_c_loop_int(c,t) { /* loop over all interior cells */
    c_face_loop(c,t,n) { /* loop over this cell's faces */
      f  = C_FACE(c,t,n);        /* get the current face */
      tf = C_FACE_THREAD(c,t,n); /* get the associated face thread */
      if (PRINCIPAL_FACE_P(f,tf)) {
        f_node_loop(f,tf,nf) {   /* loop over nodes of principal face */
          v = F_NODE(f,tf,nf);
          if (NODE_ITMP1(v) == 0) {
            NODE_ITMP1(v) = 1;
            nNodes[myid]++;
          }
        } /* end loop over face nodes */
      } /* end principal face */
    } /* end loop over cell faces */
  } end_c_loop_int(c,t);

  /* accumulate to all processes */
  PRF_GISUM(nNodes, compute_node_count, iWork);
  maxNodes = 0;
  for (i=0, totNodes=0; i<compute_node_count; i++) {
    totNodes += nNodes[i];
    maxNodes = IMAX(maxNodes, nNodes[i]);
  }
  /* allocate enough memory on the root for all nodes, otherwise only need
     for local */
  if (myid == node_zero)
    nL = totNodes;
  else
    nL = nNodes[myid];

  /* allocate the necessary storage for the nodes */
  nodeID   = (int*)   calloc(totNodes, sizeof(int));
  nodeIndx = (int*)   calloc(totNodes, sizeof(int));
  nodeX    = (float*) calloc(nL, sizeof(float));
  nodeY    = (float*) calloc(nL, sizeof(float));
  nodeZ    = (float*) calloc(nL, sizeof(float));

  /* find the offset for each node when transferring to the root */
  for (i=0; i<compute_node_count; i++) nOff[i] = 0;
  for (i=1, nOff[0]=0; i<compute_node_count; i++)
    for (j=0, nOff[i]=0; j<i; j++)
      nOff[i] += nNodes[j];

  /* fill all the node data */
  begin_c_loop_int(c,t) {
    c_node_loop(c,t,n) {
      v = C_NODE(c,t,n);
      NODE_ITMP1(v) = 0;
    }
  } end_c_loop_int(c,t);

  iN = 0;
  begin_c_loop_int(c,t) {
    c_face_loop(c,t,n) { /* loop over this cell's faces */
      f  = C_FACE(c,t,n);        /* get the current face */
      tf = C_FACE_THREAD(c,t,n); /* get the associated face thread */
      if (PRINCIPAL_FACE_P(f,tf)) {
        f_node_loop(f,tf,nf) {  /* loop over nodes of principal face */
          v = F_NODE(f,tf,nf);
          if (NODE_ITMP1(v) == 0) {
            /* store IDs and co-ordinates */
            nodeID[iN] = NODE_ID(v);
            nodeX[iN]  = NODE_X(v);
            nodeY[iN]  = NODE_Y(v);
            nodeZ[iN]  = NODE_Z(v);
            NODE_ITMP1(v) = 1;
            iN++;
          }
        }
      } /* end principal face */
    } /* end loop over cell faces */
  } end_c_loop_int(c,t);

  /* transfer all the node data to the root */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++) {
      PRF_CRECV_INT(i,  &nodeID[nOff[i]], nNodes[i], myid);
      PRF_CRECV_FLOAT(i, &nodeX[nOff[i]], nNodes[i], myid);
      PRF_CRECV_FLOAT(i, &nodeY[nOff[i]], nNodes[i], myid);
      PRF_CRECV_FLOAT(i, &nodeZ[nOff[i]], nNodes[i], myid);
    }
  }
  else {
    PRF_CSEND_INT(node_zero,  nodeID, nNodes[myid], node_zero);
    PRF_CSEND_FLOAT(node_zero, nodeX, nNodes[myid], node_zero);
    PRF_CSEND_FLOAT(node_zero, nodeY, nNodes[myid], node_zero);
    PRF_CSEND_FLOAT(node_zero, nodeZ, nNodes[myid], node_zero);
  }

  /* write the data from the root */
  if (myid == node_zero) {
    /* index array in sorted order of node ident if necessary */
    indexIntArray(totNodes, nodeID, nodeIndx);

    /* write the header information */
    writeDataStr(pFile, 80, "part");
    writeDataInt(pFile, 1, &id, 1);
    writeDataStr(pFile, 80, ex.part_name[id-1]);
    writeDataStr(pFile, 80, "coordinates");
    writeDataInt(pFile, 1, &totNodes, 1);

    /* write the node ids */
    /* writeDataInt(pFile, totNodes, nodeID, 1); */

    /* write the node coordinates */
    writeDataReal(pFile, totNodes, nodeX, 1);
    writeDataReal(pFile, totNodes, nodeY, 1);
    writeDataReal(pFile, totNodes, nodeZ, 1);

    /* send the sorted list and index out to the compute nodes */
    for (i=1; i<compute_node_count; i++) {
      PRF_CSEND_INT(i, nodeID, totNodes, myid);
      PRF_CSEND_INT(i, nodeIndx, totNodes, myid);
    }
  }
  else {
    /* receive entire sorted list from root */
    PRF_CRECV_INT(node_zero, nodeID, totNodes, node_zero);
    PRF_CRECV_INT(node_zero, nodeIndx, totNodes, node_zero);
  }

  /* set all the tmp1 values to zero */
  begin_c_loop_int(c,t) {
    c_node_loop(c,t,n) {
      v = C_NODE(c,t,n);
      NODE_ITMP1(v) = 0;
    }
  } end_c_loop_int(c,t);

  /* store the order node numbering in temporary value used later in 
     unstructured connectivity write (writeEnsightUnstructGeom) */
  begin_c_loop_int(c,t) {
    c_node_loop(c,t,n) {
      v = C_NODE(c,t,n);
      if (NODE_ITMP1(v) == 0) {
        j = gridLocateIntIndx(totNodes,nodeID,nodeIndx,NODE_ID(v));
        if (NODE_ID(v) == nodeID[nodeIndx[totNodes-1]]) j++;
        NODE_ITMP1(v) = nodeIndx[j] + 1;
      }
    }
  } end_c_loop_int(c,t);

  free(nodeID); free(nodeIndx);
  free(nodeX); free(nodeY); free(nodeZ);
  free(nNodes); free(nOff); free(iWork);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */
#else /* SERIAL */
void UDFwriteEnsightNodes(Thread *t, FILE *pFile, int id)
{
#if !RP_HOST
  cell_t c;
  face_t f;
  Node   *v;
  Domain *d = Get_Domain(1);

  int i, n, iC, iN, nCells, nNodes, iTmp = 1, *iTmpA, threadNodes;
  float fTmp, *fTmpA;

  iN = 0;
#if RampantReleaseMajor >= 15
  begin_node_loop(v,d) {
#else
  begin_node_loop(v) {
#endif /* RampantReleaseMajor >= 15 */
    iN++;
    NODE_ID(v) = iN;
  } end_node_loop(v);

  /* set flag on all nodes to zero to help avoid double-counting */
  begin_c_loop_int(c,t) {
    c_node_loop(c,t,n) {
      v = C_NODE(c,t,n);
      NODE_ITMP1(v) = 0;
    }
  } end_c_loop_int(c,t);
    
  /* get the total number of nodes and cells internal to the process */
  threadNodes = 0;
  begin_c_loop_int(c,t) {
    c_node_loop(c,t,n) {
      v = C_NODE(c,t,n);
      if (!NODE_ITMP1(v)) {
        NODE_ITMP1(v) = 1;
        threadNodes++;
      }
    }
  } end_c_loop_int(c,t);

  /* write the header information */
  writeDataStr(pFile, 80, "part");
  writeDataInt(pFile, 1, &id, 1);
  writeDataStr(pFile, 80, ex.part_name[id-1]);
  writeDataStr(pFile, 80, "coordinates");
  writeDataInt(pFile, 1, &threadNodes, 1);

  /* output the node ids and co-ordinates */
  iTmpA = (int*)   calloc(threadNodes, sizeof(int));
  fTmpA = (float*) calloc(threadNodes, sizeof(float));
  /* x co-ordinates */
  iN = 0;
#if RampantReleaseMajor >= 15
  begin_node_loop(v,d) { NODE_ITMP1(v) = 0; } end_node_loop(v)
#else
  begin_node_loop(v) { NODE_ITMP1(v) = 0; } end_node_loop(v)
#endif /* RampantReleaseMajor >= 15 */
  begin_c_loop_int(c,t) {
    c_node_loop(c,t,n) {
      v = C_NODE(c,t,n);
      if (!NODE_ITMP1(v)) {
        NODE_ITMP1(v) = 1;
        iTmpA[iN] = NODE_ID(v);
        fTmpA[iN] = NODE_X(v);
        iN++;
      }
    }
  } end_c_loop_int(c,t);
  /* fwrite(iTmpA, sizeof(int),   nNodes, pFile); */
 /*  writeDataInt(pFile, threadNodes, iTmpA); */
  writeDataReal(pFile, threadNodes, fTmpA, 1);
  /* y co-ordinates */
  iN = 0;
#if RampantReleaseMajor >= 15
  begin_node_loop(v,d) { NODE_ITMP1(v) = 0; } end_node_loop(v)
#else
  begin_node_loop(v) { NODE_ITMP1(v) = 0; } end_node_loop(v)
#endif /* RampantReleaseMajor >= 15 */
  begin_c_loop_int(c,t) {
      c_node_loop(c,t,n) {
          v = C_NODE(c,t,n);
          if (!NODE_ITMP1(v)) {
              NODE_ITMP1(v) = 1;
              fTmpA[iN] = NODE_Y(v);
              iN++;
            }
        }
  } end_c_loop_int(c,t);
  writeDataReal(pFile, threadNodes, fTmpA, 1);
  /* Z co-ordinates */
  iN = 0;
#if RampantReleaseMajor >= 15
  begin_node_loop(v,d) { NODE_ITMP1(v) = 0; } end_node_loop(v)
#else
  begin_node_loop(v) { NODE_ITMP1(v) = 0; } end_node_loop(v)
#endif /* RampantReleaseMajor >= 15 */
  begin_c_loop_int(c,t) {
      c_node_loop(c,t,n) {
          v = C_NODE(c,t,n);
          if (!NODE_ITMP1(v)) {
              NODE_ITMP1(v) = 1;
              fTmpA[iN] = NODE_Z(v);
              iN++;
            }
        }
  } end_c_loop_int(c,t);
  writeDataReal(pFile, threadNodes, fTmpA, 1);
  iN = 0;
#if RampantReleaseMajor >= 15
  begin_node_loop(v,d) { NODE_ITMP1(v) = 0; } end_node_loop(v)
#else
  begin_node_loop(v) { NODE_ITMP1(v) = 0; } end_node_loop(v)
#endif /* RampantReleaseMajor >= 15 */
  begin_c_loop_int(c,t) {
      c_node_loop(c,t,n) {
          v = C_NODE(c,t,n);
          if (NODE_ITMP1(v) == 0) {
              iN++;
              NODE_ITMP1(v) = iN;
            }
        }
  } end_c_loop_int(c,t);
  free(iTmpA); free(fTmpA);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */
#endif /* PARALLEL */

void UDFwriteEnsightUnstructGeom(Thread *t, FILE *pFile, int verts, char *cell)
{
#if !RP_HOST
  cell_t c;
  Node *v;
  int i, j, n, iC, iN, nC, totCells;

#if !PARALLEL
  int myid = 0;
#endif /* !PARALLEL */

  int *nCells = (int*) calloc(compute_node_count, sizeof(int));
  int *cOff   = (int*) calloc(compute_node_count, sizeof(int));
  int *iWork  = (int*) calloc(compute_node_count, sizeof(int));

  /* count the number of cells with the correct vertices for this thread */
  nCells[myid] = 0;
  begin_c_loop_int(c,t) {
    if (C_NNODES(c,t) == verts)
      nCells[myid]++;
  } end_c_loop_int(c,t);
  
#if PARALLEL
  PRF_GISUM(nCells, compute_node_count, iWork);
  for (i=0, totCells=0; i<compute_node_count; i++)
    totCells += nCells[i];

  nC = verts*nCells[myid];
  /* find the offset for each node when transferring to the root */
  for (i=1, cOff[0]=0; i<compute_node_count; i++)
    for (j=0, cOff[i]=0; j<i; j++)
      cOff[i] += verts*nCells[j];

#if UDF_MPI_WRITE
  /* change offset to bytes */
  cOff[myid] *= sizeof(int);

  /* advance the front material */
  cOff[myid] += 8*80*sizeof(char) + 2*sizeof(int);     /* header */
  cOff[myid] += 3*totNodes*sizeof(float);              /* coordinates */
  cOff[myid] += 80*sizeof(char) + 1*sizeof(int);       /* element type */

  MPI_File_set_view(pFile,cOff[myid],MPI_INT,MPI_INT,"native",MPI_INFO_NULL);
#else
  if (myid == node_zero) nC = verts*totCells;
#endif /* UDF_MPI_WRITE */
#else /* SERIAL */
  totCells = nCells[myid];
  nC = verts*totCells;
#endif /* PARALLEL */

#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      /* write the part cell type and number of cells */
      writeDataStr(pFile, 80, cell);
      writeDataInt(pFile, 1, &totCells, 1);
    }

  /* allocate an array for the cell connectivity */
  int *conn = (int*) calloc(nC, sizeof(int));

  /* fill the connectivity array */
  iC = 0;
  begin_c_loop_int(c,t) {
    iN = 0;
    if (C_NNODES(c,t) == verts) {
      c_node_loop(c,t,n) {
        v = C_NODE(c,t,n);
        conn[verts*iC+iN] = NODE_ITMP1(v);
        iN++;
      }
      iC++;
    }
  } end_c_loop_int(c,t);

  /* collect all the cells on the root and write to file in sequence */
#if PARALLEL
#if UDF_MPI_WRITE
  MPI_Status status;
  MPI_File_write(pFile, conn, verts*nCells[myid], MPI_INT, &status);
  MPI_File_close(&pFile);
#else /* ROOT WRITE */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++)
      PRF_CRECV_INT(i, &conn[cOff[i]], verts*nCells[i], myid);
    
    /* write the connectivity */
    writeDataInt(pFile, verts*totCells, conn, verts);
  }
  else {
    PRF_CSEND_INT(node_zero, conn, verts*nCells[myid], node_zero);
  }
#endif /* UDF_MPI_WRITE */
#else /* SERIAL */
  writeDataInt(pFile, verts*totCells, conn, verts);
#endif /* PARALLEL */
  free(conn); free(nCells); free(cOff); free(iWork);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void UDFwriteEnsightScalar(char *dir,char *prefix,char *name,char *abbrv,
                           int ind, int type)
{
#if !RP_HOST
  int i, j, id;

  char fileName[1000];
  if (ex.single_dir)
    sprintf(fileName, "%s/data/%s-%s.%05i", dir, prefix, abbrv, ex.fileNo);
  else
    sprintf(fileName, "%s/%s/%s-%s.%05i", dir, abbrv, prefix, abbrv, ex.fileNo);
  Message0("  Writing scalar file for: %20s%3s%s\n",
           name, " ", strrchr(fileName,'/')+1);

  FILE *pFile;

#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    { /* write the header */
      pFile = fopen(fileName, "wb");
      writeDataStr(pFile, 80, name);
    }

  /* write the cell geometry */
  for (i=0; i<ex.nCells; i++) {
#if PARALLEL
    if (myid == node_zero)
#endif /* PARALLEL */
      { /*  write header for this part */
        writeDataStr(pFile, 80, "part"); /* beginning of part */
        id = i+1;
        writeDataInt(pFile, 1, &id, 1);
      }
    /* write the cell connectivity for the current thread */
    switch (THREAD_ELEMENT_TYPE(ex.cell_th[i])) {
    case MIXED_CELL:
      for (j=0; j<ex.part_nv[i]; j++)
        writeEnsightScalar(ex.cell_th[i],pFile,ex.part_verts[i][j],
                           ex.part_type[i][j],ind,type);
      break;
    default:
      writeEnsightScalar(ex.cell_th[i],pFile,ex.part_verts[i][0],
                         ex.part_type[i][0],ind,type);
      break;
    } /* end switch */
    
  }
  
#if PARALLEL
  if (myid == node_zero)
#endif /*  PARALLEL */
  {
    fclose(pFile);
  }
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void writeEnsightScalar(Thread *t, FILE *pFile, int verts, char *cell, int ind,
                        int type)
{
#if !RP_HOST
  cell_t c;
  int i, j, nC, nCells, totCells, *nCellsL, *iWork;
  float *cSc;

#if PARALLEL
#if UDF_MPI_WRITE
  MPI_Status status;
  MPI_Offset *cOff=(MPI_Offset*) calloc(compute_node_count, sizeof(MPI_Offset));
#else /* PARALLEL SINGLE WRITE */
  int *cOff = (int*) calloc(compute_node_count, sizeof(int));
#endif /* UDF_MPI_WRITE */
#else /* SERIAL */

#endif /* PARALLEL */

#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      writeDataStr(pFile, 80, cell);      /* element type */
    }

  /* count the number of cells */
  nCells = 0;
  begin_c_loop_int(c,t) {
    if (C_NNODES(c,t) == verts)
      nCells++;
  } end_c_loop_int(c,t);

  nC = nCells;
#if PARALLEL
  nCellsL = (int*) calloc(compute_node_count, sizeof(int));
  iWork   = (int*) calloc(compute_node_count, sizeof(int));
  totCells = PRF_GRSUM1(nCells); /* get the total number of cells */
  /* accumulate all information about cells on every node */
  nCellsL[myid] = nCells;
  PRF_GISUM(nCellsL, compute_node_count, iWork);

  /* compute the offset */
  for (i=1, cOff[0]=0; i<compute_node_count; i++)
    for (j=0, cOff[i]=0; j<i; j++)
      cOff[i] += nCellsL[j];

#if UDF_MPI_WRITE
  cOff[myid] *= sizeof(float);     /* adjust offset to reflect byte size */
  cOff[myid] += 3*80*sizeof(char) + 1*sizeof(int); /* add the header material */

  /* set the file view */
  MPI_File_set_view(pFile,cOff[myid],MPI_INT,MPI_INT,"native",MPI_INFO_NULL);
#else /* PARALLEL SINGLE WRITE */
  if (myid == node_zero) /* make sure root has enough memory for all cells */
    nC = totCells;
#endif /* UDF_MPI_WRITE */
#endif /* PARALLEL */

  /* store the velocities in an array */
  cSc = (float*) calloc(nC, sizeof(float));
  i = 0;
  begin_c_loop_int(c,t) {
    if (C_NNODES(c,t) == verts) {
      if (type == 0)
        cSc[i] = (float)(ABS_P(C_P(c,t),op_pres)); /* absolute pressure */
      else if (type == 1)
        cSc[i] = (float)(C_STORAGE_R(c,t,ind));    /* standard */
      else if (type == 2)
        cSc[i] = (float)(C_UDMI(c,t,ind));         /* udm */
      i++;
    }
  } end_c_loop_int(c,t);

#if PARALLEL
#if UDF_MPI_WRITE
  /* write the velocities */
  MPI_File_write(pFile, cSc, nCells, MPI_INT, &status);
  MPI_File_close(&pFile);
#else /* PARALLEL SINGLE WRITE */
  /* transfer all the cell information to the root */
  if (myid == node_zero) {
    for (i=1; i<compute_node_count; i++)
      PRF_CRECV_FLOAT(i, &cSc[cOff[i]], nCellsL[i], myid);
    
      writeDataReal(pFile, totCells, cSc, 1);
  }
  else {
    PRF_CSEND_FLOAT(node_zero, cSc, nCellsL[myid], node_zero);
  }
#endif /* UDF_MPI_WRITE */
  free(nCellsL); free(iWork); free(cOff);
#else /* SERIAL */
  writeDataReal(pFile, nCells, cSc, 1);
#endif /* PARALLEL */
  free(cSc);
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void UDFwriteEnsightParticle(char *dir, char *prefix, int nStream,
                             Injection **stream)
{
#if !RP_HOST
  /* injection variables */
  Injection *I, *I_list;
  Particle  *p;
  int i, j, iP, nP, totParticles;

  /* storage variables */
  int   iTmp, *ids;
  float fTmp, *coords, *vel, *temp, *diam;

  /* currently writes all particles to a single data set */
  I_list = Get_dpm_injections();

  /* construct the file names */
  FILE *mpgFile, *velFile, *diamFile, *tempFile;
  char mpgName[1000], velName[1000], tempName[1000], diamName[1000];
  sprintf(mpgName, "%s/particle/%s-mpg.%05i",  dir, prefix, ex.fileNo);
  sprintf(velName, "%s/particle/%s-vel.%05i",  dir, prefix, ex.fileNo);
  sprintf(tempName,"%s/particle/%s-temp.%05i", dir, prefix, ex.fileNo);
  sprintf(diamName,"%s/particle/%s-diam.%05i", dir, prefix, ex.fileNo);

  /* count the number of particles*/
  nP = 0; 
  loop(I, I_list) {
    loop(p, I->p) {
      nP++; /* number of particles on this processor */
    }
  }
  totParticles = PRF_GISUM1(nP);

  Message0("  Writing particle data:   %20i%3sparticle/%s*\n", totParticles,
           " ", prefix);

  /* write the headers for each file*/
#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      /* geometry file */
      mpgFile  = fopen(mpgName, "wb");
      writeDataStr(mpgFile, 80, "C Binary");
      writeDataStr(mpgFile, 80, "Particle location for all injections");
      writeDataStr(mpgFile, 80, "particle coordinates");
      writeDataInt(mpgFile, 1, &totParticles, 1);
      
      /* velocity file */
      velFile  = fopen(velName, "wb");
      writeDataStr(velFile, 80, "C Binary");
      writeDataStr(velFile, 80, "particle velocity magnitude");
      
      /* temperature file */
      tempFile = fopen(tempName, "wb");
      writeDataStr(tempFile, 80, "C Binary");
      writeDataStr(tempFile, 80, "particle temperature");
      
      /* diameter file */
      diamFile = fopen(diamName, "wb");
      writeDataStr(diamFile, 80, "C Binary");
      writeDataStr(diamFile, 80, "particle diameter");
      fclose(mpgFile); fclose(velFile); fclose(tempFile); fclose(diamFile);
    }
  PRF_GSYNC();

  /* write the node ids*/
#if PARALLEL
  for (i=0; i<compute_node_count; i++) {
    if (i == myid)
#endif /* PARALLEL */
      {
        mpgFile  = fopen(mpgName, "ab");
        
        ids = (int*) calloc(nP, sizeof(int));
        iP = 0;
        loop(I, I_list) {
            loop(p, I->p) {
              ids[iP] = p->part_id;
              iP++;
            } /* end particle loop */
        } /* end loop over streams */
        writeDataInt(mpgFile, nP, ids, 5);
        free(ids);
        fclose(mpgFile);
      }
#if PARALLEL
    PRF_GSYNC();
  } /* end loop over compute nodes */
#endif /* PARALLEL */

  /* write the geometry as tuples and the scalar values */
#if PARALLEL
  for (i=0; i<compute_node_count; i++) {
    if (i == myid)
#endif /* PARALLEL */
      {
        mpgFile  = fopen(mpgName, "ab");
        velFile  = fopen(velName, "ab");
        tempFile = fopen(tempName, "ab");
        diamFile = fopen(diamName, "ab");
        
        coords = (float*) calloc(3*nP, sizeof(float));
        vel    = (float*) calloc(nP, sizeof(float));
        temp   = (float*) calloc(nP, sizeof(float));
        diam   = (float*) calloc(nP, sizeof(float));
        
        iP = 0;
        loop(I, I_list) {
          loop(p, I->p) {
            /* store the geometry */
            coords[3*iP]   = (float)(P_POS(p)[0]);
            coords[3*iP+1] = (float)(P_POS(p)[1]);
            coords[3*iP+2] = (float)(P_POS(p)[2]);
            
            /* store the scalars */
            fTmp=sqrt(SQR(P_VEL(p)[0])+SQR(P_VEL(p)[1])+SQR(P_VEL(p)[2]));
            vel[iP]  = (float)(fTmp);
            temp[iP] = (float)(P_T(p));
            diam[iP] = (float)(P_DIAM(p));
            iP++;
          } /* end particle loop */
        } /* end loop over streams */
        
          /* write the data */
        writeDataReal(mpgFile, 3*nP, coords, 3);
        writeDataReal(velFile, nP, vel, 5);
        writeDataReal(tempFile, nP, temp, 5);
        writeDataReal(diamFile, nP, diam, 5);
        
        /* free local memory and close files */
        free(coords); free(vel); free(temp); free(diam);
        fclose(mpgFile); fclose(velFile); fclose(tempFile); fclose(diamFile);
      }
#if PARALLEL
    PRF_GSYNC();
  } /* end loop over compute nodes */
#endif /* PARALLEL */
#endif /* !RP_HOST */
}
/* -------------------------------------------------------------------------- */

void UDFwriteEnsightInjection(char *dir, char *prefix, char *ident,
                              int nNozzles, Injection **stream)
{
  Particle *p;
  int i, count=0, perLine=6, line;
  char mpgName[1000], tempName[1000], diamName[1000];

  sprintf(mpgName,"%s/particle/%s-%s-mpg.%05i", dir, prefix, ident, ex.fileNo);
  sprintf(tempName,"%s/particle/%s-%s-temp.%05i", dir, prefix, ident,ex.fileNo);
  sprintf(diamName,"%s/particle/%s-%s-diam.%05i", dir, prefix, ident,ex.fileNo);

  /* count the number of particles*/
  for (i=0; i<nNozzles; i++)
    count += stream[i]->n_particles;

  /* write the particle geometry file */
  FILE *mpgFile  = fopen(mpgName, "w");
  FILE *tempFile = fopen(tempName, "w");
  FILE *diamFile = fopen(diamName, "w");
  fprintf(mpgFile, "Particle location for %s injection\n", ident);
  fprintf(mpgFile, "particle coordinates\n");
  fprintf(mpgFile, "%10d\n", count);
  fprintf(tempFile, "%s_Temperature\n", ident);
  fprintf(diamFile, "%s_Diameter\n", ident);

  line = 1;
  for (i=0; i<nNozzles; i++) {
    loop(p, stream[i]->p) {
      /* write particle location */
      fprintf(mpgFile, "%10d%12.5e%12.5e%12.5e\n", p->part_id,
              P_POS(p)[0], P_POS(p)[1], P_POS(p)[2]);
      fprintf(tempFile, "%12.5e", P_T(p));
      fprintf(diamFile, "%12.5e", P_DIAM(p));
      
      if (line%perLine == 0) {
        fprintf(tempFile, "\n");
        fprintf(diamFile, "\n");
      }
      line++;
    }
  }
  Message("line = %i\n", line);
  
  fclose(mpgFile); fclose(tempFile); fclose(diamFile);
}
/* -------------------------------------------------------------------------- */
