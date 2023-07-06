#include <stdlib.h>
#include <string.h>
#include "vtkIO.h"
#include "utils.h"

void VTKwriteFlamelet(char *fileName, ioStruct *io, driveStruct *driver,
                      timeStruct *time, flameStruct *flame, flowStruct *flow)
{
  int i,j,bytes;
  char x1Name[32],x2Name[32],strTmp[64];
  float *ftmp = (float*) calloc(flame->n1*flame->n2*flame->nVars,sizeof(float));

  if (flame->nDim == 1) {
    if (driver->H_DIM == 1)
      strcpy(x1Name, "H1");
    else
      strcpy(x1Name, "Z1");
    strcpy(x2Name, "y");
  }
  else if (flame->nDim == 2) {
    if (driver->H_DIM == 0) {
      if (driver->TRANS_Z2 == 0) {
        strcpy(x1Name, "Z1");
        strcpy(x2Name, "Z2");
      }
      else if (driver->TRANS_Z2 == 1) {
        strcpy(x1Name, "Z1");
        strcpy(x2Name, "Z2/(1-Z1)");
      }
      else if (driver->TRANS_Z2 == 2) {
        strcpy(x1Name, "Z1+Z2");
        strcpy(x2Name, "Z2/(Z1+Z2)");
      }
    }
    else if (driver->H_DIM == 1) {
      strcpy(x1Name, "H1");
      strcpy(x2Name, "Z2");
    }
    else if (driver->H_DIM == 2) {
      strcpy(x1Name, "Z1");
      strcpy(x2Name, "H2");
    }
  }

  FILE *dF = fopen(fileName, "w");
  
  /* ----- write the header data ----- */
  fprintf(dF, "<?xml version=\"1.0\"?>\n");
  fprintf(dF, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" ");
  if (io->byteOrder == 0)
    fprintf(dF, "byte_order=\"LittleEndian\"");
  else if (io->byteOrder == 1)
    fprintf(dF, "byte_order=\"BigEndian\"");
  fprintf(dF, ">\n");

  /* ----- mesh info ----- */
  int offset = 0;
  fprintf(dF, "%2s<RectilinearGrid WholeExtent=\"", "");
  fprintf(dF, "0 %i 0 %i 0 0\">\n",flame->n1-1,flame->n1-1);
  /* field data */
  fprintf(dF, "%4s<FieldData>\n","");
  VTKmetaDataArray(dF,6,"Float32","time_[s]",1,"ascii",offset);
  fprintf(dF, "%6s%15.6e\n","",time->t);
  fprintf(dF, "%6s</DataArray>\n","");
/*   VTKmetaDataArray(dF,6,"Float32","pressure_[Pa]",1,"ascii",offset); */
/*   fprintf(dF, "%6s%15.6e\n","",flow->P); */
/*   fprintf(dF, "%6s</DataArray>\n",""); */
  fprintf(dF, "%4s</FieldData>\n","");
  fprintf(dF, "%4s<Piece Extent=\"0 %i 0 %i 0 0\">\n","",flame->n1-1,flame->n2-1);
  /* write the co-ordinate arrays */
  fprintf(dF, "%6s<Coordinates>\n","");
  VTKmetaDataArray(dF,8,"Float32",x1Name,1,"appended",offset);
  offset += sizeof(int) + flame->n1*sizeof(float);
  VTKmetaDataArray(dF,8,"Float32",x2Name,1,"appended",offset);
  offset += sizeof(int) + flame->n2*sizeof(float);
  VTKmetaDataArray(dF,8,"Float32","x3",1,"appended",offset);
  offset += sizeof(int) + sizeof(float);
  fprintf(dF, "%6s</Coordinates>\n","");
  /* ----- write point data ----- */
  fprintf(dF, "%6s<PointData>\n","");
  /* temperature */
  VTKmetaDataArray(dF,8,"Float32","T_[K]",1,"appended",offset);
  offset += sizeof(int) + flame->n1*flame->n2*sizeof(float);
  /* species */
  for (i=0; i<io->nOutputs; i++) {
    sprintf(strTmp,"massfraction-%s",io->outputSpecies[i]);
    VTKmetaDataArray(dF,8,"Float32",strTmp,1,"appended",offset);
    offset += sizeof(int) + flame->n1*flame->n2*sizeof(float);
  }
  if (io->outputChi) {
    if (flame->chiA != NULL)
      VTKmetaDataArray(dF,8,"Float32","chiA_[1/s]",1,"appended",offset);
    if (flame->chiB != NULL)
      VTKmetaDataArray(dF,8,"Float32","chiB_[1/s]",1,"appended",offset);
    if (flame->chiC != NULL)
      VTKmetaDataArray(dF,8,"Float32","chiC_[1/s]",1,"appended",offset);
  }
  if (io->outputHtot)
    VTKmetaDataArray(dF,8,"Float32","HTOT_[J/kg]",1,"appended",offset);
  if (io->outputHR)
    VTKmetaDataArray(dF,8,"Float32","HR-rate_[J/(kg-s)]",1,"appended",offset);
  if (driver->output_tsrc) {
    VTKmetaDataArray(dF,8,"Float32","Tsrc_conv",1,"appended",offset);
    VTKmetaDataArray(dF,8,"Float32","Tsrc_chem",1,"appended",offset);
    VTKmetaDataArray(dF,8,"Float32","Tsrc_diff",1,"appended",offset);
  }
  if (io->outputConv)
    VTKmetaDataArray(dF,8,"Float32","Conv_coeff",1,"appended",offset);
  if (io->outputRho)
    VTKmetaDataArray(dF,8,"Float32","rho_[kg/m^3]",1,"appended",offset);
  if (io->outputMu)
    VTKmetaDataArray(dF,8,"Float32","mu_[Pa-s]",1,"appended",offset);
  if (io->outputCp)
    VTKmetaDataArray(dF,8,"Float32","cp_[J/(kg-K)]",1,"appended",offset);
  if (io->outputPDF)
    VTKmetaDataArray(dF,8,"Float32","pdf",1,"appended",offset);
  if (driver->Q_LOSS > 0)
    VTKmetaDataArray(dF,8,"Float32","qdot",1,"appended",offset);
  fprintf(dF, "%6s</PointData>\n","");
  fprintf(dF, "%4s</Piece>\n","");
  fprintf(dF, "%2s</RectilinearGrid>\n","");
  
  /* ----- write the appended data ----- */
  fprintf(dF, "%2s<AppendedData encoding=\"raw\">\n_", "");
  
  /* write co-ordinates */
  /* Z1 */
  VTKwriteDataArray(dF,1,flame->n1,1,0,flame->x1);
  /* Z2 */
  if (flame->nDim == 2)
    VTKwriteDataArray(dF,1,flame->n2,1,0,flame->x2);
  else {
    bytes = sizeof(int) + sizeof(float);
    ftmp[0] = 0.0;
    fwrite(ftmp,sizeof(float),1,dF);
  }
  /* dummy Z3 */
  bytes = sizeof(int) + sizeof(float);
  fwrite(&bytes,sizeof(int),1,dF);
  ftmp[0] = 0.0;
  fwrite(ftmp,sizeof(float),1,dF);
  
  /* ----- write values ----- */
  /* temperature */
  if (driver->ENERGY)
    VTKwriteDataArray(dF,flame->n2,flame->n1,1,0,flame->T);
  else
    VTKwriteDataArray(dF,flame->n2,flame->n1,flame->nVars,flame->sT,flame->Y);
  /* species */
  for (i=0; i<io->nOutputs; i++) {
    VTKwriteDataArray(dF,flame->n2,flame->n1,flame->nVars,io->outputIndex[i],
                      flame->Y);
  }
  if (io->outputChi) { /* print the scalar dissipation rate */
    if (flame->chiA != NULL) {
      if (driver->CHI_MODE == 2)
        VTKwriteDataArray(dF,flame->n2,flame->n1,flame->nVars,flame->sT,
                          flame->Y);
      else
        VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->chiA);
    }
    if (flame->nDim == 2) {
      if (flame->chiB != NULL)
        VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->chiB);
      if (flame->chiC != NULL)
        VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->chiC);
    }
  }
  if (io->outputHtot) { /* print the total enthalpy */
    if (driver->ENERGY)
      VTKwriteDataArray(dF,flame->n2,flame->n1,flame->nVars,flame->sT,flame->Y);
    else
      VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->H);
  }
  if (io->outputHR)
    VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->HRrate);
  if (driver->output_tsrc) { /* temperature equation terms */
    VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->conv_src);
    VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->chem_src);
    VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->diff_src);
  }
  if (io->outputConv)
    VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->convCoeff);
  if (io->outputRho)
    VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->rho);
  if (io->outputMu)
    VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->mu);
  if (io->outputCp)
    VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->cpmix);
  if (io->outputPDF)
    VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flow->pdf);
  if (driver->Q_LOSS > 0)
    VTKwriteDataArray(dF, flame->n2, flame->n1, 1, 0, flame->srcQ);

/*   /\* write soot source terms if present *\/ */
/*   char sootMomName[7][32] = {"number-density","volume","surface","volume-var", */
/*                                "vol-surf-var","surface-var","delta-weight"}; */
/*   char sootMomAbbrv[7][8] = {"M00", "M10", "M01", "M20", "M11", "M02", "N0"}; */
/*   char sootMomUnit[7][16] = {"1/m^3","m^3/m^3","m^2/m^3","m^6/m^3","m^5/m^3", */
/*                                 "m^4/m^3","1/m^3"}; */
/*   char sootMomSrcUnit[7][16] = {"1/(m^3-s)","m^3/(m^3-s)","m^2/(m^3-s)", */
/*                                 "m^6/(m^3-s)","m^5/(m^3-s)","m^4/(m^3-s)", */
/*                                 "1/(m^3-s)"}; */
/*   char sootSrcNames[6][32] = {"nucleation","condensation","coagulation", */
/*                               "surface-growth","oxidation","fragmentation"}; */
/*   int  sI[7] = {0, 1, 2, 3, 4, 5, 6}; */

/*   if (flame->nSootMom > 0) { */
/*     if (flame->nSootMom < 5) */
/*       sI[3] = 6; */

/*     /\* ----- moments ----- *\/ */
/*     for (i=0; i<flame->nSootMom; i++) { */
/*       fprintf(dF,"soot-%s_%s [%s]\n",sootMomName[sI[i]],sootMomAbbrv[sI[i]], */
/*               sootMomUnit[sI[i]]); */
/*       VTKwriteDataArray(dF,flame->n2,flame->n1,flame->nSootMom,sI[i],perLine, */
/*                  flame->sootMom); */
/*     } */

/*     /\* ----- moment source terms ----- *\/ */
/*     for (i=0; i<flame->nSootMom; i++) { */
/*       fprintf(dF,"soot-%s-src_%s [%s]\n",sootMomName[sI[i]], */
/*               sootMomAbbrv[sI[i]], sootMomSrcUnit[sI[i]]); */
/*       VTKwriteDataArray(dF,flame->n2,flame->n1,flame->nSootMom,sI[i],perLine, */
/*                  flame->sootMomSrc); */
/*     } */
          
/*     if (io->outputSootSrc) { /\* soot source term budget *\/ */
/*       for (j=0; j<6; j++) { */
/*         for (i=0; i<flame->nSootMom-1; i++) { */
/*           fprintf(dF,"soot-%s-%s_%s [%s]\n",sootMomName[i], */
/*                   sootSrcNames[j],sootMomAbbrv[i],sootMomSrcUnit[i]); */
/*           VTKwriteDataArray(dF,flame->n2,flame->n1,flame->nSootMom*6, */
/*                      j*flame->nSootMom+i,perLine,flame->sootMomSrcPart); */
/*         } */
/*         fprintf(dF,"soot-%s-%s_%s [%s]\n",sootMomName[6],sootSrcNames[j], */
/*                 sootMomAbbrv[6],sootMomSrcUnit[6]); */
/*         VTKwriteDataArray(dF,flame->n2,flame->n1,flame->nSootMom*6, */
/*                    j*flame->nSootMom+(flame->nSootMom-1),perLine, */
/*                    flame->sootMomSrcPart); */
/*       } */
/*     } */
/*   } */

    
  /* ----- write arrays ----- */

  fprintf(dF, "\n%2s</AppendedData>\n", "");
  fprintf(dF, "</VTKFile>\n");  
  fclose(dF);
}
/* -------------------------------------------------------------------------- */

void readVTK(char *fileName)
{

}
/* -------------------------------------------------------------------------- */

void VTKrealFieldData(FILE *pF,int tab,char *type,char *name,int n,char *format,
                      long long int offset,double *val)
{
  int i;
  fprintf(pF, "%*s<FieldData ", tab, "");
  fprintf(pF, "type=\"%s\" ", type);
  fprintf(pF, "Name=\"%s\" ", name);
  if (n > 1)
    fprintf(pF, "NumberOfTuples=\"%i\" ", n);
  fprintf(pF, "format=\"%s\" ", format);
  if (!strcmp(format, "appended")) {
    fprintf(pF, "offset=\"%Li\" ", offset);
    fprintf(pF, "/>\n");
  }
  else {
    fprintf(pF, ">\n%*s ",tab,"");
    for (i=0; i<n; i++)
      fprintf(pF, "%15.6e ", val[i]);
    fprintf(pF,"\n</FieldData>\n");
  }
}
/* -------------------------------------------------------------------------- */

void VTKmetaDataArray(FILE *pF,int tab,char *type,char *name,int n,char *format, 
                     long long int offset)
{
  fprintf(pF, "%*s<DataArray ", tab, "");
  fprintf(pF, "type=\"%s\" ", type);
  fprintf(pF, "Name=\"%s\" ", name);
  if (n > 1)
    fprintf(pF, "NumberOfComponents=\"%i\" ", n);
  fprintf(pF, "format=\"%s\" ", format);
  if (!strcmp(format, "appended")) {
    fprintf(pF, "offset=\"%Li\" ", offset);
    fprintf(pF, "/>\n");
  }
  else {
    fprintf(pF, ">\n");
  }
}
/* -------------------------------------------------------------------------- */

void VTKwriteDataArray(FILE *pF,int m,int n,int offset,int index,double *data)
{
  int j,k,bytes;
  float *ftmp = (float*) calloc(m*n, sizeof(float));

  /* write byte size of array */
  bytes = sizeof(int) + n*m*sizeof(float);
  fwrite(&bytes,sizeof(int),1,pF);

  /* write data */
  for (k=0; k<m; ++k)
    for (j=0; j<n; ++j)
      ftmp[k*n+j] = (float)data[(k*n+j)*offset+index];
  fwrite(ftmp,sizeof(float),n*m,pF);
}
/* -------------------------------------------------------------------------- */

void readFlameletVTK(char *fileName, ioStruct *io, driveStruct *driver,
                     timeStruct *time, flameStruct *flame, flowStruct *flow)
{
  int   i, j, bytes;
  int   readFlame; // flag for ending the reading of a flamelet

  // strings to store the tag for each line and comments (if any)
  char  tag[50], comment[200], tmpC;

  int nNodes;
  float *nodes;
  int nCV, *cv, nOffs, *offs, nType;
  char *types;
  
  
  FILE   *pFile;   // input stream
  pFile = fopen(fileName, "r"); // open input file

  if (!pFile) {
    printf("ERROR: cannot open  file %s\n", fileName);
    exit(1);
  }
  
  /* ----- loop until end of file is reached  ----- */
  while (!feof(pFile)) { 
    fscanf(pFile, "%s", tag); // read in the first string

    if (tag[0] == '!' || tag[0] == '#')
      fgets(comment, 200, pFile);
    else if (!strcmp(tag, "<AppendedData")) {
      fgets(comment, 200, pFile);
      fprintf(stderr, "%s\n", comment);

      /* read in the underscore */
      tmpC = fgetc(pFile);
      fprintf(stderr, "Underscore: %c\n", tmpC);

      /* get the node array */
      fread(&bytes, 1, sizeof(int), pFile);
      bytes -= 4;
      nodes = (float*) malloc(bytes);
      nNodes = (int)bytes/sizeof(float);
      fprintf(stderr, "Reading %i bytes (%i floats)\n", bytes, nNodes);
      fread(nodes, 1, bytes, pFile);

      fprintf(stderr, "%03i: ", 0);
      for (i=0; i<nNodes; i++) {
        fprintf(stderr, "%15.6e ", nodes[i]);
        if ((i+1)%3 == 0) {
          fprintf(stderr, "\n%03i: ", (i+1)/3);
        }
      }
      fprintf(stderr, "\n");

      fread(&bytes, 1, sizeof(int), pFile);
      bytes -= 4;
      cv = (int*) malloc(bytes);
      nCV = (int)bytes/sizeof(int);
      fprintf(stderr, "Reading %i bytes (%i ints)\n", bytes, nCV);
      fread(cv, 1, bytes, pFile);

      for (i=0; i<nCV; i++) {
        fprintf(stderr, "%i ", cv[i]);
        if ((i+1)%8 == 0) fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");


      fread(&bytes, 1, sizeof(int), pFile);
      bytes -= 4;
      offs = (int*) malloc(bytes);
      nOffs = (int)bytes/sizeof(int);
      fprintf(stderr, "Reading %i bytes (%i ints)\n", bytes, nOffs);
      fread(offs, 1, bytes, pFile);

      for (i=0; i<nOffs; i++) {
        fprintf(stderr, "%i ", offs[i]);
      }
      fprintf(stderr, "\n");

      fread(&bytes, 1, sizeof(int), pFile);
      bytes -= 4;
      types = (char*) malloc(bytes);
      nType = (int)bytes/sizeof(char);
      fprintf(stderr, "Reading %i bytes (%i short int)\n", bytes, nType);
      fread(types, 1, bytes, pFile);

      for (i=0; i<nType; i++) {
        fprintf(stderr, "%i ", types[i]);
      }
      fprintf(stderr, "\n");

      exit(0);

    }


  } /* end of file read */

}
/* -------------------------------------------------------------------------- */
