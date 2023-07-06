#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include "ensightIO.h"
#include "fileIO.h"
#include "parser.h"
#include "accessChem.h" // should remove dependence on this
#include "utils.h"

void initEnsight(ioStruct *io)
{
  int i;
  struct stat st;
  char buffer[100];

  // check to see if the directory exists
  if (stat(io->ensightDir, &st) != 0) mkdir(io->ensightDir, 0755);

  // species to plot
  for (i=0; i<io->nInputs; i++)
    {
      strcpy(buffer, io->ensightDir);
      strcat(buffer, "/");
      strcat(buffer, io->inputSpecies[i]);
      stripString(buffer);
      
      if (stat(buffer, &st) != 0) mkdir(buffer, 0755);
    }

  // temperature
  strcpy(buffer, io->ensightDir);
  strcat(buffer, "/T");
  stripString(buffer);
  if (stat(buffer, &st) != 0) mkdir(buffer, 0755);
}
/* -------------------------------------------------------------------------- */

void writeEnsightCaseHead(ioStruct *io, driveStruct *driver,
                          timeStruct *time, flameStruct *flame, flowStruct *flow)
{
  int i;
  char  caseFile[100], buffer[100];

  strcpy(caseFile, io->ensightDir);
  strcat(caseFile, "/mduc.case");

  FILE *pFile = fopen(caseFile, "w");

  fprintf(pFile, "FORMAT\n");
  fprintf(pFile, "type: ensight gold\n");
  fprintf(pFile, "GEOMETRY\n");
  fprintf(pFile, "model: geometry\n");
  fprintf(pFile, "VARIABLE\n");

  for (i=0; i<io->nInputs; i++)
    fprintf(pFile, "scalar per node: 1 %s %s/%s.******\n",
            io->inputSpecies[i], io->inputSpecies[i], io->inputSpecies[i]);

  fprintf(pFile, "scalar per node: 1 Temperature T/T.******\n");
  
  fprintf(pFile, "TIME\n");
  fprintf(pFile, "time set: 1\n");
  fprintf(pFile, "number of steps: %i\n", 0);
  fprintf(pFile, "filename start number: %i\n", 1);
  fprintf(pFile, "filename increment: %i\n", 1);
  fprintf(pFile, "time values: ");
  fprintf(pFile, "\n");

  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

void updateEnsightCaseTime(double newTime, ioStruct *io, timeStruct *time)
{
  int i, nt, iStart, increment;
  long int headSize;
  char  caseFile[100], tag[100], buffer[100], *headBuffer, *ident;
  double *timeVals;

  strcpy(caseFile, io->ensightDir);
  strcat(caseFile, "/mduc.case");

  FILE *pFile = fopen(caseFile, "r");

  // find the beginning of the TIME section
  while (!feof(pFile))
    { 
      fscanf(pFile, "%s", tag);  // read in string

      if (!strcmp(tag, "TIME"))
        break;
    }

  // get the total length of the header/body in bytes
  headSize = ftell(pFile) - 4;

  // read in the time info
  fgets(buffer, 80, pFile); // time set line
  fgets(buffer, 80, pFile); // time set
  fgets(buffer, 80, pFile); // number of steps
  ident = strtok(buffer, ":");
  nt = atoi(strtok(NULL,""));
  fgets(buffer, 80, pFile); // file start number
  ident = strtok(buffer, ":");
  iStart = atoi(strtok(NULL,""));
  fgets(buffer, 80, pFile); // file increment
  ident = strtok(buffer, ":");
  increment = atoi(strtok(NULL,""));

 /*  // "number of steps:" */
/*   fgets(buffer, 17, pFile); */
/*   nt = readInteger(pFile); */

/*   // "filename start number:" */
/*   fgets(buffer, 23, pFile); */
/*   iStart = readInteger(pFile); */

/*   // "filename increment:" */
/*   fgets(buffer, 20, pFile); */
/*   increment = readInteger(pFile); */

  // "time values:"
  fgets(buffer, 12, pFile);
  

  // find the beginning of the TIME section
  /* fscanf(pFile, "%s", tag); */
/*   do { */
/*     if (!strncmp(tag, "steps", 5)) */
/*       nt = readInteger(pFile); */
/*     else if (!strncmp(tag, "number", 6)) */
/*       iStart = readInteger(pFile); */
/*     else if (!strncmp(tag, "increment", 9)) */
/*       increment = readInteger(pFile); */
      
/*     fscanf(pFile, "%s", tag);  // read in string */
/*   } while (!strncmp(tag, "values", 6)); */

  // read in already existing time steps
  readRealArray(pFile, nt, &timeVals);
  
  // go back to the beginning of the file and read in the whole header
  rewind(pFile);
  headBuffer = (char*) calloc(headSize, sizeof(char));
  fread(headBuffer, 1, headSize, pFile);

  freopen(caseFile, "w", pFile);

  // re-write the header
  fwrite(headBuffer, 1, headSize, pFile);

  // update the time vectors
  nt++;
  timeVals = (double*) realloc(timeVals, nt*sizeof(double));
  timeVals[nt-1] = newTime;

  // rewrite the time info
  fprintf(pFile, "TIME\n");
  fprintf(pFile, "time set: 1\n");
  fprintf(pFile, "number of steps: %i\n", nt);
  fprintf(pFile, "filename start number: %i\n", iStart);
  fprintf(pFile, "filename increment: %i\n", increment);
  fprintf(pFile, "time values:\n");
  writeArray(pFile, 1, nt, 1, 0, 5, timeVals);

  fclose(pFile);

  free(timeVals);
  free(headBuffer);
}
/* -------------------------------------------------------------------------- */

void writeEnsightData(char *fileName, char *varName, int m, int n,
                      int offset, int index, double *data)
{
  int  iTmp = 1;
  char buffer[80];

  int j, k, row, col, curLine;
  float *tmp = (float*) calloc(m*n, sizeof(float));
  
  for (k = 0; k < m; ++k)
    {
      row = k*n*offset;    // set the row index
      for (j = 0; j < n; ++j)
        {
          col = j*(offset); // set the column index
          tmp[k*n + j] = (float)(data[row+col+index]);
        }
    }

  FILE *pFile = fopen(fileName, "wb");

  sprintf(buffer, "Variable file: %s", varName);
  fprintf(pFile, "%-80s", buffer);
  fprintf(pFile, "%-80s", "part");
  fwrite(&iTmp, sizeof(int), 1, pFile);
  fprintf(pFile, "%-80s", "block");
  fwrite(tmp, sizeof(float), n*m, pFile);

  fclose(pFile);

  free(tmp);
}
/* -------------------------------------------------------------------------- */

void writeEnsightGeom(ioStruct *io, driveStruct *driver,
                      timeStruct *time, flameStruct *flame, flowStruct *flow)
{
  int   iTmp = 1;
  char  geomFile[100];
  float fTmp = 0.0;

  strcpy(geomFile, io->ensightDir);
  strcat(geomFile, "/geometry");

  float xMin, xMax, yMin, yMax, zMin, zMax;
  xMin = 0.0;
  xMax = 1.0;
  yMin = 0.0;
  if (flame->nDim == 2)
    yMax = 1.0;
  else
    yMax = 0.0;
  zMin = 0.0;
  zMax = 0.0;

  // write geometry
  FILE *pFile = fopen(geomFile, "wb");

  fprintf(pFile, "%-80s", "C BINARY");
  fprintf(pFile, "%-80s", "Ensight Geometry File");
  fprintf(pFile, "%-80s", "Structured Geometry from MDUC");
  fprintf(pFile, "%-80s", "node id assign");
  fprintf(pFile, "%-80s", "element id assign");
  fprintf(pFile, "%-80s", "extents");
  fwrite(&xMin, sizeof(float), 1, pFile);
  fwrite(&xMax, sizeof(float), 1, pFile);
  fwrite(&yMin, sizeof(float), 1, pFile);
  fwrite(&yMax, sizeof(float), 1, pFile);
  fwrite(&zMin, sizeof(float), 1, pFile);
  fwrite(&zMax, sizeof(float), 1, pFile);
  fprintf(pFile, "%-80s", "part");
  fwrite(&iTmp, sizeof(int), 1, pFile);
  fprintf(pFile, "%-80s", "MDUC Grid");
  fprintf(pFile, "%-80s", "block rectilinear");
  fwrite(&flame->n1, sizeof(int), 1, pFile);
  fwrite(&flame->n2, sizeof(int), 1, pFile);
  fwrite(&iTmp, sizeof(int), 1, pFile);
  writeBinaryArrayReal(pFile, flame->n1, flame->x1);
  if (flame->nDim == 2)
    writeBinaryArrayReal(pFile, flame->n2, flame->x2);
  else
    fwrite(&fTmp, sizeof(float), 1, pFile);
  fwrite(&fTmp, sizeof(float), 1, pFile);

  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

void writeBinaryArrayReal(FILE *pFile, int n, double *val)
{
  int i;
  float *tmp = (float*) calloc(n, sizeof(float));
  
  for (i=0; i<n; i++) tmp[i] = (float)(val[i]);
  fwrite(tmp, sizeof(float), n, pFile);

  free(tmp);
}
/* -------------------------------------------------------------------------- */

/* void writeBinaryArrayReal(FILE *pFile, int m, int n, int offset, int index, */
/*                           double *data) */
/* { */
/*   int j, k, row, col, curLine; */
/*   float *tmp = (float*) calloc(n*m, sizeof(float)); */
  
/*   for (k = 0; k < m; ++k) */
/*     { */
/*       row = k*n*offset;    // set the row index */
/*       for (j = 0; j < n; ++j) */
/*         { */
/*           col = j*(offset); // set the column index */
/*           tmp[k*n + j] = (float)(data[row+col+index]); */
/*         } */
/*     } */
/*   fwrite(tmp, sizeof(float), n*m, pFile); */

/*   free(tmp); */
/* } */
/* -------------------------------------------------------------------------- */
