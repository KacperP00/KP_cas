#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "display.h"
#include "utils.h"

void mducHeader(char *head)
{
  int w = 80; /* default table width */
  int len = strlen(head);
  fprintf(stderr, "");
}
/* -------------------------------------------------------------------------- */

void mducInfo(char *message, char *func)
{
  int len;
  if (MDUC_PRINT_PROC) {
    len = strlen(message) + strlen(func);
    fprintf(stderr, " %s: ", func);
    if (len > 80-3) fprintf(stderr, "\n%8s","");
    fprintf(stderr, "%s\n", message);
  }
}
/* -------------------------------------------------------------------------- */

void mducErr(char *message, char *func)
{
  int len;
  if (MDUC_PRINT_PROC) {
    len = strlen(message) + strlen(func);
    fprintf(stderr, "***** MDUC Error:   %s: ", func);
    if (len > 80-24) fprintf(stderr, "\n%8s","");
    fprintf(stderr, "%s\n", message);
    exit(2);
  }
}
/* -------------------------------------------------------------------------- */

void mducWarn(char *message, char *func)
{
  int len;
  if (MDUC_PRINT_PROC) {
    len = strlen(message) + strlen(func);
    fprintf(stderr, "***** MDUC Warning: %s: ", func);
    if (len > 80-24) fprintf(stderr, "\n%8s","");
    fprintf(stderr, "%s\n", message);
  }
}
/* -------------------------------------------------------------------------- */

monStruct * monitorInit(char *id)
{
  if (MDUC_PRINT_PROC)
    fprintf(stderr, " Initializing monitor: %s\n", id);
  monStruct *mon = (monStruct*) calloc(1, sizeof(monStruct));
  strcpy(mon->id, id);
  char fileName[256];
  sprintf(fileName, "%s.out", id);
  mon->p = fopen(fileName, "w");
  mon->n     = 0;
  mon->type  = (int*)   calloc(1, sizeof(int));
  mon->names = (char**) calloc(1, sizeof(char*));
  mon->units = (char**) calloc(1, sizeof(char*));
  mon->vals  = (void**) calloc(1, sizeof(void*));
  return mon;
}
/* -------------------------------------------------------------------------- */

void monitorRegister(char *name,char *unit,void *val,int type,monStruct *mon)
{
  int i;
  mon->n++;
  i = mon->n - 1;
  mon->type  = (int*)   realloc(mon->type,  mon->n*sizeof(int));
  mon->names = (char**) realloc(mon->names, mon->n*sizeof(char*));
  mon->units = (char**) realloc(mon->units, mon->n*sizeof(char*));
  mon->vals  = (void**) realloc(mon->vals,  mon->n*sizeof(void*));
  mon->type[i] = type;
  mon->names[i] = (char*) calloc(strlen(name)+8, sizeof(char));
  mon->units[i] = (char*) calloc(strlen(unit)+8, sizeof(char));
  strcpy(mon->names[i],name);
  strcpy(mon->units[i],unit);
  mon->vals[i] = val;
}
/* -------------------------------------------------------------------------- */

void monitorWriteHead(monStruct *mon)
{
  int i;
  fprintf(mon->p,"# ");
  for (i=0; i<mon->n; i++) {
    if (strlen(mon->units[i]) > 0)
      fprintf(mon->p,"%s[%s]<%i> ",mon->names[i],mon->units[i],i+1);
    else
      fprintf(mon->p,"%s<%i> ",mon->names[i],i+1);
  }
  fprintf(mon->p,"\n");
  fflush(mon->p);
}
/* -------------------------------------------------------------------------- */

void monitorWriteLine(monStruct *mon)
{
  int i;
  for (i=0; i<mon->n; i++) {
    if (mon->type[i] == 0)
      fprintf(mon->p, "%6i ", *(int*)mon->vals[i]);
    else if (mon->type[i] < 0)
      fprintf(mon->p, "%15.6e ", *(double*)mon->vals[i]);
    else
      fprintf(mon->p, "%10.*f ", mon->type[i], *(double*)mon->vals[i]);
  }
  fprintf(mon->p,"\n");
  fflush(mon->p);
}
/* -------------------------------------------------------------------------- */

void destroyMonitor(monStruct *mon)
{
  int i;
  for (i=0; i<mon->n; i++) {
    if (mon->names[i]) free(mon->names[i]);
    if (mon->units[i]) free(mon->units[i]);
  }
  if (mon->type)  free(mon->type);
  if (mon->names) free(mon->names);
  if (mon->units) free(mon->units);
  if (mon->vals)  free(mon->vals);
  fclose(mon->p);
}
/* -------------------------------------------------------------------------- */

void printInput(ioStruct *io, driveStruct *driver, chemStruct *chem,
                timeStruct *time, flameStruct *flame)
{
  int i, j;
  
  printf("\n\n----- Data read from input file -----\n\n");
  
  printf("DRIVERS:\n");
  printf("pressure-mode:      %i\n", driver->PRES_MODE);
  printf("mixing-mode:        %i\n", driver->MIX_MODE);
  printf("chi-mode:           %i\n", driver->CHI_MODE);
  printf("diffusion-mode:     { %g %g %g }\n", 
	 driver->DIFF_MODE[0], driver->DIFF_MODE[1], driver->DIFF_MODE[2]);
  printf("convection-mode:    %g\n", driver->CONV_MODE);
  printf("time-mode:          %i\n", driver->TIME_MODE);
  printf("enthalpy-dimension: %i\n", driver->H_DIM);
  
  printf("Time Data:\n");
  printf("start-time:     %f\n", time->t0);
  printf("final-time:     %f\n", time->t_end);
  printf("iterations:     %i\n", time->nIter);
  printf("time-step:      %i\n", time->nt);
  printf("time-step-size: %f\n", time->dt);
  
  printf("IO data:\n");
  printf("write-data-file: %s\n", io->outputData);
  printf("read-data-file:  %s\n", io->inputData);
  
  printf("Flamelet Data:\n");
  printf("gridPoints:  %5i\n", flame->n1);
  printf("number-species: %5i\n", chem->nSpecies);
  printf("---------------------\n");
  printf("%5s %10s","x1","temp [k]");
  for (i=0; i < io->nInputs; i++)
    printf("%10s ", io->inputSpecies[i]);
  printf("\n");
  for (j=0; j<flame->n1; j++) 
    {
      printf("%5.3f %10.3f ", flame->x1[j], io->initT[j]);
      for (i=0; i < io->nInputs; i++)
        printf("%10.3e ", io->initY[i][j]);
      printf("\n");
    }
}
/* -------------------------------------------------------------------------- */

void printArray(int n, int index, int offset, double *y)
{
  int i, k;
  for (i=0; i<n; ++i)
    {
      k = i*offset;
      printf("%3i %13.3e\n", k+index, y[k+index]);
    }
}
/* -------------------------------------------------------------------------- */


void printSpeciesList(chemStruct *chem)
{
  int i;
  
  printf("\n----- Species in Mechanism -----\n");
  for (i=0; i<chem->nSpecies; ++i)
    {
      printf("%4i %25s\n", i+1, chem->species[i]);
    }
  printf("--------------------------------\n");
}
/* -------------------------------------------------------------------------- */

void printStreamProperties(ioStruct *io)
{	
  int i;
  
  printf("\n");
  printf("----- Mixture properties of the stream -----\n\n");
  printf("%10s %10s %10s %10s\n", "", "component", "Y", "Temp");
  printf("--------------------------------------------\n");
  for(i=0; i<io->nOx; ++i)
    {
      if (strcmp(io->oxName[i],"Temp"))
        {
          if (i == 0)
            printf("%10s %10s %10.3f %10.1f\n", "Oxidizer:",
                   io->oxName[i], io->oxY[i], io->oxT);
          else
            printf("%10s %10s %10.3f\n", "", io->oxName[i], io->oxY[i]);
        }
    }
  for(i=0; i<io->nFuel1; ++i)
    {
      if (strcmp(io->fuelName1[i],"Temp"))
        {
          if (i == 0)
            printf("%10s %10s %10.3f %10.1f\n", "Fuel 1:",
                   io->fuelName1[i], io->fuelY1[i], io->fuelT1);
          else
            printf("%10s %10s %10.3f\n", "", io->fuelName1[i], io->fuelY1[i]);
        }
    }
  for(i=0; i<io->nFuel2; ++i)
    {
      if (strcmp(io->fuelName2[i],"Temp"))
        {
          if (i == 0)
            printf("%10s %10s %10.3f %10.1f\n", "Fuel 2:",
                   io->fuelName2[i], io->fuelY2[i], io->fuelT2);
          else
            printf("%10s %10s %10.3f\n", "", io->fuelName2[i], io->fuelY2[i]);
        }
    }
  printf("--------------------------------------------\n");
}
/* -------------------------------------------------------------------------- */

void printMonitor(long int nst, double t, double dtPre, ioStruct *io, double *Y,
                  flameStruct *flame, double dt_max)
{
  int i, j, ind, jMaxT;
  double Tmax;

  if (flame->T)
    Tmax = getMaxValue(flame->n1, flame->n2, 0, 1, flame->T);
  else
    Tmax = getMaxValue(flame->n1, flame->n2, flame->sT, flame->nVars, flame->Y);
    /* jMaxT = getIndexMax(flame->n1,flame->nVars,&Y[flame->sT]); */
  // print the progress to the screen along with any monitor species
  fprintf(io->monitor, "%6li %15.6e %15.6e %15.6e %15.6e\n", nst, t, dtPre, Tmax,dt_max);

  if (io->nMonitorY > 0) {
    /* print table header */
    fprintf(io->monitor, "%5s %8s %8s", "Z", "T", "chiA");
    if (flame->nDim == 2)
      fprintf(io->monitor, "%8s %8s","chiB","chiC");
    
    for (i=0; i< io->nMonitorY; ++i)
      fprintf(io->monitor, " %13s", io->monitorSpecies[i]);
    fprintf(io->monitor, "\n");
    
    /* print table */
    for (j=0; j < io->nMonitorZ; ++j) {
      ind = floor((flame->n1-1)*j/(io->nMonitorZ-1)+0.49);
      if (flame->T)
        fprintf(io->monitor, "%5.3f %8.1f %8.3g", flame->x1[ind],
                flame->T[ind], flame->chiA[ind]);
      else
        fprintf(io->monitor, "%5.3f %8.1f %8.3g", flame->x1[ind],
                Y[ind*flame->nVars+flame->sT], flame->chiA[ind]);

      if (flame->nDim == 2)
        fprintf(io->monitor, "%8.3g %8.3g", flame->chiB[ind],
                flame->chiC[ind]);
      
      for (i=0; i < io->nMonitorY; ++i)
        fprintf(io->monitor, " %13.3e",
                Y[ind*flame->nVars + io->monitorIndex[i]]);
      
      fprintf(io->monitor, "\n");
    }
    fprintf(io->monitor, "\n");
  } /* end table */
  fflush(io->monitor);
}
/* -------------------------------------------------------------------------- */
