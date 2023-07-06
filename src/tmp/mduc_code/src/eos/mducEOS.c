#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "mducEOS.h"
#include "parser.h"

void printHelp()
{
  printf("\n");
  printf("options :\n");
  printf("  -h                           # print help\n");
  printf("  -i <input file>              # input file\n");
  printf("  -o <output file>             # output filename\n");
  printf("commands :\n");
  printf("  --PT   <P>[bar]      <T>[K]  # calculate rho from P   T \n");
  printf("  --rhoT <rho>[kg/m^3] <T>[K]  # calculate P   from rho T \n");
  printf("  --post-process <file>    # post-process specified file\n");
  exit(0);
}
/* -------------------------------------------------------------------------- */

void setDefaultOptions(optStruct *opts)
{
  /* set default file names */
  strcpy(opts->inFile,  "eos.inp");
  strcpy(opts->outFile, "none");
}
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  int i,j,k;
  optStruct opts; /* struct to store options from command line */

  /* set options */
  setDefaultOptions(&opts);
  readOptions(argc, argv, &opts);

  parseStruct *parser = parseFile(opts.inFile);
  parsePrint(parser);

  /* initialize the fluid */
  char eos[32];

  paramGetString(eos,"equation-of-state","ideal",parser);
  eosFluid *f1 = fluidInit("dodecane",eos,parser);
  eosPrintFluid(f1);

  double isoT = paramGetReal("eos-temperature",300,parser);
  fprintf(stderr,"getting density range\n");
  double *rhoRange = NULL;
  int nRho = paramGetRealArray(&rhoRange,"density-range",parser);
  fprintf(stderr,"nRho = %i\n",nRho);
  fprintf(stderr,"density range: %i %f %f %f\n",nRho,rhoRange[0],rhoRange[1],
          rhoRange[2]);
  
  /* compute an isotherm */
  double rhoMin = rhoRange[0];
  double rhoMax = rhoRange[1];
  double dRho   = rhoRange[2];
  nRho = (rhoMax-rhoMin)/dRho + 1;

  /* double dRho   = (rhoMax - rhoMin)/(nRho-1); */

  double *rho   = (double*) calloc(nRho, sizeof(double));
  double *nu    = (double*) calloc(nRho, sizeof(double));
  double *P     = (double*) calloc(nRho, sizeof(double));
  /* double isoT = 500.0; */
  for (i=0; i<nRho; i++) {
    rho[i] = rhoMin + i*dRho;
    nu[i] = f1->Mw/rho[i]/1000.0;
  }

  for (i=0; i<nRho; i++)
    P[i]  = eosCalcP_rhoT(rho[i],isoT,f1);
                 
  FILE *pF;
  char outFile[128];
  sprintf(outFile,"%s_isotherm_T%.0f.out\0",eos,isoT);
  pF = fopen(outFile, "w");

  fprintf(pF,"# nu<1>[m^3/mol] rho<1>[kg/m^3] P<3>[Pa]\n");

  for (i=0; i<nRho; i++) {
    fprintf(pF,"%15.6e %15.6e %15.6e\n",nu[i],rho[i],P[i]);
  }

  fclose(pF);

  return 0;
}
/* -------------------------------------------------------------------------- */

void readOptions(int argc, char **argv, optStruct *opts)
{
  int i, j, k, nOpt, iOpt, iEnd;
  int *optLocs;
  
  /* read through the argument list and get the option locations */
  optLocs    = (int*)  calloc(argc, sizeof(int));
  opts->args = (char*) calloc(argc, sizeof(char));
  nOpt = 0;
  for (i=1; i<argc; i++) {
    if (argv[i][0] == '-') {
      optLocs[nOpt] = i;
      opts->args[nOpt] = argv[i][1];
      nOpt++;
    }
  } // end loop over arguments
  opts->nOpts = nOpt;
  
  for (i=0; i<nOpt; i++) {
    iOpt = optLocs[i];
    
    if (i == nOpt-1)
      iEnd = argc-1;
    else
      iEnd = optLocs[i+1]-1;
    
    switch (argv[optLocs[i]][1]) {
    case 'h':
      printHelp();
      break;
    case 'i':
      strcpy(opts->inFile, argv[iOpt+1]);
      break;
    case 'o':
      strcpy(opts->outFile, argv[iOpt+1]);
      break;
    case '-': { /* has a longer word */
      if (!strcmp(argv[optLocs[i]], "--post-process")) {
        /* opts->postProcess = 1; */
        /* strcpy(opts->postProcFile, argv[iOpt+1]); */
      }
      break;
    }
    } /* end switch */
  } /* end loop over options */
  
  free(optLocs); free(opts->args);
}
/* -------------------------------------------------------------------------- */
