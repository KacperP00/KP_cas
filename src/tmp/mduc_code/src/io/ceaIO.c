#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "ceaIO.h"
#include "chemKinIO.h"
#include "parser.h"
#include "utils.h"

int CEAreadThermo(char *fileName, chemStruct *chem)
{
  int i,j,k,iSp, iL, nSpecies;
  int maxAtom = 5, maxTrange = 9;
  double Tmin[9],Tmax[9];
  char line[128],buffer[128];
  char tmpS[16][128];
  char comment[8] = "!";

  /* temporary variables for reading in */
  char   species[24],atoms[5][3];
  double nu[5];

  FILE *pF = fopen(fileName, "r");
  if (pF == NULL) {
    fprintf(stderr, "Error: cannot read CEA file %s\n", fileName);
    exit(2);
  }

  /* set up initialization species */
  chem->nSpecies = 0;
  chem->species  = (char**) calloc(1, sizeof(char*));
  chem->specs    = (molStruct*)   calloc(1, sizeof(molStruct));
  chem->elems    = (atomStruct**) calloc(1, sizeof(atomStruct*));
  molStruct *sp;

  nSpecies = 0;
  while (!feof(pF) && fgets(line,128,pF)) { 
    /* check to see if it is a comment */
    stripString(line);
    if (!strncmp("END", line, 3))
      continue;
    if (strspn(line, comment) > 0 || !strcmp(line, "\n"))
      continue;

    if (strstr(line, "thermo") || strstr(line, "THERMO")) {
      /* read in the global temperature interval */
      fgets(line,128,pF);
      sscanf(line, "%10lf%10lf%10lf%10lf%40s",&chem->Trange[0],
             &chem->Trange[1],&chem->Trange[2],&chem->Trange[3],buffer);
    }
    else { /* read in a species record */
      nSpecies++;
      /* assign to chemistry struct */
      chem->nSpecies++;
      iSp = chem->nSpecies - 1;
      chem->species = (char**) realloc(chem->species,nSpecies*sizeof(char**));
      chem->species[iSp] = (char*) calloc(32, sizeof(char));
      chem->specs = (molStruct*) realloc(chem->specs,nSpecies*sizeof(molStruct));
      sp = &chem->specs[iSp];

      /* record 1: species name and comments */
      sscanf(line,"%24s%56s",chem->species[iSp],buffer);
      strcpy(sp->name, chem->species[iSp]);
      /* record 2: number of temperature intervals, id code,
                   formula symbols and numbers,
                   condensed, molecular weight, heat of formation */
      fgets(line,128,pF);
      sreadInt(&line[1],2,1,&sp->nTrange);
      if (sp->nTrange > maxTrange) {
        fprintf(stderr,"Error: cannot have %i temperature ranges\n",
                sp->nTrange);
        exit(2);
      }
      sreadString(&line[3],6,sp->info);
      for (i=0,j=10; i<5; i++,j+=8) {
        sreadString(&line[j],2,atoms[i]);
        sreadReal(&line[j+2],6,1,&nu[i]);
      }
      sreadInt(&line[51],1,1,&sp->cond);
      sreadReal(&line[52],13,1,&sp->W);
      sreadReal(&line[65],13,1,&sp->dHf);

      for (i=0; i<sp->nTrange; i++) {
        /* record 3: */
        fgets(line,128,pF);
        /* this is a hack, but seems to work */
        sreadReal(&line[1],10,1,&Tmin[i]);
        sreadReal(&line[11],10,1,&Tmax[i]);
        sp->Trange[i] = Tmin[i];
        sreadInt(&line[22],1,1,&sp->nCpCoeff[i]);
        sreadReal(&line[23],5,8,sp->Texp[i]);
        sreadReal(&line[67],15,1,&sp->H0[i]);
        /* record 4: first 5 coeffients */
        fgets(line,128,pF);
        /* HACK: for fortran double exponential format, replace D with E */
        for (j=0; j<strlen(line); j++) if (line[j] == 'D') line[j] = 'E';
        sreadReal(line,16,5,sp->therm[i]);
        /* record 5: rest coeffients */
        fgets(line,128,pF);
        /* HACK: for fortran double exponential format, replace D with E */
        for (j=0; j<strlen(line); j++) if (line[j] == 'D') line[j] = 'E';
        sreadReal(line,16,5,&sp->therm[i][5]);
      } /* loop Trange */
      sp->Trange[sp->nTrange] = Tmax[sp->nTrange-1];

      if (sp->nTrange == 0) { /* condensed species case */
        fgets(line,128,pF);
        sreadReal(&line[1],10,1,&sp->Trange[0]);
      }

      /* allocate the atoms for this species */
      sp->nAtoms = 0;
      for (i=0; i<maxAtom; i++)
        if (nu[i] > 0.0 && strlen(atoms[i]) > 0)
          sp->nAtoms++;
      sp->nA   = (int*) calloc(sp->nAtoms, sizeof(int));
      sp->atom = (atomStruct**) calloc(sp->nAtoms, sizeof(atomStruct*));
      for (i=0,j=0; i<maxAtom; i++) {
        if (nu[i] > 0.0 && strlen(atoms[i]) > 0) {
          sp->nA[j] = nu[i];
          k = getElementIndex(atoms[i],chem);
          if (k < 0) {
            fprintf(stderr, "Warning: element |%s| not found in periodic table for species %s\n",
                    atoms[i], sp->name);
            for (k=0; k<chem->nPT; k++)
              fprintf(stderr, "|%s|\n",chem->pt[k].symb);
            exit(2);
          }
          sp->atom[j] = &chem->pt[k];
          j++;
        }
      }
      /* ----- precompute the specific coefficients ----- */
      for (i=0; i<sp->nTrange; i++) {
        /* specific heat was read directly */
        for (j=0; j<10; j++)
          sp->aCp[i][j] = sp->therm[i][j];
        /* enthalpy */
        sp->aH[i][0] = -sp->therm[i][0];     /* T^-2  */
        sp->aH[i][1] =  sp->therm[i][1];     /* ln(T) */
        sp->aH[i][2] =  sp->therm[i][8];     /* -     */
        sp->aH[i][3] =  sp->therm[i][2];     /* T     */
        sp->aH[i][4] =  sp->therm[i][3]/2.0; /* T^2   */
        sp->aH[i][5] =  sp->therm[i][4]/3.0; /* T^3   */
        sp->aH[i][6] =  sp->therm[i][5]/4.0; /* T^4   */
        sp->aH[i][7] =  sp->therm[i][6]/5.0; /* T^5   */
        /* entropy */
        sp->aS[i][0] = -sp->therm[i][0]/2.0; /* T^-2  */
        sp->aS[i][1] = -sp->therm[i][1];     /* T^-1  */
        sp->aS[i][2] =  sp->therm[i][2];     /* ln(T) */
        sp->aS[i][3] =  sp->therm[i][9];     /* -     */
        sp->aS[i][4] =  sp->therm[i][3];     /* T     */
        sp->aS[i][5] =  sp->therm[i][4]/2.0; /* T^2   */
        sp->aS[i][6] =  sp->therm[i][5]/3.0; /* T^3   */
        sp->aS[i][7] =  sp->therm[i][6]/4.0; /* T^4   */
      }
    }
  } /* end file read */

  fclose(pF);
  return nSpecies;
}
/* -------------------------------------------------------------------------- */

int CEAwriteThermo(char *fileName, chemStruct *chem)
{
  int i,j,k;
  molStruct *sp;
  char atom[3];

  FILE *pF = fopen(fileName, "w");
  if (pF == NULL) {
    fprintf(stderr, "Error: cannot read CEA file %s\n", fileName);
    exit(2);
  }

  for (i=0; i<chem->nSpecies; i++) {
    sp = &chem->specs[i];
    fprintf(pF, "%-24s\n",sp->name);
    fprintf(pF, "%2i %6s ",sp->nTrange,sp->info);
    for (j=0; j<sp->nAtoms; j++) {
      strcpy(atom,sp->atom[j]->symb);
      stripString(atom);
      for (k=0; k<strlen(atom); k++) atom[k] = toupper(atom[k]);
      fprintf(pF,"%-2s%6.2f",atom,(double)sp->nA[j]);
    }
    for (j=0; j<5-sp->nAtoms; j++)
      fprintf(pF,"%-2s%6.2f","",0.0);
    fprintf(pF," %1i%13.5f%15.3f\n",sp->cond,sp->W,sp->dHf);
    for (j=0; j<sp->nTrange; j++) {
      fprintf(pF," %10.3f%10.3f %1i",sp->Trange[j],sp->Trange[j+1],sp->nCpCoeff[j]);
      for (k=0; k<8; k++)
        fprintf(pF,"%5.1f",sp->Texp[j][k]);
      fprintf(pF,"%17.3f\n",sp->H0[j]);
      for (k=0; k<5; k++)
        fprintf(pF,"%16.9E",sp->therm[j][k]);
      fprintf(pF,"\n");
      for (k=5; k<10; k++)
        fprintf(pF,"%16.9E",sp->therm[j][k]);
      fprintf(pF,"\n");
    }
  } /* end loop over species */

  fclose(pF);
  return 0;
}
/* -------------------------------------------------------------------------- */
