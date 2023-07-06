#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "chemKinIO.h"
#include "chemUtils.h"
#include "utils.h"
#include "parser.h"
#include "display.h"

int getChemKinSpecies(char *fileName, char ***species)
{
  int nSpecies;
  FILE *pF = NULL;
  char tag[32];
 
  pF = fopen(fileName, "r");
  if (pF == NULL) {
    fprintf(stderr, "Error: cannot read ChemKin file %s\n", fileName);
    exit(2);
  }

  /* provide initial array allocation */
  species[0] = (char**) calloc(1, sizeof(char*));

  nSpecies = 0;
  while (!feof(pF)) {
    fscanf(pF, "%s", tag);

    /* read species found in file */
    if (!strcmp(tag, "SPECIES")) {
      fscanf(pF, "%s", tag);
      if (!strcmp(tag, "END")) {
        fprintf(stderr, "Error: no species in ChemKin file %s\n", fileName);
        exit(2);
      }

      /* read in species until "END" keyword found */
      do {
        nSpecies++;
        species[0] = (char**) realloc(species[0], nSpecies*sizeof(char*));
        species[0][nSpecies-1] = (char*) calloc(25, sizeof(char));
        strncpy(species[0][nSpecies-1], tag, 25);
        fscanf(pF, "%s", tag);
      } while (!strcmp(tag, "END") || feof(pF));

      fclose(pF);
      return nSpecies;
    }
  } /* end file read */
  fprintf(stderr, "Error: no species found in ChemKin file %s\n", fileName);
  exit(2);
}
/* -------------------------------------------------------------------------- */

int readChemKinMech(char *fileName, int *nAtoms, char ***atoms, int *nSpecs,
                    char ***species, rxnStruct **rxns, chemStruct *chem)
{
  int i, j, iSp, iL, iR, n, nRxn, nSpecies, nElem;
  FILE *pF = NULL;
  fpos_t pos;
  char tag[32], rxnLn[128], buffer[128], buffSp[128];
  char *pReact, *pProd, *pArr, *pS, *pM;
  char *rM;

  rxnStruct *r;

  pF = fopen(fileName, "r");
  if (pF == NULL) {
    fprintf(stderr, "Error: cannot read ChemKin file %s\n", fileName);
    exit(2);
  }

  /* pre-allocate some memory */
  atoms[0]   = (char**)     calloc(1, sizeof(char*));
  species[0] = (char**)     calloc(1, sizeof(char*));
  rxns[0]    = (rxnStruct*) calloc(1, sizeof(rxnStruct));

  nElem = nSpecies = nRxn = 0;
  while (!feof(pF)) {
    fscanf(pF, "%s\n", tag);

    if (tag[0] == '!') { /* comments */
      fgets(buffer, 128, pF);
    }

    /* ----- read in the elements ----- */
    else if (!strcmp(tag, "ELEMENTS")) {
      fscanf(pF, "%s\n", tag);
      while (strcmp(tag, "END")) {
        nElem++;
        iL = nElem-1;
        atoms[0]     = (char**) realloc(atoms[0], nElem*sizeof(char*));
        atoms[0][iL] = (char*)  calloc(3, sizeof(char));

        strcpy(atoms[0][iL], tag);
        fscanf(pF, "%s\n", tag);
      }
      *nAtoms = nElem;
    } /* end ELEMENTS */

    /* ----- read in the species ----- */
    else if (!strcmp(tag, "SPECIES")) {
      fscanf(pF, "%s\n", tag);
      while (strcmp(tag, "END")) {
        nSpecies++;
        iL = nSpecies-1;
        species[0]     = (char**) realloc(species[0], nSpecies*sizeof(char*));
        species[0][iL] = (char*)  calloc(32, sizeof(char));

        strcpy(species[0][iL], tag);
        fscanf(pF, "%s\n", tag);
      }
      *nSpecs = nSpecies;
    } /* end SPECIES */

    /* ----- read in the reactions ----- */
    else if (!strcmp(tag, "REACTIONS")) {
      /* read in the line */
      fgets(rxnLn, 128, pF);
      if (strstr(rxnLn,"MOL") != NULL) { /* read in units */
        stripString(rxnLn);
        n = strcspn(rxnLn," ");
        strncpy(chem->unitsE,rxnLn,n);
        chem->unitsE[n] = '\0';
        if (strlen(rxnLn) > n) {
          strcpy(chem->unitsA,&rxnLn[n]);
          stripString(chem->unitsA);
        }
        /* get next line */
        fgets(rxnLn, 128, pF);
      }

      while (!strstr(rxnLn, "END")) {
        /* don't operate on original line to prevent corruption */
        strcpy(buffer, rxnLn);

        /* check whether this is a reaction or not */
        rM = strchr(buffer, '=');
        if (rM) {
          nRxn++;
          rxns[0] = (rxnStruct*) realloc(rxns[0], nRxn*sizeof(rxnStruct));
          r = &rxns[0][nRxn-1]; /* tmp for readability */
          r->no = nRxn+1;

          /* ----- set defaults ----- */
          setRxnDefault(r);

          /* check whether this is a fall-off reaction */
          /* TODO: this assumes the there is no space between ( and + */
          pM = strstr(buffer,"(+");
          while (pM != NULL) {
            r->pressure = 2;
            pM[0] = ' ';
            pS = strstr(pM, ")");
            pS[0] = ' ';
            pM = strstr(pS,"(+");
          }

          /* ----- parse the reactants ----- */
          /* find the end of the reactants */
          n = strcspn(buffer, "<=");
          strncpy(buffSp, buffer, n);
          buffSp[n] = '\0';
          r->reverse = 0;
          if (buffer[n] == '<')
            r->reverse = 1;
          r->nR = parseReactionSpecies(buffSp,"<",&r->nuR,&r->spR);

          /* check whether this is a third body reaction */
          for (i=0, r->thirdBody=0; i<r->nR; i++) {
            if (!strcmp(r->spR[i], "M")) {
              r->thirdBody = 1;
              r->pressure  = 1;
            }
          }

          /* ----- parse the products ----- */
          /* find the end of the products */
          n = i = 0;
          do {
            i = strcspn(&rM[n], "+");
            if (i < strlen(rM)-n) {
              if (rM[n+i-1] == 'E' || rM[n+i-1] == 'e' && isdigit(rM[n+i+1])) {
                break;
              }
              n += i+1;
            }
            else {
              break;
            }
          } while (1 == 1);
          i = strcspn(&rM[n], " \t");
          n += i;
          j = 0;
          while (rM[j] == '=' || rM[j] == '>') j++;
          strncpy(buffSp, &rM[j], n);
          buffSp[n] = '\0';
          r->nP = parseReactionSpecies(buffSp,"{",&r->nuP,&r->spP);

          /* ----- read in the rate coefficients ----- */
          rM += n;
          sscanf(rM, "%lf%lf%lf", &r->A, &r->b, &r->E);
        }
        else { /* check for auxiliary input */
          /* NOTE: reaction struct should be from previous reaction read */
          if (r == NULL) {
            fprintf(stderr, "Error: reaction not defined for auxiliary: ");
            fprintf(stderr, " %s\n", rxnLn);
            exit(2);
          }
          strcpy(buffer, rxnLn);
          if (pM = strstr(buffer, "LOW")) {
            r->pressure = 2;
            r->lindemann = 1;
            pS = strtok(buffer, "/");
            pS = strtok(NULL, "/");
            sscanf(pS, "%lf%lf%lf", &r->A_low, &r->b_low, &r->E_low);
          }
          else if (pM = strstr(buffer, "REV")) {
            r->reverse = 3;
            pS = strtok(buffer, "/");
            pS = strtok(NULL, "/");
            sscanf(pS, "%lf%lf%lf", &r->Ar, &r->br, &r->Er);
          }
          else if (pM = strstr(buffer, "TROE")) {
            r->pressure = 2;
            r->troe = 1;
            pS = strtok(buffer, "/");
            pS = strtok(NULL, "/");
            strcpy(buffSp, pS);
            strcpy(buffer, buffSp);
            /* count the number of entries */
            n = 0;
            pS = strtok(buffer, " ");
            while (pS) {
              n++;
              pS = strtok(NULL, " ");
            }
            if (n == 4) {
              sscanf(buffSp, "%lf%lf%lf%lf", &r->troeA[1],&r->troeT[0],
                     &r->troeT[1],&r->troeT[2]);
            }
            else {
              sscanf(buffSp,"%lf%lf%lf",&r->troeA[1],&r->troeT[0],&r->troeT[1]);
              r->troeT[2] = -1.0;
            }
            r->troeA[0] = 1.0 - r->troeA[1];
          }
          else if (pM = strstr(buffer, "SRI")) {
            r->pressure = 2;
            r->sri = 1;
            pS = strtok(buffer, "/");
            pS = strtok(NULL, "/");
            strcpy(buffSp, pS);
            strcpy(buffer, buffSp);
            /* count the number of entries */
            n = 0;
            pS = strtok(buffer, " ");
            while (pS) {
              n++;
              pS = strtok(NULL, " ");
            }
            /* set defaults */
            r->sriA[3] = 1.0; /* d */
            r->sriA[4] = 0.0; /* e */
            if (n == 5)
              sscanf(buffSp,"%lf%lf%lf%lf%lf",&r->sriA[0],&r->sriA[1],
                     &r->sriA[2],&r->sriA[3],&r->sriA[4]);
            else if (n == 4)
              sscanf(buffSp, "%lf%lf%lf%lf",&r->sriA[0],&r->sriA[1],&r->sriA[2],
                     &r->sriA[3]);
            else
              sscanf(buffSp, "%lf%lf%lf", &r->sriA[0],&r->sriA[1],&r->sriA[2]);
          }
          else { /* auxiliary information */
            r->nEnhance = 0;
            r->spEnh    = (char**) calloc(1, sizeof(char*));
            r->factEnh  = (double*) calloc(1, sizeof(double));
            strcpy(buffer, rxnLn);
            stripString(buffer);
            pS = strtok(buffer, "/");
            while (pS && !strchr(pS, '\n')) {
              r->nEnhance++;
              i = r->nEnhance - 1;
              r->spEnh = (char**) realloc(r->spEnh, r->nEnhance*sizeof(char*));
              r->spEnh[i] = (char*) calloc(32, sizeof(char));
              r->factEnh = (double*) realloc(r->factEnh,
                                             r->nEnhance*sizeof(double));

              /* copy the species */
              j = 0;
              while (pS[j] == ' ') pS++;
              strcpy(r->spEnh[i], pS);

              /* get the factor */
              pS = strtok(NULL, "/");
              r->factEnh[i] = atof(pS);

              pS = strtok(NULL, "/");
            };
          }
        } /* end AUXILIARY input */

        fgets(rxnLn, 128, pF);
      }
    } /* end REACTIONS */
    else { /* blank line */

    }

  } /* end FILE */
  fclose(pF);
  
  return nRxn;
}
/* -------------------------------------------------------------------------- */

void writeChemKinMech(char *fileName, chemStruct *chem)
{
  int i, j, spMax, spLine;
  double tmp, rLine;
  char  buff[32];
  FILE *pF = NULL;

  pF = fopen(fileName, "w");
  if (!pF) {
    fprintf(stderr, "Error: cannot open file %s\n", fileName);
    exit(2);
  }

  fprintf(stderr, "Writing ChemKin mechanism: %s\n", fileName);
  /* ----- write the elements ----- */
  fprintf(pF, "ELEMENTS\n");
  for (i=0; i<chem->nElem; i++) {
    strcpy(buff, chem->elems[i]->symb);
    for (j=0; j<strlen(buff); j++)
      buff[j] = toupper(buff[j]);
    fprintf(pF, "%-2s ", buff);
  }
  fprintf(pF, "\nEND\n");

  /* ----- write the species ----- */
  spMax = 0;
  /* find the maximum size of the species */
  for (i=0; i<chem->nSpecies; i++) {
    if (strlen(chem->species[i]) > spMax)
      spMax = strlen(chem->species[i]);
  }
  tmp = modf(80.0/(spMax+1), &rLine);
  spLine = (int)(rLine);
  fprintf(pF, "SPECIES\n");
  for (i=0; i<chem->nSpecies; i++) {
    fprintf(pF, "%-*s ", spMax, chem->species[i]);
    if ((i+1)%spLine == 0)
      fprintf(pF, "\n");
  }
  fprintf(pF, "\nEND\n");

  /* ----- write the reactions ----- */
  fprintf(pF, "REACTIONS KJOULES/MOLE\n");
  for (i=0; i<chem->nRxn; i++) {
    writeChemKinRxn(pF, chem->rxns[i]);
  }
  fprintf(pF, "END\n");
  fclose(pF);

}
/* -------------------------------------------------------------------------- */

void writeChemKinRxn(FILE *pF, rxnStruct r)
{
  int i, n, nRxn = 0, nM;
  /* conversion to standard units */
  double Ecnv = 1.0e-6; /* J/kmol -> kJ/mol */
  double Acnv = pow(1.0e3,r.orderF[1]-1);  /* m^(3(n-1)) -> cm^(3(n-1) */
  char buffer[256];
  /* ----- write the reactants/products ----- */
  nRxn += writeChemKinRxnSpecies(pF, r.nR, r.nuR, r.spR, r.pressure);
  if (r.reverse) {
    fprintf(pF, "<");
    nRxn += 1;
  }
  fprintf(pF, "=>");
  nRxn += 2;
  nRxn += writeChemKinRxnSpecies(pF, r.nP, r.nuP, r.spP, r.pressure);

  /* ----- write the rate coefficients ----- */
  fprintf(pF,"%*s%10.3E%10.3f %10.2f\n",50-nRxn-1,"",r.A*Acnv,r.b,r.E*Ecnv);
  
  /* ----- write out any additional information ----- */
  if (r.reverse == 3) {
    Acnv = pow(1.0e3,r.orderR[1]-1);
    fprintf(pF,"%2sREV/%12.3E%10.3f%10.3f/\n"," ",r.Ar*Acnv,r.br,r.Er*Ecnv);
  }
  if (r.lindemann == 1 || r.troe == 1) {
    Acnv = pow(1.0e3,r.orderF[0]-1);
    fprintf(pF,"%2sLOW/%12.3E%10.3f%10.3f/\n"," ",r.A_low*Acnv,r.b_low,
            r.E_low*Ecnv);
  }
  if (r.troe == 1) {
    fprintf(pF,"%2sTROE/%10.4f%10.2f%10.3f", " ",
            r.troeA[1], r.troeT[0], r.troeT[1]);
    if (r.troeT[2] > 0.0)
      fprintf(pF, "%10.2f", r.troeT[2]);
    fprintf(pF, "/\n");
  }

  /* ----- write the enhancement factors ----- */
  if (r.nEnhance > 0) {
    for (i=0, nM=0; i<r.nEnhance; i++) {
      if (r.factEnh[i] != 1.0 && strcmp(r.spEnh[i], "OTHER")) {
        if (nM == 0)
          fprintf(pF, "%2s", " ");
        fprintf(pF, "%s/%.2f/ ", r.spEnh[i], r.factEnh[i]);
        nM++;
      }
    }
    if (nM > 0)
      fprintf(pF, "\n");
  }
}
/* -------------------------------------------------------------------------- */

int writeChemKinRxnSpecies(FILE *pF, int nSp, double *nu, char **sp, int fo)
{
  int i, len = 0;
  char buf[64];

  /* start length of the final reaction string */
  /* TODO: assumes no spaces between species */
  len += nSp - 1;       /* number of '+' required */
  if (fo > 1) len += 2; /* brackets to indicate fall-off reaction */
  
  for (i=0; i<nSp; i++) {
    if (nu[i] > 1) {
      sprintf(buf, "%.1f", nu[i]);
      if (buf[strlen(buf)-1] == '0') {
        len += strlen(buf)-2;
        buf[strlen(buf)-2] = '\0';
      }
      else
        len += strlen(buf);
      fprintf(pF, "%s", buf);
    }
    fprintf(pF, "%s", sp[i]);
    len += strlen(sp[i]);
    if (!strcmp(sp[i], "M") && fo > 1) fprintf(pF, ")");
    if (i<nSp-1) {
      if (!strcmp(sp[i+1], "M") && fo > 1) fprintf(pF, "(");
      fprintf(pF, "+");
    }
  }
  return len;
}
/* -------------------------------------------------------------------------- */

/* take the reaction starting at the first species and stop at delim */
int parseReactionSpecies(char *rxn, char *delim, double **nu, char ***sp)
{
  int i, j, iSp, nSp = 0;
  double lastNu = -1.0;
  char tag[64], nuC[16], buffer[256], *pS, splits[16];
  strcpy(buffer, rxn);

  /* allocate the memory to return species and coefficients */
  sp[0] = (char**)  calloc(1, sizeof(char*));
  nu[0] = (double*) calloc(1, sizeof(double));

  strcpy(splits, " +");
  pS = strtok(buffer, splits);
  while (pS != NULL && strcmp(pS, delim)) {
    if (strcmp(pS, "+")) {
      strcpy(tag, pS);

      /* get the coefficient for the species */
      i = 0;
      while (isdigit(tag[i]) && i < strlen(tag)) {
        i++;
      };

      if (i == strlen(tag)) { /* this was only a digit */
        lastNu = atof(tag);
        pS = strtok(NULL, splits);
        continue;
      }

      /* count the species and copy necessary values */
      nSp++;
      iSp = nSp - 1;
      sp[0]      = (char**)  realloc(sp[0], nSp*sizeof(char*));
      sp[0][iSp] = (char*)    calloc(32, sizeof(char));
      nu[0]      = (double*) realloc(nu[0], nSp*sizeof(double));

      if (lastNu > 0.0) {
        strcpy(sp[0][iSp], tag);
        nu[0][iSp] = lastNu;
        lastNu = -1.0;
      }
      else {
        if (i == 0 || (tag[i] == '-')) { /* no leading digit found */
          if      (tag[0] == '1') tag[0] = 'P';
          else if (tag[0] == '2') tag[0] = 'S';
          else if (tag[0] == '3') tag[0] = 'T';
          for (j=0; j<strlen(tag); j++) tag[j] = toupper(tag[j]);
          strcpy(sp[0][iSp], tag);
          nu[0][iSp] = 1;
        }
        else { /* leading digit found */
          strcpy(sp[0][iSp], &tag[i]);
          strncpy(nuC, tag, i);
          nuC[i] = '\0';
          nu[0][iSp] = atof(nuC);
        }
      }
    }
    pS = strtok(NULL, splits);
  }

  return nSp;
}
/* -------------------------------------------------------------------------- */

void parseMolecule(char *molID, molStruct *mol, int nElem, char **elems,
                   chemStruct *chem)
{
  /* this function assumes that preceding number coefficients have been
     removed */
  int i,j;

  strcpy(mol->name, molID);
  for (i=0; i<strlen(molID); i++) {
    /* check whether this is a two letter element */
  }
}
/* -------------------------------------------------------------------------- */

int getElementIndex(char *element, chemStruct *chem)
{
  int i;

  /* convert any all uppercase elements are */
  char atom[3];
  strcpy(atom, element);
  stripString(atom);
  for (i=1; i<3; i++) atom[i] = tolower(atom[i]);
  
  for (i=0; i<chem->nPT; i++) {
    if (!strcmp(atom, chem->pt[i].symb))
      return i;
  }
  return -1;
}
/* -------------------------------------------------------------------------- */

atomStruct * getElement(char *element, chemStruct *chem)
{
  int indx;
  indx = getElementIndex(element, chem);
  if (indx >= 0)
    return &chem->pt[indx];

  fprintf(stderr, "Error: element %s not found in periodic table\n", element);
  return NULL;
}
/* -------------------------------------------------------------------------- */

int readChemKinThermo(char *fileName, double *Trange, char ***species,
                      char ***info, char **phase, double ***spT, int **nAtoms,
                      int ***nAt, char ****atoms, double ***coeff)
{
  int i, j, iSp, iL, nSpecies;
  int maxAtom = 4;
  double Tmin, Tmid, Tmax, dum;
  FILE *pF = NULL;
  fpos_t pos;
  char tag[32], spLn[4][128], buffer[128];

  int  *nA   = (int*)   calloc(10, sizeof(int));
  char **aID = (char**) calloc(10, sizeof(char*));
  char *nTmp = (char*)  calloc(3, sizeof(char));
  for (i=0; i<10; i++)
    aID[i] = (char*) calloc(2, sizeof(char));

  pF = fopen(fileName, "r");
  if (pF == NULL) {
    fprintf(stderr, "Error: cannot read ChemKin file %s\n", fileName);
    exit(2);
  }

  /* pre-allocate memory locations */
  species[0] = (char**)   calloc(1, sizeof(char*));
  info[0]    = (char**)   calloc(1, sizeof(char*));
  phase[0]   = (char*)    calloc(1, sizeof(char));
  spT[0]     = (double**) calloc(1, sizeof(double*));
  nAtoms[0]  = (int*)     calloc(1, sizeof(int));
  nAt[0]     = (int**)    calloc(1, sizeof(int*));
  atoms[0]   = (char***)  calloc(1, sizeof(char**));
  coeff[0]   = (double**) calloc(1, sizeof(double*));

  nSpecies = 0;
  while (!feof(pF)) {
    fscanf(pF, "%s\n", tag);

    if (!strcmp(tag, "THERMO")) {
      fgetpos(pF, &pos);
      fscanf(pF, "%s", tag);
      if (!strcmp(tag, "END")) {
        fprintf(stderr, "Error: no species in ChemKin file %s\n", fileName);
        exit(2);
      }
      fsetpos(pF, &pos);

      /* otherwise, read in the temperature limits */
      fgets(buffer, 128, pF);
      sscanf(buffer, "%10lf%10lf%10lf",&Trange[0],&Trange[1],&Trange[2]);

      fgetpos(pF, &pos);
      do {
        fsetpos(pF, &pos);
        nSpecies++;
        iSp = nSpecies - 1;

        species[0] = (char**)   realloc(species[0], nSpecies*sizeof(char*));
        info[0]    = (char**)   realloc(info[0],    nSpecies*sizeof(char*));
        phase[0]   = (char*)    realloc(phase[0],   nSpecies*sizeof(char));
        spT[0]     = (double**) realloc(spT[0],     nSpecies*sizeof(double*));
        nAtoms[0]  = (int*)     realloc(nAtoms[0],  nSpecies*sizeof(int));
        nAt[0]     = (int**)    realloc(nAt[0],     nSpecies*sizeof(int*));
        atoms[0]   = (char***)  realloc(atoms[0],   nSpecies*sizeof(char**));
        coeff[0]   = (double**) realloc(coeff[0],   nSpecies*sizeof(double*));
        
        /* allocate other arrays if necessary */
        species[0][iSp] = (char*)   calloc(18, sizeof(char));
        info[0][iSp]    = (char*)   calloc(8,  sizeof(char));
        spT[0][iSp]     = (double*) calloc(3,  sizeof(double));
        coeff[0][iSp]   = (double*) calloc(14, sizeof(double));

        /* ----- read in 4 lines for each species ----- */
        for (i=0; i<4; i++) {
          /* read the line */
          fgets(spLn[i], 128, pF);
        } /* end species loop */

        /* read in the species and date */
        strncpy(species[0][iSp], spLn[0], 18);
        strncpy(info[0][iSp], &spLn[0][18], 6);
        info[0][iSp][6] = '\0';

        /* check whether there might be additional atoms in the first line */
        maxAtom = 4;
        if (strlen(spLn[0]) > 81) {
          maxAtom = (strlen(spLn[0]) - 24 - 1 - 30 - 5 - 1)/5;
          fprintf(stderr, "Warning: non-standard ChemKin Thermo format for ");
          fprintf(stderr, "species %s: checking for %i atoms\n", 
                  species[0][iSp], maxAtom);
        }

        /* read in the atomic elements */
        j = 24;
        nAtoms[0][iSp] = 0;
        for (i=0; i<maxAtom; i++) {
          strncpy(aID[i], &spLn[0][j], 2); j += 2;
          strncpy(nTmp,   &spLn[0][j], 3); j += 3;
          nA[i] = atoi(nTmp);
          if (nA[i] > 0)
            nAtoms[0][iSp]++;
        }

        /* read in the phase and the temperature range */
        sscanf(&spLn[0][24+maxAtom*5],"%c%10lf%10lf%10lf%5i",&phase[0][iSp],
               &spT[0][iSp][0],&spT[0][iSp][1],&spT[0][iSp][2],&iL);

        /* get the first 5 coefficients */
        sscanf(spLn[1],"%15lf%15lf%15lf%15lf%15lf%5i",&coeff[0][iSp][0],
               &coeff[0][iSp][1],&coeff[0][iSp][2],&coeff[0][iSp][3],
               &coeff[0][iSp][4],&iL);
        /* get the coefficients 6-10 */
        sscanf(spLn[2],"%15lf%15lf%15lf%15lf%15lf%5i",&coeff[0][iSp][5],
               &coeff[0][iSp][6],&coeff[0][iSp][7],&coeff[0][iSp][8],
               &coeff[0][iSp][9],&iL);
        /* get the coefficients 11-14 */
        sscanf(spLn[3],"%15lf%15lf%15lf%15lf%15lf%5i",&coeff[0][iSp][10],
               &coeff[0][iSp][11],&coeff[0][iSp][12],&coeff[0][iSp][13],
               &dum,&iL);

        /* assign the species appropriately */
        nAt[0][iSp]   = (int*)   calloc(nAtoms[0][iSp], sizeof(int));
        atoms[0][iSp] = (char**) calloc(nAtoms[0][iSp], sizeof(char*));
        for (i=0; i<nAtoms[0][iSp]; i++) {
          atoms[0][iSp][i] = (char*) calloc(2, sizeof(char));
        }
        for (i=0, j=0; i<maxAtom; i++) {
          if (nA[i] > 0) {
            nAt[0][iSp][j] = nA[i];
            strcpy(atoms[0][iSp][j], aID[i]);
            j++;
          }
        }

        fgetpos(pF, &pos);
        fscanf(pF, "%s", tag);
      } while (strcmp(tag, "END") && !feof(pF));
    } /* end thermo read */
  } /* end file read */

  fclose(pF);
  free(nTmp);
  for (i=0; i<10; i++) free(aID[i]);
  free(aID);
  free(nA);
  return nSpecies;
}
/* -------------------------------------------------------------------------- */

int writeChemKinThermo(char *fileName, chemStruct *chem)
{
  int   i,j,k,n;
  char  buffer[32], buff[8], *isomer, species[32];
  double coeff[14];
  FILE *pF = NULL;

  pF = fopen(fileName, "w");
  if (pF == NULL) {
    fprintf(stderr, "Error: cannot open file for writing: %s\n", fileName);
    exit(2);
  }

  fprintf(pF, "THERMO\n");
  for (i=0; i<3; i++)
    fprintf(pF, "%10.0f", chem->Trange[i]);
  fprintf(pF, "\n");

  for (i=0; i<chem->nSpecies; i++) {
    /* if the species name is longer than 18 characters, check whether it
       is an isomer and write the longest section possible */
    strcpy(buffer, chem->specs[i].name);
    strcpy(species, chem->specs[i].name);
    if (strlen(buffer) > 18) {
      isomer = strrchr(buffer, '-');
      if (!isomer) {
        fprintf(stderr, "Warning: ChemKin species longer than 18 chars (%s)\n",
                buffer);
      }
      else {
        do {
          strncpy(species, buffer, strlen(buffer)-strlen(isomer));
          species[strlen(buffer)-strlen(isomer)] = '\0';
          isomer=strrchr(isomer-1, '-');
        } while (isomer != NULL && strlen(buffer)-strlen(isomer) > 18);
        fprintf(stderr, "Warning: using isomer %s in place of %s\n", species,
                buffer);
      }
    }

    /* line 1 */
    fprintf(pF, "%-18s%-6s", species, chem->specs[i].info);
    if (chem->specs[i].nAtoms < 4)
      n = chem->specs[i].nAtoms;
    else
      n = 4;
    for (j=0; j<n; j++) {
      strcpy(buff, chem->specs[i].atom[j]->symb);
      for (k=0; k<strlen(buff); k++) buff[k] = toupper(buff[k]);
      fprintf(pF, "%-2s%3i", buff, chem->specs[i].nA[j]);
    }
    fprintf(pF, "%*s", (4-n)*5, "");
    fprintf(pF, "%c", chem->specs[i].phase);
    fprintf(pF, "%10.0f", chem->specs[i].Trange[0]);
    fprintf(pF, "%10.0f", chem->specs[i].Trange[2]);
    fprintf(pF, "%8.0f",  chem->specs[i].Trange[1]);
    if (chem->specs[j].nAtoms > 4) {
      strcpy(buff, chem->specs[i].atom[4]->symb);
      for (k=0; k<strlen(buff); k++) buff[k] = toupper(buff[k]);
      fprintf(pF, "%2s%3i 1", buff, chem->specs[i].nA[4]);
    }
    else
      fprintf(pF, "%7i\n", 1);

    /* line 2 */
    for (j=0; j<5; j++)
      fprintf(pF, "%15.8E", chem->specs[i].therm[1][j]);
    fprintf(pF, "%5i\n", 2);

    /* line 3 */
    for (j=5; j<7; j++)
      fprintf(pF, "%15.8E", chem->specs[i].therm[1][j]);
    for (j=0; j<3; j++)
      fprintf(pF, "%15.8E", chem->specs[i].therm[0][j]);
    fprintf(pF, "%5i\n", 3);

    /* line 4 */
    for (j=3; j<7; j++)
      fprintf(pF, "%15.8E", chem->specs[i].therm[0][j]);
    fprintf(pF, "%15s%5i\n", "", 4);
  }

  fclose(pF);
  return 0;
}
/* -------------------------------------------------------------------------- */

int readChemKinTrans(char *fileName, char ***species, double **eps_k,
                     double **sigma)
{
  int  i,nSpecies, nHead, empty;
  double dum;

  FILE *pF = NULL;
  /* strings to store the tag for each line and comments (if any) */
  char tag[32], comment[256], header[256];
 
  pF = fopen(fileName, "r");
  if (pF == NULL) {
    fprintf(stderr, "Error: cannot read ChemKin transport file %s\n", fileName);
    exit(2);
  }

  /* ----- read in the file ----- */
  /* get the header line and get the number of species*/
  nSpecies = nHead = 0;
  fgets(comment, 256, pF);
  while (!feof(pF)) {
    if (comment[0] == '#' || comment[0] == '!') {
      nHead++;
      strcpy (header, comment);
    }
    else {
      empty = 1;
      if (comment[0] != '\n') {
        for (i=0; i<strlen(comment); i++) {
          if (!isblank(comment[i]))
            empty = 0;
        }
      }

      if (!empty)
        nSpecies++;
    }
    fgets(comment, 256, pF);
  }
  rewind(pF);

  /* TODO: check the header */
  if (nHead > 0) {
    if (nHead > 1) {
      fprintf(stderr, "Warning: transport file has multiple header lines\n");
    }
  }

  if (MDUC_PRINT_PROC) fprintf(stderr, " Reading %i species from transport file\n", nSpecies);

  species[0] = (char**)  calloc(nSpecies, sizeof(char*));
  eps_k[0]   = (double*) calloc(nSpecies, sizeof(double));
  sigma[0]   = (double*) calloc(nSpecies, sizeof(double));
  for (i=0; i<nSpecies; i++) {
    species[0][i] = (char*) calloc(15, sizeof(char));
  }

  /* re-read the header */
  for (i=0; i<nHead; i++) {
    fgets(comment, 256, pF);
  }

  /* read the data */
  /* TODO: generalize this to take a column number, without assuming order */
  for (i=0; i<nSpecies; i++) {
    readString(pF, species[0][i]);
    dum         = readReal(pF);
    eps_k[0][i] = readReal(pF);
    sigma[0][i] = readReal(pF);
    dum         = readReal(pF);
    dum         = readReal(pF);
    dum         = readReal(pF);
  }

  fclose(pF);
  return nSpecies;
}
/* -------------------------------------------------------------------------- */

int writeChemKinTrans(char *fileName, chemStruct *chem)
{
  int    i;
  int    *iTmp = (int*)    calloc(chem->nSpecies, sizeof(int));
  double *fTmp = (double*) calloc(chem->nSpecies, sizeof(double));

  int    *i1;
  double *f1,*f2,*f3,*f4,*f5;

  /* set default to zero arrays */
  i1 = iTmp;
  f1 = f2 = f3 = f4 = f5 = fTmp;

  FILE *pF = NULL;
 
  pF = fopen(fileName, "w");
  if (pF == NULL) {
    fprintf(stderr, "Error: cannot open file for writing: %s\n", fileName);
    exit(2);
  }

  /* ----- write the header ----- */
  fprintf(pF, "%-18s", "#SpecName");
  if (chem->atomGeom == NULL) {
    fprintf(pF, "%10s", "Dummy");
  }
  else {
    fprintf(pF, "%10s", "Geometry");
    i1 = chem->atomGeom;
  }
  if (chem->eps_k == NULL) {
    fprintf(pF, "%10s", "Dummy");
  }
  else {
    fprintf(pF, "%12s", "eps/k");
    f1 = chem->eps_k;
  }
  if (chem->sigma == NULL) {
    fprintf(pF, "%12s", "Dummy");
  }
  else {
    fprintf(pF, "%12s", "sigma");
    f2 = chem->sigma;
  }
  if (chem->muDipole == NULL) {
    fprintf(pF, "%9s", "Dummy");
  }
  else {
    fprintf(pF, "%9s", "mu");
    f3 = chem->muDipole;
  }
  if (chem->alphaPolar == NULL) {
    fprintf(pF, "%9s", "Dummy");
  }
  else {
    fprintf(pF, "%9s", "alpha");
    f4 = chem->alphaPolar;
  }
  if (chem->Zrot == NULL) {
    fprintf(pF, "%9s", "Dummy");
  }
  else {
    fprintf(pF, "%9s", "Zrot");
    f5 = chem->Zrot;
  }
  fprintf(pF, "\n");
  /* ----- end header ----- */

  char buffer[32], species[32], *isomer;
  for (i=0; i<chem->nSpecies; i++) {
    if (!chem->transport[i])
      continue;

    /* if the species name is longer than 18 characters, check whether it
       is an isomer and write the longest section possible */
    strcpy(buffer, chem->specs[i].name);
    strcpy(species, chem->specs[i].name);
    if (strlen(buffer) > 18) {
      isomer = strrchr(buffer, '-');
      if (!isomer) {
        fprintf(stderr, "Warning: ChemKin species longer than 18 chars (%s)\n",
                buffer);
      }
      else {
        do {
          strncpy(species, buffer, strlen(buffer)-strlen(isomer));
          species[strlen(buffer)-strlen(isomer)] = '\0';
          isomer=strrchr(isomer-1, '-');
        } while (isomer != NULL && strlen(buffer)-strlen(isomer) > 18);
        fprintf(stderr, "Warning: using isomer %s in place of %s\n", species,
                buffer);
      }
    }

    fprintf(pF, "%-18s%10i%12.3f%12.3f%9.3f%9.3f%9.3f\n",species,
            i1[i],f1[i],f2[i],f3[i],f4[i],f5[i]);
  }

  fclose(pF);
  free(iTmp); free(fTmp);
  return 0;
}
/* -------------------------------------------------------------------------- */

void writeMechCompiledC(chemStruct *chem)
{
  int i,j,k,kk;
  char fileName[128], **spTag;
  char mechUp[128];
  int  mechFallOff = 0;
  FILE *pF;
  rxnStruct *r, *r2;

  /* construct the species tags */
  spTag = (char**) calloc(chem->nSpecies, sizeof(char*));
  for (i=0; i<chem->nSpecies; i++) {
    spTag[i] = (char*) calloc(32, sizeof(char));
    strcpy(spTag[i], chem->species[i]);
    for (j=0; j<strlen(spTag[i]); j++) {
      if (spTag[i][j] == '-')
        spTag[i][j] = 'X';
    }
  }

  strcpy(mechUp, chem->mechName);
  for (j=0; j<strlen(mechUp); j++) mechUp[j] = toupper(mechUp[j]);

  /* ----- write the header file ----- */
  sprintf(fileName, "comp_%s.h", chem->mechName);
  pF = fopen(fileName, "w");
  fprintf(pF, "#ifndef MDUC_MECH_%s_H\n", mechUp);
  fprintf(pF, "#define MDUC_MECH_%s_H\n\n", mechUp);

  /* define the species numbering */
  fprintf(pF, "typedef enum speciesLabel {\n");
  for (i=0; i<chem->nSpecies; i++)
    fprintf(pF, "  s%s = %i,\n", spTag[i], i);
  fprintf(pF, "  sEnd\n");
  fprintf(pF, "} speciesLabel;\n\n");

  /* reaction numbering */
  fprintf(pF, "typedef enum reactionLabel {\n");
  j = 0;
  for (i=0; i<chem->nRxn; i++) {
    r = &chem->rxns[i];
    if (r->reverse == 0) {
      fprintf(pF, "  r%i = %i,\n", r->no, j);
      j++;
    }
    else {
      fprintf(pF, "  r%if = %i,\n", r->no, j); j++;
      fprintf(pF, "  r%ib = %i,\n", r->no, j); j++;
    }
  }
  fprintf(pF, "  rEnd\n");
  fprintf(pF, "} reactionLabel;\n\n");

  /* determine number of third body enhancements */
  int nM  = 0;
  int *M  = (int*) calloc(1, sizeof(int));
  int *iM = (int*) calloc(chem->nRxn, sizeof(int));
  int enhMatch, noEnh = 0;

  /* check whether there are any reactions with M but no enhancement */
  /* TODO: this whole checking M stuff should be cleaned up */
  for (i=0; i<chem->nRxn; i++) {
    r = &chem->rxns[i];
    if (r->pressure > 0 && r->nEnhance == 0) {
      nM++;
      M = (int*) realloc(M, nM*sizeof(int));
      M[nM-1] = i;
      iM[i] = 0;
      break;
    }
  } /* end loop over reactions */

  int check = 0;
  for (i=0; i<chem->nRxn; i++) {
    r = &chem->rxns[i];
    /* figure out how many different 'M's need to be written */
    if (r->pressure > 0) {
      if (r->nEnhance == 0) { /* already written */
        iM[i] = 0;
      }
      else { /* check for match */
        check = 0;
        for (j=0; j<i; j++) {
          r2 = &chem->rxns[j];
          enhMatch = 0;
          if (r2->pressure > 0) {
            check++;
            if (r2->nEnhance == r->nEnhance) {
              for (k=0; k<r->nEnhance; k++) {
                for (kk=0; kk<r2->nEnhance; kk++) {
                  if (!strcmp(r->spEnh[k], r2->spEnh[kk]) &&
                      r->factEnh[k] == r2->factEnh[kk] )
                    enhMatch++;
                }
              }
            }
            if (enhMatch != r->nEnhance) { /* new M */
              nM++;
              M = (int*) realloc(M, nM*sizeof(int));
              M[nM-1] = i;
              iM[i] = nM - 1;
              break;
            }
            else {
              for (k=0; k<nM; k++)
                if (j == M[k])
                  iM[i] = k;
              break;
            }
          }
        } /* end loop over previous reactions */
        if (check == 0) { /* new M */
          nM++;
          M = (int*) realloc(M, nM*sizeof(int));
          M[nM-1] = i;
          iM[i] = nM - 1;
        }
      }
    } /* end if enhance */
  } /* end loop over reactions */

  /* third body enhancement factors */
  fprintf(pF, "typedef enum thirdBodyLabel {\n");
  for (i=0; i<nM; i++)
    fprintf(pF, "  m%i = %i,\n", i, i);
  fprintf(pF, "  mEnd\n");
  fprintf(pF, "} thirdBodyLabel;\n\n");

  /* define internal functions */
  fprintf(pF, "double getFallOffRateCoeff(double T, double logT, double RT, "
          "double P, \n%27sdouble k0, double kInf, int troe, double Fc,\n"
          "%27sdouble conc);\n", "", "");
  
  fprintf(pF, "\n#endif /* MDUC_MECH_%s_H */\n", mechUp);
  fclose(pF);

  /* ----- write the C functions ----- */
  sprintf(fileName, "comp_%s.c", chem->mechName);
  pF = fopen(fileName, "w");
  fprintf(pF, "#include <stdlib.h>\n");
  fprintf(pF, "#include <string.h>\n");
  fprintf(pF, "#include <math.h>\n");
  fprintf(pF, "#include \"comp_%s.h\"\n\n", chem->mechName);

  /* write the mechanism info */
  fprintf(pF, "int  GetNSpecies()\n");
  fprintf(pF, "{\n");
  fprintf(pF, "  return sEnd;\n");
  fprintf(pF, "}\n\n");

  fprintf(pF, "int  GetNSpecs()\n");
  fprintf(pF, "{\n");
  fprintf(pF, "  return %i;\n", chem->nSpecies);
  fprintf(pF, "}\n\n");

  fprintf(pF, "int  GetNReactions()\n");
  fprintf(pF, "{\n");
  fprintf(pF, "  return rEnd;\n");
  fprintf(pF, "}\n\n");

  fprintf(pF, "void ComputeProductionRates(double *cdot, double *w, double *k, "
          "double *c,\n                            double *M, double T, "
          "double P)\n");
  fprintf(pF, "{\n");
  
  /* write the header info */
  fprintf(pF, "  int nSpec   = %i;\n", chem->nSpecSolve);
  fprintf(pF, "  int nSpecIn = %i;\n", chem->nSpecies);
  fprintf(pF, "  double k0, kInf, fcTroe;\n");
  fprintf(pF, "  double RGAS = 8314.34;\n");
  fprintf(pF, "  double logT = log(T);\n");
  fprintf(pF, "  double RT   = RGAS*T;\n");
  
  /* write the rate coefficients */
  fprintf(pF, "\n  /* ----- rate coefficients ----- */\n");
  for (i=0; i<chem->nRxn; i++) {
    r = &chem->rxns[i];
    if (r->pressure > 1) { /* fall off reaction */
      mechFallOff = 1;
      /* write k0 and kInf */
      writeMechRateCoeff_C(pF, r->no, -2, r->A_low, r->b_low, r->E_low);
      writeMechRateCoeff_C(pF, r->no, -3, r->A, r->b, r->E);
      /* write the Troe function if necessary */
      if (r->pressure == 2) { /* fall-off reaction */
        fprintf(pF, "  fcTroe = ");
        if (r->troeA[0] > 0.0 && r->troeT[0] > 0.0)
          fprintf(pF, "%g*exp(-T/%g)", r->troeA[0], r->troeT[0]);
        if (r->troeA[1] > 0.0 && r->troeT[1] > 0.0)
          fprintf(pF, " + %g*exp(-T/%g)", r->troeA[1], r->troeT[1]);
        if (r->troeA[2] > 0.0 && r->troeT[2] >= 0.0)
          fprintf(pF, " + %g*exp(-%g/T)", r->troeA[2], r->troeT[2]);
        fprintf(pF, ";\n");
      }
      else
        fprintf(pF, "  fcTroe = -1.0\n");
      writeMechRateCoeff_C(pF, r->no, 4, 0.0, 0.0, 0.0);
    }
    else
      writeMechRateCoeff_C(pF, r->no, r->reverse, r->A, r->b, r->E);

    /* write the reverse rate */
    if (r->reverse > 0) {
      if (r->pressure > 1) { /* fall-off reaction */
        writeMechRateCoeff_C(pF, r->no, -2, r->Ar_low, r->br_low, r->Er_low);
        writeMechRateCoeff_C(pF, r->no, -3, r->Ar, r->br, r->Er);
        writeMechRateCoeff_C(pF, r->no, -4, 0.0, 0.0, 0.0);
      }
      else
        writeMechRateCoeff_C(pF, r->no, -1, r->Ar, r->br, r->Er);
    }
  } /* end loop over reactions */

  /* ----- write the enhancement factors ----- */
  int line = 0, perLine = 4;
  fprintf(pF, "\n  /* ----- enhancment factors ----- */\n");
  for (i=0; i<nM; i++) {
    fprintf(pF, "  M[m%i] = ", i);
    line = 0;
    r = &chem->rxns[M[i]];
    for (j=0; j<chem->nSpecies; j++) {
      if (r->concF[j] > 0.0) {
        fprintf(pF, " %+g*c[s%s]", r->concF[j], spTag[j]);
        line++;
        if (line != 0 && line%perLine == 0)
          fprintf(pF, "\n    ");
      }
    }
    fprintf(pF, ";\n\n");
  }

  /* ----- write the reaction rates ----- */
  fprintf(pF, "\n  /* ----- reaction rates ----- */\n");
  for (i=0; i<chem->nRxn; i++) {
    r = &chem->rxns[i];
    if (r->reverse > 0) {
      fprintf(pF, "  w[r%if] = k[r%if]", r->no, r->no);
      for (j=0; j<r->nR; j++) {
        for (k=0; k<r->nuR[j]; k++) {
          if (r->iR[j] >= 0)
            fprintf(pF, "*c[s%s]", spTag[r->iR[j]]);
        }
      }
      if (r->pressure == 1)
        fprintf(pF, "*M[m%i]", iM[i]);
      fprintf(pF, ";\n");
      fprintf(pF, "  w[r%ib] = k[r%ib]", r->no, r->no);
      for (j=0; j<r->nP; j++) {
        for (k=0; k<r->nuP[j]; k++) {
          if (r->iP[j] >= 0)
            fprintf(pF, "*c[s%s]", spTag[r->iP[j]]);
        }
      }
      if (r->pressure == 1)
        fprintf(pF, "*M[m%i]", iM[i]);
      fprintf(pF, ";\n");
    }
    else {
      fprintf(pF, "  w[r%i] = k[r%i]", r->no, r->no);
      for (j=0; j<r->nR; j++) {
        for (k=0; k<r->nuR[j]; k++) {
          if (r->iR[j] >= 0)
            fprintf(pF, "*c[s%s]", spTag[r->iR[j]]);
        }
      }
      if (r->pressure == 1)
        fprintf(pF, "*M[m%i]", iM[i]);
      fprintf(pF, ";\n");
    }
  }

  /* ----- write the concentration rates ----- */
  int dir;
  char co[8], dir0[8], dir1[8];
  fprintf(pF, "\n  /* ----- concentration rates ----- */\n");
  for (i=0; i<chem->nSpecies; i++) {
    line = 0;
    fprintf(pF, "  cdot[s%s] =", spTag[i]);

    if (chem->specs[i].nRxn == 0) {
      fprintf(pF, " 0.0");
      line = 1;
    }

    /* reactions with species present */
    for (j=0; j<chem->specs[i].nRxn; j++) {
      r = &chem->rxns[chem->specs[i].iRxn[j]];
      strcpy(co, "");
      if (fabs(chem->specs[i].nuRxn[j]) > 1.0)
        sprintf(co, "%.1f*", fabs(chem->specs[i].nuRxn[j]));
      strcpy(dir0, " - ");
      strcpy(dir1, " + ");
      if (chem->specs[i].nuRxn[j] > 0.0) {
        strcpy(dir0, " + ");
        strcpy(dir1, " - ");
      }
      fprintf(pF, "%s%sw[r%i", dir0, co, r->no);
      if (r->reverse > 0) {
        fprintf(pF, "f]");
        line++;
        if (line != 0 && line%perLine == 0 && j != chem->specs[i].nRxn-1)
          fprintf(pF, "\n    ");
        fprintf(pF, "%s%sw[r%ib]", dir1, co, r->no);
      }
      else
        fprintf(pF, "]");

      line++;
      if (line != 0 && line%perLine == 0 && j != chem->specs[i].nRxn-1)
        fprintf(pF, "\n    ");
    }

    fprintf(pF, ";\n\n");
  }
  fprintf(pF, "}\n");
  /* end compute production rate function*/

  if (mechFallOff) {
    fprintf(pF, "\n");
    writeFallOffRateFunc(pF);
    fprintf(pF, "\n");
  }

  /* ----- write the species names ----- */
  fprintf(pF, "void GetSpeciesNames(char **names)\n");
  fprintf(pF, "{\n");
  for (i=0; i<chem->nSpecies; i++) {
    fprintf(pF,"  strcpy(names[s%s],\"%-25s\");\n",spTag[i],chem->species[i]);
  }
  fprintf(pF, "}\n\n");

  /* ----- write the molecular weights ----- */
  fprintf(pF, "void GetMolarMass(double *W)\n");
  fprintf(pF, "{\n");
  for (i=0; i<chem->nSpecies; i++) {
    fprintf(pF,"  W[s%s] = %.8e;\n",spTag[i],chem->W[i]);
  }
  fprintf(pF, "}\n\n");

  /* ----- write the thermodynamic data ----- */
  double *a, *b, R_W;
  fprintf(pF, "void ComputeThermoData(double *h, double *cp, double T)\n");
  fprintf(pF, "{\n");
  fprintf(pF, "  int i;\n");
  /* write the high temperature coefficients */
  fprintf(pF, "  if (T > %.1f) {\n", chem->Trange[1]);
  for (i=0; i<chem->nSpecies; i++) {
    a = chem->specs[i].aH[1];
    b = chem->specs[i].aCp[1];
    fprintf(pF, "    h[s%s]  = %.8e*(%.8e \n",spTag[i],8314.34/chem->W[i],a[0]);
    fprintf(pF, "      + T*(%.8e + T*(%.8e \n", a[1], a[2]);
    fprintf(pF, "      + T*(%.8e + T*(%.8e + \n", a[3], a[4]);
    fprintf(pF, "      + T*(%.8e))))));\n", a[5]);
    fprintf(pF, "    cp[s%s] = %.8e*(%.8e \n",spTag[i],8314.34/chem->W[i],b[0]);
    fprintf(pF, "      + T*(%.8e + T*(%.8e \n", b[1], b[2]);
    fprintf(pF, "      + T*(%.8e + T*(%.8e)))));\n", b[3], b[4]);
  }
  fprintf(pF, "  }\n");  
  /* write the low temperature coefficients */
  fprintf(pF, "  else if (T >= %.1f) {\n", chem->Trange[0]);
  for (i=0; i<chem->nSpecies; i++) {
    a = chem->specs[i].aH[0];
    b = chem->specs[i].aCp[0];
    fprintf(pF, "    h[s%s]  = %.8e*(%.8e \n",spTag[i],8314.34/chem->W[i],a[0]);
    fprintf(pF, "      + T*(%.8e + T*(%.8e \n", a[1], a[2]);
    fprintf(pF, "      + T*(%.8e + T*(%.8e + \n", a[3], a[4]);
    fprintf(pF, "      + T*(%.8e))))));\n", a[5]);
    fprintf(pF, "    cp[s%s] = %.8e*(%.8e \n",spTag[i],8314.34/chem->W[i],b[0]);
    fprintf(pF, "      + T*(%.8e + T*(%.8e \n", b[1], b[2]);
    fprintf(pF, "      + T*(%.8e + T*(%.8e)))));\n", b[3], b[4]);
  }
  fprintf(pF, "  }\n");
  fprintf(pF, "  else {\n");
  fprintf(pF, "    ComputeThermoData(h, cp, 300.0);\n");
  fprintf(pF, "    for (i=0; i<sEnd; i++) {\n");
  fprintf(pF, "      h[i] = (T - 300.0)*cp[i] + h[i];\n");
  fprintf(pF, "    }\n");
  fprintf(pF, "  }\n");
  fprintf(pF, "}\n\n");

  fclose(pF);
}
/* -------------------------------------------------------------------------- */

void writeMechRateCoeff_C(FILE *pF, int i, int f, double A, double b, double E)
{
  if (f == 0)
    fprintf(pF, "  k[r%i] = ", i);
  else if (f > 0)
    fprintf(pF, "  k[r%if] = ", i);
  else if (f == -2)
    fprintf(pF, "  k0   = ");
  else if (f == -3)
    fprintf(pF, "  kInf = ");
  else
    fprintf(pF, "  k[r%ib] = ", i);

  if (f == -4 || f == 4) {
    fprintf(pF, "getFallOffRateCoeff(T,logT,RT,P,k0,kInf,1,fcTroe,-1.0);\n");
    return;
  }
  
  fprintf(pF, "%.10e", A);
  if (fabs(b) > 0 || fabs(E) > 0) {
    fprintf(pF, "*exp(");
    if (fabs(b) > 0)
      fprintf(pF, "%f*logT", b);
    if (fabs(E) > 0)
      fprintf(pF, "%+.9e/RT", -1.0*E);
    fprintf(pF, ")");
  }
  fprintf(pF, ";\n");
}
/* -------------------------------------------------------------------------- */

void writeFallOffRateFunc(FILE *pF)
{
  fprintf(pF, "double getFallOffRateCoeff(double T, double logT, "
          "double RT, double P,\n");
  fprintf(pF, "                           double k0, double kInf, int troe,"
          "double Fc,\n");
  fprintf(pF, "                           double conc)\n");
  fprintf(pF, "{\n");
  fprintf(pF, "  double Pr,logPr;  /* reduced pressure*/\n");
  fprintf(pF, "  double F = 1.0;   /* Troe function (default unity value is "
          "for Lindemann */\n");
  fprintf(pF, "  double c,d,n;     /* Troe or SRI coefficients */\n");
  fprintf(pF, "  double logFc;\n");
  fprintf(pF, "  double tmp;\n\n");

  fprintf(pF, "  if (conc <= 0.0) {\n");
  fprintf(pF, "    conc = P/(RT);\n");
  fprintf(pF, "  }\n\n");

  fprintf(pF, "  /* compute the reduced pressure */\n");
  fprintf(pF, "  Pr = k0*conc/kInf;\n\n");

  fprintf(pF, "  if (troe == 1) { /* use Troe form */\n");
  fprintf(pF, "    logPr = log10(Pr);\n");
  fprintf(pF, "    logFc = log10(Fc);\n");
  fprintf(pF, "    c = - 0.4  - 0.67*logFc;\n");
  fprintf(pF, "    n =   0.75 - 1.27*logFc;\n");
  fprintf(pF, "    d =   0.14;\n");
  fprintf(pF, "    tmp = (logPr+c)/(n-d*(logPr+c));\n");
  fprintf(pF, "    F = pow(Fc, 1.0/(1.0+tmp*tmp));\n");
  fprintf(pF, "  }\n");

  fprintf(pF, "  return kInf*(Pr/(1.0+Pr))*F;\n");
  fprintf(pF, "}\n");
}
/* -------------------------------------------------------------------------- */
