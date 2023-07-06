#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "flameMasterIO.h"
#include "chemKinIO.h"
#include "chemUtils.h"
#include "init.h"
#include "utils.h"
#include "parser.h"

int readFlameMasterMech(char *fileName, int *nAtoms, char ***atoms, int *nSpecs,
                        char ***species, rxnStruct **rxns, chemStruct *chem)
{
  int i, j, k, iSp, iL, iR, n, n0, nElem, nSpecies, nRxn, nLine, done;
  double *cV, tmp;
  FILE *pF = NULL;
  fpos_t pos;
  char tag[64], tag2[32], tmpF[32], rxnLn[256], buffer[256], *pBE, *pE, *pS;
  char *pReact, *pProd, *pArr, buffSp[256];
  char *rM;

  char rxnOrd[8] = {"n"};
  char Texp[8]   = {"n_k"};
  int nAddSp;
  char **addSp;

  /* space for enhancement factors */
  int nEnh = 0;
  int   *nEnhSp = (int*)     calloc(1, sizeof(int));
  double **enhF = (double**) calloc(1, sizeof(double*));
  char  **enhM  = (char**)   calloc(1, sizeof(char*));
  char ***enhSp = (char***)  calloc(1, sizeof(char**));

  rxnStruct *r;

  pF = fopen(fileName, "r");
  if (pF == NULL) {
    fprintf(stderr, "Error: cannot read FlameMaster file %s\n", fileName);
    exit(2);
  }

  /* pre-allocate some memory */
  atoms[0]   = (char**)     calloc(1, sizeof(char*));
  species[0] = (char**)     calloc(1, sizeof(char*));
  rxns[0]    = (rxnStruct*) calloc(1, sizeof(rxnStruct));

  nElem = nSpecies = nRxn = 0;
  while (!feof(pF)) {
    /* get the line */
    fgets(rxnLn, 256, pF);
    
    if (rxnLn[0] == '#') {
      /* comment, do nothing */
    }
    else if (strstr(rxnLn,"/*")) {
      do {
        fscanf(pF, "%s", tag);
      } while(strstr(tag, "*/") == NULL);
    }
    else if (!strncmp(rxnLn, "Let", 3) || !strncmp(rxnLn, "let", 3)) {
      /* get the location of the "be" */
      pBE = strstr(rxnLn, "be");

      /* copy the identifier */
      if (pBE) {
        strncpy(buffer, &rxnLn[4], strlen(rxnLn)-strlen(pBE)-5);
        buffer[strlen(rxnLn)-strlen(pBE)-5] = '\0';
        if (!strcmp(buffer, "allowed atoms")) {
          strcpy(tag, &pBE[3]);
          nElem = 0;
          
          /* get the first element */
          pE = strtok(tag, ",.");
          while (pE != NULL) {
            while (pE[0] == ' ') pE++;
            if (pE[0] != '\n') {
              nElem++;
              iL = nElem - 1;
              atoms[0]     = (char**) realloc(atoms[0], nElem*sizeof(char*));
              atoms[0][iL] = (char*)   calloc(3, sizeof(char));
              
              strcpy(atoms[0][iL], pE);
              for (i=0; i<strlen(atoms[0][iL]); i++)
                atoms[0][iL][i] = toupper(atoms[0][iL][i]);
            }
            pE = strtok(NULL, ",.");
          }
          *nAtoms = nElem;
        }
        else if (!strcmp(buffer, "additional species")) {
          strcpy(buffSp, &pBE[3]);
          addSp = (char**) calloc(1, sizeof(char*));
          nAddSp = 0;
          pE = strtok(buffSp, ",.");
          while (pE != NULL) {
            while (pE[0] == ' ') pE++;
            if (pE[0] != '\n') {
              nAddSp++;
              addSp = (char**) realloc(addSp, nAddSp*sizeof(char*));
              addSp[nAddSp-1] = calloc(32, sizeof(char));
              strcpy(addSp[nAddSp-1], pE);
            }
            pE = strtok(NULL, ",.");
          }
        }
        else if (!strcmp(buffer, "temperature exponent")) {
          strcpy(Texp, &pBE[3]);
          Texp[strlen(Texp)-2] = '\0';
        }
        else if (!strcmp(buffer, "order of reaction")) {
          strcpy(rxnOrd, &pBE[3]);
          rxnOrd[strlen(rxnOrd)-2] = '\0';
        }
        else if (!strcmp(buffer, "units for A")) {
          strcpy(chem->unitsA, &pBE[3]);
          chem->unitsA[strlen(chem->unitsA)-2] = '\0'; /* remove the period */
        }
        else if (!strcmp(buffer, "units for E")) {
          strcpy(chem->unitsE, &pBE[3]);
          chem->unitsE[strlen(chem->unitsE)-2] = '\0'; /* remove the period */
        }
        else {
          fprintf(stderr, "Did not recognize \"Let\" option: %s\n", buffer);
        }
      }
      else { /* third body reaction */
        nEnh++;
        i = nEnh - 1;
        nEnhSp   = (int*)     realloc(nEnhSp, nEnh*sizeof(int));
        enhF     = (double**) realloc(enhF,   nEnh*sizeof(double*));
        enhF[i]  = (double*)   calloc(1, sizeof(double));
        enhM     = (char**)   realloc(enhM,   nEnh*sizeof(char*));
        enhSp    = (char***)  realloc(enhSp,  nEnh*sizeof(char**));
        enhSp[i] = (char**)   calloc(1, sizeof(char*));
        strcpy(buffer, &rxnLn[4]);
        pE = strtok(buffer, "=");
        enhM[i] = (char*) calloc(32, sizeof(char));
        if (pE[0] == '[') pE++;
        strcpy(enhM[i], pE);
        for (j=0; j<strlen(enhM[i]); j++)
          if (enhM[i][j] == ']')
            enhM[i][j] = '\0';
        stripString(enhM[i]);
        nEnhSp[i] = 0;
        pE = strtok(NULL, "]");
        while (pE != NULL) {
          if (strchr(pE, '[')) {
            nEnhSp[i]++;
            j = nEnhSp[i] - 1;
            enhF[i]     = (double*) realloc(enhF[i], nEnhSp[i]*sizeof(double));
            enhSp[i]    = (char**) realloc(enhSp[i], nEnhSp[i]*sizeof(char*));
            enhSp[i][j] = (char*) calloc(32, sizeof(char));
            while (pE[0] == ' ' || pE[0] == '+') pE++;
            sscanf(pE, "%lf%s", &enhF[i][j], buffSp);
            strcpy(enhSp[i][j], &buffSp[1]);
          }
          pE = strtok(NULL, "]");
        }
      }      
    } /* end "Let" */
    else if (strlen(rxnLn) > 0 && rxnLn[0] != '\n') { /* reaction */
      /* convert any tabs to spaces */
      for (i=0; i<strlen(rxnLn); i++) {
        if (rxnLn[i] == '\t')
          rxnLn[i] = ' ';
      }
      
      /* don't operate on original line to prevent corruption */
      strcpy(buffer, rxnLn);
      /* assume that the first line in any reaction has <no.><dir>: at start */
      rM = strchr(buffer, ':');
      if (rM) {
        nRxn++;
        rxns[0] = (rxnStruct*) realloc(rxns[0], nRxn*sizeof(rxnStruct));
        r = &rxns[0][nRxn-1]; /* tmp for readability */

        /* set some defaults */
        setRxnDefault(r);
        
        /* ----- copy the ident ----- */
        n = strcspn(rxnLn, ":");
        /* get rid of any leading letters */
        /* TODO: should add reaction grouping? */
        n0 = 0;
        while (!isdigit(rxnLn[n0])) n0++;
        strncpy(tmpF, &rxnLn[n0], n-n0);
        n -= n0;
        if (isdigit(tmpF[n-1])) 
          tmpF[n] = '\0';
        else {
          /* check whether reaction is reversible */
          if (tmpF[n-1] == 'f')
            r->reverse = 1;
          else if (tmpF[n-1] == 'b')
            r->reverse = 2;
          tmpF[n-1] = '\0';
        }
        r->no = atoi(tmpF);

        /* ----- parse the reactants ----- */
        pReact = &rM[1];
        while (pReact[0] == ' ') pReact++;
        /* find the end of the reactants */
        n = strcspn(pReact, ">");
        strncpy(buffSp, pReact, n-1);
        buffSp[n-1] = '\0';
        r->nR = parseReactionSpecies(buffSp,"->",&r->nuR,&r->spR);

        /* ----- parse the products ----- */
        strcpy(buffer, rxnLn);
        pProd = strchr(buffer, '>');
        pProd++;
        while (pProd[0] == ' ') pProd++;
        /* find the end of the products (assumes at least 1 space before { ) */
        n = strcspn(pProd, "{");
        strncpy(buffSp, pProd, n-1);
        buffSp[n-1] = '\0';
        r->nP = parseReactionSpecies(buffSp,"{",&r->nuP,&r->spP);

        /* ----- parse the Arrhenius reaction coefficients ----- */
        done = 0;
        strcpy(buffer, rxnLn);
        stripString(buffer);
        pArr = strchr(buffer, '{');
        if (pArr == NULL) {
          fprintf(stderr, "Error: Arrhenius coefficients not found for rxn\n");
          fprintf(stderr, " %s", rxnLn);
        }
        pArr++;
        while (pArr[0] == ' ') pArr++;
        while (!done) {
          /* TODO: for now assume all variables have a single space before = */
          cV = NULL;
          pS = strtok(pArr, " ");
          while (pS != NULL && strcmp(pS, "}") && pS[strlen(pS)-1] != '\n') {
            if (strcmp(pS, "=")) {
              if (!strcmp(pS, "a"))
                cV = &r->A;
              else if (!strcmp(pS, "n"))
                cV = &r->b;
              else if (!strcmp(pS, "E"))
                cV = &r->E;
              else if (!strcmp(pS, "ai")) {
                r->pressure  = 2;
                r->lindemann = 1;
                cV = &r->A_low;
              }
              else if (!strcmp(pS, "ni"))
                cV = &r->b_low;
              else if (!strcmp(pS, "Ei"))
                cV = &r->E_low;
              else if (!strcmp(pS, "fca")) {
                r->lindemann = 0;
                r->troe = 1;
                cV = &r->troeA[0];
              }
              else if (!strcmp(pS, "fcta")) {
                r->troe = 1;
                cV = &r->troeT[0];
              }
              else if (!strcmp(pS, "fcb")) {
                r->troe = 1;
                cV = &r->troeA[1];
              }
              else if (!strcmp(pS, "fctb")) {
                r->troe = 1;
                cV = &r->troeT[1];
              }
              else if (!strcmp(pS, "fcc")) {
                r->troe = 1;
                cV = &r->troeA[2];
              }
              else if (!strcmp(pS, "fctc")) {
                r->troe = 1;
                cV = &r->troeT[2];
              }
              else {
                strcpy(tmpF, pS);
                cV[0] = atof(tmpF);
              }
            }
            pS = strtok(NULL, " ");
          }

          /* check whether there is more to read in for this reaction */
          pArr = strchr(rxnLn, '}');
          if (pArr == NULL) {
            fgets(rxnLn, 256, pF);
            for (i=0; i<strlen(rxnLn); i++)
              if (rxnLn[i] == '\t')
                rxnLn[i] = ' ';
            strcpy(buffer, rxnLn);
            stripString(buffer);
            pArr = &buffer[0];
            while (pArr[0] == ' ') pArr++;
          }
          else
            done = 1;
        } /* while done */

        /* if the reaction had two sets of coefficients, need to be reversed */
        if (r->pressure == 2) {
          fprintf(stderr, "Warning: swapping FlameMaster coefficients\n");
          tmp = r->A_low;  r->A_low = r->A;  r->A = tmp;
          tmp = r->E_low;  r->E_low = r->E;  r->E = tmp;
          tmp = r->b_low;  r->b_low = r->b;  r->b = tmp;
        }
      }
    }
    else {
      /* should only be empty lines here */
    }

  } /* end loop over file */
  fclose(pF);

  int newSp;
  /* ----- process the reactions ----- */
  /* accumulate a list of species */
  nSpecies = nAddSp;
  species[0] = (char**) realloc(species[0], nSpecies*sizeof(char*));
  for (i=0; i<nAddSp; i++) {
    species[0][i] = (char*) calloc(32, sizeof(char));
    strcpy(species[0][i], addSp[i]);
  }

  for (iR=0; iR<nRxn; iR++) {
    r = &rxns[0][iR]; /* tmp for readability */

    /* check the reactants */
    for (j=0; j<r->nR; j++) {
      newSp = 1;
      for (iSp=0; iSp<nSpecies; iSp++) {
        if (!strcmp(r->spR[j], species[0][iSp])) {
          newSp = 0;
          break;
        }
      }
      if (newSp &&
          !(r->spR[j][0] == 'M' && (r->spR[j][1] == '\0' ||
                                    !isalpha(r->spR[j][1])) ) ) {
        nSpecies++;
        iSp = nSpecies - 1;
        species[0]      = (char**) realloc(species[0], nSpecies*sizeof(char*));
        species[0][iSp] = (char*)   calloc(32, sizeof(char));
        strcpy(species[0][iSp], r->spR[j]);
      }
      if (r->spR[j][0] == 'M' && (r->spR[j][1] == '\0' ||
                                  !isalpha(r->spR[j][1])) ) {
        for (k=0; k<nEnh; k++) {
          if (!strcmp(r->spR[j], enhM[k])) {
            if (r->pressure == 0) r->pressure = 1;
            r->nEnhance = nEnhSp[k];
            r->factEnh = (double*) calloc(nEnhSp[k], sizeof(double));
            r->spEnh   = (char**)  calloc(nEnhSp[k], sizeof(char*));
            for (iSp=0; iSp<nEnhSp[k]; iSp++) {
              r->factEnh[iSp] = enhF[k][iSp];
              r->spEnh[iSp] = (char*)  calloc(32, sizeof(char));
              strcpy(r->spEnh[iSp], enhSp[k][iSp]);
              stripString(r->spEnh[iSp]);
            }
          }
        } /* end loop over third-body species */
        r->spR[j][1] = '\0';
      }

    } /* end loop over reactants */

    /* check the products */
    for (j=0; j<r->nP; j++) {
      newSp = 1;
      for (iSp=0; iSp<nSpecies; iSp++) {
        if (!strcmp(r->spP[j], species[0][iSp])) {
          newSp = 0;
          break;
        }
      }
      if (newSp == 1 &&
          !(r->spP[j][0] == 'M' && (r->spP[j][1] == '\0' ||
                                    !isalpha(r->spP[j][1])) ) ) {
        nSpecies++;
        iSp = nSpecies - 1;
        species[0]      = (char**) realloc(species[0], nSpecies*sizeof(char*));
        species[0][iSp] = (char*)   calloc(32, sizeof(char));
        strcpy(species[0][iSp], r->spP[j]);
      }
      if ((r->spP[j][0] == 'M' && (r->spP[j][1] == '\0' ||
                                   !isalpha(r->spP[j][1])) ) ) {
        r->spP[j][1] = '\0';
      }

    } /* end loop over products */

    /* check the Troe parameters */
    if (r->troe == 1) {
      if (r->troeA[0] == 0 && r->troeA[1] > 0.0)
        r->troeA[0] = 1.0 - r->troeA[1];
      else if (r->troeA[1] == 0 && r->troeA[0] > 0.0)
        r->troeA[1] = 1.0 - r->troeA[0];
      if ((r->troeT[0] <= 1.0E-30 || r->troeT[0] >= 1.0E30) &&
          (r->troeT[1] <= 1.0E-30 || r->troeT[1] >= 1.0E30) &&
          (r->troeT[2] <= 1.0E-30 || r->troeT[2] >= 1.0E30) )
        r->troe = 0;
      if ((r->troeT[2] <= 1.0E-30 || r->troeT[2] >= 1.0E30))
        r->troeT[2] = 0.0;
    }

  } /* end loop over reactions */
  *nSpecs = nSpecies;

  /* make sure species are upper case */
  for (i=0; i<nSpecies; i++)
    for (j=0; j<strlen(species[0][i]); j++)
      species[0][i][j] = toupper(species[0][i][j]);

  for (i=0; i<nEnh; i++) {
    for (j=0; j<nEnhSp[i]; j++) free(enhSp[i][j]);
    free(enhM[i]); free(enhF[i]); free(enhSp[i]);
  }
  free(nEnhSp); free(enhM); free(enhF); free(enhSp);

  return nRxn;
}
/* -------------------------------------------------------------------------- */
