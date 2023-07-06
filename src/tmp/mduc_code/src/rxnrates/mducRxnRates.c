#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mducRxnRates.h"
#include "fileIO.h"
#include "chemInit.h"

struct ind_val {
   int ind;
   double val;

};

int sorting_ind_val (const void * a, const void * b)
{
   double lhs = ((struct ind_val*)a)->val;
   double rhs = ((struct ind_val*)b)->val;
  
   if (fabs(lhs)<fabs(rhs))
	return 1;

   if (fabs(lhs)>fabs(rhs))
	return -1;

   return 0;
}


void printHelp()
{
  printf("\n");
  printf("options :\n");
  printf("  -h                       # print help\n");
  exit(0);
}
/* -------------------------------------------------------------------------- */

void setDefaultOptions(optStruct *opts)
{
  opts->postFileName = NULL;
  strcpy(opts->suffix, "txt");
  strcpy(opts->outFile, "reactionRates");
}
/* -------------------------------------------------------------------------- */
void getSootSurfProducts(int *products, double *sootMass,chemStruct *chem);
void getSootRxn(char *sootRxn[]);
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  int i, j;
  optStruct opts;  // struct to store options read in from command line
  char file[128], *head, *tail;  // storage for beginning and end of string
  char buffer[128], prefix[128], suffix[20];
  double *w_int;
  FILE *pFile = NULL;

  
  // declare flamelet structs so that internal functions can be used
  mducStruct *mduc = (mducStruct*) calloc(1, sizeof(mducStruct));
  /* allocate pointers to internal structs */
  mduc->io     = (ioStruct*)     calloc(1, sizeof(ioStruct));
  mduc->driver = (driveStruct*)  calloc(1, sizeof(driveStruct));
  mduc->chem   = NULL;
  mduc->time   = (timeStruct*)   calloc(1, sizeof(timeStruct));
  mduc->flow   = (flowStruct*)   calloc(1, sizeof(flowStruct));
  mduc->flame  = (flameStruct*)  calloc(1, sizeof(flameStruct));
  mduc->sdial  = (sdialStruct*)  calloc(1, sizeof(sdialStruct));
  mduc->engine = (engineStruct*) calloc(1, sizeof(engineStruct));

  /* local pointers */
  ioStruct     *io     = mduc->io;
  driveStruct  *driver = mduc->driver;
  chemStruct   *chem   = mduc->chem;
  timeStruct   *time   = mduc->time;
  flowStruct   *flow   = mduc->flow;
  flameStruct  *flame  = mduc->flame;
  sdialStruct  *sdial  = mduc->sdial;
  engineStruct *engine = mduc->engine;

  // set default options for plot script
  setDefaultOptions(&opts);
  // interpret options read passed from command line
  readOptions(argc, argv, &opts);
  if (opts.postFileName[0] == NULL) {
    printf("Error: no input file specified\n");
    exit(2);
  }

  io->parser = parseFile(opts.postFileName[0]); 
  mduc->chem = chemInit("compiled",io->parser); /* variables from mechanism */
  chem = mduc->chem;
  // generate the tabulated files for plotting
  // TODO: Include soot reactions
  w_int    = (double*) calloc(chem->nRxnTotal,sizeof(double));
  for (j=0;j<chem->nRxnTotal;j++) w_int[j]=0.0;
  fprintf(stderr, "Processing %i files:\n", opts.nFiles);
  for (i=0; i<opts.nFiles; i++) {
    fprintf(stderr, "\t%s", opts.postFileName[i]);
    
    // read in rates file
    readRates(opts.postFileName[i],chem);

    for (j=0;j<chem->nRxnTotal;j++) w_int[j]+=chem->w_int_pre[j];
    
  }
  // write the tabulated file for plotting
  for (j=0;j<chem->nRxnTotal;j++) chem->w_int[j]=w_int[j];
  
  writePostFile(opts.outFile, &opts, chem);
    
  fprintf(stderr, "\n");
  fprintf(stderr, "%i files processed\n", opts.nFiles);
}
/* -------------------------------------------------------------------------- */

void createOutputFileName(char *prefix, char *fileName, char *ext,
                          char *outputFile)
{
  strcpy(outputFile, prefix);
  strcat(outputFile, "_");
  strcat(outputFile, fileName);
  strcat(outputFile, ".");
  strcat(outputFile, ext);
}
/* -------------------------------------------------------------------------- */

/* Function takes a flamelet solution and writes a file in tabulated form */
void writePostFile(char *fileName, optStruct *opts, chemStruct *chem)
{
  int i, j, r, s, st, it, first,rs,count;
  char file_rel[100],file_abs[100], string[200];
  double sum,tmp;
  double *sumPlus,*sumMinus;
  int *products,*educts,*thirdbody,*backward,*nuF;

  // add suffix
  strcpy(file_abs, fileName);
  strcat(file_abs, ".abs");
  fprintf(stderr, " --> %s ...", file_abs);

  /* ----- write data file ----- */
  FILE *pFile_abs = fopen(file_abs, "w");
  // add suffix
  strcpy(file_rel, fileName);
  strcat(file_rel, ".rel");
  fprintf(stderr, " & %s ...", file_rel);

  /* ----- write data file ----- */
  FILE *pFile_rel = fopen(file_rel, "w");
 
  /* ----- write body ----- */ 
 
  sumPlus         = (double*) calloc(chem->nSpecSolve,sizeof(double));
  sumMinus         = (double*) calloc(chem->nSpecSolve,sizeof(double));
  educts      = (int*) calloc(chem->nSpecSolve*chem->nRxn,sizeof(int));
  products    = (int*) calloc(chem->nSpecSolve*chem->nRxn,sizeof(int));
  nuF       = (int*) calloc(chem->nSpecSolve*chem->nRxn,sizeof(int));
  thirdbody   = (int*) calloc(chem->nRxn,sizeof(int));
  backward    = (int*) calloc(chem->nRxn,sizeof(int));

  /* get reaction coefficient from mechanism */
  GetReactionNu(products);
  GetReactionNuForward(nuF);

  for (s=0;s<chem->nSpecSolve;s++) {
     for (r=0;r<chem->nRxn;r++) { 
         i = s*chem->nRxn + r;

         products[i] *= -1;

         if (products[i]<0) {
            educts[i]   = -1 * products[i];
 	    products[i] = 0;
         }

         if (nuF[i]<0 && nuF[i] != -1*products[i])  {
	 //   fprintf(stdout,"Reaction: %i Species: %i\n",r,s);
	    educts[i] -= nuF[i];
	    products[i] -= nuF[i];  
         }
     }
  }

  /* get third body reactions */
  GetReactionTB(thirdbody);

  /* get backward reactions */
  GetBackwardReactions(backward);

  /* sort reaction rates */
  struct ind_val sortedRates[chem->nRxn];

  for (r=0;r<chem->nRxn;r++) {
      sortedRates[r].ind = r;
      if (backward[r]>0)
	 sortedRates[r].val = chem->w_int[r]-chem->w_int[backward[r]];
      else
	 sortedRates[r].val = chem->w_int[r];
  }
	
  qsort(sortedRates,chem->nRxn,sizeof(struct ind_val), sorting_ind_val);

  /* species */
      sumPlus[s] = 0.0;
      sumMinus[s] = 0.0;
  for (s=0;s<chem->nSpecSolve;s++) {
      fprintf(pFile_abs,"=================== %s ===================\n",chem->species[s]);
      sum = 0.0;
      for (rs=0;rs<chem->nRxn;rs++) { 
          r = sortedRates[rs].ind;
          i = s*chem->nRxn + r;

          tmp = 0.0;
          count = 0;
          if (backward[r]>=0) {
          if (educts[i] - products[i] != 0) {
             /* print this reaction! */
             //fprintf(pFile_abs,"Reaction %5i: %13.6e kmol : ",r, sortedRates[rs].val);
             fprintf(pFile_abs,"Reaction %5i\t: ",r);
             /* print all educts first */ 
             first = 1;
             for (st=0;st<chem->nSpecSolve;st++) {
                 it = st*chem->nRxn + r;

                 if (educts[it] != 0) {
                    if (first==1) {
 			first = 0;
	                //fprintf(pFile_abs,"%1i %s",educts[it],chem->species[st]);
	                if (educts[it]>1)  
	                  count = sprintf(string,"%1i %s",educts[it],chem->species[st]);
			else
	                  count = sprintf(string,"%s",chem->species[st]);
                    }
                    else {
	                //fprintf(pFile_abs,"+ %1i %s",educts[it],chem->species[st]);
	                if (educts[it]>1)  
	                   count += sprintf(string+count," + %1i %s",educts[it],chem->species[st]);
			else
	                   count += sprintf(string+count," + %s",chem->species[st]);
                    }
                    if (s==st) {
			if (sortedRates[rs].val>0) 
    			   sumMinus[s] -= (double)educts[it]*sortedRates[rs].val;
		        else
    			   sumPlus[s] -= (double)educts[it]*sortedRates[rs].val;

 			sum -= (double)educts[it]*sortedRates[rs].val;
 			tmp -= (double)educts[it]*sortedRates[rs].val;
                    }
                 }
             }
 
             /* thirdbody ? */
             if (thirdbody[r] != 0) { 
                //fprintf(pFile_abs,"%s", "(+M)");
                count += sprintf(string+count,"%s", " (+M)");
             }

             //fprintf(pFile_abs,"<=> ");
             if (backward[r]>0) 
             count += sprintf(string + count," <=> ");
             else
                count += sprintf(string + count," => ");
             /* print all products */ 
             first = 1;
             for (st=0;st<chem->nSpecSolve;st++) {
                 it = st*chem->nRxn + r;

                 if (products[it] != 0) {
                    if (first==1) {
 			first = 0;
	                //fprintf(pFile_abs,"%1i %s",products[it],chem->species[st]);
	                if (products[it]>1) 
    	                   count += sprintf(string+count,"%1i %s",products[it],chem->species[st]);
			else
    	                   count += sprintf(string+count,"%s",chem->species[st]);
                    }
                    else {
	                //fprintf(pFile_abs,"+ %1i %s",products[it],chem->species[st]);
	                if (products[it]>1) 
	                   count += sprintf(string+count," + %1i %s",products[it],chem->species[st]);
 			else
	                   count += sprintf(string+count," + %s",chem->species[st]);
                    }
                    if (s==st) {
			if (sortedRates[rs].val>0) 
    			   sumPlus[s] += (double)products[it]*sortedRates[rs].val;
		        else
    			   sumMinus[s] += (double)products[it]*sortedRates[rs].val;

                        sum += (double)products[it]*sortedRates[rs].val;
                        tmp += (double)products[it]*sortedRates[rs].val;
                    }
                 }
             }

             /* thirdbody ? */
             if (thirdbody[r] != 0) {
             //   fprintf(pFile_abs,"%s", "(+M)");
                count += sprintf(string+count,"%s", " (+M)");
 	     }

             fprintf(pFile_abs,"%13.6e kg : %s\n",tmp*chem->W[s],string);
           } 
           } /* not a backward reaction */
      } /* end loop over reactions */
      fprintf(pFile_abs,"Overall production of %s : %15.6e kg\n",chem->species[s], sum*chem->W[s]);
      fprintf(pFile_abs,"\n");
  } /* end loop over species */
  
  
  /* ----- copy molecular weight of soot rxn species ----- */
  double *sootSpecW   = (double*) calloc(chem->nSootSp,sizeof(double));
  //for (s=0; s<chem->nSootSp; s++)
  //  for (st=0; st<chem->nSpecies; st++)
  //    if (!strcmp(chem->species[st], chem->sootRxnSpec[s]))
  //      sootSpecW[s] = chem->W[st];
  fprintf(stderr,"Warning!: PAH molecular weights hardcoded for CS paper\n");
  sootSpecW[0] = 16.0*12.0; //A4XC16H10
  sootSpecW[1] = 10.0*12.0; //A2XC10H8
  sootSpecW[2] = 14.0*12.0; //A3XC14H10
  sootSpecW[3] = 12.0*12.0; //A2R5
  sootSpecW[4] = 12.0*12.0; //P2XC12H10
  sootSpecW[5] = 16.0*12.0; //A3R5
  sootSpecW[6] = 16.0*12.0; //FLTN
  sootSpecW[7] = 18.0*12.0; //A4R5
  /* ----- write body soot PAH reactions ----- */

  fprintf(pFile_abs,"\n\n=================== Soot PAH Reactions ===================\n");
  double sum_dimer=0.0;
  double sootPlus,sootMinus=0.0;
  for (s=0; s<chem->nSootPAH; s++) {
    fprintf(pFile_abs,"\n=================== %s ===================\n", chem->sootPAHSpec[s]);
    double tmp = -2.0*chem->w_int[chem->nRxn+s]*sootSpecW[s];
    sum_dimer += tmp;
    fprintf(pFile_abs,"Reaction %5i\t: %13.6e kg : %s + %s => Dimer \n",chem->nRxn + s, tmp, chem->sootPAHSpec[s], chem->sootPAHSpec[s]);
    fprintf(pFile_abs,"Overall production of %s : \t %13.6e kg : \n", chem->sootPAHSpec[s], tmp);
  }/* end soot PAH species*/
  sootPlus = -sum_dimer;
  /*---- soot surface reactions ----*/
  int nSurfaceRxn = chem->nSootRxn-chem->nSootPAH-1;
  int sootSurfaceSp = chem->nSootSp - chem->nSootPAH; // soot rxn species without PAH
  
  /* ----- initialize the products ----- */
  int *sootProducts    = (int*) calloc(chem->nSootSp*nSurfaceRxn,sizeof(int));
  double *sootMass     = (double*) calloc(nSurfaceRxn,sizeof(double));
  getSootSurfProducts(sootProducts, sootMass, chem);
  
  char *sootRxn[nSurfaceRxn];
  getSootRxn(sootRxn);
  
  /* ----- write body soot surface reactions ----- */
  fprintf(pFile_abs,"=================== Soot Surface Reactions ===================\n");
  
  int sp, sRxn, currRxn;
  double sumSootSp, tmpSootSp;
  for (sp=0; sp<sootSurfaceSp; sp++ ){
      fprintf(pFile_abs,"\n=================== %s ===================\n",chem->sootRxnSpec[sp + chem->nSootPAH]);
      sumSootSp = 0.0;
      for(sRxn = 0; sRxn < nSurfaceRxn; sRxn++){
        tmpSootSp = 0.0;
        if(sootProducts[sp*nSurfaceRxn + sRxn] !=0){
          currRxn = sRxn + chem->nRxn + chem->nSootPAH;
          tmpSootSp = sootProducts[sp*nSurfaceRxn + sRxn]*chem->w_int[currRxn]*sootSpecW[sp + chem->nSootPAH];
          fprintf(pFile_abs,"Reaction %5i\t: %13.6e kg : %s \n", currRxn, tmpSootSp, sootRxn[sRxn]);
          sumSootSp += tmpSootSp;
        }
      }
      fprintf(pFile_abs,"Overall production of %s : %15.6e kg\n",chem->sootRxnSpec[sp + chem->nSootPAH], sumSootSp );
  }

  fprintf(pFile_abs,"\n=================== %s ===================\n","Soot & Dimer mass");

     for (s=0; s<chem->nSootPAH; s++) {
       double tmp = 2.0*chem->w_int[chem->nRxn+s]*sootSpecW[s];
       fprintf(pFile_abs,"Reaction %5i\t: %13.6e kg : %s + %s => Dimer \n",chem->nRxn + s, tmp, chem->sootPAHSpec[s], chem->sootPAHSpec[s]);
     }/* end soot PAH species*/
  for(sRxn = 0; sRxn < nSurfaceRxn; sRxn++){
     currRxn = sRxn + chem->nRxn + chem->nSootPAH;
     tmpSootSp = chem->w_int[currRxn]*sootMass[sRxn];
     if (tmpSootSp>0.0) sootPlus += tmpSootSp;
     if (tmpSootSp<0.0) sootMinus += tmpSootSp;
     fprintf(pFile_abs,"Reaction %5i\t: %13.6e kg : %s \n", currRxn, tmpSootSp, sootRxn[sRxn]);
  }
  fprintf(pFile_abs,"Overall soot and dimer production : %15.6e kg\n", sootPlus + sootMinus );
  fclose(pFile_abs);
  /* ---- write relative file ---- */
  /* species */
  for (s=0;s<chem->nSpecSolve;s++) {
      fprintf(pFile_rel,"=================== %s ===================\n",chem->species[s]);
      for (rs=0;rs<chem->nRxn;rs++) { 
          r = sortedRates[rs].ind;
          i = s*chem->nRxn + r;

          tmp = 0.0;
          count = 0;
          if (backward[r]>=0) {
          if (educts[i] - products[i] != 0) {
             /* print this reaction! */
             /* print all educts first */ 
             first = 1;
             for (st=0;st<chem->nSpecSolve;st++) {
                 it = st*chem->nRxn + r;

                 if (educts[it] != 0) {
                    if (first==1) {
 			first = 0;
	                if (educts[it]>1)  
	                  count = sprintf(string,"%1i %s",educts[it],chem->species[st]);
			else
	                  count = sprintf(string,"%s",chem->species[st]);
                    }
                    else {
	                if (educts[it]>1)  
	                   count += sprintf(string+count," + %1i %s",educts[it],chem->species[st]);
			else
	                   count += sprintf(string+count," + %s",chem->species[st]);
                    }
                    if (s==st) {
 			tmp -= (double)educts[it]*sortedRates[rs].val;
                    }
                 }
             }
 
             /* thirdbody ? */
             if (thirdbody[r] != 0) { 
                count += sprintf(string+count,"%s", " (+M)");
             }

             //count += sprintf(string + count," <=> ");
             if (backward[r]>0) 
                count += sprintf(string + count," <=> ");
             else
                count += sprintf(string + count," => ");
             /* print all products */ 
             first = 1;
             for (st=0;st<chem->nSpecSolve;st++) {
                 it = st*chem->nRxn + r;

                 if (products[it] != 0) {
                    if (first==1) {
 			first = 0;
	                if (products[it]>1) 
    	                   count += sprintf(string+count,"%1i %s",products[it],chem->species[st]);
			else
    	                   count += sprintf(string+count,"%s",chem->species[st]);
                    }
                    else {
	                if (products[it]>1) 
	                   count += sprintf(string+count," + %1i %s",products[it],chem->species[st]);
 			else
	                   count += sprintf(string+count," + %s",chem->species[st]);
                    }
                    if (s==st) {
                        tmp += (double)products[it]*sortedRates[rs].val;
                    }
                 }
             }

             /* thirdbody ? */
             if (thirdbody[r] != 0) {
                count += sprintf(string+count,"%s", " (+M)");
 	     }
             if (tmp/sumPlus[s]>=0.0 && sumPlus[s]>0.0) {
                fprintf(pFile_rel,"Reaction %5i\t: ",r);
                fprintf(pFile_rel,"%6.2f [%%] : %s\n",tmp/sumPlus[s]*100,string);
             }
             //if (tmp/sumMinus[s]>0.01 && sumMinus[s]<0.0) {
             //   fprintf(pFile_rel,"Reaction %5i\t: ",r);
             //   fprintf(pFile_rel,"%13.6e [%%] : %s\n",-tmp/sumMinus[s]*100,string);
             //}
           } 
           } /* not a backward reaction */
      } /* end loop over reactions */
      for (rs=0;rs<chem->nRxn;rs++) { 
          r = sortedRates[rs].ind;
          i = s*chem->nRxn + r;

          tmp = 0.0;
          count = 0;
          if (backward[r]>=0) {
          if (educts[i] != 0 || products[i] != 0) {
             /* print this reaction! */
             /* print all educts first */ 
             first = 1;
             for (st=0;st<chem->nSpecSolve;st++) {
                 it = st*chem->nRxn + r;

                 if (educts[it] != 0) {
                    if (first==1) {
 			first = 0;
	                if (educts[it]>1)  
	                  count = sprintf(string,"%1i %s",educts[it],chem->species[st]);
			else
	                  count = sprintf(string,"%s",chem->species[st]);
                    }
                    else {
	                if (educts[it]>1)  
	                   count += sprintf(string+count," + %1i %s",educts[it],chem->species[st]);
			else
	                   count += sprintf(string+count," + %s",chem->species[st]);
                    }
                    if (s==st) {
 			tmp -= (double)educts[it]*sortedRates[rs].val;
                    }
                 }
             }
 
             /* thirdbody ? */
             if (thirdbody[r] != 0) { 
                count += sprintf(string+count,"%s", " (+M)");
  }

             count += sprintf(string + count," <=> ");
             /* print all products */ 
             first = 1;
             for (st=0;st<chem->nSpecSolve;st++) {
                 it = st*chem->nRxn + r;

                 if (products[it] != 0) {
                    if (first==1) {
 			first = 0;
	                if (products[it]>1) 
    	                   count += sprintf(string+count,"%1i %s",products[it],chem->species[st]);
			else
    	                   count += sprintf(string+count,"%s",chem->species[st]);
                    }
                    else {
	                if (products[it]>1) 
	                   count += sprintf(string+count," + %1i %s",products[it],chem->species[st]);
 			else
	                   count += sprintf(string+count," + %s",chem->species[st]);
                    }
                    if (s==st) {
                        tmp += (double)products[it]*sortedRates[rs].val;
                    }
                 }
             }

             /* thirdbody ? */
             if (thirdbody[r] != 0) {
                count += sprintf(string+count,"%s", " (+M)");
 	     }
             //if (tmp/sumPlus[s]>0.01 && sumPlus[s]>0.0) {
             //   fprintf(pFile_rel,"Reaction %5i\t: ",r);
             //   fprintf(pFile_rel,"%13.6e [%%] : %s\n",tmp/sumPlus[s]*100,string);
             //}
             if (tmp/sumMinus[s]>=0.0 && sumMinus[s]<0.0) {
                fprintf(pFile_rel,"Reaction %5i\t: ",r);
                fprintf(pFile_rel,"%6.2f [%%] : %s\n",-tmp/sumMinus[s]*100,string);
             }
           } 
           } /* not a backward reaction */
      } /* end loop over reactions */
      fprintf(pFile_rel,"Total production :  %8.4g [mug]\n",sumPlus[s]*chem->W[s]*1.0e9);
      fprintf(pFile_rel,"Total consumption : %8.4g [mug]\n",sumMinus[s]*chem->W[s]*1.0e9);
      fprintf(pFile_rel,"\n");
  } /* end loop over species */
  /* ----- write body soot PAH reactions ----- */
  if (sum_dimer<0.0) {
     fprintf(pFile_rel,"\n=================== %s ===================\n", "Soot & Dimer Mass");

     for (s=0; s<chem->nSootPAH; s++) {
       double tmp = 2.0*chem->w_int[chem->nRxn+s]*sootSpecW[s];
       fprintf(pFile_rel,"Reaction %5i\t: %6.2f [%%] : %s + %s => Dimer \n",chem->nRxn + s, tmp/sootPlus*100.0, chem->sootPAHSpec[s], chem->sootPAHSpec[s]);
     }/* end soot PAH species*/
  }
  
  for(sRxn = 0; sRxn < nSurfaceRxn; sRxn++){
     currRxn = sRxn + chem->nRxn + chem->nSootPAH;
     tmpSootSp = chem->w_int[currRxn]*sootMass[sRxn];
     if (tmpSootSp>0.0 && sootPlus > 0.0)
        fprintf(pFile_rel,"Reaction %5i\t: %6.2f [%%] : %s \n", currRxn, tmpSootSp/sootPlus*100.0, sootRxn[sRxn]);
     if (tmpSootSp<0.0 && sootMinus < 0.0)
        fprintf(pFile_rel,"Reaction %5i\t: %6.2f [%%] : %s \n", currRxn, tmpSootSp/sootMinus*-100.0, sootRxn[sRxn]);
  }
  fprintf(pFile_rel,"Total production : %6.2e [kg]\n", sootPlus);
  fprintf(pFile_rel,"Total consumption : %6.2e [kg]\n",sootMinus);
  fclose(pFile_rel);
  free(educts); 
  free(products);
  free(thirdbody);
  free(sumPlus);
  free(sumMinus);
  /* ----- free soot var ----- */
  free(sootSpecW);
  //free(sootRxn); //?
}
/* -------------------------------------------------------------------------- */
void getSootSurfProducts(int *products, double *sootMass, chemStruct *chem){

  int nSurfaceRxn = chem->nSootRxn-chem->nSootPAH;
  int sootSurfaceSp = chem->nSootSp - chem->nSootPAH; // soot rxn species without PAH
  int sp, sRxn;
  for (sp=0; sp<sootSurfaceSp; sp++)
    for(sRxn = 0; sRxn < nSurfaceRxn; sRxn++)
      products[sp*nSurfaceRxn + sRxn] = 0;

  /* ----- OH ----- */
  products[0*nSurfaceRxn+0] = -1;
  products[0*nSurfaceRxn+5] = -1;

  /* ----- H2O ----- */
  products[1*nSurfaceRxn+0] = 1;

  /* ----- H ----- */
  products[2*nSurfaceRxn+1] = -1;
  products[2*nSurfaceRxn+2] = 1;

  /* ----- H2 ----- */
  products[3*nSurfaceRxn+1] = 1;
  
  /* ----- C2H2 ----- */
  products[4*nSurfaceRxn+3] = -1;
  
  /* ----- O2 ----- */
  products[5*nSurfaceRxn+4] = -1;

  /* ----- CO ----- */
  products[6*nSurfaceRxn+4] = 2;
  products[6*nSurfaceRxn+5] = 1;

  for(sRxn = 0; sRxn < nSurfaceRxn; sRxn++)
    sootMass[sRxn] = 0.0;

  //sootMass[0] = -1.0 * 1.0; /* soot mass loss due to H-abstraction */ 
  //sootMass[1] = -1.0 * 1.0; /* soot mass loss due to H-abstraction */ 
  //sootMass[2] = -1.0 * 1.0; /* soot mass loss due to H-abstraction */ 
  //sootMass[3] =  1.0 * (12.0*2.0+2.0*1.0); /* soot mass gain due to C2H2 addtion */
  //sootMass[4] = -2.0 * 12.0; /* soot mass los due to oxidation with O2 */
  //sootMass[5] = -1.0 * 12.0 + 1.0*1.0 ; /* soot mass loss due carbon-abstraction with H-addition */ 
  sootMass[0] = 0.0; /* soot mass loss due to H-abstraction */ 
  sootMass[1] = 0.0; /* soot mass loss due to H-abstraction */ 
  sootMass[2] = 0.0; /* soot mass loss due to H-abstraction */ 
  sootMass[3] =  1.0 * (12.0*2.0); /* soot mass gain due to C2H2 addtion */
  sootMass[4] = -2.0 * 12.0; /* soot mass los due to oxidation with O2 */
  sootMass[5] = -1.0 * 12.0; /* soot mass loss due carbon-abstraction with H-addition */ 
}
/* -------------------------------------------------------------------------- */
void getSootRxn(char *sootRxn[]){
  sootRxn[0] = "Soot-H + OH   <=> Soot* + H20";
  sootRxn[1] = "Soot-H + H    <=> Soot-* + H2";
  sootRxn[2] = "Soot-H        <=> Soot-* + H ";
  sootRxn[3] = "Soot-* + C2H2  => Soot-H";
  sootRxn[4] = "Soot-* + O2    => Soot-* + 2CO";
  sootRxn[5] = "Soot-H + OH    => Soot-H + CO";
}
/* -------------------------------------------------------------------------- */

void readList(int n, int iStart, char ***var, char **argv)
{
  int i, j;
  
  var[0] = (char**) calloc(n, sizeof(char*));
  for (i=0, j=iStart; i<n; i++, j++) {
    var[0][i] = (char*) calloc(strlen(argv[j])+10, sizeof(char));
    strcpy(var[0][i], argv[j]);
  }
}
/* -------------------------------------------------------------------------- */

void readOptions(int argc, char **argv, optStruct *opts)
{
  int i, j, k, nOpt, iOpt, iEnd;
  int *optLocs;
  
  // read through the argument list and get the option locations
  optLocs    = (int*)  calloc(argc, sizeof(int));
  opts->args = (char*) calloc(argc, sizeof(char));
  nOpt = 0;
  for (i=1; i<argc; i++)
    {
      if (argv[i][0] == '-')
        {
          optLocs[nOpt] = i;
          opts->args[nOpt] = argv[i][1];
          nOpt++;
        }
    } // end loop over arguments
  opts->nOpts = nOpt;
  
  for (i=0; i<nOpt; i++)
    {
      iOpt = optLocs[i];
      
      if (i == nOpt-1)
        iEnd = argc-1;
      else
        iEnd = optLocs[i+1]-1;
      
      switch (argv[optLocs[i]][1])
        {
        case 'h':
          printHelp();
          break;
        case 'i':
          {
            if (argv[iOpt+1][0] == '-')
              {
                printf("ERROR: definition of input file required for '-i'\n");
                exit(2);
              }
            else
              {
                opts->nFiles = iEnd - iOpt;
                readList(opts->nFiles, iOpt+1, &opts->postFileName, argv);
              }
            break;
          }
        case 's':
          {
	    /*TODO: Read in PAH species considered in soot model */            
            opts->nPAH = iEnd - iOpt;
            readList(opts->nPAH, iOpt+1, &opts->PAHlist, argv);
            break;
          }
        }
        
    } // end loop over options
  
  free(optLocs);
}
/* -------------------------------------------------------------------------- */



