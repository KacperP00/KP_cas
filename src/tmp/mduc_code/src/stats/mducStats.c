#include "mducStats.h"
#include "statsTools.h"
#include "accessMDUC.h"
#include "fileIO.h"
#include "init.h"
#include "statsUtils.h"
#include "utils.h"

void printHelp()
{
  printf("\n");
  printf("options :\n");
  printf("  -h                       # print help\n");
  printf("  -i <string>              # input file\n");
  printf("  -d <dimension>           # number of scalars\n");
  printf("  -M <string> <string> ... # method(s) to use\n");
  printf("  -t <string>              # transformation to use\n");
  printf("  -s                       # scale pdf according to grid\n");
  printf("  -T                       # run verification test cases\n");
  printf("  -E                       # run error computation\n");
  printf("  -R                       # generate regime plot\n");
  printf("  -L                       # compute mean w/out linear assumption\n");
  printf("  -I <string>              # output file identifier\n");
  printf("options :\n");
  printf("  -n <n1> <n2>             # number of grid points\n");
  printf("  -g <pivot> <perc>        # simple non-uniform grid parameters\n");
  printf("set mean(s) and variance(s):\n");
  printf("  -m <mean1> <mean2> ...   # mean(s) for each\n");
  printf("  -v <var1> <var2> ...     # variance(s)\n");
  printf("set beta parameters directly:\n");
  printf("  -b <beta1> <beta2> ...   # supply betas directly\n");
  printf("extract PDF of variable from CFD data file:\n");
  printf("  -e <variable name>       # PDF variable to generate\n");
  printf("  --mfile <string>         # moment file\n");
  printf("  --prefix <string>        # file prefix for input/ouput\n");
  printf("  --range <int> <int>      # time step range in file\n");
  printf("  --set <int>              # data set within input file\n");
  exit(0);
}
/* -------------------------------------------------------------------------- */

void setDefaultOpts(optStruct *opts)
{
  opts->readOptFile = 0;
  opts->method = NULL;
  strcpy(opts->ident, "");
  strcpy(opts->transformName, "doran");
  strcpy(opts->errorFile, "stats.err");
  opts->dim       = 1;
  opts->runTest   = 1;
  opts->errorTest = 0;
  opts->n1        = 101;
  opts->n2        = 101;
  opts->mean1     = opts->mean2 = -1;
  opts->var1      = opts->var1  = opts->var12 = -1;
  opts->beta1     = opts->beta2 = opts->beta3 = -1;
  opts->inputMean = 0;
  opts->inputBeta = 0;
  opts->transform = 1;
  opts->scaled    = 0;
  opts->verbose   = 0;
  opts->regime    = 0;
  opts->meanLin   = 1;
  opts->extract   = 0;
  opts->interpFile = 0;
  opts->range[0] = opts->range[1] = -1;
  opts->inputFile = 0;
  opts->dataSet   = 0;
  opts->nPiv1 = opts->nPiv2 = 0;
  opts->gPiv1 = opts->gPiv2 = opts->gPerc1 = opts->gPerc2 = NULL;
}
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  int i, j, k, m, f, nIter = 1;
  int n1, n2;
  int nx, ny, nz;
  char dataFile[100], interpFile[100], buffer[100];
  double tmp;

  optStruct opts;

  FILE *pFile   = NULL; // general file pointer
  FILE *errFile = NULL; // error file pointer

  // get the options from the command line
  setDefaultOpts(&opts);
  readOptions(argc, argv, &opts);

  if (opts.readOptFile)
    readOptionsFile(&opts);

  // set default method if not read or if set to all
  if (opts.method == NULL)
    {
      opts.nMethod   = 1;
      opts.method    = (char**) calloc(opts.nMethod, sizeof(char*));
      opts.method[0] = (char*)  calloc(50, sizeof(char));
      strcpy(opts.method[0], "beta");
    }

  if (opts.interpFile)
    {
      opts.nMethod++;
      opts.method    = (char**) realloc(opts.method, opts.nMethod*sizeof(char*));
      opts.method[opts.nMethod-1] = (char*)  calloc(50, sizeof(char));
      strcpy(opts.method[opts.nMethod-1], "interp");
    }

  if (!strncmp(opts.method[0], "independ", 8))
    opts.transform = 0;

  /* ----- 1D PDF ----- */
  if (opts.dim == 1)
    {
      if (opts.runTest)
        {
          printf("Running test cases for 1D pdf...\n");
          runTest1D(&opts);
          return 0;
        }

      if (opts.errorTest)
        {
          printf("Running error test for 1D pdf...\n");
          errorTest1D(&opts);
          return 0;
        }

      if (opts.regime)
        {
          printf("Computing regimes for 1D pdf...\n");
          regime(&opts);
          return 0;
        }

      int    pdfFlag, pdfBound[2];
      double *x, *dx; // mesh parameters
      double *Y, *pdf;
      double mean, var, meanPDF, varPDF;
      double dx0, eps;
      double error, errV;

      // variables required if extracting a PDF
      double ***S;
      double xmin = 0.0;
      double xmax = 1.0;

      n1 = opts.n1;

      if (opts.inputMean)
        {
          mean = opts.mean1;
          var  = opts.var1;
          tmp   = mean*(1.0-mean)/var - 1.0;
          opts.beta1 = mean*tmp;
          opts.beta2 = (1.0-mean)*tmp;
        }
      else if (opts.inputBeta)
        {
          mean = opts.beta1/(opts.beta1+opts.beta2);
          var  = opts.beta1*opts.beta2/
            ( pow(opts.beta1+opts.beta2, 2)*(opts.beta1+opts.beta2+1) );
        }
      else if (opts.extract)
        {
          printf("Extracting 1D pdf of %s from %s...\n", opts.names[0],
                 opts.file);

          strcpy(opts.method[0], "bin");

          // set pdf properties for post-processing
          pdfBound[0] = 0;
          pdfBound[1] = n1-1;

          // get data from file
          NGA_readData(opts.file, &nx, &ny, &nz, opts.names[0], &S, &tmp, &tmp);

          // get min and max values of the variable
          xmin = getArrayValueStr("min", nx, ny, nz, S);
          xmax = getArrayValueStr("max", nx, ny, nz, S);
          printf("S (min,max) = (%8.6f, %8.6f)\n", xmin, xmax);
        }
      else
        {
          printf("Error: please specify  mean/variance or beta/alpha combo\n");
          exit(2);
        }

      // allocate memory
      pdf = (double*) calloc(n1, sizeof(double));
      initGrid(n1, xmin, xmax, opts.nPiv1, opts.gPiv1, opts.gPerc1, &x, &dx);

      // generate the solution (ie Y = x), to compute the mean
      Y = (double*) calloc(n1, sizeof(double));
      for (j=0; j<n1; j++)
        Y[j] = x[j];

      strcpy(dataFile, "mduc_pdf");
      strcat(dataFile, opts.ident);
      strcat(dataFile, ".stat");

      // write header for file
      writeStatsHeader(dataFile, 1, &n1, &x, 1, &mean, &var, opts.nMethod,
                       opts.transformName);

      for (m=0; m<opts.nMethod; m++)
        {
          // compute the pdf
          if (opts.extract)
            pdfFlag = binPDF_1D(nx, ny, nz, S, n1, x, pdf);
          else
            pdfFlag = MDUCcalcPDF_1D(opts.method[m],n1,x,mean,var,pdf,pdfBound);

          sprintf(buffer, "pdf-%s", opts.method[m]);
          addVarToBody(dataFile, n1, 1, pdf, buffer);

          // compute the mean of the pdf
          if (opts.meanLin == 1)
            meanPDF = calcLinMeanY(0, 1, opts.beta1, opts.beta2, n1, x, Y,
                                   pdf, pdfBound, pdfFlag);
          else
            {
              Y[0] = Y[0] + 0.5*(x[1]-x[0]);
              Y[n1-1] = Y[n1-1] - 0.5*(x[n1-1] - x[n1-2]);
              meanPDF = calcMean(pdfBound[0], pdfBound[1], 0, 1, pdf, Y);
            }

          // compute the variance
          varPDF = 0.0;
          for (j=0; j<n1; j++)
            varPDF += (x[j]-mean)*(x[j]-mean)*pdf[j];

          // append error statistics to data file
          pFile = fopen(dataFile, "a");
          printStatsError(pFile,NULL,opts.method[m],1,&mean,&var,&meanPDF,
                          &varPDF, NULL);
          fclose(pFile);

          // print error to screen
          printStatsError(stderr,NULL,opts.method[m],1,&mean,&var,&meanPDF,
                          &varPDF, NULL);
        }
          
      free(Y);
      free(pdf);
      free(x);
      free(dx);
      return 0;
    }
  else if (opts.dim == 2) // 2D pdf
    {
      if (opts.errorTest)
        {
          errorTest2D(&opts);
          return 0;
        }

      if (opts.regime)
        {
          regime2D(&opts);
          return 0;
        }

      // grid variables
      int n[2], row, col, jj, iMom, fileExist;
      double **x, **dx, *x1, *x2;
      double x1min, x1max, x2min, x2max;

      // moment information
      double mean[3], var[4];       // moments
      double meanIn[3], varIn[4];   // moments that were input
      double meanPDF[3], varPDF[4]; // computed moments
      double **meanF, **varF;       // moments read from file
      int    *iter = NULL;
      double *time = NULL;
      double tFile;
      double fDiv[5];               // f-divergence

      // variables for dirichlet distribution
      double beta[4];  // beta parameters
      double varAn[4]; // analytic variance computed from betas
      double S, Q;     // temporary storage variables

      // pdf variables
      double **Y;
      int    **pdfBound, pdfFlag = 0;
      double **pdf, *pdf_tmp, **pdf_dns;

      // for CFD data
      double ***S1, ***S2;

      x1min = x2min = 0.0;
      x1max = x2max = 1.0;

      // process input scenarios
      if (opts.inputMean)
        {
          mean[1] = meanIn[1] = opts.mean1;
          mean[2] = meanIn[2] = opts.mean2;
          mean[0] = meanIn[0] = mean[1] + mean[2]; // total mean
          var[1]  = varIn[1]  = opts.var1;
          var[2]  = varIn[2]  = opts.var2;
          var[0]  = varIn[0]  = opts.var12;
        }
      else if (opts.inputBeta)
        {
          beta[0] = opts.beta1+opts.beta2+opts.beta3;
          // compute the means
          mean[1] = meanIn[1] = opts.beta1/beta[0];
          mean[2] = meanIn[2] = opts.beta2/beta[0];
          mean[0] = meanIn[0] = mean[1] + mean[2];

          // compute the variances
          tmp    = beta[0]*beta[0]*(beta[0]+1);
          var[1] = varIn[1] = opts.beta1*(beta[0] - opts.beta1)/tmp;
          var[2] = varIn[2] = opts.beta2*(beta[0] - opts.beta2)/tmp;
          var[0] = varIn[0] = opts.beta3*(beta[0] - opts.beta3)/tmp;
        }
      else if (opts.extract)
        {
          opts.inputFile = 0;
          printf("Extracting joint pdf of %s and %s from %s...\n",
                 opts.names[0], opts.names[1], opts.file);

          strcpy(opts.method[0], "bin");

          // get data from file
          NGA_readData(opts.file, &nx, &ny, &nz, opts.names[0], &S1, &tmp, &tmp);
          NGA_readData(opts.file, &nx, &ny, &nz, opts.names[1], &S2, &tmp, &tmp);

          // get min and max values of the variable
          x1min = getArrayValueStr("min", nx, ny, nz, S1);
          x1max = getArrayValueStr("max", nx, ny, nz, S1);
          x2min = getArrayValueStr("min", nx, ny, nz, S2);
          x2max = getArrayValueStr("max", nx, ny, nz, S2);
          printf("S1 (min,max) = (%8.6f, %8.6f)\n", x1min, x1max);
          printf("S2 (min,max) = (%8.6f, %8.6f)\n", x2min, x2max);
        }
      else if (opts.inputFile)
        {
          /* read in the moment file and return the number of moments read
             NOTE: range is also set correctly in the file */
          nIter = readMomentFile(&opts, &iter, &time, &meanF, &varF);
        }
      else
        {
          printf("Error: please specify mean/variance or betas combo\n");
          exit(2);
        }

      // store the grid size for use locally
      n[0] = n1 = opts.n1;
      n[1] = n2 = opts.n2;

      // allocate local memory
      x        = (double**) calloc(2, sizeof(double*));
      dx       = (double**) calloc(2, sizeof(double*));
      pdf      = (double**) calloc(n2, sizeof(double*));
      pdfBound = (int**)    calloc(n2, sizeof(int*));
      pdf_tmp  = (double*)  calloc(n1*n2, sizeof(double));
      pdf_dns  = (double**) calloc(n2, sizeof(double*));
      for (k=0; k<n2; k++)
        {
          pdf[k]      = (double*) calloc(n1, sizeof(double));
          pdf_dns[k]  = (double*) calloc(n1, sizeof(double));
          pdfBound[k] = (int*) calloc(2, sizeof(int));
        }

      // initialise each grid, for now using same parameters in each direction
      initGrid(n1,x1min,x1max,opts.nPiv1,opts.gPiv1,opts.gPerc1,&x[0],&dx[0]);
      initGrid(n2,x2min,x2max,opts.nPiv2,opts.gPiv2,opts.gPerc2,&x[1],&dx[1]);

      // local grid pointers
      x1 = x[0];
      x2 = x[1];
  
      // for computing the mean, compute the function Y
      Y  = (double**) calloc(3, sizeof(double*));
      for (i=0; i<3; i++) Y[i] = (double*) calloc(n2*n1, sizeof(double));
      for (k=0; k<n2; k++)
        {
          row = k*n1;
          for (j=0; j<n1; j++)
            {
              if (opts.transform == 0)
                {
                  Y[0][row+j] = x1[j] + x2[k];
                  Y[1][row+j] = x1[j];
                  Y[2][row+j] = x2[k];
                }
              else if (opts.transform == 1)
                {
                  Y[0][row+j] = x1[j] + x2[k]*(1.0-x1[j]);
                  Y[1][row+j] = x1[j];
                  Y[2][row+j] = x2[k]*(1.0-x1[j]);
                }
              else if (opts.transform == 2)
                {
                  Y[0][row+j] = x1[j];
                  Y[1][row+j] = x1[j]*(1 - x2[k]);
                  Y[2][row+j] = x2[k]*x1[j];
                }
            } // end loop over columns
        } // end loop over rows

      // clip the means/variances
      if (!opts.extract)
        clipMoments2D(&mean[1], &var[1], &mean[2], &var[2], &var[0], &var[3]);

      // write header for error file
      if (opts.inputFile)
        {
          errFile = fopen(opts.errorFile, "w");
          fprintf(errFile, "# iter time ");
          for (i=0; i<opts.nMethod; i++)
            fprintf(errFile, "%s-m v m1 v1 m2 v2 v3 KL JS HD TVD ",
                    opts.method[i]);
          fprintf(errFile, "\n");
        }

      for (f=0; f<nIter; f++)
        { 
          // write output data file name
          strcpy(dataFile, "mduc_jpdf_");
          strcat(dataFile, opts.ident);
          if (opts.inputFile)
            sprintf(buffer, "%05i", iter[f]);
          else
            sprintf(buffer, "0");
          strcat(dataFile, buffer);
          strcat(dataFile, ".stat");

          if (opts.interpFile)
            {
              // write file to interpolate from
              strcpy(interpFile, opts.prefix);
              strcat(interpFile, buffer);
              strcat(interpFile, ".flo");
              
              // check if the file exists
              pFile = fopen(interpFile, "r");
              if (!pFile)
                {
                  fprintf(stderr, "no data\n");
                  continue;
                }
              else
                {
                  fprintf(stderr, "DNS file: %s --> ", interpFile);

                  // get time from file
                  tFile = getFileTime(pFile);
                  if (tFile < 0.0)
                    {
                      fprintf(stderr, "Error: time info not present in file\n");
                      exit(2);
                    }
                  fclose(pFile);

                  // get the dns pdf to compute KL distance
                  pdfFlag = interpPDFfromFile(interpFile, n, x, pdf_dns);
                }
            }

          if (opts.inputFile)
            {
              // make sure that the time in file matches moment file
              iMom = -1;
              for (jj=0; jj<=nIter; jj++)
                {
                  if (fabs(tFile - time[jj]) < 1e-6)
                    iMom = jj;
                }

              if (iMom == -1)
                {
                  fprintf(stderr, "Error: time not found in moment file\n");
                  exit(2);
                }

              // assign moments from file according to appropriate time/row
              for (i=0; i<3; i++)
                {
                  mean[i] = meanF[iMom][i];
                  var[i]  = varF[iMom][i];
                }
              var[3] = varF[iMom][3];
            }

          // write header for file
          writeStatsHeader(dataFile, 2, n, x, tFile, mean, var, opts.nMethod,
                           opts.transformName);

          if (errFile) fprintf(errFile, "%6i%15.6e", iter[f], time[f]);

          for (m=0; m<opts.nMethod; m++)
            {
              // compute the beta parameters from the input values
              if (!strncmp(opts.method[m], "diric", 5))
                {
                  S=mean[1]*mean[1]+mean[2]*mean[2]+(1.0-mean[0])*(1.0-mean[0]);
                  Q=var[0] + var[1] + var[2];
                  tmp = (1.0 - S)/Q - 1.0;
                  beta[1] = mean[1]*tmp;
                  beta[2] = mean[2]*tmp;
                  beta[3] = (1.0-mean[0])*tmp;
                  beta[0] = beta[1] + beta[2] + beta[3];

                  // compute the analytic variances for the dirichlet
                  varAn[1] = beta[1]*(beta[0] - beta[1])/tmp;
                  varAn[2] = beta[2]*(beta[0] - beta[2])/tmp;
                  varAn[0] = beta[3]*(beta[0] - beta[3])/tmp;
                  varAn[3] = 0.5*(varAn[0] - varAn[1] - varAn[2]);
                }

              // compute the PDF
              if (opts.extract)
                pdfFlag = binPDF_2D(nx, ny, nz, S1, S2, n, x, pdf);
              else if (!strcmp(opts.method[m], "interp"))
                pdfFlag = interpPDFfromFile(interpFile, n, x, pdf);
              else
                pdfFlag = MDUCcalcPDF_2D(opts.method[m], n, x, mean[1], var[1],
                                         mean[2], var[2], var[0], pdf, pdfBound,
                                         opts.transform);

              // compute the total mean and the mean of each direction
              for (i=0; i<3; i++)
                meanPDF[i] = calcLinMeanY_2D(0, 1, n, x, Y[i], beta[1], beta[2],
                                             beta[3], pdf, pdfBound, pdfFlag,
                                             opts.transform);

              // compute all variances
              for (i=0; i<3; i++)
                varPDF[i] = calcVar(n,0,1,meanPDF[i],meanPDF[i],Y[i],Y[i],pdf);

              // compute the covariance
              varPDF[3] = calcVar(n,0,1,meanPDF[1],meanPDF[2],Y[1],Y[2],pdf);

              // unscale pdf if desired (scaled by default)
              if (opts.scaled)
                unscalePDF(1,opts.transform,n1,n2,x1,x2,pdf);

              // write the pdf to the file
              convertArrayMD_1D(n1, n2, 1, &pdf, pdf_tmp);
              if (!strcmp(opts.method[m], "interp"))
                sprintf(buffer, "pdf-dns");
              else
                sprintf(buffer, "pdf-%s", opts.method[m]);
              addVarToBody(dataFile, n1, n2, pdf_tmp, buffer);

              // compute the Kullback-Leibler distance to any DNS data
              fDiv[0] = fDiv[1] = fDiv[2] = 0.0;
              if (opts.interpFile)
                {
                  fDiv[0] = calcKullbackLeibler(n1,n2,pdf_dns,pdf);
                  fDiv[1] = calcJensonShannon(n1,n2,pdf_dns,pdf);
                  fDiv[2] = calcHellinger(n1,n1,pdf_dns,pdf);
                  fDiv[3] = calcTVD(n1,n1,pdf_dns,pdf);
                }
      
              // append the error information to the file      
              pFile = fopen(dataFile, "a");
              printStatsError(pFile,errFile,opts.method[m],2,mean,var,meanPDF,
                              varPDF, fDiv);
              fclose(pFile);

              // print results to screen
              if (!opts.inputFile)
                printStatsError(stderr,NULL,opts.method[m],2,mean,var,meanPDF,
                                varPDF, fDiv);

            } // end loop over methods
          if (errFile)
            {
              fprintf(errFile, " %15.6e %15.6e", var[1], var[2]);
              fprintf(errFile, "\n");
              fflush(errFile);
            }
        } // end loop over cases

      if (errFile) fclose(errFile);

      return 0;
    }
}
/* ------------------------------------------------------------------------- */

void readList(int n, int iStart, char ***var, char **argv)
{
  int i, j;
  
  var[0] = (char**) calloc(n, sizeof(char*));
  for (i=0, j=iStart; i<n; i++, j++)
    {
      var[0][i] = (char*) calloc(strlen(argv[j])+10, sizeof(char));
      strcpy(var[0][i], argv[j]);
    }
}
/* -------------------------------------------------------------------------- */

void readArray(int n, int iStart, double **var, char **argv)
{
  int i, j;
  var[0] = (double*) calloc(n, sizeof(double));
  for (i=0, j=iStart; i<n; i++, j++)
    var[0][i] = atof(argv[j]);
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
        case 'd':
          opts->dim = atoi(argv[iOpt+1]);
          break;
        case 's':
          opts->scaled = 1;
          break;
        case 'M':
          opts->nMethod = iEnd - iOpt;
          readList(opts->nMethod, iOpt+1, &opts->method, argv);
          break;
        case 'i':
          opts->readOptFile = 1;
          strcpy(opts->optFile, argv[iOpt+1]);
          break;
        case 'I':
          strcpy(opts->ident, argv[iOpt+1]);
          break;
        case 'T':
          opts->runTest = 1;
          break;
        case 'e':
          opts->extract = 1;
          opts->runTest = 0;
          opts->dim = iEnd - iOpt;
          if (opts->dim > 2)
             {
                printf("Error: pdf of more than 2 variables not supported\n");
                exit(2);
              }
          else 
            readList(opts->dim, iOpt+1, &opts->names, argv);
          break;
        case 'E':
          opts->errorTest = 1;
          opts->runTest = 0;
          break;
        case 'R':
          opts->regime = 1;
          opts->runTest = 0;
          break;
        case 'L':
          opts->meanLin = 0;
          break;
        case 'm':
          opts->dim = iEnd - iOpt;
          opts->mean1 = atof(argv[iOpt+1]);
          if (opts->dim == 2)
            opts->mean2 = atof(argv[iOpt+2]);
          opts->runTest   = 0;
          opts->inputMean = 1;
          break;
        case 'n':
          opts->n1 = atoi(argv[iOpt+1]);
          if (argv[iOpt+2] != NULL)
            {
              if (atoi(argv[iOpt+2]))
                opts->n2 = atoi(argv[iOpt+2]);
            }
          break;
        case 'v':
          opts->var1 = atof(argv[iOpt+1]);
          if (argv[iOpt+2] != NULL)
            {
              opts->var2 = atof(argv[iOpt+2]);
              if (argv[iOpt+3] != NULL)
                opts->var12 = atof(argv[iOpt+3]);
            }
          opts->runTest = 0;
          break;
        case 'g':
          opts->nPiv1 = iEnd - iOpt;
          readArray(opts->nPiv1, iOpt+1, &opts->gPiv2, argv);
          break;
        case 'b':
          opts->beta1 = atof(argv[iOpt+1]);
          if (argv[iOpt+2] != NULL)
            {
              if (atof(argv[iOpt+2]))
                opts->beta2 = atof(argv[iOpt+2]);
              if (argv[iOpt+3] != NULL)
                if (atof(argv[iOpt+3]))
                  opts->beta3 = atof(argv[iOpt+3]);
            }
          opts->runTest   = 0;
          opts->inputBeta = 1;
          break; 
        case 't':
          strcpy(opts->transformName, argv[iOpt+1]);
          if (!strcmp(opts->transformName, "none"))
            opts->transform = 0;
          else if (!strcmp(opts->transformName, "doran"))
            opts->transform = 1;
          else if  (!strcmp(opts->transformName, "hasse"))
            opts->transform = 2;
          break;
        case '-': // has a longer word
          {
            if (!strcmp(argv[optLocs[i]], "--prefix"))
              {
                opts->interpFile = 1;
                strcpy(opts->prefix, argv[iOpt+1]);
              }
            if (!strcmp(argv[optLocs[i]], "--range"))
              {
                opts->range[0] = atoi(argv[iOpt+1]);
                opts->range[1] = atoi(argv[iOpt+2]);
              }
            if (!strcmp(argv[optLocs[i]], "--set"))
              {
                opts->dataSet = atoi(argv[iOpt+1]);
              }
            if (!strcmp(argv[optLocs[i]], "--mfile"))
              {
                opts->inputFile = 1;
                strcpy(opts->momentFile, argv[iOpt+1]);
              }
            break;
          }
        } // end switch
    } // end loop over options
  
  free(optLocs);
}
/* -------------------------------------------------------------------------- */
