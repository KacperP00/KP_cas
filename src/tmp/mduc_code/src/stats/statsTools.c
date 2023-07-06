#include "statsTools.h"
#include "accessMDUC.h"
#include "fileIO.h"
#include "parser.h"
#include "init.h"
#include "statsUtils.h"
#include "solveUtils.h"
#include "betaPDF.h"
#include "smlPDF.h"
#include "utils.h"

void readOptionsFile(optStruct *opts)
{
  int   i, j;
  double *moments = NULL;

  // strings to store the tag for each line and comments (if any)
  char  tag[50], comment[200];
  
  FILE   *pFile;   // input stream
  pFile = fopen(opts->optFile, "r"); // open input file

  if (!pFile)
    {
      printf("ERROR: cannot open input file %s\n", opts->optFile);
      exit(1);
    }

  while (!feof(pFile))
    { 
      fscanf(pFile, "%s", tag); // read in the first string
      
      // first check to see if end of input is desired
      if (!strcmp(tag, "ENDINPUT"))
        break;
      // check to see if it is a comment
      else if (tag[0] == '!' || tag[0] == '#')
        fgets(comment, 200, pFile);

      // grid parameters
      else if (!strcmp(tag, "dimensions"))
        opts->dim = readInteger(pFile);
      else if (!strcmp(tag, "n1"))
        opts->n1 = readInteger(pFile);
      else if (!strcmp(tag, "n2"))
        opts->n2 = readInteger(pFile);
      else if (!strcmp(tag, "grid-pivot-1"))
        opts->nPiv1 = readRealArraySize(pFile, &opts->gPiv1);
      else if (!strcmp(tag, "grid-percent-1"))
        readRealArraySize(pFile, &opts->gPerc1);
      else if (!strcmp(tag, "grid-pivot-2"))
        opts->nPiv2 = readRealArraySize(pFile, &opts->gPiv2);
      else if (!strcmp(tag, "grid-percent-2"))
        readRealArraySize(pFile, &opts->gPerc2);

      // moments
      else if (!strcmp(tag, "moments"))
        {
          opts->inputMean = 1;
          opts->runTest = 0;
          readRealArray(pFile, 5, &moments);
          opts->mean1 = moments[0];
          opts->mean2 = moments[1];
          opts->var1  = moments[2];
          opts->var2  = moments[3];
          opts->var12 = moments[4];
        }
      else if (!strcmp(tag, "mean-1"))
        {
          opts->inputMean = 1;
          opts->runTest = 0;
          opts->mean1 = readReal(pFile);
        }
      else if (!strcmp(tag, "mean-2"))
        opts->mean2 = readReal(pFile);
      else if (!strcmp(tag, "variance-1"))
        opts->var1 = readReal(pFile);
      else if (!strcmp(tag, "variance-2"))
        opts->var2 = readReal(pFile);
      else if (!strcmp(tag, "variance-12"))
        opts->var12 = readReal(pFile);
      else if (!strcmp(tag, "max-1"))
        opts->max1 = readReal(pFile);
      else if (!strcmp(tag, "max-2"))
        opts->max2 = readReal(pFile);
      else if (!strcmp(tag, "data-set"))
        opts->dataSet = readReal(pFile);

      // file names
      else if (!strcmp(tag, "moment-file"))
        {
          opts->inputFile = 1;
          readString(pFile, opts->momentFile);
        }
      else if (!strcmp(tag, "file-prefix"))
        {
          opts->interpFile = 1;
          readString(pFile, opts->prefix);
        }
      else if (!strcmp(tag, "error-file"))
        readString(pFile, opts->errorFile);
      else if (!strcmp(tag, "file-ident"))
        readString(pFile, opts->ident);
      else if (!strcmp(tag, "start-step"))
        opts->range[0] = readInteger(pFile);
      else if (!strcmp(tag, "final-step"))
        opts->range[1] = readInteger(pFile);
      
      // methods and transformation
      else if (!strcmp(tag, "methods"))
        opts->nMethod = readNames(pFile, &opts->method);
      else if (!strcmp(tag, "transformation"))
        {
          readString(pFile, opts->transformName);
          if (!strcmp(opts->transformName, "none"))
            opts->transform = 0;
          else if (!strcmp(opts->transformName, "doran"))
            opts->transform = 1;
          else if  (!strcmp(opts->transformName, "hasse"))
            opts->transform = 2;
        }
      else if (!strcmp(tag, "scaling"))
        opts->scaled = readInteger(pFile);
    }

  fclose(pFile);
}
/* ------------------------------------------------------------------------- */

int readMomentFile(optStruct *opts, int **iter, double **time, double ***mean,
                   double ***var)
{
  int i, j, k, len, nt, step, nset, nhead;

  char tag[50], file[100], statFile[100], comment[500], *head;

  // read in moments for time steps specified
  FILE *pFile;
  pFile = fopen(opts->momentFile, "r");

  // get the number of data lines  
  len = 0;
  while (!feof(pFile))
    {
      strcpy(comment, "\n");
      fgets(comment, 1000, pFile);
      if (strcmp(comment, "\n"))
        len++;
    }
  len -= 1;
  rewind(pFile);

  // read in the header to see how many data sets there are
  fgets(comment, 1000, pFile);
  head = strtok(comment, " ");
  nhead = 0;
  while (head != NULL) {
    head = strtok(NULL, " ");
    nhead++;
  }
  nset = (nhead - 2)/5;
  fprintf(stderr, "%i sets found\n", nset);

  // set the range to whole file if not explicitly set
  if (opts->range[0] < 0)
    {
      opts->range[0] = 0;
      opts->range[1] = len-1;
    }

  nt = opts->range[1] - opts->range[0] + 1;

  // allocate the memory
  iter[0] = (int*)     calloc(nt, sizeof(int));
  time[0] = (double*)  calloc(nt, sizeof(double));
  mean[0] = (double**) calloc(nt, sizeof(double*));
  var[0]  = (double**) calloc(nt, sizeof(double*));
  for (i=0; i<nt; i++)
    {
      mean[0][i] = (double*) calloc(3, sizeof(double));
      var[0][i]  = (double*) calloc(4, sizeof(double));
    }

  // read in the moments
  fprintf(stderr, "%s\n", comment);
  for (i=0, j=0; i<len; i++)
    {
      step = readInteger(pFile);
      step = i;
      if (step >= opts->range[0] && step <= opts->range[1])
        {
          iter[0][j]    = step;
          time[0][j]    = readReal(pFile);

          for (k=0; k < 5*(opts->dataSet-1); k++) readReal(pFile);

          mean[0][j][1] = readReal(pFile);
          var[0][j][1]  = readReal(pFile);
          mean[0][j][2] = readReal(pFile);
          var[0][j][2]  = readReal(pFile);
          var[0][j][3]  = readReal(pFile);
          mean[0][j][0] = mean[0][j][1] + mean[0][j][2];
          var[0][j][0]  = 2*var[0][j][3] + var[0][j][1] + var[0][j][2];

          if (opts->dataSet < nset) fgets(comment, 1000, pFile);
          fprintf(stderr, "%12i %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                  step, time[0][j], mean[0][j][1], var[0][j][1], mean[0][j][2],
                  var[0][j][2], var[0][j][3]);
          j++;
        }
      else
        fgets(comment, 1000, pFile); // read in rest of line
    }

  fclose(pFile);

  return nt;
}
/* ------------------------------------------------------------------------- */

double getFileTime(FILE *pFile)
{
  // strings to store the tag for each line and comments (if any)
  char tag[50], comment[500];

  while (!feof(pFile))
    { 
      /* read in string */
      fscanf(pFile, "%s", tag);

      if (!strcmp(tag, "time"))
        return readReal(pFile);
      else
        fgets(comment, 200, pFile);
    }
  return -1;
}
/* ------------------------------------------------------------------------- */

int interpPDFfromFile(char *fileName, int *n, double **x, double **pdf)
{
  int     i, j, k, row, col;
  int     n1, n2, n1F, n2F;
  double *x1, *x2, *x1F, *x2F, dx1, dx2, sumPDF;
  double *tmp, *pdfF, *pdf_tmp;
  x1F = x2F = pdfF = NULL;

  // strings to store the tag for each line and comments (if any)
  char tag[50], comment[500];

  FILE *pFile = fopen(fileName, "r");

  if (pFile)
    {
      n1 = n[0];
      n2 = n[1];
      x1 = x[0];
      x2 = x[1];

      pdf_tmp = (double*) calloc(n1*n2, sizeof(double));
     
      while (!feof(pFile))
        { 
          /* read in string */
          fscanf(pFile, "%s", tag);

          if (tag[0] == '!')
            fgets(comment, 200, pFile);
          /* ----- read in grid data ----- */
          else if (!strcmp(tag, "gridPoints") || !strcmp(tag, "gridPoints1"))
            {
              n1F = readInteger(pFile);
              n2F = 1;
            }
          else if (!strcmp(tag, "gridPoints2"))
            {
              n2F = readInteger(pFile);
            }
          else if (!strcmp(tag, "Z1") || !strcmp(tag, "H1"))
            readRealArray(pFile, n1F, &x1F);
          else if (!strcmp(tag, "Z2") || !strcmp(tag, "H2"))
            readRealArray(pFile, n2F, &x2F);
          /* ----- read in grid data ----- */
          else if (!strcmp(tag, "pdf"))
            readRealArray(pFile, n1F*n2F, &pdfF);
        }
      fclose(pFile);

      tmp = (double*) calloc(n1F*n2F, sizeof(double));

      if (n1F == n1 && n2F == n2)
        {
          // convert to structured array
          convertArray1D_MD(n1, n2, 1, pdfF, &pdf);
        }
      else
        {
          fprintf(stderr, "Interpolating DNS data to grid...\n");
          sumPDF = 0.0;
          // compute the pdf without the grid
          for (k=0; k<n2F; k++)
            {
              if (k == 0)
                dx2 = 0.5*(x2F[1]-x2F[0]);
              else if (k == n2F-1)
                dx2 = 0.5*(x2F[k]-x2F[k-1]);
              else
                dx2 = 0.5*(x2F[k+1]-x2F[k-1]);
            
              for (j=0; j<n1F; j++)
                {
                  if (j == 0)
                    dx1 = 0.5*(x1F[1]-x1F[0]);
                  else if (j == n1F-1)
                    dx1 = 0.5*(x1F[j]-x1F[j-1]);
                  else
                    dx1 = 0.5*(x1F[j+1]-x1F[j-1]);

                  tmp[k*n1F+j] = pdfF[k*n1F + j]/(dx1*dx2);

                  sumPDF += pdfF[k*n1F + j];
                }
            }

          // interpolate solution onto new grid between 0 and 1
          /* NOTE: should consider here to simply add pdf contribution that is
             out of bounds to the boundary points */
          interpSolution2D(n1F, n2F, x1F, x2F, n1, n2, x1, x2, 1, tmp, pdf_tmp);

          // convert to structured array
          convertArray1D_MD(n1, n2, 1, pdf_tmp, &pdf);

          // normalise the new pdf
          sumPDF = 0.0;
          for (k=0; k<n2; k++)
            {
              if (k == 0)
                dx2 = 0.5*(x2[1]-x2[0]);
              else if (k == n2-1)
                dx2 = 0.5*(x2[k]-x2[k-1]);
              else
                dx2 = 0.5*(x2[k+1]-x2[k-1]);
            
              for (j=0; j<n1; j++)
                {
                  if (j == 0)
                    dx1 = 0.5*(x1[1]-x1[0]);
                  else if (j == n1-1)
                    dx1 = 0.5*(x1[j]-x1[j-1]);
                  else
                    dx1 = 0.5*(x1[j+1]-x1[j-1]);

                  pdf[k][j] = pdf[k][j]*dx1*dx2;

                  sumPDF += pdf[k][j];
                }
            }
        }

      free(x1F);
      free(x2F);
      free(tmp);
      free(pdf_tmp);
    }
  else
    {
      fprintf(stderr, "Warning: couldn't open file %s\n", fileName);
    }

  return -1;
}
/* ------------------------------------------------------------------------- */

void writeStatsHeader(char *fileName, int dim, int *n, double **x, double t,
                      double *mean, double *var, int nPDF, char *transform)
{
  time_t rawtime;
  struct tm *timeinfo;
  char date[100];

  // grid info
  int n1, n2;
  double *x1, *x2; 
  n1 = n[0];
  x1 = x[0];
  if (dim == 2)
    {
      n2 = n[1];
      x2 = x[1];
    }
  double covar, rho;
  
  // get time for run
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(date, 100, "%c", timeinfo);

  FILE *pFile = fopen(fileName, "w");
      
  /* ----- write header ----- */
  fprintf(pFile, "header\n");
  fprintf(pFile, "title = \"Statistics Data\"\n");
  fprintf(pFile, "date  = \"%s\"\n\n", date);
  
  fprintf(pFile, "time   = %13.6e [s]\n", t);
  fprintf(pFile, "Z1_min = %13.6e\n", x1[0]);
  fprintf(pFile, "Z1_max = %13.6e\n", x1[n1-1]);
  if (dim == 2)
    {
      fprintf(pFile, "Z2_min = %13.6e\n", x2[0]);
      fprintf(pFile, "Z2_max = %13.6e\n", x2[n2-1]); 
    }
  fprintf(pFile, "\n");

  // compute covariance and correlation
  covar = 0.5*(var[2]-var[1]-var[0]);
  if (var[1] > 0.0 && var[2] > 0.0)
    rho = var[3]/pow(var[1]*var[2],0.5);
  else
    rho = 0.0;

  if (dim == 1)
    {
      fprintf(pFile, "Z1_mean     = %13.6e\n", mean[dim-1]);
      fprintf(pFile, "Z1_variance = %13.6e\n", var[dim-1]);
    }
  else if (dim == 2)
    {
      fprintf(pFile, "Z1_mean           = %13.6e\n", mean[dim-1]);
      fprintf(pFile, "Z1_variance       = %13.6e\n", var[dim-1]);
      fprintf(pFile, "Z2_mean           = %13.6e\n", mean[2]);
      fprintf(pFile, "Z2_variance       = %13.6e\n", var[2]);
      fprintf(pFile, "Z1_Z2_covariance  = %13.6e\n", covar);
      fprintf(pFile, "correlation coeff = %7.4f\n", rho);
    }
  fprintf(pFile, "\n");

  fprintf(pFile, "distributions  = %i\n", nPDF);
  fprintf(pFile, "gridPoints     = %i\n", n1);
  if (dim == 2)
    {
      fprintf(pFile, "gridPoints2    = %i\n", n2);
      fprintf(pFile, "transformation = %s\n", transform);
    }
  fprintf(pFile, "\n");

  /* ----- write body ----- */
  fprintf(pFile, "body\n");
  // print the mesh
  if (!strcmp(transform, "hasse"))
    fprintf(pFile, "Z1+Z2\n");
  else
    fprintf(pFile, "Z1\n");
  writeArray(pFile, 1, n1, 1, 0, 5, x1);
  if (dim == 2)
    {
      if (!strcmp(transform, "doran"))
        fprintf(pFile, "Z2/(1-Z1)\n");
      else if (!strcmp(transform, "hasse"))
        fprintf(pFile, "Z2/(Z1+Z2)\n");
      else
        fprintf(pFile, "Z2\n");
      writeArray(pFile, 1, n2, 1, 0, 5, x2);
    }

  fprintf(pFile, "trailer\n");

  fclose(pFile);
}
/* ------------------------------------------------------------------------- */

void printStatsError(FILE *pFile, FILE *err, char *method, int dim, double *mean,
                     double *var, double *meanPDF, double *varPDF, double *fDiv)
{
  int i, nMom = 1;
  char *scalar[] = {"Z  -  ", "Z1 -  ", "Z2 -  "};

  if (dim == 1)
    nMom = 1;
  else if (dim == 2)
    nMom = 3;

  fprintf(pFile, "\nResults for %s pdf:\n", method);
  fprintf(pFile, "----------------- input ----- computed ----- Error (%%) -----");
  fprintf(pFile, "\n");
  for (i=0; i<nMom; i++)
    {
      fprintf(pFile, "%s%8s: %9.6f %12.6f %14.4f\n", scalar[i], "mean",
              mean[i], meanPDF[i], percError(meanPDF[i], mean[i]));
      fprintf(pFile, "%s%8s: %9.6f %12.6f %14.4f\n", "      ", "variance",
              var[i], varPDF[i], percError(varPDF[i], var[i]));

      if (err)
        fprintf(err, "%15.6e%15.6e", percError(meanPDF[i], mean[i]),
                percError(varPDF[i], var[i]));
    }

  if (dim == 2)
    {
      double rho, rhoPDF;
      fprintf(pFile, "%s%8s: %9.6f %12.6f %14.4f\n", "Z12 co", "variance",
              var[3], varPDF[3], percError(varPDF[3], var[3]));

      if (err) fprintf(err, "%15.6e", percError(varPDF[3], var[3]));

      // compute the correlation coefficient based on the computed co-variance
      if (var[1] > 0.0 && var[2] > 0.0)
        rho = var[3]/pow(var[1]*var[2],0.5);
      else
        rho = 0.0;

      if (varPDF[1] > 0.0 && varPDF[2] > 0.0)
        rhoPDF = varPDF[3]/pow(varPDF[1]*varPDF[2],0.5);
      else
        rhoPDF = 0.0;

      fprintf(pFile, "\n   correlation: %9.6f %12.6f %14.4f\n", rho, rhoPDF,
              percError(rhoPDF, rho));
    }

  // print Kullback-Leibler distance
  if (fDiv != NULL)
    {
      fprintf(pFile, "\n f-divergence values:\n");
      fprintf(pFile, "    Kullback-Leibler: %9.6f\n", fDiv[0]);
      fprintf(pFile, "    Jenson-Shannon:   %9.6f\n", fDiv[1]);
      fprintf(pFile, "    Hellinger:        %9.6f\n", fDiv[2]);
      fprintf(pFile, "    TVD:              %9.6f\n", fDiv[3]);
      if (err)
        fprintf(err, "%15.6e %15.6e %15.6e %15.6e",
                fDiv[0],fDiv[1],fDiv[2],fDiv[3]);
    }

  if (!strcmp(method, "beta"))
    {
      double alpha, beta, alphaPDF, betaPDF, tmp;
      tmp      = mean[0]*(1.0-mean[0])/var[0] - 1.0;
      alpha    = mean[0]*tmp;
      beta     = (1.0-mean[0])*tmp;
      tmp      = meanPDF[0]*(1.0-meanPDF[0])/varPDF[0] - 1.0;
      alphaPDF = meanPDF[0]*tmp;
      betaPDF  = (1.0-meanPDF[0])*tmp;
      fprintf(pFile, "\n%14s: %9.6f %12.6f %14.4f\n", "alpha", alpha, alphaPDF,
              percError(alphaPDF, alpha));
      fprintf(pFile, "%14s: %9.6f %12.6f %14.4f\n", "beta", beta, betaPDF,
              percError(betaPDF, beta));
    }

  if (!strncmp(method, "diric", 5))
    {
      double S, Q, beta[4], betaPDF[4], tmp;
      // compute mean and variance parameters
      S = mean[1]*mean[1] + mean[2]*mean[2] + (1-mean[0])*(1-mean[0]);
      Q = var[1] + var[2] + var[0];

      // compute the shape parameters for the beta distribution
      tmp   = (1.0 - S)/Q - 1.0;
      beta[1] = mean[1]*tmp;
      beta[2] = mean[2]*tmp;
      beta[3] = (1.0-mean[0])*tmp;
      beta[0] = beta[1] + beta[2] + beta[3];

      S = meanPDF[1]*meanPDF[1] + meanPDF[2]*meanPDF[2]
        + (1-meanPDF[0])*(1-meanPDF[0]);
      Q = varPDF[1] + varPDF[2] + varPDF[0];

      tmp   = (1.0 - S)/Q - 1.0;
      betaPDF[1] = meanPDF[1]*tmp;
      betaPDF[2] = meanPDF[2]*tmp;
      betaPDF[3] = (1.0-meanPDF[0])*tmp;
      betaPDF[0] = betaPDF[1] + betaPDF[2] + betaPDF[3];
      
      fprintf(pFile, "\n");
      for (i=0; i<4; i++)
        fprintf(pFile, "%13s%i: %9.6f %12.6f %14.4f\n", "beta", i, beta[i],
                betaPDF[i], percError(betaPDF[i], beta[i]));

    }

  fprintf(pFile, "------------------------------------------------------------");
  fprintf(pFile, "\n");
}
/* ------------------------------------------------------------------------- */

void runTest1D(optStruct *opts)
{
  int    i, j, k;
  int    n, nTest;
  int    pdfBound[11][2], pdfFlag;
  char   dataFile[200], tmp[50];
  double mean, var, meanPDF;
  double *x, *dx, *Y;
  double alpha[11], beta[11];
  double **pdf;
  double eps0, eps1;
  double error[11];

  nTest = 11;
  n = opts->n1;

  // allocate memory
  pdf  = (double**) calloc(nTest, sizeof(double*));
  for (i=0; i<nTest; i++)
    pdf[i] = (double*) calloc(n, sizeof(double));

  initGrid(n, 0.0, 1.0, opts->nPiv1, opts->gPiv1, opts->gPerc1, &x, &dx);

  // generate the solution (ie Y = x), to compute the mean
  Y = (double*) calloc(n, sizeof(double));
  for (j=0; j<n; j++)
    Y[j] = x[j];

  // create file name
  strcpy(dataFile, "PDFverify1D_");
  sprintf(tmp, "%ipt.dat", n);
  strcat(dataFile, tmp);

  FILE *pFile = fopen(dataFile, "w");

  // calculate a range of alpha/beta to cover different shapes
  // uniform distribution
  alpha[0]  = 1;
  beta[0]   = 1;
  // u-shaped distribution
  alpha[1]  = 0.5;
  beta[1]   = 0.5;
  // decreasing distributions
  alpha[2]  = 1;
  beta[2]   = 2.5; // convex
  alpha[3]  = 1;
  beta[3]   = 2;   // straight line
  alpha[4]  = 1;
  beta[4]   = 1.5; // concave
  // increasing distributions
  alpha[5]  = 2.5;
  beta[5]   = 1;   // convex
  alpha[6]  = 2;
  beta[6]   = 1;   // straight line
  alpha[7]  = 1.5;
  beta[7]   = 1;   // concave
  // unimodal distributions
  alpha[8]  = 2;
  beta[8]   = 2;
  alpha[9]  = 2;
  beta[9]   = 5;
  alpha[10] = 5;
  beta[10]  = 2;

  for (i=0; i<nTest; i++)
    {
      mean = alpha[i]/(alpha[i]+beta[i]);
      var  = alpha[i]*beta[i]/(pow(alpha[i]+beta[i],2)*(alpha[i]+beta[i]+1));

      pdfFlag = MDUCcalcPDF_1D(opts->method[0],n,x, mean, var, pdf[i], pdfBound[i]);
 
      // compute the error
      if (opts->meanLin)
        meanPDF = calcLinMeanY(0, 1, alpha[i], beta[i], n, x, Y,
                               pdf[i], pdfBound[i], pdfFlag);
      else
        {
          Y[0] = Y[0] + 0.5*(x[1]-x[0]);
          Y[n-1] = Y[n-1] - 0.5*(x[n-1] - x[n-2]);
          meanPDF = calcMean(pdfBound[i][0], pdfBound[i][1], 0,1, pdf[i], Y);
        }
          
      error[i] = fabs(mean - meanPDF)/mean;

    }

  // write the data set to a file
  for (j=0; j<n; ++j)
    {
      fprintf(pFile, "%6.4f ", x[j]);
      for (i=0; i<nTest; ++i)
        fprintf(pFile, "%15.6e ", pdf[i][j]);
      fprintf(pFile, "\n");
    }

  fclose(pFile);

  // print the results to screen
  printf("\nError of standard 1D cases:\n");
  printf("------------------------------------------------------------------\n");
  printf("%6s %6s %13s\n","alpha","beta","Error");
  for (i=0; i<nTest; ++i)
    printf("%6.2f %6.2f %10.4f [%2i,%3i]\n", alpha[i], beta[i],
           error[i]*100, pdfBound[i][0], pdfBound[i][1]);
  printf("------------------------------------------------------------------\n");

  for (i=0; i<nTest; i++) free(pdf[i]);
  free(pdf);
  free(Y);
  free(x);
  free(dx);
}
/* ------------------------------------------------------------------------- */

void errorTest1D(optStruct *opts)
{
  int i, j, k;
  int nMean, nVar, nPtsG1 = 0;
  double dMean, dVar;
  double varMax, varCur;
  double alpha, beta, tmp, meanPDF;
  double *Zmean, *Zvar; // arrays for all possible means/variances

  double *x, *dx, *Y;
  int    pdfBound[2], pdfFlag;
  double *pdf;
  double **error, **absDiff, maxErrVar, maxErrMean, maxErr = 0.0;

  nMean = 201;
  nVar  = 201;

  dMean = 1.0/(nMean-1);
  dVar = 1.0/(nVar-1);

  Zmean   = (double*)  calloc(nMean, sizeof(double));
  Zvar    = (double*)  calloc(nVar, sizeof(double));
  pdf     = (double*)  calloc(opts->n1, sizeof(double));
  error   = (double**) calloc(nMean, sizeof(double*));
  absDiff = (double**) calloc(nMean, sizeof(double*));
  for (k=0; k<nMean; k++)
    {
      error[k]   = (double*) calloc(nVar, sizeof(double));
      absDiff[k] = (double*) calloc(nVar, sizeof(double));
    }

  initGrid(opts->n1, 0.0, 1.0, opts->nPiv1, opts->gPiv1, opts->gPerc1, &x, &dx);

  // generate the solution (ie Y = x), to compute the mean
  Y = (double*) calloc(opts->n1, sizeof(double));
  for (j=0; j<opts->n1; j++)
    Y[j] = x[j];

  Zmean[0] = 0.0;
  for (k=1; k<nMean; k++) Zmean[k] = Zmean[k-1] + dMean;

  Zvar[0] = 0.0;
  for (j=1; j<nVar; j++) Zvar[j] = Zvar[j-1] + dVar;
  
  for (k=0; k<nMean; k++)
    {
      varMax = Zmean[k]*(1.0-Zmean[k]);
      for (j=0; j<nVar; j++)
        {
          varCur = Zvar[j]*varMax;
          tmp   = Zmean[k]*(1.0-Zmean[k])/varCur - 1.0;
          alpha = Zmean[k]*tmp;
          beta  = (1.0-Zmean[k])*tmp;

          pdfFlag = MDUCcalcPDF_1D(opts->method[0], opts->n1, x, Zmean[k], varCur,
                                   pdf, pdfBound);

          meanPDF = calcLinMeanY(0, 1, alpha, beta, opts->n1, x, Y,
                                 pdf, pdfBound, pdfFlag);

          absDiff[k][j] = fabs(Zmean[k] - meanPDF); 

          if (Zmean[k] != 0.0)
            error[k][j] = fabs(Zmean[k] - meanPDF)/Zmean[k];
          else
            error[k][j] = 0.0;

          if (error[k][j] > maxErr)
            {
              maxErr     = error[k][j];
              maxErrMean = Zmean[k];
              maxErrVar  = varCur;
            }

          if (error[k][j]*100 > 1.0)
            nPtsG1++;
        } // end loop over variances
    } // end loop over means

  printf("Error greater that 1.0 in %i of %i domain (%6.4f percent)\n", nPtsG1,
         nMean*nVar, (double)(nPtsG1)/nMean/nVar*100);
  printf("Maximum error computed: %6.4f %%\n", maxErr*100);
  printf("Zmean = %8.6f, Zvar = %8.6f\n", maxErrMean, maxErrVar);

  // write the solution to a file
  FILE *pFile = fopen("errorTest1D.dat", "w");
  
  fprintf(pFile, "# Zmean Zvar*Zmean*(1-Zmean)  error(%%) \n");

  for (k=0; k<nMean; k++)
    {
      for (j=0; j<nVar; j++)
        {
          fprintf(pFile, "%8.6f %8.6f %15.6e %15.6e\n", Zmean[k], Zvar[j],
                  absDiff[k][j], error[k][j]);
        } // end loop over variances
      fprintf(pFile, "\n");
    } // end loop over means

  fclose(pFile);

  free(Zmean);
  free(Zvar);
  free(pdf);
  for (k=0; k<nMean; k++)
    {
      free(error[k]);
      free(absDiff[k]);
    }
  free(error);
  free(absDiff);
}
/* ------------------------------------------------------------------------- */

void errorTest2D(optStruct *opts)
{
  int i, j, k, q, r, s, row, col;
  int n1, n2, nZ1, nV1, nZ2, nV2, nV3, nPtsG1 = 0;
  double dZ1, dZ2, dV1, dV2, dV3;
  double mean, Z1max, Z2max, V1max, V2max, V3max, Z1, Z2, V1, V2, V3, maxQ;
  double S, Q, tmp;

  // arrays for all possible means/variances
  double *Z1mean, *Z1var, *Z2mean, *Z2var, *Zvar;

  double beta1, beta2, beta3;   // beta parameters

  // arrays for holding error and absolute difference
  double mean1PDF, mean2PDF, meanPDF;
  double maxErr[3], mErr_Z1[3], mErr_Z2[3], mErr_V1[3], mErr_V2[3], mErr_V3[3];
  double *****error,*****error1, *****error2, *****absD, *****absD1, *****absD2;

  // for computing the pdf
  int n[2];
  double **x, **dx, *x1, *x2, *Y, *Y1, *Y2; // for pdf grid
  int pdfFlag, **pdfBound;
  double **pdf;

  n[0] = n1 = opts->n1;
  n[1] = n2 = opts->n2;

  for (i=0; i<3; i++) maxErr[i] = 0.0;

  nZ1 = 5;
  nV1 = 5;
  nZ2 = 5;
  nV2 = 5;
  nV3 = 5;

  dZ1 = 1.0/(nZ1-1);
  dV1 = 1.0/(nV1-1);
  dZ2 = 1.0/(nZ2-1);
  dV2 = 1.0/(nV2-1);
  dV3 = 1.0/(nV3-1);

  // allocate memory for sampling space
  Z1mean  = (double*) calloc(nZ1, sizeof(double));
  Z1var   = (double*) calloc(nV1, sizeof(double));
  Z2mean  = (double*) calloc(nZ2, sizeof(double));
  Z2var   = (double*) calloc(nV2, sizeof(double));
  Zvar    = (double*) calloc(nV3, sizeof(double));

  // allocate memory for error storage
  error  = (double*****) calloc(nZ1, sizeof(double****));
  error1 = (double*****) calloc(nZ1, sizeof(double****));
  error2 = (double*****) calloc(nZ1, sizeof(double****));
  absD   = (double*****) calloc(nZ1, sizeof(double****));
  absD1  = (double*****) calloc(nZ1, sizeof(double****));
  absD2  = (double*****) calloc(nZ1, sizeof(double****));

  for (k=0; k<nZ1; k++)
    {
      error[k]  = (double****) calloc(nV1, sizeof(double***));
      error1[k] = (double****) calloc(nV1, sizeof(double***));
      error2[k] = (double****) calloc(nV1, sizeof(double***));
      absD[k]   = (double****) calloc(nV1, sizeof(double***));
      absD1[k]  = (double****) calloc(nV1, sizeof(double***));
      absD2[k]  = (double****) calloc(nV1, sizeof(double***));

      for (j=0; j<nV1; j++)
        {
          error[k][j]  = (double***) calloc(nZ2, sizeof(double**));
          error1[k][j] = (double***) calloc(nZ2, sizeof(double**));
          error2[k][j] = (double***) calloc(nZ2, sizeof(double**));
          absD[k][j]   = (double***) calloc(nZ2, sizeof(double**));
          absD1[k][j]  = (double***) calloc(nZ2, sizeof(double**));
          absD2[k][j]  = (double***) calloc(nZ2, sizeof(double**));

          for (r=0; r<nZ2; r++)
            {
              error[k][j][r]  = (double**) calloc(nV2, sizeof(double*));
              error1[k][j][r] = (double**) calloc(nV2, sizeof(double*));
              error2[k][j][r] = (double**) calloc(nV2, sizeof(double*));
              absD[k][j][r]   = (double**) calloc(nV2, sizeof(double*));
              absD1[k][j][r]  = (double**) calloc(nV2, sizeof(double*));
              absD2[k][j][r]  = (double**) calloc(nV2, sizeof(double*));

              for (q=0; q<nV2; q++)
                {
                  error[k][j][r][q]  = (double*) calloc(nV3, sizeof(double));
                  error1[k][j][r][q] = (double*) calloc(nV3, sizeof(double));
                  error2[k][j][r][q] = (double*) calloc(nV3, sizeof(double));
                  absD[k][j][r][q]   = (double*) calloc(nV3, sizeof(double));
                  absD1[k][j][r][q]  = (double*) calloc(nV3, sizeof(double));
                  absD2[k][j][r][q]  = (double*) calloc(nV3, sizeof(double));
                }
            }
        }
    }

  // allocate memory for pdf computation
  x        = (double**) calloc(2, sizeof(double*));
  dx       = (double**) calloc(2, sizeof(double*));
  pdf      = (double**) calloc(n[1], sizeof(double*));
  pdfBound = (int**) calloc(n[1], sizeof(int*));
  for (k=0; k<n[1]; k++)
    {
      pdf[k]      = (double*) calloc(n[0], sizeof(double));
      pdfBound[k] = (int*) calloc(2, sizeof(int));
    }

  // intialise pdf grid
  initGrid(n[0], 0.0, 1.0, opts->nPiv1,opts->gPiv1, opts->gPerc1, &x[0], &dx[0]);
  initGrid(n[1], 0.0, 1.0, opts->nPiv2,opts->gPiv2, opts->gPerc2, &x[1], &dx[1]);

  x1 = x[0];
  x2 = x[1];

  // for computing the mean, compute the function Y
  fprintf(stderr, "getting function for mean calc...\n");
  Y  = (double*) calloc(n2*n1, sizeof(double));
  Y1 = (double*) calloc(n2*n1, sizeof(double));
  Y2 = (double*) calloc(n2*n1, sizeof(double));
  for (k=0; k<n2; k++)
    {
      row = k*n1;
      for (j=0; j<n1; j++)
        {
          if (opts->transform == 0)
            {
              Y[row+j] = x1[j] + x2[k];
            }
          else if (opts->transform == 1)
            {
              Y[row+j]  = x1[j] + x2[k]*(1.0-x1[j]);
              Y1[row+j] = x1[j];
              Y2[row+j] = x2[k]*(1.0-x1[j]);
            }
          else if (opts->transform == 2)
            {
              Y[row+j] = x1[j] + x2[k] + x2[k]*(1.0-x1[j]-x2[k]);
            }
        }
    }
  

  // initialise sampling space grids
  Z1mean[0] = 0.0;
  for (k=1; k<nZ1; k++) Z1mean[k] = Z1mean[k-1] + dZ1;

  Z1var[0] = 0.0;
  for (j=1; j<nV1; j++) Z1var[j] = Z1var[j-1] + dV1;

  Z2mean[0] = 0.0;
  for (r=1; r<nZ2; r++) Z2mean[r] = Z2mean[r-1] + dZ2;

  Z2var[0] = 0.0;
  for (q=1; q<nV2; q++) Z2var[q] = Z2var[q-1] + dV2;

  Zvar[0] = 0.0;
  for (s=1; s<nV3; s++) Zvar[s] = Zvar[s-1] + dV3;

  for (k=0; k<nZ1; k++)
    {
      Z1 = Z1mean[k];
      fprintf(stderr, "Z1mean = %6.4f ", Z1);
      Z2max = 1.0-Z1;
      V1max = Z1*(1.0-Z1);
      for (j=0; j<nV1; j++)
        {
          // compute the current variance of Z1 (non-normalised)
          V1 = Z1var[j]*V1max;
          for (r=0; r<nZ2; r++)
            {
              Z2 = Z2max*Z2mean[r];
              // compute the maximum possible variance of Z1
              V2max = Z2*(1.0-Z2);

              /* since both means will not change in subsequent looping, compute
                 the parameters only dependent on mean here */
              mean = Z1 + Z2;
              S    = Z1*Z1 + Z2*Z2 + (1.0-mean)*(1.0-mean);

              for (q=0; q<nV2; q++)
                {
                  // compute the current variance of Z2 (non-normalised)
                  V2 = Z2var[q]*V2max;

                  // can compute here the maximum possible variance 3
                  //V3max = 1.0 - S - V1 - V2;
                  V3max = V1max + V2max + 2*pow(V1max*V2max, 0.5);

                  for (s=0; s<nV3; s++)
                    {
                      V3 = Zvar[s]*V3max;
                      
                      Q = V1 + V2 + V3; // variance parameter

                      tmp = (1.0 - S)/Q - 1.0;

                      beta1 = Z1*tmp;
                      beta2 = Z2*tmp;
                      beta3 = (1.0-mean)*tmp;
                      
                      pdfFlag = MDUCcalcPDF_2D(opts->method[0], n, x, Z1, V1,
                                               Z2, V2, V3, pdf, pdfBound,
                                               opts->transform);

                      // compute the mean
                      meanPDF = calcLinMeanY_2D(0, 1, n, x, Y,beta1,beta2,beta3,
                                                pdf, pdfBound, pdfFlag,
                                                opts->transform);

                      mean1PDF = calcLinMeanY_2D(0,1,n,x,Y1,beta1,beta2,beta3,
                                                 pdf, pdfBound, pdfFlag,
                                                 opts->transform);

                      mean2PDF = calcLinMeanY_2D(0,1,n,x,Y2,beta1,beta2,beta3,
                                                 pdf, pdfBound, pdfFlag,
                                                 opts->transform);

                      absD[k][j][r][q][s]  = fabs(mean-meanPDF);
                      absD1[k][j][r][q][s] = fabs(Z1-mean1PDF);
                      absD2[k][j][r][q][s] = fabs(Z2-mean2PDF);


                      if (mean > 0.0)
                        error[k][j][r][q][s]  = absD[k][j][r][q][s]/mean;
                      else 
                        error[k][j][r][q][s]  = 0.0;

                      if (Z1 > 0.0)
                        error1[k][j][r][q][s] = absD1[k][j][r][q][s]/Z1;
                      else
                        error1[k][j][r][q][s] = 0.0;

                      if (Z2 > 0.0)
                        error2[k][j][r][q][s]= absD2[k][j][r][q][s]/Z2;
                      else
                        error2[k][j][r][q][s] = 0.0;

                      if (error[k][j][r][q][s] > maxErr[0])
                        {
                          maxErr[0]  = error[k][j][r][q][s];
                          mErr_Z1[0] = Z1;
                          mErr_Z2[0] = Z2;
                          mErr_V1[0] = V1;
                          mErr_V2[0] = V2;
                          mErr_V3[0] = V3;
                        }

                      if (error1[k][j][r][q][s] > maxErr[1])
                        {
                          maxErr[1]  = error1[k][j][r][q][s];
                          mErr_Z1[1] = Z1;
                          mErr_Z2[1] = Z2;
                          mErr_V1[1] = V1;
                          mErr_V2[1] = V2;
                          mErr_V3[1] = V3;
                        }

                      if (error2[k][j][r][q][s] > maxErr[2])
                        {
                          maxErr[2]  = error2[k][j][r][q][s];
                          mErr_Z1[2] = Z1;
                          mErr_Z2[2] = Z2;
                          mErr_V1[2] = V1;
                          mErr_V2[2] = V2;
                          mErr_V3[2] = V3;
                        }


                    } // end loop over variance 3
                } // end loop over variance 2
            } // end loop over mean 2
          fprintf(stderr,".");
        } // end loop over variance 1
      fprintf(stderr,"\n");
    } // end loop over mean 1

  fprintf(stderr, "\n\n");

  fprintf(stderr,"Maximum error statistics for 2D pdf:\n");
  fprintf(stderr,"------------------------------------------------------------");
  fprintf(stderr,"---------------\n");
  fprintf(stderr,"Error in mean 1: %6.4f %%\n", maxErr[1]*100);
  fprintf(stderr,"  Z1mean = %6.4f, Z1var = %8.6f\n", mErr_Z1[1], mErr_V1[1]);
  fprintf(stderr,"  Z2mean = %6.4f, Z2var = %8.6f\n", mErr_Z2[1], mErr_V2[1]);
  fprintf(stderr,"                   Zvar = %8.6f\n", mErr_V3[1]);

  fprintf(stderr,"\nError in mean 2: %6.4f %%\n", maxErr[2]*100);
  fprintf(stderr,"  Z1mean = %6.4f, Z1var = %8.6f\n", mErr_Z1[2], mErr_V1[2]);
  fprintf(stderr,"  Z2mean = %6.4f, Z2var = %8.6f\n", mErr_Z2[2], mErr_V2[2]);
  fprintf(stderr,"                   Zvar = %8.6f\n", mErr_V3[2]);

  fprintf(stderr,"\nError in total mean: %6.4f %%\n", maxErr[0]*100);
  fprintf(stderr,"  Z1mean = %6.4f, Z1var = %8.6f\n", mErr_Z1[0], mErr_V1[0]);
  fprintf(stderr,"  Z2mean = %6.4f, Z2var = %8.6f\n", mErr_Z2[0], mErr_V2[0]);
  fprintf(stderr,"                   Zvar = %8.6f\n", mErr_V3[0]);



  // write the solution to a file
 /*  FILE *pFile = fopen("pdfRegime2D.dat", "w"); */
  
/*   fprintf(pFile, "# Z1mean Z1var Z2mean Z2var Zvar  beta1  beta2 beta3 regime \n"); */

/*   for (k=0; k<nZ1; k++) */
/*     { */
/*       for (j=0; j<nV1; j++) */
/*         { */
/*           for (r=0; r<nZ2; r++) */
/*             { */
/*               for (q=0; q<nV2; q++) */
/*                 { */
/*                   for (s=0; s<nV3; s++) */
/*                     { */
/*                       fprintf(pFile, "%8.6f %8.6f %8.6f %8.6f %8.6f %15.6e %15.6e %15.6e\n", Z1mean[k], Z1var[j], Z2mean[r], Z2var[q], Zvar[s], beta1[k][j][r][q][s], beta2[k][j][r][q][s], beta3[k][j][r][q][s]); */
/*                     } // end loop over variance 3 */
/*                   fprintf(pFile, "\n"); */
/*                 } // end loop over variance 2 */
/*               fprintf(pFile, "\n"); */
/*             } // end loop over mean 2 */
/*           fprintf(pFile, "\n"); */
/*         } // end loop over variance 1 */
/*       fprintf(pFile, "\n"); */
/*     } // end loop over mean 1 */

  free(Y);
  free(Y1);
  free(Y2);
  free(Z1mean);
  free(Z1var);
  free(Z2mean);
  free(Z2var);
  free(Zvar);
}
/* ------------------------------------------------------------------------- */

void regime(optStruct *opts)
{
  int i, j, k;
  int nMean, nVar, nPtsG1 = 0;
  double dMean, dVar;
  double varMax, varCur;
  double **alpha, **beta, tmp;
  double *Zmean, *Zvar; // arrays for all possible means/variances
  int    **regime;

  nMean = opts->n1;
  nVar  = opts->n2;

  dMean = 1.0/(nMean-1);
  dVar = 1.0/(nVar-1);

  Zmean   = (double*) calloc(nMean, sizeof(double));
  Zvar    = (double*) calloc(nVar, sizeof(double));
  alpha   = (double**) calloc(nMean, sizeof(double*));
  beta    = (double**) calloc(nMean, sizeof(double*));
  regime  = (int**) calloc(nMean, sizeof(int*));
  for (k=0; k<nMean; k++)
    {
      alpha[k]   = (double*) calloc(nVar, sizeof(double));
      beta[k]    = (double*) calloc(nVar, sizeof(double));
      regime[k]  = (int*) calloc(nVar, sizeof(int));
    }

  Zmean[0] = 0.0;
  for (k=1; k<nMean; k++) Zmean[k] = Zmean[k-1] + dMean;

  Zvar[0] = 0.0;
  for (j=1; j<nVar; j++) Zvar[j] = Zvar[j-1] + dVar;
  
  for (k=0; k<nMean; k++)
    {
      varMax = Zmean[k]*(1.0-Zmean[k]);
      for (j=0; j<nVar; j++)
        {
          varCur = Zvar[j]*varMax;
          tmp   = Zmean[k]*(1.0-Zmean[k])/varCur - 1.0;
          alpha[k][j] = Zmean[k]*tmp;
          beta[k][j]  = (1.0-Zmean[k])*tmp;

          if (Zmean[k] == 0.0)
            {
              alpha[k][j] = 0.0;
              beta[k][j]  = 0.0;
            }

          //if (fabs(alpha[k][j] - 1) < 1e-9) alpha[k][j] = 1.0;
          //if (fabs(beta[k][j] - 1) < 1e-9) beta[k][j] = 1.0;

          if (alpha[k][j] <= 0.0 || beta[k][j] <= 0.0)
            regime[k][j] = -3;
          else if ( (fabs(alpha[k][j] - 1) < 1e-9) && (fabs(beta[k][j] - 1) < 1e-9) )
            regime[k][j] = 0;
          else if ( (alpha[k][j] < 1.0 && beta[k][j] < 1.0) )
            regime[k][j] = -2;
          else if ( (alpha[k][j] < 1.0 && beta[k][j] >= 1.0) || (alpha[k][j] == 1.0 && beta[k][j] > 1.0) )
            regime[k][j] = 1;
          else if ( (alpha[k][j] > 1.0 && beta[k][j] <= 1.0) || (alpha[k][j] == 1.0 && beta[k][j] < 1.0) )
            regime[k][j] = 1;
          else
            regime[k][j] = 2;

        } // end loop over variances
    } // end loop over means

  // write the solution to a file
  FILE *pFile = fopen("pdfRegime1D.dat", "w");
  
  fprintf(pFile, "# Zmean Zvar*Zmean*(1-Zmean)  alpha  beta  regime \n");

  for (k=0; k<nMean; k++)
    {
      for (j=0; j<nVar; j++)
        {
          fprintf(pFile, "%8.6f %8.6f %15.6e %15.6e %i\n", Zmean[k], Zvar[j],
                  alpha[k][j], beta[k][j], regime[k][j]);
        } // end loop over variances
      fprintf(pFile, "\n");
    } // end loop over means

  fclose(pFile);

  free(Zmean);
  free(Zvar);
}
/* ------------------------------------------------------------------------- */

void regime2D(optStruct *opts)
{
  int i, j, k, q, r, s;
  int nZ1, nV1, nZ2, nV2, nV3, nPtsG1 = 0;
  double dZ1, dZ2, dV1, dV2, dV3;
  double mean, Z1max, Z2max, V1max, V2max, V3max, Z1, Z2, V1, V2, V3, maxQ;
  double S, Q, tmp;

  // arrays for all possible means/variances
  double *Z1mean, *Z1var, *Z2mean, *Z2var, *Zvar;

  // arrays for holding beta parameters
  double *****beta1, *****beta2, *****beta3;

  nZ1 = 11;
  nV1 = 11;
  nZ2 = 11;
  nV2 = 11;
  nV3 = 11;

  dZ1 = 1.0/(nZ1-1);
  dV1 = 1.0/(nV1-1);
  dZ2 = 1.0/(nZ2-1);
  dV2 = 1.0/(nV2-1);
  dV3 = 1.0/(nV3-1);

  Z1mean  = (double*) calloc(nZ1, sizeof(double));
  Z1var   = (double*) calloc(nV1, sizeof(double));
  Z2mean  = (double*) calloc(nZ2, sizeof(double));
  Z2var   = (double*) calloc(nV2, sizeof(double));
  Zvar    = (double*) calloc(nV3, sizeof(double));

  beta1 = (double*****) calloc(nZ1, sizeof(double****));
  beta2 = (double*****) calloc(nZ1, sizeof(double****));
  beta3 = (double*****) calloc(nZ1, sizeof(double****));

  for (k=0; k<nZ1; k++)
    {
      beta1[k] = (double****) calloc(nV1, sizeof(double***));
      beta2[k] = (double****) calloc(nV1, sizeof(double***));
      beta3[k] = (double****) calloc(nV1, sizeof(double***));

      for (j=0; j<nV1; j++)
        {
          beta1[k][j] = (double***) calloc(nZ2, sizeof(double**));
          beta2[k][j] = (double***) calloc(nZ2, sizeof(double**));
          beta3[k][j] = (double***) calloc(nZ2, sizeof(double**));

          for (r=0; r<nZ2; r++)
            {
              beta1[k][j][r] = (double**) calloc(nV2, sizeof(double*));
              beta2[k][j][r] = (double**) calloc(nV2, sizeof(double*));
              beta3[k][j][r] = (double**) calloc(nV2, sizeof(double*));

              for (q=0; q<nV2; q++)
                {
                  beta1[k][j][r][q] = (double*) calloc(nV3, sizeof(double));
                  beta2[k][j][r][q] = (double*) calloc(nV3, sizeof(double));
                  beta3[k][j][r][q] = (double*) calloc(nV3, sizeof(double));
                }
            }
        }
    }

  Z1mean[0] = 0.0;
  for (k=1; k<nZ1; k++) Z1mean[k] = Z1mean[k-1] + dZ1;

  Z1var[0] = 0.0;
  for (j=1; j<nV1; j++) Z1var[j] = Z1var[j-1] + dV1;

  Z2mean[0] = 0.0;
  for (r=1; r<nZ2; r++) Z2mean[r] = Z2mean[r-1] + dZ2;

  Z2var[0] = 0.0;
  for (q=1; q<nV2; q++) Z1var[q] = Z1var[q-1] + dV2;

  Zvar[0] = 0.0;
  for (s=1; s<nV3; s++) Zvar[s] = Zvar[s-1] + dV3;
  
  for (k=0; k<nZ1; k++)
    {
      // compute the maximum possible variance of Z1
      Z1 = Z1mean[k];
      fprintf(stderr, "Z1mean = %6.4f ", Z1);
      Z2max = 1.0-Z1;
      V1max = Z1*(1.0-Z1);
      for (j=0; j<nV1; j++)
        {
          // compute the current variance of Z1 (non-normalised)
          V1 = Z1var[j]*V1max;
          for (r=0; r<nZ2; r++)
            {
              Z2 = Z2max*Z2mean[r];
              // compute the maximum possible variance of Z1
              V2max = Z2*(1.0-Z2);

              /* since both means will not change in subsequent looping, compute
                 the parameters only dependent on mean here */
              mean = Z1 + Z2;
              S = Z1*Z1 + Z2*Z2 + (1.0-mean)*(1.0-mean);

              for (q=0; q<nV2; q++)
                {
                  // compute the current variance of Z2 (non-normalised)
                  V2 = Z2var[q]*V2max;

                  // can compute here the maximum possible variance 3
                  V3max = 1.0 - S - V1 - V2;
                  for (s=0; s<nV3; s++)
                    {
                      V3 = Zvar[s]*V3max;
                      
                      Q = V1 + V2 + V3; // variance parameter

                      tmp = (1.0 - S)/Q - 1.0;

                      beta1[k][j][r][q][s] = Z1*tmp;
                      beta2[k][j][r][q][s] = Z2*tmp;
                      beta3[k][j][r][q][s] = (1.0-mean)*tmp;

                    } // end loop over variance 3
                } // end loop over variance 2
            } // end loop over mean 2
          fprintf(stderr,".");
        } // end loop over variance 1
      fprintf(stderr,"\n");
    } // end loop over mean 1

  fprintf(stderr,"writing data to file...\n");
  // write the solution to a file
  FILE *pFile = fopen("pdfRegime2D.dat", "w");
  
  fprintf(pFile, "# Z1mean Z1var Z2mean Z2var Zvar  beta1  beta2 beta3 regime \n");

  for (k=0; k<nZ1; k++)
    {
      for (j=0; j<nV1; j++)
        {
          for (r=0; r<nZ2; r++)
            {
              for (q=0; q<nV2; q++)
                {
                  for (s=0; s<nV3; s++)
                    {
                      fprintf(pFile, "%8.6f %8.6f %8.6f %8.6f %8.6f %15.6e %15.6e %15.6e\n", Z1mean[k], Z1var[j], Z2mean[r], Z2var[q], Zvar[s], beta1[k][j][r][q][s], beta2[k][j][r][q][s], beta3[k][j][r][q][s]);
                    } // end loop over variance 3
                  fprintf(pFile, "\n");
                } // end loop over variance 2
              fprintf(pFile, "\n");
            } // end loop over mean 2
          fprintf(pFile, "\n");
        } // end loop over variance 1
      fprintf(pFile, "\n");
    } // end loop over mean 1


  free(Z1mean);
  free(Z1var);
  free(Z2mean);
  free(Z2var);
  free(Zvar);
  fprintf(stderr, "goodbye\n");
}
/* ------------------------------------------------------------------------- */
