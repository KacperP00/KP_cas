#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mducPlot.h"
#include "gnuplotScript.h"
#include "fileIO.h"
#include "parser.h"
#include "init.h"
#include "destroy.h"
#include "statsUtils.h"
#include "betaPDF.h"
#include "utils.h"

void printHelp()
{
  printf("\n");
  printf("options :\n");
  printf("  -h                       # print help\n");
  printf("  -i <input files>         # flamelet file to plot\n");
  printf("  -o <output suffix>       # suffix for output file (def: .tab)\n");
  printf("  -n <min> <max>           # specify a range of file numbers\n");
  printf("  -g                       # generate gnuplot script\n");
  printf("  -f                       # use flow file format\n");
  printf("  -p                       # stats file\n");
  printf("  -M                       # generate plot scripts for movie\n");
  printf("\nplot script options :\n");
  printf("  -t <terminal>            # terminal type for gnuplot (def: X11)\n");
  printf("  -s <species list>        # list species to plot (def: maj prod)\n");
  printf("  -r <row>                 # plot row of 2D solution\n");
  printf("  -c <col>                 # plot column of 2D solution\n");
  printf("  -H                       # plot total enthalpy\n");
  printf("  -X                       # plot scalar dissipation rate\n");
  printf("  -A                       # plot all variables in file\n");
  printf("  -L                       # use semi-log scale\n");
  printf("  --RPM <n>                # engine speed\n");
  printf("  --smooth <n>             # use gnuplot interpolate function\n");
  printf("  --all-species            # plot all species!\n");
  printf("  --no-species\n");
  printf("  --progress-var <val> <species list>\n");
  printf("  --x1-range  <min> <max>  # specify the x-axis range\n");
  printf("  --x2-range  <min> <max>  # specify the x-axis range\n");
  printf("  --Y-range   <min> <max>  # specify the mass fraction plot range\n");
  printf("  --T-range   <min> <max>  # specify the temperature plot range\n");
  printf("  --chi-range <min> <max>  # specify the chi plot range\n");
  exit(0);
}
/* -------------------------------------------------------------------------- */

void setDefaultOptions(optStruct *opts)
{
  opts->plotFileName = NULL;
  strcpy(opts->termType, "X11");
  strcpy(opts->suffix, "tab");
  strcpy(opts->wLine, " ");
  opts->nPlotSpecies = -1;
  
  opts->outTerm = 0;       // plot to screen rather than file
  opts->plotDim = 1;       // 1D plot
  opts->xIndex  = 1;       // x axis 
  opts->plotRow = -2;      // default to plotting 2D field
  opts->plotCol = -2;      //     : same

  // use autoscaling
  opts->Trange    = 0; 
  opts->Yrange    = 0;
  opts->chiRange[1] = opts->chiRange[3] = opts->chiRange[2] = 0;
  opts->fileRange = opts->pdfRange  = 0;

  // set min and max values
  opts->x1min = 1e10;
  opts->x1max = 0.0;
  opts->x2min = 1e10;
  opts->x2max = 1.0;
  opts->Ymin = opts->Ymax = NULL;
  opts->miscMin = opts->miscMax = NULL;
  opts->chiMin[1] = opts->chiMin[2] = opts->chiMin[3] = 0.0;
  opts->chiMax[1] = opts->chiMax[2] = opts->chiMax[3] = -1;
  opts->Tmin = 300;
  opts->Tmax = -1;
  opts->pdfMin[0] = opts->pdfMin[1] = opts->pdfMin[2] = 0.0;
  opts->pdfMax[0] = opts->pdfMax[1] = opts->pdfMax[2] = 1.0;

  opts->readFlow  = 0;  // not a flow file
  opts->readStats = 0;  // not a stats file
  opts->readChi   = 0;  // not a chi file
  opts->nInterp   = -1; // don't interpolate anything
  opts->plotLine  = 0;  // don't use line plots
  opts->plotH     = 0;  // don't plot enthalpy
  opts->plotChi   = 0;  // don't plot scalar dissipation
  opts->plotY     = 1;  // plot major species
  opts->plotProg  = 0;  // don't plot progress variable
  opts->plotMisc  = 0;  // don't plot miscellaneous variables
  opts->pdfScript = 0;  // don't generate pdf script
  opts->chiScript = 0;  // don't generate scalar dissipation
  opts->genScript = 0;  // don't generate plot script
  opts->genMovie  = 0;  // don't do a movie
  opts->plotMulti = 1;
  opts->margBeta  = 0;
  opts->smooth    = 0;
  opts->addVar    = 0;
  opts->progVar   = NULL;
  opts->logScale  = 0;
  opts->unscaleStat = 0;
  opts->compFile = 0;
  opts->sae      = 0;
  opts->transpose = 0;
  opts->RPM      = -1;
}
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  int i, j, k, m, n, r, row, col, nTmp, iEnd, iFile, count;
  double minVal, maxVal, threshold = 1e-20;
  optStruct opts;  // struct to store options read in from command line
  plotStruct *gp;  // struct to store plotting options for each file

  char file[128], *head, *tail;  // storage for beginning and end of string
  char chiFileName[128], buffer[128], prefix[128], suffix[20];
  char cFile[256], cPre[128], cSuf[128];

  FILE *pFile = NULL;
  
  // declare flamelet structs so that internal functions can be used
  driveStruct  driver;
  ioStruct     io;
  timeStruct   time;
  flameStruct  flame;
  flowStruct   flow;
  engineStruct engine;

  int iCO2, iH2O, iCO, iFuel;
  double Yprog = 0.2807314; // progress variable based on CO2 + H2O  

  // set default options for plot script
  setDefaultOptions(&opts);
  // interpret options read passed from command line
  readOptions(argc, argv, &opts);
  if (opts.plotFileName[0] == NULL) {
    printf("Error: no input file specified\n");
    exit(2);
  }

  // if only a range of files was specified
  if (opts.fileRange) {
    // store prefix/suffix
    head = strtok(opts.plotFileName[0], ".");
    strcpy(prefix, head);
    head = strtok(NULL, "\n");
    strcpy(suffix, head);
    
    if (opts.compFile) {
      head = strtok(opts.compFileName[0], ".");
      strcpy(cPre, head);
      head = strtok(NULL, "\n");
      strcpy(cSuf, head);
    }
    
    opts.nFiles = 0;
    free(opts.plotFileName[0]);
    fprintf(stderr, "checking file range:\n");
    // loop over files to see if they exist
    for (i=opts.fileMin,iFile=0,count=0; i<=opts.fileMax; i++) {	  
      // current file name
      sprintf(buffer, "%05i", i);
      strcpy(file, prefix);
      strcat(file, buffer);
      strcat(file, ".");
      strcat(file, suffix);
      /* fprintf(stderr, "%s\n", file); */
      
      if (opts.compFile) {
        strcpy(cFile, cPre);
        strcat(cFile, buffer);
        strcat(cFile, ".");
        strcat(cFile, opts.suffix);
      }
      
      pFile = fopen(file, "r");
      if (pFile) {
        opts.nFiles++;
        fprintf(stderr, ".");
        opts.plotFileName = (char**) realloc(opts.plotFileName,
                                             opts.nFiles*sizeof(char*));
        opts.plotFileName[iFile] = (char*) calloc(200, sizeof(char));
        strcpy(opts.plotFileName[iFile], file);
        
        if (opts.compFile) {
          opts.compFileName=(char**) realloc(opts.compFileName,
                                             opts.nFiles*sizeof(char*));
          opts.compFileName[iFile] = (char*) calloc(200, sizeof(char));
          strcpy(opts.compFileName[iFile], cFile);
        }
        
        iFile++;
        fclose(pFile);
      }
      else {
        fprintf(stderr, "-");
      }
      if (count%80 == 0)
        fprintf(stderr, "\n");
      count++;
    } // end file loop
    fprintf(stderr, "\n %i files found\n\n", opts.nFiles);
  } // end check for files

  // allocate memory for plotting options
  gp = (plotStruct*) calloc(opts.nFiles, sizeof(plotStruct));
  
  // set global plot options for all files if plotting
  if (opts.genScript) {
    if (opts.genMovie)
      strcpy(opts.termType, "png");
    
    // set options for X11 terminal type so window stays open
    if (!strcmp(opts.termType, "X11"))
      strcpy(opts.postTerm, " persist");
    else
      strcpy(opts.postTerm, "");
    
    // set whether the terminal prints to a file
    if (!strcmp(opts.termType, "png") || !strcmp(opts.termType, "pdf")
        || !strcmp(opts.termType, "jpeg")
        || !strcmp(opts.termType, "png")
        || !strcmp(opts.termType, "postscript") )
      opts.outTerm = 1;
    
    // set the size of the window if a jpeg or png
    if (!strcmp(opts.termType, "png") || !strcmp(opts.termType, "jpeg"))
      /* strcat(opts.postTerm, " size 1280,800 enhanced font \"Times, 14\""); */
      strcat(opts.postTerm, " size 1280,800 enhanced");

    // set to enhanced for pdf and postscript
    if (!strcmp(opts.termType, "pdf") || !strcmp(opts.termType, "postscript"))
      strcat(opts.postTerm, " enhanced");
    
    if (opts.plotLine) strcpy(opts.wLine, " w l ");
    
    // set the row if plotting a single row
    strcpy(opts.rcIndex, "");
    if (opts.plotRow >= 0)
      sprintf(opts.rcIndex, "index %i", opts.plotRow-1);
    
    // settings if plotting a single column
    if (opts.plotCol >= 0) {
      if (opts.plotRow >= 0) {
        printf("Error: plot either rows or columns, not both\n");
        exit(2);
      }
      opts.xIndex = 2;
      sprintf(opts.rcIndex,"every ::%i::%i",opts.plotCol-1,opts.plotCol-1);
    }
  } // end set of global plot script values
  
  // simple script to plot a PDF; no need to process a file
  if (opts.pdfScript) {
    for (i=0; i<opts.nFiles; i++) {
      strcpy(gp[i].xLabel, "x1");
      strcpy(gp[i].yLabel, "x2");
      
      if (opts.plotDim == 1)
        gp[i].iY = 2;
      else if (opts.plotDim == 2) {
        if (opts.plotRow >= 0 || opts.plotCol >= 0) {
          opts.plotDim = 1;
          strcpy(gp[i].xLabel, gp[i].yLabel);
        }
        gp[i].iY = 3;
      }
      genPDFscript(opts.plotFileName[i], &opts, &gp[i]);
    }
    
    return 0;
  }
  
  // generate the tabulated files for plotting
  fprintf(stderr, "Processing %i files:\n", opts.nFiles);
  for (i=0; i<opts.nFiles; i++) {
    fprintf(stderr, "\t%s", opts.plotFileName[i]);
    
    // set flamelet defaults
    setDefaults(&io, &driver, &time, &flame, &flow, NULL);
    
    // read in flamelet file
    readFlamelet(opts.plotFileName[i],&io,&driver,&time,&flame,&flow);
    
    // read flow data if specified
    if (opts.readFlow) {
      driver.CHI_MODE = -1;
      readFlowData(opts.plotFileName[i], &driver, &flow);
    }
    
    // check to see if no transformation for stats
    if (!strcmp(io.transZ2, "none")) driver.TRANS_Z2 = 0;
    
    // check to see if beta marginal should be plotted
    if (io.mean[1] > 0.0 && io.mean[2] > 0.0) opts.margBeta = 1;
    
    // read in the leading file name and extension based on the first .
    strcpy(buffer, opts.plotFileName[i]);
    head = strtok(buffer, ".");
    strcpy(file, head);
    
    // unscale the pdfs if desired
    if (opts.unscaleStat) {
      if (opts.readFlow) {
        unscalePDF(0, driver.TRANS_Z2, flame.n1, flame.n2, flame.x1,
                   flame.x2, &flow.pdf);
        if (!opts.pdfRange) {
          maxVal = getMaxValue(flame.n1, flame.n2, 0, 1, flow.pdf);
          if (maxVal > opts.pdfMax[0]) opts.pdfMax[0] = maxVal;
        }
      }
      else if (opts.readStats) {
        for (k=0; k<io.nMisc; k++) {
          unscalePDF(0, driver.TRANS_Z2, flame.n1, flame.n2, flame.x1,
                     flame.x2, &io.miscIn[k]);
          // track the max value
          if (!opts.pdfRange) {
            maxVal = getMaxValue(flame.n1,flame.n2,0,1,io.miscIn[k]);
            if (maxVal > opts.pdfMax[0]) opts.pdfMax[0] = maxVal;
          }
        }
      } // end readStats
    } // end unscale pdf

    if (opts.plotProg) {
      // compute the progress variable
      opts.progVar = (double*) calloc(flame.n1*flame.n2, sizeof(double));
      
      for (k=0; k<flame.n2; k++) {
        r = k*flame.n1;
        for (j=0; j<flame.n1; j++) {
          for (m=0; m<opts.nProgSpecies; m++) {
            nTmp = getSpeciesIndex(io.nInputs, io.inputSpecies,
                                   opts.progSpecies[m]);
            opts.progVar[r+j] += io.initY[nTmp][r+j];
          } // end loop over progress variable species
          
          opts.progVar[r+j] /= opts.progY;
        } // end loop over columns
      } // end loop over rows
    } // end compute progress variable
    
    // write the tabulated file for plotting
    writeTabFile(file, &opts, &io, &driver, &time, &flame, &flow);
    
    // restart loop here if no plot script requested
    if (!opts.genScript)
      continue;

    /* ----- generate a gnuplot plotting script if specified ----- */
    // set axis labels depending on dimensions
    if (driver.H_DIM == 0) { // no enthalpy direction
      if (driver.TRANS_Z2 == 1) {
        strcpy(gp[i].xLabel, "Z1");
        strcpy(gp[i].yLabel, "Z2/(1-Z1)");
      }
      else if (driver.TRANS_Z2 == 2) {
        strcpy(gp[i].xLabel, "Z1+Z2");
        strcpy(gp[i].yLabel, "Z2/(Z1+Z2)");
      }
      else {                
        strcpy(gp[i].xLabel, "Z1");
        strcpy(gp[i].yLabel, "Z2");
      }
    }
    else if (driver.H_DIM == 1) {
      strcpy(gp[i].xLabel, "H");
      strcpy(gp[i].yLabel, "Z");
    }
    else if (driver.H_DIM == 2) {
      strcpy(gp[i].xLabel, "Z");
      strcpy(gp[i].yLabel, "H");
    }
    
    // set the dimension to plot
    if (flame.nDim == 1) {
      gp[i].plotDim = 1;
      gp[i].iY = 2;
    }
    else if (flame.nDim == 2) {
      opts.plotDim = 2;
      gp[i].iY = 3;
      if (opts.plotRow >= 0 || opts.plotCol >= 0) {
        gp[i].plotDim = 1;
        strcpy(gp[i].xLabel, gp[i].yLabel);
      }
      else
        gp[i].plotDim = 2;
    }

    // default to plotting major products
    if (opts.plotY == 1) {
      opts.nPlotSpecies = 3;
      gp[i].plotSpecies=(char**) calloc(opts.nPlotSpecies,sizeof(char*));
      for (m=0; m < opts.nPlotSpecies; m++)
        gp[i].plotSpecies[m] = (char*) calloc(25, sizeof(char));
      
      strcpy(gp[i].plotSpecies[0], "CO2");
      strcpy(gp[i].plotSpecies[1], "H2O");
      strcpy(gp[i].plotSpecies[2], "CO");
    }
    else if (opts.plotY == 2) { // plot all non-zero species
      opts.nPlotSpecies = 0;
      gp[i].plotSpecies = (char**) calloc(io.nInputs, sizeof(char*));
      for (m=0; m < io.nInputs; m++)
        gp[i].plotSpecies[m] = (char*) calloc(25, sizeof(char));
      opts.nonZero = (int*) calloc(io.nInputs, sizeof(int));
      
      // check for zero species
      for (m=0; m < io.nInputs; m++) {
        for (k=0; k < flame.n2; k++) {
          r = k*flame.n1;
          for (j=0; j < flame.n1; j++) {
            if (io.initY[m][r+j] > threshold)
              opts.nonZero[m] = 1;
          } // end loop over columns
        } // end loop over rows
      } // end loop over inputs
      
      for (m=0, k=0; m < io.nInputs; m++) {
        if (opts.nonZero[m] == 1) {
          // copy name of species across
          strcpy(gp[i].plotSpecies[k], io.inputSpecies[m]);
          opts.nPlotSpecies++;
          k++;
        }
      } // end loop over inputs
      free(opts.nonZero);
    } // end set default to plot all species
    else {
      gp[i].plotSpecies=(char**) calloc(opts.nPlotSpecies,sizeof(char*));
      for (m=0; m < opts.nPlotSpecies; m++) {
        gp[i].plotSpecies[m] = (char*) calloc(25, sizeof(char));
        strcpy(gp[i].plotSpecies[m], opts.plotSpecies[m]);
      }
    }

    // keep track of min/max of axes
    if (flow.phi1min < opts.x1min) opts.x1min = flow.phi1min;
    if (flow.phi1max > opts.x1max) opts.x1max = flow.phi1max;
    if (flow.phi2min < opts.x2min) opts.x2min = flow.phi2min;
    if (flow.phi2max > opts.x2max) opts.x2max = flow.phi2max;
    
    // allocate memory for min/max if necessary
    if (opts.Ymin == NULL)
      opts.Ymin = (double*) calloc(opts.nPlotSpecies, sizeof(double));
    if (opts.Ymax == NULL)
      opts.Ymax = (double*) calloc(opts.nPlotSpecies, sizeof(double));
    
    // get the index of the species to plot
    gp[i].plotIndex = (int*) calloc(opts.nPlotSpecies, sizeof(int));
    for (m=0; m<opts.nPlotSpecies; ++m) {
      // get the index of the species to plot
      gp[i].plotIndex[m] = getSpeciesIndex(io.nInputs, io.inputSpecies, 
                                           gp[i].plotSpecies[m]);

      if (!opts.Yrange) {
        opts.Ymin[m]=0.0;
        // get the max value of the species
        maxVal=getMaxValue(flame.n1,flame.n2,0,1,io.initY[gp[i].plotIndex[m]]);
        
        // track the global maximum of this species/all files
        if (maxVal > opts.Ymax[m])
          opts.Ymax[m] = maxVal;
      }
    } // end loop over plot species

    // keep track of the maximum values of temperature and chi
    if (!opts.Trange && io.initT != NULL) {
      maxVal = getMaxValue(flame.n1, flame.n2, 0, 1, io.initT);
      if (maxVal > opts.Tmax)
        opts.Tmax = maxVal;
      minVal = getMinValue(flame.n1, flame.n2, 0, 1, io.initT);
      if (minVal < opts.Tmin)
        opts.Tmin = minVal;
    }
    
    // location of temperature
    if (opts.readFlow || opts.readStats || opts.readChi)
      gp[i].iT = gp[i].iY - 1;
    else
      gp[i].iT = gp[i].iY + io.nInputs;
    iEnd = gp[i].iT;
    
    // location of total enthalpy
    gp[i].HPlot = 0;
    if (flame.H != NULL) {
      gp[i].HPlot = 1;
      gp[i].iHTOT = iEnd = gp[i].iT + 1;
    }
    
    // location of transformed scalar dissipation rates
    gp[i].chiPlot[1] = gp[i].chiPlot[2] = gp[i].chiPlot[3] = 0;
    if (flame.chiA != NULL) {
      gp[i].chiPlot[1] = 1;
      if (!opts.chiRange[1]) {
        maxVal = getMaxValue(flame.n1, flame.n2, 0, 1, flame.chiA);
        if (maxVal > opts.chiMax[1])
          opts.chiMax[1] = maxVal;
      }
      gp[i].iChi[1] = iEnd = iEnd + 1;
    }
    
    if (flame.chiB != NULL) {
      gp[i].chiPlot[3] = 1;
      if (!opts.chiRange[3]) {
        maxVal = getMaxValue(flame.n1, flame.n2, 0, 1, flame.chiB);
        if (maxVal > opts.chiMax[3])
                opts.chiMax[3] = maxVal;
      }
      gp[i].iChi[3] = iEnd = iEnd + 1;
    }
    
    if (flame.chiC != NULL) {
      gp[i].chiPlot[2] = 1;
      if (!opts.chiRange[2]) {
        maxVal = getMaxValue(flame.n1, flame.n2, 0, 1, flame.chiC);
        if (maxVal > opts.chiMax[2])
          opts.chiMax[2] = maxVal;
      }
      gp[i].iChi[2] = iEnd = iEnd + 1;
    }
    
    // location of scalar dissipation rates          
    if (flow.chi1 != NULL) {
      gp[i].chiPlot[1] = 1;
      if (!opts.chiRange[1]) {
        maxVal = getMaxValue(flame.n1, flame.n2, 0, 1, flow.chi1);
        if (maxVal > opts.chiMax[1])
          opts.chiMax[1] = maxVal;
      }
      gp[i].iChi[1] = iEnd = iEnd + 1;
    }
    
    if (flow.chi12 != NULL) {
      gp[i].chiPlot[3] = 1;
      if (!opts.chiRange[3]) {
        minVal = getMinValue(flame.n1, flame.n2, 0, 1, flow.chi12);
        if (maxVal < opts.chiMin[3])
          opts.chiMin[3] = minVal;
        
        maxVal = getMaxValue(flame.n1, flame.n2, 0, 1, flow.chi12);
        if (maxVal > opts.chiMax[3])
          opts.chiMax[3] = maxVal;
      }
      gp[i].iChi[3] = iEnd = iEnd + 1;
    }
    
    if (flow.chi2 != NULL) {
      gp[i].chiPlot[2] = 1;
      if (!opts.chiRange[2]) {
        maxVal = getMaxValue(flame.n1, flame.n2, 0, 1, flow.chi2);
        if (maxVal > opts.chiMax[2])
          opts.chiMax[2] = maxVal;
      }
      gp[i].iChi[2] = iEnd = iEnd + 1;
    }
    
    // if flow file, read in pdf
    gp[i].pdfPlot = 0;
    if (flow.pdf != NULL) {
      gp[i].pdfPlot = 1;
      maxVal = getMaxValue(flame.n1, flame.n2, 0, 1, flow.pdf);
      if (maxVal > opts.pdfMax[0])
        opts.pdfMax[0] = maxVal;
      gp[i].iPDF = iEnd = iEnd + 1;
    }
    
    // find the location of the start of the miscellaneous variables
    gp[i].iMisc = iEnd = iEnd + 1;
    
    // allocate the memory for variable names
    gp[i].nMisc = io.nMisc;
    gp[i].miscNames = (char**) calloc(io.nMisc, sizeof(char*));
    if (opts.miscMin == NULL)
      opts.miscMin = (double*) calloc(io.nMisc, sizeof(double));
    if (opts.miscMax == NULL)
      opts.miscMax = (double*) calloc(io.nMisc, sizeof(double));
    
    // store the names of the miscellaneous variables
    for (j=0; j < io.nMisc; j++) {
      gp[i].miscNames[j] = (char*) calloc(50, sizeof(char));
      strcpy(gp[i].miscNames[j], io.miscNames[j]);
      
      if (opts.readStats) stripHeader(gp[i].miscNames[j]);
      
      // track the min/max of each variable
      minVal = getMinValue(flame.n1, flame.n2, 0, 1, io.miscIn[j]);
      maxVal = getMaxValue(flame.n1, flame.n2, 0, 1, io.miscIn[j]);
      if (minVal < opts.miscMin[j]) opts.miscMin[j] = minVal;
      if (maxVal > opts.miscMax[j]) opts.miscMax[j] = maxVal;
    }
    
    gp[i].iProg = gp[i].iMisc + io.nMisc;
    
    // set the time variables
    if (opts.readStats)
      sprintf(gp[i].timeLabel, "time = %9.6f", time.t);
    else if (opts.sae)
      sprintf(gp[i].timeLabel, "CAD = %6.1f aTDC", time.t*9000+420-720);
    else if (opts.RPM > 0)
      sprintf(gp[i].timeLabel, "CAD = %6.1f aTDC", time.t*opts.RPM*6+420-720);
    else
      sprintf(gp[i].timeLabel, "time = %8.4f ms", time.t*1e3);
    
    // destroy the memory from this flamelet
    destroyFlame(&flame);
    destroyFlow(&flow);
    destroyIO(&io);
    free(opts.progVar);
    fprintf(stderr, "\n");
  } // end loop over files for reading 
  
  // generate gnuplot scripts if specified
  if (opts.genScript) {
    fprintf(stderr, "\ngenerating plot script(s) ");
    
    // set min/max of axes to nearest decimal
    if (opts.x1min < 0.0) {
      sprintf(buffer, "%3.1f", opts.x1min - 0.05);
      opts.x1min = atof(buffer);
    }
    else
      opts.x1min = 0.0;
    
    if (opts.x1max > 1.0) {
      sprintf(buffer, "%3.1f", opts.x1max + 0.05);
      opts.x1max = atof(buffer);
    }
    else
      opts.x1max = 1.0;
 
    if (opts.x2min < 0.0) {
      sprintf(buffer, "%3.1f", opts.x2min - 0.05);
      opts.x2min = atof(buffer);
    }
    else
      opts.x2min = 0.0;

    if (opts.x2max > 1.0) {
      sprintf(buffer, "%3.1f", opts.x2max + 0.05);
      opts.x2max = atof(buffer);
    }
    else
      opts.x2max = 1.0;

    // set up specifically for movie
    if (opts.genMovie) {
      fprintf(stderr, "for movie ");

      if (!opts.Yrange) {
        opts.Yrange = 1;
        for (m=0; m<opts.nPlotSpecies; m++) {
          opts.Ymin[m] = 0.0; // always set to zero for now
          opts.Ymax[m] = fmin(1.0, 1.1*opts.Ymax[m]);
          if (fabs(opts.Ymin[m] - opts.Ymax[m]) < 1.0e-6)
            opts.Ymax[m] += 0.1;
        }

        /* if one dimensional plot, get the max of all species plotted
           and store in the first value of the array */
        if (opts.plotDim == 1) {
          for (m=1; m<opts.nPlotSpecies; m++)
            if (opts.Ymax[m] > opts.Ymax[0])
              opts.Ymax[0] = opts.Ymax[m];
        }
      }

      if (!opts.Trange) {
        opts.Trange = 1;
        opts.Tmin = fmax(300.0, floor(0.9*opts.Tmin));
        opts.Tmax = fmin(3000.0, ceil(1.1*opts.Tmax));
      }

      if (!opts.chiRange[1]) {
        opts.chiRange[1] = 1;
        opts.chiMin[1] = 0.0;
        if (opts.chiMax[1] <= 1.0)
          opts.chiMax[1] = 1.1*opts.chiMax[1];
        else
          opts.chiMax[1] = ceil(1.1*opts.chiMax[1]);

        if (opts.chiMin[1] == opts.chiMax[1]) opts.chiMax[1] += 1.0;
      }

      if (!opts.chiRange[3]) {
        opts.chiRange[3] = 1;
        opts.chiMin[3] = 0.0;
        if (opts.chiMax[3] <= 1.0)
          opts.chiMax[3] = 1.1*opts.chiMax[3];
        else
          opts.chiMax[3] = ceil(1.1*opts.chiMax[3]);

        if (opts.chiMin[3] == opts.chiMax[3]) opts.chiMax[3] += 1.0;
      }

      if (!opts.chiRange[2]) {
        opts.chiRange[2] = 1;
        opts.chiMin[2] = 0.0;
        if (opts.chiMax[2] <= 1.0)
          opts.chiMax[2] = 1.1*opts.chiMax[2];
        else
          opts.chiMax[2] = ceil(1.1*opts.chiMax[2]);

        if (opts.chiMin[2] == opts.chiMax[2]) opts.chiMax[2] += 1.0;
      }

      // set the time variables
      if (opts.readStats)
        sprintf(gp[i].timeLabel, "time = %9.6f", time.t);
      else if (opts.sae)
        sprintf(gp[i].timeLabel, "CAD = %6.1f aTDC", time.t*9000+420-720);
      else if (opts.RPM > 0)
        sprintf(gp[i].timeLabel, "CAD = %6.1f aTDC", time.t*opts.RPM*6+420-720);
      else
        sprintf(gp[i].timeLabel, "time = %8.4f ms", time.t*1e3);

      if (!opts.pdfRange) {
        opts.pdfRange = 1;
        for (i=0; i<3; i++) {
          if (opts.pdfMax[i] <= 1.0)
            opts.pdfMax[i] = 1.1*opts.pdfMax[i];
          else
            opts.pdfMax[i] = ceil(1.1*opts.pdfMax[i]);
        }
      }
    }
    fprintf(stderr, "...\n");

    // make sure the min value for log scale is not zero
    if (opts.logScale) {
      opts.pdfRange = 1;
      if (opts.unscaleStat) {
        opts.pdfMin[0] = fmax(1.0e-1, opts.pdfMin[0]);
        opts.pdfMin[1] = fmax(1.0e-1, opts.pdfMin[1]);
        opts.pdfMin[2] = fmax(1.0e-1, opts.pdfMin[2]);
        
        if (opts.readStats) {
          for (k=0; k<io.nMisc; k++)
            opts.miscMin[k] = fmax(1.0e-1, opts.miscMin[k]);
        }
      }
      else {
        opts.pdfMin[0] = fmax(1.0e-7, opts.pdfMin[0]);
        opts.pdfMin[1] = fmax(1.0e-6, opts.pdfMin[1]);
        opts.pdfMin[2] = fmax(1.0e-6, opts.pdfMin[2]);
      }
      
      
      // set up specifically for movie
      if (opts.genMovie) {
        fprintf(stderr, "for movie ");
        
        if (!opts.Yrange) {
          opts.Yrange = 1;
          for (m=0; m<opts.nPlotSpecies; m++) {
            opts.Ymin[m] = 0.0; // always set to zero for now
            opts.Ymax[m] = fmin(1.0, 1.1*opts.Ymax[m]);
            if (fabs(opts.Ymin[m] - opts.Ymax[m]) < 1.0e-6)
              opts.Ymax[m] += 0.1;
          }
          
          /* if one dimensional plot, get the max of all species plotted
             and store in the first value of the array */
          if (opts.plotDim == 1) {
            for (m=1; m<opts.nPlotSpecies; m++)
              if (opts.Ymax[m] > opts.Ymax[0])
                opts.Ymax[0] = opts.Ymax[m];
          }
        }
        
        if (!opts.Trange) {
          opts.Trange = 1;
          opts.Tmin = fmax(300.0, floor(0.9*opts.Tmin));
          opts.Tmax = fmin(3000.0, ceil(1.1*opts.Tmax));
        }
        
        if (!opts.chiRange[1]) {
          opts.chiRange[1] = 1;
          opts.chiMin[1] = 0.0;
          if (opts.chiMax[1] <= 1.0)
            opts.chiMax[1] = 1.1*opts.chiMax[1];
          else
            opts.chiMax[1] = ceil(1.1*opts.chiMax[1]);
          
          if (opts.chiMin[1] == opts.chiMax[1]) opts.chiMax[1] += 1.0;
        }
        
        if (!opts.chiRange[3]) {
          opts.chiRange[3] = 1;
          opts.chiMin[3] = 0.0;
          if (opts.chiMax[3] <= 1.0)
            opts.chiMax[3] = 1.1*opts.chiMax[3];
          else
            opts.chiMax[3] = ceil(1.1*opts.chiMax[3]);
          
          if (opts.chiMin[3] == opts.chiMax[3]) opts.chiMax[3] += 1.0;
        }
        
        if (!opts.chiRange[2])  {
          opts.chiRange[2] = 1;
          opts.chiMin[2] = 0.0;
          if (opts.chiMax[2] <= 1.0)
            opts.chiMax[2] = 1.1*opts.chiMax[2];
          else
            opts.chiMax[2] = ceil(1.1*opts.chiMax[2]);
          
          if (opts.chiMin[2] == opts.chiMax[2]) opts.chiMax[2] += 1.0;
        }
        
        if (!opts.pdfRange) {
          opts.pdfRange = 1;
          for (i=0; i<3; i++)
            {
              if (opts.pdfMax[i] <= 1.0)
                opts.pdfMax[i] = 1.1*opts.pdfMax[i];
              else
                opts.pdfMax[i] = ceil(1.1*opts.pdfMax[i]);
            }
        }
        
        if (opts.readStats) {
          for (k=0; k<io.nMisc; k++)
            opts.miscMin[k] = fmax(1.0e-8, opts.miscMin[k]);
        }
      }
    }
    
    for (i=0; i<opts.nFiles; i++) {
      // read in the leading file name and extension based on the first .
      strcpy(file, opts.plotFileName[i]);
      head = strtok(file, ".");
      strcpy(file, head);
      strcat(file, ".");
      strcat(file, opts.suffix);

      if (opts.compFile)
        strcpy(cFile, opts.compFileName[i]);
      else
        cFile == NULL;

      if (opts.plotMulti) {
        if (opts.readFlow)
          genFlowScript(file, &opts, &gp[i]);
        else if (opts.readStats)
          genPDFscript(file, &opts, &gp[i]);
        else if (opts.readChi)
          genChiScript(file, &opts, &gp[i]);
        else if (opts.sae)
          genSAE_Script2D(file, cFile, &opts, &gp[i]);
        else
          genMultiPlotScript(file, cFile, &opts, &gp[i]);
      }
      else
        genPlotScript(file, &opts, &gp[i]);
    }
  }

  if (opts.Ymin) free(opts.Ymin);
  if (opts.Ymax) free(opts.Ymax);
  if (opts.miscMin) free(opts.miscMin);
  if (opts.miscMax) free(opts.miscMax);

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
void writeTabFile(char *fileName, optStruct *opts, ioStruct *io,
                  driveStruct *driver, timeStruct *time, flameStruct *flame,
                  flowStruct *flow)
{
  int i, j, k, iC;
  int r, c, n1, n2;
  char file[100], *buffer;
  double maxVal;

  // add suffix
  strcpy(file, fileName);
  strcat(file, ".");
  strcat(file, opts->suffix);

  fprintf(stderr, " --> %s ...", file);
  r = (int)strlen(file) + 2;
  buffer = (char*) calloc(r, sizeof(char));
  for (i=0; i<r-2; i++) buffer[i] = ' ';
  buffer[r-1] = '\0';

  /* ----- write data file ----- */
  FILE *pFile = fopen(file, "w");
  
  /* ##### header */
  iC = 1;
  fprintf(pFile, "# ");
  if (driver->H_DIM == 1)
    fprintf(pFile, "%5s<%i> ","H", iC);
  else {
    if (driver->TRANS_Z2 == 2)
      fprintf(pFile, "%10s<%i> ","Z1+Z2", iC);
    else
      fprintf(pFile, "%5s<%i> ","Z1", iC);
  }
  iC++;
  
  if (flame->nDim == 2) {
    if (driver->H_DIM == 2)
      fprintf(pFile, "%5s<%i> ","H", iC);
    else if (driver->TRANS_Z2 == 1)
      fprintf(pFile, "%12s<%i> ","Z2/(1-Z1)", iC);
    else if (driver->TRANS_Z2 == 2)
      fprintf(pFile, "%12s<%i> ","Z2/(Z1-Z2)", iC);
    else
      fprintf(pFile, "%5s<%i> ","Z2", iC);
    iC++;
  }
  
  if (!opts->readFlow) {
    for (i=0; i < io->nInputs; i++, iC++)
      fprintf(pFile, "Y_%s<%i> ", io->inputSpecies[i], iC);
  }
  if (io->initT != NULL && !opts->readFlow) {
    fprintf(pFile, "Temp[K]<%i> ", iC); iC++;
  }
  if (flame->H != NULL) {
    fprintf(pFile, "TotalEnthalpy[J/kg]<%i> ", iC); iC++;
  }
  if (flame->rho != NULL) {
    fprintf(pFile, "Density[kg/m3]<%i> ", iC); iC++;
  }
  if (flame->chiA != NULL) { fprintf(pFile, "chiA[1/s]<%i> ",  iC); iC++; }
  if (flame->chiB != NULL) { fprintf(pFile, "chiB[1/s]<%i> ",  iC); iC++; }
  if (flame->chiC != NULL) { fprintf(pFile, "chiC[1/s]<%i> ",  iC); iC++; }
  if (flow->chi1  != NULL) { fprintf(pFile, "chi1[1/s]<%i> ",  iC); iC++; }
  if (flow->chi12 != NULL) { fprintf(pFile, "chi12[1/s]<%i> ", iC); iC++; }
  if (flow->chi2  != NULL) { fprintf(pFile, "chi2[1/s]<%i> ",  iC); iC++; }
  if (flow->pdf   != NULL) { fprintf(pFile, "pdf<%i> ",        iC); iC++; }
  if (flame->nSootMom > 0) {
    fprintf(pFile, "soot-number-density_M00[1/m^3]<%i> ", iC); iC++;
    fprintf(pFile, "soot-volume_M10[m^3/m^3]<%i> ", iC);       iC++;
    fprintf(pFile, "soot-surface_M01[m^2/m^3]<%i> ", iC);      iC++;
    if (flame->nSootMom > 4) {
      fprintf(pFile, "soot-vol-var_M20[m^6/m^3]<%i> ", iC);        iC++;
      fprintf(pFile, "soot-vol-surf-covar_M11[m^5/m^3]<%i> ", iC); iC++;
      fprintf(pFile, "soot-surf-var_M02[m^4/m^3]<%i> ", iC);       iC++;
    }
    fprintf(pFile, "soot-delta-func_N0[1/m^3]<%i> ", iC); iC++;
  }
  if (!opts->readFlow) {
    for (i=0; i < io->nMisc; i++) {
      fprintf(pFile, "%s<%i> ", io->miscNames[i], iC);
      iC++;
    }
  }
  if (opts->progVar != NULL) {
    fprintf(pFile, "%s<%i> ", "progressVariable", iC); iC++;
  }
  fprintf(pFile, "\n");
  /* ##### end header */

  // this should be the default, but this will ensure that the loop below works
  if (flame->nDim == 1)
    flame->n2 = 1;
  
  if (opts->transpose) {
    fprintf(stderr, "transposing...\n");
    n1 = flame->n2;
    n2 = flame->n1;
  }
  else {
    n1 = flame->n1;
    n2 = flame->n2;
  }
  
  for (k=0; k<n2; k++) {
      
    /*  if (k%2 == 0) */
    /*         { */

    for (j=0; j<n1; j++) {
      if (opts->transpose) {
        r   = j*flame->n1;
        c   = k;
      }
      else {
        r   = k*flame->n1;
        c   = j;
      }

         

      /* if (opts->transpose) */
      /*             fprintf(pFile, "%6.4f ", flame->x2[k]); */
      /*           else */
      fprintf(pFile, "%6.4f ", flame->x1[j]);
      if (flame->nDim == 2) {
        /*  if (opts->transpose) */
        /*                 fprintf(pFile, "%6.4f ", flame->x1[j]); */
        /*               else */
        fprintf(pFile, "%6.4f ", flame->x2[k]);	  
      }
      if (!opts->readFlow) {
        for (i=0; i<io->nInputs; i++) {
          if (io->initY[i][r+c] >= 0.0)
            fprintf(pFile, "%15.6e ", io->initY[i][r+c]);
          else
            fprintf(pFile, "%15.6e ", 0.0);
        }
      }
      if (io->initT != NULL && !opts->readFlow)
        fprintf(pFile, "%15.6e ", io->initT[r+c]);
      if (flame->H != NULL)
        fprintf(pFile, "%15.6e ", flame->H[r+c]);
      if (flame->rho != NULL)
        fprintf(pFile, "%15.6e ", flame->rho[r+c]);
      if (flame->chiA != NULL)
        fprintf(pFile, "%15.6e ", flame->chiA[r+c]);
      if (flame->chiB != NULL)
        fprintf(pFile, "%15.6e ", flame->chiB[r+c]);
      if (flame->chiC != NULL)
        fprintf(pFile, "%15.6e ", flame->chiC[r+c]);
      if (flow->chi1 != NULL) {
        if (driver->TRANS_Z2 == 0 && (flame->x1[j] + flame->x2[j] > 1.0))
          fprintf(pFile, "%15.6e ", 1.0/0.0);
        else
          fprintf(pFile, "%15.6e ", flow->chi1[r+c]);
      }
      if (flow->chi12 != NULL) {
        if (driver->TRANS_Z2 == 0 && (flame->x1[j] + flame->x2[j] > 1.0))
          fprintf(pFile, "%15.6e ", 1.0/0.0);
        else
          fprintf(pFile, "%15.6e ", flow->chi12[r+c]);
      }
      if (flow->chi2 != NULL) {
        if (driver->TRANS_Z2 == 0 && (flame->x1[j] + flame->x2[j] > 1.0))
          fprintf(pFile, "%15.6e ", 1.0/0.0);
        else
          fprintf(pFile, "%15.6e ", flow->chi2[r+c]);
      }
      if (flow->pdf != NULL) {
        if (driver->TRANS_Z2 == 0 && (flame->x1[j] + flame->x2[j] > 1.0))
          fprintf(pFile, "%15.6e ", 1.0/0.0);
        else
          fprintf(pFile, "%15.6e ", flow->pdf[r+c]);
      }
      if (flame->nSootMom > 0) {
        fprintf(pFile, "%15.6e %15.6e %15.6e ", io->initSoot[0][j],
                io->initSoot[1][j], io->initSoot[2][j]);
        if (flame->nSootMom > 4) {
          fprintf(pFile, "%15.6e %15.6e %15.6e ", io->initSoot[3][j],
                  io->initSoot[4][j], io->initSoot[5][j]);
        }
        fprintf(pFile, "%15.6e ", io->initSoot[6][j]);
      }

      if (!opts->readFlow) {
        for (i=0; i<io->nMisc; i++) {
          if (opts->readStats || opts->readChi) {
            fprintf(pFile, "%15.6e ", io->miscIn[i][r+c]);
          }
          else
            fprintf(pFile, "%15.6e ", io->miscIn[i][r+c]);
        }
      }

      if (opts->progVar != NULL)
        fprintf(pFile, "%15.6e ", opts->progVar[r+c]);

      fprintf(pFile, "\n");

    } /* end loop over x1 */
    fprintf(pFile, "\n");
    if ( (opts->plotRow >= 0) || (opts->plotCol >= 0) )
      fprintf(pFile, "\n");
    /*             } */
  } /* end loop over x2 */
  fclose(pFile);

  /* ----- for flow/statistics files, write a file containing marginals ----- */
  if ((opts->readFlow || opts->readStats || opts->readChi)
      && flame->nDim == 2) {
    int n[2], nMarg = 1, iC, nChi = 1, *iChi1, *iChi2;
    double **x;
    int pdfBound[2];
    double **pdf, **pdfMarg1, **pdfMarg2, *beta1, *beta2;
    double **chiCond1, **chiCond2, *Z1condMean, *Z2condMean;
    char **chiNames;

    // if no pdf in flow file, return
    if (opts->readFlow && flow->pdf == NULL)
      return;

    x = (double**) calloc(2, sizeof(double*));
    n[0] = flame->n1;
    n[1] = flame->n2;
    x[0] = flame->x1;
    x[1] = flame->x2;

    if (opts->readStats) nMarg = io->nMisc;

    /* get the number of chis */
    nChi = 0;
    if (flow->chi1 != NULL) nChi = 1;
    /* find number of chi sets, assuming full set for every chi 1*/
    for (i=0; i<io->nMisc; i++)
      if (!strncmp(io->miscNames[i], "chi1", 4) &&
          strncmp(io->miscNames[i], "chi12", 5) ) 
        nChi++;

    double *chi1min, *chi1max;
    beta1    = (double*)  calloc(flame->n1, sizeof(double));
    beta2    = (double*)  calloc(flame->n2, sizeof(double));
    chi1min  = (double*)  calloc(flame->n1, sizeof(double));
    chi1max  = (double*)  calloc(flame->n1, sizeof(double));
    beta2    = (double*)  calloc(flame->n2, sizeof(double));
    pdfMarg1 = (double**) calloc(nMarg, sizeof(double*));
    pdfMarg2 = (double**) calloc(nMarg, sizeof(double*));
    for (i=0; i<nMarg; i++) {
      pdfMarg1[i] = (double*) calloc(flame->n1, sizeof(double));
      pdfMarg2[i] = (double*) calloc(flame->n2, sizeof(double));
    }
    chiCond1 = (double**) calloc(nChi, sizeof(double*));
    chiCond2 = (double**) calloc(nChi, sizeof(double*));
    iChi1    = (int*)     calloc(nChi, sizeof(int));
    iChi2    = (int*)     calloc(nChi, sizeof(int));
    /* chiNames = (char**)   calloc(nChi, sizeof(char*)); */
    for (i=0; i<nChi; i++) {
      /* chiNames[i] = (char*)   calloc(25, sizeof(char)); */
      chiCond1[i] = (double*) calloc(flame->n1, sizeof(double));
      chiCond2[i] = (double*) calloc(flame->n2, sizeof(double));
    }
    Z1condMean = (double*) calloc(flame->n2, sizeof(double));
    Z2condMean = (double*) calloc(flame->n1, sizeof(double));

    /* search for the chi indices */
    for (i=0, iC=1; i<io->nMisc; i++) {
      if (!strncmp(io->miscNames[i], "chi1", 4) &&
          strncmp(io->miscNames[i], "chi12", 5) ) {
        iChi1[iC] = i;
        iChi2[iC] = i+2; /* assumes chi12 in between, maybe bad */
        /* strncpy(chiNames[i], &io->miscNames[i][5], 10); */
        iC++;
      }
    }

    pdf = (double**) calloc(flame->n2, sizeof(double*));
    for (k=0; k<flame->n2; k++)
      pdf[k] = (double*) calloc(flame->n1, sizeof(double));

    if (opts->readFlow || opts->readChi) {
      convertArray1D_MD(flame->n1, flame->n2, 1, flow->pdf, &pdf);
      calcMarginalPDF(1, n, x, pdfMarg1[0], pdf);
      calcMarginalPDF(2, n, x, pdfMarg2[0], pdf);
      if (!opts->pdfRange) {
        maxVal = getMaxValue(flame->n1, 1, 0, 1, pdfMarg1[0]);
        if (maxVal > opts->pdfMax[1])
          opts->pdfMax[1] = maxVal;
        maxVal = getMaxValue(flame->n2, 1, 0, 1, pdfMarg2[0]);
        if (maxVal > opts->pdfMax[2])
          opts->pdfMax[2] = maxVal;
      }

      if (flow->chi1 != NULL) {
        for (j=0; j<flame->n1; j++) {
          chi1min[j] =  1e20;
          chi1max[j] = -1e20;
        }

        for (k=0; k<flame->n2; k++) {
          r = k*flame->n1;
          for (j=0; j<flame->n1; j++) {
            if (pdfMarg2[0][k] > 0.0)
              chiCond1[0][j]+=flow->chi1[r+j]*flow->pdf[r+j]/pdfMarg2[0][k];
            if (pdfMarg1[0][j] > 0.0)
              chiCond2[0][k]+=flow->chi2[r+j]*flow->pdf[r+j]/pdfMarg1[0][j];

            if (flame->x1[j] + flame->x1[k] < 1.0) {
              if (flow->chi1[r+j] < chi1min[j])
                chi1min[j] = flow->chi1[r+j];
              if (flow->chi1[r+j] > chi1max[j])
                chi1max[j] = flow->chi1[r+j];
            }
          }
        }

        for (j=0; j<flame->n1; j++) {
          if (chi1min[j] > 1e10)
            chi1min[j] = 0.0;
          if (chi1max[j] < -1e10)
            chi1max[j] = 0.0;
        }
      }
      if (nChi > 1) {
        for (i=1; i<nChi; i++) {
          for (k=0; k<flame->n2; k++) {
            r = k*flame->n1;
            for (j=0; j<flame->n1; j++) {
              if (!isnan(io->miscIn[iChi1[i]][r+j]) && pdfMarg2[0][k] > 0.0)
                chiCond1[i][j] += io->miscIn[iChi1[i]][r+j]*flow->pdf[r+j]/pdfMarg2[0][k];
              if (!isnan(io->miscIn[iChi2[i]][r+j]) && pdfMarg1[0][j] > 0.0)
                chiCond2[i][k] += io->miscIn[iChi2[i]][r+j]*flow->pdf[r+j]/pdfMarg1[0][j];
            }
          }

        }
      }
      /* compute the conditional mean that will be used as a reference */
      for (k=0; k<flame->n2; k++) {
        r = k*flame->n1;
        for (j=0; j<flame->n1; j++) {
          Z1condMean[k] += flame->x1[j]*flow->pdf[r+j]/pdfMarg2[0][k];
          Z2condMean[j] += flame->x2[k]*flow->pdf[r+j]/pdfMarg1[0][j];
        }
      }
    }
    if (opts->readStats) {
      for (i=0; i<io->nMisc; i++) {
        convertArray1D_MD(flame->n1, flame->n2, 1, io->miscIn[i], &pdf);
        calcMarginalPDF(1, n, x, pdfMarg1[i], pdf);
        calcMarginalPDF(2, n, x, pdfMarg2[i], pdf);
        if (!opts->pdfRange) {
          maxVal = getMaxValue(flame->n1, 1, 0, 1, pdfMarg1[i]);
          if (maxVal > opts->pdfMax[1])
            opts->pdfMax[1] = maxVal;
          maxVal = getMaxValue(flame->n2, 1, 0, 1, pdfMarg2[i]);
          if (maxVal > opts->pdfMax[2])
            opts->pdfMax[2] = maxVal;
        }
      }
    }

    // compute 1D betas for marginals
    if (opts->margBeta) {
      calcBetaPDF(flame->n1, flame->x1, io->mean[1], io->var[1], beta1,
                  pdfBound);
      calcBetaPDF(flame->n2, flame->x2, io->mean[2], io->var[2], beta2,
                  pdfBound);

      if (opts->unscaleStat) {
        unscalePDF(0,1,flame->n1,1,flame->x1,NULL,&beta1);
        unscalePDF(0,1,flame->n1,1,flame->x1,NULL,&beta2);
      }

      if (!opts->pdfRange) {
        maxVal = getMaxValue(flame->n1, 1, 0, 1, beta1);
        if (maxVal > opts->pdfMax[1])
          opts->pdfMax[1] = maxVal;
        maxVal = getMaxValue(flame->n2, 1, 0, 1, beta2);
        if (maxVal > opts->pdfMax[2])
          opts->pdfMax[2] = maxVal;
      }
    } // end 1D betas for marginals

      /* ----- write marginal data file ----- */
    strcpy(file, fileName);
    strcat(file, ".marg");
    pFile  = fopen(file, "w");
    fprintf(stderr, "\n\t%s --> %s ...", buffer, file);
      
    // write the header
    fprintf(pFile, "# Z1/Z2 ");
    for (i=0; i<nMarg; i++)
      fprintf(pFile, "%s ", io->miscNames[i]);
    fprintf(pFile, "\n");

    // write the marginal for Z1
    for (j=0; j<flame->n1; j++) {
      fprintf(pFile, "%6.4f ", flame->x1[j]);
      for (i=0; i<nMarg; i++)
        fprintf(pFile, "%15.6e ", pdfMarg1[i][j]);
      if (opts->margBeta)
        fprintf(pFile, "%15.6e ", beta1[j]);
      if (flow->chi1 != NULL)
        fprintf(pFile, "%15.6e ", chiCond1[0][j]);
      for (i=1; i<nChi; i++)
        fprintf(pFile, "%15.6e ", chiCond1[i][j]);
      if (flow->chi1 != NULL) {
        fprintf(pFile, "%15.6e ", chi1min[j]);
        fprintf(pFile, "%15.6e ", chi1max[j]);
      }
      fprintf(pFile, "\n");
    }

    // write another set for marginal Z2
    fprintf(pFile, "\n\n");
    for (k=0; k<flame->n2; k++) {
      fprintf(pFile, "%6.4f ", flame->x2[k]);
      for (i=0; i<nMarg; i++)
        fprintf(pFile, "%15.6e ", pdfMarg2[i][k]);
      if (opts->margBeta)
        fprintf(pFile, "%15.6e ", beta2[k]);
      if (flow->chi2 != NULL)
        fprintf(pFile, "%15.6e ", chiCond2[0][k]);
      for (i=1; i<nChi; i++)
        fprintf(pFile, "%15.6e ", chiCond2[i][k]);
      fprintf(pFile, "\n");
    }
    fclose(pFile);
      
    // free up local memory
    for (k=0; k<flame->n2; k++) free(pdf[k]);
    for (i=0; i<nMarg; i++) { free(pdfMarg1[i]); free(pdfMarg2[i]); }
    free(pdfMarg1); free(pdfMarg2);
    for (i=0; i<nChi; i++) { free(chiCond1[i]); free(chiCond2[i]); }
    free(chiCond1); free(chiCond2); free(iChi1); free(iChi2);
    free(Z1condMean); free(Z2condMean);
    free(x); free(beta1); free(beta2);
    free(pdf);
  }

  free(buffer);
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
  int i, j, k, nOpt, iOpt, iEnd, nTmp;
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
        case 't':
          strcpy(opts->termType, argv[iOpt+1]);
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
                readList(opts->nFiles, iOpt+1, &opts->plotFileName, argv);
              }
            break;
          }
        case 'C':
          {
            opts->compFile = 1;
            nTmp = iEnd - iOpt;
            readList(nTmp, iOpt+1, &opts->compFileName, argv);
            break;
          }
        case 'T':
          {
            opts->transpose = 1;
            break;
          }
        case 'o':
          {
            if (argv[iOpt+1][0] == '-')
              {
                printf("ERROR: definition of suffix required for '-o'\n");
                exit(2);
              }
            else
              {
                strcpy(opts->suffix, argv[iOpt+1]);
              }
            break;
          }
        case 'f':
          {
            opts->readFlow = 1;
            opts->plotY = 0;
            break; 
          }
        case 'g':
          {
            opts->genScript = 1;
            break; 
          }
        case 'S':
          {
            opts->plotMulti = 0;
            break; 
          }
        case 'l':
          {
            opts->plotLine = 1;
            break; 
          }
        case 'L':
          {
            opts->logScale = 1;
            break; 
          }
        case '2':
          {
            opts->plotDim = 2;
            break; 
          }
        case 'H':
          {
            opts->plotH = 1;
            break; 
          }
        case 'X':
          {
            opts->plotChi = 1;
            break; 
          }
        case 'r':
          {
            if (argv[iOpt+1] != NULL)
              {
                if (argv[iOpt+1][0] == '-')
                  opts->plotRow = -1;
                else
                  opts->plotRow = atoi(argv[iOpt+1]);
              }
            else
              opts->plotRow = -1;
            break;
          }
        case 'c':
          {
            if (argv[iOpt+1] != NULL)
              {
                if (argv[iOpt+1][0] == '-')
                  opts->plotCol = -1;
                else
                  opts->plotCol = atoi(argv[iOpt+1]);
              }
            else
              opts->plotCol = -1;
            break;
          }
        case 's':
          {
            opts->nPlotSpecies = iEnd - iOpt;
            readList(opts->nPlotSpecies, iOpt+1, &opts->plotSpecies, argv);
            opts->plotY = 3;
            break;
          }
        case 'M':
          {
            opts->genMovie  = 1;
            opts->genScript = 1;
            break; 
          }
        case 'A':
          {
            opts->plotY    = 1;
            opts->plotH    = 1;
            opts->plotChi  = 1;
            opts->plotMisc = 1;
            break; 
          }
        case 'P':
          {
            opts->pdfScript = 1;
            opts->genScript = 1;
            break; 
          }
        case 'p':
          {
            opts->readStats = 1;
            opts->plotY = 0;
            break; 
          }
        case 'n':
          {
            opts->fileMin = atoi(argv[iOpt+1]);
            opts->fileMax = atoi(argv[iOpt+2]);
            opts->fileRange = 1;
            break;
          }
        case '-': // has a longer word
          {
            if (!strcmp(argv[optLocs[i]], "--T-range")) {
              opts->Trange = 1;
              printf("%i %i\n", atoi(argv[iOpt+1]), atoi(argv[iOpt+2]));
              if (!isdigit(atoi(argv[iOpt+1])) && !isdigit(atoi(argv[iOpt+2])))
                {
                  opts->Tmin = atof(argv[iOpt+1]);
                  opts->Tmax = atof(argv[iOpt+2]);
                }
              else
                {
                  printf("Error: --Trange requires two values > zero\n");
                  exit(2);
                }
            }
            if (!strcmp(argv[optLocs[i]], "--Y-range")) {
              printf("%i %i\n", atoi(argv[iOpt+1]), atoi(argv[iOpt+2]));
              if (!isdigit(atoi(argv[iOpt+1])) && !isdigit(atoi(argv[iOpt+2])))
                {
                  opts->Ymin = (double*) calloc (1,sizeof(double));
                  opts->Ymax = (double*) calloc (1,sizeof(double));
                  opts->Ymin[0] = atof(argv[iOpt+1]);
                  opts->Ymax[0] = atof(argv[iOpt+2]);
                }
              else
                {
                  printf("Error: --Yrange requires two values\n");
                  exit(2);
                }
              opts->Yrange = 1;
            }
            
            if (!strcmp(argv[optLocs[i]], "--x1-range"))
              {
                
                opts->x1min = atof(argv[iOpt+1]);
                opts->x1max = atof(argv[iOpt+2]);
                opts->x1Range = 1;
              }
            if (!strcmp(argv[optLocs[i]], "--x2-range"))
              {
                opts->x2min = atof(argv[iOpt+1]);
                opts->x2max = atof(argv[iOpt+2]);
                opts->x2Range = 1;
              }
            if (!strcmp(argv[optLocs[i]], "--chi-range")) {
              printf("%i %i\n", atoi(argv[iOpt+1]), atoi(argv[iOpt+2]));
              if (!isdigit(atoi(argv[iOpt+1])) && !isdigit(atoi(argv[iOpt+2])))
                {
                  opts->chiMin[1] = atof(argv[iOpt+1]);
                  opts->chiMax[1] = atof(argv[iOpt+2]);
                  opts->chiMin[2] = atof(argv[iOpt+1]);
                  opts->chiMax[2] = atof(argv[iOpt+2]);
                }
              else
                {
                  printf("Error: --chirange requires two values > zero\n");
                  exit(2);
                }
              opts->chiRange[1] = 1;
              opts->chiRange[2] = 2;
            }
            if (!strcmp(argv[optLocs[i]], "--pdf-range"))
              {
                opts->pdfMin[0] = opts->pdfMin[1] = opts->pdfMin[2] = 
                  atof(argv[iOpt+1]);
                opts->pdfMax[0] = opts->pdfMax[1] = opts->pdfMax[2] =
                  atof(argv[iOpt+2]);
                opts->pdfRange = 1;
              }
            if (!strcmp(argv[optLocs[i]], "--interp"))
              {
                if (!isdigit(atoi(argv[iOpt+1])))
                  opts->nInterp = atoi(argv[iOpt+1]);
                else
                  {
                    printf("Error: --interp requires number of points\n");
                    exit(2);
                  }
              }
            if (!strcmp(argv[optLocs[i]], "--smooth"))
              {
                opts->smooth = atoi(argv[iOpt+1]);
              }
            if (!strcmp(argv[optLocs[i]], "--all-species"))
              {
                //opts->nPlotSpecies = -2;
                opts->plotY = 2;
              }
            if (!strcmp(argv[optLocs[i]], "--progress-var"))
              {
                opts->nProgSpecies = iEnd - iOpt - 1;
                opts->progY = atof(argv[iOpt+1]);
                readList(opts->nProgSpecies, iOpt+2, &opts->progSpecies, argv);
                opts->plotProg = 1;
                break;
              }
            if (!strcmp(argv[optLocs[i]], "--no-species"))
              {
                opts->nPlotSpecies = 0;
                opts->plotY = 0;
              }
            if (!strcmp(argv[optLocs[i]], "--add-variable"))
              {
                opts->addVar = 1;
              }
            if (!strcmp(argv[optLocs[i]], "--unscale-pdf"))
              {
                opts->unscaleStat = 1;
              }
            if (!strcmp(argv[optLocs[i]], "--sae"))
              {
                opts->sae = 1;
                opts->genMovie = opts->genScript =  1;
                opts->smooth = 5;
              }
            if (!strcmp(argv[optLocs[i]], "--chi"))
              {
                opts->readChi = 1;
                opts->plotY   = 0;
              }
            if (!strcmp(argv[optLocs[i]], "--RPM"))
              {
                opts->RPM = opts->chiMin[1] = atof(argv[iOpt+1]);
              }
            break;
          }
          
          
        }
    } // end loop over options
  
  free(optLocs);
}
/* -------------------------------------------------------------------------- */



