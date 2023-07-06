#ifndef MDUC_GNUPLOTSCRIPT_H
#define MDUC_GNUPLOTSCRIPT_H

#include "mducPlot.h"

/* struct to hold plot options */
struct PLOTSTRUCT {
  /* ----- flags for plotting ----- */
  int  plotDim;    // number of dimensions to be plotted
  
  /* ----- information for plot formatting ----- */  
  int    *plotIndex;     // stores the tab file index of each species to plot
  char   **plotSpecies;  // stores the name of species for titles
  int    nMisc;                         // number of additional variables
  char   **miscNames;                   // additional variable names
  int    chiPlot[4];                    // whether chi was read in
  int    pdfPlot;                       // whether pdf was read in
  int    HPlot;                         // whether total enthalpy was read
  int    iY, iT, iHTOT, iPDF, iMisc, iProg;
  int    iChi[4], iTheta[4];
  int    iFuel, iCO2, iCO, iH2O; // location of plots

  /* ----- storage for various labels, options ----- */
  char   xLabel[20], yLabel[20];
  char   rcIndex[50], postTerm[100];
  char   timeLabel[200];
};
typedef struct PLOTSTRUCT plotStruct;


void genPlotScript(char *fileName, optStruct *opts, plotStruct *gp);

void genMultiPlotScript(char *fileName, char *cFile, optStruct *opts,
                        plotStruct *gp);

void genSAE_Script2D(char *fileName, char *cFile, optStruct *opts,
                     plotStruct *gp);

void genFlowScript(char *fileName, optStruct *opts, plotStruct *gp);

void genPDFscript(char *fileName, optStruct *opts, plotStruct *gp);

void genChiScript(char *fileName, optStruct *opts, plotStruct *gp);

#endif /* MDUC_GNUPLOTSCRIPT_H */
