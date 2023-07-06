#ifndef MDUC_PLOT_H
#define MDUC_PLOT_H

#include "mducStructs.h"

/* struct to hold option read from command line */
struct OPTSTRUCT {
  /* ----- for getting command line options and storing ----- */
  int  nOpts;          // total number of options
  int  nFiles;         // number of files to process
  int  nPlotSpecies;   // number of species to plot
  int  nInterp;        // number of points on grid to interpolate to
  int  readFlow;       // whether file is flow file
  int  readStats;      // whether file is stats file
  int  readChi;        // whether file is chi file
  char *args;          // command line arguments
  char **plotFileName; // list of filenames to process
  char **plotSpecies;  // optional input for species to plot
  char suffix[10];     // file suffix to use
  char termType[10];   // type of terminal for gnuplot
  char postTerm[100];  // to write after terminal
  char **compFileName; 
  char wLine[10];
  int  outTerm;          // whether or not to set output

  int  *plotIndex;     // stores the tab file index of each species to plot
  int  *nonZero;       // temporary storage for determining non-zero species

  /* ----- plot flags ----- */
  int  genScript;  // generate gnuplot script for plotting solution file
  int  pdfScript;  // generate gnuplot script for plotting PDF
  int  chiScript;  // generate gnuplot script for plotting scalar dissipation
  int  genMovie;   // use settings for getting frames for movie (terminal, etc.)
  int  compFile;
  int  condStat;
  int  plotLine;   // use line plots (1D) instead of just points
  int  plotDim;    // number of dimensions to be plotted
  int  plotMulti;  // generate multiplot
  int  plotRow;    // row to plot if single row from 2D solution is desired
  int  plotCol;    // column to plot if single row from 2D solution is desired
  int  plotH;      // plot the enthalpy
  int  plotChi;    // plot the scalar dissipation rate
  int  plotPDF;    // plot the global pdf
  int  margBeta;   // plot the marginal beta pdfs
  int  plotY;      // plot mass fractions
  int  plotProg;   // plot the progress variable, as defined by species
  int  plotMisc;   // plot miscellaneous variables
  int  logScale;   // use semi-log scale for pdf
  int  smooth;     // use gnuplot interpolate function
  int  addVar; 
  int  unscaleStat;
  int  fileRange, fileMin, fileMax;
  int transpose;  
  
  /* ----- global plot formatting options ----- */  
  int xIndex;     // for specifying the x co-ordinate in the tabulated file
  char rcIndex[50]; // info for plotting row/column of 2D solution
  double x1min, x1max, x2min, x2max;
  double *Ymin, *Ymax, Tmin, Tmax, pdfMin[3], pdfMax[3];
  double chiMin[4], chiMax[4];
  double *miscMin, *miscMax;
  double RPM;
  int  iT, iHTOT; // location of plots
  
  int x1Range, x2Range, Trange, Yrange, range;
  int chiRange[4], pdfRange, miscRange;

  /* ----- post-processing variables ----- */
  int    nProgSpecies;
  char   **progSpecies;
  double *progVar, progY;
  int sae;
};
typedef struct OPTSTRUCT optStruct;

// function write the flamelet solution read in to a tabulated format
void writeTabFile(char *fileName, optStruct *opts, ioStruct *io,
                  driveStruct *driver, timeStruct *time, flameStruct *flame,
                  flowStruct *flow);

// read options from the command line
void readOptions(int argc, char **argv, optStruct *opts);
void readList(int n, int iStart, char ***var, char **argv);

// generate a file name with different suffix
void createOutputFileName(char *prefix, char *fileName, char *ext,
                          char *outputFile);


#endif /* MDUC_PLOT_H */
