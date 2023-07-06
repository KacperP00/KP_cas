#ifndef MDUC_IO_STRUCT_H
#define MDUC_IO_STRUCT_H

#include <stdio.h>
#include "parser.h"

/* struct required for standalone calculations*/
typedef struct CHISTRUCT {
  int	 n1;
  int    n2;
  double P;
  double *x1;
  double tChi;
  double *chiRaw;
  double *pdfRaw;
  double *hSrcRaw;
  double *rhoDotRaw;
  double *TDotRaw;
  double *chiVal;
  double *pdfVal;
  double *hSrcVal;
  double *rhoDotVal;
  double *TDotVal;
  
}chiStruct;
//typedef struct CHISTRUCT chiStruct;

struct IOSTRUCT {
  /* ----- system settings ----- */
  int byteOrder;
  parseStruct *parser;
  int nOutputs;  // number of species to output
  int nInputs;   // number of species read in
  int freqType;  // frequency type
  int writeFreq; // time steps between write commands 
  int printFreq; // time steps between print to screen commands
  double tFreq;  // for time-based output frequencies
  double tWrite; // last time data file was written
  int debug;     // debug comment level
  
  /* ----- header data ----- */
  char title[128];      // title of test case
  char author[128];     // name, affiliation, etc.
  char date[128];       // date and time of run
  char mechanism[128];  // name of chemical mechanism used

  /* ----- solver options ----- */
  char pressureMode[64];
  int  pdfMode;        // for post-processing, pdf read from file or not
  char chiMode[64];        // for post-processing, chi read from file or not
  char timeMode[64] ;      // for post-processing, how to compute time series
  char transZ2[32];        // transform for 2D mixture fraction
  char sdialSolver[10];    // equation solver from SUNDIALS to use
  char chiInterp[256];     // interpolation method for scalar dissipation
  char heatLoss[32];       // specification of heat loss model to use
  char boundFilter[32];    // filtering treatment for physical domain bounds
  char viscUpdate[32];     // implicit/explicit viscosity treatment

  /* ----- file names and options ----- */
  int  readData;        // flag for restart
  int  readRates;	// read rates for restart
  int  readFlow;        // flag for reading flow data from a restart
  int  writeFlow;       // flag to define whether flow file is written
  int  writeEnsight;    // write Ensight files
  int  dumpErrorFile;  // dump a rif file when SUNDIALS error is detected
  int  errorFileWritten;// Prevents to dump same file multiple times
  int  initScalar;      // flag whether scalar file is already initialized
  char dataDir[128];    // optional specification of root data directory
  char inputFile[256];  // file read during initialisation 
  char inputData[256];  // file to read data from a restart 
  char inputRates[256];  // file to read data from a restart 
  char chiDataDir[256]; // directory name to chi data for standalone calculation 
  char engineFile[256]; // file to specify engine operating parameters
  char engineType[32];  // engine mode: HCCI, SI, SACI
  char inputFlow[256];  // file to read flow data for restart
  char outputData[256]; // base name of output file 
  char monitorFile[256];// file to write monitor
  char errorFile[256];  // file to write error comparison to
  char scalars[256];    // file to write scalar variables (also for restart) 
  char flowData[256];   // file to write flow variables to
  char statsData[256];  // file to write PDF and chi data to
                        // this is really for the flow solver
  char numStyle[64];    // define how output files should be numbered
  int  fileNumPrec;     // how many digits to use in file numbering
  char postFile[256];   // file (prefix) for post-processing
  int  fileMin, fileMax;
  FILE *monitor;        // pointer to monitor file or stderr

  /* premixed options */
  char premixTab[256];
  char premixTabP[256];

  // options for ensight
  char ensightCas[256]; // ensight case file
  char ensightDir[256]; // directory for Ensight file output

  // optional outputs for data/scalars file
  int outputViscSrc;    // output viscous heat production
  int outputGPDF;    // output viscous heat production
  int outputHtot;       // output total enthalpy
  int outputHR;         // output heat release in the scalars file 
  int outputDryNOX;     // output Dry NOX 6 percent (only if species present)
  int outputChi;        // output average scalar dissipation rate to scalars
  int outputRho;        // output average density to scalars
  int outputRhoSigma;        // output average density to scalars
  int outputMu;         // output average viscosity
  int outputPDF;        // output pdf (transformed) to flamelet file
  int outputCFL;        // dump the CFL of the grid to the file
  int outputCp;         // output specific heat
  int outputSumYi;
  int outputQdot;       // output heat losses
  int outputGamma;      // output ratio of specific heats
  int outputConv;       // output convection coefficient
  int outputSootSrc;    // output soot source term components
  int outputSootMxmyp;
  int outputConvVel;    // output convection velocities
  /* ----- display variables for monitoring solution ----- */
  int   nMonitorY, nMonitorZ, *monitorIndex;
  char  **monitorSpecies;

  /* ----- storage for reading from input file ----- */
  int   *nDim;
  int   *n1, *n2;
  int   nOx, nFuel1, nFuel2; // number of components in ox and fuel streams
  int   nChemFuel, chemFuelisY;
  char  **chemFuel;
  double *chemFuelY;
  char  bcType[3];        // specifies type of boundary condition input
  
  char   **inputSpecies;  // character array containing names of species to
                          // read in from input file
  int    *inputIndex;     // indices of inputSpecies from mechanism
  double **initY;         // initial profiles of species
  double *initT;          // initial temperature profile
  double *gPDF;          // initial temperature profile
  char   **outputSpecies; // character array containing names of species to 
                          // output for each flamelet
  int    *outputIndex;    // indices of outputSpecies from mechanism
  int    nOutputVar;
  char   **outputVar;    // additional variables to output
  int    initSootMom;
  double **initSoot;     

  /* ----- oxidizer properties ----- */
  char   **oxName;   // oxidizer species identifiers in stream 1 
  int    *oxIndex;   // oxidizer in stream 1 mechanism index  
  double *oxY;       // mass fractions of oxidizer species in stream 1
  double oxH;        // Enthalpy of oxidizer in stream 1 
  double oxT;        // temperature of oxidizer in stream 1
  double oxP;        // pressure of ambient
  double Mox;        // molecular weight of the oxidizer
  
  /* ----- fuel stream 1 properties ----- */
  char   **fuelName1; // fuel species identifiers in stream 1 
  int    *fuelIndex1; // fuel in stream 1 mechanism index 
  double *fuelY1;     // mass fractions of fuel species in stream 1 
  double fuelH1;      // Enthalpy of fuel in stream 1 
  double fuelT1;      // temperature of fuel in stream 1
  double fuelTref1;   // reference temperature for liquid fuel stream 1
  double hVap1;       // enthalpy of vaporisation for fuel 1
  double cpLiq1;      // liquid specific heat of fuel 1
  double Mfuel1;      // molecular weight of the fuel
  char   satData[256];// file containing saturation data for liquid fuel 1

  /* ----- fuel stream 2 properties ----- */
  char   **fuelName2; // fuel species identifiers in stream 2 
  int    *fuelIndex2; // fuel in stream 2 mechanism index 
  double *fuelY2;     // mass fractions of fuel species in stream 2 
  double fuelH2;      // Enthalpy of fuel in stream 2 
  double fuelT2;      // temperature of fuel in stream 2
  double fuelTref2;   // reference temperature for liquid fuel stream 2
  double hVap2;       // enthalpy of vaporisation for fuel 2
  double cpLiq2;      // liquid specific heat of fuel 2

  /* ----- statistics data ----- */
  double mean[3], var[4]; // mean and variance read in from header
  
  /* ----- miscellaneous input to store ----- */
  int    nMisc;       // number of undefined variables
  double **miscIn;    // storage for undefined variables
  char   **miscNames; // name for undefined variables

  /* ----- variable grid parameters ----- */
  // these may change in the future
  int     nPivot1, nPivot2;  // number of pivot points
  double *gPivot1, *gPivot2; // pivot points
  double *gPerc1, *gPerc2;   // percentage of grid points below pivot

  /* ----- stand alone calculation ----- */
  chiStruct *chiData; 
  int nFiles;
};
typedef struct IOSTRUCT ioStruct;

#endif /* MDUC_IO_STRUCT_H */
