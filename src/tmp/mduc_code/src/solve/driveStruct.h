#ifndef MDUC_DRIVER_STRUCT_H
#define MDUC_DRIVER_STRUCT_H

struct DRIVESTRUCT
{
  /* !!!!!!! potentially not necessary ???????? */
  int nFlamelets;       // number of flamelets to solve
  int maxFlamelets;     // limit maximum number of flamelets allowed

  int maxNum1, maxNum2; // Max Grid points in Z1 and Z2 directions
  
  int FLNO;        // This is the FLAMELETNUMBER, which is a place holder for
                   // the array of flamelet structures in FLAME_STRUC
  
  /* ----- solver flags ----- */
  // TODO: !!!!! These need to be cleaned up
  int PHYS_COORD;  // 0: mixture fraction space
                   // 1: physical space
  int ENGINE;
  int SDIAL;       // 0: use CVODE as equation solver
                   // 1: use IDA as equation solver
  int TIME_MODE;   // 0: time variables based on CFD solver
                   // 1: base dt_final on t_end and nIter (from input file)
  int PRES_MODE;   // pressure treatment: 
                   //  P: isobaric
                   //  V: isochoric 
                   //  M: IC-engine-style with V(t)
  int MIX_MODE;    // flamelet equation transformation
                   //  1: jet ignition ?????
                   //  2: Hasse Formulation (x1=Z1+Z2, x2=Z2/x1) 
                   //  3: x1 = enthalpy, x2 = mixture fraction
  int CHI_MODE;    // scalar dissipation rate:
                   //  1: Chi(Z) based on ChiRef[0] and ZRef[0],
                   //  2: Chi(Z) given by DRIVER
  int CHI_INTERP;  // scalar dissipation rate interpolation:
                   //  0: None
                   //  1: Linear
  int CHI_SCALE;   // model scalar dissipation rate scaling
                   //  0: None
                   //  1: scale based on Zi_max
  int PDF_MODE;    // pdf treatment for post-processing
  
  int ENERGY;      // 1: solve enthalpy equation
                   // 0: solve temperature equation
  int VISC_MODE;   // 1: calculate mixture viscosity implicitly
                   // 0: calculate mixture viscosity explicitly
  int Q_LOSS;      // 0: adiabatic
  int H_DIM;       // 0: Enthalpy not a dimension
                   // 1: x1 is enthalpy dimension
                   // 2: x2 is enthalpy dimension
  int TRANS_Z2;    // 0: standard transformation
                   // 1: Hasse transformation
  int PHI_HAT;     // 0: filter the value
                   // 1: filter the derivative
  int DIM;         // flag for current direction in 2D ADI solver

  double ABSTOL;       // tolerance for convergence criteria
  double TABSTOL;      // tolerance for convergence criteria of temperature
  double SABSTOL;      // tolerance for convergence criteria of soot moments
  double QABSTOL;      // tolerance for convergence criteria for quadrature 
  double RELTOL;       // relative tolerance
  double QRELTOL;       // relative tolerance for quadrature
  double DIFF_MODE[3]; // 1: Diffusion term is solved for specified direction
                       // 0: chi is set to 0 & diffusion term is not solved
  double CONV_MODE;    // 1: Convection term is solved
                       // 0: chi is set to 0 & convection term is not solved

  int output_tsrc; // output temperature source terms
  int output_w;    // output reaction rates
  int sootMode;
  int thermophoresis; // thermophoresis enabled?
  int sootDiff; // soot diffusion enabled?
  double sootDiffusivity; // soot diffusion enabled?
  int detailedDiff; // Detailed diffusion? 
  int radiation; // thermophoresis enabled?
  int ignite; // thermophoresis enabled?
  int readBody;    // whether to read body of flamelet/flow file
  int DNS;         // Check if DNS is used

  /* ----- grid variables ----- */
  double gridPivot;    // point about which to centre finest grid
  double gridMaxDelta; // largest permissble grid spacing
};
typedef struct DRIVESTRUCT driveStruct;

#endif /* MDUC_DRIVER_STRUCT_H */

