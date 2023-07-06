#ifndef MDUC_SDIAL_STRUCT_H
#define MDUC_SDIAL_STRUCT_H

#ifdef SUNDIALS_FROM_4_0
#define SUNDIALS_FROM_3_2
#endif

#ifdef SUNDIALS_FROM_3_2
#define SUNDIALS_FROM_3_0
#endif

/* ----- sundials library headers ----- */
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_diag.h>
#include <ida/ida.h>
#include <cvodes/cvodes_bandpre.h>
#include <cvodes/cvodes_bbdpre.h>
#include <kinsol/kinsol.h>
#include <kinsol/kinsol_bbdpre.h>

#ifdef SUNDIALS_FROM_3_0

#include <cvodes/cvodes_direct.h>
#include <cvodes/cvodes_spils.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>

#else /* SUNDIALS lower than 3.0 */

#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <ida/ida_dense.h>
#include <ida/ida_band.h>
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_band.h>
#include <kinsol/kinsol_spbcgs.h>
#include <kinsol/kinsol_spgmr.h>
#include <kinsol/kinsol_sptfqmr.h>

#endif /* SUNDIALS_FROM_3_0 */

#include <nvector/nvector_serial.h>
#ifdef ACCESS_MPI
#include <nvector/nvector_parallel.h>
#endif /* ACCESS_MPI */

struct SDIALSTRUCT {
  /* ----- solver settings ----- */
  char solver[16];   /* ODE vs. DAE system */
  char method[16];   /* linear algebra method */
  char precnd[16];   /* preconditioner type */
  int  iSolve;       /* 0: use CVODE, 1: use IDA, 2: use KINSOL */
  int  iMethod;      /* 0: dense, 1: banded, 2: gmres, 3: bicgstab, 4: spqmr */
  int  iPreCnd;      /* 1: left, 2: right, 3: both */
  int  ADI;          /* 0: use implicit, 1: use ADI for 2D solver */

  /* ----- settings ----- */
#ifdef SUNDIALS_FROM_3_0
  sunindextype nEq;  /* number of equations to be solved */
  sunindextype nEqL; /* number of equations to be solved */
  sunindextype nQEq;  /* number of equations to be solved */
  sunindextype nQEqL; /* number of equations to be solved */
  sunindextype sunIndTypeTmp; /* work variable for sundials 3 interface */
#else
  int nQEq;           /* number of equations to be solved */
  int nQEqL;          /* local equations for parallel runs */
  int nEq;           /* number of equations to be solved */
  int nEqL;          /* local equations for parallel runs */
#endif /* SUNDIALS_FROM_3_0 */

  double abstol;     /* absolute tolerance */
  double reltol;     /* relative tolerance */

  /* ----- sundial memory ----- */
#ifdef SUNDIALS_FROM_3_0
  SUNMatrix A;
  SUNLinearSolver LS;
#endif /* SUNDIALS_FROM_3_0 */

  void     *mem;            /* workspace memory */
  
  N_Vector Y;               /* solution vector space */
  N_Vector Yp;              /* solution derivative vector */
  N_Vector Yscale, Fscale;  /* scaling arrays for non-linear Newton solve */
  N_Vector Q;               /* quadrature vector space */
  N_Vector tmp;             /* temporary storage */
  int      itask;           /* flag indicating the job of the solver */
  int      iret;            /* return flag of the solver */
  long int totalIt;         /* total iterations (including reinits) */

  /* assign pointers to access memory in N_Vector objects */
  realtype *Ydata, *YpData, *resData, *YscaleData, *FscaleData;
  realtype *Qdata;

  /* variables for ADI scheme */
  int dir, row, col; /* implicit direction and current row/column */

  /* array of pointers to pass data to solver functions */
  void **args;
};
typedef struct SDIALSTRUCT sdialStruct;

/* define pointer functions for RHS functions */
int (*cvodeRHS)(realtype, N_Vector, N_Vector, void *);
#ifdef SUNDIALS_FROM_2_7
int (*cvodeRHS_L)(long int, realtype, N_Vector, N_Vector, void *);
#else 
#ifdef SUNDIALS_FROM_3_0
int (*cvodeRHS_L)(sunindextype, realtype, N_Vector, N_Vector, void *);
#else /* SUNDIALS LOWER THAN 2.7 */
int (*cvodeRHS_L)(int, realtype, N_Vector, N_Vector, void *);
#endif /* SUNDIALS_FROM_3_0 */
#endif /* SUNDIALS_FROM_2_7 */
int (*idaRHS)(realtype t, N_Vector y, N_Vector yp, N_Vector r, void *data);

/* define general functions for functions that are similar between solvers */
int  (*SDialSetStopTime)(void *, double);
int  (*SDialGetNumSteps)(void *, long int *);
int  (*SDialGetNumRhsEvals)(void *, long int *);
int  (*SDialGetNumLinSolvSetups)(void *, long int *);
int  (*SDialSetInitStep)(void *, double);
int  (*SDialGetCurrentStep)(void *, realtype *);
int  (*SDialGetCurrentTime)(void *, realtype *);
int  (*SDialGetLastStep)(void *, double *);
void (*SDialFree)(void **);
int  (*SDialGetNumJacEvals)(void *, long int *);

#endif /* MDUC_SDIAL_STRUCT_H */
