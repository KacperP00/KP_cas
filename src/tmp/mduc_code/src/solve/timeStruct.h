#ifndef MDUC_TIME_STRUCT_H
#define MDUC_TIME_STRUCT_H

typedef struct timeStruct_ {
  /* ----- iteration parameters ----- */
  int    nt;         /* current time step number */
  int    nt0;        /* initial time step number */
  int    nIter;      /* total number of time steps (if known) */

  double t;          /* current time */
  double t_n;        /* time at end of step */
  double t_m1;       /* previous time */
  double t0;         /* initial time */
  double dt0;        /* initial time-step size */
  double t_end;      /* time at end of simulation */
  double dt;         /* current time-step size */
  double dt_pre;     /* previous time-step size */
  double dt_diff;    /* Diffusion Time-Scale */
  double dt_cvode_pre; /* previous cvode time step */
  double dt_max_;    /* maximum time step (local)*/
  double dt_max;     /* maximum time step (global)*/
  
  /* ----- for numerics ----- */
  double CFL;        /* For Explicit Time-step size */
  double tN1;        /* time at n+1 for use in ADI scheme */
  double tstar;      /* intermediate time */

  /* ----- for engines ----- */
  double RPM;        /* engine speed */
  double CAD_0;      /* starting crank angle */

  /* ----- define delays ----- */
  double start2D;    /* time when 2D solver begins, starting from 1D flamelet */
  double delay2D;    /* delay after 2D solver start before evaluating Z2max */
} timeStruct;
/* typedef struct TIMESTRUCT timeStruct; */

#endif /* MDUC_TIME_STRUCT_H */
