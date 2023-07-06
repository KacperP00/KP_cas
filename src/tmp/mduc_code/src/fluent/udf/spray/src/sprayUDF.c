/* ============================================================================
                   UDF for spray injection

 This UDF takes the velocity and/or mass flow rate profile for injection.
 If only velocity profile is provided, then the mass flow rate profile is
 assumed to be Cd*rho*velocity*Anozzle. Because of the above assumption,
 the amount of fuel injected, which is critical, may not match the value that
 the user wants. In order to match the amount of fuel injected, one has to tune
 the Cd above.  The udf provides an easy way to calculate Cd. Execute the
 profile_integral udf and follow the instructions and you will be able to
 calculate the Cd easily. If mass flow rate profile is provided, then there is
 no need to calculate Cd since the amount of fuel injected is already well
 defined by the mass flow rate profile.

    How to use the udf:

    - Determine what type of profiles will be used. The udf supports
      three types:
        1) velocity profile only, and the mass flow rate is assumed to be
           related to the velocity profile
        2) mass flow rate profile only
        3) both velocity and mass flow rate profiles.
    - Create an injection profile that has the following format.
        Header lines prefixed by #
        First column is time from start of injection. 
        Second column is velocity magnitude or mass flow rate.

         A sample of velocity profile:
           # time[s] velocity[m/s]
           2   0
           4   300
           15  300
           17  0

         A sample of mass flow rate profile:
           # time[s] mdot[kg/s]
           2   0
           4   0.2e-3
           10  0.3e-3
           15  0.4e-3
           17  0

    - Modify the discharge coefficient in the udf. 
      With velocity profile only, the mass flow rate is assumed to be
        Cd*rho*vel*Anozzle. 
      meaning that Cd is the parameter to tune to make sure that the amount of
      fuel injected is correct.  Cd can be calcuated by the profile_integral
      on_demand udf.
    - Set up dpm.  (Use only 1 dpm iteration per flow time step if breakup or
      collision is used.)
    - Set up injection.
    - Set up some positive velocity magnitude, which will then be changed
      by the udf if velocity profile is used. 
      YOU CAN NOT USE ZERO VELOCITY MAGNITUDE HERE.
    - Hook up the udf to dpm initialization

  Written by : Xiao HU (Fluent Inc.)
  Last updated : 12/22/2006

        changed to time instead of Crank angle for common rail systems
        last updated : Godehard Nentwig 01/17/2007

        changed discharge coefficient variable (DC) to Cd to avoid conflict in
        Fluent V12.0.15 : Eric Doran 04/30/2009

        Heavily modified to extend capability to arbitrary droplet
        distributions and compute various outputs.  Is now heavily dependent
        on appropriate Scheme variables: Eric Doran 05/31/2011

============================================================================ */

#include "sprayUDF.h"

/* global library variables */
static int  nInj, *nMdot, *nVel;
static int  DPM_writeHeader = 1;
static real **t_mdot, **mdot, **t_vel, **vel, *dropSize, *dropPDF, sumPD3;
static char **injID;
static FILE *sprayFile;

DEFINE_ON_DEMAND(loadProfile)
{
#if !RP_HOST
  int  i, j;
  char **massFile, **velFile, *dropFile;
  real   Cd, D_nozz, A_nozz, rho_f, sumPDF;
  double PI = 3.1415926535897931;
  char buffer[512];

  /* new spray profiles, make sure stats file will write header */
  DPM_writeHeader = 1;

  /* get the injection ids */
  int len, cStart, cEnd;
  char noz1[64];
  nInj = UDFsprayParseFluentString(RP_Get_String("spray/injection-ids"),&injID);
  if (nInj == 0)
    Error("UDF Spray: injection must be specified to use spray library!\n");
  /* injection names will likely be read as a list */
  for (j=0; j<nInj; j++)
    {
      len    = strlen(injID[j]);
      cStart = injID[j][0];
      cEnd   = injID[j][len-1];

      /* removing leading or trailing brackets*/
      if (cStart == '(')
        {
          strcpy(buffer, &injID[j][1]);
          strcpy(injID[j], buffer);
          len    = strlen(injID[j]);
          cEnd   = injID[j][len-1];
        }
      if (cEnd == ')')
        {
          strncpy(buffer, injID[j], len-1);
          buffer[len-1] = '\0';
          strcpy(injID[j], buffer);
        }
    }

  nMdot  = (int*) calloc(nInj, sizeof(int));
  nVel   = (int*) calloc(nInj, sizeof(int));
  t_mdot = (double**) calloc(nInj, sizeof(double*));
  mdot   = (double**) calloc(nInj, sizeof(double*));
  t_vel  = (double**) calloc(nInj, sizeof(double*));
  vel    = (double**) calloc(nInj, sizeof(double*));

  Message0("\n----- Loading spray profiles for %i injections -----\n", nInj);
  int nFile;
  nFile = UDFsprayParseFluentString(RP_Get_String("spray/injection-massflow"),
                                    &massFile);
  if (nFile != nInj)
    Message0("Warning: % injections and only %i mass flow profiles\n",
             nFile, nInj);
  nFile = UDFsprayParseFluentString(RP_Get_String("spray/injection-velocity"),
                                    &velFile);
  if (nFile != nInj)
    Message0("Warning: % injections and only %i velocity profiles\n",
             nFile, nInj);

  /* read in the profiles for each injection */
  for (j=0; j<nInj; j++)
    {
      /* prepend the injection data directory*/
      UDFsprayPrefixString(RP_Get_String("spray/injection-data"), massFile[j]);
      UDFsprayPrefixString(RP_Get_String("spray/injection-data"), velFile[j]);

      Message0("  Profiles for injection %i (%s)\n", nInj, injID[j]);
      Message0("  Mass flow profile: %s\n", massFile[j]);
      Message0("  Velocity profile:  %s\n", velFile[j]);

      if (!strncmp(velFile, "compute", 5))
        { /* specify mass flow rate only, compute velocity based on 
             u = mdot/(Cd*rho*A) 
          */
          /* read in mass flow rate profile */
          nMdot[j] = DPM_read_profile(massFile[j], &t_mdot[j], &mdot[j]);
          
          /* allocate velocity */
          nVel[j] = nMdot[j];
          vel[j]  = (real*) calloc(nVel[j], sizeof(real));
          
          /* retrieve nozzle variables */
          D_nozz = RP_Get_Real("spray/nozzle-diameter");
          Cd     = RP_Get_Real("spray/discharge-coeff");
          rho_f  = RP_Get_Real("spray/fuel-density");
          A_nozz = 0.25*PI*D_nozz*D_nozz;
          
          /* compute the velocity */
          for (i=0; i<nVel[j]; i++)
            vel[j][i] = mdot[j][i]/(Cd*rho_f*A_nozz);
        }
      else
        {
          /* read in both mass flow and velocity profiles */
          nMdot[j] = DPM_read_profile(massFile[j], &t_mdot[j], &mdot[j]);
          nVel[j]  = DPM_read_profile(velFile[j], &t_vel[j], &vel[j]);
          
          if (nMdot[j] != nVel[j])
            {
              Message0("DPM: Warning: mass flow and velocity profiles ");
              Message0("have different number of data points\n");
            }
        }
    } /* end loop over injections */
  Message0("Finished reading profiles\n");

  /* ----- initialize the droplet distribution ----- */
  /* NOTE: This will assume the same distribution for all injections/nozzles.
           Could potentially generalize to an individual injection basis, for
           now take the information from the first injection, first nozzle 
           (one nozzle is assumed to be present and numbered 'id-1' */
  dropFile = RP_Get_String("spray/injection-dropsize");
  UDFsprayPrefixString(RP_Get_String("spray/injection-data"), dropFile);

  sprintf(noz1, "%s-1", injID[0]);
  Injection *I = Pick_Injection(noz1);
  if (!I)
    Message0("UDF Spray: injection %s not found\n", noz1);

  Particle  *p;
  int    iDiam;
  int    nDrop, nDiam, nTot;
  double tmp, tmp2, delD, *tmp_dsize, *tmp_dpdf;
  double Dq0, q, alpha; /* parameters for Dumouchel */

  /* set the droplet distribution */
  nDiam = I->number_diameters;
  delD  = (I->max_diam - I->min_diam)/(nDiam-1);
  dropSize = (double*) calloc(nDiam, sizeof(double));
  dropPDF  = (double*) calloc(nDiam, sizeof(double));

  if (!strcmp(strrchr(dropFile, '/')+1, "dumouchel"))
    { /* use Dumouchel distribution */
      Message0("  Computing Dumouchel distribution\n");

      /* specify droplet sizes */
      dropSize[0] = I->min_diam;
      for (i=1; i<nDiam; i++)
        dropSize[i] = dropSize[i-1] + delD;

      /* compute distribution */
      alpha = RP_Get_Real("spray/dumouchel-alpha");
      q     = I->spread_parameter;
      Dq0   = I->mean_diam;

      for (i=0; i<nDiam; i++)
        {
          tmp = log(q) - UDF_lnGamma(alpha/q) + (alpha/q)*log(alpha/q)
            + (alpha-1)*log(dropSize[i]) - alpha*log(Dq0)
            - alpha/q*pow(dropSize[i]/Dq0,q);
          dropPDF[i] = exp(tmp);
        }
    }
  else if (strncmp(strrchr(dropFile, '/')+1, "fluent", 3))
    { 
      Message0("  Reading distribution from: %s\n", dropFile);
      /* read in from profile */
      nDrop = DPM_read_profile(dropFile, &tmp_dsize, &tmp_dpdf);

      /* normalize the distribution */
      sumPDF = 0.0;
      for (i=0; i<nDrop; i++) sumPDF     += tmp_dpdf[i];
      for (i=0; i<nDrop; i++) tmp_dpdf[i] /= sumPDF;

      /* either assign or interpolate to fluent diameter grid */
      if (nDrop == nDiam)
        {
          Message0("assigning directly\n");
          for (i=0; i<nDiam; i++)
            {
              dropSize[i] = tmp_dsize[i];
              dropPDF[i]  = tmp_dpdf[i];
            }
        }
      else
        {
          Message0("Interpolating droplet distribution to %i bins\n", nDiam);
          delD = (tmp_dsize[nDrop-1]-tmp_dsize[0])/nDiam;
          dropSize[0] = tmp_dsize[0];
          for (i=1; i<nDiam; i++)
            dropSize[i] = dropSize[i-1] + delD;

          /* interpolate the size distribution */
          for (i=0; i<nDiam; i++)
            {
              j = 1;
              while (dropSize[i] > tmp_dsize[j]) { j++; }
              dropPDF[i] = tmp_dpdf[j-1] + (tmp_dpdf[j]-tmp_dpdf[j-1])/
                (tmp_dsize[j]-tmp_dsize[j-1])*(dropSize[i] - tmp_dsize[j-1]);
            }
        }

      /* ensure the interpolated distribution is still normalized correctly */
      sumPDF = 0.0;
      for (i=0; i<nDiam; i++) sumPDF     += dropPDF[i];
      for (i=0; i<nDiam; i++) dropPDF[i] /= sumPDF;

      /* unscale the pdf */
      dropPDF[0]       = dropPDF[0]/(dropSize[1]-dropSize[0]);
      dropPDF[nDiam-1] = dropPDF[nDiam-1]/(dropSize[nDiam-1]-dropSize[nDiam-2]);
      for (i=1; i<nDiam-1; i++)
        dropPDF[i] = dropPDF[i]*2/(dropSize[i+1]-dropSize[i-1]);

      free(tmp_dsize); free(tmp_dpdf);
    }
  else /* use Fluent distribution */
    {
      Message0("  Using fluent distribution with %i bins\n", nDiam);

      /* specify droplet sizes */
      dropSize[0] = I->min_diam;
      for (i=1; i<nDiam; i++)
        dropSize[i] = dropSize[i-1] + delD;
      
      /* compute distribution */
      q     = I->spread_parameter;
      Dq0   = I->mean_diam;

      for (i=0; i<nDiam; i++)
        {
          tmp2 = dropSize[i]/Dq0;
          tmp = log(q/Dq0) + (q-1)*log(tmp2) - pow(tmp2,q);
          dropPDF[i] = exp(tmp);
        }
    }
  
  /* compute and store the term for computing total number */
  for (i=0, sumPD3=0.0; i<nDiam; i++)
    sumPD3 += dropPDF[i]*dropSize[i]*dropSize[i]*dropSize[i];

  /* integrate the distribution */
  sumPDF = dropPDF[0]*(dropSize[1]-dropSize[0])
    + dropPDF[nDiam-1]*(dropSize[nDiam-1]-dropSize[nDiam-2]);
  for (i=1; i<nDiam-1; i++)
    sumPDF += dropPDF[i]*0.5*(dropSize[i+1]-dropSize[i-1]);
  
  Message0(" ====== Droplet Distribution ======\n");
  Message0("   # diam.[m] pdf[-]\n");
  for (i=0; i<nDiam; i++)
    Message0("%15.6e %15.6e\n", dropSize[i], dropPDF[i]);
  Message0("\n15s %15.6e\n", "sum:", sumPDF);
  Message0(" ==================================\n");

  /* ----- initialize the stats file ----- */
  int  nT, cur_n;
  char fileLine[2056], dummy[2056], *prev;
  long nPrev;

  Injection *I_list = Get_dpm_injections();
  int  nNozzle;
  char **I_ids;

  nNozzle = 0;
  loop(I, I_list) { nNozzle++; }

  /* allocate memory */
  I_ids = (char**) calloc(nNozzle, sizeof(char*));
  for (i=0; i<nNozzle; i++)
    I_ids[i] = (char*) calloc(512, sizeof(char));

  i = 0;
  loop(I, I_list) {
    strcpy(I_ids[i], I->name);
    i++;
  }

#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      /* check to see if the file exists */
      sprayFile = NULL;
      sprayFile = fopen(RP_Get_String("spray/monitor-file"), "r");
      if (sprayFile) {
        Message0("  monitor %s exists\n", RP_Get_String("spray/monitor-file"));
        
        /* read in the file and search for the time step */
        cur_n = -1;
        while (!feof(sprayFile)) {
          nPrev = ftell(sprayFile);
          fgets(fileLine, 2056, sprayFile);
          
          if (fileLine[0] != '#' && fileLine[0] != '\n')
            sscanf(fileLine, "%i %*s", &cur_n);
          
          if (cur_n >= N_TIME)
            break;
        } /* end loop over file */
        
        /* go to beginning of file */
        fseek (sprayFile, 0, SEEK_SET);
        
        /* read in everything up until the current time step */
        prev = (char*) calloc(nPrev, sizeof(char));
        fread(prev,1,nPrev,sprayFile);
        
        /* reopen the file and write the previous information */
        freopen(RP_Get_String("spray/monitor-file"), "w", sprayFile);
        fwrite(prev,1,nPrev,sprayFile);
        
        free(prev);
      }
      else { /* create the file new and write the header */
        sprayFile = fopen(RP_Get_String("spray/monitor-file"), "w");
        
        j = 1;
        fprintf(sprayFile, "# nt<%i> ", j);    j++;
        fprintf(sprayFile, "time[s]<%i> ", j); j++;
        if (in_cylinder_p) {
          fprintf(sprayFile, "CAD[deg]<%i> ", j);
          j++;
        }
        fprintf(sprayFile, "total_mass[kg]<%i> ", j); j++;
        for (i=0; i<nInj; i++, j++)
          fprintf(sprayFile, "mass_%s[kg]<%i> ", injID[i], j);

        for (i=0; i<nNozzle; i++, j++)
          fprintf(sprayFile, "pen_01_%s[m]<%i> ", I_ids[i], j);

        for (i=0; i<nNozzle; i++, j++)
          fprintf(sprayFile, "pen_05_%s[m]<%i> ", I_ids[i], j);

        for (i=0; i<nNozzle; i++, j++)
          fprintf(sprayFile, "pen_max_%s[m]<%i> ", I_ids[i], j);
        
        for (i=0; i<nNozzle; i++, j++)
          fprintf(sprayFile, "SMD_%s[m]<%i> ", I_ids[i], j);

        for (i=0; i<nNozzle; i++, j++)
          fprintf(sprayFile, "film_%s[kg]<%i> ", I_ids[i], j);
        
        fprintf(sprayFile, "\n");
      }
      fflush(sprayFile);
    }

  for (i=0; i<nNozzle; i++) free(I_ids[i]);
  free(I_ids);

#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

DEFINE_DPM_INJECTION_INIT(initDPM, I)
{
#if !RP_HOST
  Particle *p;

  int  i, j, nTot, iInj, len;
  real t_start, time = RP_Get_Real("flow-time");
  real CAD;
  real vel_n, mdot_n, mdot_part, axis[ND_ND];
  real t_velm, t_velp, t_mdotm, t_mdotp;
  real tmp;
  double pi_over_6 = M_PI/6.0;
  char buffer[DPM_NAME_LENGTH];

  /* figure out which injection profile this uses */
  /* NOTE: assumes injection is <id>-[1-9].  Can't handle more than 9 holes */
  iInj = -1;
  for (j=0; j<nInj; j++)
    {
      len = strlen(I->name);
      strncpy(buffer, I->name, len-2);
      buffer[len-2] = '\0';
      if (!strcmp(buffer, injID[j]))
        iInj = j;
    }
  if (iInj < 0)
    Error("spray: couldn't find injection id %s\n", I->name);

  Message0("***** DPM: Initializing particles for: %s\n", I->name);
  /* Message0("start time: %15.6e\n", I->start_time); */
  /* Message0("unsteady start time: %15.6e\n", I->unsteady_start); */
  /* Message0("unsteady start CA:   %15.6e\n", I->unsteady_ca_start); */

  if (in_cylinder_p) {
    t_start = (I->unsteady_ca_start - in_cyn_ca_start)/(6*in_cyn_ca_rpm);
    /* account for periodicity */
    CAD = RP_Get_Real("dynamesh/in-cyn/crank-angle");
    t_start += floor(CAD/(in_cyn_ca_start+720.0))*120/in_cyn_ca_rpm;
  }
  else
    t_start = I->unsteady_start;

  /* return if current time is outside of specified profile
     mass flow and velocity should have same start/end times */
  if (time < (t_mdot[iInj][0]             + t_start) ||
      time > (t_mdot[iInj][nMdot[iInj]-1] + t_start))
    {
      loop(p, I->p_init) { p->stream_index = -1; }
      return;
    }
  
  /* find the time in profile just past current time */
  i = j = 1;
  while (((time - t_start) > t_vel[iInj][i]) && i < nVel[iInj]-1 )  { i++; }
  while (((time - t_start) > t_mdot[iInj][j]) && j < nMdot[iInj]-1) { j++; }

  /* store time either side of current time, adjusted to flow time */
  t_velm  = t_vel[iInj][i-1]  + t_start;  
  t_velp  = t_vel[iInj][i]    + t_start;
  t_mdotm = t_mdot[iInj][j-1] + t_start;
  t_mdotp = t_mdot[iInj][j]   + t_start;

  /* linearly interpolate quantities to current time */
  vel_n  = vel[iInj][i-1] 
    + (vel[iInj][i]-vel[iInj][i-1])/(t_velp-t_velm)*(time-t_velm);
  mdot_n = mdot[iInj][j-1]
    + (mdot[iInj][j]-mdot[iInj][j-1])/(t_mdotp-t_mdotm)*(time-t_mdotm);

  /* set Fluent mass flow to the profile */
  I->total_flow_rate = mdot_n;

  /* take the droplet diameters to be equispaced */
  double deltaDiam = (I->max_diam-I->min_diam)/I->number_diameters;

  if (dropSize == NULL)
    {
      Message0("Initializing particle distribution from fluent\n");
      int iDiam;
      double sumPDF, *mdot_pdf, mdot_stream, *rr_dist;
      Message0("Using fluent distribution with %i bins\n", I->number_diameters);

      dropSize = (double*) calloc(I->number_diameters, sizeof(double));
      dropPDF  = (double*) calloc(I->number_diameters, sizeof(double));
      mdot_pdf = (double*) calloc(I->number_diameters, sizeof(double));
      rr_dist  = (double*) calloc(I->number_diameters, sizeof(double));

      iDiam = nTot = 0;
      mdot_stream = 0.0;
      loop(p, I->p_init)
        {
          dropSize[iDiam] = P_DIAM(p);
          dropPDF[iDiam]  = p->number_in_parcel;
          mdot_pdf[iDiam] = p->number_in_parcel*P_MASS(p); /* p->flow_rate; */
          nTot += p->number_in_parcel;
          mdot_stream += p->number_in_parcel*P_MASS(p);/*  p->flow_rate; */
          iDiam++;
          if (iDiam == I->number_diameters)
            break;
        }

      Message0("ntotal init: %i\n", nTot);
      sumPDF = 0.0;
      for (i=0; i<I->number_diameters; i++)
        {
          dropPDF[i] /= nTot;
          
          sumPDF += dropPDF[i];
        }
      Message0("sumPDF = %15.6e\n", sumPDF);
      for (i=0; i<I->number_diameters; i++) dropPDF[i] /= sumPDF;

      /* unscale the distribution */
      for (i=0; i<I->number_diameters; i++) dropPDF[i] /= 0.5*deltaDiam;

      /* write the distribution to file to check */
      FILE *dFile = fopen("fluent_particle_dist.dat", "w");
      
      fprintf(dFile, "# diam P\n");
      for (i=0; i<I->number_diameters; i++)
        fprintf(dFile, "%15.6e %15.6e %15.6e %15.6e\n", dropSize[i], dropPDF[i],
                mdot_pdf[i], mdot_pdf[i]/mdot_stream);

      fclose(dFile);

      for (i=0, sumPD3 = 0.0; i<I->number_diameters; i++)
        sumPD3 += dropPDF[i]*dropSize[i]*dropSize[i]*dropSize[i];
      free(rr_dist); free(mdot_pdf);
    }

  double rhof, D3;
  /* TODO: figure out better way to get fuel density from injection */
  loop(p, I->p_init)
    {
      rhof = P_RHO(p);
      break;
    }

  /* if droplet distribution specified, set particle diameters */
  int iDiam;
  double mdot_stream = I->total_flow_rate/I->n_pts;

  /* compute the total number of particles required */
  nTot = (int)(floor(mdot_stream/(pi_over_6*rhof*sumPD3)));

  /* construct the distribution */
  iDiam = 0;
  loop(p, I->p_init)
    {
      P_DIAM(p) = dropSize[iDiam];
      D3        = P_DIAM(p)*P_DIAM(p)*P_DIAM(p);
      P_MASS(p) = pi_over_6*P_RHO(p)*D3;
      iDiam++;
      if (iDiam == I->number_diameters)
        iDiam = 0;
    }

  /* p: struct for particles. I: struct for injection */
  /* loop over all particles */
  iDiam = 0;
  loop(p, I->p_init)
    {
      /* Update velocity */
      tmp = NV_MAG(p->state.V);
      NV_VS(axis, =, p->state.V, /, tmp);
      NV_VS(p->state.V, =, axis, * ,vel_n);

      /* set the flow rate based on the particle number distribution */
      p->flow_rate = nTot*dropPDF[iDiam]*P_MASS(p);
          
      /* determine number of particles in parcel */
      p->number_in_parcel = p->flow_rate*CURRENT_TIMESTEP/P_MASS(p);

      iDiam++;
      if (iDiam == I->number_diameters)
        iDiam = 0;
    }
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

int DPM_read_profile(char *fileName, double **t, double **var)
{
  int i, len, nHead;
  long double tmp;
  char comment[1000];
 
  FILE *pFile = fopen(fileName, "r");

  if (!pFile)
    Error("Spray profile %s not found\n", fileName);

  /* get the number of header lines */
  nHead = -1;
  do {
    nHead++;
    strcpy(comment, "\n");
    fgets(comment, 1000, pFile);
  } while(comment[0] == '#');
  rewind(pFile);

  /* get the number of data lines */
  len = 0;
  while (!feof(pFile))
    {
      strcpy(comment, "\n");
      fgets(comment, 1000, pFile);
      if (strcmp(comment, "\n"))
        len++;
    }
  len -= nHead;
  rewind(pFile);

  Message0("Found %i data lines in %s\n", len, fileName);

  /* get the header line */
  for (i=0; i<nHead; i++) fgets(comment,1000,pFile);

  /* allocate profile memory */
  t[0]   = (real*) calloc(len, sizeof(real));
  var[0] = (real*) calloc(len, sizeof(real));

  for (i=0; i<len; i++)
    {
      /* read in time */
      fscanf(pFile, "%Lf", &tmp);
      t[0][i] = (real)(tmp);
      /* read in variable */
      fscanf(pFile, "%Lf", &tmp);
      var[0][i] = (real)(tmp);
    }

  fclose(pFile);
  return len;
}
/* ------------------------------------------------------------------------- */

/* Calculate the integral of a profile */
real DPM_integrate_profile(int n, real *t, real *data)
{
  int i;
  real integral;

  integral = 0;
  for(i=0; i<n-1; i++)
    integral += 0.5*(data[i+1]+data[i])*(t[i+1]-t[i]);

  return integral;
}
/* ------------------------------------------------------------------------- */

void DPM_print_profile(char *name, int n, real *t, real *data)
{
 int i;

 Message0("%d Points\n", n);
 Message0("# time %s\n", name);
 for(i=0; i<n; i++)
   Message0("%15.6e %15.6e\n", t[i], data[i]);
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(DPM_integrateProfile)
{
  real integral_vel, integral_massflow;

  integral_vel      = DPM_integrate_profile(nVel, t_vel, vel);
  integral_massflow = DPM_integrate_profile(nMdot, t_mdot, mdot);

  Message0("\n\nThe number of points in the vel curve is : %d", nVel);
  Message0("\nThe integral of the velocity curve is : %-12.4e\n",
           integral_vel);

  Message0("\n\nThe number of points in the massflow curve is : %d", nMdot);
  Message0("\nThe integral of the mass flow rate curve is : %-12.4e\n",
           integral_massflow);

}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(DPM_printProfile)
{
  int j;
  for (j=0; j<nInj; j++)
    {
      Message0("Profile for injection %i --> %s\n", j, injID[j]);
      if (fabs(t_mdot[j][0]-t_vel[j][0]) > 1e-5 ||
          fabs(t_mdot[j][nMdot[j]-1]-t_vel[j][nVel[j]-1]) > 1e-5)
        {
          Message0("\n\nThe two profiles do not have the same starting and ending time!!\n");
          Message0("starting time: (%-10.3e %10.3e)   ending time: (%-10.3e %10.3e)\n", t_mdot[j][0], t_vel[j][0], t_mdot[j][nMdot[j]-1], t_vel[j][nVel[j]-1]);
        }
      
      DPM_print_profile("Velocity", nVel[j], t_vel[j], vel[j]);
      DPM_print_profile("Mass flow", nMdot[j], t_mdot[j], mdot[j]);
    }
}
/* ------------------------------------------------------------------------- */

real DPM_penetration_length(Injection *I, real *pen_05, real *pen_max)
{
#if !RP_HOST
  Particle *p;

  int  i, j;
  int  nR;                      /* number of radial grid cells */
  real dr;                      /* radial grid cell size */
  real r;                       /* radial distance to current particle */
  real penetration = 0.0;       /* penetration length based on 1% mass [m] */
  real max_pen;                 /* maximum penetraion of single particle [m] */
  real one_sixth   = 1.0/6.0;   /* constant */
  real D_cubed;                 /* cube of particle diameter [m^3] */
  real part_mass;               /* mass of particles in parcel [kg] */
  real total_mass;              /* mass of all particles       [kg] */
  real accum_mass;              /* accumulated mass in radial grid [kg] */

  /* get penetration length based on spherical radius */
  total_mass = 0.0;
  max_pen    = 0.0;
  *pen_05    = 0.0;
  loop(p, I->p) {
    if (!p->on_wall) {
      total_mass += P_MASS(p)*P_N(p);
      r = sqrt(  (P_POS(p)[0] - I->pos1[0])*(P_POS(p)[0] - I->pos1[0])
               + (P_POS(p)[1] - I->pos1[1])*(P_POS(p)[1] - I->pos1[1])
               + (P_POS(p)[2] - I->pos1[2])*(P_POS(p)[2] - I->pos1[2]));
      if (r > max_pen) max_pen = r;
    }
  }

#if RP_NODE
  /* Get the total average quantities from all compute nodes */
  max_pen    = PRF_GRHIGH1(max_pen);
  total_mass = PRF_GRSUM1(total_mass);
#endif /* RP_NODE */

  /* Divide the max length into several sections */
  nR = 100;
  dr = max_pen/nR;

  /* Compute the penetration length to 99% mass */
  int pen1_done = 0;
  for(i=nR; i>0; i--) {
    loop(p, I->p) {
      if (!p->on_wall) {
        r = sqrt(  (P_POS(p)[0] - I->pos1[0])*(P_POS(p)[0] - I->pos1[0])
                 + (P_POS(p)[1] - I->pos1[1])*(P_POS(p)[1] - I->pos1[1])
                 + (P_POS(p)[2] - I->pos1[2])*(P_POS(p)[2] - I->pos1[2]));
        if((r > (i-1)*dr) && (r < i*dr))
          accum_mass += P_MASS(p)*P_N(p);
      }
    }

#if RP_NODE
    accum_mass = PRF_GRSUM1(accum_mass);
#endif /* RP_NODE */

    if (accum_mass > 0.01*total_mass && !pen1_done) {
      penetration = (i - 0.5)*dr;
      pen1_done = 1;
    }

    if (accum_mass > 0.05*total_mass) {
      *pen_05 = (i - 0.5)*dr;
      break;
    }
  }
  
  *pen_max = max_pen;
  return penetration;
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(DPM_getPenetration)
{
  Injection *I;
  real pen_1, pen_5, pen_max, D32;

  /* TODO: generalize injection name */
  I = Pick_Injection("injection-1");

  /* compute the penetration length */
  pen_1 = DPM_penetration_length(I,&pen_5,&pen_max);
  D32   = DPM_get_SauterD(I);

  Message0("Injection name: %s\n", I->name);
  Message0("Injection location: %15.6e %15.6e %15.6e\n", I->pos1[0], I->pos1[1],I->pos1[2]);

  Message0("Penetration of injection 1: %12.4e\n", pen_1);
  Message0("SMD of injection 1: %12.4e\n", D32);
}
/* ------------------------------------------------------------------------- */

real DPM_get_SauterD(Injection *I)
{
#if !RP_HOST
  Particle *p;
  real D32, den, d2, d3;

  D32 = den = 0.0;
  loop (p, I->p)
    {
      d2 = P_DIAM(p)*P_DIAM(p);
      d3 = d2*P_DIAM(p);
      D32 += p->number_in_parcel*d3;
      den += p->number_in_parcel*d2;
    }

#if RP_NODE
  /* Get the total average quantities from all compute nodes */
  D32 = PRF_GRSUM1(D32);
  den = PRF_GRSUM1(den);
#endif /* RP_NODE */

  /* finish normalization */
  if (den > 0.0)
    D32 /= den;
  else
    D32 = 0.0;

  return D32;
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

DEFINE_EXECUTE_AT_END(writeStats) 
{
#if !RP_HOST
  Message0("Writing spray statistics... ");
  Injection *I, *I_list = Get_dpm_injections();
  Particle  *p;
  int  i, j, nNozzle, lengthID;
  real *penetration, *pen_05, *pen_max, *D32, totMass;
  real *filmMass,*injMass,*rwork;

  nNozzle = 0;
  loop(I, I_list) { nNozzle++; }

  /* allocate memory */
  penetration = (real*)  calloc(nNozzle, sizeof(real));
  pen_05      = (real*)  calloc(nNozzle, sizeof(real));
  pen_max     = (real*)  calloc(nNozzle, sizeof(real));
  D32         = (real*)  calloc(nNozzle, sizeof(real));
  filmMass    = (real*)  calloc(nNozzle, sizeof(real));
  injMass     = (real*)  calloc(nInj,    sizeof(real));
  rwork       = (real*)  calloc(nInj,    sizeof(real));

  /* compute desired parameters for each injection */
  i=0;
  for (j=0; j<nInj; j++) injMass[j] = 0.0;
  loop(I, I_list) {
    penetration[i] = DPM_penetration_length(I,&pen_05[i],&pen_max[i]);
    D32[i]         = DPM_get_SauterD(I);

    /* compute liquid from nozzle i that is on wall */
    filmMass[i]    = 0.0;
    loop(p, I->p) {
      if (p->on_wall)
        filmMass[i] += P_MASS(p)*P_N(p);
    }

    lengthID = strcspn(I->name, "-"); /* NOTE: format must be ID-NOZZLE */
    for (j=0; j<nInj; j++) {
      if (!strncmp(I->name, injID[j], lengthID)) {
        loop(p, I->p) {
          injMass[j] += P_MASS(p)*P_N(p);
        } /* end loop over particles */
      }
    } /* end loop over inj ids*/
    i++;
  } /* end loop over injections */

#if RP_NODE
  /* Get the total average quantities from all compute nodes */
  PRF_GRSUM(injMass,nInj,rwork);
#endif /* RP_NODE */

  for (j=0, totMass=0.0; j<nInj; j++)
    totMass += injMass[j];

#if PARALLEL
  if (myid == node_zero)
#endif /* PARALLEL */
    {
      if (sprayFile == NULL)
        Message0("Warning: sprayFile is NULL");

      /* write time */      
      fprintf(sprayFile, "%8i %15.5e ", N_TIME, CURRENT_TIME);
      if (in_cylinder_p)
        fprintf(sprayFile,"%15.5e ",RP_Get_Real("dynamesh/in-cyn/crank-angle"));

      /* write fuel masses */
      fprintf(sprayFile, "%15.5e ", totMass);
      for (i=0; i<nInj; i++)
        fprintf(sprayFile, "%15.5e ", injMass[i]);

      /* write penetration length (1%) */
      for (i=0; i<nNozzle; i++)
        fprintf(sprayFile, "%15.6e ", penetration[i]);

      /* write penetration length (5%) */
      for (i=0; i<nNozzle; i++)
        fprintf(sprayFile, "%15.6e ", pen_05[i]);

      /* write penetration length max */
      for (i=0; i<nNozzle; i++)
        fprintf(sprayFile, "%15.6e ", pen_max[i]);

      /* write Sauter mean diameter */
      for (i=0; i<nNozzle; i++)
        fprintf(sprayFile, "%15.6e ", D32[i]);

      /* write wall film mass */
      for (i=0; i<nNozzle; i++)
        fprintf(sprayFile, "%15.6e ", filmMass[i]);

      fprintf(sprayFile, "\n");
      fflush(sprayFile);
    }
  free(penetration); free(pen_05); free(pen_max);
  free(D32); free(injMass); free(rwork);
  Message0("done.\n");
#endif /* !RP_HOST */
}
/* ------------------------------------------------------------------------- */

int UDFsprayParseFluentString(char *input, char ***var)
{
  int i, n;
  char *name, *sWork, s[1024];

  if (input == NULL)
    return 0;

  /* copy the input string */
  sWork = (char*) calloc(1024, sizeof(char));
  strcpy(sWork, input);

  var[0] = (char**) calloc(1, sizeof(char*));
  i = 0;
  while(sWork[0])
    {
      /* read in the first zone id and increment */
      sscanf(sWork,"%s%n",s,&n); sWork += n;
      var[0]    = realloc(var[0], (i+1)*sizeof(char*));
      var[0][i] = calloc(512, sizeof(char));
      strcpy(var[0][i], s);
      i++;
    }
/*   free(sWork); */
  return i;
}
/* ------------------------------------------------------------------------- */

void UDFsprayPrefixString(char *pre, char *str)
{
  char buffer[512];
  strcpy(buffer, str);
  strcpy(str, pre);
  strcat(str, buffer);
}
/* -------------------------------------------------------------------------- */

double MAX_D(double x, double y)
{
  if (x > y) return x;
  else return y;
}
/* ------------------------------------------------------------------------- */

double MIN_D(double x, double y)
{
  if (x < y) return x;
  else return y;
}
/* ------------------------------------------------------------------------- */

double UDF_lnGamma(double xx)
{
  int j;
  double x, tmp, y, ser;
  double cof[14] = {57.1562356658629235,     -59.5979603554754912,
                    14.1360979747417471,      -0.491913816097620199,
                     0.339946499848118887e-4,  0.465236289270485756e-4,
                    -0.983744753048795646e-4,  0.158088703224912494e-3,
                    -0.210264441724104883e-3,  0.217439618115212643e-3,
                    -0.164318106536763890e-3,  0.844182239838527433e-4,
                    -0.261908384015814087e-4,  0.368991826595316234e-5};
  y = x = xx;
  tmp = x + 5.24218750000000000;
  tmp = (x + 0.5)*log(tmp) - tmp;
  ser = 0.999999999999997092;
  for (j=0; j<14; j++) ser += cof[j]/++y;
  return tmp + log(2.5066282746310005*ser/x);
}
/* ------------------------------------------------------------------------- */
