#include <time.h>
#include "mducUDF.h"

int UDFtimerInit(udfTimer **timer)
{
  /* by default, first timer will be set to total step */
  timer[0] = (udfTimer*) calloc(1, sizeof(udfTimer));
  timer[0][0].nTimer = 1;
  strcpy(timer[0][0].name, "Total");
  timer[0][0].started = 0;
  timer[0][0].now = timer[0][0].in = timer[0][0].out = 0.0;
  return timer[0][0].nTimer;
}
/* -------------------------------------------------------------------------- */

int UDFtimerCreate(char *name, udfTimer **timer)
{
  if (timer[0] == NULL)
    Error("Timer has not been initialised!");

  int i;
  int iT     = timer[0][0].nTimer;
  int nTimer = timer[0][0].nTimer+1;
  timer[0] = (udfTimer*) realloc(timer[0], nTimer*sizeof(udfTimer));
  
  strcpy(timer[0][iT].name, name);
  /* set the times to be zero*/
  timer[0][iT].in = timer[0][iT].out = 0.0;
  /* make sure not started */
  timer[0][iT].started = 0;
  /* make sure all timers have information about the overall */
  for (i=0; i<nTimer; i++)
    timer[0][i].nTimer = nTimer;
  return nTimer;
}
/* -------------------------------------------------------------------------- */

void UDFtimerMonitor(char *file, udfTimer *timer)
{
  int i;
  FILE *pFile = fopen(file, "a");
  double rest;

  /* get the overall time */
#if PARALLEL
  timer[0].dt = timer[0].out - timer[0].in;
#else /* SERIAL */
  timer[0].dt = difftime(timer[0].out, timer[0].in);
#endif /* PARALLEL */

  /* compute the delta for each timer */
  for (i=1, rest=0.0; i<timer[0].nTimer; i++) {
    if (timer[i].dt > 0.0) {
      /* timer was in adjust function, do nothing */
    }
    else {
#if PARALLEL
      timer[i].dt = timer[i].out - timer[i].in;
#else /* SERIAL */
      timer[i].dt = difftime(timer[i].out, timer[i].in);
#endif /* PARALLEL */
    }
    rest += timer[i].dt;
  }

  /* write the iteration and overall time*/
  fprintf(pFile, "%6i %15.6e ", N_TIME, timer[0].dt);

  /* write each timer and its percentage */
  for (i=1; i<timer[0].nTimer; i++) {
    fprintf(pFile, "%15.6e %15.6e ", timer[i].dt,
            timer[i].dt/timer[0].dt*100.0 );
  }
  /* write the rest */
  fprintf(pFile, "%15.6e %15.6e\n", timer[0].dt - rest,
          (timer[0].dt - rest)/timer[0].dt*100.0 );

  fclose(pFile);
}
/* -------------------------------------------------------------------------- */

void UDFtimerStart(char *name, udfTimer *timer)
{
  int i = UDFtimerFind(name, timer);
  timer[i].dt = 0.0;
  timer[i].started = 1;
#if PARALLEL
  timer[i].in  = MPI_Wtime();
  timer[i].now = timer[i].in;
#else /* SERIAL */
  time(&timer[i].in);
  timer[i].now = timer[i].in;
#endif /* PARALLEL */
}
/* -------------------------------------------------------------------------- */

void UDFtimerAdd(char *name, udfTimer *timer)
{
  int i = UDFtimerFind(name, timer);
#if PARALLEL
  double now;
  now = MPI_Wtime();
  timer[i].dt += now - timer[i].now;
  timer[i].now = MPI_Wtime();
#else /* SERIAL */
  time_t now;
  time(&now);
  timer[i].dt += difftime(now, timer[i].now);
  time(&timer[i].now);
#endif /* PARALLEL */
}
/* -------------------------------------------------------------------------- */

void UDFtimerStop(char *name, udfTimer *timer)
{
  int i = UDFtimerFind(name, timer);
  timer[i].started = 0;
#if PARALLEL
  timer[i].out = MPI_Wtime();
#else /* SERIAL */
  time(&timer[i].out);
#endif /* PARALLEL */
}
/* -------------------------------------------------------------------------- */

void UDFtimerSet(char *name, udfTimer *timer)
{
  int i = UDFtimerFind(name, timer);
#if PARALLEL
  timer[i].now  = MPI_Wtime();
#else /* SERIAL */
  time(&timer[i].now);
#endif /* PARALLEL */
}
/* -------------------------------------------------------------------------- */

int UDFtimerStatus(char *name, udfTimer *timer)
{
  int i = UDFtimerFind(name, timer);
  return timer[i].started;
}
/* -------------------------------------------------------------------------- */

int UDFtimerFind(char *name, udfTimer *timer)
{
  int i;
  /* Message0("timer names: "); */
  for (i=0; i<timer[0].nTimer; i++) {
    /* Message0("%s ", timer[i].name); */
    if (!strcmp(timer[i].name, name))
      return i;
  }

  Error("Timer %s not found!", name);
}
/* -------------------------------------------------------------------------- */

