#ifndef MDUC_FLUENT_UDF_TIMING_H
#define MDUC_FLUENT_UDF_TIMING_H

#include <time.h>

typedef struct MDUC_UDF_TIMER_STRUCT {
  int    nTimer;
  char   name[25];
  int    started;
  double dt;
#if PARALLEL
  double now, in, out;
#else /* SERIAL */
  time_t now, in, out;
#endif /* PARALLEL */
};
typedef struct MDUC_UDF_TIMER_STRUCT udfTimer;

/* initialise timers */
int UDFtimerInit(udfTimer **timer);
int UDFtimerCreate(char *name, udfTimer **timer);

/* for writing monitor file */
void UDFtimerMonitor(char *file, udfTimer *timer);

/* for updating timers */
void UDFtimerStart(char *name, udfTimer *timer);
void UDFtimerAdd(char *name, udfTimer *timer);
void UDFtimerStop(char *name, udfTimer *timer);
void UDFtimerSet(char *name, udfTimer *timer);

/* find out about a timer */
int UDFtimerStatus(char *name, udfTimer *timer);
int UDFtimerFind(char *name, udfTimer *timer);

#endif /* MDUC_FLUENT_UDF_TIMING_H */
