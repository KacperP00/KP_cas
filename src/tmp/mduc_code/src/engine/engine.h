#ifndef MDUC_ENGINE_H
#define MDUC_ENGINE_H

#include "mducStructs.h"

/* initialization functions */
void defaultEngine(engineStruct *engine);
void initEngine(ioStruct *io, chemStruct *chem, engineStruct *engine);

void initEngineMass(void *mem);
double getEngineMass(char *type, void *mem);
double calcEnginePressure(double t, engineStruct *engine, flowStruct *flow, 
                          flameStruct *flame);

void writeEngineManifold(double t, engineStruct *engine);
void writeEngineMonitorHeader(ioStruct *io, engineStruct *engine);
void writeEngineMonitor(double t, ioStruct *io, chemStruct *chem, 
                        engineStruct *engine, flameStruct *flame);

void calcEngineRHS_0D(double t, timeStruct *time, chemStruct *chem,
                      flowStruct *flow, flameStruct *flame,engineStruct *engine,
                      double *Y, double *RHS);

void calcEngineSpark(double t, driveStruct *driver, timeStruct *time,
                     chemStruct *chem, flowStruct *flow, flameStruct *flame,
                     engineStruct *engine, double *Y, double *RHS);

double engineCalcFlameArea(double t, engineStruct *engine);

void destroyEngine(engineStruct *engine);
void destroyValve(valveStruct *valve);
void destroyInjection(injStruct *inj);

#endif /* MDUC_ENGINE_H */
