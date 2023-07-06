#ifndef MDUC_DESTROY_H
#define MDUC_DESTROY_H

#include "mducStructs.h"

/* functions to deallocate structs individually */
void destroyFlame(flameStruct *flame);
void destroyGrad(gradStruct *grad);
void destroyFlow(flowStruct *flow);
void destroyIO(ioStruct *io);
void destroySDial(int solver, sdialStruct *sdial);

#endif /* MDUC_DESTROY_H */
