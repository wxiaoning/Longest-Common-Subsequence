#ifndef NODE_COMPUTE_H_
#define NODE_COMPUTE_H_

#include <sys/time.h>
#include "global.h"


void node_compute(int group_rank,int rank, FILE* data,FILE* result,FILE* ref_data);
void mic_start(void* args);
void cpu_start(void* args);


#endif/*NODE_COMPUTE_H_*/
