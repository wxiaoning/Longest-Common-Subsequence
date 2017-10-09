#ifndef LCS_H_
#define LCS_H_

#include <omp.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include "global.h"

#ifdef LCS_KNC
__ONMIC__ void mic_kernel(int ref_len, int ref_num, int read32_len, int reads_num, mic_read_t* read1, mic_read_t* read2, mic_write_t* align_results, uint16_t *results, int* ref_content);

__ONMIC__ void mic_post_process(mic_read_t *align_results, uint16_t *results, int ref_len, int ref_num, int read32_len, int mic_num);
#endif

#endif/*LCS_H_*/

//GTTGCGGAGATAACCCTGACGGACCGCCATAA
//00100000000000010000010000000110
