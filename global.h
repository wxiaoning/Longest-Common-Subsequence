#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <omp.h>

typedef struct _seq_t{
    int len;
    int64_t size;
    int64_t num;
    char* content;//read content
    int* seq;//ref content
    int read32_len;
}seq_t;


#define mic_read_t uint32_t
#define mic_write_t uint32_t
#define sse_read_t uint32_t
#define sse_write_t uint32_t
#define MAX_REF 5000
#define BIT_PARA_MIC 32
#define MIC_PARA_NUM 16
#define CHUNK_SIZE 602403360
#define READ_BLOCK_SIZE 200801120	//大约190M 200801120  100400560
#define READ_TITLE_MAX_LEN 50

extern sem_t compute_a;
extern sem_t compute_b;
extern sem_t read_a;
extern sem_t read_b;
extern sem_t write_a;
extern sem_t write_b;
extern sem_t output_compute_a;
extern sem_t output_compute_b;

#ifdef LCS_KNC
#define __ONMIC__ __attribute__((target(mic)))
extern __ONMIC__ mic_read_t *read1a, *read2a, *read3a, *read4a, *read5a;//其中一块数据缓冲区，处理的数据放入其中
extern __ONMIC__ mic_read_t *read1b, *read2b, *read3b, *read4b, *read5b;//另一块
extern __ONMIC__ mic_write_t * align_results_a;
extern __ONMIC__ mic_write_t * align_results_b;
extern __ONMIC__ uint16_t *results_a, *results_b;

extern __ONMIC__ int64_t device_read_counts_a;
extern __ONMIC__ int64_t device_read_counts_b;
extern __ONMIC__ int reads_num;
extern __ONMIC__ int ref_num;
extern __ONMIC__ int ref_len;
extern __ONMIC__ int read32_len;
extern __ONMIC__ int read31_len;
extern __ONMIC__ seq_t reference;  //读取进来的ref序列内容
#else
extern mic_read_t *read1a, *read2a, *read3a, *read4a, *read5a;//其中一块数据缓冲区，处理的数据放入其中
extern mic_read_t *read1b, *read2b, *read3b, *read4b, *read5b;//另一块
extern mic_write_t * align_results_a;
extern mic_write_t * align_results_b;
extern uint16_t *results_a, *results_b;

extern int64_t device_read_counts_a;
extern int64_t device_read_counts_b;
extern int reads_num;
extern int ref_num;
extern int ref_len;
extern int read32_len;
extern int read31_len;
extern seq_t reference;  //读取进来的ref序列内容
#endif


//在这里打开read文件，第一次读取也用这里的
extern int64_t remain_size;//文件未读取大小
extern int bucket_num ;
extern seq_t read_seq_a, read_seq_b;
extern int64_t *mic_cpu_counts;//记录每块mic计算的序列数目
extern int read_bucket_num;
extern int block_num ;
extern int read_fill_len;
extern int seq_len;
extern int64_t read_actual_size;
extern int block_max_seq_num;
extern uint64_t read_file_size;


extern void print_binary( uint32_t t,int bit_len);

#endif /*_GLOBAL_H_*/
