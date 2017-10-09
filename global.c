#include "global.h"


sem_t compute_a;
sem_t compute_b;
sem_t read_a;
sem_t read_b;
sem_t write_a;
sem_t write_b;
sem_t output_compute_a;
sem_t output_compute_b;


#ifdef LCS_KNC
__ONMIC__ mic_read_t *read1a, *read2a, *read3a, *read4a, *read5a;//其中一块数据缓冲区，处理的数据放入其中
__ONMIC__ mic_read_t *read1b, *read2b, *read3b, *read4b, *read5b;//另一块
__ONMIC__ mic_write_t * align_results_a;
__ONMIC__ mic_write_t * align_results_b;
__ONMIC__ uint16_t *results_a, *results_b;

__ONMIC__ int64_t device_read_counts_a;
__ONMIC__ int64_t device_read_counts_b;
__ONMIC__ int reads_num;
__ONMIC__ int ref_num;
__ONMIC__ int ref_len;
__ONMIC__ int read32_len;
__ONMIC__ int read31_len;
__ONMIC__ seq_t reference;  //读取进来的read序列内容
#else
mic_read_t *read1a, *read2a, *read3a, *read4a, *read5a;//其中一块数据缓冲区，处理的数据放入其中
mic_read_t *read1b, *read2b, *read3b, *read4b, *read5b;//另一块
mic_write_t * align_results_a;
mic_write_t * align_results_b;
uint16_t *results_a, *results_b;

int64_t device_read_counts_a;
int64_t device_read_counts_b;
int reads_num;
int ref_num;
int ref_len;
int read32_len;
int read31_len;
seq_t reference;  //读取进来的read序列内容
#endif

//在这里打开read文件，第一次读取也用这里的
int64_t remain_size;//文件未读取大小
int bucket_num ;
seq_t read_seq_a, read_seq_b;
int64_t *mic_cpu_counts;//记录每块mic计算的序列数目
int read_bucket_num;
int block_num ;
int read_fill_len;
int seq_len;
int64_t read_actual_size;
int block_max_seq_num;
uint64_t read_file_size;


void print_binary( uint32_t t,int bit_len);
