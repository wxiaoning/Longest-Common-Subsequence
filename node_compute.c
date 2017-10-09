#include "node_compute.h"
#include <sys/time.h>
#include "timer.h"
#include "file.h"
#include "seqs_handle.h"

#include <immintrin.h>
#ifdef LCS_KNC
	#include <offload.h>
	#include "lcs.h"
#else 
#ifdef LCS_KNL
	#include "lcs_knl.h"
#endif
#endif

#include <math.h>
#include <mpi.h>

void node_compute(int group_rank,int rank, FILE* data,FILE* result,FILE* ref_data){

	double read_ref_start, read_ref_end;
	GET_TIME(read_ref_start);
	
	sem_init(&compute_a,0,1);
	sem_init(&compute_b,0,0);
	sem_init(&read_a,0,0);
	sem_init(&read_b,0,1);
	sem_init(&write_a,0,0);
	sem_init(&write_b,0,0);
	sem_init(&output_compute_a,0,1);
	sem_init(&output_compute_b,0,1);

	init(data);
	int64_t ref_size = get_ref_file_size(ref_data);

	reference.seq = (int* )_mm_malloc(sizeof(int) * ref_size,64);
	get_ref_from_file(&reference, ref_data,ref_size);

	int *ref = (int* )_mm_malloc(sizeof(int) * ref_size,64);
	ref = reference.seq;
	GET_TIME(read_ref_end);
	//printf("read ref file time is %.2fs\n\n",read_ref_end - read_ref_start);

	/*if(group_rank == 0){
		printf("reads_num is %d\n",reads_num);	
		printf("seq_len is %d\n",seq_len);
		printf("ref_num is %d\n",ref_num);
	}
	*/	
	read1a = (mic_read_t* )_mm_malloc( block_max_seq_num * read32_len * sizeof(mic_read_t),64);
	read2a = (mic_read_t* ) _mm_malloc( block_max_seq_num * read32_len * sizeof(mic_read_t),64);
	
	read1b = (mic_read_t* )_mm_malloc( block_max_seq_num * read32_len * sizeof(mic_read_t),64);
	read2b = (mic_read_t* ) _mm_malloc( block_max_seq_num * read32_len * sizeof(mic_read_t),64);

	align_results_a = (mic_write_t *)_mm_malloc(sizeof(mic_write_t) * block_max_seq_num * read32_len * ref_num,64);
	align_results_b = (mic_write_t *)_mm_malloc(sizeof(mic_write_t) * block_max_seq_num * read32_len * ref_num,64);

    results_a = (uint16_t *)_mm_malloc(sizeof(uint16_t) * block_max_seq_num * ref_num, 32);
    results_b = (uint16_t *)_mm_malloc(sizeof(uint16_t) * block_max_seq_num * ref_num, 32);

	double time_start, time_end;
	GET_TIME(time_start);

	//分配内存
   	if(read_file_size > READ_BLOCK_SIZE){
		//若总序列数目大于READ_BLOCK_SIZE,则申请2块读序列内存
		read_seq_a.content = (char*) _mm_malloc (sizeof(char) * block_max_seq_num * read_fill_len,64);
		memset(read_seq_a.content, 0 ,block_max_seq_num * read_fill_len);
	   	read_seq_b.content = (char*) _mm_malloc (sizeof(char) * block_max_seq_num * read_fill_len,64);
		memset(read_seq_b.content, 0 ,block_max_seq_num * read_fill_len);
	}else{
		//若总序列数目小于READ_BLOCK_SIZE,则申请1块read文件大小内存
		read_seq_a.content = (char*) _mm_malloc (sizeof(char) * (block_max_seq_num * read_fill_len),64);
		memset(read_seq_a.content, 0 ,block_max_seq_num * read_fill_len);
	   	read_seq_b.content = (char*) _mm_malloc (sizeof(char) * (block_max_seq_num * read_fill_len),64);
		memset(read_seq_b.content, 0 ,block_max_seq_num * read_fill_len);
		read_bucket_num = 1;
	}

	//读一块read序列
	get_read_from_file(&read_seq_a, data, READ_BLOCK_SIZE);
	read_seq_a.len = seq_len;

	if(read_file_size > READ_BLOCK_SIZE){
		read_bucket_num = (read_file_size + read_seq_a.size - 1) /  read_seq_a.size;//需要分几块读
	}

	mic_cpu_counts = (int64_t *)_mm_malloc(read_bucket_num * sizeof(int64_t),64);//记录每block处理的序列数目
	mic_cpu_counts[0] = block_num;
	device_read_counts_a = block_num;

	//处理第一块序列，将其转换成int
	mic_seqs_handle(read1a,read2a,read_seq_a.content,device_read_counts_a);
	int i = 0;
	
	GET_TIME(time_end);
   	//printf("the first block handle time  is %.2fs\n", time_end - time_start);
#ifdef LCS_KNC	
	double mic_alloc_start, mic_alloc_end;
	GET_TIME(mic_alloc_start);
	
	//在MIC上分配内存

	#pragma offload target(mic:group_rank)\
		nocopy(ref:length(ref_size) alloc_if(1) free_if(0))\
		nocopy(read1a:length(block_max_seq_num * read32_len) alloc_if(1) free_if(0))\
		nocopy(read2a:length(block_max_seq_num * read32_len) alloc_if(1) free_if(0))\
		nocopy(read1b:length(block_max_seq_num * read32_len) alloc_if(1) free_if(0))\
		nocopy(read2b:length(block_max_seq_num * read32_len) alloc_if(1) free_if(0))\
		nocopy(align_results_a:length(block_max_seq_num * read32_len * ref_num) alloc_if(1) free_if(0))\
		nocopy(align_results_b:length(block_max_seq_num * read32_len * ref_num) alloc_if(1) free_if(0))\
        nocopy(results_a:length(block_max_seq_num * ref_num) alloc_if(1) free_if(0))\
        nocopy(results_b:length(block_max_seq_num * ref_num) alloc_if(1) free_if(0))
		{}
	

	GET_TIME(mic_alloc_end);
   	printf("the malloc on MIC time  is %.2fs\n\n", mic_alloc_end - mic_alloc_start);
	
	
	//传第一块到MIC
	#pragma offload target(mic:group_rank) \
		in(read32_len,ref_len,ref_num) \
		in(ref:length(ref_size) alloc_if(0) free_if(0)) \
		in(read1a:length(block_num * read32_len) alloc_if(0) free_if(0)) \
		in(read2a:length(block_num * read32_len) alloc_if(0) free_if(0)) 
		{}
#endif
double cpu_start_time,cpu_end_time, mic_start_time, mic_end_time, offload_start_time,offload_end_time,output_start_time,output_end_time;
	omp_set_nested(1);

	double start, end;
	GET_TIME(start);
    
    double local_start, local_end, local_elapsed, elapsed;
    MPI_Barrier(MPI_COMM_WORLD);
    local_start = MPI_Wtime();

#pragma omp parallel for num_threads(3)//开启3个线程，分别mic计算、handle+offload、output
	for(int thread_id = 0; thread_id < 3; thread_id ++){
		
		int iter = 0;
		int read_bucket_index = 0;
		int pid = omp_get_thread_num();

		if(pid == 0){//mic compute thread
			while(1){
				if( iter % 2 == 0){  //a读完了，计算A

					sem_wait(&compute_a);
					sem_wait(&output_compute_a);
					GET_TIME(mic_start_time);
					#ifdef LCS_KNC
					
					#pragma offload target(mic:group_rank) \
					nocopy(ref:length(0) alloc_if(0) free_if(0)) \
					nocopy(read1a:length(0) alloc_if(0) free_if(0)) \
					nocopy(read2a:length(0) alloc_if(0) free_if(0)) \
					nocopy(align_results_a:length(0) alloc_if(0) free_if(0))\
					out(results_a:length(device_read_counts_a * ref_num) alloc_if(0) free_if(0))
					{
					    mic_kernel(ref_len,ref_num,read32_len,device_read_counts_a,read1a,read2a,align_results_a,results_a,ref);
					}
					#else
					#ifdef LCS_KNL
					
						knl_kernel(ref_len, ref_num, read32_len, device_read_counts_a, read1a, read2a, align_results_a, results_a, reference.seq);
				
					#endif
					#endif
			
					read_bucket_index++;
					iter ++;

					sem_post(&read_a);
					sem_post(&write_a);

					GET_TIME(mic_end_time);
					printf("the mic_computeA time is %.2fs\n\n",mic_end_time - mic_start_time);
					if(read_bucket_index > read_bucket_num - 1) {//计算完成
						break;
        			}
				}else if(iter % 2 == 1){//计算B

					sem_wait(&compute_b);
					sem_wait(&output_compute_b);
					GET_TIME(mic_start_time);
					#ifdef LCS_KNC
						
					#pragma offload target(mic:group_rank) \
					nocopy(ref:length(0) alloc_if(0) free_if(0)) \
					nocopy(read1b:length(0) alloc_if(0) free_if(0)) \
					nocopy(read2b:length(0) alloc_if(0) free_if(0)) \
					nocopy(align_results_b:length(0) alloc_if(0) free_if(0)) \
					out(results_b:length(device_read_counts_b * ref_num) alloc_if(0) free_if(0))
					{
    					mic_kernel(ref_len,ref_num,read32_len,device_read_counts_b,read1b,read2b,align_results_b,results_b,ref);
					}
					#else
					#ifdef LCS_KNL
					
					knl_kernel(ref_len, ref_num, read32_len, device_read_counts_b, read1b, read2b, align_results_b, results_b, reference.seq);
					
					#endif
					#endif
											
					read_bucket_index++;
					iter ++;

					sem_post(&read_b);
					sem_post(&write_b);

					GET_TIME(mic_end_time);
					printf("the mic_computeB time is %.2fs\n\n",mic_end_time - mic_start_time);
					if(read_bucket_index > read_bucket_num - 1) {//计算完成
						break;
					}
				}
			}

		}else if(pid == 1){//offload thread

			int bucket_index = 1;
			int i = 0, j = 0;
			if(read_bucket_num != 1){
    			while(1){
				if(iter % 2 == 0){//读入B

					sem_wait(&read_b);

					//GET_TIME(offload_start_time);

					get_read_from_file(&read_seq_b, data, READ_BLOCK_SIZE);//读序列
					device_read_counts_b = block_num;
					mic_cpu_counts[bucket_index] = block_num;//mic_cpu_counts记录每block处理的序列数目
										
					mic_seqs_handle(read1b,read2b,read_seq_b.content,device_read_counts_b);
					#ifdef LCS_KNC
					#pragma offload target(mic:group_rank) \
					in(read1b:length(block_num * read32_len) alloc_if(0) free_if(0)) \
					in(read2b:length(block_num * read32_len) alloc_if(0) free_if(0)) 
					{}
					#endif
			        printf("bucket_index: %d, group_rank: %d, device_read_counts_b: %d\n",bucket_index,group_rank, device_read_counts_b);	
					//GET_TIME(offload_end_time);

					sem_post(&compute_b);

					bucket_index ++;
					iter ++;
					if(bucket_index > read_bucket_num -1) {
						break;
					}
				}else if(iter % 2 == 1){//读入A

					sem_wait(&read_a);

					//GET_TIME(offload_start_time);
				
					get_read_from_file(&read_seq_a, data, READ_BLOCK_SIZE);//读序列
					mic_cpu_counts[bucket_index] = block_num;//mic_cpu_counts记录每block处理的序列数目
					device_read_counts_a = block_num;

					mic_seqs_handle(read1a,read2a,read_seq_a.content,device_read_counts_a);
					#ifdef LCS_KNC
					#pragma offload target(mic:group_rank) \
					in(read1a:length(block_num * read32_len) alloc_if(0) free_if(0)) \
					in(read2a:length(block_num * read32_len) alloc_if(0) free_if(0)) 
					{}
					#endif
					//GET_TIME(offload_end_time);
			        printf("bucket_index: %d, group_rank: %d, device_read_counts_b: %d\n",bucket_index, group_rank, device_read_counts_a);	

					sem_post(&compute_a);

					bucket_index ++;
					iter ++;
					if(bucket_index > read_bucket_num - 1) {
						break;
					}
				}
			}
			}
		}else if(pid == 2){//output thread
			int result_counts = 0;

			while(1){
				if(iter % 2 == 0){

					sem_wait(&write_a);
					//GET_TIME(output_start_time);
					result_counts = mic_cpu_counts[read_bucket_index] * ref_num;
                    if(read_bucket_index == (read_bucket_num - 1)){
                        //计算TCPUS
                        local_end = MPI_Wtime();
                        local_elapsed = local_end - local_start;
                        
                        printf("rank is %d, local_elpased is %.2f\n",rank, local_elapsed);
                        MPI_Reduce(&local_elapsed,&elapsed,1,MPI_DOUBLE,MPI_MAX,0 ,MPI_COMM_WORLD);
                        long total_counts = 0;
                        for(int i = 0; i < read_bucket_num; ++i){
                            total_counts += mic_cpu_counts[i];
                        }    

		    #ifdef LCS_KNC
                    if(rank == 0){
		    #endif
                        int node_num = 0;
                        MPI_Comm_size(MPI_COMM_WORLD,&node_num);
                        float Tcups = 0.00;                     
                        Tcups = reference.len * reference.num * total_counts * seq_len * node_num / elapsed / pow(10,12);
                        printf("total_counts is %d\n",total_counts);
                        printf("Elapsed time is %.2fs\n", elapsed);
					    printf("Tcups is %.2fT\n",Tcups);
                    #ifdef LCS_KNC    
		    }
		    #endif
                    }
                    
                    //printf("result_counts is %d\n",result_counts); 
                    //for(int i = 0; i < result_counts; ++i){
                    //    fprintf(result,"%d:%u\n", i, (unsigned)results_a[i]);
                    //}

                    fwrite(&result_counts, sizeof(uint32_t), 1, result);
                    fwrite(results_a, sizeof(uint16_t), result_counts, result);   
                    
					fflush(result);

					iter++;
					read_bucket_index ++;

					//GET_TIME(output_end_time);
				    //printf("\n result_count is %d\n",result_counts);
				    //printf("read_bucket_index is %d\n",read_bucket_index);
				    //printf("the Block A output time is %.2fs\n",output_end_time - output_start_time);

					sem_post(&output_compute_a);
					if(read_bucket_index > read_bucket_num - 1) {
						break;
					}
				}else if(iter % 2 == 1){

					sem_wait(&write_b);
					GET_TIME(output_start_time);
					result_counts = mic_cpu_counts[read_bucket_index] * ref_num;
					
                    if(read_bucket_index == (read_bucket_num - 1)){
                        //计算TCPUS
                        local_end = MPI_Wtime();
                        local_elapsed = local_end - local_start;
                        
                        printf("rank is %d, local_elpased is %.2f\n",rank, local_elapsed);
                        MPI_Reduce(&local_elapsed,&elapsed,1,MPI_DOUBLE,MPI_MAX,0 , MPI_COMM_WORLD);
                        long total_counts = 0;
                        for(int i = 0; i < read_bucket_num; ++i){
                            total_counts += mic_cpu_counts[i];
                        }  
		    #ifdef LCS_KNC
                    if(rank == 0){
		    #endif
                        int node_num = 0;
                        MPI_Comm_size(MPI_COMM_WORLD,&node_num);
                        float Tcups = 0.00; 
                        Tcups = reference.len * reference.num * total_counts * seq_len * node_num / elapsed / pow(10,12);
                        printf("total_counts is %d\n",total_counts);
                        printf("group_rank 0 Elapsed time is %.2fs\n", elapsed);
	        	        printf("Tcups is %.2fT\n",Tcups);
                    #ifdef LCS_KNC    
		    }
		    #endif
                    }
                    
                    //printf("result_counts is %d\n",result_counts); 
                    //for(int i = 0; i < result_counts; ++i){
                    //    fprintf(result,"%d:%u\n", i , (unsigned)results_b[i]);
                    //}

                    fwrite(&result_counts, sizeof(uint32_t), 1, result);
                    fwrite(results_b, sizeof(uint16_t), result_counts, result);

					fflush(result);

					iter++;
					read_bucket_index++;

					GET_TIME(output_end_time);
					//printf("the Block B output time is %.2fs\n",output_end_time - output_start_time);

					sem_post(&output_compute_b);
					if(read_bucket_index > read_bucket_num - 1) {
						break;
					}
				}
			}
		}

	}

	sem_destroy(&compute_a);
	sem_destroy(&compute_b);
	sem_destroy(&read_a);
	sem_destroy(&read_b);
	sem_destroy(&write_a);
	sem_destroy(&write_b);
	sem_destroy(&output_compute_a);
	sem_destroy(&output_compute_b);

	GET_TIME(end);
	printf("the total compute time is %.2fs\n", end - start);
#ifdef LCS_KNC
#pragma offload target(mic:group_rank) \
	nocopy(ref:length(ref_size) alloc_if(0) free_if(1))\
	nocopy(read1a:length(block_max_seq_num * read32_len) alloc_if(0) free_if(1))\
	nocopy(read2a:length(block_max_seq_num * read32_len) alloc_if(0) free_if(1))\
	nocopy(read1b:length(block_max_seq_num * read32_len) alloc_if(0) free_if(1))\
	nocopy(read2b:length(block_max_seq_num * read32_len) alloc_if(0) free_if(1))\
	nocopy(align_results_a:length(block_max_seq_num * read32_len * ref_num) alloc_if(0) free_if(1))\
	nocopy(align_results_b:length(block_max_seq_num * read32_len * ref_num) alloc_if(0) free_if(1))\
	nocopy(results_a:length(block_max_seq_num * ref_num) alloc_if(0) free_if(1))\
	nocopy(results_b:length(block_max_seq_num * ref_num) alloc_if(0) free_if(1))
	{}
#endif

    _mm_free(read1a);
    _mm_free(read2a);
    _mm_free(read1b);
    _mm_free(read2b);
    
    _mm_free(read_seq_a.content);
    _mm_free(read_seq_b.content);
    _mm_free(reference.seq);
    _mm_free(align_results_a);
    _mm_free(align_results_b);
    _mm_free(results_a);
    _mm_free(results_b);

    read1a = NULL;
    read2a = NULL;

    read1b = NULL;
    read2b = NULL;

    read_seq_a.content = NULL;
    read_seq_b.content = NULL;

    reference.seq = NULL;
    align_results_a = NULL;
    align_results_b = NULL;

    results_a = NULL;
    results_b = NULL;

}














