#include "lcs_knl.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <immintrin.h>

/*
2个bit,处理alphabet size 为4
*/
void knl_kernel(int ref_len, int ref_num, int read32_len, int mic_num, mic_read_t* read1, mic_read_t* read2, mic_write_t* align_results, uint16_t *result, int* ref_content){
	const int EVERY_LOOP = 5;
   	__m512i max, min, one, set_zero, set_one;
   	max = _mm512_set1_epi32(0xffffffff);
   	min = _mm512_set1_epi32(0);
   	one = _mm512_set1_epi32(1);
	set_zero = _mm512_set1_epi32(0x7fffffff);//最高位为0，其余为1
	set_one = _mm512_set1_epi32(0x80000000);//最高位为1，其余为0
	int i = 0;
#pragma omp parallel for num_threads(252)
	for(i = 0 ;i < mic_num ; i += MIC_PARA_NUM){//mic_num即mic_num应该为16的倍数,对于16条序列，循环ref_num次

    	__m512i tem1, tem2, tem3, tem4, tem5, tem6, tem7, tem8, tem9, notem1, notem2, notem3, notem4;
    	__m512i LP[read32_len + 1],LL[read32_len];//LP=Lj-1,LL=Lj
    	__m512i match[4][EVERY_LOOP], notmatch[4][EVERY_LOOP];
    	__m512i read1v[EVERY_LOOP],read2v[EVERY_LOOP];
    	__mmask16 carry[ref_len];
    	int current ;

	int loop = read32_len / EVERY_LOOP;
	int remainder = read32_len % EVERY_LOOP;
	int t = 0, l = 0, k = 0, j = 0, r = 0;

	for(int ref_i = 0; ref_i < ref_num; ref_i ++){
		//初始化LP
	    	for(t=0;t < read32_len; t++){
			LP[t] = _mm512_set1_epi32(0xffffffff);
			LL[t] = _mm512_set1_epi32(0);
	     	}

	    	//初始化carry
	    	for(t=0;t < ref_len; t++){
	        	carry[t] = _mm512_int2mask(0);
	     	}
		for(l = 0; l < loop ; l++){  //几次循环,当不能整除，有余数时，最后一次应该单独计算
      			for(k=0; k < EVERY_LOOP; k++){ //每次循环几个word
      				read1v[k] = _mm512_load_epi32(read1 + i * read32_len + MIC_PARA_NUM * k + l * EVERY_LOOP * MIC_PARA_NUM);
      				read2v[k] = _mm512_load_epi32(read2 + i * read32_len + MIC_PARA_NUM * k + l * EVERY_LOOP * MIC_PARA_NUM);      			
			
				tem1 = _mm512_xor_epi32(read1v[k],min);
        			tem2 = _mm512_xor_epi32(read2v[k],min);
        			tem3 = _mm512_xor_epi32(read1v[k],max);
        			tem4 = _mm512_xor_epi32(read2v[k],max);
        			notem1 = _mm512_xor_epi32(tem1,max);
        			notem2 = _mm512_xor_epi32(tem2,max);
        			notem3 = _mm512_xor_epi32(tem3,max);
        			notem4 = _mm512_xor_epi32(tem4,max);
	
        			//计算matchA、notmatchA  00
        			match[0][k] =  _mm512_and_epi32(notem1,notem2);
        			notmatch[0][k] = _mm512_xor_epi32(match[0][k],max);

				//计算matchG、notmatchG  01
        			match[1][k] =  _mm512_and_epi32(notem1,notem4);
        			notmatch[1][k] = _mm512_xor_epi32(match[1][k],max);

        			//计算matchC、notmatchC  10
        			match[2][k] =  _mm512_and_epi32(notem2,notem3);
        			notmatch[2][k] = _mm512_xor_epi32(match[2][k],max);
        		
        			//计算matchT、notmatchT  11  
        			match[3][k] =  _mm512_and_epi32(notem3,notem4);
        			notmatch[3][k] = _mm512_xor_epi32(match[3][k],max);
        		}//对应k

			for(j=0; j < ref_len; j++){
        			current = ref_content[j + ref_i * (ref_len + 1)];//ref为不变的那个序列，只有一条
				for(r=0; r < EVERY_LOOP; r++ ){
      					tem5 = _mm512_and_epi32(LP[r + l * EVERY_LOOP],match[current][r]);
     					//tem6 = _mm512_and_epi32(LP[r + l * EVERY_LOOP],notmatch[current][r]);

					tem7 = _mm512_add_epi32(LP[r + l * EVERY_LOOP],tem5);
					tem7 = _mm512_mask_add_epi32(tem7,carry[j],tem7,one);
					//carry[j] = _mm512_cmplt_epu32_mask(tem7, set_one);		
					carry[j] = _mm512_cmplt_epu32_mask(tem7,tem5);
					//LL[r + l * EVERY_LOOP] = _mm512_or_epi32(tem7, tem6);	
					//LP[r + l * EVERY_LOOP] = LL[r + l * EVERY_LOOP];
					//LP[r + l * EVERY_LOOP] = _mm512_and_epi32(LP[r + l * EVERY_LOOP],set_zero); 
					//LP[r + l * EVERY_LOOP] = _mm512_ternarylogic_epi32(set_zero, tem7, tem6, 0xe0);

					LP[r + l * EVERY_LOOP] = _mm512_ternarylogic_epi32(tem7, LP[r + l * EVERY_LOOP], notmatch[current][r], 0xf8);
				}//r
    			}//j
		//_mm512_store_epi32(align_results + read32_len * i * ref_num + MIC_PARA_NUM * l + ref_i * read32_len * MIC_PARA_NUM, LL[0]);
 		}//l

		//最后一个循环，计算remainder个word
		for(k=0; k < remainder; k++){
      			read1v[k] = _mm512_load_epi32(read1 + i * read32_len + MIC_PARA_NUM * k + loop * EVERY_LOOP * MIC_PARA_NUM);
      			read2v[k] = _mm512_load_epi32(read2 + i * read32_len + MIC_PARA_NUM * k + loop * EVERY_LOOP * MIC_PARA_NUM);
	
			tem1 = _mm512_xor_epi32(read1v[k],min);
        		tem2 = _mm512_xor_epi32(read2v[k],min);
        		tem3 = _mm512_xor_epi32(read1v[k],max);
        		tem4 = _mm512_xor_epi32(read2v[k],max);
        		notem1 = _mm512_xor_epi32(tem1,max);
        		notem2 = _mm512_xor_epi32(tem2,max);
        		notem3 = _mm512_xor_epi32(tem3,max);
        		notem4 = _mm512_xor_epi32(tem4,max);

        		//计算matchA、notmatchA  00
        		match[0][k] =  _mm512_and_epi32(notem1,notem2);
        		notmatch[0][k] = _mm512_xor_epi32(match[0][k],max);

			//计算matchG、notmatchG  01
        		match[1][k] =  _mm512_and_epi32(notem1,notem4);
        		notmatch[1][k] = _mm512_xor_epi32(match[1][k],max);

        		//计算matchC、notmatchC  10
        		match[2][k] =  _mm512_and_epi32(notem2,notem3);
        		notmatch[2][k] = _mm512_xor_epi32(match[2][k],max);
        		
        		//计算matchT、notmatchT  11  
        		match[3][k] =  _mm512_and_epi32(notem3,notem4);
        		notmatch[3][k] = _mm512_xor_epi32(match[3][k],max);

    		}//对应k


		for(j=0; j < ref_len; j++){
        		current = ref_content[j + ref_i * (ref_len + 1)];//ref为不变的那个序列，只有一条
			for(r=0; r < remainder; r++ ){
      			
				tem5 = _mm512_and_epi32(LP[r + loop * EVERY_LOOP],match[current][r]);
     				//tem6 = _mm512_and_epi32(LP[r + loop * EVERY_LOOP],notmatch[current][r]);


				tem7 = _mm512_add_epi32(LP[r + loop * EVERY_LOOP],tem5);
				tem7 = _mm512_mask_add_epi32(tem7,carry[j],tem7,one);
				//carry[j] = _mm512_cmplt_epu32_mask(tem7, set_one);				
				carry[j] = _mm512_cmplt_epu32_mask(tem7, tem5);
				
				//LP[r + loop * EVERY_LOOP] = _mm512_ternarylogic_epi32(set_zero, tem7, tem6, 0xe0);
				
				//LL[r + loop * EVERY_LOOP] = _mm512_or_epi32(tem7, tem6);	
				//LP[r + loop * EVERY_LOOP] = LL[r + loop * EVERY_LOOP];
				//LP[r + loop * EVERY_LOOP] = _mm512_and_epi32(LP[r + loop * EVERY_LOOP],set_zero); 
	    		
				LP[r + loop * EVERY_LOOP] = _mm512_ternarylogic_epi32(tem7, LP[r + loop * EVERY_LOOP], notmatch[current][r], 0xf8);	    		
			}//r
//			printf("\n");
		}//j
		//将最终结果保存到result
		for(t=0; t < read32_len; t++){
			//LP[t] = _mm512_and_epi32(LP[t],set_zero);
   			_mm512_store_epi32(align_results + read32_len * i * ref_num + MIC_PARA_NUM * t + ref_i * read32_len * MIC_PARA_NUM, LP[t]);
   		}
	}//ref_i
}//i
	//printf("align_results[0] is %x\n",align_results[0]);
	//printf("align_results[16] is %x\n",align_results[16]);
	post_process(align_results, result, ref_len, ref_num, read32_len, mic_num);

}



void post_process(mic_write_t *align_results, uint16_t *result, int ref_len, int ref_num, int read32_len, int mic_num){
	uint16_t shift_bit = ref_len % BIT_PARA_MIC;
#pragma omp parallel for
	for(int i = 0; i < mic_num; i += MIC_PARA_NUM)
		for(int j = 0; j < MIC_PARA_NUM * ref_num; j += MIC_PARA_NUM)
			for(int r = 0; r < MIC_PARA_NUM; r++){
				int bitCount = 0, bCount = 0;
				for(int k = 0; k < read32_len - 1; k++){
					bCount = _mm_popcnt_u32(align_results[i * read32_len * ref_num + j * read32_len + k * MIC_PARA_NUM + r]);
					bitCount += bCount;
				}
				mic_write_t last_word = align_results[i * read32_len * ref_num + j * read32_len + (read32_len - 1) * MIC_PARA_NUM + r];
				if(shift_bit != 0){
	        	             mic_write_t tmp_word = (1 << shift_bit) - 1;
 		                     last_word &= tmp_word;
                 		}
				bitCount += _mm_popcnt_u32(last_word);
				result[i * ref_num + j + r] = ref_len - bitCount;
			}


}




