#include <stdio.h>
#include "seqs_handle.h"
#include <omp.h>



/*
alphabet size 为4
*/

void mic_seqs_handle(uint32_t* read1,uint32_t* read2,char *seqs,int nums){
    //读取序列，转换，存储
	//printf("nums is %d\n\n",nums);
memset(read1,0,sizeof(read1));
memset(read2,0,sizeof(read2));

#pragma omp parallel for 
	for(int k=0; k < nums ; k++){
		int j = 0,i = 0;
		char current;//当前字符
		uint32_t bitmask;
		int r[4][2] = {{0,0},{0,-1},{-1,0},{-1,-1}};

		uint32_t readc1[read32_len], readc2[read32_len];
		int quotients = 0,remainder = 0;//商、余数
		//将序列转换为int
		quotients = k / MIC_PARA_NUM;
		remainder = k % MIC_PARA_NUM;

		for(i=0; i< read32_len ;i++){
        		bitmask = 0x00000001;
			readc1[i] = 0;
			readc2[i] = 0;
			
			for(j=0; j< BIT_PARA_MIC;j++){
				current = seqs[k * read_fill_len + i * BIT_PARA_MIC + j] ;
				if(current == 'A') {
					readc1[i] |= (r[0][0] & bitmask);
					readc2[i] |= (r[0][1] & bitmask);
				}else if(current == 'G') {
					readc1[i] |= (r[1][0] & bitmask);
					readc2[i] |= (r[1][1] & bitmask);
				}else if(current == 'C') {
					readc1[i] |= (r[2][0] & bitmask);
					readc2[i] |= (r[2][1] & bitmask);
				}else if(current == 'T') {
					readc1[i] |= (r[3][0] & bitmask);
					readc2[i] |= (r[3][1] & bitmask);
				}		

				bitmask <<= 1;
			}

			read1[quotients * read32_len * MIC_PARA_NUM  + i * MIC_PARA_NUM + remainder] = readc1[i];
			read2[quotients * read32_len * MIC_PARA_NUM  + i * MIC_PARA_NUM + remainder] = readc2[i];
			
		}

	}
	


}



