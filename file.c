#include <stdio.h>
#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef LCS_KNC
#include <offload.h>
#endif

#include <immintrin.h>
#include "file.h"

/*
**初始化变量，通过statistics_samples()得到序列数目，最大序列长度等
*/
void init(FILE * fp){

	statistics_samples(fp);
	read_file_size = get_file_size(fp);
	block_max_seq_num = READ_BLOCK_SIZE / (seq_len + 1);//假设全是序列，标题长度为0，READ_BLOCK_SIZE最多几条序列，分配这么大的空间
    	read32_len = (seq_len + BIT_PARA_MIC - 1) / BIT_PARA_MIC ;
   	read_fill_len = read32_len * BIT_PARA_MIC ;
}



/*
**以某种模式打开文件
*/
FILE * open_file(const char * filename, const char * mode) {
    FILE * fp;
    fp = fopen(filename, mode);
    if(fp == NULL) {
        printf("Error - can't open or create file: %s\n", filename);
        exit(1);
    }
    return fp;
}



/*
**得到read文件大小
*/
int64_t get_file_size(FILE *data){
	fseek(data, 0, SEEK_END);
	int64_t size;
	size = ftell(data);
	fseek(data, 0, SEEK_SET);
    return size;
}



/*
**得到ref文件大小
*/
int64_t get_ref_file_size(FILE *ref_data){
	fseek(ref_data, 0, SEEK_END);
	int64_t size;
	size = ftell(ref_data);
	fseek(ref_data, 0, SEEK_SET);
    return size;
}




/*
**统计fasta文件，得到序列数目，最大序列长度
*/
void statistics_samples(FILE * fp){
    	//初始化read，读取文件内容，再处理
	char *read = (char* ) malloc (CHUNK_SIZE * sizeof(char));
	memset(read,0,CHUNK_SIZE);
	int size = 1;
	int count = 0;
	int offset = 0;
	int title_count = 0;
	char c;//当前读取到的字符
	int flag = 0;//标志是否为新的序列
	int k =0;
	while(size){
		size = fread(read, sizeof(char) , CHUNK_SIZE , fp); 
		for(k = 0 ; k < size ; k ++){
			c = read[k];
			//遇到了新序列
			if(c == '>'){
				count = 0;
				flag++;
				title_count = 0;
				k++;
				while(read[k]!='\n'){
					title_count++;
					k++;
				}
			}else {//序列内容
				if(read[k]!='\n'){
					count++;
				}
			}
		}
	}

	reads_num = flag;
	seq_len = count;
	free(read);
	fseek(fp,0,SEEK_SET);//回到文件头
	
}



/*
**读fasta文件，读取section_size大小
*/
void get_read_from_file(seq_t * seq, FILE * fp, int64_t section_size){
	char* r_read = (char* ) malloc ((section_size + 2) * sizeof(char));
	memset(r_read,0,section_size + 2);
	char c;//当前读取到的字符
	int flag = 0;//标志是否为新的序列
	int64_t count = 0;
	int offset = 0;
	memset(seq->content,0,block_max_seq_num * read_fill_len);
	int position[block_max_seq_num];//第几条序列的开始位置
	int i = 0;
	read_actual_size = fread(r_read,sizeof(char),section_size,fp);

	for(i = 0; i< read_actual_size;i++){
		//遇到了新序列
		if(r_read[i] == '>' ){
			position[flag] = i;
			flag ++;
			count = 0;
			i++;
			while(r_read[i]!='\n'){
				i++;
			}
		}else {//序列内容
			if(r_read[i]!='\0' & r_read[i]!='\n'){
				seq->content[(flag - 1) * read_fill_len + count] = r_read[i];
				count++;
			}
		}
	}//对应i

	int remainder= 0;
	//移动指针头
	if(read_actual_size >= READ_BLOCK_SIZE){
		remainder = (flag - 1) % 16;
		block_num = flag - 1 - remainder;	//将序列数目变为16的倍数
	}else{
		block_num = flag;	//最后一次不管几条序列，是否16的倍数，都收着
	}
	offset = read_actual_size - position[block_num] + 1;
	//printf("block_num is %d\n",block_num);
    	seq->size = read_actual_size;
	seq->num = block_num;
	fseek(fp, -offset + 1, SEEK_CUR);
	free(r_read);
}



/*
**读取ref文件，转化为int*
*/
void get_ref_from_file(seq_t * ref, FILE * ref_data,int ref_size){
	int i = 0;
	ref_len = 0;
	char *tmp_ref = (char* )_mm_malloc(sizeof(char) * ref_size, 64);
	memset(tmp_ref,0,ref_size);
	memset(ref->seq,0,ref_size);
	fread(tmp_ref,sizeof(char),ref_size,ref_data);
	for(i = 0; ;i++) {
        if(tmp_ref[i] == '\n') {
            break;
        }
        ref_len++;
    }
	ref->len = ref_len;
	ref_num = ref_size / (ref_len + 1);
	ref->num  = ref_num;
	ref->size = ref_size;

	for(i = 0; i < ref_size; i ++){
		if(tmp_ref[i] != '\0' & tmp_ref[i] != '\n'){
			if(tmp_ref[i] == 'A'){
				ref->seq[i] = 0;
			}else if(tmp_ref[i] == 'G') {
				ref->seq[i] = 1;
			}else if(tmp_ref[i] == 'C') {
				ref->seq[i] = 2;
			}else if(tmp_ref[i] == 'T') {
				ref->seq[i] = 3;
			}
		}
	}

	_mm_free(tmp_ref);
}









