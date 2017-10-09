#include <stdio.h>
#include <string.h>
#include <stdlib.h>



void get_rand_str(int seqlen,int seqnum,FILE* file){//产生一行长度为seqlen的字符串，写入到fasta

	char* str="AGCT";
	int flag;//flag为对应的整数值
	int i = 0,k = 0;
	srand((unsigned)time(NULL));
	//char *ch = malloc(sizeof(char) * (seqlen+1));
	char ch[seqlen+2];
    char *read = "GTTGCGGAGATAACCCTGACGGACCGCCATAA";
	//写100 0000条read
	for(k=0;k<seqnum;k++)  {
		fprintf(file,"%s",">adsafra\n");
       		fprintf(file,"%s",read);
 		for(i=0;i<seqlen;i++){//生成seqlen个String
  			flag=rand() % 4;
  			ch[i]=str[flag];

   		}	
		ch[seqlen] = '\n';
		ch[seqlen+1] = '\0';
   		fprintf(file,"%s",ch);
 	}
	//free(ch);
}


int main(int argc,char* argv[]){
	int seqlen = 4000 - 32 ,seqnum = 188000 * 64;  //subject序列长度、序列数目，atoi为windows下将字符串转换为整数的函数
	const char *file_name = "./database/seq.fasta";
	FILE *fasta = fopen(file_name,"w+");
	if(!fasta){
		printf("error:can't create file\n");
		exit(0);
	}

	get_rand_str(seqlen,seqnum,fasta);
	fclose(fasta);
	
    return 0;
}
