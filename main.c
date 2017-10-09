#include <omp.h>
#include <immintrin.h>

#ifdef LCS_KNC
#include "offload.h"
#endif

#include "node_compute.h"
#include "global.h"
#include "file.h"
#include <mpi.h>

#define INFO_LEN 20
#define SEQUENCE_MAX_LEN 200000

void split_database(int nodes){
    const char *fasta_file = "./newdatabase/seq.fasta";
    FILE* data=fopen(fasta_file,"r");
    if(!data){
        printf("error:can't open input file\n");
        exit(0);
    }
    statistics_samples(data);
    int64_t reads_num_pnode = reads_num / nodes;
    //构造要切分的文件名
    const char *s1 = "./newdatabase/data64/input_seq";
    char s2[10];
    const char *s3 = ".fasta";
    char read_file[50];
    int n;
    for(int i = 0; i < nodes; ++i){
        n = reads_num_pnode;
        memset(s2,0,sizeof(s2));
        memset(read_file,0,sizeof(read_file));
        sprintf(s2,"%d",i);
        sprintf(read_file,"%s%s",s1,s2);
        sprintf(read_file,"%s%s",read_file,s3);
        //打开文件
        FILE* sub_data=fopen(read_file,"w+");
        if(!sub_data){
            printf("error:can't open input file on rank %d\n",i);
            exit(0);
        }
        //读取一定条数的序列
        char *r_read = (char*) malloc ((seq_len + INFO_LEN) * sizeof(char));
        memset(r_read,0,(seq_len + INFO_LEN));
        while(n){
            fgets(r_read,INFO_LEN,data);
            fputs(r_read,sub_data);
            memset(r_read,0,(seq_len + INFO_LEN));
            fgets(r_read,SEQUENCE_MAX_LEN,data);
            fputs(r_read,sub_data);
            memset(r_read,0,(seq_len + INFO_LEN));
            n--;
        }
        fclose(sub_data);
        free(r_read);
    }
}


int main(int argc, char* argv[]){

	int rank;
	int node_num;

    MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &node_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef LCS_KNC
    //一个节点内的rank
	MPI_Comm group_comm;
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &group_comm);
	int group_rank;
	MPI_Comm_rank(group_comm, &group_rank);

    printf("group_rank %d rank %d\n", group_rank, rank);
#else
	int group_rank = 0;
#endif

/*
    if(rank == 0){
        split_database(2);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
*/
    char read_file[50];
	char ref_file[] = "./newdatabase/ref.fasta";
	char result_file[50];

	const char *s1 = "./newdatabase/data64/input_seq";
	char s2[10];
	const char *s3 = ".fasta";
	const char *s4 = "./newdatabase/data64/result";

	//构造要读取的的文件名input_seq.fasta
	sprintf(s2,"%d", rank); 
	sprintf(read_file,"%s%s",s1,s2);
	sprintf(read_file,"%s%s",read_file,s3);	

	//打开read文件
	FILE* data=fopen(read_file,"r");
	if(!data){
		printf("error:can't open input file on rank %d\n",rank);
		exit(0);
	}

	//构造要写入的中间文件名result.fasta
	sprintf(result_file,"%s%s",s4,s2);
	sprintf(result_file,"%s%s",result_file,s3);	

	//打开result文件
	FILE* result = fopen(result_file,"w+");
	if(!result){
		printf("error:can't open file\n");
		exit(0);
	}

	//打开ref文件
	FILE* ref_data=fopen(ref_file,"r");
	if(!ref_data){
	  printf("error:can't open input_ref file\n");
	  exit(0);
	}

	//计算
	node_compute(group_rank,rank, data, result, ref_data);


    	fclose(data);
    	fclose(result);
	fclose(ref_data);


	MPI_Finalize();
	return 0;

}











