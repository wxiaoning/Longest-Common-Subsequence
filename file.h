#ifndef _FILE_H_
#define _FILE_H_


#include "global.h"


FILE * open_file(const char * filename, const char* mode);
void statistics_samples(FILE * fp);

void get_read_from_file(seq_t * seq, FILE * fp, int64_t section_size);
int64_t get_file_size(FILE * fp);

void get_ref_from_file(seq_t * ref, FILE * fp,int ref_size);
int64_t get_ref_file_size(FILE * fp);

void init( FILE * fp);

#endif /*_FILE_H_*/
