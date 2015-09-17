
#ifndef __READ_BETA_H__
#define __READ_BETA_H__
int beta_size_cnt(char* file_name,int *beta_size);
int read_beta(char*file_name,double** Map);
int normal_beta(double** Map ,int nclass,int nlex);
#endif
