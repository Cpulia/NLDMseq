#ifndef GENE_EXPRESSION_H_INCLUDED
#define GENE_EXPRESSION_H_INCLUDED
#include "feature.h"
int create_sampling(double *alpha,double *theta,int nclass);
int read_mapping(double **isoform,document*sum_exon,double *theta,double **Map,double nclass,double nlex);
double sum_gene_exon(double *p_exon,int len);
void free_sum_exon(document* sum_exon);
void expr_sig(double ** isoform,int nclass,int nlex, double* isolen,double *iso_expre,int SEQ_DEPTH);
int read_isolen(char *isolen_name,double* isolen,char t_lex);
#endif // GENE_EXPRESSION_H_INCLUDED
