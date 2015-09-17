/*
    writer.h
    a header file of vector/matrix writer.
    $Id: writer.h,v 1.2 2004/11/10 04:23:06 dmochiha Exp $

*/
#ifndef LDA_WRITER_H
#define LDA_WRITER_H
#include <stdio.h>

extern void lda_write (FILE *ap, FILE *bp, double *alpha, double **beta,
		       int nclass, int nlex);
void write_vector (FILE *fp, double *vector, int n);
void write_matrix (FILE *fp, double **matrix, int rows, int cols);
void write_log(char *exp_name,double*isoexpre,int nclass,char* gene_name);
void write_OUTCOME(char *exp_name,double*iso_expre,double*iso_var,double gene_expression,double gene_variance,int nclass,char* gene_name);
#endif
