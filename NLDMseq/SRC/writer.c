/*
    writer.c
    an implementation of vector/matrix writer.
    $Id: writer.c,v 1.3 2004/11/06 09:44:42 dmochiha Exp $

*/
#include <stdio.h>
#include "writer.h"

void
lda_write (FILE *ap, FILE *bp, double *alpha, double **beta,
	   int nclass, int nlex)
{
	//printf("writing model..\n"); fflush(stdout);
	write_vector(ap, alpha, nclass+1);
	write_matrix(bp, beta, nlex, nclass+1);
	//printf("done.\n"); fflush(stdout);
}

void
write_vector (FILE *fp, double *vector, int n)
{
	int i;
	for (i = 0; i < n; i++)
		fprintf(fp, "%.7e%s", vector[i], (i == n - 1) ? "\n" : "   ");
}

void
write_matrix (FILE *fp, double **matrix, int rows, int cols)
{
	int i, j;
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++)
			fprintf(fp, "%.7e%s", matrix[i][j],
				(j == cols - 1) ? "\n" : "   ");
}
void
write_log(char *exp_name,double*isoexpre,int nclass,char* gene_name)
{
    int cnt=0;
    FILE* input =fopen(exp_name,"a");
    fprintf(input,"%s\t",gene_name);
    for(cnt=0;cnt<nclass;cnt++)
    {
        fprintf(input,"%.7e\t",isoexpre[cnt]);
    }
    fprintf(input,"%s","\n");
    fclose(input);
}
void
write_OUTCOME(char *exp_name,double*iso_expre,double*iso_var,double gene_expression,double gene_variance,int nclass,char* gene_name)
{
    int cnt=0;
    FILE* input =fopen(exp_name,"w");
    fprintf(input,"%s\t",gene_name);
    for(cnt=0;cnt<nclass;cnt++)
    {
        fprintf(input,"%.7e\t%.7e\t",iso_expre[cnt],iso_var[cnt]);
    }
    fprintf(input,"%.7e\t%.7e",gene_expression,gene_variance);
    fprintf(input,"%s","\n");
    fclose(input);
}
