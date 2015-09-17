/*
    learn.h
    $Id: learn.h,v 1.5 2004/11/06 08:26:19 dmochiha Exp $

*/
#ifndef LDA_LEARN_H
#define LDA_LEARN_H
#include <stdlib.h>
#include "feature.h"

#define CONVERGENCE 1.0e-2
#define RANDOM ((double)rand()/(double)RAND_MAX)
extern double round (double x);

extern void lda_learn (document *data, double *alpha, double **beta,
		       int nclass, int nlex, int dlenmax,
		       int emmax, int demmax, double epsilon);

void accum_gammas (double **gammas, double *gamma, int n, int K);
void accum_betas (double **betas, double **q, int K, document *dp);
void normalize_matrix_col (double **dst, double **src, int rows, int cols);


#endif
