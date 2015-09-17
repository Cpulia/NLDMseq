/*
    likelihood.h
    $Id: likelihood.h,v 1.2 2013/01/13 14:24:05 daichi Exp $

*/
#ifndef LDA_LIKELIHOOD_H
#define LDA_LIKELIHOOD_H
#include "feature.h"

extern double lda_lik (document *data, double **beta, double **gammas,
		       int m, int nclass);
extern double lda_ppl (document *data, double **beta, double **gammas,
		       int m, int nclass);

#endif
