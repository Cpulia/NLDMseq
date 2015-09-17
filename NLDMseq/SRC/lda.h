/*
    lda.h
    header file for LDA.
    $Id: lda.h,v 1.6 2004/11/10 04:23:06 dmochiha Exp $

*/
#ifndef LDA_LDA_H
#define LDA_LDA_H
#include <stdlib.h>
#include "feature.h"

#define LDA_COPYRIGHT "lda, a Latent Dirichlet Allocation package.\nCopyright (C) 2004 Daichi Mochihashi, All rights reserved.\n"
#define CLASS_DEFAULT		50
#define EMMAX_DEFAULT		100
#define DEMMAX_DEFAULT		20
#define EPSILON_DEFAULT		1.0e-4
int gene_expre ( char * data_path,char* gene_name,int SEQ_DEPTH,char *output_path);

#endif
