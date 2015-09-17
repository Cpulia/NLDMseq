/*
    vbem.h
    a header file of VB-EM for a document.
    $Id: vbem.h,v 1.4 2004/11/06 08:26:19 dmochiha Exp $

*/
#ifndef LDA_VBEM_H
#define LDA_VBEM_H
#include "feature.h"

extern void
vbem (const document *d, double *gamma, double **q,
      double *nt, double *pnt, double *ap,
      const double *alpha, const double **beta,
      int L, int K, int emmax);

#endif
