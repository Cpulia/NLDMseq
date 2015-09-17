/*
    vbem.c
    Latent Dirichlet Allocation, VB-EM for a document.
    $Id: vbem.c,v 1.6 2004/11/06 08:26:19 dmochiha Exp $

*/
#include <math.h>
#include "lda.h"
#include "gamma.h"
#include "feature.h"
#include "util.h"

void
vbem (const document *d, double *gamma, double **q,
      double *nt, double *pnt, double *ap,
      const double *alpha, const double **beta,
      int L, int K, int emmax)
{
	int j, k, l;
	double z;

	for (k = 0; k < K; k++)
		nt[k] = (double) L / K;
	
	for (j = 0; j < emmax; j++)
	{
		/* vb-estep */
		for (k = 0; k < K; k++)
			ap[k] = exp(psi(alpha[k] + nt[k]));
		/* accumulate q */
		for (l = 0; l < L; l++)
			for (k = 0; k < K; k++)
				q[l][k] = beta[d->id[l]][k] * ap[k];
		/* normalize q */
		for (l = 0; l < L; l++) {
			z = 0;
			for (k = 0; k < K; k++)
				z += q[l][k];
			for (k = 0; k < K; k++)
				q[l][k] /= z;
		}
		/* vb-mstep */
		for (k = 0; k < K; k++) {
			z = 0;
			for (l = 0; l < L; l++)
				z += q[l][k] * d->cnt[l];
			nt[k] = z;
		}
		/* converge? */
		if ((j > 0) && converged(nt, pnt, K, 1.0e-2))
			break;
		for (k = 0; k < K; k++)
			pnt[k] = nt[k];
	}
	for (k = 0; k < K; k++)
		gamma[k] = alpha[k] + nt[k];
	
	return;
}

