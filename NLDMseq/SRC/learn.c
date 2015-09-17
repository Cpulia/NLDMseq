/*
    learn.c
    $Id: learn.c,v 1.8 2013/01/13 14:23:27 daichi Exp $

*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "learn.h"
#include "vbem.h"
#include "newton.h"
#include "likelihood.h"
#include "feature.h"
#include "dmatrix.h"
#include "util.h"

void
lda_learn (document *data, double *alpha, double **beta,
	   int nclass, int nlex, int dlenmax,
	   int emmax, int demmax, double epsilon)
{
	double *gamma, **q, **gammas, **betas, *nt, *pnt, *ap;
	double ppl, pppl = 0;
	double z;
	document *dp;
	int i, j, t, n;
	int start, elapsed;

	/*
	 *  randomize a seed
	 *
	 */
	srand(time(NULL));

	/*
	 *    count data length
	 *
	 */
	for (dp = data, n = 0; (dp->len) != -1; dp++, n++)
		;

	/*
	 *  initialize parameters
	 *
	 */
	for (i = 0; i < nclass; i++)
		alpha[i] = RANDOM;
	/*for (i = 0, z = 0; i < nclass; i++)
		z += alpha[i];
	for (i = 0; i < nclass; i++)
		alpha[i] = alpha[i] / z;
	qsort(alpha, nclass, sizeof(double), // sort alpha initially
	      (int (*)(const void *, const void *))doublecmp);

	for (j = 0; j < nclass; j++) {
		for (i = 0, z = 0; i < nlex; i++) {
			beta[i][j] = RANDOM * 10;
			z += beta[i][j];
		}
		for (i = 0; i < nlex; i++) {
			beta[i][j] = beta[i][j] / z;
		}
	}*/

	/*
	 *  initialize posteriors
	 *
	 */
	if ((gammas = dmatrix(n, nclass)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate gammas.\n");
		return;
	}
	if ((betas = dmatrix(nlex, nclass)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate betas.\n");
		return;
	}
	/*
	 *  initialize buffers
	 *
	 */
	if ((q = dmatrix(dlenmax, nclass)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate q.\n");
		return;
	}
	if ((gamma = (double *)calloc(nclass, sizeof(double))) == NULL)
	{
		fprintf(stderr, "lda_learn:: cannot allocate gamma.\n");
		return;
	}
	if ((ap = (double *)calloc(nclass, sizeof(double))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate ap.\n");
		return;
	}
	if ((nt = (double *)calloc(nclass, sizeof(double))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate nt.\n");
		return;
	}
	if ((pnt = (double*)calloc(nclass, sizeof(double))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate pnt.\n");
		return;
	}

	/*printf("Number of documents          = %d\n", n);
	printf("Number of words              = %d\n", nlex);
	printf("Number of latent classes     = %d\n", nclass);
	printf("Number of outer EM iteration = %d\n", emmax);
	printf("Number of inner EM iteration = %d\n", demmax);
	printf("Convergence threshold        = %g\n", epsilon);
	*/

	/*
	 *  learn main
	 *
	 */
	//start = myclock();
	for (t = 0; t < emmax; t++)
	{
		//printf("iteration %2d/%3d..\t", t + 1, emmax);
		fflush(stdout);
		/*
		 *  VB-E step
		 *
		 */
		/* iterate for data */
		for (dp = data, i = 0; (dp->len) != -1; dp++, i++)
		{
			vbem(dp, gamma, q, nt, pnt, ap,
			     alpha, (const double **)beta, dp->len, nclass, demmax);
			accum_gammas(gammas, gamma, i, nclass);
			accum_betas(betas, q, nclass, dp);
		}
		/*
		 *  VB-M step
		 *
		 */
		/* Newton-Raphson for alpha */
		newton_alpha(alpha, gammas, n, nclass, 0);
		/* MLE for beta */
		normalize_matrix_col(beta, betas, nlex, nclass);
		/* clean buffer */
		for (i = 0; i < nlex; i++)
			for (j = 0; j < nclass; j++)
				betas[i][j] = 0;
		/*
		 *  converge?
		 *
		 */
		ppl = lda_ppl (data, beta, gammas, n, nclass);
		elapsed = myclock() - start;
		//printf("PPL = %.04f\t", ppl); fflush(stdout);
		if ((t > 1) && (fabs((ppl - pppl)/ppl) < epsilon)) {
			if (t < 2) {
				free_dmatrix(gammas, n);
				free_dmatrix(betas, nlex);
				free_dmatrix(q, dlenmax);
				free(gamma);
				free(ap);
				free(nt);
				free(pnt);
				//printf("\nearly convergence. restarting..\n");
				lda_learn (data, alpha, beta, nclass, nlex,
					   dlenmax, emmax, demmax, epsilon);
				return;
			} else {
				//printf("\nconverged. [%s]\n", rtime(elapsed));
				break;
			}
		}
		pppl = ppl;
		/*
		 * ETA
		 *
		 */
		/*printf("ETA:%s (%d sec/step)\r",
		       rtime(elapsed * ((double) emmax / (t + 1) - 1)),
		       (int)((double) elapsed / (t + 1) + 0.5));*/
	}

	free_dmatrix(gammas, n);
	free_dmatrix(betas, nlex);
	free_dmatrix(q, dlenmax);
	free(gamma);
	free(ap);
	free(nt);
	free(pnt);

	return;
}

void
accum_gammas (double **gammas, double *gamma, int n, int K)
{
	/* gammas(n,:) = gamma for Newton-Raphson of alpha */
	int k;
	for (k = 0; k < K; k++)
		gammas[n][k] = gamma[k];
	return;
}

void
accum_betas (double **betas, double **q, int K, document *dp)
{
	int i, k;
	int n = dp->len;

	for (i = 0; i < n; i++)
		for (k = 0; k < K; k++)
			betas[dp->id[i]][k] += q[i][k] * dp->cnt[i];
}


