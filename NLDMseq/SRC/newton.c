/*
    newton.c
    Newton-Raphson iteration for alpha of LDA
    $Id: newton.c,v 1.3 2004/11/04 13:43:43 dmochiha Exp $

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "newton.h"
#include "gamma.h"
#include "util.h"

void
newton_alpha (double *alpha, double **gammas, int M, int K, int level)
{
	int i, j, t;
	double *g, *h, *pg, *palpha;
	double z, sh, hgz;
	double psg, spg, gs;
	double alpha0, palpha0;

	/* allocate arrays */
	if ((g = (double *)calloc(K, sizeof(double))) == NULL) {
		fprintf(stderr, "newton:: cannot allocate g.\n");
		return;
	}
	if ((h = (double *)calloc(K, sizeof(double))) == NULL) {
		fprintf(stderr, "newton:: cannot allocate h.\n");
		return;
	}
	if ((pg = (double *)calloc(K, sizeof(double))) == NULL) {
		fprintf(stderr, "newton:: cannot allocate pg.\n");
		return;
	}
	if ((palpha = (double *)calloc(K, sizeof(double))) == NULL) {
		fprintf(stderr, "newton:: cannot allocate palpha.\n");
		return;
	}

	/* initialize */
	if (level == 0)
	{
		for (i = 0; i < K; i++) {
			for (j = 0, z = 0; j < M; j++)
				z += gammas[j][i];
			alpha[i] = z / (M * K);
		}
	} else {
		for (i = 0; i < K; i++) {
			for (j = 0, z = 0; j < M; j++)
				z += gammas[j][i];
			alpha[i] = z / (M * K * pow(10, level));
		}
	}

	psg = 0;
	for (i = 0; i < M; i++) {
		for (j = 0, gs = 0; j < K; j++)
			gs += gammas[i][j];
		psg += psi(gs);
	}
	for (i = 0; i < K; i++) {
		for (j = 0, spg = 0; j < M; j++)
			spg += psi(gammas[j][i]);
		pg[i] = spg - psg;
	}

	/* main iteration */
	for (t = 0; t < MAX_NEWTON_ITERATION; t++)
	{
		for (i = 0, alpha0 = 0; i < K; i++)
			alpha0 += alpha[i];
		palpha0 = psi(alpha0);

		for (i = 0; i < K; i++)
			g[i] = M * (palpha0 - psi(alpha[i])) + pg[i];
		for (i = 0; i < K; i++)
			h[i] = - 1 / ppsi(alpha[i]);
		for (i = 0, sh = 0; i < K; i++)
			sh += h[i];

		for (i = 0, hgz = 0; i < K; i++)
			hgz += g[i] * h[i];
		hgz /= (1 / ppsi(alpha0) + sh);

		for (i = 0; i < K; i++)
			alpha[i] = alpha[i] - h[i] * (g[i] - hgz) / M;

		for (i = 0; i < K; i++)
			if (alpha[i] < 0) {
				if (level >= MAX_RECURSION_LIMIT) {
					fprintf(stderr, "newton:: maximum recursion limit reached.\n");
					return;
				} else {
					free(g);
					free(h);
					free(pg);
					free(palpha);
					return newton_alpha(alpha, gammas, M, K, 1 + level);
				}
			}

		if ((t > 0) && converged(alpha, palpha, K, 1.0e-4)) {
			free(g);
			free(h);
			free(pg);
			free(palpha);
			return;
		} else
			for (i = 0; i < K; i++)
				palpha[i] = alpha[i];

	}
	//fprintf(stderr, "newton:: maximum iteration reached. t = %d\n", t);

	return;

}
