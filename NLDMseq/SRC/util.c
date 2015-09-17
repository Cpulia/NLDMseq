/*
    util.c
    $Id: util.c,v 1.7 2004/11/10 04:23:06 dmochiha Exp $

*/
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "util.h"

int myclock ()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec;
}

char *
rtime (double t)
{
	int hour, min, sec;
	static char buf[BUFSIZ];

	hour = (int)floor((int)t / 60 / 60);
	min  = (int)floor(((int)t % (60 * 60)) / 60);
	sec  = (int)floor((int)t % 60);
	sprintf(buf, "%2d:%02d:%02d", hour, min, sec);
	
	return buf;
}

int
doublecmp (double *x, double *y)
{
	return (*x == *y) ? 0 : ((*x < *y) ? 1 : -1);
}

char *
strconcat (const char *s, const char *t)
{
	static char z[BUFSIZ];
	strcpy(z, s);
	strcat(z, t);
	return z;
}

int
converged (double *u, double *v, int n, double threshold)
{
	/* return 1 if |a - b|/|a| < threshold */
	double us = 0;
	double ds = 0;
	double d;
	int i;
	
	for (i = 0; i < n; i++)
		us += u[i] * u[i];

	for (i = 0; i < n; i++) {
		d = u[i] - v[i];
		ds += d * d;
	}

	if (sqrt(ds / us) < threshold)
		return 1;
	else
		return 0;

}

void
normalize_matrix_col (double **dst, double **src, int rows, int cols)
{
	/* column-wise normalize from src -> dst */
	double z;
	int i, j;

	for (j = 0; j < cols; j++) {
		for (i = 0, z = 0; i < rows; i++)
			z += src[i][j];
		for (i = 0; i < rows; i++)
			dst[i][j] = src[i][j] / z;
	}
}

void
normalize_matrix_row (double **dst, double **src, int rows, int cols)
{
	/* row-wise normalize from src -> dst */
	int i, j;
	double z;

	for (i = 0; i < rows; i++) {
		for (j = 0, z = 0; j < cols; j++)
			z += src[i][j];
		for (j = 0; j < cols; j++)
			dst[i][j] = src[i][j] / z;
	}
}


		
