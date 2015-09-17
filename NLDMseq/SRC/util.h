/*
    util.h
    $Id: util.h,v 1.6 2004/11/10 04:23:06 dmochiha Exp $

*/
#ifndef LDA_UTIL_H
#define LDA_UTIL_H

extern int  myclock (void);
extern char *rtime (double t);
extern char *strconcat(const char *s, const char *t);
extern int  converged (double *u, double *v, int n, double threshold);
extern int  doublecmp (double *x, double *y);
extern void normalize_matrix_row (double **dst, double **src, int rows, int cols);
extern void normalize_matrix_col (double **dst, double **src, int rows, int cols);

#endif
