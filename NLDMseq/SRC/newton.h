/*
    newton.h
    a header file for Newton-Raphson.
    $Id: newton.h,v 1.2 2004/11/04 13:43:43 dmochiha Exp $

*/
#ifndef LDA_NEWTON_H
#define LDA_NEWTON_H
#define MAX_RECURSION_LIMIT  20
#define MAX_NEWTON_ITERATION 20

extern void newton_alpha (double *alpha, double **gammas,
			  int nlen, int nclass, int level);


#endif
