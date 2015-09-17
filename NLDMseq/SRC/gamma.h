/*
    gamma.h
    header file for digamma and trigamma functions.
    $Id: gamma.h,v 1.1 2004/11/01 04:52:02 dmochiha Exp $

*/
#ifndef __GAMMA_H__
#define __GAMMA_H__

double digamma(double x);
double trigamma(double x);

#define psi(x)	digamma(x)
#define ppsi(x)	trigamma(x)

#endif
