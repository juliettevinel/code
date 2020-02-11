/*
 * linsys_solvers.h
 *
 *  Created on: Jan 27, 2014
 *      Author: paant
 */

#ifndef LINSYS_SOLVERS_H_
#define LINSYS_SOLVERS_H_



#include <iostream>
#include <stdio.h>
//#include <string.h>
#include <stdlib.h>
//#include <assert.h>
#include <math.h>
#include <errno.h>

using namespace std;


void BCSG(double *sL, int *ijL, double *X, double *D, double *Pre,  double error, int L, int N, int& flag);


//Functions to be used throughout the computations//


double* array_n (int n) ;
void delarray_n (double* mat) ;

double* array_np1 (int n) ;
void delarray_np1 (double* mat) ;


void tridag (double* a, double* b, double* c, double* r, double* u, int n,
		double* gam, int& flag, int k1) ;


#endif /* LINSYS_SOLVERS_H_ */
