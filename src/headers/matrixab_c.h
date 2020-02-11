/*
 * matrixab_b.h
 *
 *  Created on: Jan 11, 2011
 *      Author: paant
 */

#ifndef MATRIXAB_C_H_
#define MATRIXAB_C_H_



#include <math.h>

#include <iostream>

//#include "global.h"
#include "trivialfuncs.h"
#include "matrixab.h"

using namespace std;

void constr_a_c (double* a, int* pos, double* precond, double** phi, int n, int m, grid cells ) ;

void constr_a_b (double* a, int* pos, double* precond, double** phi, int n, int m, grid cells ) ;

void constr_b_b (double* output, double** utilde_i, double** utilde_j, double** rho1, double** rho0, double** rho_new, double** phi1, double** phi0, double** phi_new, double** M_re, int n, int m, double* y_val, grid cells, double dt) ;

#endif /* MATRIXAB_C_H_ */
