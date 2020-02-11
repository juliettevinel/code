/*
 * matrixab_b.h
 *
 *  Created on: Jan 11, 2011
 *      Author: paant
 */

#ifndef MATRIXAB_B_H_
#define MATRIXAB_B_H_



#include <math.h>
#include <iostream>


//#include "global.h"
#include "trivialfuncs.h"
#include "matrixab.h"

using namespace std;


void constr_a_b (double* a, int* pos, double* precond, double** phi, int n, int m, grid** cells );
void constr_b_b (double* output, double** utild_i, double** utild_j, double** phi, double** rho, int n , int m ,grid ** cells);
//void constr_b_b_bot (double* output, double** utild_i, double** utild_j, double** phi, double** rho, int n , int m ,grid ** cells);

#endif /* MATRIXAB_B_H_ */
