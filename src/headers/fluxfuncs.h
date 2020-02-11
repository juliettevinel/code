/*
 * fluxfuncs.h
 *
 *  Created on: Nov 30, 2010
 *      Author: paant
 */

#ifndef FLUXFUNCS_H_
#define FLUXFUNCS_H_


#include "global.h"
#include "trivialfuncs.h"
#include <math.h>
#include <iostream>

using namespace std;


void convcompu (vec& output,double** u, double** v, double** f1, double** f2, cell pt,  grid cells) ;

void fluxcompu (double** f1, double** f2, double** uti, double**vti, double** rho, double** phi, double** p, double* y_val, grid cells, int n, int m, double dt);

void fluxinit (double** f1, double**f2, double** u, double** v, double** rho, double** phi, int n, int m, grid cells) ;


#endif /* FLUXFUNCS_H_ */
