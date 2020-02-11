/*
 * basicfuncs.h
 *
 *  Created on: Nov 30, 2010
 *      Author: paant
 */

#ifndef BASICFUNCS_H_
#define BASICFUNCS_H_

#include <iostream>
#include <math.h>


#include "trivialfuncs.h"
#include "fluxfuncs.h"
#include "deformfuncs.h"


using namespace std;



void rescompu (double** resu, double** resv, double** u, double** v, double** f1, double** f2, double** phi, double** rho, double** mu,
				double F, int n, int m, grid cells) ;


#endif /* BASICFUNCS_H_ */
