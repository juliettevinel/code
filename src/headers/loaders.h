/*
 * loaders.h
 *
 *  Created on: Jan 16, 2011
 *      Author: panais
 */

#ifndef LOADERS_H_
#define LOADERS_H_



#include <math.h>
#include <fstream>
#include <cstdio>
#include "global.h"

#include "output.h"
#include "o_funcs.h"

using namespace std;


void cellcreate (grid& cells, int n, int m) ;
void celldestroy (grid& cells) ;

void fillval (double** mat, int n, int m, string file);

void fillval_bug (double** mat, int n, int m, string file) ;


void fillu_bound (double* mat, int m, string file);


#endif /* LOADERS_H_ */
