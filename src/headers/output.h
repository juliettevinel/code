/*
 * output.h
 *
 *  Created on: Dec 8, 2010
 *      Author: paant
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <math.h>

#include"global.h"


using namespace std;

void flushdata (double** mat, int n, int m, string file, int t);

void flushbug (double** mat, int init, int n, int m, string file, int t) ;
void flushbug (double* mat, int init, int n, string file, int t);
void flushbug (int* mat, int init, int n, string file, int t);
void flushbug_phi(double** phi, int n, int m);
void flushbug_rij(int n, int m, grid** cells);
void storetime (double time, double* mass, double step,  string file, int t) ;

#endif /* OUTPUT_H_ */
