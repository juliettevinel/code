/*
 * temp_s_calc.h
 *
 *  Created on: Apr 2, 2014
 *      Author: paant
 */

#ifndef TEMP_S_CALC_H_
#define TEMP_S_CALC_H_

#include <math.h>

#include "trivialfuncs.h"

#include "linsys_solvers.h"

void calc_ts_pred (double** temps1, double** temps0, double** res_ts1, double** res_ts0,
		   double** ksy,  double** phi_s, int i_min, int i_max, int j_min, int j_max, grid cells, double dt) ;

void calc_ts_corr (double** temps1, double** temps0, double** res_ts1, double** res_ts0,
		   double** ksy, double** phi_s, int i_min, int i_max, int j_min, int j_max, grid cells, double dt) ;

void res_ts_compu(double** res_t, double** temp_s, double** temp, double** phi_s, double** kx, double** M_re, double** h, grid cells, int i_min, int i_max, int j_min, int j_max) ;

void temp_s_constr_abc ( double* a, double* b, double* c, int j_min, int j_max, double** k, double** phi_s, grid cells, int i, double dt) ;

void temp_s_constr_rhs_pred (double* r, int j_min, int j_max, double** temp_s, double** res1, double** res0, double** k, double** phi_s, grid cells, int i, double dt) ;

void temp_s_constr_rhs_corr (double* r, int j_min, int j_max, double** temp_s, double** res1, double** res0, double** k, double** phi_s, grid cells, int i, double dt) ;


void ks_compu (double** ks_x, double** ks_y, double** temps, int n, int m) ;


#endif /* TEMP_S_CALC_H_ */
