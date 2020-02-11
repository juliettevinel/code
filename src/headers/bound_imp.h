/*
 * bound_imp.h
 *
 *  Created on: Jan 11, 2011
 *      Author: paant
 */

#ifndef BOUND_IMP_H_
#define BOUND_IMP_H_


#include <iostream>
#include <math.h>
#include "global.h"
#include "boundaries.h"
#include "trivialfuncs.h"

using namespace std;


void imp_bound_p(double** mat, double** rho, int n, int m, double* y_val, grid cells);

void imp_bound_u (char type_top, char type_bot, double** mat, double* left_val, double top_val, double bot_val, int n, int m, grid cells) ;

void imp_bound_u_conv (char type_top, char type_bot, double** mat, double** mat0, double* left_val, double top_val, double bot_val, double* Uc, int n, int m, grid cells, double dt) ;

void imp_bound_u_tild (char type_top, char type_bot, double** p, double** rho, double* u_infl, double** mat, int n, int m, double y_top, grid cells, double dt) ;

void imp_bound_u_tild_conv (char type_top, char type_bot, double** p, double** rho, double* u_infl, double** mat, double** mat0, double* Uc, int n, int m, double y_top, grid cells, double dt) ;

void imp_bound_v (char type_top, char type_bot, double** mat, double* left_val, double top_val, double bot_val, int n, int m, grid cells) ;

void imp_bound_v_conv (char type_top, char type_bot, double** mat, double** mat0, double* left_val, double top_val, double bot_val, double* Uc, int n, int m, grid cells, double dt) ;

void imp_bound_v_tild (char type_top, char type_bot, double** mat, int n, int m) ;

void imp_bound_v_tild_conv (char type_top, char type_bot, double* u_infl, double** mat, double** mat0, double* Uc, int n, int m,  grid cells, double dt) ;

void imp_bound_temp (double** temp, double val_left, double val_top, double val_bot, int n, int m) ;

void imp_bound_temp_conv (double** temp, double** temp0, double val_left, double val_top, double val_bot, double *Uc, int n, int m, grid cells, double dt) ;

void imp_bound_temp_s (double** temp, double val_bot, int i_min, int i_max, int j_min, int j_max) ;

void imp_bound_rho (double** rho, double p0, double** temp, int n, int m) ;

void imp_bound_cells ( grid cells, int n, int m);

void imp_bound_mat (double ** mat, int n , int m);

#endif /* BOUND_IMP_H_ */
