/*
 * utilde_calc.h
 *
 *  Created on: Jan 27, 2014
 *      Author: paant
 */

#ifndef UTILDE_CALC_H_
#define UTILDE_CALC_H_

#include "global.h"
#include "linsys_solvers.h"
#include "trivialfuncs.h"

void calc_utilde_pred (double** utilde_i, double** utilde_j, double** u, double** v, double** resu1, double** resv1, double** resu0, double** resv0,
		       double** p, double** rho1, double** rho0, double** mu, double** phist, double** phin, double** delta_x, double** delta_y, double* y_val, grid cells, int n, int m, double dt, int& flagc) ;

void calc_utilde_corr (double** utilde_i, double** utilde_j, double** u_n, double** v_n, double** resu1, double** resv1, double** resu0, double** resv0,
		       double** p, double** rhonp1, double** rho1, double** mu, double** phinp1, double** phin, double** delta_x, double** delta_y, double* y_val, grid cells, int n, int m, double dt, int& flagc) ;

void ut_constr_abc (double* a, double* b, double* c, int m, double** rho, double** phi, double** mu, double** delta, grid cells, int i, double dt) ;
void ut_constr_rhs_pred ( double* r, int m, double** u, double** res1, double** res0, double** p, double** rho, double** phi, double**  delta, double** mu, grid cells, double y_top, int i, double dt ) ;
void ut_constr_rhs_corr ( double* r, int m, double** u_n, double** res1, double** res0, double** p, double** rho, double** phi, double**  delta, double** mu, grid cells, double y_top, int i, double dt ) ;

void vt_constr_abc (double* a, double* b, double* c, int m, double** rho, double** phi, double** mu, double** delta, grid cells, int i, double dt) ;
void vt_constr_rhs_pred ( double* r, int m, double** v, double** res1, double** res0, double** rho, double** phi, double**  delta, double** mu, grid cells, int i, double dt ) ;
void vt_constr_rhs_corr ( double* r, int m, double** v_n, double** res1, double** res0, double** rho, double** phi, double**  delta, double** mu, grid cells, int i, double dt ) ;

#endif /* UTILDE_CALC_H_ */
