/*
 * temp_calc.h
 *
 *  Created on: Jan 17, 2014
 *      Author: paant
 */

#ifndef Y_CALC_H_
#define Y_CALC_H_

#include <math.h>

#include "trivialfuncs.h"

#include "linsys_solvers.h"


void calc_Y_pred (double** Y1, double** Y0, double** temps, double** rho, double** rho_s, double** res_Y1, double** res_Y0,
		  double* p0, double M_zero, double** D, double** phi, int n, int m, int pm_res, grid cells, double dt) ;

void calc_Y_corr (double** Y1, double** Y0, double** rho, double** rho_s, double** res_Y1, double** res_Y0,
		  double* p0, double M_zero, double** D, double** phinp1, double** phist, double** phin, int n, int m, int pm_res, grid cells, double dt) ;

void res_Y_compu(double** res_Y, double** f1, double** f2, double** Y, double** rho, double** phi,double** D, double** R, grid cells, int n, int m) ;

void Y_constr_abc ( double* a, double* b, double* c, int m, int pm_res, double** rho, double** D, double** phi, grid cells, int i, double dt) ;

void Y_constr_rhs_pred (double* r, int m, int pm_res, double** Y, double** res1, double** res0, double* p0, double** rho, double** D, double** phi, grid cells, int i, double dt) ;

void Y_constr_rhs_corr (double* r, int m, int pm_res, double** Y, double** res1, double** res0, double* p0, double** rho,double** D, double** phist, double** phin, grid cells, int i, double dt) ;

void h_compu (double** h, double** u, double** v, double** rho, double** mu, double** k, double** phi, int n, int m) ;

void table_vals_cm (double& c, double& m, double Re) ;

void renew_mu_k (double** mu, double** k, double** temp, int init, int n, int m) ;

void M_re_compu (double** M_re, double** temp, double** phi_s, int n, int m) ;

void update_ign (double** ignition, double t) ;

#endif /* Y_CALC_H_ */
