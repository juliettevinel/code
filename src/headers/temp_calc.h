/*
 * temp_calc.h
 *
 *  Created on: Jan 17, 2014
 *      Author: paant
 */

#ifndef TEMP_CALC_H_
#define TEMP_CALC_H_

#include <math.h>

#include "trivialfuncs.h"

#include "linsys_solvers.h"


void calc_t_pred (double** temp1, double** temp0, double** temps, double** rho1, double** rho_st, double** res_t1, double** res_t0,
		  double* p0, double M_zero, double** cp, double** kf, double** phi, int n, int m, int pm_res, grid cells, double dt) ;

void calc_t_corr (double** temp1, double** temp0, double** temps, double** rhonp1, double** rho_st, double** res_t1, double** res_t0,
		  double* p0, double M_zero, double** cp, double** kf,  double** phinp1, double** phist, double** phin, int n, int m, int pm_res, grid cells, double dt) ;

void res_t_compu(double** res_t, double** f1, double** f2, double** temp, double** temps, double** rho, double** phi, double** cp, double** M_re,double p0, double** h, double** k, grid cells, int n, int m) ;

void temp_constr_abc ( double* a, double* b, double* c, int m, int pm_res, double** rho, double** cp, double** k, double** phi, grid cells, int i, double dt) ;

void temp_constr_rhs_pred (double* r, int m, int pm_res, double** temp, double** res1, double** res0, double* p0, double** rho, double** cp, double** k, double** phi, grid cells, int i, double dt) ;

void temp_constr_rhs_corr (double* r, int m, int pm_res, double** temp, double** res1, double** res0, double* p0, double** rho,  double** cp, double** k, double** phist, double** phin, grid cells, int i, double dt) ;

void h_compu (double** h, double** u, double** v, double** rho, double** mu, double** k, double** phi, int n, int m) ;

void table_vals_cm (double& c, double& m, double Re) ;

void renew_mu_k (double** mu, double** k, double** temp, int init, int n, int m) ;

void M_re_compu (double** M_re, double** temp, double** phi_s, int n, int m) ;

void update_ign (double** ignition, double t) ;

#endif /* TEMP_CALC_H_ */
