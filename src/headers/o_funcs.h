/*
 * o_funcs.h
 *
 *  Created on: Jan 24, 2014
 *      Author: paant
 */

#ifndef O_FUNCS_H_
#define O_FUNCS_H_

#include <math.h>

#include "matrixab_c.h"
#include "linsys_solvers.h"


void calc_p (double** p2D, double* p, double* matA, int* matA_pos, double* precond, double* matB,
	     double** utilde_i, double** utilde_j, double** rho1, double** rho0, double** rho_new, double** phi1, double** phi0, double** phi_new, double** M_re, double* y_val, grid cells, int n, int m, double dt, int& flagc, int ti) ;

void calc_mass_eq(double* mass, double** rho1, double** rho0, double** phi1, double** phi0, double** u1, double** v1, double** M_re, int n, int m, grid cells, double dt) ;

void calc_vel (double** u, double** v, double** utilde_i, double** utilde_j, double** p, double** rho, int n, int m, double* y_val, grid cells, double dt) ;

void calc_Uc(double& Uc, double** M_re, double** rho1, double** rho0, double** phi1, double** phi0, double* u_in, int n, int m, grid cells, double dt) ;

void calc_phi( double** phi1, double** phi0, double** M_re, int n, int m, double dt) ;
void calc_phi_s( double** phi_s, double** phi, int init, int n, int m) ;

void filldelta (double** mat, int dir, int n, int m, double** u, double** phi, double** rho, double** mu ) ;

void filldelta_s (double** mat, double** u, double** v, double** rho, double** mu, double** phi, int n , int m) ;

double cd_cyl ( double Rep );
double cd_sph ( double Rep );


double** matrix_n(int n, int m) ;
void delmat_n(double** mat, int n, int m) ;

double** matrix_np(int n, int m) ;
void delmat_np (double** mat, int n, int m) ;

#endif /* O_FUNCS_H_ */
