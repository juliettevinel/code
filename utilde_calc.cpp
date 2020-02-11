/*
 * utilde_calc.cpp
 *
 *  Created on: Jan 27, 2014
 *      Author: paant
 */

#include "./headers/utilde_calc.h"
#include "./headers/output.h"


void calc_utilde_pred (double** utilde_i, double** utilde_j, double** u, double** v, double** resu1, double** resv1, double** resu0, double** resv0,
		       double** p, double** rho_st, double** rho_n, double** mu, double** phist, double** phin, double** delta_x, double** delta_y, double* y_val, grid cells, int n, int m, double dt, int& flagc) {

	double *matUt_a, *matUt_b, *matUt_c, *rhs, *sol, *gam ;

	matUt_a = new double[m];  matUt_a--;
	matUt_b = new double[m];  matUt_b--;
	matUt_c = new double[m];  matUt_c--;
	rhs = new double[m]; rhs--;
	sol = new double[m]; sol--;
	gam = new double[m]; gam--;

	// flag detecting errors in tri_dag
	int flag = 0;

	double y_top;
	y_top = y_val[m] + cells.rj[m] ;


	for (int i=1 ; i<=n ; i++) {

	    // for u_tilde
	    ut_constr_abc ( matUt_a, matUt_b, matUt_c, m, rho_st, phist, mu, delta_x, cells, i, dt) ;
	    ut_constr_rhs_pred( rhs, m, u, resu1, resu0, p, rho_n, phin, delta_x, mu, cells, y_top, i,dt);
	    tridag(matUt_a, matUt_b, matUt_c, rhs, sol, m, gam, flag, i);

	    for (int j=1 ; j<=m ; j++)
	      utilde_i[i][j] = sol[j] ;

	    // for v_tilde
	    vt_constr_abc ( matUt_a, matUt_b, matUt_c, m, rho_st, phist, mu, delta_y, cells, i,dt) ;
	    vt_constr_rhs_pred( rhs, m, v, resv1, resv0, rho_n, phin, delta_y, mu, cells, i, dt);
	    tridag(matUt_a, matUt_b, matUt_c, rhs, sol, m, gam, flag, i);

	   
	    
	    for (int j=1 ; j<=m ; j++)
	      utilde_j[i][j] = sol[j] ;

	}

	if (flag ==1)
		flagc = 1;

	matUt_a++ ; delete [] matUt_a;
	matUt_b++ ; delete [] matUt_b;
	matUt_c++ ; delete [] matUt_c;
	rhs++ ; delete [] rhs ;
	sol++ ; delete [] sol ;
	gam++ ; delete [] gam;

}



void calc_utilde_corr (double** utilde_i, double** utilde_j, double** u_n, double** v_n, double** resu1, double** resv1, double** resu0, double** resv0,
		       double** p, double** rhonp1, double** rho1, double** mu, double** phinp1, double** phin, double** delta_x, double** delta_y, double* y_val, grid cells, int n, int m, double dt, int& flagc) {

	double *matUt_a, *matUt_b, *matUt_c, *rhs, *sol, *gam ;

	matUt_a = new double[m];  matUt_a--;
	matUt_b = new double[m];  matUt_b--;
	matUt_c = new double[m];  matUt_c--;
	rhs = new double[m]; rhs--;
	sol = new double[m]; sol--;
	gam = new double[m]; gam--;

	// flag detecting errors in tri_dag
	int flag = 0;

	double y_top;
	y_top = y_val[m] + cells.rj[m] ;


	for (int i=1 ; i<=n ; i++) {

	    // for u_tilde
	    ut_constr_abc ( matUt_a, matUt_b, matUt_c, m, rhonp1, phinp1, mu, delta_x, cells, i, dt) ;
	    ut_constr_rhs_corr( rhs, m, u_n, resu1, resu0, p, rho1, phin, delta_x, mu, cells, y_top, i,dt);
	    tridag(matUt_a, matUt_b, matUt_c, rhs, sol, m, gam, flag, i);

	    for (int j=1 ; j<=m ; j++)
	      utilde_i[i][j] = sol[j] ;

	    // for v_tilde
	    vt_constr_abc ( matUt_a, matUt_b, matUt_c, m, rhonp1, phinp1, mu, delta_y, cells, i,dt) ;
	    vt_constr_rhs_corr( rhs, m, v_n, resv1, resv0, rho1, phin, delta_y, mu, cells, i, dt);
	    tridag(matUt_a, matUt_b, matUt_c, rhs, sol, m, gam, flag, i);

	    for (int j=1 ; j<=m ; j++)
	      utilde_j[i][j] = sol[j] ;

	}

	if (flag ==1)
		flagc = 1;

	matUt_a++ ; delete [] matUt_a;
	matUt_b++ ; delete [] matUt_b;
	matUt_c++ ; delete [] matUt_c;
	rhs++ ; delete [] rhs ;
	sol++ ; delete [] sol ;
	gam++ ; delete [] gam;

}





void ut_constr_abc (double* a, double* b, double* c, int m, double** rho, double** phi, double** mu, double** delta, grid cells, int i, double dt) {

  double alpha, betta;

  cell pt ;

  pt.i = i ;

  for (int j=1; j<=m ; j++ ) {

	  pt.j = j ;

    alpha =  0.5 *dt /Re * valtop(phi,mu, pt,cells) /2. /cells.rj[j] /( cells.rj[j] + cells.rj[j+1] ) ;

    betta =  0.5 *dt /Re * valbot(phi,mu, pt,cells) /2. /cells.rj[j] /( cells.rj[j] + cells.rj[j-1] ) ;


    a[j] = - betta;

    c[j] = - alpha ;

    b[j] = phi[i][j]*rho[i][j] + alpha + betta + 0.5*dt *delta[i][j] ;   // 1./2./Re is included in alpha and betta
    
  }


  // Lower boundary condition, zero Dirichlet for u ( add 0.5*betta )
  if (Bcomb_ub == 'w')
	  b[1] -= a[1] ;
  else if (Bcomb_ub == 'o')  // zero Neumann
	  b[1] += a[1] ;
  else
	  cout << "Error in bottom b-cond. implementation for u_tilde." << endl;


  // Upper boundary condition, zero Dirichlet for u
  if (Bcomb_ut == 'w')
	  b[m] -= c[m];
  else if (Bcomb_ut == 'o') // zero Neumann for u
	  b[m] += c[m];
  else
	  cout << "Error in top b-cond. implementation for u_tilde." << endl;


}  // endof function constr_a




void ut_constr_rhs_pred ( double* r, int m, double** u, double** res1, double** res0, double** p, double** rho, double** phi, double**  delta, double** mu, grid cells, double y_top, int i, double dt ) {

  double term_n, term_s ;

  double dummy ;

  cell pt ;

  pt.i = i ;

  for (int j=1 ; j<=m ; j++) {

	  pt.j = j ;

    term_n = valtop(phi,mu, pt,cells) * derytop(u, pt,cells) /2. /cells.rj[j] ;

    term_s = valbot(phi,mu, pt,cells) * derybot(u, pt,cells) /2. /cells.rj[j] ;


    r[j] = rho[i][j] * phi[i][j] * u[i][j] + 1.5 *dt *res1[i][j] - 0.5 *dt *res0[i][j] - 0.5 *dt *delta[i][j] * u[i][j]
      + 0.5 *dt /Re *( term_n - term_s ) ;
  }

  /*
  if (Bcomb_ub == 'w')
    {
      pt.j = 1 ;

      dummy = - 0.5 *dt /Re * valbot(phi,mu, pt,cells) /2. /cells.rj[1] /( cells.rj[1] + cells.rj[0] ) ;

      r[1] -= 2. *dt /valbot(rho,pt,cells) * derxbot(p, pt,cells) * dummy ;
    }
  
  if (Bcomb_ut == 'w')
    {
      pt.j = m;

      dummy = - 0.5 *dt /Re * valtop(phi,mu, pt,cells) /2. /cells.rj[m] /( cells.rj[m] + cells.rj[m+1] ) ;

      r[m] -= 2. *dt /valtop(rho,pt,cells) * (derxtop(p,pt,cells) - Ri * y_top * derxtop(rho,pt,cells) ) * dummy ;
    }
  */
}



void ut_constr_rhs_corr ( double* r, int m, double** u, double** res1, double** res0, double** p, double** rho, double** phi, double**  delta, double** mu, grid cells, double y_top, int i, double dt ) {

  double term_n, term_s ;

  double dummy ;

  cell pt ;

  pt.i = i;

  for (int j=1 ; j<=m ; j++) {

	pt.j = j ;

    term_n = valtop(phi,mu, pt,cells) * derytop(u, pt,cells) /2. /cells.rj[j] ;

    term_s = valbot(phi,mu, pt,cells) * derybot(u, pt,cells) /2. /cells.rj[j] ;

    
    r[j] = rho[i][j] * phi[i][j] * u[i][j] + 0.5 *dt *res1[i][j] + 0.5 *dt *res0[i][j] - 0.5 *dt *delta[i][j] * u[i][j]
      + 0.5 *dt /Re *( term_n - term_s ) ;

  }

  /*
  if  (Bcomb_ub == 'w')
    {
      pt.j = 1 ;

      dummy = - 0.5 *dt /Re * valbot(phi,mu, pt,cells) /2. /cells.rj[1] /( cells.rj[1] + cells.rj[0] ) ;

      r[1] -= 2. *dt /valbot(rho,pt,cells) * derxbot(p, pt,cells) * dummy ;
    }
  
  if (Bcomb_ut == 'w')
    {
      pt.j = m;

      dummy = - 0.5 *dt /Re * valtop(phi,mu, pt,cells) /2. /cells.rj[m] /( cells.rj[m] + cells.rj[m+1] ) ;

      r[m] -= 2. *dt /valtop(rho,pt,cells) * (derxtop(p,pt,cells) - Ri * y_top * derxtop(rho,pt,cells) ) * dummy ;
    }
  */
}




void vt_constr_abc (double* a, double* b, double* c, int m, double** rho, double** phi, double** mu, double** delta, grid cells, int i, double dt) {

  double alpha, betta;

  cell pt;

  pt.i = i ;

  for (int j=1; j<=m ; j++ ) {

	pt.j = j ;

    alpha =  2. /3. *dt /Re * valtop(phi,mu, pt,cells) /2. /cells.rj[j] /( cells.rj[j] + cells.rj[j+1] ) ;

    betta =  2. /3. *dt /Re * valbot(phi,mu, pt,cells) /2. /cells.rj[j] /( cells.rj[j] + cells.rj[j-1] ) ;


    a[j] = - betta;

    c[j] = - alpha ;

    b[j] = phi[i][j]*rho[i][j]  + alpha + betta + 0.5*dt* delta[i][j] ;   // 2./3./Re is included in alpha and betta
  }


  // Lower boundary condition, zero dirichlet for v ( add 0.5*betta )
  if (Bcomb_vb == 'w')
    b[1] -= a[1] ;
  else if (Bcomb_vb == 'o') // zero Neumann
    b[1] += a[1] ;
  else
    cout << "Error in bottom b-cond. implementation for v_tilde." << endl;


  // Upper boundary condition, zero Dirichlet for v
  if (Bcomb_vt == 'w')
    b[m] -= c[m];
  else if (Bcomb_vt == 'o') // zero Neumann for v
    b[m] += c[m];
  else
    cout << "Error in top b-cond. implementation for v_tilde." << endl;

}  // endof function constr_a




void vt_constr_rhs_pred ( double* r, int m, double** v, double** res1, double** res0, double** rho, double** phi, double**  delta, double** mu, grid cells, int i, double dt ) {

  double term_n, term_s ;

  cell pt ;

  pt.i = i;

  for (int j=1 ; j<=m ; j++) {

	  pt.j = j ;

    term_n = valtop(phi,mu, pt,cells) * derytop(v, pt,cells) /2. /cells.rj[j] ;

    term_s = valbot(phi,mu, pt,cells) * derybot(v, pt,cells) /2. /cells.rj[j] ;


    r[j] = rho[i][j] * phi[i][j] * v[i][j] + 1.5 *dt *res1[i][j] - 0.5 *dt *res0[i][j] - 0.5 *dt *delta[i][j] * v[i][j]
      + 2./3. *dt /Re * ( term_n - term_s )  ;

  }

}


void vt_constr_rhs_corr ( double* r, int m, double** v, double** res1, double** res0, double** rho, double** phi, double**  delta, double** mu, grid cells, int i, double dt ) {

  double term_n, term_s ;

  cell pt ;

  pt.i = i ;

  for (int j=1 ; j<=m ; j++) {

	  pt.j = j ;

    term_n = valtop(phi,mu, pt,cells) * derytop(v, pt,cells) /2. /cells.rj[j] ;

    term_s = valbot(phi,mu, pt,cells) * derybot(v, pt,cells) /2. /cells.rj[j] ;


    r[j] = rho[i][j] * phi[i][j] * v[i][j] + 0.5 *dt *res1[i][j] + 0.5 *dt *res0[i][j] - 0.5 *dt *delta[i][j] * v[i][j]
      + 2./3. *dt /Re * ( term_n - term_s )  ;

  }

}


