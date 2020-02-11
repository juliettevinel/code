
#include "./headers/temp_s_calc.h"




void calc_ts_pred (double** temps1, double** temps0, double** res_ts1, double** res_ts0,
		   double** ksy,  double** phi_s, int i_min, int i_max, int j_min, int j_max, grid cells, double dt) {



	double *mat_a, *mat_b, *mat_c, *rhs, *sol, *gam ;

	int pm_res = j_max - j_min + 1 ;

	mat_a = new double[pm_res];  mat_a--;
	mat_b = new double[pm_res];  mat_b--;
	mat_c = new double[pm_res];  mat_c--;
	rhs = new double[pm_res]; rhs--;
	sol = new double[pm_res]; sol--;
	gam = new double[pm_res]; gam--;

	// flag detecting errors in tri_dag
	int flag = 0;



	//------------------
	// start
	//------------------


	for (int i=i_min ; i<=i_max ; i++) {

	  temp_s_constr_abc (mat_a, mat_b, mat_c, j_min, j_max, ksy, phi_s, cells, i, dt) ;
	  temp_s_constr_rhs_pred ( rhs, j_min, j_max, temps0, res_ts1, res_ts0, ksy, phi_s, cells, i , dt) ;

	  tridag (mat_a, mat_b, mat_c, rhs, sol, pm_res, gam, flag, i) ;

	 
	  for (int j=j_min ; j<=j_max ; j++)
	    temps1[i][j] = sol[j - j_min +1] ;

	}



	if (flag == 1)
	  cout << "ERROR IN calc_ts_pred" << endl;



	mat_a++ ; delete [] mat_a;
	mat_b++ ; delete [] mat_b;
	mat_c++ ; delete [] mat_c;
	rhs++ ; delete [] rhs ;
	sol++ ; delete [] sol ;
	gam++ ; delete [] gam;



	//boundary conditions for temperature solid.


	/*

	for (int i=1 ; i<=n ; i++)
	  for (int j=1 ; j<=pmedia_res ; j++) 
	    temps1[i][j] = temps0[i][j]
	      + dt /cp_s /(1.-phi[j]) /rho_s * (3./2. *res_ts1[i][j] - 1./2. *res_ts0[i][j] - h[i][j] * (tempf0[i][j] - temps0[i][j]) ) ;

	for (int i=1 ; i<=n ; i++)
	  for (int j=pmedia_res+1 ; j<=m ; j++) 
	    temps1[i][j] = (kf[i][j] + kf[i][j-1])/2. * (tempf1[i][j] - tempf1[i][j-1]) /(ks[i][j] + ks[i][j-1]) *2. + temps1[i][j-1] ;
	    
	*/
	    

}



void calc_ts_corr (double** temps1, double** temps0, double** res_ts1, double** res_ts0,
		   double** ksy, double** phi_s, int i_min, int i_max, int j_min, int j_max, grid cells, double dt) {



	double *mat_a, *mat_b, *mat_c, *rhs, *sol, *gam ;

	int pm_res = j_max - j_min + 1 ;

	mat_a = new double[pm_res];  mat_a--;
	mat_b = new double[pm_res];  mat_b--;
	mat_c = new double[pm_res];  mat_c--;
	rhs = new double[pm_res]; rhs--;
	sol = new double[pm_res]; sol--;
	gam = new double[pm_res]; gam--;

	// flag detecting errors in tri_dag
	int flag = 0;



	for (int i=i_min ; i<=i_max ; i++) {

	  temp_s_constr_abc (mat_a, mat_b, mat_c, j_min, j_max, ksy, phi_s, cells, i, dt) ;
	  temp_s_constr_rhs_corr ( rhs, j_min, j_max, temps0, res_ts1, res_ts0, ksy, phi_s, cells, i , dt) ;

	  tridag (mat_a, mat_b, mat_c, rhs, sol, pm_res, gam, flag, i) ;

	  for (int j=j_min ; j<=j_max ; j++)
	    temps1[i][j] = sol[j - j_min +1] ;

	}



	if (flag == 1)
	  cout << "ERROR IN calc_ts_pred" << endl;



	mat_a++ ; delete [] mat_a;
	mat_b++ ; delete [] mat_b;
	mat_c++ ; delete [] mat_c;
	rhs++ ; delete [] rhs ;
	sol++ ; delete [] sol ;
	gam++ ; delete [] gam;


	//boundary conditions for temperature solid.


	/*
	// calculate rho* T*
	for (int i=1 ; i<=n ; i++)
	  for (int j=1 ; j<=pmedia_res ; j++) 
	    temps1[i][j] = temps0[i][j]
	      + dt /cp_s /(1.-phi[j]) /rho_s * (1./2. *res_ts1[i][j] + 1./2. *res_ts0[i][j] - h[i][j] * (tempf0[i][j] - temps0[i][j]) ) ;

	for (int i=1 ; i<=n ; i++)
	  for (int j=pmedia_res +1 ; j<=m ; j++) 
	    	    temps1[i][j] = (kf[i][j] + kf[i][j-1])/2. * (tempf1[i][j] - tempf1[i][j-1]) /(ks[i][j] + ks[i][j-1]) *2. + temps1[i][j-1] ;
	    
	*/

	  

}







void res_ts_compu(double** res_t, double** temp_s, double** temp, double** phi_s, double** kx, double** M_re, double** h, grid cells, int i_min, int i_max, int j_min, int j_max) {

  double diffa  ;

  double diffb ;

  cell pt ;

  for (int i=i_min ; i<=i_max ; i++) {
    for (int j=j_min ; j<=j_max ; j++) {

    	pt.i = i;
    	pt.j = j;



      	  // diffusive term

      diffa = valfor(phi_s,kx, pt,cells) * derxfor(temp_s, pt,cells) /2. /cells.ri[i] ;

      diffb = valback(phi_s,kx, pt,cells) * derxback(temp_s, pt,cells) /2. /cells.ri[i] ;

      
      /*    THIS PART IS INTEGRATED IN TIME IMPLICITLY, SO IT'S LEFT OUT OF THE RESIDUE
      diffb = (ky[i][j+1]*(1.-phi[j+1]) + ky[i][j]*(1.-phi[j])) /2. * derytop(temp_s, pt, cells)
	-  (ky[i][j-1]*(1.-phi[j-1]) + ky[i][j]*(1.-phi[j])) /2. * derybot(temp_s, pt, cells) ;
      //diffb = 0. ;

      diffb *= 1. /2. /cells[i][j].rj ;
      */

      res_t[i][j] = 1./Re /Pr * (diffa - diffb) - h[i][j] * (temp[i][j] - temp_s[i][j]) + M_re[i][j] * cp_s * (temp_s[i][j] - 1.) ;

    }
  }

}


// The 'k' in this function is the ky, since the term treated implicitly is the one with derivatives along the y
void temp_s_constr_abc ( double* a, double* b, double* c, int j_min, int j_max, double** k, double** phi_s, grid cells, int i, double dt) {

  double alpha, betta;

  cell pt ;

  pt.i = i ;

  for (int j=j_min ; j<=j_max ; j++) {

	  pt.j = j ;

    if (j == j_max) {

      alpha = 0. ;
      betta = 0.5 /Re /Pr * dt * valbot(phi_s,k, pt,cells) /2. /cells.rj[j] /(cells.rj[j] + cells.rj[j-1]) ;
    
    }

    else if ( j == j_min ) {

      alpha = 0.5 /Re /Pr * dt * valtop(phi_s,k, pt,cells) /2. /cells.rj[j] /(cells.rj[j] + cells.rj[j+1]) ;
      betta = 0.;

    } 

    else {

      alpha = 0.5 /Re /Pr * dt * valtop(phi_s,k, pt,cells) /2. /cells.rj[j] /(cells.rj[j] + cells.rj[j+1]) ;
      betta = 0.5 /Re /Pr * dt * valbot(phi_s,k, pt,cells) /2. /cells.rj[j] /(cells.rj[j] + cells.rj[j-1]) ;

    }

    a[j-j_min+1] = - betta ;

    c[j-j_min+1] = - alpha ;

    b[j-j_min+1] = cp_s * rho_s * phi_s[i][j] + alpha + betta ;

  }

  /*
  if (Bcomb_Tt == 'o'){
	  // lower boundary condition
	  b[1] += a[1];
  }
  else if (Bcomb_Tt == 'd') {
	  // lower boundary condition
	  b[1] -= a[1];
  }
  else
  	  cout << endl << "ERROR in the choice of bcomb_Tt, in file temp_s_calc.cpp" << endl;
  */
}



void temp_s_constr_rhs_pred (double* r, int j_min, int j_max, double** temp_s, double** res1, double** res0, double** k, double** phi_s, grid cells, int i, double dt) {

  double term_n, term_s ;

  cell pt ;
  pt.i = i ;

  for (int j=j_min ; j<=j_max ; j++) {

	pt.j = j ;

    if (j == j_max) {
      term_n = 0. ;
      term_s = valbot(phi_s,k, pt,cells) * derybot(temp_s, pt,cells) /2. /cells.rj[j] ;
    }   
    else if (j==j_min){
      term_n = valtop(phi_s,k, pt,cells) * derytop(temp_s, pt,cells) /2. /cells.rj[j] ;
      term_s = 0. ;
    }
    else {
      term_n = valtop(phi_s,k, pt,cells) * derytop(temp_s, pt,cells) /2. /cells.rj[j] ;
      term_s = valbot(phi_s,k, pt,cells) * derybot(temp_s, pt,cells) /2. /cells.rj[j] ;
    }

    r[j-j_min +1] = cp_s * rho_s * phi_s[i][j] * temp_s[i][j] + 1.5*dt*res1[i][j] - 0.5*dt*res0[i][j]
      + dt*0.5 /Re /Pr * (term_n - term_s) ;

  }

  /*
  if (Bcomb_Tt == 'd'){
	  pt.j = 1;
	  r[1] += 2. * Tval_b * 0.5 *dt /Re /Pr * valbot(phi_s,k, pt,cells) /2. /cells.rj[1] /(cells.rj[1] + cells.rj[0]) ;
  }
  */

}



void temp_s_constr_rhs_corr (double* r, int j_min, int j_max, double** temp_s, double** res1, double** res0, double** k, double** phi_s, grid cells, int i, double dt) {

  double term_n, term_s ;

  cell pt ;
  pt.i = i ;

  for (int j=j_min ; j<=j_max ; j++) {

	pt.j = j ;

    if (j == j_max) {
      term_n = 0. ;
      term_s = valbot(phi_s,k, pt,cells) * derybot(temp_s, pt,cells) /2. /cells.rj[j] ;
    }   
    else if (j ==j_min) {
      term_n = valtop(phi_s,k, pt,cells) * derytop(temp_s, pt,cells) /2. /cells.rj[j] ;
      term_s =  0.;
    }
    else {
      term_n = valtop(phi_s,k, pt,cells) * derytop(temp_s, pt,cells) /2. /cells.rj[j] ;
      term_s = valbot(phi_s,k, pt,cells) * derybot(temp_s, pt,cells) /2. /cells.rj[j] ;
    }

    r[j - j_min +1] = cp_s * rho_s * phi_s[i][j] * temp_s[i][j] + 0.5*dt*res1[i][j] + 0.5*dt*res0[i][j]
      + dt*0.5 /Re /Pr * (term_n - term_s) ;

  }

  /*
  if (Bcomb_Tt == 'd'){
	  pt.j = 1;
	  r[1] += 2. * Tval_b * 0.5 *dt /Re /Pr * valbot(phi_s,k, pt,cells) /2. /cells.rj[1] /(cells.rj[1] + cells.rj[0]) ;
  }
  */

}




void ks_compu (double** ks_x, double** ks_y, double** temps, int n, int m) {

  // b_hat (dimensional from experimental data) = 0.00037

  double betta = 0.00037* T_ref /k_ref ;
  double alpha = 0.0103345 /k_ref ;

  for (int i=0 ; i<= n+1 ; i++)
    for (int j=0 ; j<=m+1 ; j++) {

      ks_x[i][j] = alpha + betta * temps[i][j] ;
      
      ks_y[i][j] = 1.8 * ks_x[i][j] ;

    }
  

}
