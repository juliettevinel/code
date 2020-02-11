
#include "./headers/Y_calc.h"


void calc_Y_pred (double** Y1, double** Y0, double** rho, double** rho_s, double** res_Y1, double** res_Y0,
		  double* p0, double M_zero, double** D,  double** phi, int n, int m, int pm_res, grid cells, double dt) {




	double *mat_a, *mat_b, *mat_c, *rhs, *sol, *gam ;

	mat_a = new double[m];  mat_a--;
	mat_b = new double[m];  mat_b--;
	mat_c = new double[m];  mat_c--;
	rhs = new double[m]; rhs--;
	sol = new double[m]; sol--;
	gam = new double[m]; gam--;

	// flag detecting errors in tri_dag
	int flag = 0;

	//------------------
	// start
	//------------------

    

	for (int i=1 ; i<=n ; i++) {

	  Y_constr_abc (mat_a, mat_b, mat_c, m, pm_res, rho, D, phi, cells, i, dt) ;
	  Y_constr_rhs_pred ( rhs, m, pm_res, Y0, res_Y1, res_Y0, p0, rho,D, phi, cells, i , dt) ;

	  tridag (mat_a, mat_b, mat_c, rhs, sol, m, gam, flag, i) ;

	  for (int j=1 ; j<=m ; j++)
	    Y1[i][j] = sol[j] ;

	}


	if (flag == 1)
	  cout << "ERROR IN calc_t_pred" << endl;

// POUR CLOSED DOMAIN JE COMPRENDS PAS
	 if ( open_d == 0 ) {
		// integration dummy holding the cumulative sum
		double integ_d = 0.;

		// calculate a first estimate for t*
		for (int i=1 ; i<=n ; i++)
			for (int j=1 ; j<=m ; j++) // integ (1/T)dv
				integ_d += phi[i][j] * (2.*cells.ri[i]) * (2.*cells.rj[j]) / temp1[i][j] ;


		// we will be replacing p0[1], so we move the value to the obsolete p0[0]
		//	p0[0] = p0[1] ;

		// a first estimate for the p0^*
		p0[2] = M_zero / integ_d ;
		//	p0[2] = 1. ;
	}

	//for (int i=1 ; i<=n ; i++)
	  //for (int j=1 ; j<=m ; j++)
	   // rho_st[i][j] = p0[2] /temp1[i][j] ;


	mat_a++ ; delete [] mat_a;
	mat_b++ ; delete [] mat_b;
	mat_c++ ; delete [] mat_c;
	rhs++ ; delete [] rhs ;
	sol++ ; delete [] sol ;
	gam++ ; delete [] gam;


	//boundary conditions for temperature.

}




void calc_Y_corr (double** Y1, double** Y0, double** rho, double** rho_s, double** res_Y1, double** res_Y0,
		  double* p0, double M_zero, double** D, double** phinp1, double** phist, double** phin, int n, int m, int pm_res, grid cells, double dt) {



	double *mat_a, *mat_b, *mat_c, *rhs, *sol, *gam ;

	mat_a = new double[m];  mat_a--;
	mat_b = new double[m];  mat_b--;
	mat_c = new double[m];  mat_c--;
	rhs = new double[m]; rhs--;
	sol = new double[m]; sol--;
	gam = new double[m]; gam--;

	// flag detecting errors in tri_dag
	int flag = 0;


	//------------------
	// start
	//------------------


	for (int i=1 ; i<=n ; i++) {

	  Y_constr_abc (mat_a, mat_b, mat_c, m, pm_res, rho_s,D, phist, cells, i, dt) ;
	  temp_constr_rhs_corr ( rhs, m, pm_res, Y0, res_Y1, res_Y0, p0, rho_s,D, phist, phin, cells, i , dt) ;

	  tridag (mat_a, mat_b, mat_c, rhs, sol, m, gam, flag, i) ;

	  for (int j=1 ; j<=m ; j++)
	    Y1[i][j] = sol[j] ;

	}


	if (flag == 1)
	  cout << "ERROR IN calc_t_corr" << endl;



	if (open_d == 0) {
		// integration dummy holding the cumulative sum
		double integ_d = 0.;

        // JE COMPRENDS PAS CE QUE C'EST
		// calculate a first estimate for t*
		for (int i=1 ; i<=n ; i++)
			for (int j=1 ; j<=m ; j++) // integ (1/T)dv
				integ_d += phinp1[i][j] * (2.*cells.ri[i]) * (2.*cells.rj[j]) / temp1[i][j] ;

		// This is the corrector stage.
		// Therefore, we don't put the (*) values at the place of the (n) values.
		// p0[0] = p0[1] ;

		// a first estimate for the new p0
		p0[2] = M_zero / integ_d ;
		//	p0[2] = 1. ;
	}


	// Now we are about to calculate the rho(n+1) values.
	
	//for (int i=1 ; i<=n ; i++)
	//  for (int j=1 ; j<=m ; j++)
	  //  rhonp1[i][j] = p0[2] /temp1[i][j] ;
	

	//boundary conditions for temperature.


	mat_a++ ; delete [] mat_a;
	mat_b++ ; delete [] mat_b;
	mat_c++ ; delete [] mat_c;
	rhs++ ; delete [] rhs ;
	sol++ ; delete [] sol ;
	gam++ ; delete [] gam;


}




void res_Y_compu(double** res_Y, double** f1, double** f2, double** Y, double** rho, double** phi,double** D, double** R, grid cells, int n, int m) {

  double conva, convb, conv_term ;
  double diffa, diffb,  diff_term ;

  //  double diffb ;

  cell pt ;

  for (int i=1 ; i<=n ; i++) {
    for (int j=1 ; j<=m ; j++) {

    	pt.i = i;
    	pt.j = j;

    	// convective term

	conva = ( f1[i][j] * derxfor(Y, pt, cells) + f1[i-1][j] * derxback(Y, pt, cells) ) /2. ;

	convb = ( f2[i][j] * derytop(Y, pt, cells) + f2[i][j-1] * derybot(Y, pt, cells) ) /2. ;

	conv_term = -(conva + convb) ;


      	  // diffusive term

      diffa = valfor(phi,D, pt,cells) * derxfor(Y, pt, cells) /2. /cells.ri[i]  ;

      diffb = valback(phi,D, pt,cells) * derxback(Y, pt,cells) /2. /cells.ri[i] ;

      
      diff_term = rho[i][j]*(diffa - diffb) ;
      
      
        res_Y[i][j] = conv_term + diff_term + rho[i][j]*phi[i][j] * Y[i][j] * R[i][j];

    }
  }

}



void Y_constr_abc ( double* a, double* b, double* c, int m, int pm_res, double** rho, double** D, double** phi, grid cells, int i, double dt) {

  double alpha, betta;

  cell pt;

  pt.i = i;

  for (int j=1 ; j<=m ; j++) {

	pt.j = j ;

	/*
    if (j == pm_res) {
    	alpha = 0.5 /Re /Pr * dt * valtop(k, pt,cells) /2. /cells.rj[j] /(cells.rj[j] + cells.rj[j+1]) ;

    	betta = 0.5 /Re /Pr * dt * valbot(phi,k, pt,cells) /2. /cells.rj[j] /(cells.rj[j] + cells.rj[j-1]) ;
     }
    else if (j == pm_res+1) {
    	alpha = 0.5 /Re /Pr * dt * valtop(phi,k, pt,cells) /2. /cells.rj[j] /(cells.rj[j] + cells.rj[j+1]) ;

    	betta = 0.5 /Re /Pr * dt * valbot(k, pt,cells) /2. /cells.rj[j] /(cells.rj[j] + cells.rj[j-1]) ;
    }
	*/
	//    else {
    	alpha = 0.5 * dt * valtop(phi,D, pt,cells) /2. /cells.rj[j] /(cells.rj[j] + cells.rj[j+1]) ;

    	betta = 0.5 * dt * valbot(phi,D, pt,cells) /2. /cells.rj[j] /(cells.rj[j] + cells.rj[j-1]) ;
	//    }

    a[j] = - betta*rho[i][j] ;

    c[j] = - alpha *rho[i][j];

    b[j] = phi[i][j] *rho[i][j] + alpha + betta ;

  }

  if (Bcomb_Yt == 'o'){
	  // lower boundary
	  b[1] += a[1];
	  // upper boundary
	  b[m] += c[m];
  }
  else if (Bcomb_Yt == 'd') {
	  // lower boundary
	  b[1] -= a[1];
	  // upper boundary
	  b[m] -= c[m];
  }
  else
	  cout << endl << "ERROR in the choice of bcomb_Tt, in file Y_calc.cpp" << endl;

}



void Y_constr_rhs_pred (double* r, int m, int pm_res, double** Y, double** res1, double** res0, double* p0, double** rho, double** D, double** phi, grid cells, int i, double dt) {

  double term_n, term_s ;

  cell pt;

  pt.i = i ;

  for (int j=1 ; j<=m ; j++) {

	  pt.j = j;

	  /*
    if (j == pm_res) {
    	term_n = valtop(k, pt,cells) * derytop(temp, pt,cells) /2. /cells.rj[j] ;

    	term_s = valbot(phi,k, pt,cells) * derybot(temp, pt,cells) /2. /cells.rj[j] ;
    }   
    else if (j == pm_res+1) {
    	term_n = valtop(phi,k, pt,cells) * derytop(temp, pt,cells) /2. /cells.rj[j] ;

    	term_s = valbot(phi,k, pt,cells) * derybot(temp, pt,cells) /2. /cells.rj[j] ;
    }
	  */
	  //    else {
      term_n = valtop(phi,D, pt,cells) * derytop(Y, pt,cells) /2. /cells.rj[j] ;

      term_s = valbot(D, pt,cells) * derybot(Y, pt,cells) /2. /cells.rj[j] ;
      //    }

    r[j] =  rho[i][j]*phi[i][j] * Y[i][j]  + 1.5*dt*res1[i][j] - 0.5*dt*res0[i][j] + dt*0.5 * rho[i][j]*(term_n - term_s) ;

  }

  
// C'est quoi ce truc ????
  if (Bcomb_Yt == 'd') {

	pt.j = 1;
    r[1] += 2. * Yval_b * 0.5 *dt * valbot(phi,D, pt,cells)*rho[i][j] /2. /cells.rj[1] /(cells.rj[1] + cells.rj[0]) ;

    pt.j = m ;
    r[m] += 2. * Yval_t * 0.5 *dt * valtop(phi,D, pt,cells)*rho[i][j] /2. /cells.rj[m] /(cells.rj[m] + cells.rj[m+1]);

  }

}



void Y_constr_rhs_corr (double* r, int m, int pm_res, double** Y, double** res1, double** res0, double* p0, double** rho, double** D, double** phist, double** phin, grid cells, int i, double dt) {

  double term_n, term_s ;

  cell pt;

  pt.i = i ;

  for (int j=1 ; j<=m ; j++) {

	  pt.j = j;

	  /*
    if (j == pm_res) {
    	term_n = valtop(k, pt,cells) * derytop(temp, pt,cells) /2. /cells.rj[j] ;

    	term_s = valbot(phist,k, pt,cells) * derybot(temp, pt,cells) /2. /cells.rj[j] ;
    }   
    else if (j == pm_res+1) {
    	term_n = valtop(phist,k, pt,cells) * derytop(temp, pt,cells) /2. /cells.rj[j] ;

    	term_s = valbot(phist,k, pt,cells) * derybot(temp, pt,cells) /2. /cells.rj[j] ;
    }
	  */
	  //    else {
      term_n = valtop(phist,D, pt,cells) * derytop(Y, pt,cells) /2. /cells.rj[j] ;

      term_s = valbot(D, pt,cells) * derybot(Y, pt,cells) /2. /cells.rj[j] ;
      //    }

    r[j] = rho[i][j] * phist[i][j] * Y[i][j] + 0.5*dt*res1[i][j] + 0.5*dt*res0[i][j]
      + dt*0.5 * rho[i][j]*(term_n - term_s) ;

  }

  
// c'est quoi ce truc
  if (Bcomb_Yt == 'd') {

	pt.j = 1;
    r[1] += 2. * Yval_b * 0.5 *dt * valbot(phist,D, pt,cells)*rho[i][j] /2. /cells.rj[1] /(cells.rj[1] + cells.rj[0]) ;

    pt.j = m ;
    r[m] += 2. * Yval_t * 0.5 *dt *valtop(phist,D, pt,cells)*rho[i][j] /2. /cells.rj[m] /(cells.rj[m] + cells.rj[m+1]);

  }

}





void h_compu (double** h, double** u, double** v, double** rho, double** mu, double** k, double** phi, int n, int m) {


	double Nu_dp, Re_dp, C, mm ;


	for (int i=1 ; i<=n ; i++)
		for (int j=1 ; j<=m ; j++) {

			Re_dp = rho[i][j] * sqrt(u[i][j]*u[i][j] + v[i][j]*v[i][j]) * dp / mu[i][j] * Re ;

			table_vals_cm (C, mm, Re_dp) ;


			Nu_dp = C * pow(Re_dp, mm) * pow (Pr, 1./3.) ;

			h[i][j] = -4. * (1 - phi[i][j]) /dp /dp * Nu_dp * pow(Pr, -2./3.) /Re * k[i][j] ;

		}

}


void table_vals_cm (double& c, double& m, double Re) {

	if (Re < 4.) {
		c = 0.989 ; m = 0.33 ;
	}
	else if ((Re >= 4. ) && (Re < 40.)) {
		c = 0.911 ; m = 0.385 ;
	}
	else if ((Re >=40. ) && (Re < 4000.)) {
		c = 0.683 ; m = 0.466 ;
	}
	else if ((Re >= 4000.) && (Re < 40000.)) {
		c = 0.193 ; m = 0.618 ;
	}
	else if (Re >= 40000) {
		c = 0.027 ; m = 0.805 ;
	}

}

void renew_mu_k (double** mu, double** k, double** temp, int init, int n, int m) {

	for (int i=init ; i<=n ; i++)
		for (int j=init ; j<=m ; j++) {

		  // mu[i][j] = pow(temp[i][j], 0.7) ;
			mu[i][j] = 0.10607 + 1.23749 * temp[i][j] - 0.571147 * temp[i][j] * temp[i][j]
			          +0.322822 * temp[i][j] * temp[i][j] * temp[i][j]
			          -0.112499 * temp[i][j] * temp[i][j] * temp[i][j] * temp[i][j]
					  +0.0172672* temp[i][j] * temp[i][j] * temp[i][j] * temp[i][j] * temp[i][j] ;
		  //		  mu[i][j] = 1. ;

			k[i][j] = mu[i][j] ;

		}

}


void M_re_compu (double** M_re, double** temp, double** phi_s, int n, int m) {

	for (int i=1 ; i<=n ; i++)
		for (int j=1 ; j<=m ; j++)
			M_re[i][j] = K_reaction * exp(-e_activation /temp[i][j]) *4. * phi_s[i][j] /dp ;


}


void update_ign (double** ignition, double t) {

  double val, val_max ;

  val_max = exp(1.51);

  val = t>1. ? val_max : exp(t*1.51) ;

  if (t<8.) {

    for (int i = 250 ; i<= 265 ; i++)
      for (int j = 70 ; j<=85 ; j++)
	ignition[i][j] = ignition[i][j] > val_max ? ignition[i][j] : val ;

  }


}

