/*
 * o_funcs.cpp
 *
 *  Created on: Jan 24, 2014
 *      Author: paant
 */

#include "./headers/o_funcs.h"


// calculate pressure

void calc_p (double** p2D, double* p, double* matA, int* matA_pos, double* precond, double* matB,
	     double** utilde_i, double** utilde_j, double** rho1, double** rho0, double** rho_new, double** phi1, double** phi0, double** phi_new, double** M_re, double* y_val, grid cells, int n, int m, double dt, int& flagc, int ti) {


	int flag = 0 ;

	double mean = 0.;

	/* ********************************************************* *
	 * CALCULATE THE SOURCE TERM 1D MATRIX OF THE LINEAR SYSTEM* *
	 *********************************************************** */

	if ( (Bcomb_p =='a') || (Bcomb_p == 'c') )
	  constr_b(matB,utilde_i,utilde_j, rho1, rho0, rho_new, phi1, phi0, phi_new, M_re, n,m, y_val, cells, dt) ;
	else if ( Bcomb_p =='b' )
	  constr_b_b(matB,utilde_i,utilde_j, rho1, rho0, rho_new, phi1, phi0, phi_new, M_re, n,m, y_val, cells, dt) ;
	else {
		cout << endl << "ERROR IN calc_p : these boundary conditions for pressure have not been taken into account" << endl;
//		return 0;
	}


	/***************************
	 *  SOLVE THE LINEAR SYSTEM
	 ***************************/


	if (ti > 100)
	  BCSG(matA, matA_pos, p , matB ,precond, 1e-13 ,2500, n*m ,flag);
	else
	  BCSG(matA, matA_pos, p , matB ,precond, 1e-14 ,4000, n*m ,flag);


	if (flag == 1) {
		cout << "ERROR IN THE BCGC SOLVER !" << endl;
		flagc = 1 ;
	}
	/*
	for (int i=1 ; i<=n*m ; i++)
	  mean += p[i] ;

	mean /= n*m ;

	for (int i=1 ; i<=n*m ; i++)
	  p[i] -= mean ;
	*/

	//convert 1d-p to 2d-p
	convert_to_2d(p2D,p,n,m) ;

}



// calculate mass balance

void calc_mass_eq(double* mass, double** rho1, double** rho0, double** phi1, double** phi0, double** u1, double** v1, double** M_re, int n, int m, grid cells, double dt) {

  // 1: volume integral
  // 2: surface integral left
  // 3:   -//-  right
  // 4: surface top
  
  double val_dummy = 0.;

  for (int i=1 ; i<=5 ; i++)
    mass[i] = 0.;
  

  for (int i=1 ; i<=n ; i++)
    for (int j=1 ; j<=m ; j++) {
      
      val_dummy =(phi1[i][j]*rho1[i][j] - phi0[i][j]*rho0[i][j]) /dt  ;
      mass[1] += val_dummy *4. *cells.ri[i] *cells.rj[j] ; 

      mass[2] -= M_re[i][j] *4. *cells.ri[i] *cells.rj[j] ;

    }
  
  
  for (int j=1 ; j<=m ; j++) {

    val_dummy = (phi1[1][j]*rho1[1][j]*u1[1][j] + phi1[0][j]*rho1[0][j]*u1[0][j]) /2. ;
    mass[3] -= val_dummy *2. *cells.rj[j] ;

    val_dummy = (phi1[n][j]*rho1[n][j]*u1[n][j] + phi1[n+1][j]*rho1[n+1][j]*u1[n+1][j]) /2. ;
    mass[4] += val_dummy *2. *cells.rj[j] ;

  }
    
  
  for (int i=1 ; i<=n ; i++) {

    val_dummy = (phi1[i][m]*rho1[i][m]*v1[i][m] + phi1[i][m+1]*rho1[i][m+1]*v1[i][m+1]) /2. ;
    mass[5] += val_dummy *2. *cells.ri[i] ;

  }


}


// calculate valocities

void calc_vel (double** u, double** v, double** utilde_i, double** utilde_j, double** p, double** rho, int n, int m, double* y_val, grid cells, double dt) {

	cell pt;

	// calculate values of u* and store on the u1 matrix
	for (int i=1 ; i<=n; i++) {
		for (int j=1 ; j<=m ; j++) {

			pt.i = i;
			pt.j = j;

			u[i][j] = utilde_i[i][j] - dt / rho[i][j]  * ( derx(p,pt,cells) - Ri * y_val[j] * derx(rho,pt,cells) )  ;

			v[i][j] = utilde_j[i][j] - dt / rho[i][j]  * ( dery(p,pt,cells) - Ri * y_val[j] * dery(rho,pt,cells) )  ;

		}
	}

}


// calculate convective velocity

void calc_Uc(double& Uc, double** M_re, double** rho1, double** rho0, double** phi1, double** phi0, double* u_in, int n, int m, grid cells, double dt) {

  double vol_integ = 0.;
  double surf_integ_in = 0.;
  double surf_integ_out = 0. ;


  double val_dummy ;

  for (int i=1 ; i<=n ; i++)
    for (int j=1 ; j<=m ; j++) {
      val_dummy = M_re[i][j] - (phi1[i][j]*rho1[i][j] - phi0[i][j]*rho0[i][j]) /dt ;
      vol_integ += val_dummy *4. *cells.ri[i] *cells.rj[j] ;
    }
  

  for (int j=1 ; j<=m ; j++) {

    surf_integ_in += u_in[j] *2. *cells.rj[j] ;

    val_dummy = (rho1[n][j] + rho1[n+1][j]) /2. ;
    surf_integ_out += val_dummy *2. *cells.rj[j] ;

  }

  Uc = (vol_integ + surf_integ_in) /surf_integ_out ;
  
}


// calculate phi

void calc_phi( double** phi1, double** phi0, double** M_re, int n, int m, double dt) {

	for (int i=1 ; i<=n ; i++)
		for (int j=1 ; j<=m ; j++)
			phi1[i][j] = phi0[i][j] + dt /rho_s * M_re[i][j]  ;

}



void calc_phi_s( double** phi_s, double** phi, int init, int n, int m) {

	for (int i=init ; i<=n ; i++)
		for (int j=init ; j<=m ; j++)
			phi_s[i][j] = 1. - phi[i][j] ;

}


// fill the matrix with the interphasial drag parameter values
// interphasial drag parameter for cylinders

void filldelta (double** mat, int dir,int n, int m, double** u, double** phi, double** rho, double** mu) {

	double cd = 0. ;  //  cd : drag coefficient //  dp : particle diameter
			// dir for direction : considering that the porous m. is not isotropic we have preference in direction.
			// Values calculated for gereric isotropic porous medium
			// Reference :  Papalexandris "Numerical simulation of detonations in mixtures of gases and solid particles"
			//  JFM 2004

	double u_loc, Rep ;  // scalar of the local velocity , Reynolds of the particle


	if (dir ==1) {


		for (int i=1 ; i <= n ; i++) {
			for (int j=1 ; j<=m ; j++ ) {
	//			u_loc = sqrt ( pow(u[i][j],2.) + pow(v[i][j],2.) ) ;

				u_loc = u[i][j];

				Rep = rho[i][j] * fabs(u_loc) * dp / mu[i][j] * Re ;

				cd = cd_cyl(Rep);

				mat[i][j] = 2./M_PI * cd * rho[i][j] * fabs(u_loc) * (1 - phi[i][j]) / dp ;

				if(mat[i][j] != mat[i][j]) { cout << "delta = " <<  mat[i][j] << "  busted at i = " << i << "\t j= " << j << endl;

					cout << "uloc =" << u_loc << endl;
					cout << "Rep = " << Rep << endl << "cd = " << cd << endl ;

				}

			}
		}
	}  // end if its the x-direction delta matrix



	else if (dir==2) {

		for (int i=1 ; i <= n ; i++) {
			for (int j=1 ; j<=m ; j++ ) {


			  //				mat[i][j] = 2.656 * (1. - phi[j]) / dp  * sqrt( mu[i][j] * rho[i][j] * fabs( u[i][j] ) / p_height ) / sqrt (Re) ;
//for debugging purposes, simulations with no p.m. or non-standard domain configuration
			  mat[i][j] = 2.656 * (1. - phi[i][j]) / dp  * sqrt( mu[i][j] * rho[i][j] * fabs( u[i][j] ) ) / sqrt (Re) ;
				// attention! : sqrt includes height of the porous medium which is l_ref and therefore 1

			}
		}

	}//  endof else, aka, the y-direction delta matrix contruction.

	else cout << endl << "direction for delta not expected, file 'rescalc.cpp' " << endl;

}



// interphasial drag parameter for spheres

void filldelta_s (double** mat, double** u, double** v, double** rho, double** mu, double** phi, int n, int m) {

  double cd = 0.;
  double u_loc,Rep;

  for (int j=1 ; j<=m ; j++ ) {
    for (int i=1 ; i <= n ; i++) {

      u_loc = sqrt (u[i][j] * u[i][j] + v[i][j] * v[i][j]) ;

      Rep = rho[i][j] * u_loc * dp / mu[i][j] * Re;

      cd = cd_sph(Rep) ;

      mat[i][j] = 3./4. * cd * (1 - phi[i][j]) * rho[i][j] * u_loc / dp ;

    }
  }

}


double cd_cyl ( double Rep ) {

	double cd;

	if (Rep < .1)
		cd = 8. * M_PI / (Rep + 1e-16) / (1./2. - 0.577215665 - log( (Rep + 1e-16) /8.)) ;
	else if ((Rep >= .1) && (Rep < 4.))
		cd = pow(Rep, -0.726705485) * pow(10., 1.039580289) ;
	else if ((Rep >= 4.) && (Rep < 40.))
		cd = pow(Rep, -0.301029996) * pow(10., 0.783298108) ;
	else if ((Rep >= 40.) && (Rep < 600.))
		cd = pow(Rep, -0.255958025) * pow(10., 0.711090107) ;
	else if ((Rep >= 600.))
		cd = 1.;
	else {
		cout << endl << "drag coefficient value not expected" << endl;
		cd = 1./0.;
	}

	return cd;

}



double cd_sph ( double Rep ) {

		double cd;
		if (Rep < 0.1)
			cd = 24. / (Rep + 1e-16) ;
		else if ((0.1 <= Rep) && (Rep < 1))
			cd = pow(Rep, -0.982271233) * pow(10., 1.397940009) ;
	    else if ((1. <= Rep) && (Rep < 20.))
	    	cd = pow(Rep, -0.707761356) * pow(10., 1.397940009) ;
	    else if ((20. <= Rep) && (Rep < 1000.))
	    	cd = pow(Rep, -0.51505398) * pow(10., 1.147221933) ;
	    else if (1000. <= Rep)
	    	cd = 0.4 ;
		else {
			cout << endl << "drag coefficient value not expected" << endl;
			cd = 1./0.;
		}

	    return cd;
}







//
//
//  FUNCTIONS THAT CONSTRUCT AND DESTROY ARRAYS
//
//




double** matrix_n(int n, int m) {

double **mat;

mat = new double*[n];
for (int i=0 ; i<n ; i++) {
	mat[i] = new double[m];
	mat[i]--;
}
mat--;

return mat;

}


void delmat_n(double** mat, int n, int m) {

mat++;

for (int i=0 ; i<n ; i++) {
	mat[i]++;
	delete [] mat[i] ;
}
delete [] mat ;

}



double** matrix_np(int n, int m) {

	double **mat;

	mat = new double*[n];
	for (int i=0 ; i<n ; i++)
		mat[i] = new double[m] ;

	return mat;
}

void delmat_np (double** mat, int n, int m) {

	for (int i=0 ; i<n ; i++)
		delete [] mat[i] ;

	delete [] mat ;

}
