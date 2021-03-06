/*
 * main.cpp
 *
 *  Created on: Nov 29, 2010
 *      Author: Juliette
 */

#include "./headers/LFOPM.h"


int main () {

	cout.precision(16);


	//-------GRID RESOLUTION---------

	int xres, yres ;
	double dimx, dimy ;

	int pmedia_jmin, pmedia_jmax, pmedia_imin, pmedia_imax ;

	ifstream griddata;
	griddata.open("../data/grid.dat");  if (! griddata.is_open()) cout << "input file error, grid.data" ;

	griddata >> xres >> yres >> dimx >> dimy  >> pmedia_jmin >> pmedia_jmax >> pmedia_imin >> pmedia_imax ;
	griddata.close();

	//resolution of grid n x m
	const int n = xres;   // length
	const int m = yres;	// height

	double F = 0.;

	double Rep = dp * Re ;
	double cd;


	if ((forcing == 1)&&(delta_rank == 'c')) {
		cd = cd_cyl(Rep);
		F = - 2. /M_PI * cd * (1. - vf_fluid) / vf_fluid / dp;
	}
	// complete expression : gradP = Cd * rho * (1 - phi) * u^2 / phi / dp
	else if ((forcing == 1)&&(delta_rank == 's')) {
		cd = cd_sph(Rep);
		F = - 3./4. * cd * (1 - vf_fluid) / vf_fluid / dp ;
	}
	else if (forcing==0)
		F = 0.;
	else
		cout << "invalid value in F.dat" << endl;

	// test-case, laminar flow in between infinite parallel plates, with distance 2
	//F = -2. /Re;
	// ... with distance 3
	//F = -8. /9. /Re ;

	//F = -1./2./Re ;

	ifstream tdata;
	tdata.open("../data/time.dat");  if (! tdata.is_open()) cout << "input file error, time.dat" ;

	double t, dt;
	int step_init;

	tdata >> t >> step_init >> dt ;
	tdata.close();

// to read things from tstep_data  (converts chars to int)
	ifstream tstep_data;
	tstep_data.open("../data/tsteps.dat");  if (! tstep_data.is_open()) cout << "input file error, tsteps.dat" ;

	int tsteps, ti_min, sample_out ;

	tstep_data >> tsteps >> ti_min >> sample_out ;
	tstep_data.close();



	double uval_t, vval_t, uval_b, vval_b;

	uval_t = readdouble("../data/uval_t.dat");
	vval_t = readdouble("../data/vval_t.dat");

	uval_b = readdouble("../data/uval_b.dat");
	vval_b = readdouble("../data/vval_b.dat");

// same logic as in the previous use of <<  and  >> 
	cout << "RUN for " << tsteps << " time steps, dt = " << dt << endl;

	cout << endl << "Re no :" <<  Re << "\t Pr no : " << Pr << "\t Ri no: " << Ri << endl;

	cout << endl << "Domain Resolution and Dimensions :" << endl ;
	cout << "n=" << n << "\t m=" << m << "\t \t Pmedia resolution : " << pmedia_jmax << endl;
	cout << "x * y = " << dimx << "x" << dimy << "\t \t pmedia height :"
			<< endl;

	cout << endl << "pressure gradient forcing F=" << F << endl;


	//************
	//DEFINITIONS
	//*************

	//CELLS matrix holding data of grid dimensions and cell distribution

	grid cells ;

	double *y_val;
	y_val = new double[m+2] ;

	double **phinm1;
	phinm1 = matrix_np(n+2,m+2) ;

	double **phin;
	phin = matrix_np(n+2,m+2) ;

	double **phist;
	phist = matrix_np(n+2,m+2) ;

	double **phinp1;
	phinp1 = matrix_np(n+2,m+2) ;


	double **phi_sn;
	phi_sn = matrix_np(n+2,m+2) ;

	double **phi_sst;
	phi_sst = matrix_np(n+2,m+2) ;

	double **phi_snp1;
	phi_snp1 = matrix_np(n+2,m+2) ;

	double **delta_x ;
	delta_x = matrix_n(n,m) ;
	double **delta_y ;
	delta_y = matrix_n(n,m) ;

	double **mu ;
	mu = matrix_np(n+2, m+2) ;

	initialize(mu, 0, n+1, m+1) ;


	// U AND V FOR  t = -1 and t = 0   (matrices u0 and u1 respectively)
	//Dynamically construct u and v two dimensional matrices destined to hold velocity values
	//attention : the size of these matrices holds as well the values of the ghost cells.

	double **u0 ;
	u0 = matrix_np(n+2,m+2);
	double **v0 ;
	v0 = matrix_np(n+2,m+2);

	double **u1 ;
	u1 = matrix_np(n+2,m+2);
	double **v1 ;
	v1 = matrix_np(n+2,m+2);

	double **f1;
	f1 = new double*[n+1];
	for (int i=0 ; i<n+1 ; i++){
		f1[i] = new double[m];
		f1[i] --;
	}

	double **f2;
	f2 = new double*[n];
	for (int i=0 ; i<n ; i++)
		f2[i] = new double[m+1];
	f2--;

	double **resu0 ;
	resu0 = matrix_n(n,m) ;
	double **resv0 ;
	resv0 = matrix_n(n,m) ;

	double **resu1 ;
	resu1 = matrix_n(n,m) ;
	double **resv1 ;
	resv1 = matrix_n(n,m) ;


	// U and V matrices for the inflow profiles

	double * u_infl, * v_infl;
	u_infl = new double[m];
	v_infl = new double[m];

	// preserve consistency with other u matrices : internal cells 1->m (1 is the index of the 1st interior cell)
	u_infl--;
	v_infl--;



	//MATRIX P//

	double *p ;
		p = new double[n*m] ;
		p--;

	double **p2D;
	p2D = matrix_np(n+2,m+2) ;

	double** mb2d ;
	mb2d = matrix_n(n,m) ;

	double p0[3] ;

	double **utilde_i ;
	utilde_i = matrix_np(n+2,m+2);

	double **utilde_i0 ;
	utilde_i0 = matrix_np(n+2,m+2);

	double **utilde_j ;
	utilde_j = matrix_np(n+2,m+2);


	double **utilde_j0 ;
	utilde_j0 = matrix_np(n+2,m+2);



	double **rho0;
	rho0 = matrix_np(n+2,m+2);
	double **rho1;
	rho1 = matrix_np(n+2,m+2);
	double **rhonp1;
	rhonp1 = matrix_np(n+2,m+2) ;
	double **rho_st;
	rho_st = matrix_np(n+2,m+2);

	double **temp0 ;
	temp0 = matrix_np(n+2,m+2);
	double **temp1 ;
	temp1 = matrix_np(n+2,m+2);

	initialize(temp1, 0, n+1, m+1) ;

	double **temp_s0 ;
	temp_s0 = matrix_np(n+2,m+2);

	double **temp_s1 ;
	temp_s1 = matrix_np(n+2,m+2);

	double **res_ts0;
	res_ts0 = matrix_n(n,m);
	double **res_ts1;
	res_ts1 = matrix_n(n,m);

	initialize (res_ts0, 1,n,m) ;
	initialize (res_ts1, 1,n,m) ;

	double **res_t0;
	res_t0 = matrix_n(n,m);
	double **res_t1;
	res_t1 = matrix_n(n,m);

	double **h;
	h = matrix_n(n,m) ;

	double **kf ;
	kf = matrix_np(n+2,m+2) ;

	double **ks_x ;
	ks_x = matrix_np(n+2,m+2) ;

	double **ks_y ;
	ks_y = matrix_np(n+2,m+2) ;

	initialize(kf, 0, n+1, m+1) ;

//	double **ks1 ;

//	double **ks2 ;

	double **cp ;
	cp = matrix_n(n,m) ;

	double **M_re;
	M_re = matrix_n(n,m) ;


	double *mass_balance ;
	mass_balance = new double[5] ;
	mass_balance--;


	/* A AND B MATRICES OF THE LINEAR SYSTEM****/

	//sparce storage of matrix A, consists of 1D matrix A, + the position matrix

		int sizeofmatA ;  // the  size of the A matrix (sparce storaged) depends on the boundary conditions for p

		if (Bcomb_p == 'c')
			sizeofmatA = 5*n*m - 2*n + 1;
		else if ((Bcomb_p == 'a')||(Bcomb_p == 'b'))
			sizeofmatA = 5*n*m - 2*(n+m) + 1 ;
		else {
			cout << endl << "Error in the choice of boundary conditions for pressure !!!" << endl;
			return 0;
		}



	double *matA ;
		matA = new double[sizeofmatA] ;
		matA-- ;
	int *matA_pos;
		matA_pos = new int[sizeofmatA] ;
		matA_pos--;

	double *precond ;        // precondition matrix - necessity for the solver-function of the linear system
		precond = new double [n*m] ;
		precond-- ;

	double *matB ;
		matB = new double[n*m] ;
		matB --;







		/************
		 * Construct the domain
		 ************/

		cellcreate(cells,n,m);
		imp_bound_cells(cells,n,m);

		y_val[0] = -cells.rj[0] ;
		y_val[1] = cells.rj[1] ;
		for (int j=2 ; j<=m+1 ; j++) y_val[j] = y_val[j-1] + cells.rj[j] + cells.rj[j-1] ;

		fillval(phin,n,m,"../init/centphi.txt");
		imp_bound_mat(phin, n,m);

		replace(phist, phin, 0, n+1, m+1) ;
		replace(phinp1, phin, 0, n+1, m+1) ;
		replace(phinm1, phin, 0, n+1, m+1) ;

		calc_phi_s( phi_sn, phin, 0, n+1, m+1) ;
		replace(phi_sst, phi_sn, 0, n+1, m+1) ;
		replace(phi_snp1, phi_sn, 0, n+1, m+1) ;


		/*******************************
		 * FILLING IN MATRIX A AND POS *
		 *******************************/





		//meet dvector equivalence (equivalence with linear system solver)


/*****************************************
 *         LOAD INITIAL DATA
 ******************************************/

	for (int i=1 ; i<=n ; i++)
		for (int j=1 ; j<=m ; j++)
			cp[i][j] = 1.;


	fillval(temp1, n,m, "../init/centtemp.txt") ;
	imp_bound_temp(temp1, 1., Tval_t, Tval_b, n,m) ;

	replace(temp0, temp1, 0, n+1, m+1) ;
	


	fillval(temp_s1, n,m, "../init/centtemp_s.txt") ;
	imp_bound_temp_s(temp_s1, Tval_b, pmedia_imin, pmedia_imax, pmedia_jmin,pmedia_jmax) ;

	replace(temp_s0, temp_s1, 0, n+1, m+1) ;


	ks_compu (ks_x, ks_y, temp_s0, n,m) ;
	


	// integration dummy holding the cumulative sum
	double integ_d = 0.;

	// integ (1/T)dv
	for (int i=1 ; i<=n ; i++)
		for (int j=1 ; j<=m ; j++)
			integ_d += phin[i][j] * (2.*cells.ri[i]) * (2.*cells.rj[j]) / temp1[i][j] ;


	double M_zero = 0.;
	double M_zero_test ;
	double M_zero_o = 0.;

	p0[0] = 1.;
	p0[1] = 1.;
	p0[2] = 1.;

	
	for( int i=1 ; i<=n ; i++)
	  for (int j=1 ; j<=m ; j++)
	    rho1[i][j] = p0[1] / temp1[i][j] ;
	    


	//	fillval(rho1, n,m, "../init/centrho.txt") ;

	imp_bound_rho(rho1, p0[1], temp1, n,m) ;

	replace(rho0,rho1, 0, n+1, m+1) ;

	//	fillval(rho0, n,m, "../init/centrhom1.txt") ;
	//imp_bound_rho(rho0, p0[1], temp1, n,m) ;


	replace(rho_st, rho1, 0, n+1, m+1) ;
	replace(rhonp1, rho1, 0, n+1, m+1) ;



	  for (int i=1 ; i<=n ; i++)
	    for (int j=1 ; j<=m ; j++)
	      M_zero += rho1[i][j] * phin[i][j] * 4. *cells.ri[i] *cells.rj[j] ;
	

	 M_re_compu(M_re, temp1, phi_sn, n,m) ;




	renew_mu_k(mu, kf, temp1, 0, n+1, m+1) ;





	//vels for t = n-1
	// fillval_bug(u0,n,m,"../init/centu_db.txt");
	// fillval_bug(v0,n,m,"../init/centv_db.txt");

	fillval(u0,n,m,"../init/centu.txt");
	fillval(v0,n,m,"../init/centv.txt");


	//fill matrices for u, v inflow
	//	fillu_bound(u_infl,m,"../init/u_dir.txt");
	//	fillu_bound(v_infl,m,"../init/v_dir.txt");

	for (int j=1 ; j<=m ; j++) {
	  u_infl[j] = u0[1][j];
	  v_infl[j] = v0[1][j] ;
	}



	// THE FOLLOWING COMMENTS APPLY PERTURBATIONS ON THE VELOCITIES

/*  
	srand(194585);
	double dumrand1, dumrand2 ;


	  for (int i=1 ; i<=n ; i++){
	    for (int j=1 ; j<=m ; j++) {
	      dumrand1 = rand()%100000 / 100000. - 0.5 ; dumrand1 *= 1e-9 ;
	      dumrand2 = rand()%100000 / 100000. - 0.5 ; dumrand2 *= 1e-9 ;
	    	  u0[i][j] += dumrand1 ;
	    	  v0[i][j] += dumrand2 ;

	    }
	  }
*/



	/*
	double pert_amp = 1e-9 ;
	double gauss_steep = 20.;
	double num_per1 = 2.;
	//	double num_per2 = 2.;


	for (int i=1 ; i<=n ; i++) {
	  for (int j=1 ; j<=m ; j++) {
	    
	    	    u0[i][j] +=  pert_amp * exp ( - gauss_steep * pow ( (2.*j+1)*cells[i][j].rj - 1., 2. ) ) * cos( 2.*M_PI / dimx * (2*i+1)*cells[i][j].ri * num_per1) ;
	    //    u0[i][j] +=  phi[j] * pert_amp * sin( 2.*M_PI / dimx * (2*i-1)*cells[i][j].ri * num_per1) *  cos( 2.*M_PI / dimy * (2*j-1)*cells[i][j].rj * num_per1)  ;
	    	    v0[i][j] +=  pert_amp * exp ( - gauss_steep * pow ( (2.*j+1)*cells[i][j].rj - 1., 2. ) ) * sin( 2.*M_PI / dimx * (2*i+1)*cells[i][j].ri * num_per1) ;
	    //v0[i][j] -=  phi[j] * pert_amp * cos( 2.*M_PI / dimx * (2*i-1)*cells[i][j].ri * num_per1) *  sin( 2.*M_PI / dimy * (2*j-1)*cells[i][j].rj * num_per1)  ;


	  }
	}
	*/
	


	imp_bound_u(Bcomb_ut , Bcomb_ub , u0 ,u_infl, uval_t, uval_b, n , m ,cells );
	imp_bound_v(Bcomb_vt , Bcomb_vb , v0 ,v_infl, vval_t, vval_b, n , m ,cells );

	replace (u1, u0, 0, n+1, m+1) ;
	replace (v1, v0, 0, n+1, m+1) ;

	
	double phiUc ;

	calc_Uc (phiUc, M_re, rhonp1, rho1, phinp1, phin, u_infl, n,m, cells, dt) ;
	
	double* Uc;
	Uc = new double[m] ;
	Uc-- ;

	for (int j=1 ; j<=m ; j++)
	  Uc[j] = phiUc /phinp1[n][j] ;
	  


	if (delta_rank == 'c') {
		filldelta(delta_x,1,n,m,u1,phin, rho1,mu);
		filldelta(delta_y,2,n,m,v1,phin, rho1,mu);
	}
	else if (delta_rank == 's') {
	  filldelta_s(delta_x,u1,v1,rho_st,mu,phin,n,m);
	  filldelta_s(delta_y,u1,v1,rho_st,mu,phin,n,m);
	}




	fluxinit(f1, f2, u1, v1, rho1, phin, n, m, cells) ;

	rescompu(resu1,resv1, u1,v1, f1,f2, phin, rho1, mu, F, n,m,cells) ;

	replace(resu0, resu1, 1, n,m ) ;
	replace(resv0, resv1, 1, n,m ) ;

	h_compu(h,u1, v1, rho1, mu, kf, phin, n,m );

	res_t_compu(res_t1, f1,f2, temp1, temp_s1, rho1, phin, cp, M_re, p0[1], h, kf, cells, n,m) ;
	replace (res_t0, res_t1, 1, n,m) ;

	
	ks_compu (ks_x, ks_y, temp_s1, n,m) ;
	res_ts_compu (res_ts1, temp_s1, temp1, phi_sn, ks_x, M_re, h, cells, pmedia_imin, pmedia_imax, pmedia_jmin , pmedia_jmax ) ;

	

	replace (res_ts0, res_ts1) ;
		      
		     

	//fill matrix p for solver efficiency
	//fillval_bug(p2D, n,m, "../init/centp_db.txt") ;
	fillval(p2D, n,m, "../init/centp.txt") ;

	imp_bound_p(p2D , rho1, n , m , y_val, cells) ;

	convert_to_1d(p, p2D, n,m) ;

	int flag = 0 ; // It takes the value 1, when the linear solver doesn't converge, leading to the termination of the execution


	// This is needed only so that we can apply the convective boundary condition later on. 
	// Therefore, only the utilde_i is of interest
	  calc_utilde_pred(utilde_i,utilde_j, u1,v1, resu1,resv1,resu0,resv0, p2D, rho_st,rho1,mu,phist, phin, delta_x, delta_y, y_val, cells,n,m,dt, flag ) ;
	  // boundary conditions for utilde
	  imp_bound_u_tild(Bcomb_ut, Bcomb_ub, p2D, rho1, u_infl, utilde_i, n,m, dimy, cells, dt) ;
	  imp_bound_v_tild(Bcomb_vt, Bcomb_vb, utilde_j, n,m) ;

	  replace(utilde_i0, utilde_i, 0, n+1, m+1) ;
	  replace(utilde_j0, utilde_j, 0, n+1, m+1) ;
/*****************************************
 *         LOADING INITIAL DATA - END
 ******************************************/




    //**************************************
	/***************************************
	 * BEGINNING OF THE ITERATING PROCEDURE*
	 **************************************/
	//**************************************




	// ADAPTIVE DT VARS
	double umax, vmax ; // auxiliary var. for adaptive dt routine.
	double dt_cand_u, dt_cand_v ;  // candidate for next dt value;


	double temp_inflow = 1.; 


	for (int ti = 1 + step_init ; ti <= tsteps ; ti++) {  // ti for t_i, integer number of step
	  
	  if (t<4.)
	    temp_inflow = 1.5/4. *t + 1.;
	  else
	    temp_inflow = 2.5 ;
	  


		//calculate phi*
		calc_phi (phist, phin, M_re, n,m, dt) ;
		imp_bound_mat(phist,n,m) ;

		calc_phi_s (phi_sst, phist, 0, n+1, m+1) ;

		//temp0 -> (n)
		//calculate T*, T_s*
	  replace(temp0, temp1) ;
	  calc_t_pred(temp1,temp0, temp_s1, rho1,rho_st, res_t1,res_t0, p0,M_zero,cp, kf, phin, n,m, pmedia_jmax, cells, dt);

	  //	  update_ign(temp1, t) ;

	  imp_bound_temp_conv(temp1, temp0, temp_inflow, Tval_t, Tval_b, Uc, n,m, cells, dt) ;
	  imp_bound_rho(rho_st, p0[2], temp1, n,m) ;


	  renew_mu_k(mu, kf, temp1, 0, n+1, m+1) ;

	  
	  replace(temp_s0,temp_s1) ;
	  calc_ts_pred(temp_s1, temp_s0, res_ts1, res_ts0, ks_y, phi_sn, pmedia_imin, pmedia_imax, pmedia_jmin , pmedia_jmax , cells, dt) ;

	  //	  update_ign(temp_s1, t) ;

	  imp_bound_temp_s(temp_s1, Tval_b, pmedia_imin, pmedia_imax, pmedia_jmin,pmedia_jmax) ;
	  

	  ks_compu(ks_x, ks_y, temp_s1, n,m) ;


	  replace(utilde_i0, utilde_i, 0, n+1, m+1) ;
	  replace(utilde_j0, utilde_j, 0, n+1, m+1) ;

	  // ATTENTION : UTILDE CALCULATION CONFIGURED FOR WALL ON TOP
	  calc_utilde_pred(utilde_i,utilde_j, u1,v1, resu1,resv1,resu0,resv0, p2D, rho_st,rho1,mu,phist, phin, delta_x, delta_y, y_val, cells,n,m,dt, flag ) ;
	  // boundary conditions for utilde
	  imp_bound_u_tild_conv(Bcomb_ut, Bcomb_ub, p2D, rho1, u_infl, utilde_i, utilde_i0, Uc, n,m, dimy, cells, dt) ;
	  imp_bound_v_tild_conv(Bcomb_vt, Bcomb_vb, v_infl, utilde_j, utilde_j0, Uc, n,m, cells, dt) ;


	  if (Bcomb_p == 'c')
	  		constr_a_c (matA, matA_pos, precond, phist, n, m, cells) ;
	  else if (Bcomb_p == 'a')
			constr_a (matA, matA_pos, precond, phist, n, m, cells) ;
	  else if (Bcomb_p == 'b')
			constr_a_b (matA, matA_pos, precond, phist, n, m, cells) ;
	  else {
		  cout << endl << "ERROR IN calc_p : these boundary conditions for pressure have not been taken into account" << endl;
	  		return 0;
	  }


	  calc_p (p2D,p, matA,matA_pos,precond,matB, utilde_i, utilde_j, rho1,rho0, rho_st, phin, phinm1, phist, M_re, y_val, cells, n,m,dt, flag, ti) ;
	  // boundary conditions for pressure
	  imp_bound_p(p2D ,rho_st, n , m , y_val, cells) ;


	  /*
	  //  ELLIPTIC EQUATION SOLVER DEBUGGING PACKAGE, predictor
	    if (( ti >= ti_min)&&( ti%sample_out == 0 )) {  // output
	    //flushdata(res_t1, n,m, "rtpr",ti) ;
	    convert_to_2d(mb2d, matB, n,m) ;
	    flushdata( mb2d, n,m, "mbpr", ti) ;
	    //constr_b_testa (testB, utilde_i, utilde_j, rho1, rho0, rho_st, phi, n,m,cells, dt) ;
	    //dimensionalize (mb2d, testB, n,m) ;
	    //flushdata (mb2d, n,m, "tbpr", ti) ;
	    //flushdata (rho_st,n,m, "rst", ti) ;
	    //constr_b_testb (testB, utilde_i, utilde_j, rho1, rho0, rho_st, phi, n,m,cells, dt) ;
	    //dimensionalize (mb2d, testB, n,m) ;
	    //flushdata (mb2d, n,m, "tb2pr", ti) ;
	    }	   
	  */


	  // the u(-1) values considered obsolete and switched with u(0)
	  // u0 -> (n)
	  // u1 -> (*)
	  replace(u0,u1) ;
	  replace(v0,v1) ;

	  calc_vel (u1,v1, utilde_i, utilde_j, p2D, rho_st, n,m, y_val, cells, dt) ;
	  //impose boundary conditions to u*

	  if (Bcomb_vel == 2)
	    calc_Uc (phiUc, M_re, rho_st, rho1, phist, phin, u_infl, n,m, cells, dt) ;

	  for (int j=1 ; j<=m ; j++)
	    Uc[j] = phiUc /phist[n][j] ;


	  imp_bound_u_conv(Bcomb_ut, Bcomb_ub, u1 ,u0, u_infl, uval_t, uval_b, Uc, n , m , cells, dt );
	  imp_bound_v_conv(Bcomb_vt, Bcomb_vb, v1 ,v0, v_infl, vval_t, vval_b, Uc, n , m , cells, dt );


	  //
	  // PREDICTOR STEP COMPLETED
	  // Preparing all (*) values
	  //


	  fluxcompu(f1, f2, utilde_i, utilde_j, rho_st, phist, p2D, y_val, cells, n, m, dt) ;


	  if (delta_rank == 'c') {
	    filldelta(delta_x,1,n,m,u1,phist, rho_st,mu);
	    filldelta(delta_y,2,n,m,v1,phist, rho_st,mu);
	  }
	  else if (delta_rank == 's') {
	    filldelta_s(delta_x,u1,v1,rho_st,mu,phist,n,m);
	    filldelta_s(delta_y,u1,v1,rho_st,mu,phist,n,m);
	  }





	  // res0 -> (n)
	  // calculate (*)
	  replace(resu0,resu1);
	  replace(resv0,resv1);

	  rescompu(resu1, resv1, u1,v1, f1,f2, phist, rho_st, mu, F, n,m,cells) ;

	  h_compu(h,u1, v1, rho_st, mu, kf, phist, n,m );

	  M_re_compu(M_re, temp1, phi_sst, n,m) ;

	  // the same for res_t
	  replace(res_t0, res_t1) ;

	  res_t_compu(res_t1, f1,f2, temp1, temp_s1, rho_st, phist, cp, M_re, p0[2], h, kf, cells, n,m) ;

	  // the same for res_ts
	  replace(res_ts0, res_ts1) ;
	  
	  res_ts_compu(res_ts1, temp_s1, temp1, phi_sst, ks_x, M_re, h, cells, pmedia_imin, pmedia_imax, pmedia_jmin , pmedia_jmax ) ;


	  //ENDOF ALL U* PROCEDURES



	  // predictor output
	  /*
	  if (( ti >= ti_min)&&( ti%sample_out == 0 )) {  // output
	    flushbug (u1,0,n+1,m+1,"upr",ti);
	    flushbug (v1,0,n+1,m+1,"vpr",ti);
	    flushbug (p2D,0,n+1,m+1,"ppr",ti);

	    flushbug (temp1, 0,n+1,m+1,"tpr",ti) ;
	    flushbug (temp_s1,0,n+1,m+1,"tspr", ti) ;

	    convert_to_2d(mb2d, matB, n,m) ;
	    flushdata(mb2d, n,m, "mbpr", ti) ;

	    flushbug(rho_st, 0,n+1,m+1,"rpr",ti) ;

	    flushdata (resu1, n,m, "ru",ti) ;
	    flushdata (resv1, n,m, "rv", ti);

	    flushbug (utilde_i, 0,n+1, m+1, "utipr", ti) ;
	    flushbug (utilde_j, 0,n+1, m+1, "utjpr", ti) ;

	    flushbug (f1[10], 1, m, "f1", ti) ;
	    flushbug (f2[10], 0, m, "f2", ti) ;

	    flushdata (res_ts1, n,pmedia_jmax, "rtsp", ti) ;
	  }
	  */



	  //****************************************************
	  //*************************************************
	  /***************************************************
	   * 			CORRECTOR STAGE
	   ***************************************************/
	  //***************************************************




	  //calculate phi(n+1)
	  calc_phi(phinp1,phin, M_re, n,m,dt) ;
	  imp_bound_mat(phinp1, n,m) ;

	  calc_phi_s(phi_snp1,phinp1, 0, n+1, m+1) ;


	  //temp0 -> (n)
	  //temp1 -> (n+1)
	  calc_t_corr(temp1,temp0, temp_s0, rhonp1, rho_st, res_t1,res_t0, p0,M_zero, cp, kf, phinp1, phist, phin, n,m, pmedia_jmax, cells, dt);

	  //	  update_ign(temp1, t) ;

	  imp_bound_temp_conv(temp1, temp0, temp_inflow, Tval_t, Tval_b, Uc, n,m, cells, dt) ;
	  imp_bound_rho (rhonp1, p0[2], temp1,  n,m) ;


	  
	  renew_mu_k(mu, kf, temp1, 0, n+1, m+1) ;

	  calc_ts_corr(temp_s1, temp_s0, res_ts1, res_ts0, ks_y, phi_sst, pmedia_imin, pmedia_imax, pmedia_jmin , pmedia_jmax , cells, dt) ;

	  //	  update_ign(temp_s1, t) ;

	  imp_bound_temp_s(temp_s1, Tval_b, pmedia_imin, pmedia_imax, pmedia_jmin,pmedia_jmax) ;

	  

	  ks_compu(ks_x, ks_y, temp_s1, n,m) ;





	  // ATTENTION : UTILDE CALCULATION CONFIGURED FOR WALL ON TOP
	  calc_utilde_corr(utilde_i,utilde_j, u0,v0, resu1,resv1,resu0,resv0, p2D, rhonp1,rho1,mu,phinp1, phin, delta_x, delta_y, y_val, cells,n,m,dt, flag ) ;
	  // boundary conditions for utilde
	  imp_bound_u_tild_conv(Bcomb_ut, Bcomb_ub, p2D, rho1, u_infl, utilde_i, utilde_i0, Uc, n,m, dimy, cells, dt) ;
	  imp_bound_v_tild_conv(Bcomb_vt, Bcomb_vb, v_infl, utilde_j, utilde_j0, Uc, n,m, cells, dt) ;


	  if (Bcomb_p == 'c')
	  		constr_a_c (matA, matA_pos, precond, phinp1, n , m , cells) ;
	  else if (Bcomb_p == 'a')
			constr_a (matA, matA_pos, precond, phinp1, n , m , cells) ;
	  else if (Bcomb_p == 'b')
			constr_a_b (matA, matA_pos, precond, phist, n, m, cells) ;
	  else {
		  cout << endl << "ERROR IN calc_p : these boundary conditions for pressure have not been taken into account" << endl;
	  		return 0;
	  }

	  calc_p (p2D,p, matA,matA_pos,precond,matB, utilde_i, utilde_j, rho1,rho0, rhonp1, phin, phinm1, phinp1, M_re, y_val, cells, n,m,dt, flag, ti) ;
	  // boundary conditions for pressure
	  imp_bound_p(p2D ,rhonp1, n , m , y_val, cells) ;

	  // It's the corrector stage, so we don't want to replace the (n) values with the (*) values
	  // replace(u0,u1, 0, n+1,m+1);
	  // replace(v0,v1, 0, n+1,m+1);

	  calc_vel (u1,v1, utilde_i, utilde_j, p2D, rhonp1, n,m, y_val, cells, dt) ;
	  //impose boundary conditions to u1

	  if (Bcomb_vel == 2)
	    calc_Uc (phiUc, M_re, rhonp1, rho1, phinp1, phin, u_infl, n,m, cells, dt) ;

	  for (int j=1 ; j<=m ; j++)
	    Uc[j] = phiUc /phinp1[n][j] ;

	  imp_bound_u_conv(Bcomb_ut, Bcomb_ub, u1 ,u0, u_infl, uval_t, uval_b, Uc, n , m , cells, dt );
	  imp_bound_v_conv(Bcomb_vt, Bcomb_vb, v1 ,v0, v_infl, vval_t, vval_b, Uc, n , m , cells, dt );



	  //
	  // END OF CORRECTOR STAGE.
	  // calculate new (n+1) values to be used as (n) values at the next predictor step


	  fluxcompu(f1, f2, utilde_i, utilde_j, rhonp1, phinp1, p2D, y_val, cells, n, m, dt) ;


	  if (delta_rank == 'c') {
	    filldelta(delta_x,1,n,m,u1,phinp1, rhonp1,mu);
	    filldelta(delta_y,2,n,m,v1,phinp1, rhonp1,mu);
	  }
	  else if (delta_rank == 's') {
	    filldelta_s(delta_x,u1,v1,rhonp1,mu,phinp1,n,m);
	    filldelta_s(delta_y,u1,v1,rhonp1,mu,phinp1,n,m);
	  }


	  // SOLVER DEBUGGING PACKAGE, corrector
	      /*
	  if (( ti >= ti_min)&&( ti%sample_out == 0 )) {  // output

	    convert_to_2d(mb2d, matB, n,m) ;
	    flushdata( mb2d, n,m, "mb", ti) ;


	    flushdata(res_t1, n,m, "rt",ti) ;
	    constr_b_testa (testB, utilde_i, utilde_j, rho1, rho0, rhonp1, phi, n,m,cells, dt) ;		
	    dimensionalize (mb2d, testB, n,m) ;
	    flushdata (mb2d, n,m, "tb", ti) ;
	    constr_b_testb (testB, utilde_i, utilde_j, rho1, rho0, rhonp1, phi, n,m,cells, dt) ;		
	    flushdata (rho1, n,m, "r1", ti) ;
	    flushdata (rho0, n,m, "r0", ti) ;
	    flushdata (rhonp1, n,m, "rn", ti) ;
	      
	    dimensionalize (mb2d, testB, n,m) ;
	    flushdata (mb2d, n,m, "tb2", ti) ;

	    }
	      */	  


	  // It's the corrector stage, so we don't want to replace the (n) values with the (*) values
	  //replace(resu0, resu1) ;
	  //replace(resv0, resv1) ;

	  rescompu(resu1, resv1, u1,v1, f1,f2, phinp1, rhonp1, mu, F, n,m,cells) ;


	  //
	  // *************
	  // re-interation of the temperatures

	  // first calculate res_t with the new velocities and densities.
	  res_t_compu(res_t1, f1,f2, temp1, temp_s1, rhonp1 ,phinp1, cp, M_re, p0[2], h, kf, cells, n,m) ;

	  res_ts_compu(res_ts1, temp_s1, temp1, phi_snp1, ks_x, M_re, h, cells, pmedia_imin, pmedia_imax, pmedia_jmin , pmedia_jmax ) ;


	  calc_t_corr(temp1,temp0, temp_s0, rhonp1, rhonp1, res_t1,res_t0, p0,M_zero, cp, kf, phinp1, phinp1, phinp1, n,m, pmedia_jmax, cells, dt);

	  //	  update_ign(temp1, t) ;

	  imp_bound_temp_conv(temp1, temp0, temp_inflow,Tval_t, Tval_b, Uc, n,m, cells, dt) ;
	  imp_bound_rho (rhonp1, p0[2], temp1,  n,m) ;


	  
	  renew_mu_k(mu, kf, temp1, 0, n+1, m+1) ;

	  calc_ts_corr(temp_s1, temp_s0, res_ts1, res_ts0, ks_y, phi_snp1, pmedia_imin, pmedia_imax, pmedia_jmin , pmedia_jmax , cells, dt) ;

	  //	  update_ign(temp_s1, t) ;

	  imp_bound_temp_s(temp_s1, Tval_b, pmedia_imin, pmedia_imax, pmedia_jmin,pmedia_jmax) ;

	  

	  ks_compu(ks_x, ks_y, temp_s1, n,m) ;

	  //
	  //************
	  //


	  calc_mass_eq (mass_balance, rhonp1, rho1, phinp1, phin, u1, v1, M_re, n,m, cells, dt) ;


	  h_compu(h,u1, v1, rhonp1, mu, kf, phinp1, n,m );

	  M_re_compu(M_re, temp1, phi_snp1, n,m) ;


	  // Same as the comment above for res_T and res_Ts
	  // replace(res_t0,res_t1,1, n,m);

	  res_t_compu(res_t1, f1,f2, temp1, temp_s1, rhonp1 ,phinp1, cp, M_re, p0[2], h, kf, cells, n,m) ;

	  res_ts_compu(res_ts1, temp_s1, temp1, phi_snp1, ks_x, M_re, h, cells, pmedia_imin, pmedia_imax, pmedia_jmin , pmedia_jmax ) ;

	  // Re-arrange rho matrices
	  replace(rho0, rho1);
	  replace(rho1, rhonp1) ;

	  // Re-arrange phi matrices
	  replace(phinm1, phin) ;
	  replace(phin, phinp1) ;

	  replace(phi_sn, phi_snp1) ;


	  p0[0] = p0[1] ;
	  p0[1] = p0[2] ;


	  


	  // update time

	  t += dt ;

	  //*******
	  //OUTPUT*
	  //*******

	  if (ti % 10 == 0 ) cout << endl << "t=" << t << "\t dt=" << dt << "\t tstep:" << ti << endl;

	  if (( ti >= ti_min)&&( (ti+1)%sample_out == 0 )) // output rho (n-1)
	    flushdata (rho1, n,m, "rho", ti) ;

	  if (( ti >= ti_min)&&( ti%sample_out == 0 )) {  // output

	    cout << endl << endl<< "t=" << t << "\t step :" << ti << endl;

	    //	    flushdata (u1,n,m,"u",ti);
	    //	    flushdata (v1,n,m,"v",ti);
	    //	    flushdata (p2D,n,m,"p",ti);


	    flushbug (u1,0,n+1,m+1,"u",ti);
	    flushbug (v1,0,n+1,m+1,"v",ti);
	    flushbug (p2D,0,n+1,m+1,"p",ti);
 
	    flushbug (temp1, 0,n+1,m+1, "temp", ti) ;
	    flushdata (temp_s1,n,m,"temp_s", ti) ;
	    //flushdata (res_ts1,n,m,"res_s",ti);
	    //flushdata (h, n,m, "h", ti) ;
	    //flushdata (rho1, n,m, "rho", ti) ;
	    flushbug (rho1, 0,n+1,m+1, "rho", ti) ;

	    flushdata (phin, n,m, "phi", ti) ;

	    //flushdata (res_ts1, n,m, "rts", ti) ;

	    //	    flushbug (utilde_i, 0,n+1, m+1, "uti", ti) ;
	    //	    flushbug (utilde_j, 0,n+1, m+1, "utj", ti) ;
						
	    //			flushdata (resu1, n,m, "rucor",ti) ;
	    //			flushdata (resv1, n,m, "rvcor", ti);

	 
	    storetime(t,mass_balance,dt,"t",ti);


	  }


	  // calculate dt for the next cycle

	  
	  umax = 0.;
	  vmax = 0.;

	  for (int i=1 ; i<=n ; i++)
	    for(int j=1 ; j<=m ; j++) {
	      umax = umax < fabs(u1[i][j]) ? fabs(u1[i][j]) : umax ;
	      vmax = vmax < fabs(v1[i][j]) ? fabs(v1[i][j]) : vmax ;
	    }

	  dt_cand_u = CFL * dimx /n /umax ;
	  dt_cand_v = CFL * dimy /m /vmax ;

	  dt_cand_u = dt_cand_u < dt_cand_v ? dt_cand_u : dt_cand_v ;

	  dt = dt_cand_u > 0.005 ? 0.005 : dt_cand_u ;
	  
	  //dt = dt_cand;


	  // Conditional terminate in case of error.
	  if (flag == 1) {


	    ofstream crash;
	    crash.open("crash.out");  if (! crash.is_open()) cout << "input file error, crash.out" ;

	    crash << endl << "Crashing at t=" << t << " and step :" << ti << endl;

	    crash.close();

	    break;
	  }



	}  //  end iterating   t=0 -> tsteps  procedure


	mass_balance++;
	delete [] mass_balance ;

	Uc ++ ;
	delete [] Uc ;

	matB ++;
	delete [] matB ;

	matA++;
	matA_pos++;
	precond++;
	delete [] matA ;
	delete [] matA_pos ;
	delete [] precond ;

	delmat_n(M_re, n,m) ;

	delmat_n(cp, n,m) ;
	delmat_np(kf, n+2, m+2);
	delmat_n(h, n,m) ;

	delmat_np(ks_x, n+2, m+2);
	delmat_np(ks_y, n+2, m+2);

	delmat_n(res_t1, n,m);
	delmat_n(res_t0, n,m);

	delmat_n(res_ts1, n,pmedia_jmax);
	delmat_n(res_ts0, n,pmedia_jmax);

	delmat_np(temp_s0, n+2,m+2) ;
	delmat_np(temp_s1, n+2,m+2) ;

	delmat_np(temp1, n+2,m+2) ;
	delmat_np(temp0, n+2,m+2) ;

	delmat_np(rho0, n+2,m+2) ;
	delmat_np(rho1, n+2,m+2) ;
	delmat_np(rhonp1,n+2,m+2) ;
	delmat_np(rho_st, n+2,m+2) ;

	delmat_np(utilde_i, n+2,m+2) ;
	delmat_np(utilde_i0, n+2,m+2) ;
	delmat_np(utilde_j, n+2,m+2) ;
	delmat_np(utilde_j0, n+2,m+2) ;

	delmat_np(p2D, n+2,m+2) ;

	p++ ;
	delete [] p;

	delmat_n(mb2d,n,m) ;

	u_infl++;
	delete [] u_infl ;
	v_infl++;
	delete [] v_infl ;
       
	delmat_n(resu1, n,m);
	delmat_n(resv1, n,m);
	delmat_n(resu0, n,m);
	delmat_n(resv0, n,m);

	f2++;
	for (int i=0 ; i<n ; i++)
	  delete [] f2[i];
	delete [] f2 ;

	for (int i=0 ; i<n+1 ; i++) {
	  f1[i] ++;
	  delete [] f1[i] ;
	}
	delete [] f1 ;

	delmat_np(u1, n+2,m+2) ;
	delmat_np(v1, n+2,m+2) ;
	delmat_np(u0, n+2,m+2) ;
	delmat_np(v0, n+2,m+2) ;

	delmat_np(mu, n+2,m+2) ;
	
	delmat_n(delta_x, n,m) ;
	delmat_n(delta_y, n,m) ;

	delmat_np(phinm1, n+2,m+2) ;
	delmat_np(phin, n+2,m+2) ;
	delmat_np(phist, n+2,m+2) ;
	delmat_np(phinp1, n+2,m+2) ;

	delmat_np(phi_sn, n+2,m+2) ;
	delmat_np(phi_sst, n+2,m+2) ;
	delmat_np(phi_snp1, n+2,m+2) ;

	delete [] y_val ;

	celldestroy(cells) ;

	return 0;
}




