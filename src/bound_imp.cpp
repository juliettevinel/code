/*
 * bound_imp.cpp
 *
 *  Created on: Jan 11, 2011
 *      Author: paant
 */

#include "./headers/bound_imp.h"



//REMINDER : ALWAYS, FIRST IMPOSE TOP & BOTTOM CONDITIONS, THEN LEFT AND RIGHT - see boundaries.cpp




void imp_bound_p(double** mat, double** rho, int n, int m, double* y_val, grid cells) {

	if (Bcomb_p == 'a') {
		left_bound_neu(mat,n,m,cells);
		right_bound_neu(mat,n,m,cells);
	}

	else if (Bcomb_p =='b') {
	  //		left_bound_dir(mat,Pval_l,n,m,cells);
	  //right_bound_dir(mat,Pval_r,n,m,cells);

	  for (int j=1 ; j<=m ; j++) {
	    mat[0][j] = 2.*Pval_l - mat[1][j] ;

	    mat[n+1][j] = 2.*Pval_r - mat[n][j] ;
	  }

	}

	else if (Bcomb_p == 'c') {
		left_bound_perio(mat,n,m,cells);
		right_bound_perio(mat,n,m,cells);
	}
	else cout << endl << "THERE WAS AN ERROR IN THE imp_bound_p FUNCTION" <<endl ;

	top_bound_neu(mat,rho,n,m,y_val,cells);
	bot_bound_neu(mat,n,m,cells);

}




void imp_bound_u (char type_top, char type_bot, double** mat, double* left_val, double top_val, double bot_val, int n, int m, grid cells) {

	// LEFT -- RIGHT

	if (Bcomb_vel == 1) {
		left_bound_u_perio(mat,n,m);
		right_bound_u_perio(mat,n,m);
	}
	else if (Bcomb_vel == 2) {
	  left_bound_u_dir(mat, left_val, n,m,cells);
		right_bound_u_neu(mat,n,m);
	}
	else if (Bcomb_vel == 3) {
		left_bound_u_neu(mat,n,m);
		right_bound_u_neu(mat,n,m);
	}
	else cout << endl << "THERE WAS AN ERROR IN THE Bcomb_vel VALUE" <<endl ;


	// TOP
	if (type_top == 'w')
		top_bound_u_wall(mat,n,m);

	else if (type_top == 'o')
		top_bound_u_neu(mat,n,m);

	else if (type_top == 'd')
		top_bound_u_dir(mat,top_val,n,m,cells);

	else
		cout << endl << "THERE WAS AN ERROR IN THE type_top VALUE IN imp_bound_u" <<endl ;


	// BOT

	if (type_bot == 'w')
		bot_bound_u_wall(mat,n,m);

	else if (type_bot == 'o')
		bot_bound_u_neu(mat,n,m);

	else if (type_bot == 'd')
		bot_bound_u_dir(mat,bot_val,n,m,cells);

	else
		cout << endl << "THERE WAS AN ERROR IN THE type_bot VALUE IN imp_bound_u" <<endl ;


}


void imp_bound_u_conv (char type_top, char type_bot, double** mat, double** mat0, double* left_val, double top_val, double bot_val, double* Uc, int n, int m, grid cells, double dt) {

	// LEFT -- RIGHT

	if (Bcomb_vel == 1) {
		left_bound_u_perio(mat,n,m);
		right_bound_u_perio(mat,n,m);
	}
	else if (Bcomb_vel == 2) {
	  left_bound_u_dir(mat,left_val,n,m, cells);

	}
	else if (Bcomb_vel == 3) {
		left_bound_u_neu(mat,n,m);
		right_bound_conv(mat,mat0, Uc, n,m, cells, dt);
	}
	else cout << endl << "THERE WAS AN ERROR IN THE Bcomb_vel VALUE" <<endl ;


	// TOP
	if (type_top == 'w')
		top_bound_u_wall(mat,n,m);

	else if (type_top == 'o')
		top_bound_u_neu(mat,n,m);

	else if (type_top == 'd')
		top_bound_u_dir(mat,top_val,n,m,cells);

	else
		cout << endl << "THERE WAS AN ERROR IN THE type_top VALUE IN imp_bound_u" <<endl ;


	// BOT

	if (type_bot == 'w')
		bot_bound_u_wall(mat,n,m);

	else if (type_bot == 'o')
		bot_bound_u_neu(mat,n,m);

	else if (type_bot == 'd')
		bot_bound_u_dir(mat,bot_val,n,m,cells);

	else
		cout << endl << "THERE WAS AN ERROR IN THE type_bot VALUE IN imp_bound_u" <<endl ;


}




void imp_bound_u_tild (char type_top, char type_bot, double** p, double** rho, double* u_infl, double** mat, int n, int m, double y_top, grid cells, double dt) {
  
	cell pt ;


	// LEFT AND RIGHT

	if (Bcomb_vel == 1) {
		//PERIODIC
		for (int j=1 ; j<=m ; j++)
			mat[0][j] = mat[n][j];

		for (int j=1; j <= m ; j++ )
			mat[n+1][j] = mat[1][j];
	}
	else if (Bcomb_vel == 2) {
		//NEUMANN 0 BOTH LEFT AND RIGHT
		for (int j=1 ; j<=m ; j++)
		  mat[0][j] = 2.* u_infl[j] - mat[1][j];

		for (int j=1; j <= m ; j++)
			mat[n+1][j] = mat[n][j];
	}
	else if (Bcomb_vel == 3) {
		//NEUMANN 0 BOTH LEFT AND RIGHT
		for (int j=1 ; j<=m ; j++)
		  mat[0][j] = mat[1][j];

		for (int j=1; j <= m ; j++)
		  mat[n+1][j] = mat[n][j];
	}
	else cout << endl << "THERE WAS AN ERROR IN THE Bcomb_vel VALUE, U_TILDE B-COND IMPLEMENTATION" <<endl ;


	// TOP -- NEUMAN ZERO OR DIRICHLET ZERO
	if (type_top == 'o') {
		for (int i = 0 ; i <= n+1 ; i++)
			mat [i][m+1] = mat[i][m];
	}
	else if (type_top == 'w') {

	  pt.j = m ;

	  //for (int i = 1 ; i <= n ; i++)
	  for (int i = 0 ; i <= n+1 ; i++)
		  {
		    pt.i = i ;

		    //mat [i][m+1] = 2. *dt /valtop(rho,pt,cells) * (derxtop(p, pt,cells) - Ri * y_top * derxtop(rho,pt,cells)) - mat[i][m];
		    mat [i][m+1] = - mat[i][m];
		  }

	}
	else {
		cout << "ERROR APPLYING THE U_TILDE BOUNDARY CONDITION - TOP B" ;
//		return 0;
	}


	// BOTTOM - WALL
	if (type_bot == 'w') {

	  pt.j = 1 ;

	  //		for (int i = 1 ; i <= n ; i++)
	  for (int i = 0 ; i <= n+1 ; i++)
		  {
		    pt.i = i ;
		    mat[i][0] =  - mat[i][1] ;
		  }

		
	}
	else if (type_bot == 'o') { // NEUMANN 0
		for (int i = 0 ; i <= n+1 ; i++)
			mat[i][0] =  mat[i][1] ;
	}
	else {
			cout << "ERROR APPLYING THE U_TILDE BOUNDARY CONDITION - BOTTOM B" ;
	//		return 0;
		}

}

void imp_bound_u_tild_conv (char type_top, char type_bot, double** p, double** rho, double* u_infl, double** mat, double** mat0, double* Uc, int n, int m, double y_top, grid cells, double dt) {
  
	cell pt ;


	// LEFT AND RIGHT

	if (Bcomb_vel == 1) {
		//PERIODIC
		for (int j=1 ; j<=m ; j++)
			mat[0][j] = mat[n][j];

		for (int j=1; j <= m ; j++ )
			mat[n+1][j] = mat[1][j];
	}
	else if (Bcomb_vel == 2) {
		//NEUMANN 0 BOTH LEFT AND RIGHT
		for (int j=1 ; j<=m ; j++)
		  mat[0][j] = 2.* u_infl[j] - mat[1][j];

		right_bound_conv(mat,mat0, Uc, n,m, cells, dt);
		//for (int j=1; j <= m ; j++)
		//	mat[n+1][j] = mat[n][j];
	}
	else if (Bcomb_vel == 3) {
		//NEUMANN 0 BOTH LEFT AND RIGHT
		for (int j=1 ; j<=m ; j++)
		  mat[0][j] = mat[1][j];
		// mat[0][j] = 3.*mat[1][j] - 3.*mat[2][j] + mat[3][j];
		
		right_bound_conv(mat,mat0, Uc, n,m, cells, dt);
		//for (int j=1; j <= m ; j++)
		  //mat[n+1][j] = 3.*mat[n][j] - 3.*mat[n-1][j] + mat[n-2][j];
			
	}
	else cout << endl << "THERE WAS AN ERROR IN THE Bcomb_vel VALUE, U_TILDE B-COND IMPLEMENTATION" <<endl ;


	// TOP -- NEUMAN ZERO OR DIRICHLET ZERO
	if (type_top == 'o') {
		for (int i = 0 ; i <= n+1 ; i++)
			mat [i][m+1] = mat[i][m];
	}
	else if (type_top == 'w') {

	  pt.j = m ;

	  //	for (int i = 1 ; i <= n ; i++)
	  for (int i = 0 ; i <= n+1 ; i++)
		  {
		    pt.i = i ;

		    //		    mat [i][m+1] = 2. *dt /valtop(rho,pt,cells) * (derxtop(p, pt,cells) - Ri * y_top * derxtop(rho,pt,cells)) - mat[i][m];
		    mat [i][m+1] = - mat[i][m];
		  }
	}
	else {
		cout << "ERROR APPLYING THE U_TILDE BOUNDARY CONDITION - TOP B" ;
//		return 0;
	}


	// BOTTOM - WALL
	if (type_bot == 'w') {

	  pt.j = 1 ;

	  //	for (int i = 1 ; i <= n ; i++)
	  for (int i = 0 ; i <= n+1 ; i++)
		  {
		    pt.i = i ;
		    //		    mat[i][0] = 2. *dt /valbot(rho,pt,cells) * derxbot(p, pt,cells) - mat[i][1] ;
		    mat[i][0] = - mat[i][1] ;
		  }
	}
	else if (type_bot == 'o') { // NEUMANN 0
		for (int i = 0 ; i <= n+1 ; i++)
			mat[i][0] =  mat[i][1] ;
	}
	else {
			cout << "ERROR APPLYING THE U_TILDE BOUNDARY CONDITION - BOTTOM B" ;
	//		return 0;
		}

}

void imp_bound_v (char type_top, char type_bot, double** mat, double* left_val, double top_val, double bot_val, int n, int m, grid cells) {

	// LEFT -- RIGHT

	if (Bcomb_vel == 1) {
		left_bound_u_perio(mat,n,m);
		right_bound_u_perio(mat,n,m);
	}
	else if (Bcomb_vel == 2) {
	  left_bound_u_dir(mat,left_val,n,m, cells);
	  right_bound_u_neu(mat,n,m);
	}
	else if (Bcomb_vel == 3) {
		left_bound_u_neu(mat,n,m);
		right_bound_u_neu(mat,n,m);
	}
	else cout << endl << "THERE WAS AN ERROR IN THE Bcomb_vel VALUE" <<endl ;


	// TOP
	if (type_top == 'w')
		top_bound_u_wall(mat,n,m);

	else if (type_top == 'o')
		top_bound_u_neu(mat,n,m);

	else if (type_top == 'd')
		top_bound_u_dir(mat,top_val,n,m,cells);

	else
		cout << endl << "THERE WAS AN ERROR IN THE type_top VALUE IN imp_bound_u" <<endl ;


	// BOT

	if (type_bot == 'w')
		bot_bound_u_wall(mat,n,m);

	else if (type_bot == 'o')
		bot_bound_u_neu(mat,n,m);

	else if (type_bot == 'd')
		bot_bound_u_dir(mat,bot_val,n,m,cells);

	else
		cout << endl << "THERE WAS AN ERROR IN THE type_bot VALUE IN imp_bound_u" <<endl ;


}

void imp_bound_v_conv (char type_top, char type_bot, double** mat, double** mat0, double* left_val, double top_val, double bot_val, double* Uc, int n, int m, grid cells, double dt) {

	// LEFT -- RIGHT

	if (Bcomb_vel == 1) {
		left_bound_u_perio(mat,n,m);
		right_bound_u_perio(mat,n,m);
	}
	else if (Bcomb_vel == 2) {

	  left_bound_u_dir(mat, left_val, n,m, cells);
	  
	  right_bound_conv(mat,mat0, Uc, n,m, cells, dt);

	}
	else if (Bcomb_vel == 3) {
	  
	  left_bound_u_neu(mat,n,m);
	  right_bound_conv(mat,mat0, Uc, n,m, cells, dt);
	}
	else cout << endl << "THERE WAS AN ERROR IN THE Bcomb_vel VALUE" <<endl ;


	// TOP
	if (type_top == 'w')
		top_bound_u_wall(mat,n,m);

	else if (type_top == 'o')
		top_bound_u_neu(mat,n,m);

	else if (type_top == 'd')
		top_bound_u_dir(mat,top_val,n,m,cells);

	else
		cout << endl << "THERE WAS AN ERROR IN THE type_top VALUE IN imp_bound_u" <<endl ;


	// BOT

	if (type_bot == 'w')
		bot_bound_u_wall(mat,n,m);

	else if (type_bot == 'o')
		bot_bound_u_neu(mat,n,m);

	else if (type_bot == 'd')
		bot_bound_u_dir(mat,bot_val,n,m,cells);

	else
		cout << endl << "THERE WAS AN ERROR IN THE type_bot VALUE IN imp_bound_u" <<endl ;


}


void imp_bound_v_tild (char type_top, char type_bot, double** mat, int n, int m) {

	// LEFT AND RIGHT

	if (Bcomb_vel == 1) {
		//PERIODIC
		for (int j=1 ; j<=m ; j++)
			mat[0][j] = mat[n][j];

		for (int j=1; j <= m ; j++ )
			mat[n+1][j] = mat[1][j];
	}
	else if (Bcomb_vel == 2) {
		//NEUMANN 0 BOTH LEFT AND RIGHT
		for (int j=1 ; j<=m ; j++)
		  mat[0][j] = -mat[1][j];

		for (int j=1; j <= m ; j++)
		  mat[n+1][j] = mat[n][j];
	}
	else if (Bcomb_vel == 3) {
		//NEUMANN 0 BOTH LEFT AND RIGHT
		for (int j=1 ; j<=m ; j++)
		  mat[0][j] = mat[1][j];
		
		for (int j=1; j <= m ; j++)
		  mat[n+1][j] = mat[n][j];
	}
	else cout << endl << "THERE WAS AN ERROR IN THE Bcomb_vel VALUE, V_TILDE B-COND IMPLEMENTATION" <<endl ;



	// TOP -- NEUMAN ZERO OR DIRICHLET ZERO
	if (type_top == 'o') {
		for (int i = 0 ; i <= n+1 ; i++)
			mat [i][m+1] = mat[i][m];
	}
	else if (type_top == 'w') {
		for (int i = 0 ; i <= n+1 ; i++)
			mat [i][m+1] = -mat[i][m];
	}
	else {
		cout << "ERROR APPLYING THE U_TILDE BOUNDARY CONDITION - TOP B" ;
//		return 0;
	}


	// BOTTOM - WALL
	if (type_bot == 'w'){
		for (int i = 0 ; i <= n+1 ; i++)
			mat[i][0] = - mat[i][1] ;
	}
	else if (type_bot == 'o') { // NEUMANN 0
		for (int i = 0 ; i <= n+1 ; i++)
			mat[i][0] =  mat[i][1] ;
	}
	else {
			cout << "ERROR APPLYING THE U_TILDE BOUNDARY CONDITION - BOTTOM B" ;
	//		return 0;
		}

}



void imp_bound_v_tild_conv (char type_top, char type_bot, double* u_infl, double** mat, double** mat0, double* Uc, int n, int m,  grid cells, double dt) {
  
	cell pt ;


	// LEFT AND RIGHT

	if (Bcomb_vel == 1) {
		//PERIODIC
		for (int j=1 ; j<=m ; j++)
			mat[0][j] = mat[n][j];

		for (int j=1; j <= m ; j++ )
			mat[n+1][j] = mat[1][j];
	}
	else if (Bcomb_vel == 2) {
		//NEUMANN 0 BOTH LEFT AND RIGHT
		for (int j=1 ; j<=m ; j++)
		  mat[0][j] = -mat[1][j];

		right_bound_conv(mat,mat0, Uc, n,m, cells, dt);
		//for (int j=1; j <= m ; j++)
		//	mat[n+1][j] = mat[n][j];
	}
	else if (Bcomb_vel == 3) {
		//NEUMANN 0 BOTH LEFT AND RIGHT
		for (int j=1 ; j<=m ; j++)
		  mat[0][j] = mat[1][j];
		//t[0][j] = 3.*mat[1][j] - 3.*mat[2][j] + mat[3][j];


		for (int j=1; j <= m ; j++)
		  right_bound_conv(mat,mat0, Uc, n,m, cells, dt);
		//t[n+1][j] = 3.*mat[n][j] - 3.*mat[n-1][j] + mat[n-2][j];
			//mat[n+1][j] = mat[n][j];
	}
	else cout << endl << "THERE WAS AN ERROR IN THE Bcomb_vel VALUE, U_TILDE B-COND IMPLEMENTATION" <<endl ;


	// TOP -- NEUMAN ZERO OR DIRICHLET ZERO
	if (type_top == 'o') {
		for (int i = 0 ; i <= n+1 ; i++)
			mat [i][m+1] = mat[i][m];
	}
	else if (type_top == 'w') {
		for (int i = 0 ; i <= n+1 ; i++)
			mat [i][m+1] = -mat[i][m];
	}
	else {
		cout << "ERROR APPLYING THE U_TILDE BOUNDARY CONDITION - TOP B" ;
//		return 0;
	}


	// BOTTOM - WALL
	if (type_bot == 'w') {
		for (int i = 0 ; i <= n+1 ; i++)
			mat[i][0] = - mat[i][1] ;
	}
	else if (type_bot == 'o') { // NEUMANN 0
		for (int i = 0 ; i <= n+1 ; i++)
			mat[i][0] =  mat[i][1] ;
	}
	else {
			cout << "ERROR APPLYING THE U_TILDE BOUNDARY CONDITION - BOTTOM B" ;
	//		return 0;
		}

}

 


/*
//
//  Extrapolates to obtain the value at the ghost cell
//
void imp_bound_v_tild (int type, char type_top, double** mat, int n, int m) {

	// LEFT AND RIGHT PERIODIC
	for (int j=1 ; j<=m ; j++)
		mat[0][j] = mat[n][j];



	for (int j=1; j <= m ; j++ )
		mat[n+1][j] = mat[1][j];



	// TOP -- NEUMAN ZERO OR DIRICHLET ZERO

	for (int i = 0 ; i <= n+1 ; i++)
	  mat [i][m+1] = 3.*mat[i][m] -3.*mat[i][m-1] + mat[i][m-2] ;


	// BOTTOM - WALL

	for (int i = 0 ; i <= n+1 ; i++)
	  mat[i][0] = 3.* mat[i][1] -3.*mat[i][2] + mat[i][3] ;

}
*/



void imp_bound_Y (double** Y, double val_left, double val_top, double val_bot, int n, int m) {
// Temperature boundary condition. It is called only in the first step.

	if (Bcomb_Y == 1) {
		// LEFT AND RIGHT PERIODIC
		for (int j=1 ; j<=m ; j++)
			Y[0][j] = Y[n][j];

		for (int j=1; j <= m ; j++ )
			Y[n+1][j] = Y[1][j];
	}
	else if (Bcomb_Y == 2) {
		// LEFT DIRICHLET AND RIGHT NEUMANN ZERO
		for (int j=1 ; j<=m ; j++)
            Y[0][j] = 2.- Y[1][j];// Temperature boundary conditions. It is called only in the first step.*val_left - Y[1][j];
		  //Y[0][j] = Y[1][j];
		  

		for (int j=1; j <= m ; j++ )
			Y[n+1][j] = Y[n][j];
	}




	if (Bcomb_Yt == 'o') {
		//top - Neumann zero
		for (int i=0 ; i<=n+1 ; i++)
			Y[i][m+1] = Y[i][m] ;

		//bottom - Neumann zero
		for (int i=0 ; i<=n+1 ; i++)
			Y[i][0] = Y[i][1] ;
	}
	else if (Bcomb_Yt == 'd') {
		//top - Dirichlet
		for (int i=0 ; i<=n+1 ; i++)
			Y[i][m+1] = 2.*val_top - Y[i][m] ;

		//bottom - Dirichlet
		for (int i=0 ; i<=n+1 ; i++)
			Y[i][0] = 2.*val_bot - Y[i][1] ;
	}
	else
		cout << "ERROR APPLYING THE Y BOUNDARY CONDITION - TOP-BOTTOM B" ;

}


void imp_bound_Y_conv (double** Y, double** Y0, double val_left, double val_top, double val_bot, double *Uc, int n, int m, grid cells, double dt) {
// Temperature boundary condition. It is called from the second step onwards

	if (Bcomb_Y == 1) {
		// LEFT AND RIGHT PERIODIC
		for (int j=1 ; j<=m ; j++)
			Y[0][j] = Y[n][j];

		for (int j=1; j <= m ; j++ )
			Y[n+1][j] = Y[1][j];
	}
	else if (Bcomb_Y == 2) {
		// LEFT DIRICHLET
		for (int j=1 ; j<=m ; j++)
		  Y[0][j] = 2.*val_left - Y[1][j];
		  //Y[0][j] = Y[1][j];
		  
		// RIGHT CONVECTIVE
		right_bound_conv(Y,Y0, Uc, n,m, cells, dt);
		//for (int j=1; j <= m ; j++ )
		//			temp[n+1][j] = temp[n][j];
	}




	if (Bcomb_Yt == 'o') {
		//top - Neumann zero
		for (int i=0 ; i<=n+1 ; i++)
			Y[i][m+1] = Y[i][m] ;

		//bottom - Neumann zero
		for (int i=0 ; i<=n+1 ; i++)
			Y[i][0] = Y[i][1] ;
	}
	else if (Bcomb_Yt == 'd') {
		//top - Dirichlet
		for (int i=0 ; i<=n+1 ; i++)
			Y[i][m+1] = 2.*val_top - Y[i][m] ;

		//bottom - Dirichlet
		for (int i=0 ; i<=n+1 ; i++)
			Y[i][0] = 2.*val_bot - Y[i][1] ;
	}
	else
		cout << "ERROR APPLYING THE Y BOUNDARY CONDITION - TOP-BOTTOM B" ;

}


void imp_bound_temp_s (double** temp, double val_bot, int i_min, int i_max, int j_min, int j_max) {


		// LEFT AND RIGHT NEUMANN ZERO
		for (int j=j_min ; j<=j_max ; j++)
			temp[i_min-1][j] = temp[i_min][j];

		for (int j=j_min; j <=j_max ; j++ )
			temp[i_max+1][j] = temp[i_max][j];


	for (int i=i_min ; i<=i_max ; i++)
	  temp[i][j_max+1] = temp[i][j_max] ;
	
	for (int i=i_min ; i<=i_max ; i++)
	  temp[i][j_min-1] = temp[i][j_min] ;

	
	if (Bcomb_Tt == 'o'){
		//bottom - Neumann zero
		for (int i=0 ; i<=n+1 ; i++)
			temp[i][0] = temp[i][1] ;
	}
	else if (Bcomb_Tt == 'd') {
		//bottom - Dirichlet
		for (int i=0 ; i<=n+1 ; i++)
			temp[i][0] = 2.*val_bot - temp[i][1] ;
	}
	else
		cout << "ERROR APPLYING THE TEMP_S BOUNDARY CONDITION - BOTTOM" ;
	

}

void imp_bound_rho (double** rho, double p0, double** temp, int n, int m) {

	if (Bcomb_T == 1) {
		// LEFT AND RIGHT PERIODIC
		for (int j=1 ; j<=m ; j++)
			rho[0][j] = rho[n][j];

		for (int j=1; j <= m ; j++ )
			rho[n+1][j] = rho[1][j];
	}
	else if (Bcomb_T == 2) {
		// LEFT AND RIGHT NEUMANN ZERO
		for (int j=1 ; j<=m ; j++)
		  rho[0][j] = p0 / temp[0][j];
		  //rho[0][j] = rho[1][j];
		  

		for (int j=1; j <= m ; j++ )
		  rho[n+1][j] = p0 / temp[n+1][j];
		  //rho[n+1][j] = rho[n][j];
	}
	else
		cout << "ERROR APPLYING THE RHO BOUNDARY CONDITION - LEFT-RIGHT B" ;

	if (Bcomb_Tt == 'o') {
		//top - NEUMANN ZERO
		for (int i=0 ; i<=n+1 ; i++)
			rho[i][m+1] = rho[i][m] ;

		//bottom - Dirichlet
		for (int i=0 ; i<=n+1 ; i++)
			rho[i][0] = rho[i][1] ;
	}
	else if (Bcomb_Tt == 'd') {
		//top - calculated on the ghost cell, based on p0 and temperature
		for (int i=0 ; i<=n+1 ; i++)
			rho[i][m+1] = p0 /temp[i][m+1] ;
		//	  rho[i][m+1] = 3.*rho[i][m] - 3.*rho[i][m-1] + rho[i][m-2];
	  
		//bottom - Dirichlet
		for (int i=0 ; i<=n+1 ; i++)
			rho[i][0] = p0 /temp[i][0] ;
		//	  rho[i][0] = 3.*rho[i][1] -3.*rho[i][2] + rho[i][3];
	}

	else
		cout << "ERROR APPLYING THE RHO BOUNDARY CONDITION - TOP-BOTTOM B" ;

}


void imp_bound_cells (grid cells, int n, int m) {

  // LEFT BOUND 
  cells.ri[0] = cells.ri[1] ;

  // RIGHT BOUND 
  cells.ri[n+1] = cells.ri[n] ;

  // BOTTOM BOUND 
  cells.rj[0] = cells.rj[1] ;

  // TOP BOUND 
  cells.rj[m+1] = cells.rj[m] ;

}


void imp_bound_mat (double** mat, int n , int m) {


	// PERIODIC
	if (Bcomb_vel == 1) {

		// LEFT BOUND - implementation DOES NOT reach the edges
		for (int j = 1; j <= m ; j++)
			mat[0][j] =mat[n][j] ;


		// RIGHT BOUND - implementation DOES NOT reach the edges
		for (int j = 1 ; j <= m ; j++)
			mat[n+1][j] = mat[1][j] ;

	}
	else {

		// LEFT BOUND - implementation DOES NOT reach the edges
		for (int j = 1; j <= m ; j++)
			mat[0][j] =mat[1][j] ;


		// RIGHT BOUND - implementation DOES NOT reach the edges
		for (int j = 1 ; j <= m ; j++)
		mat[n+1][j] = mat[n][j] ;

	}


	// TOP BOUND - implementation DOES reach the edges
	for (int i = 0 ; i <= n+1 ; i++)
		mat[i][m+1] = mat[i][m] ;


	// BOTTOM BOUND - implementation DOES reach the edges
	for (int i = 0 ; i <= n+1 ; i++)
		mat[i][0] = mat[i][1] ;



}



