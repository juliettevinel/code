/*
 * boundaries.cpp
 *
 *  Created on: Dec 4, 2010
 *      Author: panais
 */

#include "./headers/boundaries.h"






/********************************************
 * For thermodynamic variables OTHER THAN U *
 *******************************************/

//TOP

//top - implementing Neumann , the routine does reach the edge
void top_bound_neu (double** mat, double** rho, int n, int m, double* y_val, grid cells) {

  double der, rho_der, y_top ;

  y_top = y_val[m] + cells.rj[m] ;

  if (Bcomb_vt == 'w')
	  for (int i = 0 ; i <= n+1 ; i++) {
	  
		  rho_der = ( rho[i][m+1] - rho[i][m] ) / (cells.rj[m] + cells.rj[m+1]) ;
	  
		  der = Ri * y_top * rho_der  ;
	  
		  mat [i][m+1] = mat[i][m] + der * (cells.rj[m+1] + cells.rj[m]);
	  
	  }
  else
	  for (int i = 0 ; i <= n+1 ; i++)
		  mat [i][m+1] = mat[i][m] ;

}

//BOTTOM

void bot_bound_neu (double** mat, int n , int m, grid cells) {

  //	double rhoval_b, der ; 

  // zero neumann, because y=0 at the bottom boundary
	for (int i = 0 ; i <= n+1 ; i++) {
	  
	  mat[i][0] = mat[i][1] ;

	}

}


// right - implementing Neuman, the routine doesn't reach the edge (corner ghost cells)
void left_bound_neu (double** mat, int n, int m, grid cells) {

	for (int j=1; j <= m ; j++ )
		mat[0][j] = mat[1][j] ;

}

//left - implementing periodic, the routine doesn't reach the edge (corner ghost cells)
void left_bound_perio (double** mat, int n, int m, grid cells) {

	for (int j=1 ; j<=m ; j++)
		mat[0][j] = mat[n][j] ;

}





// right - implementing Neuman, the routine doesn't reach the edge (corner ghost cells)
void right_bound_neu (double** mat, int n, int m, grid cells) {

		for (int j=1; j <= m ; j++ )
			mat[n+1][j] =  mat[n][j] ;

}


// right - implementing periodic, the routine doesn't reach the edge (corner ghost cells)
void right_bound_perio (double** mat, int n, int m, grid cells) {

		for (int j=1; j <= m ; j++ )
			mat[n+1][j] =  mat[1][j];

}



/********************************************
 * ENDFOR thermodynamic variables OTHER THAN U *
 *******************************************/








/*****************************
 * **********FOR U************
 *****************************/


//TOP -- reach the edge
//top u- implementing wall
void top_bound_u_wall (double** mat, int n, int m) {

	int j = m + 1;

	for (int i = 0 ; i <= n+1 ; i++) {
		mat [i][j] = - mat[i][j-1];
	}

}

//top u- implementing outflow
void top_bound_u_neu (double** mat, int n, int m) {

	int j = m + 1;

	for (int i = 0 ; i <= n+1 ; i++) {
		mat [i][j] = mat[i][j-1];
	}

}



void top_bound_u_dir (double** mat, double val, int n, int m, grid cells) {

	int j = m + 1;

	for (int i = 0 ; i <= n+1 ; i++) {
		mat [i][j] = - mat[i][j-1] * cells.rj[j] /cells.rj[j-1] + val * (cells.rj[j] + cells.rj[j-1]) / cells.rj[j-1];
	}

}



//BOTTOM -- reach the edge
//bottom u- implementing wall
void bot_bound_u_wall (double** mat, int n , int m) {

	int j = 0;

	for (int i = 0 ; i <= n+1 ; i++) {
		mat[i][j] = - mat[i][j+1];
	}

}


void bot_bound_u_dir (double** mat, double val, int n, int m, grid cells) {

	int j = 0;

	for (int i = 0 ; i <= n+1 ; i++) {
		mat[i][j] = - mat[i][j+1] * cells.rj[j] / cells.rj[j+1] + val * (cells.rj[j] + cells.rj[j+1]) / cells.rj[j+1];
	}

}


void bot_bound_u_neu (double** mat, int n , int m) {

	int j = 0;

	for (int i = 0 ; i <= n+1 ; i++) {
		mat[i][j] = mat[i][j+1];
	}

}


//LEFT -- don't reach the edge
//left u- implementing periodic
void left_bound_u_perio (double** mat,int n, int m) {  // val the value on the boundary (dirichlet)

	int i = 0;

	for (int j=1 ; j<=m ; j++) {
		mat[i][j] = mat[n][j];
	}

}


//left u- implementing inflow (dirichlet)
void left_bound_u_dir (double** mat, double* val, int n, int m, grid cells) {  // val the value on the boundary (dirichlet)

	int i = 0;

	for (int j=1 ; j<=m ; j++) {
		mat[i][j] = - mat[i+1][j] * cells.ri[i] / cells.ri[i+1] + val[j] * (cells.ri[i] + cells.ri[i+1]) / cells.ri[i+1] ;
	}

}


//LEFT -- don't reach the edge
//left u- implementing Neumann 0
void left_bound_u_neu (double** mat,int n, int m) {

	for (int j=1 ; j<=m ; j++) {
		mat[0][j] = mat[1][j];
	}

}


//RIGHT -- don't reach the edge
//right u- implementing periodic
void right_bound_u_perio (double** mat, int n, int m)  {      // n and m the actual mat dimensions

	int	i = n+1;

	for (int j=1; j <= m ; j++ ) {
		mat[i][j] = mat[1][j];
	}

}

void right_bound_conv (double** mat1, double** mat0, double* Uc, int n, int m, grid cells, double dt) {
  
  double frac1 ;
  double frac2 = 1. /2. /dt ;


  for (int j=1 ; j<=m ; j++) {
    frac1 = Uc[j] /2. /cells.ri[n] ;

    mat1[n+1][j] = mat1[n][j] * (frac1 - frac2) + frac2 * (mat0[n][j] + mat0[n+1][j]) ;
    mat1[n+1][j] /= frac1 + frac2 ;

  }

}

//right u- implementing outflow (Neumann)
void right_bound_u_neu (double** mat, int n, int m)  {      // n and m the actual mat dimensions

	int	i = n+1;

	for (int j=1; j <= m ; j++ ) {
		mat[i][j] = mat[i-1][j];
	}

}





