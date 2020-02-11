/* trivialfuncs.cpp
 * Created on: Nov 29, 2010
 *      Author: paant
 *
 *  trivial funcs : computation of 1st and 2nd order derivatives as well as mixed derivatives
 */



//DER** :	x-derivative y-derivative  x^2-derivative y^2-derivative
//arguements : 	matrix a: the variable who's derivative is computed
//		pt the cell at which the der. is evaluated
// 		cells : grid type matrix containing grid data

#include "./headers/trivialfuncs.h"


//first order x derivative - evaluated at the cell center
double derx (double** a, cell pt, grid cells) {

  //return ( a[pt.i+1][pt.j] - a[pt.i-1][pt.j] ) / (cells[pt.i+1][pt.j].ri + 2.*cells[pt.i][pt.j].ri + cells[pt.i-1][pt.j].ri);

  return ( derxfor(a,pt,cells) + derxback(a,pt,cells) ) /2. ;

}


//first order y derivative - evaluated at the cell center
double dery (double** a, cell pt, grid cells) {

  //return  ( a[pt.i][pt.j+1] - a[pt.i][pt.j-1] ) / (cells[pt.i][pt.j+1].rj + 2.*cells[pt.i][pt.j].rj + cells[pt.i][pt.j-1].rj);

  return ( derytop(a,pt,cells) + derybot(a,pt,cells) ) /2. ;

}







//*****************************************************************************
//  used to calculate trivial functions when the cell is next to a boundary


//DER*FOR/BACK/TOP/BOT : compute x, y derivatives on all walls
//			parameters as defined before
//arguements :		a - the matrix of the variable
//			pt the cell at which the der. is evaluated
// 			cells : cell type matrix containing grid data

double derxfor (double** a, cell pt, grid cells) {

	return  (  a[pt.i+1][pt.j] - a[pt.i][pt.j]  ) / (  cells.ri[pt.i+1] + cells.ri[pt.i]  )  ;

}


double derxback (double** a, cell pt, grid cells) {

	return (  a[pt.i][pt.j] - a[pt.i-1][pt.j]  ) / (  cells.ri[pt.i-1] + cells.ri[pt.i]  )  ;

}

double derxtop (double** a, cell pt, grid cells) {

	cell pt_n = pt;
	pt_n.j++;

	return   (derx(a,pt,cells) * cells.rj[pt_n.j] + derx(a,pt_n,cells) * cells.rj[pt.j])
			/ ( cells.rj[pt.j] + cells.rj[pt_n.j] ) ;

}


double derxbot(double** a, cell pt, grid cells) {

	cell pt_n = pt;
	pt_n.j--;

	return   (derx(a,pt,cells) * cells.rj[pt_n.j] + derx(a,pt_n,cells) * cells.rj[pt.j])
			/ ( cells.rj[pt.j] + cells.rj[pt_n.j] ) ;

}


double deryfor (double** a, cell pt, grid cells) {

	cell pt_n = pt;
	pt_n.i++;

	return (dery(a,pt,cells) * cells.ri[pt_n.i] + dery(a,pt_n,cells) * cells.ri[pt.i])
			/ ( cells.ri[pt.i] + cells.ri[pt_n.i] ) ;

}

double deryback (double** a, cell pt, grid cells) {

	cell pt_n = pt;
	pt_n.i--;

	return (dery(a,pt,cells) * cells.ri[pt_n.i] + dery(a,pt_n,cells) * cells.ri[pt.i])
			/ ( cells.ri[pt.i] + cells.ri[pt_n.i] ) ;

}



double derytop  (double** a, cell pt, grid cells) {

	return  (  a[pt.i][pt.j+1] - a[pt.i][pt.j]  ) / (  cells.rj[pt.j+1] + cells.rj[pt.j]  )  ;

}


double derybot (double** a, cell pt, grid cells) {

	return (  a[pt.i][pt.j] - a[pt.i][pt.j-1]  ) / (  cells.rj[pt.j-1] + cells.rj[pt.j]  )  ;

}


//
//   Values at the interfaces
//   FOR/BACK  TOP/BOT
//
//

double valfor(double** a, cell pt, grid cells) {
	return ( a[pt.i][pt.j] * cells.ri[pt.i+1] + a[pt.i+1][pt.j] * cells.ri[pt.i] ) / ( cells.ri[pt.i] + cells.ri[pt.i+1] ) ;
}

double valback(double** a, cell pt, grid cells) {
	return ( a[pt.i][pt.j] * cells.ri[pt.i-1] + a[pt.i-1][pt.j] * cells.ri[pt.i] ) / ( cells.ri[pt.i] + cells.ri[pt.i-1] ) ;
}

double valtop(double** a, cell pt, grid cells) {
	return ( a[pt.i][pt.j] * cells.rj[pt.j+1] + a[pt.i][pt.j+1] * cells.rj[pt.j] ) / ( cells.rj[pt.j] + cells.rj[pt.j+1] ) ;
}

double valtop(double* a, cell pt, grid cells) {
	return ( a[pt.j] * cells.rj[pt.j+1] + a[pt.j+1] * cells.rj[pt.j] ) / ( cells.rj[pt.j] + cells.rj[pt.j+1] ) ;
}

double valbot(double** a, cell pt, grid cells) {
	return ( a[pt.i][pt.j] * cells.rj[pt.j-1] + a[pt.i][pt.j-1] * cells.rj[pt.j] ) / ( cells.rj[pt.j] + cells.rj[pt.j-1] ) ;
}

double valbot(double* a, cell pt, grid cells) {
	return ( a[pt.j] * cells.rj[pt.j-1] + a[pt.j-1] * cells.rj[pt.j] ) / ( cells.rj[pt.j] + cells.rj[pt.j-1] ) ;
}



//  Interpolate at the interface a product of 2 functions

double valfor (double** a, double** b, cell pt, grid cells) {

	return ( a[pt.i+1][pt.j] * b[pt.i+1][pt.j] * cells.ri[pt.i] + a[pt.i][pt.j] * b[pt.i][pt.j] * cells.ri[pt.i+1] )
			/ ( cells.ri[pt.i] + cells.ri[pt.i+1] ) ;

}

double valback (double** a, double** b, cell pt, grid cells) {

	return ( a[pt.i-1][pt.j] * b[pt.i-1][pt.j] * cells.ri[pt.i] + a[pt.i][pt.j] * b[pt.i][pt.j] * cells.ri[pt.i-1] )
			/ ( cells.ri[pt.i] + cells.ri[pt.i-1] ) ;

}

double valtop (double** a, double** b, cell pt, grid cells) {

	return ( a[pt.i][pt.j+1] * b[pt.i][pt.j+1] * cells.rj[pt.j] + a[pt.i][pt.j] * b[pt.i][pt.j] * cells.rj[pt.j+1] )
			/ ( cells.rj[pt.j] + cells.rj[pt.j+1] ) ;

}

double valbot (double** a, double** b, cell pt, grid cells) {

	return ( a[pt.i][pt.j-1] * b[pt.i][pt.j-1] * cells.rj[pt.j] + a[pt.i][pt.j] * b[pt.i][pt.j] * cells.rj[pt.j-1] )
			/ ( cells.rj[pt.j] + cells.rj[pt.j-1] ) ;

}

double valtop(double** a, double* b, cell pt, grid cells) {

	return ( a[pt.i][pt.j+1] * b[pt.j+1] * cells.rj[pt.j] + a[pt.i][pt.j] * b[pt.j] * cells.rj[pt.j+1] )
			/ ( cells.rj[pt.j] + cells.rj[pt.j+1] ) ;

}


double valbot(double** a, double* b, cell pt, grid cells) {

	return ( a[pt.i][pt.j-1] * b[pt.j-1] * cells.rj[pt.j] + a[pt.i][pt.j] * b[pt.j] * cells.rj[pt.j-1] )
			/ ( cells.rj[pt.j] + cells.rj[pt.j-1] ) ;

}


//  Interpolate at the interface a product of 3 functions

double valfor (double** a, double** b, double** c, cell pt, grid cells) {

	return ( a[pt.i+1][pt.j]*b[pt.i+1][pt.j]*c[pt.i+1][pt.j] * cells.ri[pt.i] + a[pt.i][pt.j]*b[pt.i][pt.j]*c[pt.i][pt.j] * cells.ri[pt.i+1] )
			/ ( cells.ri[pt.i] + cells.ri[pt.i+1] ) ;

}

double valback (double** a, double** b, double** c, cell pt, grid cells) {

	return ( a[pt.i-1][pt.j]*b[pt.i-1][pt.j]*c[pt.i-1][pt.j] * cells.ri[pt.i] + a[pt.i][pt.j]*b[pt.i][pt.j]*c[pt.i][pt.j] * cells.ri[pt.i-1] )
			/ ( cells.ri[pt.i] + cells.ri[pt.i-1] ) ;

}

double valtop (double** a, double** b, double** c, cell pt, grid cells) {

	return ( a[pt.i][pt.j+1]*b[pt.i][pt.j+1]*c[pt.i][pt.j+1] * cells.rj[pt.j] + a[pt.i][pt.j]*b[pt.i][pt.j]*c[pt.i][pt.j] * cells.rj[pt.j+1] )
			/ ( cells.rj[pt.j] + cells.rj[pt.j+1] ) ;

}

double valbot (double** a, double** b, double** c, cell pt, grid cells) {

	return ( a[pt.i][pt.j-1]*b[pt.i][pt.j-1]*c[pt.i][pt.j-1] * cells.rj[pt.j] + a[pt.i][pt.j]*b[pt.i][pt.j]*c[pt.i][pt.j] * cells.rj[pt.j-1] )
			/ ( cells.rj[pt.j] + cells.rj[pt.j-1] ) ;

}

//**************************************************************
// overloaded replace function performs   a = b  ///  i.e. prepares the matrices for 1 step advancement
//*************************************************************

void replace (double* a, double* b, int dim) {

	for (int i=0 ; i<= dim ; i++) {
		a[i] = b[i];
	}

}


void replace (double** a, double** b, int init, int n, int m) {

	for (int j = init; j<=m ;j++) {
		for (int i = init ; i<=n; i++) {
			a[i][j] = b[i][j] ;
		}
	}

}



void replace (double** &a, double** &b) {

  double** temp ;

  temp = a;

  a = b;

  b = temp ;

}

void replace (double* &a, double* &b) {

  double* temp ;

  temp = a;

  a = b;

  b = temp ;

}

// transform the 1D result of the p matrix to the 2D equivalent, corresponding to the grid
// the new 2D matrix INCLUDES POSITIONS FOR GHOST CELLS
// n,m are the dimensions merged at 1D in mat. mat2D has 2 more elements for each dimention,
// including that way the ghost elements

void convert_to_2d (double** mat2d, double* mat, int n, int m) {

	for (int i=1 ; i<=n ; i++)
		for (int j=1 ; j<=m ; j++)
			mat2d[i][j] = mat[(j-1) * n + i] ;

}

void convert_to_1d ( double* mat, double** mat2d, int n, int m) {

	for (int i=1 ; i<=n ; i++)
		for (int j=1 ; j<=m ; j++)
			mat[(j-1) * n + i] = mat2d[i][j] ;

}


void initialize (double** mat, int init, int n, int m) {

	for (int i=init ; i<=n ; i++)
		for (int j=init ; j<=m ; j++)
			mat[i][j] = 0.;

}




