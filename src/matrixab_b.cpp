/*
 * matrixa.cpp
 *
 *  Created on: Dec 1, 2010
 *      Author: paant
 */


#include "./headers/matrixab_b.h"


//function constructing the matrix A of the linear system
// "row-indexed sparse storage mode" , numerical recipes pg. 78

//suplementary objective - fill in the preconditioner matrix

// a - the A matrix (output)
//pos - the matrix storing the positions of the non-zero values of A
//precond - preconditioner output
//n,m dimentions of the physical grid
// ... the cells matrix



void constr_a_b (double* a, int* pos, double* precond, double** phi, int n, int m, grid** cells ) {

	int alpha;  //  mapping the 2D natural grid to 1D (which is then mapped again for the 1D A matrix)
	double centval;  // the value of the corresponding central cell, which is augmented if close to boundary
	int index = 0;  // overall index to keep right numbering at the a matrix
	int indexpos = 0; // indexing for the pos matrix
	int dummy;      // flag to detect if value recorded is the first non-diagonal item on the row
	int N = n*m + 1 ; //the starting index of non-diagonal elements of A

	for (int j=0; j<m ; j++) {

		for (int i=0; i<n ; i++) {

			alpha = j * n + i ;

			centval = 0.;
			dummy = 0;


			if (j==0) {  // close to the bottom BOUNDARY   //  wall

				centval += sthcell_a(phi,i,j,cells) ;

			}
			else {

				a[N + index] = sthcell_a(phi,i,j,cells);
				pos[alpha] = N + index ;
				pos[N + indexpos] = alpha - n ;
				index++;
				indexpos++;
				dummy = 1;
			}

			if (i==0) { // close to the left BOUNDARY   // dirichlet

				centval += - cells[i][j+1].ri/cells[i+1][j+1].ri * wstcell_a(phi,i,j,cells) ;  // note that the i,j used directly at the cells matr. are augmented by one
																				 // as is automaticaly done within the ***cell_a functions
			}
			else{
				a[N + index] = wstcell_a(phi,i,j,cells) ;
				if (dummy == 0) {  // if it is the FIRST non-diagonal in the row

					pos[alpha] = N + index;
					dummy = 1;
					pos[N + indexpos] = alpha - 1;
					indexpos++;

				}
				else     {// if it is a regular non-diagonal
					pos[N + indexpos] = alpha - 1;
					indexpos++ ;
				}

				index ++;

			}


			if (i==n-1) {// close to the right BOUNDARY

				centval += estcell_a(phi,i,j,cells) ;

			}
			else{
				a[N + index] = estcell_a(phi,i,j,cells) ;
				if (dummy ==0 ) {  // if the FIRST non-diagonal in the row
					pos[alpha] = N + index;
					dummy = 1;
					pos[N + indexpos] =  alpha + 1 ;
					indexpos ++;

				}
				else{
					pos[N + indexpos] =  alpha + 1 ;
					indexpos ++;
				}
				index++;

			}

			if (j==m-1) {// close to the top BOUNDARY   // Neumann
				centval += nthcell_a(phi,i,j,cells) ;
			}
			else{
				a[N + index] = nthcell_a(phi,i,j,cells) ;
				if (dummy ==0 ){ // if the FIRST non-diagonal in the row
					pos[alpha] = N + index;
					dummy = 1;
					pos[N + indexpos] = alpha + n ;
					indexpos ++;
				}
				else {
					pos[N + indexpos] = alpha + n ;
					indexpos ++;
				}
				index++;

			}

			centval += ctcell_a(phi,i,j,cells) ;
			a[alpha] = centval;

			precond[alpha] = centval;   // preconditioner matrix loaded with the values of the diagonal of matrix A

		}  // for (:) rows of the grid

	}  //  for(:) columns of the grid





}  // endof function constr_a


// function that constructs the source matrix B.
// output the resulting matrix
// utild the Utilda matrix who's divergence is calculated
// n,m dimensions and the classic cells




void constr_b_b (double* output, double** utild_i, double** utild_j, double** phi, double** rho, int n , int m ,grid ** cells) {

	int pos;
	cell pt;

	for (int j=0; j < m ; j++) {
		for (int i=0 ; i<n ; i++) {

			pt.i = i+1;
			pt.j = j+1;
			pos = j * n + i ;

			output[pos] = rho[i+1][j+1] * ( derx(utild_i,pt,cells) + dery(utild_j,pt,cells)   ) ;

			if( i ==0) {   // if close to the left boundary, subtract the correcting term  // dirichlet
				output[pos] -=  Pval_l * ( cells[pt.i][pt.j].ri + cells[pt.i-1][pt.j].ri) /  cells[pt.i][pt.j].ri * wstcell_a(phi,i,j,cells) ;
			}
			else if (i == n-1) {   // if close to the right boundary, bubtract the correcting term  // Neumann
				output[pos] -=  Pder_r * ( cells[pt.i][pt.j].ri + cells[pt.i+1][pt.j].ri ) * estcell_a(phi,i,j,cells) ;
			}

//assume for now that Pder_t = 0
//			if ( j == m-1) {  // if close to the top boundary, subtract the correcting term   //   Neumann
//				output[pos] -=  Pder_t * (  cells[pt.i][pt.j].rj + cells[pt.i][pt.j+1].rj) * nthcell_a(phi,i,j,cells) ;
//			}
		}
	}  //  endfor calculate the source term 1D matrix of the linear system

} // endof function constr_b




