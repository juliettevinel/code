/*
 * matrixa.cpp
 *
 *  Created on: Dec 1, 2010
 *      Author: paant
 */


#include "./headers/matrixab_c.h"


//function constructing the matrix A of the linear system
// "row-indexed sparse storage mode" , numerical recipes pg. 78

//suplementary objective - fill in the preconditioner matrix

// a - the A matrix (output)
//pos - the matrix storing the positions of the non-zero values of A
//precond - preconditioner output
//n,m dimentions of the physical grid
// ... the cells matrix



void constr_a_b (double* a, int* pos, double* precond, double** phi, int n, int m, grid cells ) {

	double centval = 0.;  // the value of the corresponding central cell, which is augmented if close to boundary

	double nthval, sthval, estval, wstval;

	int alpha = 0;  //  mapping the 2D natural grid to 1D (which is then mapped again for the 1D A matrix)
	int index = 0;  // overall index to keep right numbering at the a matrix

	//flag
	int firstel = 1;      // flag to detect if value recorded is the first non-diagonal item on the row
	int N = n*m + 2 ; //the starting index of non-diagonal elements of A

	cell pt;

	a[N-1] = 0. ;
	pos[N-1] = 5*n*m - 2*(n+m) + 2 ;

	for (int j=1; j<=m ; j++) {

		pt.j = j ;

		for (int i=1; i<=n ; i++) {

			pt.i = i ;

			nthval = nthcell_a(phi,pt,cells);
			sthval = sthcell_a(phi,pt,cells);
			estval = estcell_a(phi,pt,cells);
			wstval = wstcell_a(phi,pt,cells);

			alpha = position_oned(n,m,i,j) ;

			centval = 0.;
			firstel = 1;



			if (j==1) {  // close to the bottom BOUNDARY   // wall

			  centval += sthval ;

			}
			else {

			  a[N + index] = sthval ;
				pos[alpha] = N + index ;
				pos[N + index] = position_oned(n,m,i,j-1) ;
				index++;

				firstel = 0;
			}




			if (i == 1) {  // if it's close to the LEFT boundary
				centval -= wstval;
			}
			else {

			  a[N + index] = wstval ;
				if (firstel == 1) {  // if it is the FIRST non-diagonal in the row

					pos[alpha] = N + index;
					firstel = 0;

				}
					// whether it is the FIRST or NOT
				pos[N + index] = alpha - 1;


				index ++;

			}



			if (i == n) { // if it's close to the RIGHT boundary
				centval -= estval;
			}
			else {

			  a[N + index] = estval;
				if (firstel == 1 ) {  // if the FIRST non-diagonal in the row

					pos[alpha] = N + index;
					firstel = 0;
				}
					// whether it is the FIRST or NOT
				pos[N + index] =  alpha + 1 ;


				index++;

			}




			if (j==m) {// close to the top BOUNDARY   // outflow
			  //				if ( (i<50)||(i>52) )
			  //					centval += nthval;
			  //				else
					centval += nthval;
			}
			else{
			  a[N + index] = nthval ;

				// This coefficient cannot be the first one of the row, either the (i+1) or the (i-1) elements will be written first

				pos[N + index] = position_oned(n,m,i,j+1) ;

				index++;

			}

			centval += -nthval -sthval -estval -wstval ;
			a[alpha] = centval;

			precond[alpha] = centval;   // preconditioner matrix loaded with the values of the diagonal of matrix A

		}  // for (:) rows of the grid

	}  //  for(:) columns of the grid


	 if (a[N-1] != 0.) {
	    cout << "Error in the stensil construction function, matrix a" << endl ;

	    for (int i=1 ; i<= n*m ; i++ )
	      a[i] = 1./0. ;
	  }

	  if (pos[N-1] != 5*n*m - 2*(n+m) + 2) {
	    cout << "Error in the stensil construction function, matrix pos[N-1]" << endl ;

	    for (int i=1 ; i<= n*m ; i++ )
	      a[i] = 1./0. ;
	  }

	  if (pos[1] != n*m + 2) {
		    cout << "Error in the stensil construction function, matrix pos[1]" << endl ;

		    for (int i=1 ; i<= n*m ; i++ )
		      a[i] = 1./0. ;
		  }


}  // endof function constr_a


// function that constructs the source matrix B.
// output the resulting matrix
// utild: the Utilda matrix who's divergence is calculated
// n,m dimensions and the classic cells




void constr_b_b (double* output, double** utilde_i, double** utilde_j, double** rho1, double** rho0, double** rho_new, double** phi1, double** phi0, double** phi_new, double** M_re, int n, int m, double* y_val, grid cells, double dt) {

	int pos;
	cell pt;

	double term_a, term_b;

	for (int j=1; j <= m ; j++) {
		for (int i=1 ; i <= n ; i++) {

			pt.i = i;
			pt.j = j;

			pos = position_oned(n,m, i,j) ;


			// Term from continuity

			output[pos] = (3.*rho_new[i][j]*phi_new[i][j] - 4.*rho1[i][j]*phi1[i][j] + rho0[i][j]*phi0[i][j]) /2. /dt /dt
						  - M_re[i][j] /dt * (rho_s - rho_new[i][j]) /rho_s ;

			//			output[pos] = 0.;

					
			// Divergence rho*phi*u_tilde
			// x-component
			term_a = valfor(phi_new,rho_new,utilde_i, pt,cells) ;
			term_b = valback(phi_new,rho_new,utilde_i, pt,cells) ;

			output[pos] += (term_a - term_b) /2. /cells.ri[pt.i] /dt;

			// y-component
			term_a = valtop(phi_new,rho_new,utilde_j, pt,cells);
			term_b = valbot(phi_new,rho_new,utilde_j, pt,cells);

			output[pos] += (term_a - term_b) /2. /cells.rj[pt.j] /dt ;
			

			// gravity term
			// x - component
			
			term_a = valfor(phi_new, pt,cells) * derxfor(rho_new, pt,cells) ;
			term_b = valback(phi_new, pt,cells) * derxback(rho_new, pt,cells) ;

			output[pos] += Ri *y_val[pt.j] *(term_a - term_b) /2. /cells.ri[pt.i];

			// gravity term
			// y - component

			term_a = valtop(phi_new,y_val, pt,cells) * derytop(rho_new, pt,cells) ;
			term_b = valbot(phi_new,y_val, pt,cells) * derybot(rho_new, pt,cells) ;

			output[pos] += Ri * (term_a - term_b) /2. /cells.rj[pt.j] ;
			


			// Bottom boundary condition for pressure
			//			if (pt.j == 1)
			//			  output[pos] += 0. ! (y at the bottom boundary is zero


			if (i==1)
			  output[pos] -= 2.* Pval_l * wstcell_a(phi_new, pt,cells) ;

			// Top boundary condition
			if ( (j==m) && (Bcomb_vt=='w') ) {

			  output[pos] -= Ri * (y_val[j] + cells.rj[j]) * derytop(rho_new, pt,cells)
			    * nthcell_a(phi_new, pt,cells) * (cells.rj[j] + cells.rj[j+1]) ;

			}

		}
	}  //  endfor calculate the source term 1D matrix of the linear system


} // endof function constr_b




void constr_a_c (double* a, int* pos, double* precond, double** phi, int n, int m, grid cells ) {

	double centval = 0.;  // the value of the corresponding central cell, which is augmented if close to boundary

	double nthval, sthval, estval, wstval;

	int alpha = 0;  //  mapping the 2D natural grid to 1D (which is then mapped again for the 1D A matrix)
	int index = 0;  // overall index to keep right numbering at the a matrix

	//flag
	int firstel = 1;      // flag to detect if value recorded is the first non-diagonal item on the row
	int N = n*m + 2 ; //the starting index of non-diagonal elements of A

	cell pt;

	a[N-1] = 0. ;
	pos[N-1] = 5*n*m - 2*n + 2 ;


	for (int j=1; j<=m ; j++) {

		pt.j = j ;

		for (int i=1; i<=n ; i++) {

			pt.i = i ;

			nthval = nthcell_a(phi,pt,cells);
			sthval = sthcell_a(phi,pt,cells);
			estval = estcell_a(phi,pt,cells);
			wstval = wstcell_a(phi,pt,cells);

			alpha = position_oned(n,m,i,j) ;

			centval = 0.;
			firstel = 1;


			if (j==1) {  // close to the bottom BOUNDARY   // wall

			  centval += sthval ;

			}
			else {

			  a[N + index] = sthval ;
				pos[alpha] = N + index ;
				pos[N + index] = position_oned(n,m,i,j-1) ;
				index++;

				firstel = 0;
			}



			if (i == n) {   // close to the right boundary  ---  periodic

			  a[N + index] = estval;

				if (firstel == 1) {  // if it is the FIRST non-diagonal in the row

					pos[alpha] = N + index ;
					firstel = 0;

				}
					// whether it is the FIRST or NOT
					pos[N + index] = position_oned(n,m,1,j) ;

					index ++ ;

			}



			if (i != 1) {  // if it's NOT close to the LEFT boundary

			  a[N + index] = wstval ;
				if (firstel == 1) {  // if it is the FIRST non-diagonal in the row

					pos[alpha] = N + index;
					firstel = 0;

				}
					// whether it is the FIRST or NOT
				pos[N + index] = alpha - 1;


				index ++;

			}



			if (i != n) { // if it's NOT close to the RIGHT boundary

			  a[N + index] = estval;
				if (firstel == 1 ) {  // if the FIRST non-diagonal in the row

					pos[alpha] = N + index;
					firstel = 0;
				}
					// whether it is the FIRST or NOT
				pos[N + index] =  alpha + 1 ;


				index++;

			}

			if ( i==1 ) {  // if close to the LEFT boundary

			  a[N + index] = wstval ;

				// If close to the left boundary, then p(i+1,j,k) exists, so it cannot be the first element of the row

				pos[N + index] = position_oned(n,m,n,j) ;

				index ++ ;
			}



			if (j==m) {// close to the top BOUNDARY   // outflow
			  centval += nthval;
			}
			else{
			  a[N + index] = nthval ;

				// This coefficient cannot be the first one of the row, either the (i+1) or the (i-1) elements will be written first

				pos[N + index] = position_oned(n,m,i,j+1) ;

				index++;

			}

			centval += -nthval -sthval -estval -wstval ;
			a[alpha] = centval;

			precond[alpha] = centval;   // preconditioner matrix loaded with the values of the diagonal of matrix A

		}  // for (:) rows of the grid

	}  //  for(:) columns of the grid


	 if (a[N-1] != 0.) {
	    cout << "Error in the stensil construction function, matrix a" << endl ;

	    for (int i=1 ; i<= n*m ; i++ )
	      a[i] = 1./0. ;
	  }

	  if (pos[N-1] != 5*n*m - 2*n + 2) {
	    cout << "Error in the stensil construction function, matrix pos[N-1]" << endl ;

	    for (int i=1 ; i<= n*m ; i++ )
	      a[i] = 1./0. ;
	  }

	  if (pos[1] != n*m + 2) {
		    cout << "Error in the stensil construction function, matrix pos[1]" << endl ;

		    for (int i=1 ; i<= n*m ; i++ )
		      a[i] = 1./0. ;
		  }


}  // endof function constr_a


// function that constructs the source matrix B.
// output the resulting matrix
// utild: the Utilda matrix who's divergence is calculated
// n,m dimensions and the classic cells






