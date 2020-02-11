/*
 * fluxfuncs.cpp
 *
 *  Created on: Nov 30, 2010
 *      Author: paant
 */

//
// Functions concerning the FLUX term inside the residue.
//
//

#include "./headers/fluxfuncs.h"


//FLUXCOMPU : 	calculate the (nabla) F u , F the aux. flux, at the cell center
//	      	via the F values calculated at the cell walls
//arguments: 	u : matrix of the velocities which are used in the calculations
//		pt : the cell for which the calculations are made
//		cells : pointer to a "cell" type array // aka containing grid / cell info

void convcompu (vec& output,double** u, double** v, double** f1, double** f2, cell pt,  grid cells) {

	double term1a, term1b, term2a, term2b ;

	//x - component ;

	term1a = f1[pt.i][pt.j] * valfor(u, pt,cells) ;
	term1b = f1[pt.i-1][pt.j] * valback(u, pt,cells) ;

	term2a = f2[pt.i][pt.j] * valtop(u, pt, cells) ;
	term2b = f2[pt.i][pt.j-1] * valbot(u, pt, cells) ;

	output.i = (term1a - term1b) /2. /cells.ri[pt.i] + (term2a - term2b) /2. /cells.rj[pt.j] ;


	//y - component ;

	term1a = f1[pt.i][pt.j] * valfor(v, pt,cells) ;
	term1b = f1[pt.i-1][pt.j] * valback(v, pt,cells) ;

	term2a = f2[pt.i][pt.j] * valtop(v, pt,cells) ;
	term2b = f2[pt.i][pt.j-1] * valbot(v, pt,cells) ;

	output.j = (term1a - term1b) /2./cells.ri[pt.i] + (term2a - term2b) /2. /cells.rj[pt.j] ;


}    //endof function fluxcompu




//*************************************************************





void fluxcompu (double** f1, double** f2, double** uti, double**vti, double** rho, double** phi, double** p, double* y_val, grid cells, int n, int m, double dt) {

	cell pt;

	for (int i=0 ; i<=n ; i++)
		for (int j=1 ; j<=m ; j++) {

			pt.i = i;
			pt.j = j;

			// utilde_term
			f1[i][j] = valfor(phi,rho,uti, pt,cells) - valfor(phi, pt,cells) * derxfor(p, pt,cells) *dt
					+ Ri * dt * y_val[j] * valfor(phi, pt,cells) * derxfor(rho, pt,cells) ;

		}

	for (int i=1 ; i<=n ; i++)
		for (int j=0 ; j<=m ; j++) {

			pt.i = i;
			pt.j = j;

			//utilde_term
			f2[i][j] = valtop(phi,rho,vti, pt,cells) - valtop(phi, pt,cells) * derytop(p, pt,cells) *dt
					+ Ri *dt * (phi[i][j+1]*y_val[j+1]*cells.rj[j] + phi[i][j]*y_val[j]*cells.rj[j+1])
					/(cells.rj[j]+cells.rj[j+1]) * derytop(rho, pt,cells) ;

		}


}



// For the first step
void fluxinit (double** f1, double**f2, double** u, double** v, double** rho, double** phi, int n, int m, grid cells) {

	cell pt;

	for (int i=0 ; i<=n ; i++)
		for (int j=1 ; j<=m ; j++) {

			pt.i = i;
			pt.j = j;

			f1[i][j] = valfor(phi,rho,u, pt,cells) ;

		}

	for (int i=1 ; i<=n ; i++)
		for (int j=0 ; j<=m ; j++) {

			pt.i = i;
			pt.j = j;

			f2[i][j] = valtop(phi,rho,v, pt,cells) ;

		}

}



