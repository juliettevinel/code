/*
 * basicfuncs.cpp
 *
 *  Created on: Nov 30, 2010
 *      Author: paant
 */

#include "./headers/rescalc.h"




//


void rescompu (double** resu, double** resv, double** u, double** v, double** f1, double** f2, double** phi, double** rho, double** mu,
				double F, int n, int m, grid cells)  {


	cell pt; // the center cell

	vec output ;


	for (int i=1 ; i<=n ; i++)
		for (int j=1 ; j<=m ; j++) {

			pt.i = i ;
			pt.j = j ;

			deformcompu (output, u, v, phi, mu, pt, cells) ;

			resu[i][j] = 1. /Re *output.i ;
			//resu[pt.i][pt.j] = 0. ;

			resv[i][j] = 1. /Re *output.j ;
			//resv[pt.i][pt.j] = 0. ;

			


			convcompu(output,u,v,f1, f2, pt,cells) ;

			resu[i][j] -= output.i ;

			resv[i][j] -= output.j ;


			resu[i][j] -= phi[i][j] * F ;

			//resv[pt.i][pt.j] -= rho[pt.i][pt.j] * phi[pt.j] * g ;

		}




} //  end function rescompu












