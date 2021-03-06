/*
 * deformfuncs.cpp
 *
 *  Created on: Jan 29, 2014
 *      Author: paant
 */


#include "headers/deformfuncs.h"


// computing the deformation tensor term

void deformcompu (vec& tensor, double** u , double** v, double** phi, double** mu, cell pt, grid cells) {

	double terma_e, terma_w;  // term "a", east and west wall values in order to compute x-derivative on cell-center
	double termb_n, termb_s;  // term "b", north and south wall values in order to compute y-derivative on cell-center

	double phixmu_e = valfor(phi,mu, pt,cells) ;
	double phixmu_w = valback(phi,mu, pt,cells) ;
	double phixmu_n = valtop(phi,mu, pt,cells) ;
	double phixmu_s = valbot(phi,mu, pt,cells) ;

	terma_e = 4./3. * derxfor(u,pt,cells) - 2./3. *deryfor(v,pt,cells)  ;
	terma_w = 4./3. * derxback(u,pt,cells) - 2./3. *deryback(v,pt,cells) ;

	termb_n = derxtop(v,pt,cells) ;
	termb_s = derxbot(v,pt,cells) ;

	tensor.i =	(terma_e * phixmu_e - terma_w * phixmu_w) /2. /cells.ri[pt.i] +
				(termb_n * phixmu_n - termb_s * phixmu_s) /2. /cells.rj[pt.j]
			; // the outmost parenthesis includes the divergence

	//the same follows of the 2nd component of the  div(T) matix

	terma_e = deryfor(u,pt,cells) + derxfor(v,pt,cells) ;
	terma_w = deryback(u,pt,cells) + derxback(v,pt,cells) ;

	termb_n = derxtop(u,pt,cells)  ;
	termb_s = derxbot(u,pt,cells)  ;


	tensor.j =	(terma_e * phixmu_e - terma_w * phixmu_w) /2. /cells.ri[pt.i]
			- 2./3. * (termb_n * phixmu_n - termb_s * phixmu_s) /2./cells.rj[pt.j]
			; // the outmost parenthesis includes the divergence

} // endof function deformcompu


