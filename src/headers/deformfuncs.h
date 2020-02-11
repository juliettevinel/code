/*
 * deformfuncs.h
 *
 *  Created on: Jan 29, 2014
 *      Author: paant
 */

#ifndef DEFORMFUNCS_H_
#define DEFORMFUNCS_H_

#include "global.h"

#include "trivialfuncs.h"

void deformcompu (vec& tensor, double** u , double** v, double** phi, double** mu, cell pt, grid cells) ;

#endif /* DEFORMFUNCS_H_ */
