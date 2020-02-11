/*
 * global.h
 *
 *  Created on: Nov 29, 2010
 *      Author: paant
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_


#include "readfiles.h"



/****************
  COEFFICIENTS
****************/


const double CFL = readdouble("../data/CFL.dat");


const double Re = readdouble("../data/Re.dat"); ;

const double Pr = readdouble("../data/Pr.dat"); ;

const double Ri = readdouble("../data/Ri.dat");;

// For the on-the-fly non-D of k_s
const double k_ref = readdouble("../data/k_ref.dat");
const double Y_ref = readdouble("../data/Y_ref.dat") ;

//dp : particle diameter for a generic isotropic medium
const double dp= readdouble("../data/dp.dat");

//Solid non-D values
const double rho_s = readdouble ("../data/rho_s.dat") ;
const double cp_s = readdouble ("../data/cp_s.dat") ;


const double gama = 1.4 ;


// Chemical kinetics parameters
const double e_activation = readdouble ("../data/e_activation.dat") ;
const double h_formation = readdouble ("../data/h_formation.dat") ;
const double K_reaction = readdouble ("../data/K_reaction.dat") ;



/***************************
 **BOUNDARY CONDITION SETS**
 ***************************/

const int open_d = readint("../data/open_domain.dat") ;

// defining the P boundary condition (top-right) combination
// a: top=wall, right=Dirichlet
// b: top=outflow, right=outflow
const char Bcomb_p = readchar("../data/bcomb_p.dat");

// defining the U boundary condition (left-right) combination
// applies to both u and v
// 1: left,right = periodic
// 2: left=inflow, right=outflow
// 3: left=outflow, right = outflow
const int Bcomb_vel = readint("../data/bcomb_vel_l-r.dat");


// defining the U boundary condition at the top
// w: top = wall,
// o: top = outflow
// d: top = dirichlet
const char Bcomb_ut = readchar("../data/bcomb_ut.dat");
//for V
const char Bcomb_vt = readchar("../data/bcomb_vt.dat");

// defining the U boundary condition at the bottom
// char values, same as in 'top'
const char Bcomb_ub = readchar("../data/bcomb_ub.dat");
//for V
const char Bcomb_vb = readchar("../data/bcomb_vb.dat");


//BOUNDARY CONDITIONS - values
//bottom always wall


// defining the T boundary conditions
// left - right
const int Bcomb_Y = readint("../data/bcomb_Y.dat") ;
// top and bottom
const char Bcomb_Yt = readchar("../data/bcomb_Yt.dat") ;



// values for P matrix
	//left - always inflow

	const double vf_fluid = readdouble("../data/phi.dat");

	const char delta_rank = readchar("../data/delta.dat");

	const int forcing = readint("../data/F.dat");

	const double Pval_l = readdouble("../data/pval_l.dat");


	//right - either Dirichlet, or Neumann
	const double Pval_r = readdouble("../data/pval_r.dat");
	const double Pder_r = readdouble("../data/pder_r.dat");

	const double Tval_t = readdouble("../data/Yval_t.dat");
	const double Tval_b = readdouble("../data/Yval_b.dat");

// ENDOF BOUNDARY CONDITION SETS & VALUES







/**************************
//  STRUCTURE DEFINITIONS
***************************/

// structure grid : grid parameters & scalar quantities


struct grid {

//ri, rj HALF of the rectangular's width or hight
double *ri;
double *rj;

} ;




// structure vector : used to help out with the two components of vector equations

struct vec {

double i , j;

};


// structure cell : containing the coordinates of a cell in the grid

struct cell {
	int i  , j;
};


// structure wally : contain the 4 values, one for each wall direction (N,S,E,W)

struct wally {

double n  ,s  ,e , w ;

};





#endif /* GLOBAL_H_ */
