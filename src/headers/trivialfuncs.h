/*
 * trivialfuncs.h
 *
 *  Created on: Nov 29, 2010
 *      Author: paant
 */

#ifndef TRIVIALFUNCS_H_
#define TRIVIALFUNCS_H_


# include <math.h>
# include "global.h"


using namespace std;


//***********************
// FUNCTION DECLARATIONS
//***********************



double derx (double** a, cell pt, grid cells);
double dery (double** a, cell pt, grid cells);

double derxfor (double** a, cell pt, grid cells);
double derxback (double** a, cell pt, grid cells);
double derxtop (double** a, cell pt, grid cells);
double derxbot(double** a, cell pt, grid cells) ;

double deryfor (double** a, cell pt, grid cells);
double deryback (double** a, cell pt, grid cells) ;
double derytop  (double** a, cell pt, grid cells);
double derybot (double** a, cell pt, grid cells);

double valfor(double** a, cell pt, grid cells);
double valback(double** a, cell pt, grid cells);

double valtop(double** a, cell pt, grid cells) ;
double valtop(double* a, cell pt, grid cells) ;

double valbot(double** a, cell pt, grid cells);
double valbot(double* a, cell pt, grid cells) ;

double valfor(double** a, double** b, cell pt, grid cells);
double valback(double** a, double** b, cell pt, grid cells);
double valtop(double** a, double** b, cell pt, grid cells) ;
double valbot(double** a, double** b, cell pt, grid cells);

double valtop(double** a, double* b, cell pt, grid cells) ;
double valbot(double** a, double* b, cell pt, grid cells);

double valfor(double** a, double** b, double** c, cell pt, grid cells);
double valback(double** a, double** b, double** c, cell pt, grid cells);
double valtop(double** a, double** b, double** c, cell pt, grid cells) ;
double valbot(double** a, double** b, double** c, cell pt, grid cells);

void replace (double* a, double*b, int dim) ;
void replace (double** a, double** b, int init, int n, int m);

void replace (double** &a, double** &b) ;
void replace (double* &a, double* &b) ;

void convert_to_2d (double** mat2d, double* mat, int n, int m);
void convert_to_1d (double* mat, double** mat2d, int n, int m) ;

void initialize (double** mat, int init, int n, int m) ;

#endif /* TRIVIALFUNCS_H_ */
