/*
 * global.cpp
 *
 *  Created on: Jan 17, 2011
 *      Author: paant
 */

# include "./headers/readfiles.h"



double readdouble(string file) {

	double a;

	ifstream input;
	input.open( file.c_str() );  if (! input.is_open()) cout << "input file error, function readdouble" << file ;
	input >> a;
	input.close();

	return a;

}


int readint (string file) {

	int a;

	ifstream input;
	input.open( file.c_str() );  if (! input.is_open()) cout << "input file error, function readint" << file ;
	input >> a;
	input.close();

	return a;

}


char readchar (string file) {

	char a;

	ifstream input;
	input.open( file.c_str() );  if (! input.is_open()) cout << "input file error, function readchar" << file ;
	input >> a;
	input.close();

	return a;

}
