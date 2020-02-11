/*
 * output.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: paant
 */


#include "./headers/output.h"



// output to file - only the values of the interior cells of the domain

void flushdata (double** mat, int n, int m, string file, int t) {


	std::stringstream ss;
	ss << t << file << ".txt" ;

	string output2;
	ss >> output2;


	FILE *output= fopen(output2.c_str() , "wb") ;

	for (int  j=1; j<=m ; j++) {
		for (int i=1;i<=n ; i++) {
			fwrite( &mat[i][j], sizeof(double), 1, output) ;
		}
	}


	fclose(output);

}


// output to file - values of both the interior and the ghost cells of the domain

void flushbug (double** mat, int init, int n, int m, string file, int t) {

	std::stringstream ss;
	ss << t << file << "_db.txt" ;

	string output2;
		ss >> output2;

	FILE *output= fopen(output2.c_str() , "wb") ;

	for (int  j=init; j<=m ; j++) {
		for (int i=init;i<=n ; i++) {
			fwrite( &mat[i][j], sizeof(double), 1, output) ;
		}
	}

	fclose(output);


}



// (same as before, for 1D array)

void flushbug (double* mat, int init, int n, string file, int t) {

	std::stringstream ss;
	ss << t << file << "_db.txt" ;

	string output2;
		ss >> output2;

	FILE *output= fopen(output2.c_str() , "wb") ;

	for (int i=init;i<=n ; i++) {
		fwrite( &mat[i], sizeof(double), 1, output) ;
	}


	fclose(output);

}

// (same as before, for 1D integer array)

void flushbug (int* mat, int init, int n, string file, int t) {

	std::stringstream ss;
	ss << t << file << "_db.txt" ;

	string output2;
		ss >> output2;


	FILE *output= fopen(output2.c_str() , "wb") ;

	for (int i=init;i<=n ; i++) {
		fwrite( &mat[i], sizeof(double), 1, output) ;
	}

	fclose(output);

}



// outputs the time file

void storetime (double time, double* mass, double step,  string file, int t) {


	std::stringstream ss;
	ss << t << file << ".txt" ;


	string output2;
	ss >> output2;

	ofstream output;

	output.open( output2.c_str() ) ;

	output.precision(16);

	output << "t = " << time << "\t with dt =" << step << endl ;
	output << "vol_integ_d(rho phi)/dt = " << mass[1] << "\t vol_integ_M = " << mass[2] << endl;
	output << "surf_integ_L = " << mass[3] << "\t surf_integ_R = " << mass[4] << endl;
	output << "surf_integ_T = " << mass[5] << endl;
	output << "mass balance = " << mass[1] + mass[2] + mass[3] + mass[4] + mass[5] << endl;

	output.close();


}



// LOW LEVEL debugging funcs for phi and rij matrices

void flushbug_phi(double** phi, int n, int m) {

	ofstream output;

	output.open("../grid/deb_phi");


	for (int j = -1 ; j <= m+2 ; j++) {
		for (int i = -1 ; i <= n+2 ; i++) {

			output << phi[i][j] << "\t" ;

		}
		output << endl;
	}

	output.close();

}


void flushbug_rij(int n, int m, grid** cells) {

	ofstream outputi, outputj;

	outputi.open("../grid/deb_ri");
	outputj.open("../grid/deb_rj");


	for (int j = -1 ; j <= m+2 ; j++) {
		for (int i = -1 ; i <= n+2 ; i++) {

			outputi << cells[i][j].ri << "\t" ;
			outputj << cells[i][j].rj << "\t" ;

		}
		outputi << endl;
		outputj << endl;
	}

	outputi.close();
	outputj.close();

}

