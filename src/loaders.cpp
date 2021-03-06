/*
 * fillers.cpp
 *
 *  Created on: Dec 5, 2010
 *      Author: panais
 */


#include "./headers/loaders.h"



void cellcreate (grid& cells, int n, int m) {

  cells.ri = new double[n+2];
  cells.rj = new double[m+2];


	double ** x_cor ;
			x_cor = new double*[n+1] ;
			for (int i=0; i<n+1 ; i++)
						x_cor[i] = new double[m+1];

	double ** y_cor ;
			y_cor = new double*[n+1] ;
			for (int i=0; i<n+1 ; i++)
						y_cor[i] = new double[m+1];


	FILE *xin = fopen("../grid/cornx.txt", "rb") ;  if (xin == NULL) cout << "input file error, cornx.txt // loaders.cpp" ;
	FILE *yin = fopen("../grid/corny.txt", "rb") ;  if (yin == NULL) cout << "input file error, corny.txt // loaders.cpp" ;


	for (int j=0 ; j < m+1 ; j++) {
		for (int i = 0 ; i < n+1 ; i++) {

			fread ( &x_cor[i][j], sizeof(double), 1, xin) ;
			fread ( &y_cor[i][j], sizeof(double), 1, yin) ;

		}
	}

	fclose(xin);
	fclose(yin);


	FILE *xcent, *ycent;

	xcent = fopen("../grid/centx.txt", "wb") ; if (xcent == NULL) cout << "output file error, centx.txt  //  loaders.cpp" ;
	ycent = fopen("../grid/centy.txt", "wb") ; if (ycent == NULL) cout << "output file error, centy.txt  //  loaders.cpp" ;


	FILE *ari, *arj;

	ari = fopen ("../grid/ri.txt", "wb") ; if (ari == NULL) cout << "output file error, ri.txt  //  loaders.cpp" ;
	arj = fopen ("../grid/rj.txt", "wb") ; if (arj == NULL) cout << "output file error, rj.txt  //  loaders.cpp" ;


	double xdummy, ydummy, ridummy, rjdummy;

	for (int j=1 ; j <= m ; j++) {
		for (int i=1 ; i <= n ; i++) {

			ridummy = (x_cor[i][j] - x_cor[i-1][j]) / 2. ;
			rjdummy = (y_cor[i][j] - y_cor[i][j-1]) / 2. ;

			cells.ri[i] = ridummy;
			cells.rj[j] = rjdummy;

			fwrite( &ridummy, sizeof(double), 1, ari) ;
			fwrite( &rjdummy, sizeof(double), 1, arj) ;

			xdummy = (x_cor[i][j] + x_cor[i-1][j]) /2. ;
			ydummy = (y_cor[i][j] + y_cor[i][j-1]) /2. ;

			fwrite( &xdummy, sizeof(double), 1, xcent) ;
			fwrite( &ydummy, sizeof(double), 1, ycent) ;

		}
	}

	fclose(ari);
	fclose(arj);
	fclose(xcent);
	fclose(ycent);


	// FREE CORNER MATRICES FROM MEMORY

	for (int i=0; i<n+1 ; i++) {
		delete [] x_cor[i];
		delete [] y_cor[i];
	}

	delete [] x_cor;
	delete [] y_cor;


}  // endof cellinfo func.



void celldestroy (grid& cells) {

  delete [] cells.ri ;
  delete [] cells.rj ;

}




// fills only interior, not ghost => requires application of boundary condition functions

void fillval (double** mat, int n, int m, string file) {


	FILE *values;

	values = fopen(file.c_str(),"rb") ; if (values == NULL)  cout << "input file error at fillval, " << file << " // loaders.cpp" ;

	for (int j=1 ; j<=m ; j++) {
		for (int i=1 ; i<=n ; i++) {
			fread( &mat[i][j], sizeof(double) ,1, values);
		}
	}

//	if (! uvels.eof()) cout << endl << "ERROR WITH READING 'centers_u.txt' - dimensions inconsistency - see loaders.cpp" << endl;
//	if (! vvels.eof()) cout << endl << "ERROR WITH READING 'centers_v.txt' - dimensions inconsistency - see loaders.cpp" << endl;

	fclose(values);

} // endof filling u,v function


// Fills matrices, including values at the ghost cells

void fillval_bug (double** mat, int n, int m, string file) {


	FILE *values;

	values = fopen(file.c_str(),"rb") ; if (values == NULL)  cout << "input file error at fillval, " << file << " // loaders.cpp" ;

	for (int j=0 ; j<=m+1 ; j++) {
		for (int i=0 ; i<=n+1 ; i++) {
			fread( &mat[i][j], sizeof(double) ,1, values);
		}
	}

//	if (! uvels.eof()) cout << endl << "ERROR WITH READING 'centers_u.txt' - dimensions inconsistency - see loaders.cpp" << endl;
//	if (! vvels.eof()) cout << endl << "ERROR WITH READING 'centers_v.txt' - dimensions inconsistency - see loaders.cpp" << endl;

	fclose(values);

} // endof filling u,v function





void fillu_bound (double* mat, int m, string file) {

	FILE *values;

	values = fopen(file.c_str(),"rb") ; if (values == NULL)  cout << "input file error at fillval, " << file << " // loaders.cpp" ;

	for (int j=1 ; j<=m ; j++) {
			fread( &mat[j], sizeof(double) ,1, values);
	}

//	if (! uvels.eof()) cout << endl << "ERROR WITH READING 'centers_u.txt' - dimensions inconsistency - see loaders.cpp" << endl;
//	if (! vvels.eof()) cout << endl << "ERROR WITH READING 'centers_v.txt' - dimensions inconsistency - see loaders.cpp" << endl;

	fclose(values);

} // endof filling u,v function




// fills the 1D p matrix : ONLY THE INTERIOR POINTS
// this filler assumes UNIFORM MESH IN THE X-DIRECTION
