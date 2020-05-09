/*
 * init.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: paant
 */


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#include "./headers/readfiles.h"


using namespace std;

double tanfunc (double a, int res, int j , double height);
double tanfunc_sym (double a, int res, int j , double height);
double profile (double* y_cor, int j);


// !!!!  WARNING - THE PHI INTEGRATION IS NOT CORRECT - THE STRUCTURE CORRESPONDS TO CORNERS AND NOT CELL CENTERS





int main () {

	double Re = readdouble("../data/Re.dat") ;

	double F = -8. /9. /Re ;

	ifstream dimensions;
//	ofstream xpos, ypos;  // at corners of the cells
	FILE *xpos, *ypos;

//	ofstream u, v, p, phi, rho, mu ;	// at the cell centers
	FILE *u, *v, *p, *phi, *temp, *temps, *Y;   // at the cell centers
//	ofstream u_lb, v_lb; // u and v left boundary profiles
	FILE *u_lb, *v_lb; // u and v left boundary profiles
//	ofstream u_up_conv, v_up_conv, p_up_conv; // u, v, p boundary matrices helping with the top convective boundary condition
	FILE *u_up_conv, *v_up_conv, *p_up_conv; // u, v, p boundary matrices helping with the top convective boundary condition
//	ofstream u_rt_conv, v_rt_conv, p_rt_conv; // u, v, p boundary matrices helping with the right convective boundary condition
	FILE *u_rt_conv, *v_rt_conv, *p_rt_conv; // u, v, p boundary matrices helping with the right convective boundary condition

	dimensions.open ("../data/grid.dat") ;  if (! dimensions.is_open()) cout << "input file error, grid.dat" ;

	double dimx, dimy ; // x, y dimensions
	int resx, resy ;  // x, y resolution

	int pmedia_jmin, pmedia_jmax, pmedia_imin, pmedia_imax ; // the y-resolution dedicated to the porous medium


	dimensions >> resx >> resy >> dimx >> dimy  >> pmedia_jmin >> pmedia_jmax >> pmedia_imin >> pmedia_imax ;
	dimensions.close();

	if (pmedia_jmax > resy) {
		cout << "resolusion dedicated for the porus medium cannot be higher than that of the whole domain" << endl;
		return 0;
	}

	cout  << "dimensions \t" << dimx << "\t" << dimy << "\nresolusions\t" << resx << "\t" << resy ;
	cout << endl << "pmedia resolution \t" << pmedia_jmax ;

	xpos = fopen("../grid/cornx.txt","wb") ; if (xpos == NULL) cout << "output file error, corners_x.txt // init.cpp" ;
	ypos = fopen("../grid/corny.txt","wb") ; if (ypos == NULL) cout << "output file error, corners_y.txt // init.cpp" ;





	const double Tval_t = readdouble("../data/tval_t.dat");
	const double Tval_b = readdouble("../data/tval_b.dat");

	const double Ri = readdouble("../data/Ri.dat") ;



	double dx = dimx / resx;
	double dy;

	double * y_corns;
	y_corns = new double [resy + 1]; //  assistive array to hold y values of the corners of the cells.
	int dummy = 0;

	//	if (pmedia_jmax == 0) dy = 0. ;
//	else dy = pmedia_lim / pmedia_jmax;
//	else dy = dimy / resy;
	dy = dimy / resy;



	//****************************************************
	// DOMAIN CONSTRUCTION - AERIA COVERED BY POROUS MEDIUM
	//****************************************************

	double xdummy, ydummy;

	for (int j = 0 ; j < pmedia_jmax + 1 ; j++) {
		for (int i = 0 ; i < resx + 1 ; i++) {

			xdummy = i*dx ;   ydummy = j*dy ;
			fwrite(&xdummy, sizeof(double), 1, xpos);
			fwrite(&ydummy, sizeof(double), 1, ypos);

		}

		y_corns[dummy] = j*dy;
		dummy++;
	}


	//**********************************************
	// DOMAIN CONSTRUCTION - AREA OF PURE-FLUID MEDIUM
	//**********************************************

	int air_resy = resy - pmedia_jmax ;
	//	double air_height = dimy - pmedia_lim ;

	//	cout << endl <<"free flow y- domain " << air_height << endl;
	cout << "air resolution " << air_resy << endl;


	for (int j = 1 ; j <= air_resy ; j++) {
		for (int i = 0 ; i < resx + 1 ; i++) {


			xdummy = i*dx;
			//			ypos << pmedia_lim + tanfunc_sym( .98346, air_resy, j, air_height ) << "\t" ;
			//			ypos << pmedia_lim + tanfunc( .978, air_resy, j, air_height ) << "\t" ;
			ydummy = (j + pmedia_jmax)*dy;

			fwrite(&xdummy, sizeof(double), 1, xpos);
			fwrite(&ydummy, sizeof(double), 1, ypos);

		}
//		y_corns[dummy] = pmedia_lim + tanfunc( .978, air_resy, j, air_height )  ;   // assistive matrix for the blausius profile implementation
		y_corns[dummy] = (j + pmedia_jmax)*dy;
		dummy++;

	}

	fclose(xpos);
	fclose(ypos);




	/**********************
	***********************
	// variables initiation
	***********************
	***********************/

	u = fopen("../init/centu.txt","wb") ;  if (u == NULL) cout << "output file error, centu.txt // init.cpp" ;
	v = fopen("../init/centv.txt","wb") ;  if (v == NULL) cout << "output file error, centv.txt // init.cpp" ;
	p = fopen("../init/centp.txt","wb") ;  if (p == NULL) cout << "output file error, centp.txt // init.cpp" ;

	phi = fopen("../init/centphi.txt","wb") ;  if (phi == NULL) cout << "output file error, centphi.txt // init.cpp" ;
	temp = fopen("../init/centtemp.txt","wb"); if (temp == NULL) cout << "output file error, centtemp.txt // init.cpp" ;
	temps = fopen("../init/centtemp_s.txt","wb"); if (temp == NULL) cout << "output file error, centtemps.txt // init.cpp" ;
    Y = fopen("../init/centY.txt","wb"); if (temp == NULL) cout << "output file error, centY.txt // init.cpp" ;
    
	u_lb = fopen("../init/u_dir.txt","wb") ;  if (u_lb == NULL) cout << "output file error, u_dir.txt // init.cpp" ;
	v_lb = fopen("../init/v_dir.txt","wb") ;  if (v_lb == NULL) cout << "output file error, v_dir.txt // init.cpp" ;

	u_up_conv = fopen("../init/u_up_conv.txt","wb") ;  if (u_up_conv == NULL) cout << "output file error, u_up_conv.txt // init.cpp" ;
	v_up_conv = fopen("../init/v_up_conv.txt","wb") ;  if (v_up_conv == NULL) cout << "output file error, v_up_conv.txt // init.cpp" ;
	p_up_conv = fopen("../init/p_up_conv.txt","wb") ;  if (p_up_conv == NULL) cout << "output file error, p_up_conv.txt // init.cpp" ;

	u_rt_conv = fopen("../init/u_rt_conv.txt","wb") ;  if (u_rt_conv == NULL) cout << "output file error, u_rt_conv.txt // init.cpp" ;
	v_rt_conv = fopen("../init/v_rt_conv.txt","wb") ;  if (v_rt_conv == NULL) cout << "output file error, v_rt_conv.txt // init.cpp" ;
	p_rt_conv = fopen("../init/p_rt_conv.txt","wb") ;  if (p_rt_conv == NULL) cout << "output file error, p_rt_conv.txt // init.cpp" ;

	
	/********************
	input data initiation	
	********************/
	
	double Pval_l, Pval_r,  slope;	// data needed for p gradient construction (initial condition)

	ifstream bcond, value;

	bcond.open("../data/pval_l.dat");  if (! bcond.is_open()) cout << "input file error, pval_l.dat" ;
	bcond >> Pval_l ;
	bcond.close();

	bcond.open("../data/pval_r.dat"); if (! bcond.is_open()) cout << "input file error, pval_r.dat" ;
	bcond >> Pval_r;
	bcond.close();



	double phi_val;
	value.open("../data/phi.dat"); if (! value.is_open()) cout << "input file error, phi.dat" ;
	value >> phi_val;
	value.close();


	slope = (Pval_r - Pval_l)/dimx ;
//	double a,b; // auxiliary vars, to help with the profile of u and v in case of dirichlet condition on the left

	


	double phidummy, udummy, vdummy, pdummy, tempdummy, dispdummy;
    
    
    ifstream fichier("../prehand_new/src/u385_0.48.txt");
           if (fichier)
           {
               string ligne;
               while(getline(fichier,ligne))
               {
                    for (int i=1 ; i<= resx; i++) {
                   //cout << ligne << endl;
                   udummy= atof(ligne.c_str());
                   fwrite(&udummy, sizeof(double), 1, u);
                    }
               }
               
           }
          
          else
          { cout << "Erreur: impossible d'ouvrir le fichier udummy.txt" << endl;
          }

        fichier.close();
    
    
    


	// FOR TEMPERATURE PROFILE
	double k_betta = 0.00037* 273.15 /0.0243 ;
	double k_alpha = 0.0103345 /0.0243;
	double ks = 1.8 *(k_alpha + k_betta * 1.1);
	double kf = 1.068993045;
	
	
	double y_val, x_val;

	for (int j =1 ; j <= resy ; j++) {

//		a = profile (y_corns,j) ;
//		b = profile (y_corns,j) ;
		for (int i=1 ; i <= resx ; i++) {

		  //		  y_val = (j+0.5)*dy - pmedia_lim ;
		  y_val = (j-0.5)*dy ;
		  x_val = (i-0.5)*dx ;

		  phidummy = 1. ;


		  
		  /*
		  if ((i>=pmedia_imin)&&(i<=pmedia_imax)&&(j<=pmedia_jmax)&&(j>=pmedia_jmin))
		    phidummy = (1.-phi_val)/4.*(-1) + 1. ;

		  if ((i>=pmedia_imin+1)&&(i<=pmedia_imax-1)&&(j<=pmedia_jmax-1)&&(j>=pmedia_jmin+1))
		    phidummy = (1.-phi_val)/4.*(-2) + 1. ;

		  if ((i>=pmedia_imin+2)&&(i<=pmedia_imax-2)&&(j<=pmedia_jmax-2)&&(j>=pmedia_jmin+2))
		    phidummy = (1.-phi_val)/4.*(- 3) + 1. ;

		  if ((i>=pmedia_imin+3)&&(i<=pmedia_imax-3)&&(j<=pmedia_jmax-3)&&(j>=pmedia_jmin+3))
		    phidummy = phi_val ;
		  */


		  if ((i>=pmedia_imin)&&(i<=pmedia_imax)&&(j<=pmedia_jmax)&&(j>=pmedia_jmin))
		    phidummy = (1.-phi_val)/2.*(-1) + 1. ;

		  if ((i>=pmedia_imin)&&(i<=pmedia_imax)&&(j<=pmedia_jmax-1)&&(j>=pmedia_jmin))
		    phidummy = phi_val ;
		    

		  

		  /*

		  if ((i==202)||(j==pmedia_jmax -3)||(i==302))
		    phidummy = (1.-phi_val)/4.*(-3) + 1. ;

		    if ((i==201)||(j==pmedia_jmax -2)||(i==301))
		    phidummy = (1.-phi_val)/4.*(-2) + 1. ;
		  
		    if ((i==200)||(j==pmedia_jmax -1)||(i==300))
		    phidummy = (1.-phi_val)/4.*(- 1) + 1. ;

		  */


		  // phidummy = (1. - phi_val)/4. * j - (1. - phi_val)/4.*pmedia_jmax + 1. ;

			fwrite(&phidummy, sizeof(double),1 , phi) ;

			//udummy = -dimy * dimy /2. * 1./2. * ( y_val*y_val /dimy /dimy - y_val /dimy   )  ;
			//udummy = 0.;

			//tempdummy = atan(50. * y_val ) * (Tval_t - Tval_b)/2. / M_PI * 2. + (Tval_t + Tval_b)/2.;
			
			tempdummy = 1. ;
            dispdummy= 0.;
			//tempdummy = -tanh( (x_val - 1.)*4. )/2. + 1.5 ;
			//if (x_val < 1.8)
			//tempdummy = -1. /5. * x_val + 2. ;
			//else
			//  tempdummy = 1. ;
			
			/*
			if (i>197 && i<208 && j>17 && j<28) 
			  tempdummy = 2. ;

			if (i>198 && i<207 && j>18 && j<27) 
			  tempdummy = 3. ;

			if (i>199 && i<206 && j>19 && j<26) 
			  tempdummy = 4. ;

			if (i>200 && i<205 && j>20 && j<25) 
			  tempdummy = 5. ;
			 */

			fwrite(&tempdummy, sizeof(double), 1, temp) ;
			fwrite(&tempdummy, sizeof(double), 1, temps) ;
            fwrite(&dispdummy, sizeof(double), 1, Y) ;

// U profile all over the domain (uniform)
			//udummy is determined higher, now it depends on y
			//udummy = 0.;
			//fwrite(&udummy, sizeof(double), 1, u);
//			u << 1. << "\t";
//			u << a << "\t";
			vdummy = 0.;
			fwrite(&vdummy, sizeof(double), 1, v);
//			v << b << "\t";




			//pdummy = slope * (i+.5)*dx + Pval_l ;
			//pdummy = -Ri * (j+.5)*dy ;
			pdummy = 0. ;
			fwrite(&pdummy, sizeof(double), 1, p);


		}



// U profile at left boundary (dirichlet condition on the left)
//
		fwrite(&udummy, sizeof(double), 1, u_lb);

		fwrite(&vdummy, sizeof(double), 1, v_lb);





//		constructing auxiliary matrices for convective boundary condition on the right;
//		OLD FORMATING, TEXT
//		u_rt_conv << 0. << "\t" << 0. << endl;
//		v_rt_conv << 0. << "\t" << 0. << endl;
//		p_rt_conv << 0. << "\t" << 0. << endl;
//		NEW FORMATING, BIN
		udummy = 0.;
		fwrite(&udummy, sizeof(double), 1, u_rt_conv);
		fwrite(&udummy, sizeof(double), 1, u_rt_conv);
		vdummy = 0.;
		fwrite(&vdummy, sizeof(double), 1, v_rt_conv);
		fwrite(&vdummy, sizeof(double), 1, v_rt_conv);
		pdummy = 0.;
		fwrite(&pdummy, sizeof(double), 1, p_rt_conv);
		fwrite(&pdummy, sizeof(double), 1, p_rt_conv);
	}


// top convective boundaries and first internal cells

	for (int j = 0 ; j<=1 ; j++) {

		for (int i = -1 ; i<= resx+2 ; i++) { // velocities with two layers of ghost cells
			udummy = 0.;
			fwrite(&udummy, sizeof(double), 1, u_up_conv);
			vdummy = 0.;
			fwrite(&vdummy, sizeof(double), 1, v_up_conv);
		}

		for (int i = 0 ; i <= resx+1 ; i++) { // pressure with only one layer of ghost cells

			pdummy = 0.;
			fwrite(&pdummy, sizeof(double), 1, p_up_conv);

		}

	}


	fclose(phi);
	fclose(u);
	fclose(v);
	fclose(p);

	fclose(temp);
	fclose(temps);

	fclose(u_lb);
	fclose(v_lb);

	fclose(u_up_conv);
	fclose(v_up_conv);
	fclose(p_up_conv);

	fclose(u_rt_conv);
	fclose(v_rt_conv);
	fclose(p_rt_conv);


	delete [] y_corns;

	return 0;
}









/************************************
 * **********************************
 * **ASSISTIVE - CONSTRUCTIVE FUNCS**
 * **********************************
 * *********************************/

// modified tanh func, devided by half so that it's not symmetric, and also to adapt to desirable height

double tanfunc (double a, int res, int j , double height) {

	double atanha = .5 * log((1+a)/(1-a)) ;

	double ksi = -1. + 2.*j/2./res ;

	return ((1./a) * tanh (ksi*atanha) + 1.) * height ;

}



// original tanh func, only adjusted to cover desirable height

double tanfunc_sym (double a, int res, int j , double height) {

	double atanha = .5 * log((1+a)/(1-a)) ;

	double ksi = -1. + 2.*j/res ;

	return ((1./a) * tanh (ksi*atanha) + 1.) * height /2. ;

}




//implementing the boundary layer u profile (integral analysis)
//for a certain Rex
//ATTENTION : NOT EFFICIENT IMPLEMENTATION - the Rex could be calculated
//out of this function // this way it is calculated for every j which is 
//useless- other than that, it doesn't affect the overall execution time
//considerably, and contributes to the better looking of the main func.
/*
double profile (double* y_cor, int j) {

	double x = 0.05;
	double mu = readdouble("../data/mu_ref.dat");
	double U = readdouble("../data/u_ref.dat");
	double l_ref = readdouble("../data/l_ref.dat");



	double delta, y;
	// delta is the boundary layer thickness

	delta = sqrt(30. * mu * x / rho / U );	
	y = (y_cor[j+1] + y_cor[j]) / 2. * l_ref;   //  multiply times (Lref) because this must be the dimensional value

	if (delta > y)
		return 2. * y / delta - pow ( y / delta , 2. );
	else 
		return 1. ;	

}
*/




