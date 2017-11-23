#ifndef _DERIVATIVESH_
#define _DERIVATIVESH_

#include <math.h>
#include <cstring>
#include <iostream>
#include "Macros.hpp"
#include "Domain.hpp"
#include "BC.hpp"
#include "Utils.hpp"

class Derivatives{


    public:

	int Nx, Ny, Nz, N;
	int i, j, k;
	double dx, dy, dz, dd;

        //First Derivative

        // Interior Coefficients	
	double alpha_1D;
	double a_1D, b_1D;

	//T6 Dirichlet Coefficients w/ T4 at edge
	double alpha11_1D, alpha21_1D, alpha22_1D;
	double a1_1D, b1_1D, c1_1D, d1_1D, e1_1D, f1_1D;
	double a2_1D, b2_1D, c2_1D, d2_1D, e2_1D; 

	double *diag_1D, *offlower_1D, *offupper_1D;


        //Second Derivative

        // Interior Coefficients	
	double alpha_2D;
	double a_2D, b_2D;

	//T6 Dirichlet Coefficients w/ T4 at edge
	double alpha11_2D, alpha21_2D, alpha22_2D;
	double a1_2D, b1_2D, c1_2D, d1_2D, e1_2D;
	double a2_2D, b2_2D, c2_2D, d2_2D, e2_2D, f2_2D; 

	double *diag_2D, *offlower_2D, *offupper_2D;
	enum Direct {DIRX, DIRY, DIRZ};
	Direct currentDir;

	BC::BCType bcType;

	//Constructor
        Derivatives(Domain *dom, BC::BCType bcType, Direct currentDir){

   	    this->Nx = dom->Nx;
	    this->Ny = dom->Ny;
	    this->Nz = dom->Nz;
	    this->dx = dom->dx; 
	    this->dy = dom->dy;
	    this->dz = dom->dz;

	    this->currentDir = currentDir;

	    this->bcType = bcType;

	    if(currentDir == DIRX){
		N  = Nx;
		dd = dx;
	    }else if(currentDir == DIRY){
		N  = Ny;
		dd = dy;
	    }else if(currentDir == DIRZ){
		N  = Nz;
		dd = dz; 
	    }

	    //1st Derivative coefficients
	    alpha_1D = 1.0/3.0;
	    a_1D     = 14.0/9.0;
	    b_1D     = 1.0/9.0;

   	    alpha11_1D = 5.0;
	    a1_1D = -197.0/60.0;
	    b1_1D =   -5.0/12.0;
	    c1_1D =    5.0;
	    d1_1D =   -5.0/3.0;
 	    e1_1D =    5.0/12.0;
	    f1_1D =   -1.0/20.0;

	    alpha21_1D = 1.0/8.0;
	    alpha22_1D = 3.0/4.0;
	    a2_1D = -43.0/96.0;
	    b2_1D = -5.0/6.0;
	    c2_1D =  9.0/8.0;
	    d2_1D =  1.0/6.0;
	    e2_1D = -1.0/96.0;

	    //2nd Derivative coefficients
	    alpha_2D = 2.0/11.0;
	    a_2D     = 12.0/11.0;
	    b_2D     = 3.0/11.0;

	    //4th order here
   	    alpha11_2D = 10.0;
	    a1_2D =  145.0/12.0;
	    b1_2D =  -76.0/3.0;
	    c1_2D =   29.0/2.0;
	    d1_2D =   -4.0/3.0;
 	    e1_2D =    1.0/12.0;

	    //6th order here...
	    alpha21_2D = 2.0/11.0;
	    alpha22_2D = -131.0/22.0;
	    a2_2D =  177.0/88.0;
	    b2_2D = -507.0/44.0;
	    c2_2D =  783.0/44.0;
	    d2_2D = -201.0/22.0;
	    e2_2D =  81.0/88.0;
	    f2_2D =  -3.0/44.0;


	    diag_1D     = new double[N]; 
	    offlower_1D = new double[N];
	    offupper_1D = new double[N];

	    diag_2D     = new double[N]; 
	    offlower_2D = new double[N];
	    offupper_2D = new double[N];

	//Left off here!     
  
	for(int ip = 0; ip < N; ip++){ //Not correct...
	    diag_1D[ip] = 1.0;
	    offlower_1D[ip]  = alpha_1D;
	    offupper_1D[ip]  = alpha_1D;
	}

	for(int ip = 0; ip < N; ip++){
	    diag_2D[ip] = 1.0;
	    offlower_2D[ip]  = alpha_2D;
	    offupper_2D[ip]  = alpha_2D;
	}

	if(bcType == BC::DIRICHLET_SOLVE){
	    offupper_1D[0] = alpha11_1D;  
	    offupper_1D[1] = alpha22_1D;
	    offlower_1D[1] = alpha21_1D;

	    offlower_1D[N-1] = alpha11_1D;
	    offlower_1D[N-2] = alpha22_1D;
	    offupper_1D[N-2] = alpha21_1D;

	    offupper_2D[0] = alpha11_2D;  
	    offupper_2D[1] = alpha22_2D;
	    offlower_2D[1] = alpha21_2D;

	    offlower_2D[N-1] = alpha11_2D;
	    offlower_2D[N-2] = alpha22_2D;
	    offupper_2D[N-2] = alpha21_2D;
	}	


    }

    //Function's to call...
    void calc1stDerivField(double *dataIn, double *dataOut);
    void calc2ndDerivField(double *dataIn, double *dataOut);

    void calc1stDeriv(double *phi, double *dphi);
    void calc2ndDeriv(double *phi, double *dphi);

    //Need a cleaner way of passing these things...
    void multRHS1stDerivPeriodic(double dh, double *phi, int N, double *RHSvec);
    void multRHS2ndDerivPeriodic(double dh, double *phi, int N, double *RHSvec);
    void multRHS1stDerivDirichlet(double dh, double *phi, int N, double *RHSvec);
    void multRHS2ndDerivDirichlet(double dh, double *phi, int N, double *RHSvec);


    void Compact1stPeriodic(double *phi, double *dphidx);
    void Compact2ndPeriodic(double *phi, double *dphidx);
    void Compact1stDirichlet(double *phi, double *dphidx);
    void Compact2ndDirichlet(double *phi, double *dphidx);
};

#endif
