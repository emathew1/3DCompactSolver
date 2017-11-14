#include <math.h>
#include <cstring>
#include <iostream>
#include "Utils.hpp"

class Filter{

    public:
	int Nx, Ny;

 	double alphaF;	
	double *alphaFx;
	double *alphaFy;
	double a0, a1, a2, a3, a4, a5;
	double b0, b1, b2, b3, b4, b5, b6, b7, b8;
	double c0, c1, c2, c3, c4, c5, c6, c7, c8;
	double d0, d1, d2, d3, d4, d5, d6, d7, d8;
	double e0, e1, e2, e3, e4, e5, e6, e7, e8;


	double *diagFx, *offFx1, *offFx2;
	double *diagFy, *offFy1, *offFy2;

    Filter(){
	alphaF  = 0.49;
	alphaFx = NULL;
	alphaFy = NULL;

	a0 = 0.0; 
        a1 = 0.0;
        a2 = 0.0;
        a3 = 0.0;
        a4 = 0.0;
        a5 = 0.0;

	b0 = 0.0; c0 = 0.0; d0 = 0.0; e0 = 0.0;
	b1 = 0.0; c1 = 0.0; d1 = 0.0; e1 = 0.0;
	b2 = 0.0; c2 = 0.0; d2 = 0.0; e2 = 0.0;
	b3 = 0.0; c3 = 0.0; d3 = 0.0; e3 = 0.0;
	b4 = 0.0; c4 = 0.0; d4 = 0.0; e4 = 0.0;
	b5 = 0.0; c5 = 0.0; d5 = 0.0; e5 = 0.0;
	b6 = 0.0; c6 = 0.0; d6 = 0.0; e6 = 0.0;
	b7 = 0.0; c7 = 0.0; d7 = 0.0; e7 = 0.0;
	b8 = 0.0; c8 = 0.0; d8 = 0.0; e8 = 0.0;

	diagFx = NULL; offFx1 = NULL; offFx2 = NULL;
	diagFy = NULL; offFy1 = NULL; offFy2 = NULL;

	Nx = 0; Ny = 0;
    }


    Filter(int NX, int NY){

	Nx = NX; Ny = NY;

	alphaFx = new double[Nx];
  	alphaFy = new double[Ny];

	alphaF = 0.49;
	for(int ip = 0; ip < Nx; ip++){
	    alphaFx[ip] = alphaF;
	}

	for(int ip = 0; ip < Ny; ip++){
	    alphaFy[ip] = alphaF;
	}

	a0     = (93.0 + 70.0*alphaF)/128.0;
        a1     = ( 7.0 + 18.0*alphaF)/16.0;
        a2     = (-7.0 + 14.0*alphaF)/32.0;
        a3     = ( 1.0 -  2.0*alphaF)/16.0;
        a4     = (-1.0 +  2.0*alphaF)/128.0;
        a5     = 0.0;

	b0 = 0.0; c0 = 0.0; d0 = 0.0; e0 = 0.0;
	b1 = 0.0; c1 = 0.0; d1 = 0.0; e1 = 0.0;
	b2 = 0.0; c2 = 0.0; d2 = 0.0; e2 = 0.0;
	b3 = 0.0; c3 = 0.0; d3 = 0.0; e3 = 0.0;
	b4 = 0.0; c4 = 0.0; d4 = 0.0; e4 = 0.0;
	b5 = 0.0; c5 = 0.0; d5 = 0.0; e5 = 0.0;
	b6 = 0.0; c6 = 0.0; d6 = 0.0; e6 = 0.0;
	b7 = 0.0; c7 = 0.0; d7 = 0.0; e7 = 0.0;
	b8 = 0.0; c8 = 0.0; d8 = 0.0; e8 = 0.0;

	diagFx = new double[NX]; 
	offFx1  = new double[NX];
	offFx2  = new double[NX];
	diagFy = new double[NY]; 
	offFy1  = new double[NY];
	offFy2  = new double[NY];

        for(int ip = 0; ip < Nx; ip++){
            diagFx[ip] = 1.0;
            offFx1[ip]  = alphaFx[ip];
            offFx2[ip]  = alphaFx[ip];
        }
        
        for(int ip = 0; ip < Ny; ip++){
            diagFy[ip] = 1.0;
            offFy1[ip]  = alphaFy[ip];
            offFy2[ip]  = alphaFy[ip];
        }
	
    }

    void calcFilterCoefficients(double *alphaFdir);

    void multRHSFilter(double *phi, int N, double *RHSvec);
    void FilterPeriodicY(double *phi, double *phiF);
    void FilterPeriodicX(double *phi, double *phiF);

    void multRHSFilterFiniteDomain(double *phi, int N, double *RHSvec);
    void FilterFiniteDomainY(double *phi, double *phiF);
    void FilterFiniteDomainX(double *phi, double *phiF);

};
