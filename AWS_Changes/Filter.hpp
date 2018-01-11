#ifndef _FILTERH_
#define _FILTERH_

#include <math.h>
#include <cstring>
#include <iostream>
#include "Derivatives.hpp"
#include "BC.hpp"
#include "Utils.hpp"

class Filter{

    public:
	int Nx, Ny, Nz, N;

 	double alphaF;	
	double a0_8, a1_8, a2_8, a3_8, a4_8;
	double a0_6, a1_6, a2_6, a3_6;

	//explicit filtering coefficients at boundaries
	double a00, a01, a02, a03, a04;
	double b00, b01, b02, b03, b04;
	double c00, c01, c02, c03, c04;

	double *diagF, *offlowerF, *offupperF;
	
	Derivatives::Direct currentDir;
	BC::BCType bcType;


    Filter(double alphaF, Domain *dom, BC::BCType bcType, Derivatives::Direct currentDir){

	this->alphaF = alphaF;
	this->Nx = dom->Nx;
	this->Ny = dom->Ny;
	this->Nz = dom->Nz;

	this->currentDir = currentDir;

	this->bcType = bcType;

	if(currentDir == Derivatives::DIRX){
	    N = Nx;
	}else if(currentDir == Derivatives::DIRY){
	    N = Ny;
	}else if(currentDir == Derivatives::DIRZ){
	    N = Nz;
	}

	diagF = new double[N];
	offlowerF = new double[N];
	offupperF = new double[N];

	a0_8     = (93.0 + 70.0*alphaF)/128.0;
        a1_8     = ( 7.0 + 18.0*alphaF)/16.0;
        a2_8     = (-7.0 + 14.0*alphaF)/32.0;
        a3_8     = ( 1.0 -  2.0*alphaF)/16.0;
        a4_8     = (-1.0 +  2.0*alphaF)/128.0;

	a0_6     = (11.0 + 10.0*alphaF)/16.0;
	a1_6     = (15.0 + 34.0*alphaF)/32.0;
	a2_6     = (-3.0 +  6.0*alphaF)/16.0;
	a3_6	 = ( 1.0 -  2.0*alphaF)/32.0;	

	a00 = 15.0/16.0;	    
	a01 =  4.0/16.0;
	a02 = -6.0/16.0;
	a03 =  4.0/16.0;
	a04 = -1.0/16.0;

	b00 =  1.0/16.0;
	b01 =  3.0/4.0;
	b02 =  6.0/16.0;
	b03 = -4.0/16.0;
	b04 =  1.0/16.0;

	c00 = -1.0/16.0;
	c01 =  4.0/16.0;
	c02 =  5.0/8.0;
	c03 =  4.0/16.0;
	c04 = -1.0/16.0;

        for(int ip = 0; ip < N; ip++){
            diagF[ip] = 1.0;
            offlowerF[ip]  = alphaF;
            offupperF[ip]  = alphaF;
        }
        
	if(bcType == BC::DIRICHLET_SOLVE){
	    //explicitly filter the boundary points (Lele)
	    offlowerF[0] = 0.0;
	    offlowerF[1] = 0.0;
	    offlowerF[2] = 0.0;
	    offupperF[0] = 0.0;
	    offupperF[1] = 0.0;
	    offupperF[2] = 0.0;

	    offlowerF[N-1] = 0.0;
	    offlowerF[N-2] = 0.0;
	    offlowerF[N-3] = 0.0;
	    offupperF[N-1] = 0.0;
	    offupperF[N-2] = 0.0;
	    offupperF[N-3] = 0.0;
	}

    }

    void multRHSPeriodicFilter(double *phi, double *RHSvec);
    void multRHSDirichletFilter(double *phi, double *RHSvec);

    void FilterPeriodic(double *phi, double *phiF);
    void FilterDirichlet(double *phi, double *phiF);

    void compactFilter(double *phi, double *phiF);
    void filterField(double *dataIn, double *dataOut);

};

#endif
