#ifndef _AEROOPTICSH_
#define _AEROOPTICSH_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include "Macros.hpp"
#include "Utils.hpp"
#include "Domain.hpp"
#include "IdealGas.hpp"

//This class is for operations in computing the Aero-Optic distortion in the flow
//Right now, hard-coded for integrating in the Y direction of the domain
//and rotation of the beam grid along the Z-axis

using namespace std;

class AeroOptics{

    public:

	Domain *domain; //domain object from solver
	IdealGas *idealGas; //ideal gas object from solver

	int Nx, Ny, Nz, N; //Dimensions of the grid
	double dx, dy, dz, ds;

	int Nline;
 	double dxLine, dyLine, dsLine;

	double *rho, *p, *T; //density, pressure, and temperature pointers from solver

	int numAngles;
	double maxAngle, lowerIBound, upperIBound;
	double totalDX, totalDY, Ldiag;

	double *OPL;

	AeroOptics(Domain *domain, IdealGas *idealGas, double *rho, double *p, double *T, int numAngles, double maxAngle){

	    this->domain   = domain;
	    this->idealGas = idealGas;

	    this->Nx = domain->Nx;
	    this->Ny = domain->Ny;
	    this->Nz = domain->Nz;
	    N = Nx*Ny*Nz;

	    this->dx = domain->dx;
	    this->dy = domain->dy;
	    this->dz = domain->dz;
	    ds = fmin(dx,fmin(dy,dz));

	    this->rho = rho;
	    this->p   = p;
	    this->T   = T;

	    this->numAngles = numAngles;
	    this->maxAngle  = maxAngle;

	    //These should be outside of the sponge zone of the domain
	    lowerIBound = 0.2*domain->Ly;
	    upperIBound = 0.8*domain->Ly;
	    totalDY = upperIBound - lowerIBound;

	    //Allocate the OPL container...
	    OPL = new double[numAngles*Nx*Nz]; 

	}

        void generateBaseBeamDiff(int angleNumber);
	void writeOPLFiles();
	void computeAO();
};

#endif
