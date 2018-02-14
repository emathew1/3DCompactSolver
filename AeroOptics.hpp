#ifndef _AEROOPTICSH_
#define _AEROOPTICSH_

#include <iostream>
#include "Macros.hpp"
#include "Domain.hpp"
#include "IdealGas.hpp"

//This class is for operations in computing the Aero-Optic distortion in the flow
//Right now, hard-coded for integrating in the Y direction of the domain
//and rotation of the beam grid along the Z-axis

class AeroOptics{

    public:

	Domain *domain; //domain object from solver
	IdealGas *idealGas; //ideal gas object from solver

	int Nx, Ny, Nz, N; //Dimensions of the grid
	double dx, dy, dz, ds;

	double *rho, *p, *T; //density, pressure, and temperature pointers from solver

	int numAngles;
	double maxAngle, lowerIBound, upperIBound, totalDY;

	AeroOptics(Domain *domain, IdealGas *idealGas, *rho, *p, *T, numAngles, maxAngle){

	    this->domain   = domain;
	    this->idealGas = idealGas;

	    this->Nx = domain->Nx;
	    this->Ny = domain->Ny;
	    this->Nz = domain->Nz;
	    N = Nx*Ny*Nz;

	    this->dx = domain->dx;
	    this->dy = domain->dy;
	    this->dz = domain->dz;

	    this->rho = rho;
	    this->p   = p;
	    this->T   = T;

	    this->numAngles = numAngles;
	    this->maxAngle  = maxAngle;

	    //These should be outside of the sponge zone of the domain
	    lowerIBound = 0.15*domain->Ly;
	    upperIBound = 0.85*domain->Ly;
	    totalDY = upperIBound - lowerIBound;

	}

};

#endif
