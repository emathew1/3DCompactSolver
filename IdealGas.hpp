#ifndef _IDEALGASH_
#define _IDEALGASH_

#include "Macros.hpp"
#include "Domain.hpp"
#include <math.h>

class IdealGas{

    public:
	
	int Nx, Ny, Nz;
	Domain *domain;
	double gamma;
	double Pr;
	double R_gas;
	double cp;
	double p_ref;
	double rho_ref;
	double T_ref;
	double mu_ref;

    IdealGas(){
	Nx = 0; Ny = 0; Nz = 0;
	gamma = 1.4;
	Pr    = 0.7;
	p_ref = 1.0/gamma;
	rho_ref = 1.0;
	T_ref = 1.0;
	mu_ref = 1.0;
	R_gas = p_ref/rho_ref/T_ref;
	cp = R_gas*gamma/(gamma - 1.0);
    }

    IdealGas(Domain *domain){
	this->Nx = domain->Nx; 
	this->Ny = domain->Ny;
	this->Nz = domain->Nz;
	gamma = 1.4;
	Pr    = 0.7;
	p_ref = 1.0/gamma;
	rho_ref = 1.0;
	T_ref = 1.0;
	mu_ref = 1.0;
	R_gas = p_ref/rho_ref/T_ref;
	cp = R_gas*gamma/(gamma - 1.0);
    } 

    void solveMu(double *T, double *mu);
    void solveAmu(double *T, double *Amu);
    void solveU(double *rho, double *rhoU, double *U);
    void solvep(double *rho, double *rhoE, double *U, double *V, double *W, double *p);
    void solveT(double *rho, double *p, double *T);
    void solveSOS(double *rho, double *p, double *SOS); 
    void solverhoE(double *rho, double *p, double *U, double *V, double *W, double *rhoE); 

};

#endif
