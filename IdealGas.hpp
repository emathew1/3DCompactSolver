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

    inline double solveMu(double T){
        return mu_ref*pow(T/T_ref, 0.76);
    }

    inline double solveAmu(double T){
        return 0.76*(mu_ref/pow(T_ref, 0.76))*pow(T, 0.76-1.0);
    } 

    inline double solveU(double rho, double rhoU){
        return rhoU/rho;
    }

    inline double solvep(double rho, double rhoE, double U, double V, double W){
        return (gamma-1)*(rhoE - 0.5 * rho*(U*U + V*V + W*W));
    }

    inline double solveT(double rho, double p){
        return  p/(rho*R_gas);
    }

    inline double solveSOS(double rho, double p){
        return  sqrt(gamma*p/rho);
    }

    inline double solverhoE(double rho, double p, double U, double V, double W){
        return p/(gamma-1.0) + 0.5*rho*(U*U + V*V + W*W);
    }

};

#endif
