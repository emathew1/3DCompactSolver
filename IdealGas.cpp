#include "IdealGas.hpp"

void IdealGas::solveMu(double *T, double *mu){
    #pragma omp parallel for
    FOR_XYZ{
        mu[ip] = mu_ref*pow(T[ip]/T_ref, 0.76);
    }
}

void IdealGas::solveAmu(double *T, double *Amu){
    #pragma omp parallel for
    FOR_XYZ{
	Amu[ip] = 0.76*(mu_ref/pow(T_ref, 0.76))*pow(T[ip], 0.76-1.0);
    }
}

void IdealGas::solveU(double *rho, double *rhoU, double *U){
    #pragma omp parallel for
    FOR_XYZ{
        U[ip] = rhoU[ip]/rho[ip];
    }
}

void IdealGas::solvep(double *rho, double *rhoE, double *U, double *V, double *W, double *p){
    #pragma omp parallel for
    FOR_XYZ{
        p[ip] = (gamma-1)*(rhoE[ip] - 0.5 * rho[ip]*(U[ip]*U[ip] + V[ip]*V[ip] + W[ip]*W[ip]));
    }
}

void IdealGas::solveT(double *rho, double *p, double *T){
    #pragma omp parallel for
    FOR_XYZ{
        T[ip] = p[ip]/(rho[ip]*R_gas);
    }
}

void IdealGas::solveSOS(double *rho, double *p, double *SOS){
    #pragma omp parallel for
    FOR_XYZ{
        SOS[ip] = sqrt(gamma*p[ip]/rho[ip]);
    }
}

void IdealGas::solverhoE(double *rho, double *p, double *U, double *V, double *W, double *rhoE){
    #pragma omp parallel for
    FOR_XYZ{
        rhoE[ip] = p[ip]/(gamma-1.0) + 0.5*rho[ip]*(U[ip]*U[ip] + V[ip]*V[ip] + W[ip]*W[ip]);
    }
}

