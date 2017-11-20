#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <chrono>

#include "Macros.hpp"
#include "Utils.hpp"
#include "BC.hpp"
#include "TimeStepping.hpp"
#include "CSolver.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]){

    cout << endl;
    cout << "----------------------------------" << endl;
    cout << " 3D 6th-Order Compressible Solver " << endl;
    cout << "----------------------------------" << endl;
    cout << endl;


    /////////////////////////
    //Initialize the Domain//
    /////////////////////////
    int    Nx = 32, 
	   Ny = 32, 
	   Nz = 32;
    double Lx = 2.0*M_PI*((double)Nx - 1.0)/(double(Nx)), 
	   Ly = 2.0*M_PI*((double)Ny - 1.0)/(double(Ny)), 
	   Lz = 2.0*M_PI*((double)Nz - 1.0)/(double(Nz));
    Domain *dom = new Domain(Nx, Ny, Nz, Lx, Ly, Lz);



    ////////////////////////////////////
    //Time Stepping info intialization//
    ////////////////////////////////////
    TimeStepping::TimeSteppingType timeSteppingType = TimeStepping::CONST_CFL;
    double CFL 	     = 0.25;
    int maxTimeStep  = 10;
    double maxTime   = 10.0;
    int filterStep   = 5;
    TimeStepping *ts = new TimeStepping(timeSteppingType, CFL, maxTimeStep, maxTime, filterStep);



    ///////////////////////////
    //Boundary Condition Info//
    ///////////////////////////
    BC::BCType bcXType = BC::PERIODIC_SOLVE; 
    BC::BCType bcYType = BC::DIRICHLET_SOLVE; 
    BC::BCType bcZType = BC::PERIODIC_SOLVE; 

    BC::BCKind bcX0 = BC::PERIODIC;
    BC::BCKind bcX1 = BC::PERIODIC;
    BC::BCKind bcY0 = BC::SPONGE;
    BC::BCKind bcY1 = BC::SPONGE;
    BC::BCKind bcZ0 = BC::PERIODIC;
    BC::BCKind bcZ1 = BC::PERIODIC;

    BC *bc = new BC(bcXType, bcX0, bcX1,
		    bcYType, bcY0, bcY1,
		    bcZType, bcZ0, bcZ1);



    /////////////////////////
    //Initialize the Solver//
    /////////////////////////
    double alphaF = 0.48;
    double mu_ref = 0.0001;
    CSolver *cs   = new CSolver(dom, bc, ts, alphaF, mu_ref); 

    cs->initializeSolverData(); 


    ///////////////////////////////
    //Set flow initial conditions//
    ///////////////////////////////
    FOR_XYZ{
	cs->rho0[ip] = 1.0;
	cs->U0[ip]   = 0.1;
	cs->V0[ip]   = 0.1;
	cs->W0[ip]   = 0.1;
	cs->p0[ip]   = 1.0/cs->ig->gamma;
    }

    cs->setInitialConditions();
 
    cs->calcDtFromCFL();

    double *test = new double[Nx];
    double *dtest = new double[Nx];
    FOR_X{
	test[i] = sin(cs->dom->x[i]);
	dtest[i] = 0.0;
    }
  
    cs->filtX->compactFilter(test, dtest);

    FOR_X{
	cout << test[i] << " " << dtest[i] << endl;
    }

    FOR_Y{
	test[j] = sin(cs->dom->y[j]);
	dtest[j] = 0.0;
    }
  
    cs->filtY->compactFilter(test, dtest);
    cout << endl;
    FOR_Y{
	cout << test[j] << " " << dtest[j] << endl;
    }



    return 0;
}









