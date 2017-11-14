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
    int    Nx = 128, 
	   Ny = 128, 
	   Nz = 128;
    double Lx = 1.0, 
	   Ly = 1.0, 
	   Lz = 1.0;
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

    return 0;
}









