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

    //Initialize the Domain
    int Nx = 10, 
	Ny = 10, 
	Nz = 10;
    double Lx = 1.0, 
	   Ly = 1.0, 
	   Lz = 1.0;
    Domain *domain = new Domain(Nx, Ny, Nz, Lx, Ly, Lz);

    //Time Stepping info intialization
    TimeStepping::TimeSteppingType timeSteppingType = TimeStepping::CONST_CFL;
    double CFL 	    = 0.25;
    int maxTimeStep = 10;
    double maxTime  = 10.0;
    int filterStep  = 5;
    TimeStepping *timeStepping = new TimeStepping(timeSteppingType, CFL, maxTimeStep, maxTime, filterStep);

    //Boundary Condition Info
    BC::BCType bcXType = BC::PERIODIC_SOLVE; 
    BC::BCType bcYType = BC::PERIODIC_SOLVE; 
    BC::BCType bcZType = BC::PERIODIC_SOLVE; 

    BC::BCKind bcX0 = BC::PERIODIC;
    BC::BCKind bcX1 = BC::PERIODIC;
    BC::BCKind bcY0 = BC::PERIODIC;
    BC::BCKind bcY1 = BC::PERIODIC;
    BC::BCKind bcZ0 = BC::PERIODIC;
    BC::BCKind bcZ1 = BC::PERIODIC;

    BC *bc = new BC(bcXType, bcX0, bcX1,
		    bcYType, bcY0, bcY1,
		    bcZType, bcZ0, bcZ1);

    //Initialize the Solver
    double alphaF = 0.48;
    double mu_ref = 0.0001;
    CSolver *cSolver = new CSolver(domain, bc, timeStepping, alphaF, mu_ref); 
 
    return 0;
}
