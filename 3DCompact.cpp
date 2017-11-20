#include <omp.h>
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

    #pragma omp parallel
    {
      int threadCount = omp_get_num_threads();
      int threadID    = omp_get_thread_num();
      if(threadCount > 1 && threadID == 0){
          cout << "Running with OpenMP using " << threadCount << " threads" << endl;
      }
    }


    /////////////////////////
    //Initialize the Domain//
    /////////////////////////
    int    Nx = 100, 
	   Ny = 100, 
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
    double alphaF = 0.49;
    double mu_ref = 0.0001;
    CSolver *cs   = new CSolver(dom, bc, ts, alphaF, mu_ref); 



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

    double test[Nx*1000];
    double dtest[Nx*1000];
    FOR_X{
	for(int ip = 0; ip < 1000; ip++){
	    test[ip*Nx+i] = sin(cs->dom->x[i]);
	    dtest[ip*Nx+i] = 0.0;
	}
    }

    auto  t1 = std::chrono::system_clock::now();
    #pragma omp parallel for 
    for(int jp = 0; jp < 1000; jp++){
	double *testpoint = &test[jp*Nx];
	double *dtestpoint = &dtest[jp*Nx];
        cs->derivX->calc1stDeriv(testpoint, dtestpoint);
    }
    auto  t2 = std::chrono::system_clock::now();
    auto  t3 = t2-t1;

    cout << "1st derivative:  " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000/(double)1000 << endl;

    t1 = std::chrono::system_clock::now();
    for(int jp = 0; jp < 1000; jp++){
	double *testpoint = &test[jp*Nx];
	double *dtestpoint = &dtest[jp*Nx];
        cs->derivX->calc1stDeriv(testpoint, dtestpoint);
    }
    t2 = std::chrono::system_clock::now();
    cout << "1st derivative:  " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000/(double)1000 << endl;


/*
    FOR_X{
	cout << test[i] << " " << dtest[i] << endl;
    }
*/
    double testy[Ny*1000];
    double dtesty[Ny*1000];

    FOR_Y{
	for(int ip = 0; ip < 1000; ip++){
	    testy[ip*Ny + j] = sin(cs->dom->y[j]);
	    dtesty[ip*Ny + j] = 0.0;
	}
    }

    t1 = std::chrono::system_clock::now();
    #pragma omp parallel for 
    for(int jp = 0; jp < 1000; jp++){ 
	double *testypoint = &testy[jp*Ny]; 
	double *dtestypoint = &dtesty[jp*Ny]; 
        cs->derivY->calc2ndDeriv(testypoint, dtestypoint);
    }
    t2 = std::chrono::system_clock::now();
    cout << "2nd Derivative: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000/(double)1000 << endl;
/*
    cout << endl;
    FOR_Y{
	cout << testy[j] << " " << dtesty[j] << endl;
    }
*/

    return 0;
}









