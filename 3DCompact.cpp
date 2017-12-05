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
    int    Nx = 32, 
	   Ny = 32, 
	   Nz = 32;
    double Lx = 0.1, 
	   Ly = 0.1, 
	   Lz = 0.1;
    Domain *dom = new Domain(Nx, Ny, Nz, Lx, Ly, Lz);

    ////////////////////////////////////
    //Time Stepping info intialization//
    ////////////////////////////////////
    TimeStepping::TimeSteppingType timeSteppingType = TimeStepping::CONST_CFL;
    double CFL 	     = 0.1;
    int maxTimeStep  = 100;
    double maxTime   = 10.0;
    int filterStep   = 1000;
    int checkStep    = 1;
    int dumpStep     = 10;
    TimeStepping *ts = new TimeStepping(timeSteppingType, CFL, maxTimeStep, maxTime, filterStep, checkStep, dumpStep);


    ///////////////////////////
    //Boundary Condition Info//
    ///////////////////////////
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



    /////////////////////////
    //Initialize the Solver//
    /////////////////////////
    double alphaF = 0.40;
    double mu_ref = 1.0000;
    CSolver *cs   = new CSolver(dom, bc, ts, alphaF, mu_ref); 


    ///////////////////////////////
    //Set flow initial conditions//
    ///////////////////////////////
    FOR_Z{
	FOR_Y{
	    FOR_X{
		int ii = GET3DINDEX_XYZ;
		cs->U0[ii]   = 0.00;
		cs->V0[ii]   = 0.00;
		cs->W0[ii]   = 0.00;
	
//		if(cs->dom->x[i] > 0.5){
//		    cs->rho0[ii] = 0.125;
//		    cs->p0[ii]   = 0.1/cs->ig->gamma;
//		}else{
	  	    cs->rho0[ii] = 1.0;
		    cs->p0[ii]   = 1.0/cs->ig->gamma;
//		}
	    }
	}
    }

    cs->setInitialConditions();

    
 
    while(cs->endFlag == false){

	//Get the dt for this time step
        cs->calcDtFromCFL();

	//start rkStep 1
        cs->rkStep = 1;

	cs->preStepBCHandling();
	
	cs->preStepDerivatives();

	getRange(cs->momYEulerX, "momYEulerX", Nx, Ny, Nz);
	getRange(cs->momYEulerY, "momYEulerY", Nx, Ny, Nz);
	getRange(cs->momYEulerZ, "momYEulerZ", Nx, Ny, Nz);
	getRange(cs->momXEulerX, "momXEulerX", Nx, Ny, Nz);
	getRange(cs->momXEulerY, "momXEulerY", Nx, Ny, Nz);
	getRange(cs->momXEulerZ, "momXEulerZ", Nx, Ny, Nz);


	cs->solveContinuity();
	cs->solveXMomentum();
	cs->solveYMomentum();
	cs->solveZMomentum();
	cs->solveEnergy();

	cs->postStepBCHandling();

	cs->updateConservedData();
	cs->updateNonConservedData();


	//start rkStep 2
        cs->rkStep = 2;

	cs->preStepBCHandling();
	
	cs->preStepDerivatives();

	cs->solveContinuity();
	cs->solveXMomentum();
	cs->solveYMomentum();
	cs->solveZMomentum();
	cs->solveEnergy();

	cs->postStepBCHandling();

	cs->updateConservedData();
	cs->updateNonConservedData();

	getRange(cs->rhoE2, "rhoE2", Nx, Ny, Nz);
	//start rkStep 3
        cs->rkStep = 3;

	cs->preStepBCHandling();
	
	cs->preStepDerivatives();

	cs->solveContinuity();
	cs->solveXMomentum();
	cs->solveYMomentum();
	cs->solveZMomentum();
	cs->solveEnergy();

	cs->postStepBCHandling();

	cs->updateConservedData();
	cs->updateNonConservedData();

	getRange(cs->rhoE2, "rhoE2", Nx, Ny, Nz);
	//start rkStep 4
        cs->rkStep = 4;

	cs->preStepBCHandling();
	
	cs->preStepDerivatives();

	cs->solveContinuity();
	cs->solveXMomentum();
	cs->solveYMomentum();
	cs->solveZMomentum();
	cs->solveEnergy();

	cs->postStepBCHandling();

	cs->updateConservedData();

	getRange(cs->rhoE2, "rhoE2", Nx, Ny, Nz);
	getRange(cs->rhoE1, "rhoE1", Nx, Ny, Nz);
	//on rk step 4 we'll filter the data if need be first...
	cs->filterConservedData();
	getRange(cs->rhoE1, "rhoE1", Nx, Ny, Nz);

	//Then we'll update the nonconserved variables...
	cs->updateNonConservedData();

	//Update the sponge if using it...
	cs->updateSponge();

	//Now lets output the solution info to the console...
	cs->checkSolution();

	//Dump the conserved variables if its the time to do so...
	cs->dumpSolution();

	//Check if we've met our end conditions yet...
	cs->checkEnd();

    }




    return 0;
}









