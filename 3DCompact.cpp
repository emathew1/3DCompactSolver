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
    double Lx = 1.0, 
	   Ly = 1.0, 
	   Lz = 1.0;
    Domain *dom = new Domain(Nx, Ny, Nz, Lx, Ly, Lz);

    ////////////////////////////////////
    //Time Stepping info intialization//
    ////////////////////////////////////
    TimeStepping::TimeSteppingType timeSteppingType = TimeStepping::CONST_CFL;
    double CFL 	     = 0.25;
    int maxTimeStep  = 10000;
    double maxTime   = 1.0;
    int filterStep   = 5;
    int checkStep    = 1;
    int dumpStep     = 1000;
    TimeStepping *ts = new TimeStepping(timeSteppingType, CFL, maxTimeStep, maxTime, filterStep, checkStep, dumpStep);


    ///////////////////////////
    //Boundary Condition Info//
    ///////////////////////////
    BC::BCType bcXType = BC::DIRICHLET_SOLVE; 
    BC::BCType bcYType = BC::DIRICHLET_SOLVE; 
    BC::BCType bcZType = BC::DIRICHLET_SOLVE; 

    BC::BCKind bcX0 = BC::ADIABATIC_WALL;
    BC::BCKind bcX1 = BC::ADIABATIC_WALL;
    BC::BCKind bcY0 = BC::ADIABATIC_WALL;
    BC::BCKind bcY1 = BC::ADIABATIC_WALL;
    BC::BCKind bcZ0 = BC::ADIABATIC_WALL;
    BC::BCKind bcZ1 = BC::ADIABATIC_WALL;

    BC *bc = new BC(bcXType, bcX0, bcX1,
		    bcYType, bcY0, bcY1,
		    bcZType, bcZ0, bcZ1);



    /////////////////////////
    //Initialize the Solver//
    /////////////////////////
    double alphaF = 0.48;
    double mu_ref = 0.00001;
    int blocksize = 16;
    CSolver *cs   = new CSolver(dom, bc, ts, alphaF, mu_ref, blocksize); 


    ///////////////////////////////
    //Set flow initial conditions//
    ///////////////////////////////
    FOR_Z{
	FOR_Y{
	    FOR_X{
		int ii = GET3DINDEX_XYZ;
		cs->U0[ii]   = 0.0;//sin(cs->dom->x[i]);
		cs->V0[ii]   = 0.0;//sin(cs->dom->y[j]);
		cs->W0[ii]   = 0.0;//sin(cs->dom->z[k]);
	
//		if(cs->dom->x[i] > 0.5){
//		    cs->rho0[ii] = 0.125;
// 		    cs->p0[ii]   = 0.1/cs->ig->gamma;
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

	//on rk step 4 we'll filter the data if need be first...
	cs->filterConservedData();

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









