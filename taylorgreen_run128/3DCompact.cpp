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
    int    Nx = 128, 
	   Ny = 128, 
	   Nz = 128;
    double Lx = 2.0*M_PI - (2.0*M_PI/((double)(Nx+1.0))), 
	   Ly = 2.0*M_PI - (2.0*M_PI/((double)(Ny+1.0))), 
	   Lz = 2.0*M_PI - (2.0*M_PI/((double)(Nz+1.0)));
    Domain *dom = new Domain(Nx, Ny, Nz, Lx, Ly, Lz);

    ////////////////////////////////////
    //Time Stepping info intialization//
    ////////////////////////////////////
    TimeStepping::TimeSteppingType timeSteppingType = TimeStepping::CONST_CFL;
    double CFL 	     = 0.25;
    int maxTimeStep  = 25000;
    double maxTime   = 25.0;
    int filterStep   = 50;
    int checkStep    = 1;
    int dumpStep     = 500;
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
    double alphaF = 0.495;
    double mu_ref = 0.0000188745;
    int blocksize = 16;
    CSolver *cs   = new CSolver(dom, bc, ts, alphaF, mu_ref, blocksize); 


    /////////////////////////////
    //load in turbulence output//
    /////////////////////////////
/*    ifstream uFile, vFile, wFile;
    uFile.open("U_Mt0p3_N128_k8.dat",ifstream::in);
    vFile.open("V_Mt0p3_N128_k8.dat",ifstream::in);
    wFile.open("W_Mt0p3_N128_k8.dat",ifstream::in);

    double *u_temp = new double[Nx*Ny*Nz];
    double *v_temp = new double[Nx*Ny*Nz];
    double *w_temp = new double[Nx*Ny*Nz];

    FOR_Z{
	FOR_Y{
	    FOR_X{
		int ii = GET3DINDEX_XYZ;
		uFile >> u_temp[ii];
		vFile >> v_temp[ii];
		wFile >> w_temp[ii];
	    }
	}
    }

    uFile.close();
    vFile.close();
    wFile.close();
*/

    ///////////////////////////////
    //Set flow initial conditions//
    ///////////////////////////////
    double u0 = 1.0;
    FOR_Z{
	FOR_Y{
	    FOR_X{
		int ii = GET3DINDEX_XYZ;

		cs->U0[ii]   =  u0*sin(cs->dom->x[i])*cos(cs->dom->y[j])*cos(cs->dom->z[k]);
		cs->V0[ii]   = -u0*cos(cs->dom->x[i])*sin(cs->dom->y[j])*cos(cs->dom->z[k]);
		cs->W0[ii]   =  0.0; 
		cs->rho0[ii] = 1.0;
		cs->p0[ii]   = (100.0 + (1.0/16.0)*(cos(cs->dom->x[i]*2.0) + cos(cs->dom->y[j]*2.0))*(cos(cs->dom->z[k]*2.0) + 2.0))/cs->ig->gamma;
		
/*
		cs->rho0[ii] = 1.0;
		cs->p0[ii]   = 100.0/cs->ig->gamma;
		cs->U0[ii]   = u_temp[ii];
		cs->V0[ii]   = v_temp[ii];
		cs->W0[ii]   = w_temp[ii];
*/
/*
		if(cs->dom->y[j] > 0.75){

		    cs->rho0[ii] = 1.05;
 		    cs->p0[ii]   = 1.0/cs->ig->gamma;
	 	    cs->U0[ii]   = -0.5;//sin(cs->dom->x[i]);
		    cs->V0[ii]   = 0.0;//sin(cs->dom->y[j]);
		    cs->W0[ii]   = 0.0;//sin(cs->dom->z[k]);
		}else{
		    cs->rho0[ii] = 0.95;
 		    cs->p0[ii]   = 1.0/cs->ig->gamma;
		    cs->U0[ii]   = 0.5;//sin(cs->dom->x[i]);
		    cs->V0[ii]   = 0.0;//sin(cs->dom->y[j]);
		    cs->W0[ii]   = 0.0;//sin(cs->dom->z[k]);
		}

		cs->V0[ii] += fRand(-0.3,0.3)*exp(-(cs->dom->y[j]-0.75)*(cs->dom->y[j]-0.75)*1000.0);
*/
	    }
	}
    }
/*
    delete[] u_temp;
    delete[] v_temp;
    delete[] w_temp;
*/
    cs->setInitialConditions();
    cs->dumpSolution();

    double initMach = 0.0, initMu = 0.0;;
    FOR_XYZ{
	initMach += cs->sos[ip]/((double)(Nx*Ny*Nz));
	initMu   += cs->mu[ip]/((double)(Nx*Ny*Nz));
    }
    initMach = u0/initMach;
    cout << " > Initial Mach # = " << initMach << endl;
    cout << " > Re = r*u0*L/mu = " << 1.0*u0*1.0/initMu;

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


	//Do some extra calculations if we need too...
	//cs->calcTurbulenceQuantities();
	cs->calcTaylorGreenQuantities();

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









