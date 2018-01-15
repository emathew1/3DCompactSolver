#include <omp.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

#include "Macros.hpp"
#include "Utils.hpp"
#include "BC.hpp"
#include "TimeStepping.hpp"
#include "Domain.hpp"
#include "AbstractCSolver.hpp"
#include "AbstractRK.hpp"
#include "TVDRK3.hpp"
#include "RK4.hpp"
#include "CSolver.hpp"
#include "CSolver_AWS.hpp"

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
          cout << " > Running with OpenMP using " << threadCount << " threads" << endl;
	  cout << " > Thread limit: " << omp_get_thread_limit() << endl;
      }
    }
    omp_set_max_active_levels(2);

    /////////////////////////
    //Initialize the Domain//
    /////////////////////////
    int    Nx = 128, 
	   Ny = 128, 
	   Nz = 128;
    double Lx = 2.0*M_PI, 
	   Ly = 2.0*M_PI, 
	   Lz = 2.0*M_PI;
    Domain *dom = new Domain(Nx, Ny, Nz, Lx, Ly, Lz);

    ////////////////////////////////////
    //Time Stepping info intialization//
    ////////////////////////////////////
    TimeStepping::TimeSteppingType timeSteppingType = TimeStepping::CONST_CFL;
    double CFL 	     = 1.0;
    int maxTimeStep  = 10000;
    double maxTime   = 10.0;
    int filterStep   = 1;
    int checkStep    = 1;
    int dumpStep     = 2500;
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
    double alphaF  = 0.495;
    double mu_ref  = 0.000704;
    int blocksize  = 16;
    bool useTiming = false;
    AbstractCSolver *cs;
    cs = new CSolver_AWS(dom, bc, ts, alphaF, mu_ref, blocksize, useTiming); 

    ///////////////////////////////////////////
    //Initialize Execution Loop and RK Method//
    ///////////////////////////////////////////
    AbstractRK *rk;
    rk = new TVDRK3(cs);

    /////////////////////////////
    //load in turbulence output//
    /////////////////////////////
    ifstream uFile, vFile, wFile;
    uFile.open("homogenous_turbulence/U_Mt0p3_N128_k8.dat",ifstream::in);
    vFile.open("homogenous_turbulence/V_Mt0p3_N128_k8.dat",ifstream::in);
    wFile.open("homogenous_turbulence/W_Mt0p3_N128_k8.dat",ifstream::in);

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

    ///////////////////////////////
    //Set flow initial conditions//
    ///////////////////////////////
    double u0 = 1.0;
    FOR_Z{
	FOR_Y{
	    FOR_X{
		int ii = GET3DINDEX_XYZ;

		cs->rho0[ii] = 1.0;
		cs->p0[ii]   = 4.0/cs->ig->gamma;
		cs->U0[ii]   = u_temp[ii];
		cs->V0[ii]   = v_temp[ii];
		cs->W0[ii]   = w_temp[ii];
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
    delete[] u_temp;
    delete[] v_temp;
    delete[] w_temp;

    //Run the simulation!
    rk->executeSolverLoop();


    return 0;
}









