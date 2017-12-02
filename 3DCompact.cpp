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
    double Lx = 9,//2.0*M_PI*((double)Nx - 1.0)/(double(Nx)), 
	   Ly = 9,//2.0*M_PI*((double)Ny - 1.0)/(double(Ny)), 
	   Lz = 9;//2.0*M_PI*((double)Nz - 1.0)/(double(Nz));
    Domain *dom = new Domain(Nx, Ny, Nz, Lx, Ly, Lz);



    ////////////////////////////////////
    //Time Stepping info intialization//
    ////////////////////////////////////
    TimeStepping::TimeSteppingType timeSteppingType = TimeStepping::CONST_CFL;
    double CFL 	     = 0.25;
    int maxTimeStep  = 10;
    double maxTime   = 10.0;
    int filterStep   = 5;
    int checkStep    = 1;
    int dumpStep     = 100;
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
    double alphaF = 0.49;
    double mu_ref = 0.0001;
    CSolver *cs   = new CSolver(dom, bc, ts, alphaF, mu_ref); 



    ///////////////////////////////
    //Set flow initial conditions//
    ///////////////////////////////
    FOR_Z{
	FOR_Y{
	    FOR_X{
		int ii = GET3DINDEX_XYZ;
		cs->rho0[ii] = 1.0;
		cs->U0[ii]   = (double)j*j;
		cs->V0[ii]   = (double)k*k;
		cs->W0[ii]   = (double)i*i;
		cs->p0[ii]   = 1.0/cs->ig->gamma;
	    }
	}
    }

    cs->setInitialConditions();
 
    cs->calcDtFromCFL();

    cs->rkStep = 1;

    auto t1 = std::chrono::system_clock::now();
    auto t2 = std::chrono::system_clock::now();
 
    double *test = new double[Nx*Ny*Nz];
    double *testTrans = new double[Nx*Ny*Nz];

    FOR_Z{
	FOR_Y{
	    FOR_X{
		int ii = GET3DINDEX_XYZ;
		test[ii] = (double)ii;
		testTrans[ii] = 0.0;
	    }
	}
    }
/*
    FOR_Z{
	FOR_Y{
	    FOR_X{
		int ii = GET3DINDEX_XYZ;
		cout << test[ii] << " ";
	    }
	    cout << endl;
	}
	cout << endl;
    }
*/
 
    t1 = std::chrono::system_clock::now();
    transposeXYZtoYZX(test, Nx, Ny, Nz, testTrans);
    t2 = std::chrono::system_clock::now();
    cout << "transposeXYZtoYXZ: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    transposeXYZtoYZX_Fast(test, Nx, Ny, Nz, testTrans, 1);
    t2 = std::chrono::system_clock::now();
    cout << "transposeXYZtoYXZ_Fast Block 1: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    transposeXYZtoYZX_Fast(test, Nx, Ny, Nz, testTrans, 2);
    t2 = std::chrono::system_clock::now();
    cout << "transposeXYZtoYXZ_Fast Block 2: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    transposeXYZtoYZX_Fast(test, Nx, Ny, Nz, testTrans, 4);
    t2 = std::chrono::system_clock::now();
    cout << "transposeXYZtoYXZ_Fast Block 4: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    transposeXYZtoYZX_Fast(test, Nx, Ny, Nz, testTrans, 8);
    t2 = std::chrono::system_clock::now();
    cout << "transposeXYZtoYXZ_Fast Block 8: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;
    t1 = std::chrono::system_clock::now();
    transposeXYZtoYZX_Fast(test, Nx, Ny, Nz, testTrans, 16);
    t2 = std::chrono::system_clock::now();
    cout << "transposeXYZtoYXZ_Fast Block 16: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    transposeXYZtoYZX_Fast(test, Nx, Ny, Nz, testTrans, 24);
    t2 = std::chrono::system_clock::now();
    cout << "transposeXYZtoYXZ_Fast Block 24: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    transposeXYZtoYZX_Fast(test, Nx, Ny, Nz, testTrans, 32);
    t2 = std::chrono::system_clock::now();
    cout << "transposeXYZtoYXZ_Fast Block 32: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;


    t1 = std::chrono::system_clock::now();
    cs->preStepDerivatives();
    t2 = std::chrono::system_clock::now();
    cout << "preStepDerivatives: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;


    t1 = std::chrono::system_clock::now();
    cs->solveContinuity();  
    t2 = std::chrono::system_clock::now();
    cout << "Continuity:" << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;


    t1 = std::chrono::system_clock::now();
    cs->solveXMomentum();
    t2 = std::chrono::system_clock::now();
    cout << "XMom:" << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    cs->solveYMomentum();
    t2 = std::chrono::system_clock::now();
    cout << "YMom:" << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    cs->solveZMomentum();
    t2 = std::chrono::system_clock::now();
    cout << "ZMom:" << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    cs->solveEnergy();
    t2 = std::chrono::system_clock::now();
    cout << "Energy:" << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    cs->updateConservedData();
    t2 = std::chrono::system_clock::now();
    cout << "updateConserved:" << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    cs->updateNonConservedData();
    t2 = std::chrono::system_clock::now();
    cout << "updateNonConserved:" << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    cs->filterConservedData();
    t2 = std::chrono::system_clock::now();
    cout << "filterConserved:" << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    cs->updateSponge();
    t2 = std::chrono::system_clock::now();
    cout << "updateSponge:" << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    cs->checkSolution();
    t2 = std::chrono::system_clock::now();
    cout << "checkSolution:" << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;




    return 0;
}









