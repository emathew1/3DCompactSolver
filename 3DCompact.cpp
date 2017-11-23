//#include <omp.h>
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
/*
    #pragma omp parallel
    {
      int threadCount = omp_get_num_threads();
      int threadID    = omp_get_thread_num();
      if(threadCount > 1 && threadID == 0){
          cout << "Running with OpenMP using " << threadCount << " threads" << endl;
      }
    }
*/

    /////////////////////////
    //Initialize the Domain//
    /////////////////////////
    int    Nx = 64, 
	   Ny = 64, 
	   Nz = 64;
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
    FOR_XYZ{
	cs->rho0[ip] = 1.0;
	cs->U0[ip]   = 0.1;
	cs->V0[ip]   = 0.1;
	cs->W0[ip]   = 0.1;
	cs->p0[ip]   = 1.0/cs->ig->gamma;
    }

    cs->setInitialConditions();
 
    cs->calcDtFromCFL();


    double *Test1 = new double[Nx*Ny*Nz];
    double *Test2 = new double[Nx*Ny*Nz];
    double *Test3 = new double[Nx*Ny*Nz];

    int i, j, k;
    FOR_XYZ{
	Test1[GET3DINDEX_XYZ] = 1.0;
	Test2[GET3DINDEX_XYZ] = 1.0;
	Test3[GET3DINDEX_XYZ] = 1.0;
    }


    auto t1 = std::chrono::system_clock::now();
    FOR_XYZ{
	Test1[GET3DINDEX_XYZ] = Test2[GET3DINDEX_XYZ]/Test3[GET3DINDEX_XYZ];
    }
    auto t2 = std::chrono::system_clock::now();
    cout << "SIMPLE MATH TEST: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    transposeXYZtoYZX(Test1, Nx, Ny, Nz, Test2);
    t2 = std::chrono::system_clock::now();
    cout << "Trans XYZ to YZX: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    transposeYZXtoZXY(Test2, Nx, Ny, Nz, Test3);
    t2 = std::chrono::system_clock::now();
    cout << "Trans YZX to ZXY: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    transposeZXYtoXYZ(Test3, Nx, Ny, Nz, Test1);
    t2 = std::chrono::system_clock::now();
    cout << "Trans ZXY to XYZ: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    transposeYZXtoXYZ(Test3, Nx, Ny, Nz, Test1);
    t2 = std::chrono::system_clock::now();
    cout << "Trans YZX to XYZ: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    cs->preStepDerivatives(0);
    t2 = std::chrono::system_clock::now();
    cout << "preStepDerivatives: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;

    t1 = std::chrono::system_clock::now();
    cs->preStepDerivatives2(0);
    t2 = std::chrono::system_clock::now();
    cout << "preStepDerivatives2: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;




    return 0;
}









