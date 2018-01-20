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
    int    Nx = 256, 
	   Ny = 192, 
	   Nz = 128;
    double Lx = 172.0, 
	   Ly = 129.0, 
	   Lz = 86.0;;
    Domain *dom = new Domain(Nx, Ny, Nz, Lx, Ly, Lz);

    ////////////////////////////////////
    //Time Stepping info intialization//
    ////////////////////////////////////
    TimeStepping::TimeSteppingType timeSteppingType = TimeStepping::CONST_CFL;
    double CFL 	     = 0.75;
    int maxTimeStep  = 25000;
    double maxTime   = 3000.0;
    int filterStep   = 5;
    int checkStep    = 1;
    int dumpStep     = 500;
    TimeStepping *ts = new TimeStepping(timeSteppingType, CFL, maxTimeStep, maxTime, filterStep, checkStep, dumpStep);


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
    double alphaF  = 0.495;
    double mu_ref  = 0.00375;
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
    ifstream uFile, vFile, wFile, pFile;
    uFile.open("ShearLayer/U_uprime1_N128_k12_256x128x42.dat",ifstream::in);
    vFile.open("ShearLayer/V_uprime1_N128_k12_256x128x42.dat",ifstream::in);
    wFile.open("ShearLayer/W_uprime1_N128_k12_256x128x42.dat",ifstream::in);
    pFile.open("ShearLayer/P_uprime1_N128_k12_256x128x42.dat",ifstream::in);

    int inNx = 256;
    int inNy = 128;
    int inNz = 42;

    double *u_temp = new double[inNx*inNy*inNz];
    double *v_temp = new double[inNx*inNy*inNz];
    double *w_temp = new double[inNx*inNy*inNz];
    double *p_temp = new double[inNx*inNy*inNz];

    cout << " > Reading in perturbations...";
    for(int kp = 0; kp < inNz; kp++){
	for(int jp = 0; jp < inNy; jp++){
	    for(int ip = 0; ip < inNx; ip++){
		int ii =  kp*inNy*inNx + jp*inNx + ip;
		uFile >> u_temp[ii];
		vFile >> v_temp[ii];
		wFile >> w_temp[ii];
		pFile >> p_temp[ii];	
	    }
	}
    }

    uFile.close();
    vFile.close();
    wFile.close();
    pFile.close();
    cout << "done!" << endl;

    //Only thing is we need to switch Y and Z directions...
    //Divergence Free condition should hold whether we're mirroring
    //or rotating...

    //Now we need to scale the data for this simulation...
    //Data read in has been normalized so that u'=1
    double delta_u = 0.6;
    double intensity = 0.1*delta_u; //10% turbulence intensity?
    cout << " > Scaling the perturbations...";
    #pragma omp parallel for
    for(int ip = 0; ip < inNz*inNy*inNx; ip++){
	u_temp[ip] *= intensity;
	v_temp[ip] *= intensity;
	w_temp[ip] *= intensity;
	p_temp[ip] *= intensity*intensity; //scales as u^2
    }
    cout << "done!" << endl;;

    //Create density perturbuations assuming constant temperature...for now?
    double *r_temp = new double[inNx*inNy*inNz];
    #pragma omp parallel for
    for(int ip = 0; ip < inNz*inNy*inNx; ip++){
	r_temp[ip] = p_temp[ip];
    }

    //Lets plug these into a bigger fluc matrix to make it easier (lazy?)
    double *uFluc = new double[Nx*Ny*Nz];
    double *vFluc = new double[Nx*Ny*Nz];
    double *wFluc = new double[Nx*Ny*Nz];
    double *pFluc = new double[Nx*Ny*Nz];
    double *rFluc = new double[Nx*Ny*Nz];

    int midNum = Ny/2;
    int startYind = midNum - inNz/2;

    #pragma omp parallel for 
    FOR_Z{
	FOR_Y{
	    FOR_X{
                int ii = GET3DINDEX_XYZ;
	        uFluc[ii] = 0.0;
	        vFluc[ii] = 0.0;
	        wFluc[ii] = 0.0;
	        pFluc[ii] = 0.0;
	        rFluc[ii] = 0.0;
	    }
	}
    }

    #pragma omp parallel for 
    for(int kp = 0; kp < Nz; kp++){
	for(int jp = startYind; jp < startYind+inNz; jp++){
	    for(int ip = 0; ip < Nx; ip++){
		int ii  = kp*Ny*Nx + jp*Nx + ip;
		int iip = (jp-startYind)*Nz*Nx + kp*Nx + ip; 

	        double yTemp = dom->y[jp]-Ly/2.0;
		double scalingFactor = exp(-pow(yTemp/2.0,2.0));

		uFluc[ii] = u_temp[iip]*scalingFactor; 
		vFluc[ii] = v_temp[iip]*scalingFactor; 
		wFluc[ii] = w_temp[iip]*scalingFactor; 
		pFluc[ii] = p_temp[iip]*scalingFactor; 
		rFluc[ii] = r_temp[iip]*scalingFactor; 
	    }
	}
    }

    ///////////////////////////////
    //Set flow initial conditions//
    ///////////////////////////////
    #pragma omp parallel for
    FOR_Z{
	FOR_Y{
	    FOR_X{
		int ii = GET3DINDEX_XYZ;

		cs->rho0[ii] = 1.0 + rFluc[ii];
		cs->p0[ii]   = (1.0 + pFluc[ii])/cs->ig->gamma;
		cs->U0[ii]   = (delta_u/2.0)*tanh(-(dom->y[j]-Ly/2.0)/2.0) + uFluc[ii];
		cs->V0[ii]   = 0.0 + vFluc[ii];
		cs->W0[ii]   = 0.0 + wFluc[ii];
	    }
	}
    }
    delete[] u_temp;
    delete[] v_temp;
    delete[] w_temp;
    delete[] p_temp;
    delete[] r_temp;

    //Run the simulation!
    rk->executeSolverLoop();


    return 0;
}









