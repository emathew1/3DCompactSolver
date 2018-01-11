#ifndef _CSOLVERAWSH_
#define _CSOLVERAWSH_

#include <omp.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include "Macros.hpp"
#include "Utils.hpp"
#include "BC.hpp"
#include "TimeStepping.hpp"
#include "IdealGas.hpp"
#include "SpongeBC.hpp"
#include "Derivatives.hpp"
#include "Filter.hpp"

using namespace std::chrono;

class CSolver_AWS{

    public:

	Domain *dom;
	BC *bc;
	TimeStepping *ts;
	IdealGas *ig;
	double alphaF;
	double mu_ref;

	//Make local copies for macros...
	int Nx, Ny, Nz, N;
	int blocksize;

	//Track the local step in the Runge-Kutta integration...
	int rkStep;
	
	//Track the current time and timestep
        int timeStep;
        double time;
	int filterTimeStep;
	bool endFlag;

        std::chrono::system_clock::time_point t1Save, t2Save;


	//Kill solver condition
        bool done;

	//initial conditions
	double *rho0, 
	         *U0,
	         *V0,
		 *W0,
		 *p0;
	
	//non-conserved data
	double  *U, 
		*V,
		*W,
	        *T,
		*p,
	       *mu,
	      *Amu,
		*k, 
	      *sos;

	//derivatives of data
	double *Ux, *Uy, *Uz;
	double *Vx, *Vy, *Vz;
	double *Wx, *Wy, *Wz;

	double *Uxx, *Uyy, *Uzz;
	double *Vxx, *Vyy, *Vzz;
	double *Wxx, *Wyy, *Wzz;

	double *Uxy, *Uxz, *Uyz;
	double *Vxy, *Vxz, *Vyz;
	double *Wxy, *Wxz, *Wyz;

	double *Tx,  *Ty,  *Tz;
	double *Txx, *Tyy, *Tzz;

	double *contEulerX, *contEulerY, *contEulerZ;
	double *momXEulerX, *momXEulerY, *momXEulerZ;
	double *momYEulerX, *momYEulerY, *momYEulerZ;
	double *momZEulerX, *momZEulerY, *momZEulerZ;
	double *engyEulerX, *engyEulerY, *engyEulerZ;

	double *rho1,  *rhok,  *rhok2,  *rho2; 
	double *rhoU1, *rhoUk, *rhoUk2, *rhoU2; 
	double *rhoV1, *rhoVk, *rhoVk2, *rhoV2; 
	double *rhoW1, *rhoWk, *rhoWk2, *rhoW2; 
	double *rhoE1, *rhoEk, *rhoEk2, *rhoE2;

	double *temp, *temp2, *temp3, *temp4;
	double *transRho, *transRhoU, *transRhoV, *transRhoW, *transRhoE;
	double *transUx, *transVx, *transWx;
	double *transUy, *transVy, *transWy;

	//New memory allocated for AWS solver...
        double *transTempUy;
        double *transTempVy;
    	double *transTempWy;
    	double *transTempUyy;
    	double *transTempVyy;
    	double *transTempWyy;
    	double *transTempTy; 
    	double *transTempTyy; 
	double *transTempUxy;
	double *transTempVxy;
	double *transTempWxy;
	double *transTempUyz;
	double *transTempVyz;
	double *transTempWyz;
	double *transTempContEuler;
	double *transTempXEuler;
	double *transTempYEuler;
	double *transTempZEuler;
	double *transTempEngEuler;




	double *turbdiss, *uprime2, *uiprime2, *uvar, *kineticEng;

	bool spongeFlag;
	SpongeBC *spg; 

	//Moving Wall BC Velocities
	double X0WallV, X0WallW, X1WallV, X1WallW;
	double Y0WallU, Y0WallW, Y1WallU, Y1WallW;
	double Z0WallU, Z0WallV, Z1WallU, Z1WallV;

	enum Eqn {CONT, XMOM, YMOM, ZMOM, ENGY};

	Derivatives *derivX, *derivY, *derivZ;
	Filter *filtX, *filtY, *filtZ;


	//Constructor to use for this class...
	CSolver_AWS(Domain *dom, BC *bc, TimeStepping *ts, double alphaF, double mu_ref, int blocksize){

	    //Take in input information and initialize data structures...
	    this->dom = dom;
	    this->bc = bc;
	    this->ts = ts;
	    this->alphaF = alphaF;
	    this->mu_ref = mu_ref;

	    ig = new IdealGas(dom, mu_ref);

	    //give ourselves the local copies of the domain sizes
	    Nx = dom->Nx;
	    Ny = dom->Ny;
	    Nz = dom->Nz;
	    N  = dom->N;

	    //initialize time and timestep
   	    time = 0.0;
	    timeStep = 0;
	    filterTimeStep = 0;
	    endFlag = false;
            t1Save = std::chrono::system_clock::now();
            t2Save = std::chrono::system_clock::now();

	    done = false;

	    //Allocate our arrays for the solver data
	    initializeSolverData();		    	    

	    //Initialize the sponge boundary conditions if necessary
	    if(bc->bcX0 == BC::SPONGE || bc->bcX1 == BC::SPONGE || bc->bcY0 == BC::SPONGE || bc->bcY1 == BC::SPONGE || bc->bcZ0 == BC::SPONGE || bc->bcZ1 == BC::SPONGE){
		spongeFlag = true;
		spg = new SpongeBC(dom, ig, bc);
	    }else{
		spg = NULL;
	    }

	    //Initialize our derivative calculations for each direction...
	    derivX = new Derivatives(dom, bc->bcXType, Derivatives::DIRX);
	    derivY = new Derivatives(dom, bc->bcYType, Derivatives::DIRY);
	    derivZ = new Derivatives(dom, bc->bcZType, Derivatives::DIRZ);

	    //Initialize the filters we're going to use for each direction
	    filtX  = new Filter(alphaF, dom, bc->bcXType, Derivatives::DIRX);
	    filtY  = new Filter(alphaF, dom, bc->bcYType, Derivatives::DIRY);
	    filtZ  = new Filter(alphaF, dom, bc->bcZType, Derivatives::DIRZ);

	    this->blocksize = blocksize;

 	    X0WallV = 0.0; X0WallW = 0.0; X1WallV = 0.0; X1WallW = 0.0;
	    Y0WallU = 0.0; Y0WallW = 0.0; Y1WallU = 0.0; Y1WallW = 0.0;
	    Z0WallU = 0.0; Z0WallV = 0.0; Z1WallU = 0.0; Z1WallV = 0.0;

	}


	void initializeSolverData();

	void setInitialConditions();

	void calcDtFromCFL();

	inline double calcSpongeSource(double phi, double phiSpongeAvg, double sigma){
        	return sigma*(phiSpongeAvg - phi);
	};

	void preStepBCHandling();

	void preStepDerivatives();

	void solveContinuity();
	void solveXMomentum();
	void solveYMomentum();
	void solveZMomentum();
	void solveEnergy();

	void postStepBCHandling();

	void updateConservedData();
	
	void filterConservedData();

	void updateNonConservedData();

	void calcTurbulenceQuantities();
	void calcTaylorGreenQuantities();

	void updateSponge();

	void checkSolution();

	void dumpSolution();
	
	void checkEnd();

	void reportAll();

};

#endif
