#ifndef _CSOLVERH_
#define _CSOLVERH_

#include <omp.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include "Macros.hpp"
#include "Utils.hpp"
#include "SpongeBC.hpp"
#include "AbstractCSolver.hpp"

using namespace std::chrono;

class CSolver: public AbstractCSolver{

    public:

	double alphaF;
	double mu_ref;

	//Make local copies for macros...
	int Nx, Ny, Nz, N;
	int blocksize;

	//Track the current time and timestep
        int timeStep;
        double time;
	int filterTimeStep;

        std::chrono::system_clock::time_point t1Save, t2Save;


	//Kill solver condition
        bool done;

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

	double *temp, *temp2, *temp3, *temp4;
	double *transRho, *transRhoU, *transRhoV, *transRhoW, *transRhoE;
	double *transUx, *transVx, *transWx;
	double *transUy, *transVy, *transWy;

	double *turbdiss, *uprime2, *uiprime2, *uvar, *kineticEng;

	bool spongeFlag;
	SpongeBC *spg; 

	//Moving Wall BC Velocities
	double X0WallV, X0WallW, X1WallV, X1WallW;
	double Y0WallU, Y0WallW, Y1WallU, Y1WallW;
	double Z0WallU, Z0WallV, Z1WallU, Z1WallV;

	enum Eqn {CONT, XMOM, YMOM, ZMOM, ENGY};


	//Constructor to use for this class...
	CSolver(Domain *dom, BC *bc, TimeStepping *ts, double alphaF, double mu_ref, int blocksize, bool useTiming){

	    //Take in input information and initialize data structures...
	    this->dom = dom;
	    this->bc = bc;
	    this->ts = ts;
	    this->alphaF = alphaF;
	    this->mu_ref = mu_ref;
	    this->blocksize = blocksize;
	    this->useTiming = useTiming;

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


 	    X0WallV = 0.0; X0WallW = 0.0; X1WallV = 0.0; X1WallW = 0.0;
	    Y0WallU = 0.0; Y0WallW = 0.0; Y1WallU = 0.0; Y1WallW = 0.0;
	    Z0WallU = 0.0; Z0WallV = 0.0; Z1WallU = 0.0; Z1WallV = 0.0;

	}

	//Functions required from AbstractCSolver
	void initializeSolverData();
	void setInitialConditions();
 	void preStep();
	void preSubStep();
	void solveEqnSet();
	void postSubStep();
	void updateData();
	void postStep();

	//Pre Step Functions
	void calcDtFromCFL();
	
	//Pre Sub Step Functions
	void preStepBCHandling();
	void preStepDerivatives();

	//Solve Eqn Set Functions
	void solveContinuity();
	void solveXMomentum();
	void solveYMomentum();
	void solveZMomentum();
	void solveEnergy();

	//Post Sub Step Functions
	void postStepBCHandling();
	void updateConservedData();
	void filterConservedData();
	void updateNonConservedData();

	//Post Step Functions
	void calcTurbulenceQuantities();
	void calcTaylorGreenQuantities();
	void updateSponge();
	void checkSolution();
	void dumpSolution();
	void checkEnd();
	void reportAll();

	//Inline, Or general solver functions
	inline double calcSpongeSource(double phi, double phiSpongeAvg, double sigma){
        	return sigma*(phiSpongeAvg - phi);
	};


};

#endif
