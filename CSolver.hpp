#ifndef _CSOLVERH_
#define _CSOLVERH_

#include "Macros.hpp"
#include "Utils.hpp"
#include "BC.hpp"
#include "TimeStepping.hpp"
#include "IdealGas.hpp"
#include "SpongeBC.hpp"
#include "Derivatives.hpp"

class CSolver{

    public:

	Domain *dom;
	BC *bc;
	TimeStepping *ts;
	IdealGas *ig;
	double alphaF;
	double mu_ref;

	//Make local copies for macros...
	int Nx, Ny, Nz, N;

	//Track the local step in the Runge-Kutta integration...
	int rkStep;
	
	//Track the current time and timestep
        int timeStep;
        double time;

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

	bool spongeFlag;
	SpongeBC *spg; 
	enum Eqn {CONT, XMOM, YMOM, ZMOM, ENGY};

	//Constructor to use for this class...
	CSolver(Domain *dom, BC *bc, TimeStepping *ts, double alphaF, double mu_ref){

	    //Take in input information and initialize data structures...
	    this->dom = dom;
	    this->bc = bc;
	    this->ts = ts;
	    this->alphaF = alphaF;
	    this->mu_ref = mu_ref;

	    ig = new IdealGas(dom);

	    //give ourselves the local copies of the domain sizes
	    Nx = dom->Nx;
	    Ny = dom->Ny;
	    Nz = dom->Nz;
	    N  = dom->N;

	    //initialize time and timestep
   	    time = 0.0;
	    timeStep = 0;

	    done = false;

	    //Allocate our arrays for the solver data
	    initializeSolverData();		    	    

	    //Initialize the sponge boundary conditions if necessary
	    if(bc->bcX0 == BC::SPONGE || bc->bcX1 == BC::SPONGE || bc->bcY0 == BC::SPONGE || bc->bcY1 == BC::SPONGE || bc->bcZ0 == BC::SPONGE || bc->bcZ1 == BC::SPONGE){
		spongeFlag = true;
		spg = new SpongeBC(dom, ig, bc);
	    }

	}


	void initializeSolverData();

	void setInitialConditions();

	void calcDtFromCFL();

	void calcSpongeSource(double *phi, double *phiSpongeAvg, double *spongeSource);

	void preStepBCHandling(double *rho, double *rhoU, double *rhoV, double *rhoW, double *rhoE);

	void preStepDerivatives(int rkStep);

	void solveContinuity();

	void solveXMomentumEuler();
	void solveXMomentumViscous();
	void solveXMomentum();

	void solveYMomentumEuler();
	void solveYMomentumViscous();
	void solveYMomentum();

	void solveZMomentumEuler();
	void solveZMomentumViscous();
	void solveZMomentum();

	void solveEnergyEuler();
	void solveEnergyViscous();
	void solveEnergy();

	void postStepBCHandling(double *rho, double *rhoU, double *rhoV, double *rhoW, double *rhoE);

	void updateConservedData(int rkStep);
	void updateNonConservedData(int rkStep);
	
	void filterConservedData();

	void updateSponge();

	void checkSolution();

};

#endif
