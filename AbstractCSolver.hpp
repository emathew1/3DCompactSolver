#ifndef _CABSTRACTSOLVERH_
#define _CABSTRACTSOLVERH_

#include "Derivatives.hpp"
#include "Filter.hpp"
#include "IdealGas.hpp"
#include "TimeStepping.hpp"

class AbstractCSolver{

    public:
	
	int rkStep;
	bool rkLast;
	bool endFlag;
	bool useTiming;
	 
	//Initial Conditions  
	double *rho0, *p0, *U0, *V0, *W0;

	//Common Solver Data
        double *rho1,  *rhok,  *rhok2,  *rho2;
        double *rhoU1, *rhoUk, *rhoUk2, *rhoU2;
        double *rhoV1, *rhoVk, *rhoVk2, *rhoV2;
        double *rhoW1, *rhoWk, *rhoWk2, *rhoW2;
        double *rhoE1, *rhoEk, *rhoEk2, *rhoE2;       
 
	Domain *dom;
	BC *bc;
	TimeStepping *ts;
	IdealGas *ig;
        Derivatives *derivX, *derivY, *derivZ;
        Filter *filtX, *filtY, *filtZ;

	//Each Solver Class needs to have these functions to overwrite the pure virtual ones
	virtual void initializeSolverData() = 0;
	virtual void setInitialConditions() = 0;
	virtual void preStep() = 0;
	virtual void preSubStep() = 0;
	virtual void solveEqnSet() = 0;
	virtual void postSubStep() = 0;
	virtual void updateData() = 0;
	virtual void postStep() = 0;

};

#endif
