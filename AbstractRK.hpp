#ifndef _CABSTRACTRKH_
#define _CABSTRACTRKH_

#include "AbstractCSolver.hpp"

class AbstractRK{

    public:
	
	AbstractCSolver *cs;

	//Each RK Class needs to have these functions to overwrite the pure virtual ones
	virtual void executeSolverLoop() = 0;
	virtual void updateConservedData() = 0;

};

#endif
