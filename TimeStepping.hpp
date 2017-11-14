#ifndef _TIMESTEPPINGH_
#define _TIMESTEPPINGH_

#include <iostream>

class TimeStepping{

    public:

	enum TimeSteppingType {CONST_DT, CONST_CFL};
	TimeSteppingType timeSteppingType;
        double CFL;
	int maxTimeStep;
	double maxTime;
	int filterStep;  

	TimeStepping(TimeSteppingType timeSteppingType, 
		     double CFL, int maxTimeStep, double maxTime, int filterStep){
	    this->timeSteppingType = timeSteppingType;
	    this->CFL = CFL;
	    this->maxTimeStep = maxTimeStep;
	    this->maxTime = maxTime;
	    this->filterStep = filterStep;

	    std::cout << std::endl;
	    std::cout << " >Initializing time dependent options..." << std::endl;

	}
};

#endif
