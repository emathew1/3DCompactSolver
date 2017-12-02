#ifndef _TIMESTEPPINGH_
#define _TIMESTEPPINGH_

#include <iostream>

class TimeStepping{

    public:

	enum TimeSteppingType {CONST_DT, CONST_CFL};
	TimeSteppingType timeSteppingType;
        double CFL;
	double dt;
	int maxTimeStep;
	double maxTime;
	int filterStep;  
	int checkStep;
	int dumpStep;

	TimeStepping(TimeSteppingType timeSteppingType, 
		     double CFL, int maxTimeStep, double maxTime, int filterStep, int checkStep, int dumpStep){
	    this->timeSteppingType = timeSteppingType;
	    this->CFL = CFL;
	    this->maxTimeStep = maxTimeStep;
	    this->maxTime = maxTime;
	    this->filterStep = filterStep;
	    this->checkStep  = checkStep;
	    this->dumpStep   = dumpStep;

	    std::cout << std::endl;
	    std::cout << " > Initializing time dependent options..." << std::endl;
	    if(timeSteppingType == CONST_CFL){
		std::cout << " > Using constant CFL timestepping of value " << CFL << std::endl;
	    }else if(timeSteppingType == CONST_DT){
		std::cout << " > Using constant dt timestepping of value " << dt << std::endl;
	    }
	    std::cout << " > Max time steps = " << maxTimeStep << ", max Time = " << maxTime << std::endl;
	    std::cout << " > Filtering every " << filterStep << " steps" << std::endl;

	}
};

#endif
