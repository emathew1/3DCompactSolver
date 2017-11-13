#ifndef _TIMESTEPPINGH_
#define _TIMESTEPPINGH_

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

	    cout << endl;
	    cout << " >Initializing time dependent options..." << endl;
	    cout << " >Using " << endl;

	}
};

#endif
