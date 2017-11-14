#include "CSolver.hpp"

void CSolver::initializeSolverData(){

    cout << endl;
    cout << " > Allocating Solver Arrays..." << endl;
    double workSize = 0;
    workSize = 82.0 * (double)N * 8.0;
    cout << " > Need " << workSize/1024.0/1024.0/1024.0 << " Gb of memory required to allocate solver arrays " << endl;

    //3
    Ux = new double[N];
    Uy = new double[N];
    Uz = new double[N];

    //6
    Vx = new double[N];
    Vy = new double[N];
    Vz = new double[N];

    //9
    Wx = new double[N];
    Wy = new double[N];
    Wz = new double[N];

    //12
    Uxx = new double[N];
    Uyy = new double[N];
    Uzz = new double[N];

    //15
    Vxx = new double[N];
    Vyy = new double[N];
    Vzz = new double[N];

    //18
    Wxx = new double[N];
    Wyy = new double[N];
    Wzz = new double[N];

    //21
    Uxy = new double[N];
    Uxz = new double[N];
    Uyz = new double[N];

    //24
    Vxy = new double[N];
    Vxz = new double[N];
    Vyz = new double[N];

    //27
    Wxy = new double[N];
    Wxz = new double[N];
    Wyz = new double[N];

    //30
    Tx = new double[N];
    Ty = new double[N];
    Tz = new double[N];

    //33
    Txx = new double[N];
    Tyy = new double[N];
    Tzz = new double[N];

    //36
    contEulerX = new double[N];
    contEulerY = new double[N];
    contEulerZ = new double[N];

    //39
    momXEulerX = new double[N];
    momXEulerY = new double[N];
    momXEulerZ = new double[N];

    //42
    momYEulerX = new double[N];
    momYEulerY = new double[N];
    momYEulerZ = new double[N];

    //45
    momZEulerX = new double[N];
    momZEulerY = new double[N];
    momZEulerZ = new double[N];

    //48
    engyEulerX = new double[N];
    engyEulerY = new double[N];
    engyEulerZ = new double[N];

    //52
    rho1  = new double[N];
    rhok  = new double[N];
    rhok2 = new double[N];
    rho2  = new double[N];

    //56
    rhoU1  = new double[N];
    rhoUk  = new double[N];
    rhoUk2 = new double[N];
    rhoU2  = new double[N];

    //60
    rhoV1  = new double[N];
    rhoVk  = new double[N];
    rhoVk2 = new double[N];
    rhoV2  = new double[N];

    //64
    rhoW1  = new double[N];
    rhoWk  = new double[N];
    rhoWk2 = new double[N];
    rhoW2  = new double[N];
 
    //68
    rhoE1  = new double[N];
    rhoEk  = new double[N];
    rhoEk2 = new double[N];
    rhoE2  = new double[N];

    //73 these will be cleared though...
    rho0 = new double[N];
    U0   = new double[N];
    V0   = new double[N];
    W0   = new double[N];
    p0   = new double[N];

    //82 finally, the non-conserved data
    U   = new double[N];
    V   = new double[N];
    W   = new double[N];
    T   = new double[N];
    p   = new double[N];
    mu  = new double[N];
    Amu = new double[N];
    sos = new double[N];
    
}


void CSolver::setInitialConditions(){

    cout << endl;
    cout << " > Setting initial conditions..." << endl; 

    //just do the simple stuff in a loop...
    FOR_XYZ{
	rho1[ip] = rho0[ip];
	U[ip]	 = U0[ip];
	V[ip] 	 = V0[ip];
	W[ip] 	 = W0[ip];
	p[ip]	 = p0[ip];
	rhoU1[ip] = rho1[ip]*U[ip];	
	rhoV1[ip] = rho1[ip]*V[ip];	
	rhoW1[ip] = rho1[ip]*W[ip];
    }

    //Can now release initial condition data...
    delete[] rho0;
    delete[] U0;
    delete[] V0;
    delete[] W0;
    delete[] p0;

    //Call the ideal gas relations for the slightly more involved stuff..
    ig->solverhoE(rho1, p, U, V, W, rhoE1);
    ig->solveT(rho1, p, T);
    ig->solveMu(T, mu);
    ig->solveAmu(T, Amu);
    ig->solveSOS(rho1, p, sos);

    //This is where we'll do the boundary condition specific stuff...
    bool wallBCFlag = false;

    if(bc->bcX0 == BC::ADIABATIC_WALL){
    
	wallBCFlag = true;
    }

    if(bc->bcX1 == BC::ADIABATIC_WALL){
 
	wallBCFlag = true;
    }

    if(bc->bcY0 == BC::ADIABATIC_WALL){
 
	wallBCFlag = true;
    }

    if(bc->bcY1 == BC::ADIABATIC_WALL){
 
	wallBCFlag = true;
    }

    if(bc->bcZ0 == BC::ADIABATIC_WALL){
 
	wallBCFlag = true;
    }

    if(bc->bcZ1 == BC::ADIABATIC_WALL){
 
	wallBCFlag = true;
    }

    if(wallBCFlag == true){
	//Need to update the pressure, sos, and rhoE fields at the boundaries with walls...
	//Maybe need to change the way the ideal gas relation functions are programmed to be
	//"on-demand" for cells that are calculated rather than just doing all locations by default
    }

    if(spongeFlag == true){
	FOR_XYZ{
	    spg->spongeRhoAvg[ip]  = rho1[ip];
	    spg->spongeRhoUAvg[ip] = rhoU1[ip];
	    spg->spongeRhoVAvg[ip] = rhoV1[ip];
	    spg->spongeRhoWAvg[ip] = rhoW1[ip];
	    spg->spongeRhoEAvg[ip] = rhoE1[ip];
 	}
    }

    std::cout << " > Finished initialization of flow field " << std::endl;

}


void CSolver::calcDtFromCFL(){
    
    //Calculate the wave speed over the local spacings...
    double *UChar_dx = new double[N];
    FOR_XYZ{
	UChar_dx[ip] = (fabs(U[ip]) + sos[ip])/dom->dx + (fabs(V[ip])+sos[ip])/dom->dy + (fabs(W[ip]) + sos[ip])/dom->dz;
    }

    //Get the largest value in the domain
    double max_UChar_dx = -100000.0;
    FOR_XYZ{
	if(UChar_dx[ip] > max_UChar_dx){
	    max_UChar_dx = UChar_dx[ip];
	}
    }
    
    //done with UChar_dx
    delete[] UChar_dx;

    if(ts->CONST_CFL){
	ts->dt = ts->CFL/max_UChar_dx;
    }else if(ts->CONST_DT){
	ts->CFL = ts->dt*max_UChar_dx;
    }
  
    if(timeStep == 0){
	timeStep++;
	time = 0.0;
	cout << endl;
    }else{
	timeStep++;
	time += ts->dt;
    }

    cout << " > Timestep " << timeStep << ", T = " << time << ", dt = " << ts->dt << ", CFL = " << ts->CFL << endl; 

}

void CSolver::calcSpongeSource(double *phi, double *phiSpongeAvg, double *spongeSource){

    FOR_XYZ{
	spongeSource[ip] = spg->sigma[ip]*(phiSpongeAvg[ip] - phi[ip]);
    }

}

void CSolver::preStepBCHandling(double *rho, double *rhoU, double *rhoV, double *rhoW, double *rhoE){

    if(bc->bcX0 == BC::ADIABATIC_WALL){

    }

    if(bc->bcX1 == BC::ADIABATIC_WALL){

    }   

    if(bc->bcY0 == BC::ADIABATIC_WALL){

    }

    if(bc->bcY1 == BC::ADIABATIC_WALL){

    }

    if(bc->bcZ0 == BC::ADIABATIC_WALL){

    }

    if(bc->bcZ1 == BC::ADIABATIC_WALL){

    }

}


