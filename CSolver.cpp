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
   
    temp  = new double[Nx*Ny*Nz];
    temp2 = new double[Nx*Ny*Nz]; 
    temp3 = new double[Nx*Ny*Nz]; 
    temp4 = new double[Nx*Ny*Nz]; 

    transRho = new double[Nx*Ny*Nz];
    transRhoU = new double[Nx*Ny*Nz];
    transRhoV = new double[Nx*Ny*Nz];
    transRhoW = new double[Nx*Ny*Nz];
    transRhoE = new double[Nx*Ny*Nz];
    transUx   = new double[Nx*Ny*Nz];
    transVx   = new double[Nx*Ny*Nz];
    transWx   = new double[Nx*Ny*Nz];
    transUy = new double[Nx*Ny*Nz];
    transVy = new double[Nx*Ny*Nz];
    transWy = new double[Nx*Ny*Nz];



 
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

double CSolver::calcSpongeSource(double phi, double phiSpongeAvg, double sigma){

    if(spongeFlag){
        return sigma*(phiSpongeAvg - phi);
    }else{
        return 0.0;
    }

}

void CSolver::preStepBCHandling(){

    double *rhoP, *rhoUP, *rhoVP, *rhoWP, *rhoEP;
    if(rkStep == 1){
	rhoP  = rho1;
	rhoUP = rhoU1;
	rhoVP = rhoV1;
	rhoWP = rhoW1;
	rhoEP = rhoE1;
    }else if(rkStep == 2 || rkStep == 3 || rkStep == 4){
	rhoP  = rhok;
	rhoUP = rhoUk; 
	rhoVP = rhoVk;
	rhoWP = rhoWk;
	rhoEP = rhoEk;
    }

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


void CSolver::preStepDerivatives(){

    //TODO
    //WHICH RHOU! etc. NEEDS TO CHANGE WITH THE RKSTEP!!!!
    double *rhoP;
    double *rhoUP;
    double *rhoVP;
    double *rhoWP;
    double *rhoEP;

    if(rkStep == 1){
	rhoP  = rho1;
	rhoUP = rhoU1;
	rhoVP = rhoV1;
	rhoWP = rhoW1;
	rhoEP = rhoE1; 
    }else if(rkStep == 2 || rkStep == 3 || rkStep == 4){
	rhoP  = rhok;
	rhoUP = rhoUk;
	rhoVP = rhoVk;
	rhoWP = rhoWk;
	rhoEP = rhoEk; 
    }

 
    ///////////////////
    // X-DERIVATIVES //
    ///////////////////

    //First we'll do all of the X-Direction derivatives since we're in XYZ order

    //Calculate the Euler Components of the equations... 
    #pragma omp parallel for
    FOR_XYZ temp[ip]  = rhoUP[ip]*U[ip] + p[ip];

    #pragma omp parallel for
    FOR_XYZ temp2[ip] = rhoVP[ip]*U[ip];

    #pragma omp parallel for
    FOR_XYZ temp3[ip] = rhoWP[ip]*U[ip];

    #pragma omp parallel for
    FOR_XYZ temp4[ip] = rhoEP[ip]*U[ip] + U[ip]*p[ip];
    
    auto t1 = std::chrono::system_clock::now();
    omp_set_nested(1);
    #pragma omp parallel sections num_threads(8) 
    {
        //Calculate the stuff needed for viscous derivatives
		#pragma omp section
        derivX->calc1stDerivField(U, Ux);
		#pragma omp section
        derivX->calc2ndDerivField(U, Uxx);
		#pragma omp section
        derivX->calc1stDerivField(V, Vx);
		#pragma omp section
        derivX->calc2ndDerivField(V, Vxx);
		#pragma omp section
        derivX->calc1stDerivField(W, Wx);
		#pragma omp section
        derivX->calc2ndDerivField(W, Wxx);
		#pragma omp section
        derivX->calc1stDerivField(T, Tx);
		#pragma omp section
        derivX->calc2ndDerivField(T, Txx);

        //Compute the Euler Derivatives
		#pragma omp section
        derivX->calc1stDerivField(rhoUP, contEulerX);
		#pragma omp section
        derivX->calc1stDerivField(temp,  momXEulerX);
		#pragma omp section
        derivX->calc1stDerivField(temp2, momYEulerX);
		#pragma omp section
        derivX->calc1stDerivField(temp3, momZEulerX);
		#pragma omp section
        derivX->calc1stDerivField(temp4, engyEulerX);

    }
    auto t2 = std::chrono::system_clock::now();
    cout << "part 1 deriv test, sectioned: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;    



    omp_set_nested(1);
    #pragma omp parallel sections num_threads(4) 
    {
        //Do the transposes
		#pragma omp section
        transposeXYZtoYZX(rhoP,  Nx, Ny, Nz, transRho);

		#pragma omp section
        transposeXYZtoYZX(rhoUP, Nx, Ny, Nz, transRhoU);

		#pragma omp section
        transposeXYZtoYZX(rhoVP, Nx, Ny, Nz, transRhoV);

		#pragma omp section
        transposeXYZtoYZX(rhoWP, Nx, Ny, Nz, transRhoW);

		#pragma omp section
        transposeXYZtoYZX(rhoEP, Nx, Ny, Nz, transRhoE);

		#pragma omp section
        transposeXYZtoYZX(Ux,    Nx, Ny, Nz, transUx);

		#pragma omp section
        transposeXYZtoYZX(Vx,    Nx, Ny, Nz, transVx);

		#pragma omp section
        transposeXYZtoYZX(Wx,    Nx, Ny, Nz, transWx);
    }

    ///////////////////
    // Y-DERIVATIVES //
    ///////////////////

    //Now recalculate properties in the new space
    #pragma omp parallel for
    FOR_XYZ U[ip] = transRhoU[ip]/transRho[ip];
    #pragma omp parallel for
    FOR_XYZ V[ip] = transRhoV[ip]/transRho[ip];
    #pragma omp parallel for
    FOR_XYZ W[ip] = transRhoW[ip]/transRho[ip];

    ig->solvep(transRho, transRhoE, U, V, W, p);
    ig->solveT(transRho, p, T);

    //Calculate the stuff for Euler derivatives in new space
    #pragma omp parallel for
    FOR_XYZ temp[ip]  = transRhoU[ip]*V[ip];
    #pragma omp parallel for
    FOR_XYZ temp2[ip] = transRhoV[ip]*V[ip] + p[ip];
    #pragma omp parallel for
    FOR_XYZ temp3[ip] = transRhoW[ip]*V[ip];
    #pragma omp parallel for
    FOR_XYZ temp4[ip] = transRhoE[ip]*V[ip] + V[ip]*p[ip];

    omp_set_nested(1);
    #pragma omp parallel sections num_threads(8)
    {

        //Calculate Viscous Derivatives
		#pragma omp section
        derivY->calc1stDerivField(U, Uy);

		#pragma omp section
        derivY->calc2ndDerivField(U, Uyy);

		#pragma omp section
        derivY->calc1stDerivField(V, Vy);

		#pragma omp section
        derivY->calc2ndDerivField(V, Vyy);

		#pragma omp section
        derivY->calc1stDerivField(W, Wy);

		#pragma omp section
        derivY->calc2ndDerivField(W, Wyy);

		#pragma omp section
        derivY->calc1stDerivField(T, Ty);

		#pragma omp section
        derivY->calc2ndDerivField(T, Tyy);

		#pragma omp section
        derivY->calc1stDerivField(transUx, Uxy);

		#pragma omp section
        derivY->calc1stDerivField(transVx, Vxy);

		#pragma omp section
        derivY->calc1stDerivField(transWx, Wxy);

        //Calculate the Euler Derivatives

		#pragma omp section
        derivY->calc1stDerivField(transRhoV, contEulerY);

		#pragma omp section
        derivY->calc1stDerivField(temp, 	 momXEulerY);

		#pragma omp section
        derivY->calc1stDerivField(temp2,	 momYEulerY);

		#pragma omp section
        derivY->calc1stDerivField(temp3,	 momZEulerY);

		#pragma omp section
        derivY->calc1stDerivField(temp4,	 engyEulerY);
    }

    //Moving data to ZXY

    omp_set_nested(1);
    #pragma omp parallel sections num_threads(4)
    {
        //Get the original conserved data from XYZ to ZXY
		#pragma omp section
        transposeXYZtoZXY(rhoP,  Nx, Ny, Nz, transRho);
		#pragma omp section
        transposeXYZtoZXY(rhoUP, Nx, Ny, Nz, transRhoU);
		#pragma omp section
        transposeXYZtoZXY(rhoVP, Nx, Ny, Nz, transRhoV);
		#pragma omp section
        transposeXYZtoZXY(rhoWP, Nx, Ny, Nz, transRhoW);
		#pragma omp section
        transposeXYZtoZXY(rhoEP, Nx, Ny, Nz, transRhoE);

        //Move the Y Derivative data from YZX to ZXY
		#pragma omp section
        transposeYZXtoZXY(Uy, Nx, Ny, Nz, transUy);
		#pragma omp section
        transposeYZXtoZXY(Vy, Nx, Ny, Nz, transVy);
		#pragma omp section
        transposeYZXtoZXY(Wy, Nx, Ny, Nz, transWy);

        //Move the X Derivative data from XYZ to ZXY
		#pragma omp section
        transposeXYZtoZXY(Ux, Nx, Ny, Nz, transUx);
		#pragma omp section
        transposeXYZtoZXY(Vx, Nx, Ny, Nz, transVx);
		#pragma omp section
        transposeXYZtoZXY(Wx, Nx, Ny, Nz, transWx);
    }

    //Moving Data from YZX to XYZ

    #pragma omp parallel sections num_threads(4)
    {

		#pragma omp section
		{
            memcpy(temp, Uy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(temp, Nx, Ny, Nz, Uy);
		}

		#pragma omp section
		{
            memcpy(temp2, Uyy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(temp2, Nx, Ny, Nz, Uyy);
		}

		#pragma omp section
		{
            memcpy(temp3, Vy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(temp3, Nx, Ny, Nz, Vy);
		}

		#pragma omp section
		{
            memcpy(temp4, Vyy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(temp4, Nx, Ny, Nz, Vyy);
		}
	
		#pragma omp section
		{
            memcpy(U, Wy, sizeof(double)*Nx*Ny*Nz); //Starting here getting cute with reusing memory
            transposeYZXtoXYZ(U, Nx, Ny, Nz, Wy);
		}


		#pragma omp section
		{
            memcpy(V, Wyy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(V, Nx, Ny, Nz, Wyy);
		}


		#pragma omp section
		{
            memcpy(W, Ty, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(W, Nx, Ny, Nz, Ty);
		}

		#pragma omp section
		{
            memcpy(T, Tyy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(T, Nx, Ny, Nz, Tyy);
		}
    } 

    #pragma omp parallel sections num_threads(4)
    {

		#pragma omp section
		{
            memcpy(temp, Uxy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(temp, Nx, Ny, Nz, Uxy);
		}

		#pragma omp section
		{
            memcpy(temp2, Vxy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(temp2, Nx, Ny, Nz, Vxy);
		}

		#pragma omp section
		{
            memcpy(temp3, Wxy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(temp3, Nx, Ny, Nz, Wxy);
		}

		#pragma omp section
		{
            memcpy(temp4, contEulerY, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(temp4, Nx, Ny, Nz, contEulerY);
		}

		#pragma omp section
		{
            memcpy(U, momXEulerY, sizeof(double)*Nx*Ny*Nz); //Getting cute with memory here again
            transposeYZXtoXYZ(U, Nx, Ny, Nz, momXEulerY);
		}

		#pragma omp section
		{
            memcpy(V, momYEulerY, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(V, Nx, Ny, Nz, momYEulerY);
		}

		#pragma omp section
		{
            memcpy(W, momZEulerY, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(W, Nx, Ny, Nz, momZEulerY);
		}

		#pragma omp section
		{
            memcpy(T, engyEulerY, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ(T, Nx, Ny, Nz, engyEulerY);
		}
    }

    ///////////////////
    // Z-DERIVATIVES //
    ///////////////////

    //Now recalculate properties in the new space
    #pragma omp parallel for
    FOR_XYZ U[ip] = transRhoU[ip]/transRho[ip];
    #pragma omp parallel for
    FOR_XYZ V[ip] = transRhoV[ip]/transRho[ip];
    #pragma omp parallel for
    FOR_XYZ W[ip] = transRhoW[ip]/transRho[ip];

    ig->solvep(transRho, transRhoE, U, V, W, p);
    ig->solveT(transRho, p, T);

    //Calculate the stuff for the Euler Derivatives
    #pragma omp parallel for
    FOR_XYZ temp[ip]  = transRhoU[ip]*W[ip];
    #pragma omp parallel for
    FOR_XYZ temp2[ip] = transRhoV[ip]*W[ip] + p[ip];
    #pragma omp parallel for
    FOR_XYZ temp3[ip] = transRhoW[ip]*W[ip] + p[ip];
    #pragma omp parallel for
    FOR_XYZ temp4[ip] = transRhoE[ip]*W[ip] + W[ip]*p[ip];

    omp_set_nested(1);
    #pragma omp parallel sections num_threads(8)
    {
        //Calculate the viscous derivatives
		#pragma omp section
        derivZ->calc1stDerivField(U, Uz);
		#pragma omp section
        derivZ->calc2ndDerivField(U, Uzz);

		#pragma omp section
        derivZ->calc1stDerivField(V, Vz);
		#pragma omp section
        derivZ->calc2ndDerivField(V, Vzz);

		#pragma omp section
        derivZ->calc1stDerivField(W, Wz);
		#pragma omp section
        derivZ->calc2ndDerivField(W, Wzz);

		#pragma omp section
        derivZ->calc1stDerivField(T, Tz);
		#pragma omp section
        derivZ->calc2ndDerivField(T, Tzz);

		#pragma omp section
        derivZ->calc1stDerivField(transUx, Uxz);
		#pragma omp section
        derivZ->calc1stDerivField(transVx, Vxz);
		#pragma omp section
        derivZ->calc1stDerivField(transWx, Wxz);

		#pragma omp section
        derivZ->calc1stDerivField(transUy, Uyz);
		#pragma omp section
        derivZ->calc1stDerivField(transVy, Vyz);
		#pragma omp section
        derivZ->calc1stDerivField(transWy, Wyz);

        //Calculate the Euler Derivatives
		#pragma omp section
        derivZ->calc1stDerivField(transRhoW, contEulerZ);
		#pragma omp section
        derivZ->calc1stDerivField(temp,	     momXEulerZ);
		#pragma omp section
        derivZ->calc1stDerivField(temp2,     momYEulerZ);
		#pragma omp section
        derivZ->calc1stDerivField(temp3,     momZEulerZ);
		#pragma omp section
        derivZ->calc1stDerivField(temp4,     engyEulerZ);
    }


    omp_set_nested(1);
    #pragma omp parallel sections num_threads(4)
    {

		//Moving all the data back to XYZ
		#pragma omp section
		{
			memcpy(temp, Uz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(temp, Nx, Ny, Nz, Uz);
		}

		#pragma omp section
		{
			memcpy(temp2, Uzz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(temp2, Nx, Ny, Nz, Uzz);
		}

		#pragma omp section
		{
			memcpy(temp3, Vz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(temp3, Nx, Ny, Nz, Vz);
		}

		#pragma omp section
		{
			memcpy(temp4, Vzz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(temp4, Nx, Ny, Nz, Vzz);
		}

		#pragma omp section
		{
			memcpy(transUx, Wz, sizeof(double)*Nx*Ny*Nz); //Once again getting cute with memory starting here
			transposeZXYtoXYZ(transUx, Nx, Ny, Nz, Wz);
		}

		#pragma omp section
		{
			memcpy(transVx, Wzz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(transVx, Nx, Ny, Nz, Wzz);
		}

		#pragma omp section
		{
			memcpy(transWx, Tz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(transWx, Nx, Ny, Nz, Tz);
		}

		#pragma omp section
		{
			memcpy(transUy, Tzz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(transWy, Nx, Ny, Nz, Tzz);
		}

		#pragma omp section
		{
			memcpy(transVy, Uxz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(transVy, Nx, Ny, Nz, Uxz);
		}

		#pragma omp section
		{
			memcpy(transWy, Vxz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(transWy, Nx, Ny, Nz, Vxz);
		}

		#pragma omp section
		{
			memcpy(transRhoU, Wxz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(transRhoU, Nx, Ny, Nz, Wxz);
		}

		#pragma omp section
		{
			memcpy(transRhoV, Uyz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(transRhoV, Nx, Ny, Nz, Uyz);
		}

		#pragma omp section
		{
			memcpy(transRhoW, Vyz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(transRhoW, Nx, Ny, Nz, Vyz);
		}

		#pragma omp section
		{
			memcpy(transRhoE, Wyz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(transRhoE, Nx, Ny, Nz, Wyz);
		}

		#pragma omp section
		{
			memcpy(U, contEulerZ, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(U, Nx, Ny, Nz, contEulerZ);
		}

		#pragma omp section
		{
			memcpy(V, momXEulerZ, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(V, Nx, Ny, Nz, momXEulerZ);
		}

		#pragma omp section
		{
			memcpy(W, momYEulerZ, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(W, Nx, Ny, Nz, momYEulerZ);
		}

		#pragma omp section
		{
			memcpy(T, momZEulerZ, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(T, Nx, Ny, Nz, momZEulerZ);
		}
	
		#pragma omp section
		{
			memcpy(p, engyEulerZ, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ(p, Nx, Ny, Nz, engyEulerZ);
		}
    }

    //Going back to original...
    #pragma omp parallel for
    FOR_XYZ U[ip] = rhoUP[ip]/rhoP[ip];
    #pragma omp parallel for
    FOR_XYZ V[ip] = rhoVP[ip]/rhoP[ip];
    #pragma omp parallel for
    FOR_XYZ W[ip] = rhoWP[ip]/rhoP[ip];

    ig->solvep(rhoP, rhoEP, U, V, W, p);
    ig->solveT(rhoP, p, T);

}


void CSolver::solveContinuity(){
    
	#pragma omp parallel
	{

	#pragma omp for
    FOR_XYZ rhok2[ip]  = - contEulerX[ip] - contEulerY[ip] - contEulerZ[ip];

    if(spongeFlag){
        double *rhoP;
        if(rkStep == 1){
	    rhoP = rho1;
        }else if(rkStep == 2 || rkStep == 3 || rkStep == 4){
	    rhoP = rhok;
        }

	    #pragma omp for
	    FOR_XYZ{
	        rhok2[ip]  += calcSpongeSource(rhoP[ip], spg->spongeRhoAvg[ip], spg->sigma[ip]);	
	    }
    }

	#pragma omp for
    FOR_XYZ rhok2[ip] *= ts->dt;

	}
}

void CSolver::solveXMomentum(){

    #pragma omp parallel
    {

    	double MuX, MuY, MuZ;

		#pragma omp for
    	FOR_XYZ{
			//Viscous Terms
        	rhoUk2[ip]  = mu[ip]*((4.0/3.0)*Uxx[ip] + Uyy[ip] + Uzz[ip] + (1.0/3.0)*Vxy[ip] + (1.0/3.0)*Wxz[ip]);
    	}

		#pragma omp for
    	FOR_XYZ{
			MuX = Amu[ip]*Tx[ip];
			rhoUk2[ip] += (4.0/3.0)*MuX*(Ux[ip] - 0.5*Vy[ip] - 0.5*Wz[ip]);
    	}

		#pragma omp for
    	FOR_XYZ{
			MuY = Amu[ip]*Ty[ip];
			rhoUk2[ip] += MuY*(Uy[ip] + Vx[ip]); 
    	}
  
		#pragma omp for
    	FOR_XYZ{
			MuZ = Amu[ip]*Tz[ip];
			rhoUk2[ip] += MuZ*(Wx[ip] + Uz[ip]); 
    	}

		#pragma omp for
    	FOR_XYZ{
			//Euler Terms
			rhoUk2[ip] += -momXEulerX[ip] -momXEulerY[ip] -momXEulerZ[ip];
    	}

    	if(spongeFlag){
        	double *rhoUP;
        	if(rkStep == 1){
            	rhoUP = rhoU1;
        	}else if(rkStep == 2 || rkStep == 3 || rkStep == 4){
            	rhoUP = rhoUk;
        	}

			#pragma omp for
        	FOR_XYZ{
            	rhoUk2[ip]  += calcSpongeSource(rhoUP[ip], spg->spongeRhoUAvg[ip], spg->sigma[ip]);
        	}
   		}

		#pragma omp for
	    FOR_XYZ rhoUk2[ip] *= ts->dt;

    }
}

void CSolver::solveYMomentum(){

    double MuX, MuY, MuZ;
    FOR_XYZ{
        rhoVk2[ip]  = mu[ip]*((4.0/3.0)*Vyy[ip] + Vxx[ip] + Vzz[ip] + (1.0/3.0)*Uxy[ip] + (1.0/3.0)*Wyz[ip]);
    }
	//Viscous Terms
    FOR_XYZ{ 
	MuX = Amu[ip]*Tx[ip];
	rhoVk2[ip] += (4.0/3.0)*MuY*(Vy[ip] - 0.5*Ux[ip] - 0.5*Wz[ip]);
    }

    FOR_XYZ{
	MuY = Amu[ip]*Ty[ip];
	rhoVk2[ip] += MuX*(Uy[ip] + Vx[ip]); 
    }

    FOR_XYZ{
  	MuZ = Amu[ip]*Tz[ip];
	rhoVk2[ip] += MuZ*(Wy[ip] + Vz[ip]); 
    }

    //Euler Terms
    FOR_XYZ{
	rhoVk2[ip] += -momYEulerX[ip] -momYEulerY[ip] -momYEulerZ[ip];
	
    }

    if(spongeFlag){
        double *rhoVP;
        if(rkStep == 1){
            rhoVP = rhoV1;
        }else if(rkStep == 2 || rkStep == 3 || rkStep == 4){
            rhoVP = rhoVk;
        }

        FOR_XYZ{
            rhoVk2[ip]  += calcSpongeSource(rhoVP[ip], spg->spongeRhoVAvg[ip], spg->sigma[ip]);
        }
    }

    FOR_XYZ rhoVk2[ip] *= ts->dt;

}

void CSolver::solveZMomentum(){

    double MuX, MuY, MuZ;
    FOR_XYZ{
        rhoWk2[ip]  = mu[ip]*((4.0/3.0)*Wzz[ip] + Wyy[ip] + Wxx[ip] + (1.0/3.0)*Uxz[ip] + (1.0/3.0)*Vyz[ip]);
    }

	//Viscous Terms
    FOR_XYZ{
	MuX = Amu[ip]*Tx[ip];
	rhoWk2[ip] += (4.0/3.0)*MuZ*(Wz[ip] - 0.5*Ux[ip] - 0.5*Vy[ip]);
    }

    FOR_XYZ{
	MuY = Amu[ip]*Ty[ip];
	rhoWk2[ip] += MuX*(Wx[ip] + Uz[ip]); 
    }

    FOR_XYZ{
	MuZ = Amu[ip]*Tz[ip];
	rhoWk2[ip] += MuY*(Wy[ip] + Vz[ip]); 
    }

    FOR_XYZ{
	//Euler Terms
	rhoWk2[ip] += -momZEulerX[ip] -momZEulerY[ip] -momZEulerZ[ip];
	
    }

    if(spongeFlag){
        double *rhoWP;
        if(rkStep == 1){
            rhoWP = rhoW1;
        }else if(rkStep == 2 || rkStep == 3 || rkStep == 4){
            rhoWP = rhoWk;
        }

        FOR_XYZ{
            rhoWk2[ip]  += calcSpongeSource(rhoWP[ip], spg->spongeRhoWAvg[ip], spg->sigma[ip]);
        }
    }

    FOR_XYZ rhoWk2[ip] *= ts->dt;

}


void CSolver::solveEnergy(){

    double MuX, MuY, MuZ;
    double *qtemp = new double[N];
    double *vtemp1 = new double[N];
    double *vtemp2 = new double[N];
    double *engyEuler = new double[N];

    //Heat Transfer Terms
    FOR_XYZ{
	MuX = Amu[ip]*Tx[ip];
	qtemp[ip]  = MuX*Tx[ip];
    }

    FOR_XYZ{
	MuY = Amu[ip]*Ty[ip];
	qtemp[ip] += MuY*Ty[ip];
    }  

    FOR_XYZ{
	MuZ = Amu[ip]*Tz[ip];
	qtemp[ip] += MuZ*Tz[ip];
    }

    FOR_XYZ{ 
	//Heat Transfer Terms
	qtemp[ip] += mu[ip]*Txx[ip] + mu[ip]*Tyy[ip] + mu[ip]*Tzz[ip];
	qtemp[ip] *= ig->cp/ig->Pr;
    }

    //Viscous Energy terms w/o viscosity derivatives...
    FOR_XYZ vtemp1[ip]  = U[ip]*((4.0/3.0)*Uxx[ip] + Uyy[ip] + Uzz[ip]);
    FOR_XYZ vtemp1[ip] += V[ip]*(Vxx[ip] + (4.0/3.0)*Vyy[ip] + Vzz[ip]);
    FOR_XYZ vtemp1[ip] += W[ip]*(Wxx[ip] + Wyy[ip] + (4.0/3.0)*Wzz[ip]);
	
    FOR_XYZ vtemp1[ip] += (4.0/3.0)*(Ux[ip]*Ux[ip] + Vy[ip]*Vy[ip] + Wz[ip]*Wz[ip]);

    FOR_XYZ vtemp1[ip] += Uy[ip]*Uy[ip] + Uz[ip]*Uz[ip];
    FOR_XYZ vtemp1[ip] += Vx[ip]*Vx[ip] + Vz[ip]*Vz[ip];
    FOR_XYZ vtemp1[ip] += Wx[ip]*Wx[ip] + Wy[ip]*Wy[ip];

    FOR_XYZ vtemp1[ip] += -(4.0/3.0)*(Ux[ip]*Vy[ip] + Ux[ip]*Wz[ip] + Vy[ip]*Wz[ip]);
	
    FOR_XYZ vtemp1[ip] += 2.0*(Uy[ip]*Vx[ip] + Uz[ip]*Wx[ip] + Vz[ip]*Wy[ip]);

    FOR_XYZ vtemp1[ip] += (1.0/3.0)*(U[ip]*Vxy[ip] + U[ip]*Wxz[ip] + V[ip]*Uxy[ip]); 
    FOR_XYZ vtemp1[ip] += (1.0/3.0)*(V[ip]*Wyz[ip] + W[ip]*Uxz[ip] + W[ip]*Vyz[ip]); 
	
    FOR_XYZ vtemp1[ip] *= mu[ip];

    //Viscous Energy terms w/ viscosity derivatives...
    FOR_XYZ{ 
	MuX = Amu[ip]*Tx[ip];
	MuY = Amu[ip]*Ty[ip];
	MuZ = Amu[ip]*Tz[ip];
	vtemp2[ip]   = (4.0/3.0)*(U[ip]*MuX*Ux[ip] + V[ip]*MuY*Vy[ip] + W[ip]*MuZ*Wz[ip]);
    }

    FOR_XYZ{
	MuX = Amu[ip]*Tx[ip];
	vtemp2[ip]  += -(2.0/3.0)*U[ip]*MuX*(Vy[ip] + Wz[ip]);
    }

    FOR_XYZ{
	MuY = Amu[ip]*Ty[ip];
	vtemp2[ip]  += -(2.0/3.0)*V[ip]*MuY*(Ux[ip] + Wz[ip]);
    }

    FOR_XYZ{
	MuZ = Amu[ip]*Tz[ip];
	vtemp2[ip]  += -(2.0/3.0)*W[ip]*MuZ*(Ux[ip] + Vy[ip]);
    }

    FOR_XYZ{
	MuY = Amu[ip]*Ty[ip];
	vtemp2[ip]  += U[ip]*MuY*(Uy[ip] + Vx[ip]);
    }

    FOR_XYZ{
	MuZ = Amu[ip]*Tz[ip];
        vtemp2[ip]  += U[ip]*MuZ*(Uz[ip] + Wx[ip]);
    }

    FOR_XYZ{
	MuX = Amu[ip]*Tx[ip];
	vtemp2[ip]  += V[ip]*MuX*(Uy[ip] + Vx[ip]);
    }

    FOR_XYZ{
	MuZ = Amu[ip]*Tz[ip];
	vtemp2[ip]  += V[ip]*MuZ*(Vz[ip] + Wy[ip]);
    }

    FOR_XYZ{
	MuX = Amu[ip]*Tx[ip];
	vtemp2[ip]  += W[ip]*MuX*(Uz[ip] + Wx[ip]);
    }

    FOR_XYZ{
	MuY = Amu[ip]*Ty[ip];
	vtemp2[ip]  += W[ip]*MuY*(Vz[ip] + Wy[ip]);
    }

    //Euler terms
    FOR_XYZ engyEuler[ip]  = -engyEulerX[ip] - engyEulerY[ip] - engyEulerZ[ip];

    FOR_XYZ{
	double engySponge;
        if(spongeFlag){
            double *rhoEP;
            if(rkStep == 1){
                rhoEP = rhoE1;
            }else if(rkStep == 2 || rkStep == 3 || rkStep == 4){
                rhoEP = rhoEk;
            }

            engySponge  = calcSpongeSource(rhoEP[ip], spg->spongeRhoEAvg[ip], spg->sigma[ip]);
        }else{
	    engySponge = 0.0;
	}	

	rhoEk2[ip] = ts->dt*(qtemp[ip] + vtemp1[ip] + vtemp2[ip] + engyEuler[ip] + engySponge);

    }

    delete[] qtemp;
    delete[] vtemp1;
    delete[] vtemp2;
    delete[] engyEuler;

}

void CSolver::postStepBCHandling(){

    double *rhoP, *rhoUP, *rhoVP, *rhoWP, *rhoEP;
    if(rkStep == 1){
	rhoP  = rho1;
	rhoUP = rhoU1;
	rhoVP = rhoV1;
	rhoWP = rhoW1;
	rhoEP = rhoE1;
    }else if(rkStep == 2 || rkStep == 3 || rkStep == 4){
	rhoP  = rhok;
	rhoUP = rhoUk; 
	rhoVP = rhoVk;
	rhoWP = rhoWk;
	rhoEP = rhoEk;
    }

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

void CSolver::updateConservedData(){

    if(rkStep == 1){

	//Add to final solution
	FOR_XYZ rho2[ip]  = rho1[ip]  + rhok2[ip]/6.0;
	FOR_XYZ rhoU2[ip] = rhoU1[ip] + rhoUk2[ip]/6.0;
	FOR_XYZ rhoV2[ip] = rhoV1[ip] + rhoVk2[ip]/6.0;
	FOR_XYZ rhoW2[ip] = rhoW1[ip] + rhoWk2[ip]/6.0;
	FOR_XYZ rhoE2[ip] = rhoE1[ip] + rhoEk2[ip]/6.0;

	//Calculate intermediate solution
	FOR_XYZ rhok[ip]  = rho1[ip]  + rhok2[ip]/2.0; 
	FOR_XYZ rhoUk[ip] = rhoU1[ip] + rhoUk2[ip]/2.0; 
	FOR_XYZ rhoVk[ip] = rhoV1[ip] + rhoVk2[ip]/2.0; 
	FOR_XYZ rhoWk[ip] = rhoW1[ip] + rhoWk2[ip]/2.0; 
	FOR_XYZ rhoEk[ip] = rhoE1[ip] + rhoEk2[ip]/2.0; 

    }else if(rkStep == 2){

	//Add to final solution
	FOR_XYZ rho2[ip]  += rhok2[ip]/3.0;
	FOR_XYZ rhoU2[ip] += rhoUk2[ip]/3.0;
	FOR_XYZ rhoV2[ip] += rhoVk2[ip]/3.0;
	FOR_XYZ rhoW2[ip] += rhoWk2[ip]/3.0;
	FOR_XYZ rhoE2[ip] += rhoEk2[ip]/3.0;

	//Calculate intermediate solution
	FOR_XYZ rhok[ip]  = rho1[ip]  + rhok2[ip]/2.0; 
	FOR_XYZ rhoUk[ip] = rhoU1[ip] + rhoUk2[ip]/2.0; 
	FOR_XYZ rhoVk[ip] = rhoV1[ip] + rhoVk2[ip]/2.0; 
	FOR_XYZ rhoWk[ip] = rhoW1[ip] + rhoWk2[ip]/2.0; 
	FOR_XYZ rhoEk[ip] = rhoE1[ip] + rhoEk2[ip]/2.0; 

    }else if(rkStep == 3){

	//Add to final solution
	FOR_XYZ rho2[ip]  += rhok2[ip]/3.0;
	FOR_XYZ rhoU2[ip] += rhoUk2[ip]/3.0;
	FOR_XYZ rhoV2[ip] += rhoVk2[ip]/3.0;
	FOR_XYZ rhoW2[ip] += rhoWk2[ip]/3.0;
	FOR_XYZ rhoE2[ip] += rhoEk2[ip]/3.0;

	//Calculate intermediate solution
	FOR_XYZ rhok[ip]  = rho1[ip]  + rhok2[ip]; 
	FOR_XYZ rhoUk[ip] = rhoU1[ip] + rhoUk2[ip]; 
	FOR_XYZ rhoVk[ip] = rhoV1[ip] + rhoVk2[ip]; 
	FOR_XYZ rhoWk[ip] = rhoW1[ip] + rhoWk2[ip]; 
	FOR_XYZ rhoEk[ip] = rhoE1[ip] + rhoEk2[ip]; 

    }else if(rkStep == 4){

	//Add to final solution
	FOR_XYZ rho2[ip]  += rhok2[ip]/6.0;
	FOR_XYZ rhoU2[ip] += rhoUk2[ip]/6.0;
	FOR_XYZ rhoV2[ip] += rhoVk2[ip]/6.0;
	FOR_XYZ rhoW2[ip] += rhoWk2[ip]/6.0;
	FOR_XYZ rhoE2[ip] += rhoEk2[ip]/6.0;

    }

}

void CSolver::filterConservedData(){

    //Need to do round robin of filtering directions...
    if(ts->filterStep%timeStep == 0){

        //Advance the filtering time step       
        filterTimeStep++;

        //Going to try and be cute to minimize dmemory allocation
        if(filterTimeStep%3 == 1){

            //Here we'll do X->Y->Z     
            filtX->filterField(rho2,  rho1);
            filtX->filterField(rhoU2, temp);
            filtX->filterField(rhoV2, temp2);
            filtX->filterField(rhoW2, temp3);
            filtX->filterField(rhoE2, temp4);

            //Do the transpose to YZX space
            transposeXYZtoYZX(rho1,   Nx, Ny, Nz, rho2);
            transposeXYZtoYZX(temp,   Nx, Ny, Nz, rhoU2);
            transposeXYZtoYZX(temp2,  Nx, Ny, Nz, rhoV2);
            transposeXYZtoYZX(temp3,  Nx, Ny, Nz, rhoW2);
            transposeXYZtoYZX(temp4,  Nx, Ny, Nz, rhoE2);

            //filter in the Y direction
            filtY->filterField(rho2,  rho1);
            filtY->filterField(rhoU2, temp);
            filtY->filterField(rhoV2, temp2);
            filtY->filterField(rhoW2, temp3);
            filtY->filterField(rhoE2, temp4);

            //tranpose from YZX to ZXY
            transposeYZXtoZXY(rho1,   Nx, Ny, Nz, rho2);
            transposeYZXtoZXY(temp,   Nx, Ny, Nz, rhoU2);
            transposeYZXtoZXY(temp2,  Nx, Ny, Nz, rhoV2);
            transposeYZXtoZXY(temp3,  Nx, Ny, Nz, rhoW2);
            transposeYZXtoZXY(temp4,  Nx, Ny, Nz, rhoE2);

            //filter in the Y direction
            filtZ->filterField(rho2,  rho1);
            filtZ->filterField(rhoU2, temp);
            filtZ->filterField(rhoV2, temp2);
            filtZ->filterField(rhoW2, temp3);
            filtZ->filterField(rhoE2, temp4);

            //get us back to XYZ from ZXY
            transposeZXYtoXYZ(rho1,   Nx, Ny, Nz, rho2);
            transposeZXYtoXYZ(temp,   Nx, Ny, Nz, rhoU1);
            transposeZXYtoXYZ(temp2,  Nx, Ny, Nz, rhoV1);
            transposeZXYtoXYZ(temp3,  Nx, Ny, Nz, rhoW1);
            transposeZXYtoXYZ(temp4,  Nx, Ny, Nz, rhoE1);

            //from being cute with memory allocation need to copy rho2 to rho1
            memcpy(rho2, rho1, sizeof(double)*Nx*Ny*Nz);

        }else if(filterTimeStep%3 == 2){

            //Here we'll do Y->Z->X     

            //Do the transpose to YZX space first
            transposeXYZtoYZX(rho2,   Nx, Ny, Nz, rho1);
            transposeXYZtoYZX(rhoU2,  Nx, Ny, Nz, temp);
            transposeXYZtoYZX(rhoV2,  Nx, Ny, Nz, temp2);
            transposeXYZtoYZX(rhoW2,  Nx, Ny, Nz, temp3);
            transposeXYZtoYZX(rhoE2,  Nx, Ny, Nz, temp4);

            //filter in the Y direction
            filtY->filterField(rho1,  rho2);
            filtY->filterField(temp,  rhoU2);
            filtY->filterField(temp2, rhoV2);
            filtY->filterField(temp3, rhoW2);
            filtY->filterField(temp4, rhoE2);

            //Move to ZXY space next
            transposeYZXtoZXY(rho2,   Nx, Ny, Nz, rho1);
            transposeYZXtoZXY(rhoU2,  Nx, Ny, Nz, temp);
            transposeYZXtoZXY(rhoV2,  Nx, Ny, Nz, temp2);
            transposeYZXtoZXY(rhoW2,  Nx, Ny, Nz, temp3);
            transposeYZXtoZXY(rhoE2,  Nx, Ny, Nz, temp4);

            //filter in the Z direction
            filtZ->filterField(rho1,  rho2);
            filtZ->filterField(temp,  rhoU2);
            filtZ->filterField(temp2, rhoV2);
            filtZ->filterField(temp3, rhoW2);
            filtZ->filterField(temp4, rhoE2);

            //transpose from ZXY to XYZ
            transposeZXYtoXYZ(rho2,   Nx, Ny, Nz, rho1);
            transposeZXYtoXYZ(rhoU2,  Nx, Ny, Nz, temp);
            transposeZXYtoXYZ(rhoV2,  Nx, Ny, Nz, temp2);
            transposeZXYtoXYZ(rhoW2,  Nx, Ny, Nz, temp3);
            transposeZXYtoXYZ(rhoE2,  Nx, Ny, Nz, temp4);

            //filter in the X direction
            filtX->filterField(rho1,  rho2);
            filtX->filterField(temp,  rhoU1);
            filtX->filterField(temp2, rhoV1);
            filtX->filterField(temp3, rhoW1);
            filtX->filterField(temp4, rhoE1);

            //from being cute with memory allocation need to copy rho2 to rho1
            memcpy(rho2, rho1, sizeof(double)*Nx*Ny*Nz);


        }else{
            //Here we'll do Z->X->Y     
            //Do the transpose to ZXY space first
            transposeXYZtoZXY(rho2,   Nx, Ny, Nz, rho1);
            transposeXYZtoZXY(rhoU2,  Nx, Ny, Nz, temp);
            transposeXYZtoZXY(rhoV2,  Nx, Ny, Nz, temp2);
            transposeXYZtoZXY(rhoW2,  Nx, Ny, Nz, temp3);
            transposeXYZtoZXY(rhoE2,  Nx, Ny, Nz, temp4);

            //filter in the Z direction
            filtZ->filterField(rho1,  rho2);
            filtZ->filterField(temp,  rhoU2);
            filtZ->filterField(temp2, rhoV2);
            filtZ->filterField(temp3, rhoW2);
            filtZ->filterField(temp4, rhoE2);

            //Move to XYZ space next
            transposeZXYtoXYZ(rho2,   Nx, Ny, Nz, rho1);
            transposeZXYtoXYZ(rhoU2,  Nx, Ny, Nz, temp);
            transposeZXYtoXYZ(rhoV2,  Nx, Ny, Nz, temp2);
            transposeZXYtoXYZ(rhoW2,  Nx, Ny, Nz, temp3);
            transposeZXYtoXYZ(rhoE2,  Nx, Ny, Nz, temp4);

            //filter in the X direction
            filtX->filterField(rho1,  rho2);
            filtX->filterField(temp,  rhoU2);
            filtX->filterField(temp2, rhoV2);
            filtX->filterField(temp3, rhoW2);
            filtX->filterField(temp4, rhoE2);

            //transpose from XYZ to YZX
            transposeXYZtoYZX(rho2,   Nx, Ny, Nz, rho1);
            transposeXYZtoYZX(rhoU2,  Nx, Ny, Nz, temp);
            transposeXYZtoYZX(rhoV2,  Nx, Ny, Nz, temp2);
            transposeXYZtoYZX(rhoW2,  Nx, Ny, Nz, temp3);
            transposeXYZtoYZX(rhoE2,  Nx, Ny, Nz, temp4);

            //filter in the Y direction
            filtY->filterField(rho1,  rho2);
            filtY->filterField(temp,  rhoU2);
            filtY->filterField(temp2, rhoV2);
            filtY->filterField(temp3, rhoW2);
            filtY->filterField(temp4, rhoE2);

            //have an extra step here to go from YZX to XYZ
            transposeYZXtoXYZ(rho2,   Nx, Ny, Nz, rho1);
            transposeYZXtoXYZ(rhoU2,  Nx, Ny, Nz, rhoU1);
            transposeYZXtoXYZ(rhoV2,  Nx, Ny, Nz, rhoV1);
            transposeYZXtoXYZ(rhoW2,  Nx, Ny, Nz, rhoW1);
            transposeYZXtoXYZ(rhoE2,  Nx, Ny, Nz, rhoE1);

        }
 
    //If not filtering, need to copy the solution over to the *1 variables
    }else{

	memcpy(rho2,  rho1,  sizeof(double)*N);
	memcpy(rhoU2, rhoU1, sizeof(double)*N);
	memcpy(rhoV2, rhoV1, sizeof(double)*N);
	memcpy(rhoW2, rhoW1, sizeof(double)*N);
	memcpy(rhoE2, rhoE1, sizeof(double)*N);

    }
};

void CSolver::updateNonConservedData(){
    if(rkStep == 1 || rkStep == 2 || rkStep == 3){

	ig->solveU(rhok, rhoUk, U);
	ig->solveU(rhok, rhoVk, V);
	ig->solveU(rhok, rhoWk, W);
	ig->solvep(rhok, rhoEk, U, V, W, p);
	ig->solveT(rhok, p, T);
	ig->solveMu(T, mu);
	ig->solveAmu(T, Amu);
	ig->solveSOS(rhok, p, sos);

    }else if(rkStep == 4){

	ig->solveU(rho1, rhoU1, U);
	ig->solveU(rho1, rhoV1, V);
	ig->solveU(rho1, rhoW1, W);
	ig->solvep(rho1, rhoE1, U, V, W, p);
	ig->solveT(rho1, p, T);
	ig->solveMu(T, mu);
	ig->solveAmu(T, Amu);
	ig->solveSOS(rho1, p, sos);

    }
}


void CSolver::updateSponge(){
    if(spongeFlag){
	double eps = 1.0/(spg->avgT/ts->dt + 1.0);
	FOR_XYZ spg->spongeRhoAvg[ip]  += eps*(rho1[ip]  - spg->spongeRhoAvg[ip]);	
	FOR_XYZ spg->spongeRhoUAvg[ip] += eps*(rhoU1[ip] - spg->spongeRhoUAvg[ip]);	
	FOR_XYZ spg->spongeRhoVAvg[ip] += eps*(rhoV1[ip] - spg->spongeRhoVAvg[ip]);	
	FOR_XYZ spg->spongeRhoWAvg[ip] += eps*(rhoW1[ip] - spg->spongeRhoWAvg[ip]);	
	FOR_XYZ spg->spongeRhoEAvg[ip] += eps*(rhoE1[ip] - spg->spongeRhoEAvg[ip]);	
	
	FOR_XYZ spg->spongeRhoEAvg[ip] = spg->epsP*spg->spongeRhoEAvg[ip] + (1.0 -  spg->epsP)*(spg->spongeP/(ig->gamma-1.0) \
					 + 0.5*(spg->spongeRhoUAvg[ip]*spg->spongeRhoUAvg[ip] + spg->spongeRhoVAvg[ip]*spg->spongeRhoVAvg[ip] \
					 + spg->spongeRhoWAvg[ip]*spg->spongeRhoWAvg[ip])/spg->spongeRhoAvg[ip]);

	
        if(bc->bcX0 == BC::SPONGE){

        }

        if(bc->bcX1 == BC::SPONGE){

        }   

        if(bc->bcY0 == BC::SPONGE){

        }

        if(bc->bcY1 == BC::SPONGE){

        }

        if(bc->bcZ0 == BC::SPONGE){

        }

        if(bc->bcZ1 == BC::SPONGE){

        }

    }
};


void CSolver::test(){

    derivX->calc1stDerivField(U, Ux);
}






