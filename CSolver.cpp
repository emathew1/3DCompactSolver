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

void CSolver::preStepDerivatives(int rkStep){

    ///////////////////
    // X-DERIVATIVES //
    ///////////////////

    //First we'll do all of the X-Direction derivatives since we're in XYZ order
 
    FOR_Z{
	FOR_Y{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn     =   &U[k*Nx*Ny + j*Nx];
	   dataOutd1  =  &Ux[k*Nx*Ny + j*Nx];
	   dataOutd2  = &Uxx[k*Nx*Ny + j*Nx];
	   derivX->calc1stDeriv(dataIn, dataOutd1);
	   derivX->calc2ndDeriv(dataIn, dataOutd2);
	}
    }
  
    FOR_Z{
	FOR_Y{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn    =   &V[k*Nx*Ny + j*Nx];
	   dataOutd1 =  &Vx[k*Nx*Ny + j*Nx];
	   dataOutd2 = &Vxx[k*Nx*Ny + j*Nx];
	   derivX->calc1stDeriv(dataIn, dataOutd1);
	   derivX->calc2ndDeriv(dataIn, dataOutd2);
	}
    }

    FOR_Z{
	FOR_Y{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn    =   &W[k*Nx*Ny + j*Nx];
	   dataOutd1 =  &Wx[k*Nx*Ny + j*Nx];
	   dataOutd2 = &Wxx[k*Nx*Ny + j*Nx];
	   derivX->calc1stDeriv(dataIn, dataOutd1);
	   derivX->calc2ndDeriv(dataIn, dataOutd2);
	}
    }

    FOR_Z{
	FOR_Y{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn    =   &T[k*Nx*Ny + j*Nx];
	   dataOutd1 =  &Tx[k*Nx*Ny + j*Nx];
	   dataOutd2 = &Txx[k*Nx*Ny + j*Nx];
	   derivX->calc1stDeriv(dataIn, dataOutd1);
	   derivX->calc2ndDeriv(dataIn, dataOutd2);
	}
    }

    //Calculate the Euler component...

//!!!!! 
    //WHICH RHOU! etc. NEEDS TO CHANGE WITH THE RKSTEP!!!!
//!!!!!

    FOR_Z{
	FOR_Y{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn    =      &rhoU1[k*Nx*Ny + j*Nx];
	   dataOutd1 = &contEulerX[k*Nx*Ny + j*Nx];
	   derivX->calc1stDeriv(dataIn, dataOutd1);
	}
    }

    FOR_XYZ temp[ip] = rhoU1[ip]*U[ip] + p[ip];
    

    FOR_Z{
	FOR_Y{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn    =       &temp[k*Nx*Ny + j*Nx];
	   dataOutd1 = &momXEulerX[k*Nx*Ny + j*Nx];
	   derivX->calc1stDeriv(dataIn, dataOutd1);
	}
    }

    FOR_XYZ temp[ip] = rhoV1[ip]*U[ip];
    

    FOR_Z{
	FOR_Y{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn    =       &temp[k*Nx*Ny + j*Nx];
	   dataOutd1 = &momYEulerX[k*Nx*Ny + j*Nx];
	   derivX->calc1stDeriv(dataIn, dataOutd1);
	}
    }

    FOR_XYZ temp[ip] = rhoW1[ip]*U[ip];
    

    FOR_Z{
	FOR_Y{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn    =       &temp[k*Nx*Ny + j*Nx];
	   dataOutd1 = &momZEulerX[k*Nx*Ny + j*Nx];
	   derivX->calc1stDeriv(dataIn, dataOutd1);
	}
    }


    FOR_XYZ temp[ip] = rhoE1[ip]*U[ip] + U[ip]*p[ip];
    

    FOR_Z{
	FOR_Y{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn    =       &temp[k*Nx*Ny + j*Nx];
	   dataOutd1 = &engyEulerX[k*Nx*Ny + j*Nx];
	   derivX->calc1stDeriv(dataIn, dataOutd1);
	}
    }


    transposeXYZtoYZX(rho1,  Nx, Ny, Nz, transRho);
    transposeXYZtoYZX(rhoU1, Nx, Ny, Nz, transRhoU);
    transposeXYZtoYZX(rhoV1, Nx, Ny, Nz, transRhoV);
    transposeXYZtoYZX(rhoW1, Nx, Ny, Nz, transRhoW);
    transposeXYZtoYZX(rhoE1, Nx, Ny, Nz, transRhoE);
    transposeXYZtoYZX(Ux,    Nx, Ny, Nz, transUx);
    transposeXYZtoYZX(Vx,    Nx, Ny, Nz, transVx);
    transposeXYZtoYZX(Wx,    Nx, Ny, Nz, transWx);


    ///////////////////
    // Y-DERIVATIVES //
    ///////////////////

    //Now recalculate properties in the new space
    FOR_XYZ{
	U[ip] = transRhoU[ip]/transRho[ip];
	V[ip] = transRhoV[ip]/transRho[ip];
	W[ip] = transRhoW[ip]/transRho[ip];
    }
    ig->solvep(transRho, transRhoE, U, V, W, p);
    ig->solveT(transRho, p, T);

    FOR_X{
	FOR_Z{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn     =   &U[i*Nz*Ny + k*Ny];
	   dataOutd1  =  &Uy[i*Nz*Ny + k*Ny];
	   dataOutd2  = &Uyy[i*Nz*Ny + k*Ny];
	   derivY->calc1stDeriv(dataIn, dataOutd1);
	   derivY->calc2ndDeriv(dataIn, dataOutd2);
	}
    }
 
    FOR_X{
	FOR_Z{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn     =   &V[i*Nz*Ny + k*Ny];
	   dataOutd1  =  &Vy[i*Nz*Ny + k*Ny];
	   dataOutd2  = &Vyy[i*Nz*Ny + k*Ny];
	   derivY->calc1stDeriv(dataIn, dataOutd1);
	   derivY->calc2ndDeriv(dataIn, dataOutd2);
	}
    }
 
    FOR_X{
	FOR_Z{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn     =   &W[i*Nz*Ny + k*Ny];
	   dataOutd1  =  &Wy[i*Nz*Ny + k*Ny];
	   dataOutd2  = &Wyy[i*Nz*Ny + k*Ny];
	   derivY->calc1stDeriv(dataIn, dataOutd1);
	   derivY->calc2ndDeriv(dataIn, dataOutd2);
	}
    }
 
    FOR_X{
	FOR_Z{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn     =   &T[i*Nz*Ny + k*Ny];
	   dataOutd1  =  &Ty[i*Nz*Ny + k*Ny];
	   dataOutd2  = &Tyy[i*Nz*Ny + k*Ny];
	   derivY->calc1stDeriv(dataIn, dataOutd1);
	   derivY->calc2ndDeriv(dataIn, dataOutd2);
	}
    }

  
    FOR_X{
	FOR_Z{
	   double *dataIn, *dataOutd1;
	   dataIn     =  &transUx[i*Nz*Ny + k*Ny];
	   dataOutd1  =      &Uxy[i*Nz*Ny + k*Ny];
	   derivY->calc1stDeriv(dataIn, dataOutd1);
	}
    }  


    FOR_X{
	FOR_Z{
	   double *dataIn, *dataOutd1;
	   dataIn     =  &transVx[i*Nz*Ny + k*Ny];
	   dataOutd1  =      &Vxy[i*Nz*Ny + k*Ny];
	   derivY->calc1stDeriv(dataIn, dataOutd1);
	}
    }  

    FOR_X{
	FOR_Z{
	   double *dataIn, *dataOutd1;
	   dataIn     =  &transWx[i*Nz*Ny + k*Ny];
	   dataOutd1  =      &Wxy[i*Nz*Ny + k*Ny];
	   derivY->calc1stDeriv(dataIn, dataOutd1);
	}
    } 


    FOR_X{
	FOR_Z{
	   double *dataIn, *dataOutd1;
	   dataIn     =  &transRhoV[i*Nz*Ny + k*Ny];
	   dataOutd1  = &contEulerY[i*Nz*Ny + k*Ny];
	   derivY->calc1stDeriv(dataIn, dataOutd1);
	}
    } 
 
    
    FOR_XYZ temp[ip] = transRhoU[ip]*V[ip];
    

    FOR_X{
        FOR_Z{
           double *dataIn, *dataOutd1;
           dataIn     =       &temp[i*Nz*Ny + k*Ny];
           dataOutd1  = &momXEulerY[i*Nz*Ny + k*Ny];
           derivY->calc1stDeriv(dataIn, dataOutd1);
        }
    }

    FOR_XYZ temp[ip] = transRhoV[ip]*V[ip] + p[ip];
    

    FOR_X{
        FOR_Z{
           double *dataIn, *dataOutd1;
           dataIn     =       &temp[i*Nz*Ny + k*Ny];
           dataOutd1  = &momYEulerY[i*Nz*Ny + k*Ny];
           derivY->calc1stDeriv(dataIn, dataOutd1);
        }
    }


    FOR_XYZ temp[ip] = transRhoW[ip]*V[ip];
    

    FOR_X{
        FOR_Z{
           double *dataIn, *dataOutd1;
           dataIn     =       &temp[i*Nz*Ny + k*Ny];
           dataOutd1  = &momZEulerY[i*Nz*Ny + k*Ny];
           derivY->calc1stDeriv(dataIn, dataOutd1);
        }
    }

    FOR_XYZ temp[ip] = transRhoE[ip]*V[ip] + V[ip]*p[ip];
    

    FOR_X{
        FOR_Z{
           double *dataIn, *dataOutd1;
           dataIn     =       &temp[i*Nz*Ny + k*Ny];
           dataOutd1  = &engyEulerY[i*Nz*Ny + k*Ny];
           derivY->calc1stDeriv(dataIn, dataOutd1);
        }
    }
    
//!!!!!!
    //NEED TO TRANSPOSE "INPLACE" THESE EULER DERIVATIVES!!!!
    //AND THE CROSS DERIVATIVES!!!! 
//!!!!!!   

    //Moving data to ZXY

    //Get the original conserved data from XYZ to ZXY
    transposeXYZtoZXY(rho1,  Nx, Ny, Nz, transRho);
    transposeXYZtoZXY(rhoU1, Nx, Ny, Nz, transRhoU);
    transposeXYZtoZXY(rhoV1, Nx, Ny, Nz, transRhoV);
    transposeXYZtoZXY(rhoW1, Nx, Ny, Nz, transRhoW);
    transposeXYZtoZXY(rhoE1, Nx, Ny, Nz, transRhoE);

    //Move the Y Derivative data from YZX to ZXY
    transposeYZXtoZXY(Uy, Nx, Ny, Nz, transUy);
    transposeYZXtoZXY(Vy, Nx, Ny, Nz, transVy);
    transposeYZXtoZXY(Wy, Nx, Ny, Nz, transWy);

    //Move the X Derivative data from XYZ to ZXY
    transposeXYZtoZXY(Ux, Nx, Ny, Nz, transUx);
    transposeXYZtoZXY(Vx, Nx, Ny, Nz, transVx);
    transposeXYZtoZXY(Wx, Nx, Ny, Nz, transWx);



    //Moving Data from YZX to XYZ

    memcpy(temp, Uy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Uy);
    memcpy(temp, Uyy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Uyy);

    memcpy(temp, Vy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Vy);
    memcpy(temp, Vyy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Vyy);

    memcpy(temp, Wy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Wy);
    memcpy(temp, Wyy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Wyy);

    memcpy(temp, Ty, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Ty);
    memcpy(temp, Tyy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Tyy);

    memcpy(temp, Uxy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Uxy);
    memcpy(temp, Vxy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Vxy);
    memcpy(temp, Wxy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Wxy);

    memcpy(temp, contEulerY, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, contEulerY);

    memcpy(temp, momXEulerY, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, momXEulerY);

    memcpy(temp, momYEulerY, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, momYEulerY);

    memcpy(temp, momZEulerY, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, momZEulerY);

    memcpy(temp, engyEulerY, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, engyEulerY);


    ///////////////////
    // Z-DERIVATIVES //
    ///////////////////

    //Now recalculate properties in the new space
    FOR_XYZ{
	U[ip] = transRhoU[ip]/transRho[ip];
	V[ip] = transRhoV[ip]/transRho[ip];
	W[ip] = transRhoW[ip]/transRho[ip];
    }
    ig->solvep(transRho, transRhoE, U, V, W, p);
    ig->solveT(transRho, p, T);

    
    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn     =   &U[j*Nz*Nx + i*Nz];
	   dataOutd1  =  &Uz[j*Nz*Nx + i*Nz];
	   dataOutd2  = &Uzz[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	   derivZ->calc2ndDeriv(dataIn, dataOutd2);
	}
    }

    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn     =   &V[j*Nz*Nx + i*Nz];
	   dataOutd1  =  &Vz[j*Nz*Nx + i*Nz];
	   dataOutd2  = &Vzz[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	   derivZ->calc2ndDeriv(dataIn, dataOutd2);
	}
    }

    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn     =   &W[j*Nz*Nx + i*Nz];
	   dataOutd1  =  &Wz[j*Nz*Nx + i*Nz];
	   dataOutd2  = &Wzz[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	   derivZ->calc2ndDeriv(dataIn, dataOutd2);
	}
    }

    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1, *dataOutd2;
	   dataIn     =   &T[j*Nz*Nx + i*Nz];
	   dataOutd1  =  &Tz[j*Nz*Nx + i*Nz];
	   dataOutd2  = &Tzz[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	   derivZ->calc2ndDeriv(dataIn, dataOutd2);
	}
    }

    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1;
	   dataIn     =  &transUx[j*Nz*Nx + i*Nz];
	   dataOutd1  =      &Uxz[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	}
    }

    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1;
	   dataIn     =  &transVx[j*Nz*Nx + i*Nz];
	   dataOutd1  =      &Vxz[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	}
    }

    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1;
	   dataIn     =  &transWx[j*Nz*Nx + i*Nz];
	   dataOutd1  =      &Wxz[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	}
    }


    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1;
	   dataIn     =  &transUy[j*Nz*Nx + i*Nz];
	   dataOutd1  =      &Uyz[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	}
    }


    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1;
	   dataIn     =  &transVy[j*Nz*Nx + i*Nz];
	   dataOutd1  =      &Vyz[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	}
    }


    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1;
	   dataIn     =  &transWy[j*Nz*Nx + i*Nz];
	   dataOutd1  =      &Wyz[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	}
    }


    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1;
	   dataIn     =  &transRhoW[j*Nz*Nx + i*Nz];
	   dataOutd1  = &contEulerZ[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	}
    }

    
    FOR_XYZ temp[ip] = transRhoU[ip]*W[ip];

    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1;
	   dataIn     =        &temp[j*Nz*Nx + i*Nz];
	   dataOutd1  =  &momXEulerZ[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	}
    }


    FOR_XYZ temp[ip] = transRhoV[ip]*W[ip];

    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1;
	   dataIn     =        &temp[j*Nz*Nx + i*Nz];
	   dataOutd1  =  &momYEulerZ[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	}
    }


    FOR_XYZ temp[ip] = transRhoW[ip]*W[ip] + p[ip];

    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1;
	   dataIn     =        &temp[j*Nz*Nx + i*Nz];
	   dataOutd1  =  &momZEulerZ[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	}
    }



    FOR_XYZ temp[ip] = transRhoE[ip]*W[ip] + W[ip]*p[ip];

    FOR_Y{
	FOR_X{
	   double *dataIn, *dataOutd1;
	   dataIn     =        &temp[j*Nz*Nx + i*Nz];
	   dataOutd1  =  &engyEulerZ[j*Nz*Nx + i*Nz];
	   derivZ->calc1stDeriv(dataIn, dataOutd1);
	}
    }



    //Moving all the data back to XYZ
    memcpy(temp, Uz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Uz);
    memcpy(temp, Uzz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Uzz);

    memcpy(temp, Vz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Vz);
    memcpy(temp, Vzz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Vzz);

    memcpy(temp, Wz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Wz);
    memcpy(temp, Wzz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Wzz);

    memcpy(temp, Tz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Tz);
    memcpy(temp, Tzz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Tzz);

    memcpy(temp, Uxz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Uxz);
    memcpy(temp, Vxz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Vxz);
    memcpy(temp, Wxz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Wxz);

    memcpy(temp, Uyz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Uyz);
    memcpy(temp, Vyz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Vyz);
    memcpy(temp, Wyz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Wyz);

    memcpy(temp, contEulerZ, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, contEulerZ);
    memcpy(temp, momXEulerZ, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, momXEulerZ);
    memcpy(temp, momYEulerZ, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, momYEulerZ);
    memcpy(temp, momZEulerZ, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, momZEulerZ);
    memcpy(temp, engyEulerZ, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, engyEulerZ);



}


void CSolver::preStepDerivatives2(int rkStep){

    //TODO
    //WHICH RHOU! etc. NEEDS TO CHANGE WITH THE RKSTEP!!!!

 
    ///////////////////
    // X-DERIVATIVES //
    ///////////////////

    //First we'll do all of the X-Direction derivatives since we're in XYZ order

    //Calculate the Euler Components of the equations... 
    FOR_XYZ temp[ip]  = rhoU1[ip]*U[ip] + p[ip];
    FOR_XYZ temp2[ip] = rhoV1[ip]*U[ip];
    FOR_XYZ temp3[ip] = rhoW1[ip]*U[ip];
    FOR_XYZ temp4[ip] = rhoE1[ip]*U[ip] + U[ip]*p[ip];

    //Calculate the stuff needed for viscous derivatives
    derivX->calc1stDerivField(U, Ux);
    derivX->calc2ndDerivField(U, Uxx);
    derivX->calc1stDerivField(V, Vx);
    derivX->calc2ndDerivField(V, Vxx);
    derivX->calc1stDerivField(W, Wx);
    derivX->calc2ndDerivField(W, Wxx);
    derivX->calc1stDerivField(T, Tx);
    derivX->calc2ndDerivField(T, Txx);

    //Compute the Euler Derivatives
    derivX->calc1stDerivField(rhoU1, contEulerX);
    derivX->calc1stDerivField(temp,  momXEulerX);
    derivX->calc1stDerivField(temp2, momYEulerX);
    derivX->calc1stDerivField(temp3, momZEulerX);
    derivX->calc1stDerivField(temp4, engyEulerX);

    //Do the transposes
    transposeXYZtoYZX(rho1,  Nx, Ny, Nz, transRho);
    transposeXYZtoYZX(rhoU1, Nx, Ny, Nz, transRhoU);
    transposeXYZtoYZX(rhoV1, Nx, Ny, Nz, transRhoV);
    transposeXYZtoYZX(rhoW1, Nx, Ny, Nz, transRhoW);
    transposeXYZtoYZX(rhoE1, Nx, Ny, Nz, transRhoE);
    transposeXYZtoYZX(Ux,    Nx, Ny, Nz, transUx);
    transposeXYZtoYZX(Vx,    Nx, Ny, Nz, transVx);
    transposeXYZtoYZX(Wx,    Nx, Ny, Nz, transWx);

    ///////////////////
    // Y-DERIVATIVES //
    ///////////////////

    //Now recalculate properties in the new space
    FOR_XYZ U[ip] = transRhoU[ip]/transRho[ip];
    FOR_XYZ V[ip] = transRhoV[ip]/transRho[ip];
    FOR_XYZ W[ip] = transRhoW[ip]/transRho[ip];

    ig->solvep(transRho, transRhoE, U, V, W, p);
    ig->solveT(transRho, p, T);

    //Calculate the stuff for Euler derivatives in new space
    FOR_XYZ temp[ip]  = transRhoU[ip]*V[ip];
    FOR_XYZ temp2[ip] = transRhoV[ip]*V[ip] + p[ip];
    FOR_XYZ temp3[ip] = transRhoW[ip]*V[ip];
    FOR_XYZ temp4[ip] = transRhoE[ip]*V[ip] + V[ip]*p[ip];

    //Calculate Viscous Derivatives
    derivY->calc1stDerivField(U, Uy);
    derivY->calc2ndDerivField(U, Uyy);
    derivY->calc1stDerivField(V, Vy);
    derivY->calc2ndDerivField(V, Vyy);
    derivY->calc1stDerivField(W, Wy);
    derivY->calc2ndDerivField(W, Wyy);
    derivY->calc1stDerivField(T, Ty);
    derivY->calc2ndDerivField(T, Tyy);
    derivY->calc1stDerivField(transUx, Uxy);
    derivY->calc1stDerivField(transVx, Vxy);
    derivY->calc1stDerivField(transWx, Wxy);

    //Calculate the Euler Derivatives
    derivY->calc1stDerivField(transRhoV, contEulerY);
    derivY->calc1stDerivField(temp, 	 momXEulerY);
    derivY->calc1stDerivField(temp2,	 momYEulerY);
    derivY->calc1stDerivField(temp3,	 momZEulerY);
    derivY->calc1stDerivField(temp4,	 engyEulerY);

    
    //Moving data to ZXY

    //Get the original conserved data from XYZ to ZXY
    transposeXYZtoZXY(rho1,  Nx, Ny, Nz, transRho);
    transposeXYZtoZXY(rhoU1, Nx, Ny, Nz, transRhoU);
    transposeXYZtoZXY(rhoV1, Nx, Ny, Nz, transRhoV);
    transposeXYZtoZXY(rhoW1, Nx, Ny, Nz, transRhoW);
    transposeXYZtoZXY(rhoE1, Nx, Ny, Nz, transRhoE);

    //Move the Y Derivative data from YZX to ZXY
    transposeYZXtoZXY(Uy, Nx, Ny, Nz, transUy);
    transposeYZXtoZXY(Vy, Nx, Ny, Nz, transVy);
    transposeYZXtoZXY(Wy, Nx, Ny, Nz, transWy);

    //Move the X Derivative data from XYZ to ZXY
    transposeXYZtoZXY(Ux, Nx, Ny, Nz, transUx);
    transposeXYZtoZXY(Vx, Nx, Ny, Nz, transVx);
    transposeXYZtoZXY(Wx, Nx, Ny, Nz, transWx);


    //Moving Data from YZX to XYZ

    memcpy(temp, Uy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Uy);
    memcpy(temp, Uyy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Uyy);

    memcpy(temp, Vy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Vy);
    memcpy(temp, Vyy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Vyy);

    memcpy(temp, Wy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Wy);
    memcpy(temp, Wyy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Wyy);

    memcpy(temp, Ty, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Ty);
    memcpy(temp, Tyy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Tyy);

    memcpy(temp, Uxy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Uxy);
    memcpy(temp, Vxy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Vxy);
    memcpy(temp, Wxy, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, Wxy);

    memcpy(temp, contEulerY, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, contEulerY);

    memcpy(temp, momXEulerY, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, momXEulerY);

    memcpy(temp, momYEulerY, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, momYEulerY);

    memcpy(temp, momZEulerY, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, momZEulerY);

    memcpy(temp, engyEulerY, sizeof(double)*Nx*Ny*Nz);
    transposeYZXtoXYZ(temp, Nx, Ny, Nz, engyEulerY);

    ///////////////////
    // Z-DERIVATIVES //
    ///////////////////

    //Now recalculate properties in the new space
    FOR_XYZ U[ip] = transRhoU[ip]/transRho[ip];
    FOR_XYZ V[ip] = transRhoV[ip]/transRho[ip];
    FOR_XYZ W[ip] = transRhoW[ip]/transRho[ip];

    ig->solvep(transRho, transRhoE, U, V, W, p);
    ig->solveT(transRho, p, T);

    //Calculate the stuff for the Euler Derivatives
    FOR_XYZ temp[ip]  = transRhoU[ip]*W[ip];
    FOR_XYZ temp2[ip] = transRhoV[ip]*W[ip] + p[ip];
    FOR_XYZ temp3[ip] = transRhoW[ip]*W[ip] + p[ip];
    FOR_XYZ temp4[ip] = transRhoE[ip]*W[ip] + W[ip]*p[ip];

    //Calculate the viscous derivatives
    derivZ->calc1stDerivField(U, Uz);
    derivZ->calc2ndDerivField(U, Uzz);

    derivZ->calc1stDerivField(V, Vz);
    derivZ->calc2ndDerivField(V, Vzz);

    derivZ->calc1stDerivField(W, Wz);
    derivZ->calc2ndDerivField(W, Wzz);

    derivZ->calc1stDerivField(T, Tz);
    derivZ->calc2ndDerivField(T, Tzz);

    derivZ->calc1stDerivField(transUx, Uxz);
    derivZ->calc1stDerivField(transVx, Vxz);
    derivZ->calc1stDerivField(transWx, Wxz);

    derivZ->calc1stDerivField(transUy, Uyz);
    derivZ->calc1stDerivField(transVy, Vyz);
    derivZ->calc1stDerivField(transWy, Wyz);

    //Calculate the Euler Derivatives
    derivZ->calc1stDerivField(transRhoW, contEulerZ);
    derivZ->calc1stDerivField(temp,	 contEulerZ);
    derivZ->calc1stDerivField(temp2,	 contEulerZ);
    derivZ->calc1stDerivField(temp3, 	 contEulerZ);
    derivZ->calc1stDerivField(temp4,	 contEulerZ);

    //Moving all the data back to XYZ
    memcpy(temp, Uz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Uz);
    memcpy(temp, Uzz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Uzz);

    memcpy(temp, Vz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Vz);
    memcpy(temp, Vzz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Vzz);

    memcpy(temp, Wz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Wz);
    memcpy(temp, Wzz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Wzz);

    memcpy(temp, Tz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Tz);
    memcpy(temp, Tzz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Tzz);

    memcpy(temp, Uxz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Uxz);
    memcpy(temp, Vxz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Vxz);
    memcpy(temp, Wxz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Wxz);

    memcpy(temp, Uyz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Uyz);
    memcpy(temp, Vyz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Vyz);
    memcpy(temp, Wyz, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, Wyz);

    memcpy(temp, contEulerZ, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, contEulerZ);
    memcpy(temp, momXEulerZ, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, momXEulerZ);
    memcpy(temp, momYEulerZ, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, momYEulerZ);
    memcpy(temp, momZEulerZ, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, momZEulerZ);
    memcpy(temp, engyEulerZ, sizeof(double)*Nx*Ny*Nz);
    transposeZXYtoXYZ(temp, Nx, Ny, Nz, engyEulerZ);

    //Going back to original...
    FOR_XYZ U[ip] = rhoU1[ip]/rho1[ip];
    FOR_XYZ V[ip] = rhoV1[ip]/rho1[ip];
    FOR_XYZ W[ip] = rhoW1[ip]/rho1[ip];
    ig->solvep(rho1, rhoE1, U, V, W, p);
    ig->solveT(rho1, p, T);

}




