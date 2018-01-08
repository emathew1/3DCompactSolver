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

    //turbulence quantities
    turbdiss   = new double[Nx*Ny*Nz];
    uprime2    = new double[Nx*Ny*Nz];
    uiprime2    = new double[Nx*Ny*Nz];
    uvar       = new double[Nx*Ny*Nz];
    kineticEng = new double[Nx*Ny*Nz];
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
    #pragma omp parallel for
    FOR_XYZ rhoE1[ip] = ig->solverhoE(rho1[ip], p[ip], U[ip], V[ip], W[ip]);
    #pragma omp parallel for
    FOR_XYZ T[ip]     = ig->solveT(rho1[ip], p[ip]);
    #pragma omp parallel for
    FOR_XYZ mu[ip]    = ig->solveMu(T[ip]);
    #pragma omp parallel for
    FOR_XYZ Amu[ip]   = ig->solveAmu(T[ip]);
    #pragma omp parallel for
    FOR_XYZ sos[ip]   = ig->solveSOS(rho1[ip], p[ip]);

    //This is where we'll do the boundary condition specific stuff...
    bool wallBCFlag = false;

    //--------------------------------
    //No-Slip Wall Boundary Conditions
    //--------------------------------

    if(bc->bcX0 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        #pragma omp parallel for
        FOR_X0{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Xp1],
                                T[GET3DINDEX_XYZ_Xp2],
                                T[GET3DINDEX_XYZ_Xp3],
                                T[GET3DINDEX_XYZ_Xp4],
                                T[GET3DINDEX_XYZ_Xp5],
                                T[GET3DINDEX_XYZ_Xp6]);
        }END_FORX0
    }

    if(bc->bcX1 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        #pragma omp parallel for
        FOR_X1{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Xm1],
                                T[GET3DINDEX_XYZ_Xm2],
                                T[GET3DINDEX_XYZ_Xm3],
                                T[GET3DINDEX_XYZ_Xm4],
                                T[GET3DINDEX_XYZ_Xm5],
                                T[GET3DINDEX_XYZ_Xm6]);
        }END_FORX1
    }

    if(bc->bcY0 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        #pragma omp parallel for
        FOR_Y0{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Yp1],
                                T[GET3DINDEX_XYZ_Yp2],
                                T[GET3DINDEX_XYZ_Yp3],
                                T[GET3DINDEX_XYZ_Yp4],
                                T[GET3DINDEX_XYZ_Yp5],
                                T[GET3DINDEX_XYZ_Yp6]);
        }END_FORY0
    }

    if(bc->bcY1 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        #pragma omp parallel for
        FOR_Y1{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Ym1],
                                T[GET3DINDEX_XYZ_Ym2],
                                T[GET3DINDEX_XYZ_Ym3],
                                T[GET3DINDEX_XYZ_Ym4],
                                T[GET3DINDEX_XYZ_Ym5],
                                T[GET3DINDEX_XYZ_Ym6]);
        }END_FORY1
    }

    if(bc->bcZ0 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        #pragma omp parallel for
        FOR_Z0{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Zp1],
                                T[GET3DINDEX_XYZ_Zp2],
                                T[GET3DINDEX_XYZ_Zp3],
                                T[GET3DINDEX_XYZ_Zp4],
                                T[GET3DINDEX_XYZ_Zp5],
                                T[GET3DINDEX_XYZ_Zp6]);
        }END_FORZ0
    }

    if(bc->bcZ1 == BC::ADIABATIC_WALL){
	wallBCFlag = true;
        FOR_Z1{
            U[ip]  = 0.0;
            V[ip]  = 0.0;
            W[ip]  = 0.0;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = 0.0;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Zm1],
                                T[GET3DINDEX_XYZ_Zm2],
                                T[GET3DINDEX_XYZ_Zm3],
                                T[GET3DINDEX_XYZ_Zm4],
                                T[GET3DINDEX_XYZ_Zm5],
                                T[GET3DINDEX_XYZ_Zm6]);
        }END_FORZ1
    }

    //-------------------------------
    //Moving Wall Boundary Conditions
    //-------------------------------

    if(bc->bcX0 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        #pragma omp parallel for
        FOR_X0{
            U[ip]  = 0.0;
            V[ip]  = X0WallV;
            W[ip]  = X0WallW;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = rho1[ip]*X0WallV;
            rhoW1[ip] = rho1[ip]*X0WallW;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Xp1],
                                T[GET3DINDEX_XYZ_Xp2],
                                T[GET3DINDEX_XYZ_Xp3],
                                T[GET3DINDEX_XYZ_Xp4],
                                T[GET3DINDEX_XYZ_Xp5],
                                T[GET3DINDEX_XYZ_Xp6]);
        }END_FORX0
    }

    if(bc->bcX1 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        #pragma omp parallel for
        FOR_X1{
            U[ip]  = 0.0;
            V[ip]  = X1WallV;
            W[ip]  = X1WallW;
            rhoU1[ip] = 0.0;
            rhoV1[ip] = rho1[ip]*X1WallV;
            rhoW1[ip] = rho1[ip]*X1WallW;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Xm1],
                                T[GET3DINDEX_XYZ_Xm2],
                                T[GET3DINDEX_XYZ_Xm3],
                                T[GET3DINDEX_XYZ_Xm4],
                                T[GET3DINDEX_XYZ_Xm5],
                                T[GET3DINDEX_XYZ_Xm6]);
        }END_FORX1
    }

    if(bc->bcY0 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        #pragma omp parallel for
        FOR_Y0{
            U[ip]  = Y0WallU;
            V[ip]  = 0.0;
            W[ip]  = Y0WallW;
            rhoU1[ip] = rho1[ip]*Y0WallU;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = rho1[ip]*Y0WallW;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Yp1],
                                T[GET3DINDEX_XYZ_Yp2],
                                T[GET3DINDEX_XYZ_Yp3],
                                T[GET3DINDEX_XYZ_Yp4],
                                T[GET3DINDEX_XYZ_Yp5],
                                T[GET3DINDEX_XYZ_Yp6]);
        }END_FORY0
    }

    if(bc->bcY1 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        #pragma omp parallel for
        FOR_Y1{
            U[ip]  = Y1WallU;
            V[ip]  = 0.0;
            W[ip]  = Y1WallW;
            rhoU1[ip] = rho1[ip]*Y1WallU;
            rhoV1[ip] = 0.0;
            rhoW1[ip] = rho1[ip]*Y1WallW;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Ym1],
                                T[GET3DINDEX_XYZ_Ym2],
                                T[GET3DINDEX_XYZ_Ym3],
                                T[GET3DINDEX_XYZ_Ym4],
                                T[GET3DINDEX_XYZ_Ym5],
                                T[GET3DINDEX_XYZ_Ym6]);
        }END_FORY1
    }

    if(bc->bcZ0 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        #pragma omp parallel for
        FOR_Z0{
            U[ip]  = Z0WallU;
            V[ip]  = Z0WallV;
            W[ip]  = 0.0;
            rhoU1[ip] = rho1[ip]*Z0WallU;
            rhoV1[ip] = rho1[ip]*Z0WallV;
            rhoW1[ip] = 0.0;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Zp1],
                                T[GET3DINDEX_XYZ_Zp2],
                                T[GET3DINDEX_XYZ_Zp3],
                                T[GET3DINDEX_XYZ_Zp4],
                                T[GET3DINDEX_XYZ_Zp5],
                                T[GET3DINDEX_XYZ_Zp6]);
        }END_FORZ0
    }

    if(bc->bcZ1 == BC::MOVING_ADIABATIC_WALL){
	wallBCFlag = true;
        #pragma omp parallel for
        FOR_Z1{
            U[ip]  = Z1WallU;
            V[ip]  = Z1WallV;
            W[ip]  = 0.0;
            rhoU1[ip] = rho1[ip]*Z1WallU;
            rhoV1[ip] = rho1[ip]*Z1WallV;
            rhoW1[ip] = 0.0;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Zm1],
                                T[GET3DINDEX_XYZ_Zm2],
                                T[GET3DINDEX_XYZ_Zm3],
                                T[GET3DINDEX_XYZ_Zm4],
                                T[GET3DINDEX_XYZ_Zm5],
                                T[GET3DINDEX_XYZ_Zm6]);
        }END_FORZ1
    }




    if(wallBCFlag == true){
	//Need to update the pressure, sos, and rhoE fields at the boundaries with walls...
	FOR_XYZ{
	    p[ip]     = ig->solvep_idealgas(rho1[ip], T[ip]);
	    sos[ip]   = ig->solveSOS(rho1[ip],p[ip]);
	    rhoE1[ip] = ig->solverhoE(rho1[ip], p[ip], U[ip], V[ip], W[ip]);
	}
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
    #pragma omp parallel for
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

    //--------------------------------
    //No-slip wall boundary conditions
    //--------------------------------

    if(bc->bcX0 == BC::ADIABATIC_WALL){
	#pragma omp parallel for
	FOR_X0{
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Xp1],
				T[GET3DINDEX_XYZ_Xp2],
				T[GET3DINDEX_XYZ_Xp3],
				T[GET3DINDEX_XYZ_Xp4],
				T[GET3DINDEX_XYZ_Xp5],
				T[GET3DINDEX_XYZ_Xp6]);
	   p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	   rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
	}END_FORX0

    }

    if(bc->bcX1 == BC::ADIABATIC_WALL){
	#pragma omp parallel for
	FOR_X1{
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Xm1],
				T[GET3DINDEX_XYZ_Xm2],
				T[GET3DINDEX_XYZ_Xm3],
				T[GET3DINDEX_XYZ_Xm4],
				T[GET3DINDEX_XYZ_Xm5],
				T[GET3DINDEX_XYZ_Xm6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
	}END_FORX1
    }   

    if(bc->bcY0 == BC::ADIABATIC_WALL){
	#pragma omp parallel for
	FOR_Y0{
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Yp1],
				T[GET3DINDEX_XYZ_Yp2],
				T[GET3DINDEX_XYZ_Yp3],
				T[GET3DINDEX_XYZ_Yp4],
				T[GET3DINDEX_XYZ_Yp5],
				T[GET3DINDEX_XYZ_Yp6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
	}END_FORY0
    }

    if(bc->bcY1 == BC::ADIABATIC_WALL){
	#pragma omp parallel for
	FOR_Y1{
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Ym1],
				T[GET3DINDEX_XYZ_Ym2],
				T[GET3DINDEX_XYZ_Ym3],
				T[GET3DINDEX_XYZ_Ym4],
				T[GET3DINDEX_XYZ_Ym5],
				T[GET3DINDEX_XYZ_Ym6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
	}END_FORY1
    }

    if(bc->bcZ0 == BC::ADIABATIC_WALL){
	#pragma omp parallel for
	FOR_Z0{
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Zp1],
				T[GET3DINDEX_XYZ_Zp2],
				T[GET3DINDEX_XYZ_Zp3],
				T[GET3DINDEX_XYZ_Zp4],
				T[GET3DINDEX_XYZ_Zp5],
				T[GET3DINDEX_XYZ_Zp6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
	}END_FORZ0
    }

    if(bc->bcZ1 == BC::ADIABATIC_WALL){
	#pragma omp parallel for
	FOR_Z1{
	    U[ip]  = 0.0;
	    V[ip]  = 0.0;
	    W[ip]  = 0.0;
	    rhoUP[ip] = 0.0;
	    rhoVP[ip] = 0.0;
	    rhoWP[ip] = 0.0;
	    T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Zm1],
				T[GET3DINDEX_XYZ_Zm2],
				T[GET3DINDEX_XYZ_Zm3],
				T[GET3DINDEX_XYZ_Zm4],
				T[GET3DINDEX_XYZ_Zm5],
				T[GET3DINDEX_XYZ_Zm6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
	}END_FORZ1
    }

    //-------------------------------
    //Moving wall boundary conditions
    //-------------------------------

    if(bc->bcX0 == BC::MOVING_ADIABATIC_WALL){
        #pragma omp parallel for
        FOR_X0{
            U[ip]  = 0.0;
            V[ip]  = X0WallV;
            W[ip]  = X0WallW;
            rhoUP[ip] = 0.0;
            rhoVP[ip] = rhoP[ip]*X0WallV;
            rhoWP[ip] = rhoP[ip]*X0WallW;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Xp1],
                                T[GET3DINDEX_XYZ_Xp2],
                                T[GET3DINDEX_XYZ_Xp3],
                                T[GET3DINDEX_XYZ_Xp4],
                                T[GET3DINDEX_XYZ_Xp5],
                                T[GET3DINDEX_XYZ_Xp6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
        }END_FORX0
    }

    if(bc->bcX1 == BC::MOVING_ADIABATIC_WALL){
        #pragma omp parallel for
        FOR_X1{
            U[ip]  = 0.0;
            V[ip]  = X1WallV;
            W[ip]  = X1WallW;
            rhoUP[ip] = 0.0;
            rhoVP[ip] = rhoP[ip]*X1WallV;
            rhoWP[ip] = rhoP[ip]*X1WallW;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Xm1],
                                T[GET3DINDEX_XYZ_Xm2],
                                T[GET3DINDEX_XYZ_Xm3],
                                T[GET3DINDEX_XYZ_Xm4],
                                T[GET3DINDEX_XYZ_Xm5],
                                T[GET3DINDEX_XYZ_Xm6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
        }END_FORX1
    }

    if(bc->bcY0 == BC::MOVING_ADIABATIC_WALL){
        #pragma omp parallel for
        FOR_Y0{
            U[ip]  = Y0WallU;
            V[ip]  = 0.0;
            W[ip]  = Y0WallW;
            rhoUP[ip] = rhoP[ip]*Y0WallU;
            rhoVP[ip] = 0.0;
            rhoWP[ip] = rhoP[ip]*Y0WallW;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Yp1],
                                T[GET3DINDEX_XYZ_Yp2],
                                T[GET3DINDEX_XYZ_Yp3],
                                T[GET3DINDEX_XYZ_Yp4],
                                T[GET3DINDEX_XYZ_Yp5],
                                T[GET3DINDEX_XYZ_Yp6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
        }END_FORY0
    }

    if(bc->bcY1 == BC::MOVING_ADIABATIC_WALL){
        #pragma omp parallel for
        FOR_Y1{
            U[ip]  = Y1WallU;
            V[ip]  = 0.0;
            W[ip]  = Y1WallW;
            rhoUP[ip] = rhoP[ip]*Y1WallU;
            rhoVP[ip] = 0.0;
            rhoWP[ip] = rhoP[ip]*Y1WallW;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Ym1],
                                T[GET3DINDEX_XYZ_Ym2],
                                T[GET3DINDEX_XYZ_Ym3],
                                T[GET3DINDEX_XYZ_Ym4],
                                T[GET3DINDEX_XYZ_Ym5],
                                T[GET3DINDEX_XYZ_Ym6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
        }END_FORY1
    }

    if(bc->bcZ0 == BC::MOVING_ADIABATIC_WALL){
        #pragma omp parallel for
        FOR_Z0{
            U[ip]  = Z0WallU;
            V[ip]  = Z0WallV;
            W[ip]  = 0.0;
            rhoUP[ip] = rhoP[ip]*Z0WallU;
            rhoVP[ip] = rhoP[ip]*Z0WallV;
            rhoWP[ip] = 0.0;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Zp1],
                                T[GET3DINDEX_XYZ_Zp2],
                                T[GET3DINDEX_XYZ_Zp3],
                                T[GET3DINDEX_XYZ_Zp4],
                                T[GET3DINDEX_XYZ_Zp5],
                                T[GET3DINDEX_XYZ_Zp6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
        }END_FORZ0
    }

    if(bc->bcZ1 == BC::MOVING_ADIABATIC_WALL){
        #pragma omp parallel for
        FOR_Z1{
            U[ip]  = Z1WallU;
            V[ip]  = Z1WallV;
            W[ip]  = 0.0;
            rhoUP[ip] = rhoP[ip]*Z1WallU;
            rhoVP[ip] = rhoP[ip]*Z1WallV;
            rhoWP[ip] = 0.0;
            T[ip] = calcNeumann(T[GET3DINDEX_XYZ_Zm1],
                                T[GET3DINDEX_XYZ_Zm2],
                                T[GET3DINDEX_XYZ_Zm3],
                                T[GET3DINDEX_XYZ_Zm4],
                                T[GET3DINDEX_XYZ_Zm5],
                                T[GET3DINDEX_XYZ_Zm6]);
	    p[ip] = ig->solvep_idealgas(rhoP[ip], T[ip]);
	    rhoEP[ip] = ig->solverhoE(rhoP[ip], p[ip], U[ip], V[ip], W[ip]); //Not 100% about this?
        }END_FORZ1
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

    omp_set_nested(1);

    const int halfThreadCount = omp_get_num_threads()/4;

    #pragma omp parallel sections num_threads(halfThreadCount) 
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


    //auto t1 = std::chrono::system_clock::now();
   //derivX->calc1stDerivField(temp4, engyEulerX);
    //auto t2 = std::chrono::system_clock::now();
    //cout << "part 1 deriv test, sectioned: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()/(double)1000000000 << endl;    



    omp_set_nested(1);
    #pragma omp parallel sections num_threads(halfThreadCount) 
    {
        //Do the transposes
	#pragma omp section
        transposeXYZtoYZX_Fast(rhoP,  Nx, Ny, Nz, transRho, blocksize);

	#pragma omp section
        transposeXYZtoYZX_Fast(rhoUP, Nx, Ny, Nz, transRhoU, blocksize);

	#pragma omp section
        transposeXYZtoYZX_Fast(rhoVP, Nx, Ny, Nz, transRhoV, blocksize);

	#pragma omp section
        transposeXYZtoYZX_Fast(rhoWP, Nx, Ny, Nz, transRhoW, blocksize);

	#pragma omp section
        transposeXYZtoYZX_Fast(rhoEP, Nx, Ny, Nz, transRhoE, blocksize);

	#pragma omp section
        transposeXYZtoYZX_Fast(Ux,    Nx, Ny, Nz, transUx, blocksize);

	#pragma omp section
        transposeXYZtoYZX_Fast(Vx,    Nx, Ny, Nz, transVx, blocksize);

	#pragma omp section
        transposeXYZtoYZX_Fast(Wx,    Nx, Ny, Nz, transWx, blocksize);
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

    #pragma omp parallel for
    FOR_XYZ p[ip] = ig->solvep(transRho[ip], transRhoE[ip], U[ip], V[ip], W[ip]);
    #pragma omp parallel for
    FOR_XYZ T[ip] = ig->solveT(transRho[ip], p[ip]);

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
    #pragma omp parallel sections num_threads(halfThreadCount)
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
    #pragma omp parallel sections num_threads(halfThreadCount)
    {
        //Get the original conserved data from XYZ to ZXY
	#pragma omp section
        transposeXYZtoZXY_Fast(rhoP,  Nx, Ny, Nz, transRho, blocksize);
	#pragma omp section
        transposeXYZtoZXY_Fast(rhoUP, Nx, Ny, Nz, transRhoU, blocksize);
	#pragma omp section
        transposeXYZtoZXY_Fast(rhoVP, Nx, Ny, Nz, transRhoV, blocksize);
	#pragma omp section
        transposeXYZtoZXY_Fast(rhoWP, Nx, Ny, Nz, transRhoW, blocksize);
	#pragma omp section
        transposeXYZtoZXY_Fast(rhoEP, Nx, Ny, Nz, transRhoE, blocksize);

        //Move the Y Derivative data from YZX to ZXY
	#pragma omp section
        transposeYZXtoZXY_Fast(Uy, Nx, Ny, Nz, transUy, blocksize);
	#pragma omp section
        transposeYZXtoZXY_Fast(Vy, Nx, Ny, Nz, transVy, blocksize);
	#pragma omp section
        transposeYZXtoZXY_Fast(Wy, Nx, Ny, Nz, transWy, blocksize);

        //Move the X Derivative data from XYZ to ZXY
	#pragma omp section
        transposeXYZtoZXY_Fast(Ux, Nx, Ny, Nz, transUx, blocksize);
	#pragma omp section
        transposeXYZtoZXY_Fast(Vx, Nx, Ny, Nz, transVx, blocksize);
	#pragma omp section
        transposeXYZtoZXY_Fast(Wx, Nx, Ny, Nz, transWx, blocksize);
    }

    //Moving Data from YZX to XYZ

    #pragma omp parallel sections num_threads(halfThreadCount)
    {

	#pragma omp section
	{
            memcpy(temp, Uy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(temp, Nx, Ny, Nz, Uy, blocksize);
	}

	#pragma omp section
	{
            memcpy(temp2, Uyy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(temp2, Nx, Ny, Nz, Uyy, blocksize);
	}

	#pragma omp section
	{
            memcpy(temp3, Vy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(temp3, Nx, Ny, Nz, Vy, blocksize);
	}

	#pragma omp section
	{
            memcpy(temp4, Vyy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(temp4, Nx, Ny, Nz, Vyy, blocksize);
	}
	
	#pragma omp section
	{
            memcpy(U, Wy, sizeof(double)*Nx*Ny*Nz); //Starting here getting cute with reusing memory
            transposeYZXtoXYZ_Fast(U, Nx, Ny, Nz, Wy, blocksize);
	}


	#pragma omp section
	{
            memcpy(V, Wyy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(V, Nx, Ny, Nz, Wyy, blocksize);
	}


	#pragma omp section
	{
            memcpy(W, Ty, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(W, Nx, Ny, Nz, Ty, blocksize);
	}

	#pragma omp section
	{
            memcpy(T, Tyy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(T, Nx, Ny, Nz, Tyy, blocksize);
	}
    } 

    #pragma omp parallel sections num_threads(halfThreadCount)
    {

	#pragma omp section
	{
            memcpy(temp, Uxy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(temp, Nx, Ny, Nz, Uxy, blocksize);
	}

	#pragma omp section
	{
            memcpy(temp2, Vxy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(temp2, Nx, Ny, Nz, Vxy, blocksize);
	}

	#pragma omp section
	{
            memcpy(temp3, Wxy, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(temp3, Nx, Ny, Nz, Wxy, blocksize);
	}

	#pragma omp section
	{
            memcpy(temp4, contEulerY, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(temp4, Nx, Ny, Nz, contEulerY, blocksize);
	}

	#pragma omp section
	{
            memcpy(U, momXEulerY, sizeof(double)*Nx*Ny*Nz); //Getting cute with memory here again
            transposeYZXtoXYZ_Fast(U, Nx, Ny, Nz, momXEulerY, blocksize);
	}

	#pragma omp section
	{
            memcpy(V, momYEulerY, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(V, Nx, Ny, Nz, momYEulerY, blocksize);
	}

	#pragma omp section
	{
            memcpy(W, momZEulerY, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(W, Nx, Ny, Nz, momZEulerY, blocksize);
	}

	#pragma omp section
	{
            memcpy(T, engyEulerY, sizeof(double)*Nx*Ny*Nz);
            transposeYZXtoXYZ_Fast(T, Nx, Ny, Nz, engyEulerY, blocksize);
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

    #pragma omp parallel for
    FOR_XYZ p[ip] = ig->solvep(transRho[ip], transRhoE[ip], U[ip], V[ip], W[ip]);
    #pragma omp parallel for
    FOR_XYZ T[ip] = ig->solveT(transRho[ip], p[ip]);

    //Calculate the stuff for the Euler Derivatives
    #pragma omp parallel for
    FOR_XYZ temp[ip]  = transRhoU[ip]*W[ip];
    #pragma omp parallel for
    FOR_XYZ temp2[ip] = transRhoV[ip]*W[ip];
    #pragma omp parallel for
    FOR_XYZ temp3[ip] = transRhoW[ip]*W[ip] + p[ip];
    #pragma omp parallel for
    FOR_XYZ temp4[ip] = transRhoE[ip]*W[ip] + W[ip]*p[ip];

    omp_set_nested(1);
    #pragma omp parallel sections num_threads(halfThreadCount)
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
    #pragma omp parallel sections num_threads(halfThreadCount)
    {

		//Moving all the data back to XYZ
		#pragma omp section
		{
			memcpy(temp, Uz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(temp, Nx, Ny, Nz, Uz, blocksize);
		}

		#pragma omp section
		{
			memcpy(temp2, Uzz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(temp2, Nx, Ny, Nz, Uzz, blocksize);
		}

		#pragma omp section
		{
			memcpy(temp3, Vz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(temp3, Nx, Ny, Nz, Vz, blocksize);
		}

		#pragma omp section
		{
			memcpy(temp4, Vzz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(temp4, Nx, Ny, Nz, Vzz, blocksize);
		}

		#pragma omp section
		{
			memcpy(transUx, Wz, sizeof(double)*Nx*Ny*Nz); //Once again getting cute with memory starting here
			transposeZXYtoXYZ_Fast(transUx, Nx, Ny, Nz, Wz, blocksize);
		}

		#pragma omp section
		{
			memcpy(transVx, Wzz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(transVx, Nx, Ny, Nz, Wzz, blocksize);
		}

		#pragma omp section
		{
			memcpy(transWx, Tz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(transWx, Nx, Ny, Nz, Tz, blocksize);
		}

		#pragma omp section
		{
			memcpy(transUy, Tzz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(transUy, Nx, Ny, Nz, Tzz, blocksize);
		}

		#pragma omp section
		{
			memcpy(transVy, Uxz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(transVy, Nx, Ny, Nz, Uxz, blocksize);
		}

		#pragma omp section
		{
			memcpy(transWy, Vxz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(transWy, Nx, Ny, Nz, Vxz, blocksize);
		}

		#pragma omp section
		{
			memcpy(transRhoU, Wxz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(transRhoU, Nx, Ny, Nz, Wxz, blocksize);
		}

		#pragma omp section
		{
			memcpy(transRhoV, Uyz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(transRhoV, Nx, Ny, Nz, Uyz, blocksize);
		}

		#pragma omp section
		{
			memcpy(transRhoW, Vyz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(transRhoW, Nx, Ny, Nz, Vyz, blocksize);
		}

		#pragma omp section
		{
			memcpy(transRhoE, Wyz, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(transRhoE, Nx, Ny, Nz, Wyz, blocksize);
		}

		#pragma omp section
		{
			memcpy(U, contEulerZ, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(U, Nx, Ny, Nz, contEulerZ, blocksize);
		}

		#pragma omp section
		{
			memcpy(V, momXEulerZ, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(V, Nx, Ny, Nz, momXEulerZ, blocksize);
		}

		#pragma omp section
		{
			memcpy(W, momYEulerZ, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(W, Nx, Ny, Nz, momYEulerZ, blocksize);
		}

		#pragma omp section
		{
			memcpy(T, momZEulerZ, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(T, Nx, Ny, Nz, momZEulerZ, blocksize);
		}
	
		#pragma omp section
		{
			memcpy(p, engyEulerZ, sizeof(double)*Nx*Ny*Nz);
			transposeZXYtoXYZ_Fast(p, Nx, Ny, Nz, engyEulerZ, blocksize);
		}
    }

    //Going back to original...
    #pragma omp parallel for
    FOR_XYZ U[ip] = rhoUP[ip]/rhoP[ip];
    #pragma omp parallel for
    FOR_XYZ V[ip] = rhoVP[ip]/rhoP[ip];
    #pragma omp parallel for
    FOR_XYZ W[ip] = rhoWP[ip]/rhoP[ip];

    #pragma omp parallel for
    FOR_XYZ p[ip] = ig->solvep(rhoP[ip], rhoEP[ip], U[ip], V[ip], W[ip]);
    #pragma omp parallel for
    FOR_XYZ T[ip] = ig->solveT(rhoP[ip], p[ip]);

}


void CSolver::solveContinuity(){
   

    #pragma omp parallel
    {
        double *rhoP;
  	if(rkStep == 1){
            rhoP = rho1;
        }else if(rkStep == 2 || rkStep == 3 || rkStep == 4){
            rhoP = rhok;
        }
	
	double spgSource; 

	#pragma omp for
        FOR_XYZ{

	    if(spongeFlag)
		spgSource = calcSpongeSource(rhoP[ip], spg->spongeRhoAvg[ip], spg->sigma[ip]);
	    else
		spgSource = 0.0;
		
	    rhok2[ip]  = ts->dt*(-contEulerX[ip] - contEulerY[ip] - contEulerZ[ip] + spgSource);
	}
    }

}

void CSolver::solveXMomentum(){

	#pragma omp parallel
	{
            double *rhoUP;
            if(rkStep == 1){
                rhoUP = rhoU1;
            }else if(rkStep == 2 || rkStep == 3 || rkStep == 4){
                rhoUP = rhoUk;
            }

	    double MuX, MuY, MuZ, spgSource;
	    #pragma omp for
     	    FOR_XYZ{


	        if(spongeFlag)
	   	    spgSource = calcSpongeSource(rhoUP[ip], spg->spongeRhoUAvg[ip], spg->sigma[ip]);
	        else
		    spgSource = 0.0;

		//Calculate the viscosity derivatives
	        MuX = Amu[ip]*Tx[ip];
	        MuY = Amu[ip]*Ty[ip];
	        MuZ = Amu[ip]*Tz[ip];

		//Viscous Terms
                rhoUk2[ip]  = mu[ip]*((4.0/3.0)*Uxx[ip] + Uyy[ip] + Uzz[ip] + (1.0/3.0)*Vxy[ip] + (1.0/3.0)*Wxz[ip]);
	        rhoUk2[ip] += (4.0/3.0)*MuX*(Ux[ip] - 0.5*Vy[ip] - 0.5*Wz[ip]) + MuY*(Uy[ip] + Vx[ip]) + MuZ*(Wx[ip] + Uz[ip]);

		//Euler Terms
	        rhoUk2[ip] += -momXEulerX[ip] -momXEulerY[ip] -momXEulerZ[ip] + spgSource;
		rhoUk2[ip] *= ts->dt;
	    }

	}

}

void CSolver::solveYMomentum(){


    #pragma omp parallel
    {

        double *rhoVP;
        if(rkStep == 1){
            rhoVP = rhoV1;
        }else if(rkStep == 2 || rkStep == 3 || rkStep == 4){
            rhoVP = rhoVk;
        }

	double MuY, MuX, MuZ, spgSource;

        #pragma omp for
        FOR_XYZ{ 

            if(spongeFlag)
      	        spgSource = calcSpongeSource(rhoVP[ip], spg->spongeRhoVAvg[ip], spg->sigma[ip]);
	    else
	        spgSource = 0.0;

		
	    MuX = Amu[ip]*Tx[ip];
	    MuY = Amu[ip]*Ty[ip];
  	    MuZ = Amu[ip]*Tz[ip];

            //Viccous Terms 
	    rhoVk2[ip]  = mu[ip]*((4.0/3.0)*Vyy[ip] + Vxx[ip] + Vzz[ip] + (1.0/3.0)*Uxy[ip] + (1.0/3.0)*Wyz[ip]);
	    rhoVk2[ip] += (4.0/3.0)*MuY*(Vy[ip] - 0.5*Ux[ip] - 0.5*Wz[ip]) + MuX*(Uy[ip] + Vx[ip]) + MuZ*(Wy[ip] + Vz[ip]);
	    //Euler Terms
	    rhoVk2[ip] += -momYEulerX[ip] -momYEulerY[ip] -momYEulerZ[ip] + spgSource;
	
	    rhoVk2[ip] *= ts->dt;
        }

    }


}

void CSolver::solveZMomentum(){

    //Viscous Terms
    #pragma omp parallel for
    FOR_XYZ{
        rhoWk2[ip]  = mu[ip]*((4.0/3.0)*Wzz[ip] + Wyy[ip] + Wxx[ip] + (1.0/3.0)*Uxz[ip] + (1.0/3.0)*Vyz[ip]);
    }

    #pragma omp parallel for
    FOR_XYZ{
	double MuZ = Amu[ip]*Tz[ip];
	rhoWk2[ip] += (4.0/3.0)*MuZ*(Wz[ip] - 0.5*Ux[ip] - 0.5*Vy[ip]);
    }

    #pragma omp parallel for
    FOR_XYZ{
	double MuX = Amu[ip]*Tx[ip];
	rhoWk2[ip] += MuX*(Wx[ip] + Uz[ip]); 
    }

    #pragma omp parallel for
    FOR_XYZ{
	double MuY = Amu[ip]*Ty[ip];
	rhoWk2[ip] += MuY*(Wy[ip] + Vz[ip]); 
    }

    #pragma omp parallel for
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

    #pragma omp parallel for
        FOR_XYZ{
            rhoWk2[ip]  += calcSpongeSource(rhoWP[ip], spg->spongeRhoWAvg[ip], spg->sigma[ip]);
        }
    }

    #pragma omp parallel for
    FOR_XYZ rhoWk2[ip] *= ts->dt;

}


void CSolver::solveEnergy(){

    double *qtemp = new double[N];
    double *vtemp1 = new double[N];
    double *vtemp2 = new double[N];
    double *engyEuler = new double[N];

    //Heat Transfer Terms
    #pragma omp parallel for
    FOR_XYZ{
	double MuX = Amu[ip]*Tx[ip];
	qtemp[ip]  = MuX*Tx[ip];
    }

    #pragma omp parallel for 
    FOR_XYZ{
	double MuY = Amu[ip]*Ty[ip];
	qtemp[ip] += MuY*Ty[ip];
    }  

    #pragma omp parallel for
    FOR_XYZ{
	double MuZ = Amu[ip]*Tz[ip];
	qtemp[ip] += MuZ*Tz[ip];
    }

    #pragma omp parallel for
    FOR_XYZ{ 
	//Heat Transfer Terms
	qtemp[ip] += mu[ip]*Txx[ip] + mu[ip]*Tyy[ip] + mu[ip]*Tzz[ip];
    }

    #pragma omp parallel for
    FOR_XYZ qtemp[ip] *= ig->cp/ig->Pr;

    //Viscous Energy terms w/o viscosity derivatives...
    #pragma omp parallel for
    FOR_XYZ vtemp1[ip]  = U[ip]*((4.0/3.0)*Uxx[ip] + Uyy[ip] + Uzz[ip]);
    #pragma omp parallel for
    FOR_XYZ vtemp1[ip] += V[ip]*(Vxx[ip] + (4.0/3.0)*Vyy[ip] + Vzz[ip]);
    #pragma omp parallel for
    FOR_XYZ vtemp1[ip] += W[ip]*(Wxx[ip] + Wyy[ip] + (4.0/3.0)*Wzz[ip]);
	
    #pragma omp parallel for
    FOR_XYZ vtemp1[ip] += (4.0/3.0)*(Ux[ip]*Ux[ip] + Vy[ip]*Vy[ip] + Wz[ip]*Wz[ip]);

    #pragma omp parallel for
    FOR_XYZ vtemp1[ip] += Uy[ip]*Uy[ip] + Uz[ip]*Uz[ip];
    #pragma omp parallel for
    FOR_XYZ vtemp1[ip] += Vx[ip]*Vx[ip] + Vz[ip]*Vz[ip];
    #pragma omp parallel for
    FOR_XYZ vtemp1[ip] += Wx[ip]*Wx[ip] + Wy[ip]*Wy[ip];

    #pragma omp parallel for
    FOR_XYZ vtemp1[ip] += -(4.0/3.0)*(Ux[ip]*Vy[ip] + Ux[ip]*Wz[ip] + Vy[ip]*Wz[ip]);
	
    #pragma omp parallel for
    FOR_XYZ vtemp1[ip] += 2.0*(Uy[ip]*Vx[ip] + Uz[ip]*Wx[ip] + Vz[ip]*Wy[ip]);

    #pragma omp parallel for
    FOR_XYZ vtemp1[ip] += (1.0/3.0)*(U[ip]*Vxy[ip] + U[ip]*Wxz[ip] + V[ip]*Uxy[ip]); 
    #pragma omp parallel for
    FOR_XYZ vtemp1[ip] += (1.0/3.0)*(V[ip]*Wyz[ip] + W[ip]*Uxz[ip] + W[ip]*Vyz[ip]); 
	
    #pragma omp parallel for
    FOR_XYZ vtemp1[ip] *= mu[ip];

    //Viscous Energy terms w/ viscosity derivatives...
    #pragma omp parallel for
    FOR_XYZ{ 
	double MuX = Amu[ip]*Tx[ip];
	double MuY = Amu[ip]*Ty[ip];
	double MuZ = Amu[ip]*Tz[ip];
	vtemp2[ip]   = (4.0/3.0)*(U[ip]*MuX*Ux[ip] + V[ip]*MuY*Vy[ip] + W[ip]*MuZ*Wz[ip]);
    }

    #pragma omp parallel for
    FOR_XYZ{
	double MuX = Amu[ip]*Tx[ip];
	vtemp2[ip]  += -(2.0/3.0)*U[ip]*MuX*(Vy[ip] + Wz[ip]);
    }

    #pragma omp parallel for
    FOR_XYZ{
	double MuY = Amu[ip]*Ty[ip];
	vtemp2[ip]  += -(2.0/3.0)*V[ip]*MuY*(Ux[ip] + Wz[ip]);
    }

    #pragma omp parallel for
    FOR_XYZ{
	double MuZ = Amu[ip]*Tz[ip];
	vtemp2[ip]  += -(2.0/3.0)*W[ip]*MuZ*(Ux[ip] + Vy[ip]);
    }

    #pragma omp parallel for
    FOR_XYZ{
	double MuY = Amu[ip]*Ty[ip];
	vtemp2[ip]  += U[ip]*MuY*(Uy[ip] + Vx[ip]);
    }

    #pragma omp parallel for
    FOR_XYZ{
	double MuZ = Amu[ip]*Tz[ip];
        vtemp2[ip]  += U[ip]*MuZ*(Uz[ip] + Wx[ip]);
    }

    #pragma omp parallel for
    FOR_XYZ{
	double MuX = Amu[ip]*Tx[ip];
	vtemp2[ip]  += V[ip]*MuX*(Uy[ip] + Vx[ip]);
    }

    #pragma omp parallel for
    FOR_XYZ{
	double MuZ = Amu[ip]*Tz[ip];
	vtemp2[ip]  += V[ip]*MuZ*(Vz[ip] + Wy[ip]);
    }

    #pragma omp parallel for
    FOR_XYZ{
	double MuX = Amu[ip]*Tx[ip];
	vtemp2[ip]  += W[ip]*MuX*(Uz[ip] + Wx[ip]);
    }

    #pragma omp parallel for
    FOR_XYZ{
	double MuY = Amu[ip]*Ty[ip];
	vtemp2[ip]  += W[ip]*MuY*(Vz[ip] + Wy[ip]);
    }

    //Euler terms
    #pragma omp parallel for
    FOR_XYZ engyEuler[ip]  = -engyEulerX[ip] - engyEulerY[ip] - engyEulerZ[ip];

    #pragma omp parallel for
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

    ////////////////////////////////
    //ADIABATIC AND MOVING WALL BC// 
    ////////////////////////////////

    if(bc->bcX0 == BC::ADIABATIC_WALL || bc->bcX0 == BC::MOVING_ADIABATIC_WALL){
	#pragma omp parallel for
	FOR_X0{
	    rhok2[ip]  = -ts->dt*contEulerX[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    //rhoEk2[ip] = -ts->dt*(engyEulerX[ip] - (mu[ip]/ig->Pr/(ig->gamma-1.0))*Txx[ip]);
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORX0
    }

    if(bc->bcX1 == BC::ADIABATIC_WALL  || bc->bcX1 == BC::MOVING_ADIABATIC_WALL){
	#pragma omp parallel for
	FOR_X1{
	    rhok2[ip]  = -ts->dt*contEulerX[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    //rhoEk2[ip] = -ts->dt*(engyEulerX[ip] - (mu[ip]/ig->Pr/(ig->gamma-1.0))*Txx[ip]);
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORX1
    }   

    if(bc->bcY0 == BC::ADIABATIC_WALL || bc->bcY0 == BC::MOVING_ADIABATIC_WALL){
	#pragma omp parallel for
	FOR_Y0{
	    rhok2[ip]  = -ts->dt*contEulerY[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    //rhoEk2[ip] = -ts->dt*(engyEulerY[ip] - (mu[ip]/ig->Pr/(ig->gamma-1.0))*Tyy[ip]);
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORY0
    }

    if(bc->bcY1 == BC::ADIABATIC_WALL || bc->bcY1 == BC::MOVING_ADIABATIC_WALL){
	#pragma omp parallel for
	FOR_Y1{
	    rhok2[ip]  = -ts->dt*contEulerY[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    //rhoEk2[ip] = -ts->dt*(engyEulerY[ip] - (mu[ip]/ig->Pr/(ig->gamma-1.0))*Tyy[ip]);
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORY1
    }

    if(bc->bcZ0 == BC::ADIABATIC_WALL || bc->bcZ0 == BC::MOVING_ADIABATIC_WALL){
	#pragma omp parallel for
	FOR_Z0{
	    rhok2[ip]  = -ts->dt*contEulerZ[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    //rhoEk2[ip] = -ts->dt*(engyEulerZ[ip] - (mu[ip]/ig->Pr/(ig->gamma-1.0))*Tzz[ip]);
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORZ0
    }

    if(bc->bcZ1 == BC::ADIABATIC_WALL || bc->bcZ1 == BC::MOVING_ADIABATIC_WALL){
	#pragma omp parallel for
	FOR_Z1{
	    rhok2[ip]  = -ts->dt*contEulerZ[ip];
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    //rhoEk2[ip] = -ts->dt*(engyEulerZ[ip] - (mu[ip]/ig->Pr/(ig->gamma-1.0))*Tzz[ip]);
	    rhoEk2[ip] = 0.0;//Not 100% 
	}END_FORZ1
    }


    /////////////
    //SPONGE BC//
    /////////////

    if(bc->bcX0 == BC::SPONGE){
	#pragma omp parallel for
	FOR_X0{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORX0
    }

    if(bc->bcX1 == BC::SPONGE){
	#pragma omp parallel for
	FOR_X1{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORX1
    }   

    if(bc->bcY0 == BC::SPONGE){
	#pragma omp parallel for
	FOR_Y0{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORY0
    }

    if(bc->bcY1 == BC::SPONGE){
	#pragma omp parallel for
	FOR_Y1{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORY1
    }

    if(bc->bcZ0 == BC::SPONGE){
	#pragma omp parallel for
	FOR_Z0{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORZ0
    }

    if(bc->bcZ1 == BC::SPONGE){
	#pragma omp parallel for
	FOR_Z1{
	    rhok2[ip]  = 0.0;
	    rhoUk2[ip] = 0.0;
	    rhoVk2[ip] = 0.0;
	    rhoWk2[ip] = 0.0;
	    rhoEk2[ip] = 0.0;
	}END_FORZ1
    }


}

void CSolver::updateConservedData(){

    if(rkStep == 1){

	//Add to final solution
        #pragma omp parallel for
	FOR_XYZ rho2[ip]  = rho1[ip]  + rhok2[ip]/6.0;
        #pragma omp parallel for
	FOR_XYZ rhoU2[ip] = rhoU1[ip] + rhoUk2[ip]/6.0;
        #pragma omp parallel for
	FOR_XYZ rhoV2[ip] = rhoV1[ip] + rhoVk2[ip]/6.0;
        #pragma omp parallel for
	FOR_XYZ rhoW2[ip] = rhoW1[ip] + rhoWk2[ip]/6.0;
        #pragma omp parallel for
	FOR_XYZ rhoE2[ip] = rhoE1[ip] + rhoEk2[ip]/6.0;

	//Calculate intermediate solution
        #pragma omp parallel for
	FOR_XYZ rhok[ip]  = rho1[ip]  + rhok2[ip]/2.0; 
        #pragma omp parallel for
	FOR_XYZ rhoUk[ip] = rhoU1[ip] + rhoUk2[ip]/2.0; 
        #pragma omp parallel for
	FOR_XYZ rhoVk[ip] = rhoV1[ip] + rhoVk2[ip]/2.0; 
        #pragma omp parallel for
	FOR_XYZ rhoWk[ip] = rhoW1[ip] + rhoWk2[ip]/2.0; 
        #pragma omp parallel for
	FOR_XYZ rhoEk[ip] = rhoE1[ip] + rhoEk2[ip]/2.0; 

    }else if(rkStep == 2){

	//Add to final solution
        #pragma omp parallel for
	FOR_XYZ rho2[ip]  += rhok2[ip]/3.0;
        #pragma omp parallel for
	FOR_XYZ rhoU2[ip] += rhoUk2[ip]/3.0;
        #pragma omp parallel for
	FOR_XYZ rhoV2[ip] += rhoVk2[ip]/3.0;
        #pragma omp parallel for
	FOR_XYZ rhoW2[ip] += rhoWk2[ip]/3.0;
        #pragma omp parallel for
	FOR_XYZ rhoE2[ip] += rhoEk2[ip]/3.0;

	//Calculate intermediate solution
        #pragma omp parallel for
	FOR_XYZ rhok[ip]  = rho1[ip]  + rhok2[ip]/2.0; 
        #pragma omp parallel for
	FOR_XYZ rhoUk[ip] = rhoU1[ip] + rhoUk2[ip]/2.0; 
        #pragma omp parallel for
	FOR_XYZ rhoVk[ip] = rhoV1[ip] + rhoVk2[ip]/2.0; 
        #pragma omp parallel for
	FOR_XYZ rhoWk[ip] = rhoW1[ip] + rhoWk2[ip]/2.0; 
        #pragma omp parallel for
	FOR_XYZ rhoEk[ip] = rhoE1[ip] + rhoEk2[ip]/2.0; 

    }else if(rkStep == 3){

	//Add to final solution
        #pragma omp parallel for
	FOR_XYZ rho2[ip]  += rhok2[ip]/3.0;
        #pragma omp parallel for
	FOR_XYZ rhoU2[ip] += rhoUk2[ip]/3.0;
        #pragma omp parallel for
	FOR_XYZ rhoV2[ip] += rhoVk2[ip]/3.0;
        #pragma omp parallel for
	FOR_XYZ rhoW2[ip] += rhoWk2[ip]/3.0;
        #pragma omp parallel for
	FOR_XYZ rhoE2[ip] += rhoEk2[ip]/3.0;

	//Calculate intermediate solution
        #pragma omp parallel for
	FOR_XYZ rhok[ip]  = rho1[ip]  + rhok2[ip]; 
        #pragma omp parallel for
	FOR_XYZ rhoUk[ip] = rhoU1[ip] + rhoUk2[ip]; 
        #pragma omp parallel for
	FOR_XYZ rhoVk[ip] = rhoV1[ip] + rhoVk2[ip]; 
        #pragma omp parallel for
	FOR_XYZ rhoWk[ip] = rhoW1[ip] + rhoWk2[ip]; 
        #pragma omp parallel for
	FOR_XYZ rhoEk[ip] = rhoE1[ip] + rhoEk2[ip]; 

    }else if(rkStep == 4){

	//Add to final solution
        #pragma omp parallel for
	FOR_XYZ rho2[ip]  += rhok2[ip]/6.0;
        #pragma omp parallel for
	FOR_XYZ rhoU2[ip] += rhoUk2[ip]/6.0;
        #pragma omp parallel for
	FOR_XYZ rhoV2[ip] += rhoVk2[ip]/6.0;
        #pragma omp parallel for
	FOR_XYZ rhoW2[ip] += rhoWk2[ip]/6.0;
        #pragma omp parallel for
	FOR_XYZ rhoE2[ip] += rhoEk2[ip]/6.0;

    }

}

void CSolver::filterConservedData(){

    const int halfThreadCount = omp_get_num_threads()/4;
    omp_set_nested(1);


    //Need to do round robin of filtering directions...
    if(timeStep%ts->filterStep == 0){

        //Advance the filtering time step       
        filterTimeStep++;

        //Going to try and be cute to minimize dmemory allocation
        if(filterTimeStep%3 == 1){

            //Here we'll do X->Y->Z     

    	    #pragma omp parallel sections num_threads(halfThreadCount) 
 	    {
		#pragma omp section
                filtX->filterField(rho2,  rho1);
		#pragma omp section
                filtX->filterField(rhoU2, temp);
		#pragma omp section
                filtX->filterField(rhoV2, temp2);
		#pragma omp section
                filtX->filterField(rhoW2, temp3);
		#pragma omp section
                filtX->filterField(rhoE2, temp4);
	    }

            //Do the transpose to YZX space
    	    #pragma omp parallel sections num_threads(halfThreadCount) 
 	    {
		#pragma omp section
                transposeXYZtoYZX_Fast(rho1,   Nx, Ny, Nz, rho2, blocksize);
		#pragma omp section
                transposeXYZtoYZX_Fast(temp,   Nx, Ny, Nz, rhoU2, blocksize);
		#pragma omp section
                transposeXYZtoYZX_Fast(temp2,  Nx, Ny, Nz, rhoV2, blocksize);
		#pragma omp section
                transposeXYZtoYZX_Fast(temp3,  Nx, Ny, Nz, rhoW2, blocksize);
		#pragma omp section
                transposeXYZtoYZX_Fast(temp4,  Nx, Ny, Nz, rhoE2, blocksize);
	    }

            //filter in the Y direction
     	    #pragma omp parallel sections num_threads(halfThreadCount) 
 	    {
		#pragma omp section
  	        filtY->filterField(rho2,  rho1);
		#pragma omp section
                filtY->filterField(rhoU2, temp);
		#pragma omp section
                filtY->filterField(rhoV2, temp2);
		#pragma omp section
                filtY->filterField(rhoW2, temp3);
		#pragma omp section
                filtY->filterField(rhoE2, temp4);
	    }

            //tranpose from YZX to ZXY
     	    #pragma omp parallel sections num_threads(halfThreadCount) 
 	    {
		#pragma omp section
                transposeYZXtoZXY_Fast(rho1,   Nx, Ny, Nz, rho2, blocksize);
		#pragma omp section
                transposeYZXtoZXY_Fast(temp,   Nx, Ny, Nz, rhoU2, blocksize);
		#pragma omp section
                transposeYZXtoZXY_Fast(temp2,  Nx, Ny, Nz, rhoV2, blocksize);
		#pragma omp section
                transposeYZXtoZXY_Fast(temp3,  Nx, Ny, Nz, rhoW2, blocksize);
		#pragma omp section
                transposeYZXtoZXY_Fast(temp4,  Nx, Ny, Nz, rhoE2, blocksize);
	    }

            //filter in the Y direction
     	    #pragma omp parallel sections num_threads(halfThreadCount) 
 	    {
		#pragma omp section
                filtZ->filterField(rho2,  rho1);
		#pragma omp section
                filtZ->filterField(rhoU2, temp);
		#pragma omp section
                filtZ->filterField(rhoV2, temp2);
		#pragma omp section
                filtZ->filterField(rhoW2, temp3);
		#pragma omp section
                filtZ->filterField(rhoE2, temp4);
	    }

            //get us back to XYZ from ZXY
     	    #pragma omp parallel sections num_threads(halfThreadCount) 
 	    {
		#pragma omp section
                transposeZXYtoXYZ_Fast(rho1,   Nx, Ny, Nz, rho2, blocksize);
		#pragma omp section
                transposeZXYtoXYZ_Fast(temp,   Nx, Ny, Nz, rhoU1, blocksize);
		#pragma omp section
                transposeZXYtoXYZ_Fast(temp2,  Nx, Ny, Nz, rhoV1, blocksize);
		#pragma omp section
                transposeZXYtoXYZ_Fast(temp3,  Nx, Ny, Nz, rhoW1, blocksize);
		#pragma omp section
                transposeZXYtoXYZ_Fast(temp4,  Nx, Ny, Nz, rhoE1, blocksize);
	    }
            //from being cute with memory allocation need to copy rho2 to rho1
            memcpy(rho1, rho2, sizeof(double)*Nx*Ny*Nz);

        }else if(filterTimeStep%3 == 2){

            //Here we'll do Y->Z->X     

            //Do the transpose to YZX space first
      	    #pragma omp parallel sections num_threads(halfThreadCount) 
 	    { 
		#pragma omp section
                transposeXYZtoYZX_Fast(rho2,   Nx, Ny, Nz, rho1, blocksize);
		#pragma omp section
                transposeXYZtoYZX_Fast(rhoU2,  Nx, Ny, Nz, temp, blocksize);
		#pragma omp section
                transposeXYZtoYZX_Fast(rhoV2,  Nx, Ny, Nz, temp2, blocksize);
		#pragma omp section
                transposeXYZtoYZX_Fast(rhoW2,  Nx, Ny, Nz, temp3, blocksize);
		#pragma omp section
                transposeXYZtoYZX_Fast(rhoE2,  Nx, Ny, Nz, temp4, blocksize);
	    }	

            //filter in the Y direction
       	    #pragma omp parallel sections num_threads(halfThreadCount) 
 	    { 
		#pragma omp section
                filtY->filterField(rho1,  rho2);
		#pragma omp section
                filtY->filterField(temp,  rhoU2);
		#pragma omp section
                filtY->filterField(temp2, rhoV2);
		#pragma omp section
                filtY->filterField(temp3, rhoW2);
		#pragma omp section
                filtY->filterField(temp4, rhoE2);
	    }

            //Move to ZXY space next
       	    #pragma omp parallel sections num_threads(halfThreadCount) 
 	    { 
		#pragma omp section
                transposeYZXtoZXY_Fast(rho2,   Nx, Ny, Nz, rho1, blocksize);
		#pragma omp section
                transposeYZXtoZXY_Fast(rhoU2,  Nx, Ny, Nz, temp, blocksize);
		#pragma omp section
                transposeYZXtoZXY_Fast(rhoV2,  Nx, Ny, Nz, temp2, blocksize);
		#pragma omp section
                transposeYZXtoZXY_Fast(rhoW2,  Nx, Ny, Nz, temp3, blocksize);
		#pragma omp section
                transposeYZXtoZXY_Fast(rhoE2,  Nx, Ny, Nz, temp4, blocksize);
	    }

            //filter in the Z direction
       	    #pragma omp parallel sections num_threads(halfThreadCount) 
 	    { 
		#pragma omp section
                filtZ->filterField(rho1,  rho2);
		#pragma omp section
                filtZ->filterField(temp,  rhoU2);
		#pragma omp section
                filtZ->filterField(temp2, rhoV2);
		#pragma omp section
                filtZ->filterField(temp3, rhoW2);
		#pragma omp section
                filtZ->filterField(temp4, rhoE2);
	    }

            //transpose from ZXY to XYZ
       	    #pragma omp parallel sections num_threads(halfThreadCount) 
 	    { 
		#pragma omp section
                transposeZXYtoXYZ_Fast(rho2,   Nx, Ny, Nz, rho1, blocksize);
		#pragma omp section
                transposeZXYtoXYZ_Fast(rhoU2,  Nx, Ny, Nz, temp, blocksize);
		#pragma omp section
                transposeZXYtoXYZ_Fast(rhoV2,  Nx, Ny, Nz, temp2, blocksize);
		#pragma omp section
                transposeZXYtoXYZ_Fast(rhoW2,  Nx, Ny, Nz, temp3, blocksize);
		#pragma omp section
                transposeZXYtoXYZ_Fast(rhoE2,  Nx, Ny, Nz, temp4, blocksize);
	    }

            //filter in the X direction
       	    #pragma omp parallel sections num_threads(halfThreadCount) 
  	    { 
		#pragma omp section
                filtX->filterField(rho1,  rho2);
		#pragma omp section
                filtX->filterField(temp,  rhoU1);
		#pragma omp section
                filtX->filterField(temp2, rhoV1);
		#pragma omp section
                filtX->filterField(temp3, rhoW1);
		#pragma omp section
                filtX->filterField(temp4, rhoE1);
	    }

            //from being cute with memory allocation need to copy rho2 to rho1
            memcpy(rho1, rho2, sizeof(double)*Nx*Ny*Nz);


        }else{

            //Here we'll do Z->X->Y     
            //Do the transpose to ZXY space first
       	    #pragma omp parallel sections num_threads(halfThreadCount) 
  	    { 
		#pragma omp section
                transposeXYZtoZXY_Fast(rho2,   Nx, Ny, Nz, rho1, blocksize);
		#pragma omp section
                transposeXYZtoZXY_Fast(rhoU2,  Nx, Ny, Nz, temp, blocksize);
		#pragma omp section
                transposeXYZtoZXY_Fast(rhoV2,  Nx, Ny, Nz, temp2, blocksize);
		#pragma omp section
                transposeXYZtoZXY_Fast(rhoW2,  Nx, Ny, Nz, temp3, blocksize);
		#pragma omp section
                transposeXYZtoZXY_Fast(rhoE2,  Nx, Ny, Nz, temp4, blocksize);
	    }

            //filter in the Z direction
            #pragma omp parallel sections num_threads(halfThreadCount) 
  	    {
		#pragma omp section
                filtZ->filterField(rho1,  rho2);
		#pragma omp section
                filtZ->filterField(temp,  rhoU2);
		#pragma omp section
                filtZ->filterField(temp2, rhoV2);
		#pragma omp section
                filtZ->filterField(temp3, rhoW2);
		#pragma omp section
                filtZ->filterField(temp4, rhoE2);
	    }

            //Move to XYZ space next
            #pragma omp parallel sections num_threads(halfThreadCount) 
  	    {
		#pragma omp section
                transposeZXYtoXYZ_Fast(rho2,   Nx, Ny, Nz, rho1, blocksize);
		#pragma omp section
                transposeZXYtoXYZ_Fast(rhoU2,  Nx, Ny, Nz, temp, blocksize);
		#pragma omp section
                transposeZXYtoXYZ_Fast(rhoV2,  Nx, Ny, Nz, temp2, blocksize);
		#pragma omp section
                transposeZXYtoXYZ_Fast(rhoW2,  Nx, Ny, Nz, temp3, blocksize);
		#pragma omp section
                transposeZXYtoXYZ_Fast(rhoE2,  Nx, Ny, Nz, temp4, blocksize);
	    }

            //filter in the X direction
            #pragma omp parallel sections num_threads(halfThreadCount) 
  	    {
		#pragma omp section
                filtX->filterField(rho1,  rho2);
		#pragma omp section
                filtX->filterField(temp,  rhoU2);
		#pragma omp section
                filtX->filterField(temp2, rhoV2);
		#pragma omp section
                filtX->filterField(temp3, rhoW2);
		#pragma omp section
                filtX->filterField(temp4, rhoE2);
	    }

            //transpose from XYZ to YZX
            #pragma omp parallel sections num_threads(halfThreadCount) 
  	    { 
		#pragma omp section
                transposeXYZtoYZX_Fast(rho2,   Nx, Ny, Nz, rho1, blocksize);
		#pragma omp section
                transposeXYZtoYZX_Fast(rhoU2,  Nx, Ny, Nz, temp, blocksize);
		#pragma omp section
                transposeXYZtoYZX_Fast(rhoV2,  Nx, Ny, Nz, temp2, blocksize);
		#pragma omp section
                transposeXYZtoYZX_Fast(rhoW2,  Nx, Ny, Nz, temp3, blocksize);
		#pragma omp section
                transposeXYZtoYZX_Fast(rhoE2,  Nx, Ny, Nz, temp4, blocksize);
	    }

            //filter in the Y direction
            #pragma omp parallel sections num_threads(halfThreadCount) 
  	    { 
		#pragma omp section
	        filtY->filterField(rho1,  rho2);
		#pragma omp section
                filtY->filterField(temp,  rhoU2);
		#pragma omp section
                filtY->filterField(temp2, rhoV2);
		#pragma omp section
                filtY->filterField(temp3, rhoW2);
		#pragma omp section
                filtY->filterField(temp4, rhoE2);
	    }

            //have an extra step here to go from YZX to XYZ
            #pragma omp parallel sections num_threads(halfThreadCount) 
  	    { 
		#pragma omp section
	        transposeYZXtoXYZ_Fast(rho2,   Nx, Ny, Nz, rho1, blocksize);
		#pragma omp section
                transposeYZXtoXYZ_Fast(rhoU2,  Nx, Ny, Nz, rhoU1, blocksize);
		#pragma omp section
                transposeYZXtoXYZ_Fast(rhoV2,  Nx, Ny, Nz, rhoV1, blocksize);
		#pragma omp section
                transposeYZXtoXYZ_Fast(rhoW2,  Nx, Ny, Nz, rhoW1, blocksize);
		#pragma omp section
                transposeYZXtoXYZ_Fast(rhoE2,  Nx, Ny, Nz, rhoE1, blocksize);
	    }

        }
 
    //If not filtering, need to copy the solution over to the *1 variables
    }else{
        #pragma omp parallel sections  
  	{

	    #pragma omp section
	    memcpy(rho1,  rho2, sizeof(double)*Nx*Ny*Nz);
	    #pragma omp section
	    memcpy(rhoU1, rhoU2, sizeof(double)*Nx*Ny*Nz);
	    #pragma omp section
	    memcpy(rhoV1, rhoV2, sizeof(double)*Nx*Ny*Nz);
	    #pragma omp section
	    memcpy(rhoW1, rhoW2, sizeof(double)*Nx*Ny*Nz);
	    #pragma omp section
	    memcpy(rhoE1, rhoE2, sizeof(double)*Nx*Ny*Nz);
	}

    }
};

void CSolver::updateNonConservedData(){
    if(rkStep == 1 || rkStep == 2 || rkStep == 3){

        #pragma omp parallel for
	FOR_XYZ U[ip]   = ig->solveU(rhok[ip], rhoUk[ip]);
        #pragma omp parallel for
	FOR_XYZ V[ip]   = ig->solveU(rhok[ip], rhoVk[ip]);
        #pragma omp parallel for
	FOR_XYZ W[ip]   = ig->solveU(rhok[ip], rhoWk[ip]);
        #pragma omp parallel for
	FOR_XYZ p[ip]   = ig->solvep(rhok[ip], rhoEk[ip], U[ip], V[ip], W[ip]);
        #pragma omp parallel for
	FOR_XYZ T[ip]   = ig->solveT(rhok[ip], p[ip]);
        #pragma omp parallel for
	FOR_XYZ mu[ip]  = ig->solveMu(T[ip]);
        #pragma omp parallel for
	FOR_XYZ Amu[ip] = ig->solveAmu(T[ip]);
        #pragma omp parallel for
	FOR_XYZ sos[ip] = ig->solveSOS(rhok[ip], p[ip]);

    }else if(rkStep == 4){

        #pragma omp parallel for
	FOR_XYZ U[ip]   = ig->solveU(rho1[ip], rhoU1[ip]);
        #pragma omp parallel for
	FOR_XYZ V[ip]   = ig->solveU(rho1[ip], rhoV1[ip]);
        #pragma omp parallel for
	FOR_XYZ W[ip]   = ig->solveU(rho1[ip], rhoW1[ip]);
        #pragma omp parallel for
	FOR_XYZ p[ip]   = ig->solvep(rho1[ip], rhoE1[ip], U[ip], V[ip], W[ip]);
        #pragma omp parallel for
	FOR_XYZ T[ip]   = ig->solveT(rho1[ip], p[ip]);
        #pragma omp parallel for
	FOR_XYZ mu[ip]  = ig->solveMu(T[ip]);
        #pragma omp parallel for
	FOR_XYZ Amu[ip] = ig->solveAmu(T[ip]);
        #pragma omp parallel for
	FOR_XYZ sos[ip] = ig->solveSOS(rho1[ip], p[ip]);

    }

}


void CSolver::updateSponge(){
    if(spongeFlag){
	double eps = 1.0/(spg->avgT/ts->dt + 1.0);
	#pragma omp parallel for
	FOR_XYZ spg->spongeRhoAvg[ip]  += eps*(rho1[ip]  - spg->spongeRhoAvg[ip]);	
	#pragma omp parallel for
	FOR_XYZ spg->spongeRhoUAvg[ip] += eps*(rhoU1[ip] - spg->spongeRhoUAvg[ip]);	
	#pragma omp parallel for
	FOR_XYZ spg->spongeRhoVAvg[ip] += eps*(rhoV1[ip] - spg->spongeRhoVAvg[ip]);	
	#pragma omp parallel for
	FOR_XYZ spg->spongeRhoWAvg[ip] += eps*(rhoW1[ip] - spg->spongeRhoWAvg[ip]);	
	#pragma omp parallel for
	FOR_XYZ spg->spongeRhoEAvg[ip] += eps*(rhoE1[ip] - spg->spongeRhoEAvg[ip]);	
	
	#pragma omp parallel for
	FOR_XYZ spg->spongeRhoEAvg[ip] = spg->epsP*spg->spongeRhoEAvg[ip] + (1.0 -  spg->epsP)*(spg->spongeP/(ig->gamma-1.0) \
					 + 0.5*(spg->spongeRhoUAvg[ip]*spg->spongeRhoUAvg[ip] + spg->spongeRhoVAvg[ip]*spg->spongeRhoVAvg[ip] \
					 + spg->spongeRhoWAvg[ip]*spg->spongeRhoWAvg[ip])/spg->spongeRhoAvg[ip]);

	
        if(bc->bcX0 == BC::SPONGE){
	    #pragma omp parallel for
	    FOR_X0{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORX0
        }

        if(bc->bcX1 == BC::SPONGE){
	    #pragma omp parallel for
	    FOR_X1{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORX1
        }   

        if(bc->bcY0 == BC::SPONGE){
	    #pragma omp parallel for
	    FOR_Y0{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORY0
        }

        if(bc->bcY1 == BC::SPONGE){
	    #pragma omp parallel for
	    FOR_Y1{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORY1
        }

        if(bc->bcZ0 == BC::SPONGE){
	    #pragma omp parallel for
	    FOR_Z0{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORZ0

        }

        if(bc->bcZ1 == BC::SPONGE){
	    #pragma omp parallel for
	    FOR_Z1{
		rho1[ip]  = spg->spongeRhoAvg[ip];
		rhoU1[ip] = spg->spongeRhoUAvg[ip];
		rhoV1[ip] = spg->spongeRhoVAvg[ip];
		rhoW1[ip] = spg->spongeRhoWAvg[ip];
		rhoE1[ip] = spg->spongeRhoEAvg[ip];
	    }END_FORZ1

        }

    }
};

void CSolver::calcTurbulenceQuantities(){


    FOR_XYZ{
	double Sij[3][3];
	Sij[0][0] = 0.5*(Ux[ip] + Ux[ip]) - (1.0/3.0)*(Ux[ip] + Vy[ip] + Wz[ip]);
	Sij[1][1] = 0.5*(Vy[ip] + Vy[ip]) - (1.0/3.0)*(Ux[ip] + Vy[ip] + Wz[ip]);
	Sij[2][2] = 0.5*(Wz[ip] + Wz[ip]) - (1.0/3.0)*(Ux[ip] + Vy[ip] + Wz[ip]);
	Sij[0][1] = 0.5*(Uy[ip] + Vx[ip]);
	Sij[1][0] = 0.5*(Uy[ip] + Vx[ip]);
	Sij[0][2] = 0.5*(Uz[ip] + Wx[ip]);
	Sij[2][0] = 0.5*(Uz[ip] + Wx[ip]);
	Sij[1][2] = 0.5*(Vz[ip] + Wy[ip]);
	Sij[2][1] = 0.5*(Vz[ip] + Wy[ip]);

	turbdiss[ip]  = Sij[0][0]*Ux[ip];
	turbdiss[ip] += Sij[0][1]*Uy[ip];
	turbdiss[ip] += Sij[0][2]*Uz[ip];
	turbdiss[ip] += Sij[1][0]*Vx[ip];
	turbdiss[ip] += Sij[1][1]*Vy[ip];
	turbdiss[ip] += Sij[1][2]*Vz[ip];
	turbdiss[ip] += Sij[2][0]*Wx[ip];
	turbdiss[ip] += Sij[2][1]*Wy[ip];
	turbdiss[ip] += Sij[2][2]*Wz[ip];
	turbdiss[ip] *= 2*mu[ip];
	
	uprime2[ip]    = (U[ip]*U[ip] + V[ip]*V[ip] + W[ip]*W[ip])/3.0; 
	uvar[ip]       = (U[ip]*U[ip] + V[ip]*V[ip] + W[ip]*W[ip]); 
	uiprime2[ip]   = (Ux[ip]*Ux[ip] + Vy[ip]*Vy[ip] + Wz[ip]*Wz[ip])/3.0; 
	kineticEng[ip] = rho1[ip]*uvar[ip];

    }

    double taylor2_denom = 0;
    double meanRho = 0.0;
    double meanMu = 0.0;
    double meanTurbDiss = 0.0;
    double meanSOS = 0.0; 

    double turbMach = 0.0;
    double uprime = 0.0;
    double uiprime = 0.0;
    double urms = 0.0;
    double kolNu = 0.0;
    double meanKineticEng = 0.0;
 
    double rhoprime = 0.0;
    double dilprime = 0.0;
    double vortprime = 0.0;

    FOR_XYZ{
	meanRho += rho1[ip]/((double)(Nx*Ny*Nz));
	meanMu  += mu[ip]/((double)(Nx*Ny*Nz));
	meanTurbDiss += turbdiss[ip]/((double)(Nx*Ny*Nz));
	uprime += uprime2[ip]/((double)(Nx*Ny*Nz));
	uiprime += uiprime2[ip]/((double)(Nx*Ny*Nz));
	urms   += uvar[ip]/((double)(Nx*Ny*Nz));
	meanSOS += sos[ip]/((double)(Nx*Ny*Nz));
	taylor2_denom += (Ux[ip]*Ux[ip])/((double)(Nx*Ny*Nz));
	meanKineticEng += kineticEng[ip]/((double)(Nx*Ny*Nz))/2.0;
    }


    FOR_XYZ{
	rhoprime += ((rho1[ip] - meanRho)*(rho1[ip] - meanRho))/((double)(Nx*Ny*Nz));
	dilprime += ((Ux[ip] + Vy[ip] + Wz[ip])*(Ux[ip] + Vy[ip] + Wz[ip]))/((double)(Nx*Ny*Nz));

	double wx, wy, wz;
        wx = Wy[ip] - Vz[ip];
        wy = Uz[ip] - Wx[ip];
        wz = Vx[ip] - Uy[ip];
	vortprime += ((wx*wx + wy*wy + wz*wz))/((double)(Nx*Ny*Nz));
	
    }

    double enstrophy2Mu = meanMu*vortprime;

    rhoprime  = sqrt(rhoprime);
    dilprime  = sqrt(dilprime);
    vortprime = sqrt(vortprime);

    double meanNu = meanMu/meanRho; 
    meanTurbDiss /= meanRho;
    uprime = sqrt(uprime);
    uiprime = sqrt(uiprime);
    urms = sqrt(urms);

    turbMach = urms/meanSOS;

    double taylorMicro2 = uprime/uiprime; //this is the correct one for Samtaney paper
    double taylorMicro  = sqrt(15.0 * meanNu / meanTurbDiss)*urms;
    double taylorReyn  = uprime * taylorMicro / meanNu;
    double taylorReyn2 = urms*urms*sqrt(5.0/(meanNu*meanTurbDiss));
    double taylorReyn3  = uprime * taylorMicro2 / meanNu; //correct one for Samtaney paper
    double k0 = 8.0;
    double tau = sqrt(2.0*M_PI)/k0/uprime;


    //Kolmogorov Scale
    kolNu = pow(pow(meanNu,3.0)/meanTurbDiss,0.25);

    cout << " taylorReyn3: "  << taylorReyn3 << endl;
    cout << " turbMach: "  << turbMach << endl;
    cout << " taylorMicro2: " << taylorMicro2 << endl;
    cout << " meanTurbDiss: " << meanTurbDiss << endl;
    cout << " enstropy2Mu: " << enstrophy2Mu << endl;
    cout << " uprime: "      << uprime << endl;
    cout << " urms: "        << urms << endl;
    cout << " meanRho: "     << meanRho << endl;
    cout << " meanMu: "      << meanMu << endl;
    cout << " KolmogorovNu: " << kolNu << endl;
    cout << " rhoprime: " << rhoprime << endl;
    cout << " dilprime: " << dilprime << endl;
    cout << " vortprime: " << vortprime << endl;
    cout << " tau: " << tau << endl;
    cout << " t/tau: " << time/tau << endl;
    cout << " meanKineticEng: " << meanKineticEng << endl << endl;

        ofstream outfile;
        outfile.precision(17);
        string outputFileName;
        outputFileName = "turbdata.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
        outfile << time << " " << tau << " " << taylorReyn3 << " " << turbMach << " ";
	outfile << meanTurbDiss << " " << uprime << " " << kolNu << " ";
	outfile << rhoprime << " " << dilprime << " " << vortprime << " " << meanKineticEng << " ";
	outfile << endl;
        outfile.close();



};

void CSolver::calcTaylorGreenQuantities(){

    double enstrophySum = 0.0;
    double enstrophySum2Mu = 0.0;
    double turbDissSum = 0.0;
    double kineticEngSum = 0.0;
    double rhoprime = 0.0;
    double dilprime = 0.0;
    double meanMu   = 0.0;

    FOR_XYZ{
	double Sij[3][3];
	Sij[0][0] = 0.5*(Ux[ip] + Ux[ip]) - (1.0/3.0)*(Ux[ip] + Vy[ip] + Wz[ip]);
	Sij[1][1] = 0.5*(Vy[ip] + Vy[ip]) - (1.0/3.0)*(Ux[ip] + Vy[ip] + Wz[ip]);
	Sij[2][2] = 0.5*(Wz[ip] + Wz[ip]) - (1.0/3.0)*(Ux[ip] + Vy[ip] + Wz[ip]);
	Sij[0][1] = 0.5*(Uy[ip] + Vx[ip]);
	Sij[1][0] = 0.5*(Uy[ip] + Vx[ip]);
	Sij[0][2] = 0.5*(Uz[ip] + Wx[ip]);
	Sij[2][0] = 0.5*(Uz[ip] + Wx[ip]);
	Sij[1][2] = 0.5*(Vz[ip] + Wy[ip]);
	Sij[2][1] = 0.5*(Vz[ip] + Wy[ip]);

	turbdiss[ip]  = Sij[0][0]*Ux[ip];
	turbdiss[ip] += Sij[0][1]*Uy[ip];
	turbdiss[ip] += Sij[0][2]*Uz[ip];
	turbdiss[ip] += Sij[1][0]*Vx[ip];
	turbdiss[ip] += Sij[1][1]*Vy[ip];
	turbdiss[ip] += Sij[1][2]*Vz[ip];
	turbdiss[ip] += Sij[2][0]*Wx[ip];
	turbdiss[ip] += Sij[2][1]*Wy[ip];
	turbdiss[ip] += Sij[2][2]*Wz[ip];
	turbdiss[ip] *= 2*mu[ip];

        uvar[ip]       = (U[ip]*U[ip] + V[ip]*V[ip] + W[ip]*W[ip]);
	kineticEng[ip] = rho1[ip]*uvar[ip];

	double wx, wy, wz;
        wx = Wy[ip] - Vz[ip];
        wy = Uz[ip] - Wx[ip];
        wz = Vx[ip] - Uy[ip];
	enstrophySum += (1.0/2.0)*rho1[ip]*(wx*wx + wy*wy + wz*wz);
	enstrophySum2Mu += mu[ip]*rho1[ip]*(wx*wx + wy*wy + wz*wz);
	
	turbDissSum += turbdiss[ip];
	kineticEngSum += kineticEng[ip]/2.0;

	dilprime += ((Ux[ip] + Vy[ip] + Wz[ip])*(Ux[ip] + Vy[ip] + Wz[ip]))/((double)(Nx*Ny*Nz));

	meanMu += mu[ip]/((double)(Nx*Ny*Nz));
    }

    dilprime  = sqrt(dilprime);


    cout << " > kineticEngSum: " << kineticEngSum <<  endl;
    cout << " > turbDissSum:   " << turbDissSum <<  endl;
    cout << " > enstrophySum:  " << enstrophySum <<  endl;
    cout << " > enstrophySum*2mu :  " << enstrophySum2Mu <<  endl;
    cout << " > meanMu:        " << meanMu << endl;
    cout << " > dilprime:      " << dilprime << endl << endl;

        ofstream outfile;
        outfile.precision(17);
        string outputFileName;
        outputFileName = "taylorgreen.out";
        outfile.open(outputFileName, fstream::app);
        outfile.precision(17);
        outfile << time << " " << kineticEngSum << " " << turbDissSum << " " << enstrophySum << " " << enstrophySum2Mu << " " << meanMu << " " << dilprime;
	outfile << endl;
        outfile.close();



};

void CSolver::checkSolution(){
    if(timeStep%ts->checkStep == 0){
        t2Save = std::chrono::system_clock::now();
        cout << "-------------------------------------------------" << endl;
        cout << " Step = "<< timeStep << ", time = " << time << ", dt = " << ts->dt << endl;
        cout << "-------------------------------------------------" << endl;
        cout << "  Time since last timestep = " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2Save-t1Save).count()/(double)1000000000 << endl;
        getRange(rho1, "RHO", Nx, Ny, Nz);
        getRange(U, "U", Nx, Ny, Nz);
        getRange(V, "V", Nx, Ny, Nz);
        getRange(W, "W", Nx, Ny, Nz);
        getRange(p, "P", Nx, Ny, Nz);
        getRange(T, "T", Nx, Ny, Nz);
        getRange(mu, "mu", Nx, Ny, Nz);
        getRange(rhoE1, "RHOE", Nx, Ny, Nz);
        getRange(sos, "SOS", Nx, Ny, Nz);
        cout << endl;

        t1Save = std::chrono::system_clock::now();
    }
};


void CSolver::dumpSolution(){

    if(timeStep%ts->dumpStep == 0){
        cout << endl;
        cout << "> ===============" << endl;
        cout << ">  DUMPING FIELD " << endl;
        cout << "> ===============" << endl;

        ofstream outfile;
        outfile.precision(17);
        string outputFileName;

        outputFileName = "rho.out.";
        outputFileName.append(to_string(timeStep));
        outfile.open(outputFileName);
        outfile.precision(17);
        FOR_XYZ{
            outfile << rho1[ip] << " ";
        }
        outfile.close();

        outputFileName = "rhoU.out.";
        outputFileName.append(to_string(timeStep));
        outfile.open(outputFileName);
        outfile.precision(17);
	FOR_XYZ{
            outfile << rhoU1[ip] << " ";
        }
        outfile.close();

        outputFileName = "rhoV.out.";
        outputFileName.append(to_string(timeStep));
        outfile.open(outputFileName);
        outfile.precision(17);
	FOR_XYZ{
            outfile << rhoV1[ip] << " ";
        }
        outfile.close();

        outputFileName = "rhoW.out.";
        outputFileName.append(to_string(timeStep));
        outfile.open(outputFileName);
        outfile.precision(17);
	FOR_XYZ{
            outfile << rhoW1[ip] << " ";
        }
        outfile.close();

        outputFileName = "rhoE.out.";
        outputFileName.append(to_string(timeStep));
        outfile.open(outputFileName);
        outfile.precision(17);
	FOR_XYZ{
            outfile << rhoE1[ip] << " ";
        }
        outfile.close();

	int dumpSigma = 0;
	if(dumpSigma == 1){
	    outputFileName = "sigma.out";
            outputFileName.append(to_string(timeStep));
            outfile.open(outputFileName);
            outfile.precision(17);
            FOR_XYZ{
                 outfile << spg->sigma[ip] << " ";
            }
	    outfile.close();
	}

	int dumpVorticityMag = 0;
	if(dumpVorticityMag == 1){
	   double *vortMag = new double[N];
	   #pragma omp parallel for
	   FOR_XYZ{
		double wx, wy, wz;
		wx = Wy[ip] - Vz[ip];
		wy = Uz[ip] - Wx[ip];
		wz = Vx[ip] - Uy[ip];
		vortMag[ip] = sqrt(wx*wx + wy*wy + wz*wz);
	   }
	   string str = "vortMag.";
	   str.append(to_string(timeStep));
	   str.append(".out");
	   FILE *fp = fopen(str.c_str(), "wb");
	   fwrite((void *)vortMag, sizeof(double), Nx*Ny*Nx, fp);
	   fclose(fp);

           outputFileName = "vortMag.";
           outputFileName.append(to_string(timeStep));
	   outputFileName.append(".bov");
           outfile.open(outputFileName);
           outfile.precision(17);
	   outfile << "TIME: " << time << endl;	   
	   outfile << "DATA_FILE: " << str << endl;	   
	   outfile << "DATA_SIZE: " << Nx << " " << Ny << " " << Nz << endl;	   
	   outfile << "DATA_FORMAT: DOUBLE" << endl;	   
	   outfile << "VARIABLE: vortmag" << endl;	   
	   outfile << "DATA_ENDIAN: LITTLE" << endl;	   
	   outfile << "CENTERING: ZONAL" << endl;	   
	   outfile << "BRICK_ORIGIN: 0. 0. 0." << endl;	   
	   outfile << "BRICK_SIZE:" << Nx << ". " << Ny << ". " << Nz << ". " << endl;	   
	   outfile.close();

	}

    }
}

void CSolver::checkEnd(){

    if(time >= ts->maxTime){
	cout << "=================" << endl;
	cout << " HIT END OF TIME " << endl;
	cout << "=================" << endl;

	endFlag = true;
    }

    if(timeStep >= ts->maxTimeStep){
	cout << "=================" << endl;
	cout << " HIT END OF TIME " << endl;
	cout << "=================" << endl;

	endFlag = true;

    } 

    if(endFlag){
	dumpSolution();
    }

}

void CSolver::reportAll(){

   cout << "REPORT ALL" << endl;

   getRange(Ux, "Ux", Nx, Ny, Nz);
   getRange(Uxx, "Uxx", Nx, Ny, Nz);
   getRange(Uy, "Uy", Nx, Ny, Nz);
   getRange(Uyy, "Uyy", Nx, Ny, Nz);
   getRange(Uz, "Uz", Nx, Ny, Nz);
   getRange(Uzz, "Uzz", Nx, Ny, Nz);
   getRange(Uxy, "Uxy", Nx, Ny, Nz);
   getRange(Uyz, "Uyz", Nx, Ny, Nz);
   getRange(Uxz, "Uxz", Nx, Ny, Nz);
cout << " " << endl;
   getRange(Vx, "Vx", Nx, Ny, Nz);
   getRange(Vxx, "Vxx", Nx, Ny, Nz);
   getRange(Vy, "Vy", Nx, Ny, Nz);
   getRange(Vyy, "Vyy", Nx, Ny, Nz);
   getRange(Vz, "Vz", Nx, Ny, Nz);
   getRange(Vzz, "Vzz", Nx, Ny, Nz);
   getRange(Vxy, "Vxy", Nx, Ny, Nz);
   getRange(Vyz, "Vyz", Nx, Ny, Nz);
   getRange(Vxz, "Vxz", Nx, Ny, Nz);
cout << " " << endl;
   getRange(Wx, "Wx", Nx, Ny, Nz);
   getRange(Wxx, "Wxx", Nx, Ny, Nz);
   getRange(Wy, "Wy", Nx, Ny, Nz);
   getRange(Wyy, "Wyy", Nx, Ny, Nz);
   getRange(Wz, "Wz", Nx, Ny, Nz);
   getRange(Wzz, "Wzz", Nx, Ny, Nz);
   getRange(Wxy, "Wxy", Nx, Ny, Nz);
   getRange(Wyz, "Wyz", Nx, Ny, Nz);
   getRange(Wxz, "Wxz", Nx, Ny, Nz);
cout << " " << endl;
   getRange(Tx, "Ux", Nx, Ny, Nz);
   getRange(Txx, "Uxx", Nx, Ny, Nz);
   getRange(Ty, "Uy", Nx, Ny, Nz);
   getRange(Tyy, "Uyy", Nx, Ny, Nz);
   getRange(Tz, "Uz", Nx, Ny, Nz);
   getRange(Tzz, "Uzz", Nx, Ny, Nz);
cout << " " << endl;
   getRange(contEulerX, "contEulerX", Nx, Ny, Nz);
   getRange(contEulerY, "contEulerY", Nx, Ny, Nz);
   getRange(contEulerZ, "contEulerZ", Nx, Ny, Nz);
cout << " " << endl;
   getRange(momXEulerX, "momXEulerX", Nx, Ny, Nz);
   getRange(momXEulerY, "momXEulerY", Nx, Ny, Nz);
   getRange(momXEulerZ, "momXEulerZ", Nx, Ny, Nz);
cout << " " << endl;
   getRange(momYEulerX, "momYEulerX", Nx, Ny, Nz);
   getRange(momYEulerY, "momYEulerY", Nx, Ny, Nz);
   getRange(momYEulerZ, "momYEulerZ", Nx, Ny, Nz);
cout << " " << endl;
   getRange(momZEulerX, "momZEulerX", Nx, Ny, Nz);
   getRange(momZEulerY, "momZEulerY", Nx, Ny, Nz);
   getRange(momZEulerZ, "momZEulerZ", Nx, Ny, Nz);
cout << " " << endl;
   getRange(engyEulerX, "engyEulerX", Nx, Ny, Nz);
   getRange(engyEulerY, "engyEulerY", Nx, Ny, Nz);
   getRange(engyEulerZ, "engyEulerZ", Nx, Ny, Nz);
cout << " " << endl;
   getRange(rho1, "rho1", Nx, Ny, Nz);
   getRange(rhok, "rhok", Nx, Ny, Nz);
   getRange(rhok2, "rhok2", Nx, Ny, Nz);
   getRange(rho2, "rho2", Nx, Ny, Nz);
cout << " " << endl;
   getRange(rhoU1, "rhoU1", Nx, Ny, Nz);
   getRange(rhoUk, "rhoUk", Nx, Ny, Nz);
   getRange(rhoUk2, "rhoUk2", Nx, Ny, Nz);
   getRange(rhoU2, "rhoU2", Nx, Ny, Nz);
cout << " " << endl;
   getRange(rhoV1, "rhoV1", Nx, Ny, Nz);
   getRange(rhoVk, "rhoVk", Nx, Ny, Nz);
   getRange(rhoVk2, "rhoVk2", Nx, Ny, Nz);
   getRange(rhoV2, "rhoV2", Nx, Ny, Nz);
cout << " " << endl;
   getRange(rhoW1, "rhoW1", Nx, Ny, Nz);
   getRange(rhoWk, "rhoWk", Nx, Ny, Nz);
   getRange(rhoWk2, "rhoWk2", Nx, Ny, Nz);
   getRange(rhoW2, "rhoW2", Nx, Ny, Nz);
cout << " " << endl;
   getRange(rhoE1, "rhoE1", Nx, Ny, Nz);
   getRange(rhoEk, "rhoEk", Nx, Ny, Nz);
   getRange(rhoEk2, "rhoEk2", Nx, Ny, Nz);
   getRange(rhoE2, "rhoE2", Nx, Ny, Nz);
cout << " " << endl;
   getRange(p, "p", Nx, Ny, Nz);
   getRange(U, "U", Nx, Ny, Nz);
   getRange(V, "V", Nx, Ny, Nz);
   getRange(W, "W", Nx, Ny, Nz);
   getRange(T, "T", Nx, Ny, Nz);
   getRange(mu, "mu", Nx, Ny, Nz);
   getRange(Amu, "Amu", Nx, Ny, Nz);
   getRange(sos, "sos", Nx, Ny, Nz);
cout << " " << endl;


}
