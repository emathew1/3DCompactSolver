#ifndef _CRK4H_
#define _CRK4H_

#include "Utils.hpp"
#include "AbstractRK.hpp"

class RK4:public AbstractRK{

    public:

	double a1, a2, a3, a4;
	double b1, b2, b3;

	RK4(AbstractCSolver *cs){
	    this->cs = cs;
	
	    a1 = 1.0/6.0;
	    a2 = 1.0/3.0;
	    a3 = 1.0/3.0;
	    a4 = 1.0/6.0;

	    b1 = 1.0/2.0;
	    b2 = 1.0/2.0;
	    b3 = 1.0;

	}


        void executeSolverLoop();
        void updateConservedData();

};


void RK4::updateConservedData(){

    if(cs->useTiming) tic();

    int Nx, Ny, Nz;
    Nx = cs->dom->Nx;
    Ny = cs->dom->Ny;
    Nz = cs->dom->Nz;


    if(cs->rkStep == 1){

        #pragma omp parallel for
        FOR_XYZ{
            //Add to final solution 
            cs->rho2[ip]  = cs->rho1[ip]  + a1*cs->rhok2[ip];
            cs->rhoU2[ip] = cs->rhoU1[ip] + a1*cs->rhoUk2[ip];
            cs->rhoV2[ip] = cs->rhoV1[ip] + a1*cs->rhoVk2[ip];
            cs->rhoW2[ip] = cs->rhoW1[ip] + a1*cs->rhoWk2[ip];
            cs->rhoE2[ip] = cs->rhoE1[ip] + a1*cs->rhoEk2[ip];
            
            //Calculate intermediate solution
            cs->rhok[ip]  = cs->rho1[ip]  + b1*cs->rhok2[ip]; 
            cs->rhoUk[ip] = cs->rhoU1[ip] + b1*cs->rhoUk2[ip];
            cs->rhoVk[ip] = cs->rhoV1[ip] + b1*cs->rhoVk2[ip];
            cs->rhoWk[ip] = cs->rhoW1[ip] + b1*cs->rhoWk2[ip];
            cs->rhoEk[ip] = cs->rhoE1[ip] + b1*cs->rhoEk2[ip];
        }
    }else if(cs->rkStep == 2){

        #pragma omp parallel for
        FOR_XYZ{
            //Add to final solution
            cs->rho2[ip]  += a2*cs->rhok2[ip];
            cs->rhoU2[ip] += a2*cs->rhoUk2[ip];
            cs->rhoV2[ip] += a2*cs->rhoVk2[ip];
            cs->rhoW2[ip] += a2*cs->rhoWk2[ip];
            cs->rhoE2[ip] += a2*cs->rhoEk2[ip];
            
            //Calculate intermediate solution
            cs->rhok[ip]  = cs->rho1[ip]  + b2*cs->rhok2[ip]; 
            cs->rhoUk[ip] = cs->rhoU1[ip] + b2*cs->rhoUk2[ip];
            cs->rhoVk[ip] = cs->rhoV1[ip] + b2*cs->rhoVk2[ip];
            cs->rhoWk[ip] = cs->rhoW1[ip] + b2*cs->rhoWk2[ip];
            cs->rhoEk[ip] = cs->rhoE1[ip] + b2*cs->rhoEk2[ip];
        }

    }else if(cs->rkStep == 3){

        #pragma omp parallel for
        FOR_XYZ{
            //Add to final solution
            cs->rho2[ip]  += a3*cs->rhok2[ip];
            cs->rhoU2[ip] += a3*cs->rhoUk2[ip];
            cs->rhoV2[ip] += a3*cs->rhoVk2[ip];
            cs->rhoW2[ip] += a3*cs->rhoWk2[ip];
            cs->rhoE2[ip] += a3*cs->rhoEk2[ip];

            //Calculate intermediate solution
            cs->rhok[ip]  = cs->rho1[ip]  + b3*cs->rhok2[ip];
            cs->rhoUk[ip] = cs->rhoU1[ip] + b3*cs->rhoUk2[ip];
            cs->rhoVk[ip] = cs->rhoV1[ip] + b3*cs->rhoVk2[ip];
            cs->rhoWk[ip] = cs->rhoW1[ip] + b3*cs->rhoWk2[ip];
            cs->rhoEk[ip] = cs->rhoE1[ip] + b3*cs->rhoEk2[ip];
        }

    }else if(cs->rkStep == 4){

        #pragma omp parallel for
        FOR_XYZ{
            //Add to final solution
            cs->rho2[ip]  += a4*cs->rhok2[ip];
            cs->rhoU2[ip] += a4*cs->rhoUk2[ip];
            cs->rhoV2[ip] += a4*cs->rhoVk2[ip];
            cs->rhoW2[ip] += a4*cs->rhoWk2[ip];
            cs->rhoE2[ip] += a4*cs->rhoEk2[ip];
        }
    }

    if(cs->useTiming){
        cout << " > updateCons Timing: ";
        toc();
    }
}

void RK4::executeSolverLoop(){

    cs->setInitialConditions();

    while(cs->endFlag == false){

	cs->rkLast = false;
	
	cs->preStep();

	//Step 1
	cs->rkStep = 1;
	
	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();

	//Step 2
	cs->rkStep = 2;

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();

	//Step 3
	cs->rkStep = 3;

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();

	//Step 4
	cs->rkStep = 4;
	cs->rkLast = true;	

	cs->preSubStep();
	cs->solveEqnSet();
	cs->postSubStep();
	updateConservedData();
	cs->updateData();

	cs->postStep();	

    }
}
        
#endif
